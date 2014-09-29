#!/usr/bin/env python

from itertools import product
from collections import defaultdict
from redis import Redis
import sys
import json

from Bio import SeqIO

from nearpy import Engine
from nearpy.hashes import RandomBinaryProjections
from nearpy.storage import RedisStorage
from nearpy.filters import DistanceThresholdFilter

import numpy as np

# 20 amino acids
ALPHABET = [ 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K',
             'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', ]

assert len(ALPHABET) == 20

K = 2 # k-mer size for shingles
P = 100 # number of Random Binary Projections
MAX_DIST = 20 # maximum distance threshold

# The distance ratio of an ANN y is it's distance to the minimal hypersphere
# around the query vector x, that contains all exact nearest neighbours n,
# clamped to zero and normalized with this hypersphere's radius

TEST = False # test on a little bit of sequences

def shingles(s, k):
    ''' return shingles of a given string given a k-mer size k '''
    return [ s[i : i + k ] for i in range(0, len(s) - k + 1) ]


def to_features(s, f):
    ''' convert shingles to features vector '''
    d = defaultdict(lambda: 0)
    for i in s:
        d[i] += 1
    # automatically filters out features that aren't in ALPHABET!
    return np.array([ d[i] for i in f ])

# array of all possible features
features = list(''.join(i) for i in product(ALPHABET, repeat=2))

rbp = RandomBinaryProjections('rbp', P)

redis_storage = RedisStorage(Redis(host='localhost', port=6379, db=0))

lsh = Engine(len(features),
             lshashes=[rbp],
             vector_filters = [DistanceThresholdFilter(MAX_DIST)],
             storage=redis_storage
             )

# train
with open('proteins.fasta') as handle:
    records = SeqIO.parse(handle, 'fasta')

    for i, record in enumerate(records):
        seq = str(record.seq)
        feat = to_features(shingles(seq, K), features)

        lsh.store_vector(feat, record.id)

        if i % 10000 == 0:
            print >> sys.stderr, 'training: %s' % (i)

        if TEST:
            if i > 10000:
                break

# predict

neighbor_counts = []

with open('proteins.fasta') as handle:
    records = SeqIO.parse(handle, 'fasta')

    for i, record in enumerate(records):
        seq = str(record.seq)
        feat = to_features(shingles(seq, K), features)

        neighbours = lsh.neighbours(feat)

        hits = [ n[1] for n in neighbours ]

        if i % 10000 == 0:
            print >> sys.stderr, 'predicting: %s' % (i)

        print json.dumps({ 'query': record.id, 'hits': hits })

        if TEST:
            if i > 10000:
                break
