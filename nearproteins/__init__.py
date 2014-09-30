#!/usr/bin/env python

from itertools import product
from collections import defaultdict
from redis import Redis
import sys
import json

from Bio import SeqIO

from nearpy import Engine
from nearpy import hashes
from nearpy.storage import RedisStorage
from nearpy.filters import DistanceThresholdFilter

import numpy as np

# The distance ratio of an ANN y is it's distance to the minimal hypersphere
# around the query vector x, that contains all exact nearest neighbours n,
# clamped to zero and normalized with this hypersphere's radius

class FeatureGenerator:

    def __init__(self, **kwargs):
        ''' '''

        self.alphabet = [ 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K',
                     'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', ]

        assert len(self.alphabet) == 20
        self.k = kwargs['K']

        self.feature_space = list(''.join(i) for i in product(self.alphabet,
            repeat=self.k))

        self.n_features = len(self.feature_space)


    def shingles(self, s, k):
        ''' return shingles of a given string given a k-mer size k '''
        return [ s[i : i + k ] for i in range(0, len(s) - k + 1) ]


    def vectorize(self, s):
        ''' convert shingles to features vector '''

        d = defaultdict(lambda: 0)

        for i in s:
            d[i] += 1

        return np.array([ d[i] for i in self.feature_space ])


    def transform(self, str):
        return self.vectorize(self.shingles(str, self.k))




class SimilarStringStore:

    def __init__(self, **kwargs):

        defaults = { 'seed': 42, 'K': 2, 'P': 100, 'MAX_DIST': 100 }

        defaults.update(kwargs)
        self.config = defaults

        self.transformer = FeatureGenerator(K = self.config['K'])

        self.hasher = hashes.RandomBinaryProjections('rbp', self.config['P'],
                rand_seed=42)

        self.backend = RedisStorage(Redis(host='localhost', port=6379, db=0))

        self.filters = [ DistanceThresholdFilter(self.config['MAX_DIST']) ]

        self.engine = Engine(self.transformer.n_features,
                     lshashes=[self.hasher],
                     vector_filters = self.filters,
                     storage=self.backend
                     )

    def vectorize(self, s):
        return self.engine.transformer.transform(s)

    def add(self, s, id):
        ''' add a string to index '''
        vector = self.transformer.transform(s)
        self.engine.store_vector(vector, str(id))
        return vector

    def query(self, s):
        ''' query index '''
        vector = self.transformer.transform(s)
        neighbours = self.engine.neighbours(vector)
        return neighbours

    def remove(self, id):
        ''' remove a string from the index '''
        pass
