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
        self.feature_space = list(''.join(i) for i in product(self.alphabet, repeat=self.k))

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

    def __init__(self):

        self.config = {}
        self.config['K'] = 2 # k-mer size for shingles
        self.config['P'] = 100 # number of Random Binary Projections
        self.config['MAX_DIST'] = 20 # maximum distance threshold

        self.transformer = FeatureGenerator(K = self.config['K'])

        self.hasher = RandomBinaryProjections('rbp', self.config['P'])
        self.backend = RedisStorage(Redis(host='localhost', port=6379, db=0))

        self.filters = [ DistanceThresholdFilter(self.config['MAX_DIST']) ]

        self.engine = Engine(self.transformer.n_features,
                     lshashes=[self.hasher],
                     vector_filters = self.filters,
                     storage=self.backend
                     )

    def add(self, id, str):
        ''' add a string to index '''

        vector = self.transformer.transform(str)
        self.engine.store_vector(vector, id)
        return vector


    def query(self, str):
        ''' query index '''
        vector = self.transformer.transform(str)
        neighbours = self.engine.neighbours(vector)
        return neighbours
