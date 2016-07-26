#!/usr/bin/env python

from collections import defaultdict
from itertools import product
import json
import random
import sys

from annoy import AnnoyIndex

from Bio import SeqIO

import numpy as np

class FeatureGenerator:

    def __init__(self, k=2):
        ''' '''

        self.k = k

        self.alphabet = [ 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K',
                     'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', ]

        assert len(self.alphabet) == 20

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

        # convert to counts in feature space
        vector = np.array([ d[i] for i in self.feature_space ])

        return vector


    def transform(self, str):
        return self.vectorize(self.shingles(str, self.k))


class SimilarStringStore:

    def __init__(self, **kwargs):

        self.transformer = FeatureGenerator(k=1)

        print(self.transformer.n_features)

        self.store = AnnoyIndex(self.transformer.n_features)

    def vectorize(self, s):
        return self.transformer.transform(s)

    def add(self, id, s):
        ''' add a string to index '''

        vector = self.transformer.transform(s)
        self.store.add_item(int(id), vector)
        return vector

    def build(self):
        self.store.build(500)

    def save(self, filename='store.knn'):
        self.store.save(filename)

    def build_and_save(self, filename='store.knn'):
        self.build()
        self.save(filename)

    def load(self, filename='store.knn'):
        self.store.load(filename)


    def query(self, s):
        ''' query index '''
        vector = self.transformer.transform(s)
        neighbors = self.store.get_nns_by_vector(vector, 40)
        return neighbors


    def remove(self, id):
        ''' remove a string from the index '''
        pass
