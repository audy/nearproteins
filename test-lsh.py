#!/usr/bin/env python

from sklearn.neighbors import LSHForest
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.pipeline import Pipeline
from Bio import SeqIO

import sys

# build pipeline
vectorizer = TfidfVectorizer(ngram_range=(2, 2), analyzer='char')
lsh = LSHForest()

# load sequences
with open(sys.argv[1]) as handle:
    records = [ ( r.id, str(r.seq) ) for r in SeqIO.parse(handle, 'fasta') ]

print('loaded %s records' % len(records))

# transform sequences
print('transforming...')
X = vectorizer.fit_transform(i[1] for i in records)

# fit LSH
print('fitting...')
lsh.fit(X)

# do some queries
print('querying...')

neighbors = lsh.kneighbors(X[0], 10)


ids = list(neighbors[1][0])

query = records[0][0]

print('query = %s' %  query)
print('hits:')

print ids
for id_ in ids:
    print('- %s' % records[id_][0])
