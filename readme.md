# nearproteins

Store for finding similar amino acid sequences using Locality Sensitive Hashing
and Approximate Nearest Neighbors.

## Installation

Python, with dependencies:

- NearPy, using [this](https://github.com/pixelogik/NearPy/tree/2d05bf38d8dc52cb765534094cb5006c9ed622b6) ref on GitHub.
- BioPython

Redis

## Instructions

1. Load proteins into database

```sh
$ ./load-proteins < data/proteins.fasta
```

2. Query database

```sh
$ ./query-proteins < data/proteins.fasta # returns JSON for each record
```

## Python API

Very basic. I plan to add more configuration.

### Loading data into store

```python

import nearproteins

store = nearproteins.SimilarStringStore()

store.engine.clean_all_buckets()

records = SeqIO.parse(handle, 'fasta')

for record in records:
    store.add(str(record.seq), record.id)
```

### Retrieving records from store

```python
# returns array of vectors, match IDs, similarities
results = store.get(str(record.seq)) 
```

## Use as a server

You can query and add records to the database using simple sockets.

```
./server # start the server, listens on port 1234
```

In another window...

```
nc 127.0.0.1 1234 # connect
SET 1 AUSTIN
SET 2 BOSTON
GET AUSTIN
```
