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
