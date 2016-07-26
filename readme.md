# nearproteins

Store for finding similar amino acid sequences using Locality Sensitive Hashing
and Approximate Nearest Neighbors.

## Installation

Python, with dependencies:

- Annoy
- BioPython

```
pip install -r requirements.txt
```

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
from Bio import SeqIO
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

```sh
$ ./server # start the server, listens on port 1234
```

In another window...

```sh
$ nc 127.0.0.1 1234 # connect
SET 1 AUSTIN
SET 2 BOSTON
GET AUSTIN
{"1": 0.0}
GET BOSTON
{"2": 0.0}
```

You can use this to build a simple client in another language such as Ruby

```ruby
require 'dna'
require 'socket'
require 'json'

HOSTNAME = '127.0.0.1'
PORT = '1234'

socket = TCPSocket.open HOSTNAME, PORT

File.open('proteins.fasta') do |handle|
  records = Dna.new handle, :format => :fasta

  records.each do |record|
    socket.puts "GET #{record.sequence}"
    resp = JSON.parse(socket.gets)
    p resp
  end
end
```
