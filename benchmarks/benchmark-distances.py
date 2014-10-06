from Bio import SeqIO
import numpy as np

import nearproteins as nearp
import redis
import logging
import timeit

# # within distance benchmark

# cluster a group of proteins and measure the intra-family and inter-family
# distances as a function of various parameters:
#
# - minimum distance.
# - number of hashers.
# - length of hash (random projections).
# - number of proteins


def benchmark(records, params):

    all_distances = []
    neighbor_counts = []
    add_times = []
    query_times = []

    # build up database
    logging.info('clearing store')

    store = nearp.SimilarStringStore(n_hashers = params['n_hashers'],
                                     P = params['n_projections'] )
    store.engine.storage.redis_object.flushall()

    logging.info('inserting records')
    for record in records:
        start = timeit.timeit()
        store.add(str(record.seq), record.id)
        stop = timeit.timeit()

    logging.info('finding neighbors')
    for record in records:
        start = timeit.timeit()
        res = store.query(str(record.seq))
        stop = timeit.timeit()

        query_times.append(stop - start)

        neighbors = [ i[1] for i in res ]
        distances = [ i[2] for i in res ]

        all_distances.append(distances)
        neighbor_counts.append(len(neighbors))

        add_times.append(stop - start)

    return { 'distances': distances,
             'neighbor_counts': neighbor_counts,
             'add_times': add_times,
             'query_times': query_times }



def main():

    logging.basicConfig(level=logging.INFO, logfile='/dev/stderr')

    with open('proteins.fasta') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))

    for n_hashers in [1, 2, 10, 30]: # number of hashes
        for n_projections in [1, 5, 10, 20]: # hash length
            params = { 'n_hashers': n_hashers,
                       'n_projections': n_projections }

            logging.info('params: %s' % params)

            res = benchmark(records, params)

            logging.info('distances: %.2f +/- %.2f' % (np.median(res['distances']),
                np.std(res['distances'])))

            logging.info('neighbors: %.2f +/- %.2f' %
                    (np.median(res['neighbor_counts']),
                        np.std(res['neighbor_counts'])))

            logging.info('insert time: %.2g +/- %.2g' %
                    (np.median(res['add_times']), np.std(res['add_times'])))

            logging.info('query time: %.2g +/- %.2g' %
                    (np.median(res['query_times']), np.std(res['query_times'])))




if __name__ == '__main__':
    main()
