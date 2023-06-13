from itertools import product


alphabet = 'ACGT'


def gen_all_kmers(k):
    return [''.join(x) for x in product(alphabet, repeat=k)]