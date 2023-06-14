from collections import Counter
import itertools


def create_kmer_index(k):
    # Create a dictionary of all possible k-mers of length k
    kmers = [''.join(p) for p in itertools.product('ACGT', repeat=k)]
    kmer_index = {kmer: i for i, kmer in enumerate(kmers)}
    return kmer_index


def get_freqs(sequence: str) -> list:
    k = 4
    kmer_index = create_kmer_index(k)
    vector = [0] * len(kmer_index)
    kmer_counts = Counter()
    
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        
        if set(kmer).issubset({'A', 'C', 'G', 'T'}):
            kmer_counts[kmer] += 1
            index = kmer_index[kmer]
            vector[index] = kmer_counts[kmer]
    
    return vector



def normalize_sequence(sequence: str) -> list:
    k = 4
    kmer_index = create_kmer_index(k)
    return get_freqs(sequence, kmer_index, k)
