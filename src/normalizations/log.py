from collections import Counter
import math
import os
import sys

from kmer_gen import gen_all_kmers


def log_normalization(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)

    all_kmers = gen_all_kmers(k)
    
    # Get the frequencies of all possible k-mers (including those with zero counts)
    frequencies = [kmer_counts.get(kmer, 0) for kmer in all_kmers]

    # Apply log transformation to the frequencies
    log_frequencies = [math.log(freq) if freq != 0 else 0 for freq in frequencies]

    return log_frequencies


def normalize_sequence(sequence: str) -> list:
    k = os.getenv("KMER")
    if k is None:
        sys.stderr.write("Missing KMER env variable\n")
        sys.exit(1)
    else:
        k = int(k)
    return log_normalization(sequence, k)