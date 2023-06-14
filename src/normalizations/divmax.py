from collections import Counter
import os
import sys
import numpy as np
from normalizations.kmer_gen import gen_all_kmers


def normalize_kmers(sequence, k):
    # Count the occurrence of each k-mer in the sequence
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    
    # Find the maximum count of k-mers in the sequence
    max_count = max(kmer_counts.values())
    
    # Normalize the k-mer counts
    normalized_kmers = np.zeros(4**k)
    for i, kmer in enumerate(gen_all_kmers(k)):
        normalized_kmers[i] = kmer_counts[kmer] / max_count
    
    return normalized_kmers


def normalize_sequence(sequence: str) -> list:
    k = os.getenv("KMER")
    if k is None:
        sys.stderr.write("Missing KMER env variable\n")
        sys.exit(1)
    else:
        k = int(k)
        print(k)
    return normalize_kmers(sequence, k)