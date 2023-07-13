from collections import Counter
import os
import sys
from kmer_gen import gen_all_kmers

def normalize_kmers(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmers = [kmer for kmer in kmers if set(kmer).issubset({'A', 'C', 'G', 'T'})]
    kmer_counts = Counter(kmers)

    all_kmers = gen_all_kmers(k)
    # Normalize the k-mer counts
    normalized_kmers = [kmer_counts.get(kmer, 0) for kmer in all_kmers]

    return normalized_kmers


def normalize_sequence(sequence: str) -> list:
    k = os.getenv("KMER")
    if k is None:
        sys.stderr.write("Missing KMER env variable\n")
        sys.exit(1)
    else:
        k = int(k)
    return normalize_kmers(sequence, k)