from collections import Counter
import os
import sys
from kmer_gen import gen_all_kmers
import statistics

def normalize_with_zscore(data):
    mean = statistics.mean(data)
    std_dev = statistics.stdev(data)
    normalized_data = [(x - mean) / std_dev for x in data]
    return normalized_data

def normalize_kmers(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmers = [kmer for kmer in kmers if set(kmer).issubset({'A', 'C', 'G', 'T'})]
    kmer_counts = Counter(kmers)

    all_kmers = gen_all_kmers(k)

    # Normalize the k-mer counts
    all_kmer_counts = [kmer_counts.get(kmer, 0) for kmer in all_kmers]
    norm = normalize_with_zscore(all_kmer_counts)
    return norm


def normalize_sequence(sequence: str) -> list:
    k = os.getenv("KMER")
    if k is None:
        sys.stderr.write("Missing KMER env variable\n")
        sys.exit(1)
    else:
        k = int(k)
    normalized = normalize_kmers(sequence, k)
    return normalized
