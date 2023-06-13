import numpy as np
from sklearn.preprocessing import quantile_transform


def calculate_kmer_counts(sequence, k):
    kmer_counts = {}
    n = len(sequence)
    for i in range(n - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts


def normalize_kmer_counts(sequence, k):
    kmer_counts = calculate_kmer_counts(sequence, k)
    # kmers = list(kmer_counts.keys())
    counts = np.array(list(kmer_counts.values())).reshape(-1, 1)

    n_quantiles = min(counts.shape[0], 1000)
    
    # Apply quantile normalization using sklearn's quantile_transform
    normalized_counts = quantile_transform(counts, output_distribution='normal', random_state=0, n_quantiles=n_quantiles)

    return normalized_counts


def normalize_sequence(sequence: str) -> list:
    k = 4
    return normalize_kmer_counts(sequence, k)