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
    kmers = list(kmer_counts.keys())
    counts = np.array(list(kmer_counts.values())).reshape(-1, 1)

    # Apply quantile normalization using sklearn's quantile_transform
    normalized_counts = quantile_transform(counts, output_distribution='normal', random_state=0)

    normalized_kmer_counts = dict(zip(kmers, normalized_counts.flatten()))
    return normalized_kmer_counts

# Example usage:
# Create a sample genomic sequence
genomic_sequence = 'ACGTACGTACGT'

# Set the desired k-mer length
k = 3

# Normalize the k-mer counts
normalized_kmer_counts = normalize_kmer_counts(genomic_sequence, k)

# Display the original and normalized k-mer counts
print("Original k-mer counts:")
print(calculate_kmer_counts(genomic_sequence, k))
print("\nNormalized k-mer counts:")
print(normalized_kmer_counts)
