from collections import Counter
import os
import sys
from kmer_gen import gen_all_kmers


def min_max_normalization(sequence, k):
    # Generate all possible k-mers
    all_kmers = gen_all_kmers(k)
    
    # Count the occurrences of each k-mer in the sequence
    kmer_counts = Counter(sequence[i:i+k] for i in range(len(sequence)-k+1))
    
    # Get the frequencies of k-mers in the sequence
    kmer_frequencies = [kmer_counts[kmer] for kmer in all_kmers]
    
    # Normalize the frequencies
    min_freq = min(kmer_frequencies)
    max_freq = max(kmer_frequencies)
    
    normalized_frequencies = [(freq - min_freq) / (max_freq - min_freq) if (max_freq - min_freq) != 0 else 0 for freq in kmer_frequencies]
    
    return normalized_frequencies

def normalize_sequence(sequence: str) -> list:
    k = os.getenv("KMER")
    if k is None:
        sys.stderr.write("Missing KMER env variable\n")
        sys.exit(1)
    else:
        k = int(k)
        print(k)
    return min_max_normalization(sequence, k)