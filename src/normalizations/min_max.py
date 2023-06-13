from collections import Counter
from kmer_gen import gen_all_kmers


def min_max_normalization(sequence, k):
    # Generate all possible k-mers
    all_kmers = gen_all_kmers(k)
    
    # Count the occurrences of each k-mer
    kmer_counts = Counter(all_kmers)
    
    # Get the frequencies of k-mers in the sequence
    sequence_kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_frequencies = [kmer_counts[kmer] for kmer in sequence_kmers]
    
    # Normalize the frequencies
    min_freq = min(kmer_frequencies)
    max_freq = max(kmer_frequencies)
    
    normalized_frequencies = [(freq - min_freq) / (max_freq - min_freq) if (max_freq - min_freq) != 0 else 0 for freq in kmer_frequencies]
    
    return normalized_frequencies

def normalize_sequence(sequence: str) -> list:
    k = 4
    return min_max_normalization(sequence, k)