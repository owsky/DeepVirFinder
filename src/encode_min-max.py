from collections import Counter

def min_max_normalization(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    
    # Get the frequencies of k-mers
    frequencies = list(kmer_counts.values())
    
    min_freq = min(frequencies)
    max_freq = max(frequencies)
    
    normalized_frequencies = [(freq - min_freq) / (max_freq - min_freq) for freq in frequencies]
    
    return normalized_frequencies

raw_sequence = "ATCGATCGATCG"
k = 3
normalized_frequencies = min_max_normalization(raw_sequence, k)
print(normalized_frequencies)
