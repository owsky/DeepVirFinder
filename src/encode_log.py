from collections import Counter
import math

def log_normalization(sequence, k):
    # Count the occurrences of each k-mer
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    
    # Get the frequencies of k-mers
    frequencies = list(kmer_counts.values())
    
    # Apply log transformation to the frequencies
    log_frequencies = [math.log(freq) for freq in frequencies]
    
    return log_frequencies

raw_sequence = "ATCGATCGATCG"
k = 3
log_normalized_frequencies = log_normalization(raw_sequence, k)
print(log_normalized_frequencies)
