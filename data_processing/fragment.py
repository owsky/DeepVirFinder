import random
import os
import sys

def generate_random_fragment(sequence):
    fragment_size = random.randint(150, 1000)
    return sequence[:fragment_size]

def read_fasta_file(file_path):
    sequences = {}
    current_sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ""
            else:
                sequences[current_sequence] += line
    return sequences

def write_fasta_file(sequences, output_path):
    with open(output_path, 'w') as file:
        for sequence_name, sequence in sequences.items():
            file.write(">" + sequence_name + "\n")
            file.write(sequence + "\n")

def generate_random_fragment_fasta(input_file_path, output_file_path):
    sequences = read_fasta_file(input_file_path)
    fragmented_sequences = {}
    for sequence_name, sequence in sequences.items():
        fragment = generate_random_fragment(sequence)
        fragmented_sequences[sequence_name] = fragment
    write_fasta_file(fragmented_sequences, output_file_path)

# Usage example
base_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "test")
input_path = os.path.join(base_path, "host_test.fa")
output_path = os.path.join(base_path, "host_test_fragmented.fa")
generate_random_fragment_fasta(input_path, output_path)
