import os
import sys
import numpy as np
from Bio.Seq import Seq

def encodeSeq(sequence):
    # Encode each base pair of the sequence
    encoded_sequence = []
    for base in sequence:
        if base == 'A':
            encoded_sequence.append([1, 0, 0, 0])
        elif base == 'C':
            encoded_sequence.append([0, 1, 0, 0])
        elif base == 'G':
            encoded_sequence.append([0, 0, 1, 0])
        elif base == 'T':
            encoded_sequence.append([0, 0, 0, 1])
        else:
            encoded_sequence.append([1/4, 1/4, 1/4, 1/4])
    
    return encoded_sequence

def encodeFastaFile(fasta_file, direction):
    sequences = []
    with open(fasta_file, 'r') as file:
        header = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                if header is not None and sequence:
                    if direction == "bw":
                        sequence = Seq(sequence).reverse_complement()
                    sequences.append((header, encodeSeq(sequence)))  # Append previous sequence and its header
                header = line[1:]  # Extract the sequence header (without '>')
                sequence = ''
            else:  # Sequence line
                sequence += line.upper()  # Append the sequence to the existing sequence
        
        if header is not None and sequence:
            if direction == "bw":
                sequence = Seq(sequence).reverse_complement()
            sequences.append((header, encodeSeq(sequence)))  # Append the last sequence and its header
    
    return sequences

def main():
    input_file = sys.argv[1]
    input_path = os.path.dirname(input_file)
    input_file_name = os.path.basename(input_file).split(".")[0]
    os.makedirs(os.path.join(input_path, "encoded"), exist_ok=True)
    direction = sys.argv[2]
    output_path = os.path.join(input_path, "encoded", input_file_name + "_" + direction + ".npy")
    encoded_sequences = encodeFastaFile(input_file, direction)

    # Convert the encoded sequences to numpy array
    encoded_array = np.array([encoded_seq for _, encoded_seq in encoded_sequences], dtype=object)

    # Save the encoded sequences to disk using numpy.save
    np.save(output_path, encoded_array)
    print(f"Encoded sequences saved to '{output_path}'.")

if __name__ == "__main__":
    main()