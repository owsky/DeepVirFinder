from sklearn.preprocessing import StandardScaler
import numpy as np
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import os

fasta_file_path = sys.argv[1]
contig_type = sys.argv[2]
out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_z-score")
os.makedirs(out_dir, exist_ok=True)
scaler = StandardScaler()

def encodeSeq(seq) : 
    seq_code = list()
    for pos in range(len(seq)) :
        letter = seq[pos]
        if letter in ['A', 'a'] :
            code = [1,0,0,0]
        elif letter in ['C', 'c'] :
            code = [0,1,0,0]
        elif letter in ['G', 'g'] :
            code = [0,0,1,0]
        elif letter in ['T', 't'] :
            code = [0,0,0,1]
        else :
            code = [1/4, 1/4, 1/4, 1/4]
        seq_code.append(code)
    return seq_code 

encoded_sequences = []
encoded_sequencesbw = []
with open(fasta_file_path) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequencebw = Seq(sequence).reverse_complement()

        encoded_sequence = encodeSeq(sequence)
        encoded_sequencebw = encodeSeq(sequencebw)

        # Apply z-score normalization
        normalized_sequence = scaler.fit_transform(encoded_sequence)
        normalized_sequencebw = scaler.fit_transform(encoded_sequencebw)
        encoded_sequences.append(normalized_sequence)
        encoded_sequencesbw.append(normalized_sequencebw)

    encoded_sequences = np.array(encoded_sequences)
    encoded_sequencesbw = np.array(encoded_sequencesbw)
    np.save(os.path.join(out_dir, contig_type + "_0k_codefw.npy"), encoded_sequences)
    np.save(os.path.join(out_dir, contig_type + "_0k_codebw.npy"), encoded_sequencesbw)
fasta_file.close()



