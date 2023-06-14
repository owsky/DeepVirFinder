import sys
import os
from typing import List, Tuple, Callable
subdirectory_path = os.path.join(os.path.dirname(__file__), 'normalizations')
sys.path.append(subdirectory_path)
from normalizations.normalizers import normalizers
from split_sequence import split_sequence
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
import optparse


def process_sequence(sequence: str, normalize_sequence: Callable[[str], str]) -> Tuple[str, str]:
    sequencebw = Seq(sequence).reverse_complement()
    normalized_sequence = normalize_sequence(sequence)
    normalized_sequencebw = normalize_sequence(sequencebw)
    return (normalized_sequence, normalized_sequencebw)


def normalize_contigs(sequences: List[str], normalize_sequence: Callable[[str], str]) -> Tuple[List[str], List[str]]:
    normalized_sequences = []
    normalized_sequencesbw = []
    for sequence in sequences:
        (ns, nsbw) = process_sequence(sequence, normalize_sequence)
        normalized_sequences.append(ns)
        normalized_sequencesbw.append(nsbw)
    return (normalized_sequences, normalized_sequencesbw)


def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", action = "store", type = "string", dest = "fasta_file_path", help = "fasta file input path")
    parser.add_option("-n", "--norm", action = "store", type = "string", dest = "normalization", help = "normalization method")
    parser.add_option("-c", "--contig", action = "store", type = "string", dest = "contig_type", help = "contig type")
    parser.add_option("-v", "--validation", action="store_true", default=False, dest="is_validation", help="validation flag")

    (options, _) = parser.parse_args()
    if (options.contig_type is None or options.fasta_file_path is None or options.normalization is None) :
        sys.stderr.write(sys.argv[0] + ": ERROR: missing required command-line argument")
        parser.print_help()
        sys.exit(1)

    fasta_file_path = options.fasta_file_path
    normalization = options.normalization
    contig_type = options.contig_type
    is_validation = options.is_validation
    contig_lengths = [150, 300, 500, 1000]

    out_dir = os.path.join(os.path.dirname(fasta_file_path), "encode_" + normalization)
    os.makedirs(out_dir, exist_ok=True)

    if normalization in normalizers:
        normalize_sequence = normalizers[normalization]
        print(f"Normalizing with {normalization}")
    else:
        sys.stderr.write(f"Unsupported normalization method. Use any of the following: {normalizers.keys()}")
        sys.exit(1)

    encoded_sequences = []
    encoded_sequencesbw = []

    if is_validation:
        for l in contig_lengths:
            encoded_sequences = []
            encoded_sequencesbw = []
            with open(fasta_file_path) as fasta_file:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequence = str(record.seq)
                    contigs = split_sequence(sequence, l)
                    (nss, nssbw) = normalize_contigs(contigs, normalize_sequence)
                    for ns in nss:
                        encoded_sequences.append(ns)
                    for nsbw in nssbw:
                        encoded_sequencesbw.append(nsbw)

                encoded_sequences = np.array(encoded_sequences)
                encoded_sequencesbw = np.array(encoded_sequencesbw)
                if len(encoded_sequences.shape) < 3:
                    encoded_sequences = np.reshape(encoded_sequences, (*encoded_sequences.shape, 1))
                    encoded_sequencesbw = np.reshape(encoded_sequencesbw, (*encoded_sequencesbw.shape, 1))
                    
                np.save(os.path.join(out_dir, contig_type + f"_{l/1000}k_codefw.npy"), encoded_sequences)
                np.save(os.path.join(out_dir, contig_type + f"_{l/1000}k_codebw.npy"), encoded_sequencesbw)
            fasta_file.close()
    else:
        with open(fasta_file_path) as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence = str(record.seq)
                (ns, nsbw) = process_sequence(sequence, normalize_sequence)
                encoded_sequences.append(ns)
                encoded_sequencesbw.append(nsbw)

            encoded_sequences = np.array(encoded_sequences)
            encoded_sequencesbw = np.array(encoded_sequencesbw)
            if len(encoded_sequences.shape) < 3:
                encoded_sequences = np.reshape(encoded_sequences, (*encoded_sequences.shape, 1))
                encoded_sequencesbw = np.reshape(encoded_sequencesbw, (*encoded_sequencesbw.shape, 1))
                
            np.save(os.path.join(out_dir, contig_type + "_0k_codefw.npy"), encoded_sequences)
            np.save(os.path.join(out_dir, contig_type + "_0k_codebw.npy"), encoded_sequencesbw)
        fasta_file.close()

if __name__ == "__main__":
    main()