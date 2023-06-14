from typing import List


def split_sequence(sequence: str, l: int) -> List[str]:
    contigs = []
    sequence_length = len(sequence)
    num_contigs = sequence_length // l
    
    for i in range(num_contigs):
        contig = sequence[i*l:(i+1)*l]
        contigs.append(contig)
    
    return contigs
