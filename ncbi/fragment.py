from Bio import SeqIO
import random
import os
from threading import Thread


def get_subsequence(l: int, directory: str, output_files: list):
    output_dir = os.path.join(directory, f'L_{l}')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'output.fa')
    with open(output_file, 'w') as output_handle:
        for filename in os.listdir(directory):
            if filename.endswith('.fa'):
                with open(os.path.join(directory, filename), 'r') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        sequence = str(record.seq)
                        if l > len(sequence):
                            # print(f"Sequence in file {directory + '/' + filename} is shorter than l: {l}")
                            continue
                        start = random.randint(0, len(sequence)-l)
                        end = start + l
                        subsequence = sequence[start:end]
                        output_handle.write(f'>{record.id}#{l}#{start}#{end}\n{subsequence}\n')
    output_files.append(output_file)



def get_random_subsequences(directories: list, L: list) -> list:
    threads = []
    output_files = []
    for directory in directories:
        for l in L:
            t = Thread(target=get_subsequence, args=(l, directory, output_files))
            t.start()
            threads.append(t)
    for t in threads:
        t.join()
    return output_files


def main():
    directories = ["./downloads/virus-after-2015-refseqs/", "./downloads/virus-before-2014-refseqs/", "./downloads/virus-between-2014.1-2015-refseqs/", "./downloads/prokaryote-2014.1-2015-refseqs", "./downloads/prokaryote-after-2015-refseqs", "./downloads/prokaryote-before-2014-refseqs"]
    get_random_subsequences(directories, [3000])


if __name__ == "__main__":
    main()