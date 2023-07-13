from Bio import Entrez
from csv import reader
from os import listdir, system, remove
from os.path import join, splitext, isdir, isfile, getsize
from subprocess import check_output
from concurrent.futures import ThreadPoolExecutor
from progress_bar import percent_complete
from time import sleep

Entrez.email = "EMAIL"
Entrez.api_key = "API_KEY"

def line_count(filename):
    return int(check_output(['wc', '-l', filename]).split()[0])

def download_data_thread(accession: str, save_path: str):
    file_out = save_path + "/" + accession + ".fa"
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession)
        record = Entrez.read(handle)
        handle.close()
        sleep(0.4)

        id_list = record["IdList"]

        fasta_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta")
        fasta_data = fasta_handle.read()
        fasta_handle.close()
        sleep(0.4)
            
        with open(file_out, "w") as f:
            f.write(fasta_data)
    except Exception as e:
        with open(file_out, "w") as f:
            print(f"\n####################################\n{e}\nAccession number: {accession}\n##############################")
            f.write("")

i = 0
def done(num_genomes: int, title: str):
    global i
    percent_complete(step=i, total_steps=num_genomes, title=title)
    i = i + 1
    
def download_data(file_path: str, save_path: str, file_name: str, num_genomes: int):
    with open(file_path) as csvfile:
        file = reader(csvfile)
        next(file)  # skip the headings
        with ThreadPoolExecutor(max_workers=10) as executor:
            for row in file:
                accession = row[0]
                file_out = save_path + "/" + accession + ".fa"
                if not isfile(file_out) or getsize(file_out) == 0:
                    future = executor.submit(download_data_thread, accession, save_path)
                    future.add_done_callback(lambda _: done(num_genomes, file_name))
                else:
                    # print(f"Already got {accession} genome")
                    done(num_genomes, file_name)

def main():
    global i
    base_path = "./SupplementaryTable1/"
    csv_files = [f for f in listdir(base_path) if join(base_path, f)]
    
    for csv_file in csv_files:
        i = 0
        file_name = splitext(csv_file)[0]
        file_path = base_path + csv_file
        save_path = "./downloads/" + file_name + "-refseqs"

        num_genomes = line_count(file_path) - 1 # account for headings

        print(f"\nNow downloading {num_genomes} genomes from file {csv_file}")

        if not isdir(save_path):
            print(f"Creating folder {save_path}")
            system("mkdir " + save_path)

        percent_complete(step=0, total_steps=num_genomes, title=file_name)
        download_data(file_path, save_path, file_name, num_genomes)
        percent_complete(step=num_genomes, total_steps=num_genomes, title=file_name)
        print("\n")

        print("Checking for failures...")
        results = [f for f in listdir(save_path) if join(save_path, f)]
        failures = []
        for result in results:
            result_path = join(save_path, result)
            if getsize(result_path) == 0:
                failures.append(splitext(result)[0] + "\n")
                remove(result_path)
        if len(failures) > 0:
            fail_path = save_path + "/!failures.txt"
            print(f"\nFailed to retrieve {len(failures)} genomes out of {num_genomes}. See {fail_path} for the list accession numbers\n")
            with open(fail_path, "w") as f:
                f.writelines(failures)



if __name__ == "__main__":
    main()