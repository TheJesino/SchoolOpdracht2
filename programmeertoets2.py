import tkinter
from tkinter import *
from tkinter import filedialog
import pickle
import matplotlib.pyplot as plt
import re

class MyClass():
    def __init__(self, param):
        self.param = param

def parse_fasta(fasta_file):
    """Parse a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(fasta_file, 'r') as f:
        print("Opening FASTA file..")
        counter = 0
        lastentry = 0
        for line in f:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                sequences[seq_name] = ""
            else:
                sequences[seq_name] += line.strip()
                counter+=1
                if int((counter/519756)*100) is not lastentry:
                    lastentry = int((counter / 519756) * 100)
                    print(int((counter/519756)*100),"%")

    print("Werkt nog")
    return sequences


def parse_gff3(gff3_file):
    """Parse a GFF3 file and return a dictionary of gene locations."""
    genes = {}
    print("Opening GFF3 file...")
    with open(gff3_file, 'r') as f:
        for gene_name, gene_data in genes.items():
            if line.startswith("#"):
                continue
            elif not line.startswith("NC"):
                continue
            fields = line.strip().split("\t")
            print(fields)
            if fields[2] == "gene":
                gene_name = re.search("gene=(.+?);", fields[8]).group(1)
                type = fields[2]
                start, end = int(fields[3]), int(fields[4])
                strand = fields[6]
                seq_id = fields[0]
                genes[gene_name] = {"seq_id": seq_id, "start": start, "end": end, "strand": strand}
    print("GFF3 File opened!")
    return genes

def extract_gene_sequences(fasta_file, gff3_file, output_file):
    """Extract the gene sequences from a FASTA file and write them to a new FASTA file."""
    sequences = parse_fasta(fasta_file)
    genes = parse_gff3(gff3_file)
    counter = 0
    print("Begint met schrijven..")
    with open(output_file, 'w') as f:
        for gene_name, gene_data in genes.items():
            seq_id = gene_data["seq_id"]
            start = gene_data["start"]
            end = gene_data["end"]
            strand = gene_data["strand"]
            if strand == "+":
                if seq_id is not "NC_05238.1":
                    continue
                gene_seq = sequences[seq_id][start-1:end]
            elif strand == "-":
                if seq_id is not "NC_05238.1":
                    continue
                gene_seq = sequences[seq_id][start-1:end].reverse_complement()
            f.write(">" + gene_name + "\n")
            f.write(gene_seq + "\n")
            counter += 1
            if int((counter / 1057) * 100) is not lastentry:
                lastentry = int((counter / 1057) * 100)
                print(int((counter / 1057) * 100), "%")
    print("Is klaar met schrijven!")



def pickle_time1(obj):
    try:
        with open("datafasta.pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error:", ex)

def pickle_time2(obj):
    try:
        with open("datagff.pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error:", ex)

def load_object(filename):
    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error:", ex)

def fileGUI():
    def browseFasta():
        filenameFASTA = filedialog.askopenfilename(initialdir="/",
                                                   title="Selecteer een bestand",
                                                   filetypes=(("FASTA bestanden",
                                                               ".fasta"),
                                                              ("Text bestanden",
                                                                ".txt")))
        FASTA_label.configure(text="Geopend bestand: " + filenameFASTA)
        obj = MyClass(str(filenameFASTA))
        pickle_time1(obj)



    def browseGFF():
        filenameGFF = filedialog.askopenfilename(initialdir="/",
                                                 title="Selecteer een bestand",
                                                 filetypes=(("GFF3 bestanden",
                                                             ".GFF3"),
                                                            ("Text bestanden",
                                                             ".txt")))
        GFF3_label.configure(text="Geopend bestand: " + filenameGFF)
        obj = MyClass(str(filenameGFF))
        pickle_time2(obj)


    win = tkinter.Tk()
    win.title("File Explorer")
    win.geometry('900x900')
    FASTA_label = Label(win,
                           text = "Nog geen FASTA bestand geopend",
                           width = 100, height = 4,
                           fg = "blue")
    FASTA_button = Button(win,
                          text = "Browse FASTA bestanden",
                          command = browseFasta)

    GFF3_label = Label(win,
                       text = "Nog geen GFF3 bestand geopend",
                       width = 100, height = 4,
                       fg = "blue")
    GFF3_button = Button(win,
                         text = "Browse GFF3 bestanden",
                         command = browseGFF)

    exit_button = Button(win,
                         text = "Afsluiten",
                         command = exit)
    FASTA_label.grid(column = 1, row = 1)
    FASTA_button.grid(column = 1, row = 2)
    GFF3_label.grid(column = 1, row = 5)
    GFF3_button.grid(column = 1, row = 8)
    exit_button.grid(column = 1, row = 10)
    win.mainloop()

def plot():
    # Open the GFF3 file and iterate through the rows
    genes = {}
    with open("gallus_gallus.gff3", 'w') as f:
        for gene_name, gene_data in genes.items():
            if line.startswith("#"):
                continue
            elif not line.startswith("NC_052538.1"):
                continue
            fields = line.strip().split("\t")
            print(fields)
            if fields[2] == "gene":
                gene_name = re.search("gene=(.+?);", fields[8]).group(1)
                type = fields[2]
                start, end = int(fields[3]), int(fields[4])
                strand = fields[6]
                seq_id = fields[0]
                genes[gene_name] = {"seq_id": seq_id, "start": start, "end": end, "strand": strand}

    with open("gallus_gallus.gff3", "r") as gff3_file:
        for row in gff3_file:
            # Skip comment rows
            if row.startswith("#"):
                continue

            # Split the row into columns
            columns = row.strip().split("\t")

            # Increment the count for the current feature type
            if columns[2] == "gene":
                feature_type = columns[2]
                feature_counts[feature_type] += 1

            # Extract the feature types and counts as separate lists
        feature_types = list(feature_counts.keys())
        feature_counts = list(feature_counts.values())

    # Create a bar chart of the feature occurrences
    plt.bar(feature_types, feature_counts)
    plt.xlabel("Feature Type")
    plt.ylabel("Occurrences")
    plt.show()
def main():
    extract_gene_sequences("gallus_gallus.FASTA", "gallus_gallus.GFF3", "new.fasta")
    fileGUI()
    extract_gene_sequences("gallus_gallus.FASTA", "gallus_gallus.GFF3", "new.fasta")
    obj = load_object("datagff.pickle")
    obj2 = load_object("datafasta.pickle")
    print(list(obj))
    print(list(obj2))
    plot()



main()