""" This script converts fasta files to genbank files
    It calls two functions:
    fasta_to_genbank to create a genbank file from a fasta file
    genbank_to_dict to create a dictionary from the genbank file
    This could be run for example as:
    python3 fastaTogenbank.py -i example_file.fasta -o OUTPUTS/example_file.gb
    python3 fastaTogenbank.py -i example_file.fasta -o OUTPUTS/example_file.gb --print_stats
    python3 fastaTogenbank.py -i example_file.fasta -o OUTPUTS/example_file.gb -p

"""

import argparse
from Bio import SeqIO
import statistics
import re
import os

statistics_help_message = "If true prints file statistics"


def print_fasta_statistics(fasta_file): 

    """ 
    This function calculates and prints the following statistics to the command line:
        Mean read length, Median read length, Maximum read length and Minimum read length.

    Parameters:
        filename (str): The name of the file to check

    Returns:
        statisticsdictionary (dictionary): contains the calculated statistics
    """

    # Load all records from your FASTA file into a list
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Calculate total number of reads
    total_reads = len(records)
    print(f"Total reads: {total_reads}")

    # Calculate sizes of all reads
    sizes = [len(record.seq) for record in records]

    # Calculate and print various statistics
    print(f"Mean read length: {statistics.mean(sizes)}")
    print(f"Median read length: {statistics.median(sizes)}")
    print(f"Max read length: {max(sizes)}")
    print(f"Min read length: {min(sizes)}")
    statisticsdictionary = {'total_reads': total_reads, 'mean_read_length': statistics.mean(sizes),'median_read_length': statistics.median(sizes), 'max_read_length': max(sizes), 'min_read_length': min(sizes)}
    return statisticsdictionary



# Function to check the file extension
def is_fasta(filename):

    """ 
    This function checks if the file is a FASTA file based on its extension. 
    It returns true if the extension of the file filename is either .fasta or .fa.

    Parameters:
        filename (str): The name of the file to check

    Returns:
        bool: True if the file has a .fasta or a .fa extension, False otherwise
    """

    return filename.endswith('.fasta') or filename.endswith('.fa')

# Function to check the FASTA header
def has_valid_header(filename):
    
    """
    Check if a FASTA file has a valid header.

    The function reads the first line of the file and checks if it starts with '>',
    which is the standard header prefix in FASTA format.

    Parameters:
    filename (str): The name of the FASTA file to check.

    Returns:
    bool: True if the file has a valid FASTA header, False otherwise.
    """

    with open(filename, 'r') as file:
        first_line = file.readline()
        # A typical FASTA header starts with '>'
        return bool(re.match(r'^>', first_line))



def parse_file(fasta_file):    
    """
    Parse a FASTA file and annotate each sequence with its molecule type (DNA).

    This function reads a FASTA file, parses each sequence into a SeqRecord object,
    and adds an annotation indicating the molecule type as 'DNA'. A list of
    annotated SeqRecord objects is then returned.

    Parameters:
    fasta_file (str): The path to the FASTA file to be parsed.

    Returns:
    list: A list of SeqRecord objects with 'molecule_type' annotations.
    """

    # Parse the fasta file into a list
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    # Add the molecule_type annotation necessary for SeqIO.write
    for sequence in sequences:
        sequence.annotations["molecule_type"] = "DNA"
    return(sequences)

def write_genbank(parsed_list, genbank_file="genbankOutput"):   
    
    """
    This function takes in the parsed list derived from a fasta file and the output path of a genbank file.
    It will create a genbank file at that path.

    Parameters:
        parsed_list (list): parsed list object
        genbank_file (str): output path 
    Returns:
        genbank file generated from the parsed list
    """
    output = None
    ## Write the genbank file
    try:
        output = SeqIO.write(parsed_list, genbank_file, "genbank")
    except IOError as error:
        print(f"An I/O error occurred: {error}")
        raise IOError("An I/O error occurred:  ")
        return error
    except Exception as error:
        print(f"An unexpected error occurred: {error}")
        raise Exception
        return error
    else:
        print("Conversion process completed Successfully.")
        return output


def main():

    # This creates a parser object
    # the description may be accessed by adding the argument -h
    parser = argparse.ArgumentParser(description="A script to convert fasta to genbank and create a dictionary")

    # Add arguments to the parser, I set default names in the functions
    parser.add_argument("-i", "--input", help="The input fasta file", type=str, required=True)
    parser.add_argument("-o", "--output", help="The output genbank file", type=str, required=False)
    parser.add_argument("-p", "--print_stats", help=statistics_help_message, action='store_true', required=False)

    # Parse the arguments
    args = parser.parse_args()

    if is_fasta(args.input) and has_valid_header(args.input):
        print(f"{args.input} is a valid FASTA file.")
        if args.print_stats is True:
            print(print_fasta_statistics(args.input))
        write_genbank(parse_file(args.input), args.output)
    else:
        print(f"{args.input} is not a valid FASTA file.")
    pass


if __name__ == "__main__":
    main()
