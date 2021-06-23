#!/usr/bin/env python

import os
import sys
import errno
import argparse

def parse_args(args=None):
    Description = "Generate fasta file from library reference file"
    Epilog = "Example usage: python get_library_fasta <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Library reference file.")
    parser.add_argument("FILE_OUT", help="Output fasta file.")
    return parser.parse_args(args)

def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check library reference file -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please library reference file -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

# TODO nf-core: Update the get_library_fasta function
def get_library_fasta(file_in, file_out):
    """
    This function checks that the library fasta follows the following structure:

    id,sequence,gene
    s_10007,TGTTCACAGTATAGTTTGCC,CCNA1
    s_10008,TTCTCCCTAATTGCTTGCTG,CCNA1
    s_10027,ACATGTTGCTTCCCCTTGCA,CCNC
    """

    fasta_run_dict = {}
    print(file_in)
    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 3
        # TODO nf-core: Update the column names for the input samplesheet
        HEADER = ["id", "sequence", "gene"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check reference file header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            id, sequence, gene = lspl[: len(HEADER)]
            if id:
                if id.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            ## Check sequence consists only of ACTG bases
            ## Check reference library file extension
            for seq in [sequence]:
                if seq.find(" ") != -1:
                    print_error("sequence contains spaces!", "Line", line)

                for base in seq:
                    if base not in "ACTG":
                        print_error("library sequence contains wrong bases!", "Line", line)

            ## Create Dictionary with sgRNAs, genes and sequence
            # new_id = "_".join([id, gene])

            ## Create sample mapping dictionary = {sample: {replicate : [ single_end, fastq_1, fastq_2 ]}}
            if id not in fasta_run_dict:
              fasta_run_dict[id] = sequence
            else:
                print_error("Reference Library File contains duplicate rows!", "Line", line)

    ## Write fasta reference file
    if len(fasta_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            for id in sorted(fasta_run_dict.keys()):
                ## Check that replicate ids are in format 1..<NUM_REPS>
                ## Write to file
                seq = fasta_run_dict[id]
                fout.write(">" + id + "\n")
                fout.write(seq + "\n")

def main(args=None):
    args = parse_args(args)
    get_library_fasta(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
