#!/usr/bin/env python

import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Generates library input files for Mageck and PinAPL-py"
    Epilog = "Example usage: python get_library_fasta <FILE_IN> <MAGECK_OUT> <PINAPL_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Library reference file.")
    # parser.add_argument("MAGECK_OUT", help="Output file for Mageck")
    parser.add_argument("PINAPL_OUT", help="Output file in the correct format for Pinapl-py")
    parser.add_argument("MAGECK_OUT", help="Output file in the correct format for Mageck")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    try:
        read_file = pd.read_excel(args.FILE_IN, header=0)
    except Exception:
        read_file = pd.read_csv(args.FILE_IN, header=0)

    # reorder and write file

    out = read_file[["id", "sequence", "gene"]]
    out.to_csv(args.PINAPL_OUT, header=True, index=False, sep=",")
    out.to_csv(args.MAGECK_OUT, header=True, index=False, sep=",")

main()
exit()
# sgRNA library conversions - reorder columns for library file of PinAPL-Py and MaGeCK #
df = pd.read_csv ('mouse_geckov2_library_a_2_mageck.csv', header=None) # change file to file of your choice
# check which column contains the gene sequence (length = 20, only letters)
# save the original file under a new name to reflect the tool it belongs to (i.e. mageck or pinapl-py)
# re-order the columns and save under the name of the tool it belongs to
if len(df.iloc[0][1]) == 20 and df.iloc[0][1].isalpha() == True: # if true, the file is in t he correct order for Mageck
    print('Mageck', df.iloc[0]) # not necessary, just to check if the order is indeed correct
    df.to_csv('library_mageck.csv', header=False, index=False)
    df.columns=['ID', 'seq', 'gene']
    df_pinapl = df[['gene', 'ID', 'seq']]
    df_pinapl.to_csv('library_pinapl.csv', header=False, index=False)
elif len(df.iloc[0][2]) == 20 and df.iloc[0][1].isalpha() == True: # if true, the file is in the correct order for Pinapl-py
    print('Pinapl_py', df.iloc[0])  #not necessary, just to check if the order is indeed correct
    df.to_csv('library_pinapl.csv', header=False, index=False)
    df.columns=['gene', 'seq', 'ID']
    df_mageck=df[['ID', 'seq', 'gene']]
    df_mageck.to_csv('library_mageck.csv', header=False, index=False)
else:
    print('Column order does not fit mageck or pinapl-py', df.iloc[0])


exit()
