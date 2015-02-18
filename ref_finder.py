#!/usr/bin/env python
"""The main script that do everything."""


import os
import sys
import logging
import argparse

import numpy as np
import pandas as pd
from pyfaidx import Fasta


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Sylvie Provost"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}


def main():
    """The main function."""
    # Creating the option parser
    desc = ("Finds the human reference allele.".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # We run the script
    try:
        # Parsing the options
        args = parse_args(parser)
        check_args(args)

        # Adding the logging capability
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler(),
                      logging.FileHandler(args.log, mode="w")]
        )
        logging.info("Logging everything into '{}'".format(args.log))

        # Reading the BIM file
        variations = read_variations(args.i_filename)

        # Reading the reference
        reference = read_reference(args.reference)

        # Finds the reference and the alternative allele
        variations = find_ref_alt(variations, reference)

        # Saving the VCF file
        write_vcf_file(variations, args.output)

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    # Catching the ProgramError
    except ProgramError as e:
        parser.error(e.message)


def read_variations(filename):
    """Reads variations from a BIM file."""
    logging.info("Reading variations '{}'".format(filename))
    variations = pd.read_csv(filename, sep="\t", names=["chrom", "name", "cm",
                                                        "pos", "a1", "a2"])
    logging.info("  - read {:,d} variations".format(len(variations)))

    return variations


def read_reference(filename):
    """Reads the reference using pyfaidx."""
    logging.info("Reading the reference's index '{}'".format(filename))
    ref = Fasta(filename)
    logging.info("  - read {:,d} sequences".format(len(ref.keys())))

    return ref


def find_ref_alt(variations, reference):
    """Finds the reference and the alternative alleles."""
    # Finding the reference alleles
    logging.info("Finding the reference allele for all positions")
    variations["ref"] = [
        find_ref(reference, *data)
        for data in variations.loc[:, ["chrom", "pos"]].values
    ]

    # Are there any unknown reference alleles?
    unknown_ref = variations.ref.isnull()
    if unknown_ref.any():
        logging.warning("The following variations had unknown or invalid "
                        "REF alleles")
        for name in variations.name[unknown_ref]:
            logging.warning("  - {}".format(name))

    # Finding the alternative alleles
    logging.info("Finding the alternative allele for all positions")
    variations["alt"] = [
        find_alt(*data)
        for data in variations.loc[:, ["ref", "a1", "a2"]].values
    ]

    # Are there any NaN values?
    unknown_alt = variations.alt.isnull()
    if (unknown_alt & (~unknown_ref)).any():
        logging.warning("The following variations were incompatible with REF "
                        "allele")
        for index in variations.name[unknown_alt & (~unknown_ref)].index:
            name, a1, a2, ref = variations.loc[index,
                                               ["name", "a1", "a2", "ref"]]
            logging.warning(
                "  - {}: REF={}, GENO={}/{}".format(name, ref, a1, a2)
            )

    return variations


def write_vcf_file(variations, filename):
    """Write a VCF file for CADD."""
    logging.info("Writing '{}'".format(filename))

    # Finding the complete data
    missing = variations.isnull().any(axis=1)
    with open(filename, "w") as i_file:
        print("##fileformat=VCFv4.1", file=i_file)
        print("#CHROM", "POS", "ID", "REF ALT", sep="\t", file=i_file)
        variations[~missing].to_csv(
            i_file,
            sep="\t",
            header=False,
            index=False,
            columns=["chrom", "pos", "name", "ref", "alt"],
        )
    logging.info("  - wrote {:,d} variations".format((~missing).sum()))


def find_ref(ref, chrom, pos):
    """Finds the reference allele at this position."""
    if encode_chrom(chrom) not in ref:
        raise ProgramError("{}: invalid chromosome".format(chrom))
    ref_allele = ref[encode_chrom(chrom)][int(pos)-1].seq.upper()

    # Is it a valid allele?
    if ref_allele in _complement:
        return ref_allele

    return np.nan


def find_alt(ref, a1, a2):
    """Finds the alternative allele (a2 is the most common)."""
    logging.debug("ref={}, a1={}, a2={}".format(ref, a1, a2))

    # The reference allele is undetermined
    if pd.isnull(ref):
        return np.nan

    # Being sure that everything is upper case
    a1 = a1.upper()
    a2 = a2.upper()

    # Checks for undetermined alleles
    if (a1 not in _complement) or (a2 not in _complement):
        return np.nan

    # a2 is the reference allele
    if ref == a2:
        return a1

    # a2 is the complement of the reference allele
    if ref == _complement[a2]:
        return _complement[a1]

    # a1 is the reference allele
    if ref == a1:
        return a2

    # a1 is the complement of the reference allele
    if ref == _complement[a1]:
        return _complement[a2]

    return np.nan


def encode_chrom(chrom):
    """Encodes the chromosome."""
    if str(chrom) == "23":
        return "X"

    if str(chrom) == "24":
        return "Y"

    if str(chrom) == "26":
        return "M"

    return str(chrom)


def check_args(args):
    """Checks the arguments and options."""
    # Checking the input file
    if not os.path.isfile(args.i_filename):
        raise ProgramError("{}: no such file".format(args.i_filename))

    # Checking the fasta file (and index)
    if not os.path.isfile(args.reference):
        raise ProgramError("{}: no such file".format(args.reference))
    if not os.path.isfile(args.reference + ".fai"):
        raise ProgramError("{}: no such file".format(args.reference + ".fai"))

    # Checking if the output file is a VCF
    if not args.output.endswith(".vcf"):
        args.output += ".vcf"

    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__))
    parser.add_argument("--debug", action="store_true",
                        help="Set the logging to debug")
    parser.add_argument("--log", type=str, metavar="LOGFILE",
                        default="ref_finder.log",
                        help="The log file [%(default)s]")

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument("-i", "--input", type=str, metavar="FILE",
                       required=True, dest="i_filename",
                       help="The input file (BIM format).")
    group.add_argument("-r", "--reference", type=str, metavar="FASTA",
                       required=True, help=("The human reference genome in "
                                            "FASTA format."))

    # The result file
    group = parser.add_argument_group("Result File")
    group.add_argument("-o", "--output", type=str, metavar="FILE",
                       default="ref_alleles.vcf",
                       help="The name of the output file [%(default)s]")

    return parser.parse_args()


class ProgramError(Exception):
    """An Exception raised in case of a problem."""
    def __init__(self, msg):
        """Construction of the ProgramError class."""
        self.message = str(msg)

    def __str__(self):
        return self.message


if __name__ == "__main__":
    main()
