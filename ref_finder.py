#!/usr/bin/env python


# Finds the reference (and alternative) allele in the reference genome

# The MIT License (MIT)
#
# Copyright (c) 2015 Beaulieu-Saucier Pharmacogenomics Centre
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from __future__ import print_function

import os
import sys
import logging
import argparse
import unittest

import numpy as np
import pandas as pd
from pyfaidx import Fasta


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015, Beaulieu-Saucier Pharmacogenomics Centre. "
                 "All rights reserved.")
__credits__ = ["Louis-Philippe Lemieux Perreault", "Sylvie Provost"]
__license__ = "MIT"
__version__ = "0.2"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}


def main(args=None):
    """The main function."""
    # Creating the option parser
    desc = ("Finds the human reference allele "
            "(version {}).".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # The files that will need closing
    logging_fh = None
    reference = None

    # We run the script
    try:
        # Parsing the options
        args = parse_args(parser, args)
        check_args(args)

        # Adding the logging capability
        logging_fh = logging.FileHandler(args.log, mode="w")
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler(), logging_fh]
        )
        logging.info("Logging everything into '{}'".format(args.log))

        # Reading the BIM file
        variations = read_variations(args.i_filename)

        # Reading the reference
        reference = read_reference(args.reference)

        # Finds the reference and the alternative allele
        variations = find_ref_alt(variations, reference, args.prefix)

        # Saving the VCF file
        write_vcf_file(variations, args.prefix + ".vcf")

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    # Catching the ProgramError
    except ProgramError as e:
        parser.error(e.message)

    finally:
        if logging_fh is not None:
            logging_fh.close()
        if reference is not None:
            reference.close()


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


def find_ref_alt(variations, reference, prefix):
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
        logging.warning("There were variations with unknown or invalid "
                        "REF alleles (see '{}.invalid_ref' for a "
                        "list)".format(prefix))
        variations[unknown_ref].to_csv(
            prefix + ".invalid_ref",
            sep="\t",
            index=False,
            columns=["chrom", "pos", "name"],
        )

    # Finding the alternative alleles
    logging.info("Finding the alternative allele for all positions")
    variations["alt"] = [
        find_alt(*data)
        for data in variations.loc[:, ["ref", "a1", "a2"]].values
    ]

    # Are there any NaN values?
    unknown_alt = variations.alt.isnull()
    if (unknown_alt & (~unknown_ref)).any():
        logging.warning("There were variations incompatible with REF "
                        "allele (see '{}.invalid_ref_alt' for a "
                        "list)".format(prefix))
        variations[unknown_alt & (~unknown_ref)].to_csv(
            prefix + ".invalid_ref_alt",
            sep="\t",
            index=False,
            columns=["chrom", "pos", "name", "ref", "a1", "a2"],
        )

    return variations


def write_vcf_file(variations, filename):
    """Write a VCF file for CADD."""
    logging.info("Writing '{}'".format(filename))

    # Finding the complete data
    missing = variations.isnull().any(axis=1)
    with open(filename, "w") as i_file:
        print("##fileformat=VCFv4.1", file=i_file)
        print("#CHROM", "POS", "ID", "REF", "ALT", sep="\t", file=i_file)
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

    # The reference allele is undetermined
    if pd.isnull(ref):
        return np.nan

    # Being sure that everything is upper case
    a1 = a1.upper()
    a2 = a2.upper()

    # The 0 is a special case (since the variation is homozygous major in the
    # data set) (a2 will never be 0)
    if a1 == "0" and a2 in _complement:
        if (a2 != ref) and (_complement[a2] != ref):
            logging.debug("ref={}, a1={}, a2={}".format(ref, a1, a2))
            return a2
        return np.nan

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

    return True


def parse_args(parser, args=None):
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
    group.add_argument("-o", "--output-prefix", type=str, metavar="FILE",
                       dest="prefix", default="ref_alleles.vcf",
                       help="The name of the output file [%(default)s]")

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


class ProgramError(Exception):
    """An Exception raised in case of a problem."""
    def __init__(self, msg):
        """Construction of the ProgramError class."""
        self.message = str(msg)

    def __str__(self):
        return self.message


class Test(unittest.TestCase):
    def setUp(self):
        """Setting up the test cases."""
        # Tests with invalid data
        self.prefix1 = os.path.join("test", "ref_finder_test_with_invalid")
        args = [
            "--log", self.prefix1 + ".log",
            "--input", os.path.join("test", "test_invalid.bim"),
            "--reference", os.path.join("test", "test.reference"),
            "--output-prefix", self.prefix1,
        ]
        main(args=args)

        # Tests without invalid data
        self.prefix2 = os.path.join("test", "ref_finder_test")
        args = [
            "--log", self.prefix2 + ".log",
            "--input", os.path.join("test", "test.bim"),
            "--reference", os.path.join("test", "test.reference"),
            "--output-prefix", self.prefix2,
        ]
        main(args=args)

    def tearDown(self):
        """Deletes the output files."""
        import glob

        # Deleting files for the invalid data test
        for filename in glob.glob(self.prefix1 + ".*"):
            os.remove(filename)

        # Deleting files for the test
        for filename in glob.glob(self.prefix2 + ".*"):
            os.remove(filename)

    def test_output_file_present(self):
        """Check if all output files are present."""
        suffixes = [".log", ".vcf"]
        invalid_suffixes = [".invalid_ref", ".invalid_ref_alt"]

        # The first test should have all the output files
        checked_suffixes = suffixes + invalid_suffixes
        for filename in [self.prefix1 + suffix for suffix in checked_suffixes]:
            self.assertTrue(os.path.isfile(filename))

        # The second test should not have the invalid files
        for filename in [self.prefix2 + suffix for suffix in suffixes]:
            self.assertTrue(os.path.isfile(filename))
        for filename in [self.prefix2 + suffix for suffix in invalid_suffixes]:
            self.assertFalse(os.path.isfile(filename))

    def test_input_file_missing(self):
        """Tests the script raises an error if an input is missing."""
        # Missing BIM
        args = [
            "--log", self.prefix2 + ".log",
            "--input", os.path.join("test", "test_missing.bim"),
            "--reference", os.path.join("test", "test.reference"),
            "--output-prefix", self.prefix2,
        ]
        with self.assertRaises(SystemExit) as cm:
            main(args=args)
            self.assertEqual(cm.exception.code, 2)

        # Missing reference
        args = [
            "--log", self.prefix2 + ".log",
            "--input", os.path.join("test", "test.bim"),
            "--reference", os.path.join("test", "test_missing.reference"),
            "--output-prefix", self.prefix2,
        ]
        with self.assertRaises(SystemExit) as cm:
            main(args=args)
            self.assertEqual(cm.exception.code, 2)

    def test_invalid_reference(self):
        """Tests the script raises an error if there is a missing chrom."""
        args = [
            "--log", self.prefix2 + ".log",
            "--input", os.path.join("test", "test_invalid.bim"),
            "--reference", os.path.join("test", "test.invalid_reference"),
            "--output-prefix", self.prefix2,
        ]
        with self.assertRaises(SystemExit) as cm:
            main(args=args)
            self.assertEqual(cm.exception.code, 2)

    def test_missing_index(self):
        """Tests the script raises an error if the index is missing."""
        args = [
            "--log", self.prefix2 + ".log",
            "--input", os.path.join("test", "test.bim"),
            "--reference", os.path.join("test", "test.reference_no_index"),
            "--output-prefix", self.prefix2,
        ]
        with self.assertRaises(SystemExit) as cm:
            main(args=args)
            self.assertEqual(cm.exception.code, 2)

    def test_vcf(self):
        """Checks the VCF results."""
        # The valid file
        output = (
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\n"
            "1\t1\tmarker_1\tA\tG\n"
            "1\t2\tmarker_2\tC\tG\n"
            "1\t3\tmarker_3\tA\tC\n"
            "1\t4\tmarker_4\tC\tA\n"
            "1\t9\tmarker_9\tT\tG\n"
            "23\t1\tmarker_5\tA\tT\n"
            "24\t1\tmarker_6\tA\tT\n"
            "26\t1\tmarker_7\tA\tT\n"
        )
        content = None
        with open(self.prefix2 + ".vcf", "r") as i_file:
            content = i_file.read()
        self.assertEqual(output, content)

        # The invalid file
        output = (
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\n"
            "1\t1\tmarker_1\tA\tG\n"
            "1\t2\tmarker_2\tC\tG\n"
            "1\t3\tmarker_3\tA\tC\n"
            "1\t4\tmarker_4\tC\tA\n"
            "1\t5\tmarker_5\tT\tC\n"
            "1\t6\tmarker_6\tG\tT\n"
            "1\t9\tmarker_9\tT\tG\n"
            "1\t10\tmarker_10\tA\tG\n"
        )
        content = None
        with open(self.prefix1 + ".vcf", "r") as i_file:
            content = i_file.read()
        self.assertEqual(output, content)

    def test_invalid_ref(self):
        """Checks the 'invalid_ref' file."""
        output = (
            "chrom\tpos\tname\n"
            "1\t12\tmarker_12\n"
            "3\t1\tmarker_15\n"
        )
        content = None
        with open(self.prefix1 + ".invalid_ref", "r") as i_file:
            content = i_file.read()
        self.assertEqual(output, content)

    def test_invalid_ref_alt(self):
        """Checks the 'invalid_ref_alt' file."""
        output = (
            "chrom\tpos\tname\tref\ta1\ta2\n"
            "1\t7\tmarker_7\tT\t0\tA\n"
            "1\t8\tmarker_8\tA\t0\tA\n"
            "1\t11\tmarker_11\tC\t1\t2\n"
            "2\t1\tmarker_13\tA\tC\tG\n"
            "2\t2\tmarker_14\tC\tA\tT\n"
        )
        content = None
        with open(self.prefix1 + ".invalid_ref_alt", "r") as i_file:
            content = i_file.read()
        self.assertEqual(output, content)

if __name__ == "__main__":
    main()
