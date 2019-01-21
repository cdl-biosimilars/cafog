#!/usr/bin/env python3

from argparse import ArgumentParser
import logging
import os
import sys

from correction import GlycationGraph, read_clean_datasets, read_library


def setup_parser() -> ArgumentParser:
    """
    Set up options for the argument parser

    :return: a parser for handling command line arguments
    :rtype: ArgumentParser
    """
    parser = ArgumentParser(
        description="Correct glycation influence on glycoform abundances. "
                    "Results are written in CSV format to STDOUT.")

    parser.add_argument("-f", "--glycoforms",
                        action="store",
                        help="CSV file containing glycoform abundances",
                        required=True)
    parser.add_argument("-g", "--glycation",
                        action="store",
                        help="CSV file containing glycation abundances",
                        required=True)
    parser.add_argument("-l", "--glycan-library",
                        action="store",
                        help="CSV file containing a glycan library")
    parser.add_argument("-o", "--graph-output-format",
                        action="store",
                        help="graph output format, either 'dot' or 'gexf' ",
                        metavar="FORMAT",
                        choices=["dot", "gexf"])
    parser.add_argument("-v", "--version",
                        action="version",
                        help="print the version number",
                        version="%(prog)s 1.0")
    return parser


def _main() -> None:
    """
    Parse command line arguments, correct abundances and write output.

    :return: nothing
    :rtype: None
    """

    logging.basicConfig(format="%(levelname)s: %(message)s",
                        level=logging.INFO)
    args = setup_parser().parse_args()

    # read input files
    dataset_name = os.path.splitext(args.glycoforms)[0]
    try:
        glycoforms = read_clean_datasets(args.glycoforms)
        glycation = read_clean_datasets(args.glycation)
        if args.glycan_library is None:
            glycan_library = None
        else:
            glycan_library = read_library(args.glycan_library)
    except (OSError, ValueError) as e:
        logging.error(e)
        sys.exit(-1)

    # assemble the glycation graph, correct abundances and store
    logging.info("Correcting dataset '{}' …".format(args.glycoforms))
    try:
        G = GlycationGraph(glycan_library, glycoforms, glycation)
        G.correct_abundances()
        G.to_dataframe().to_csv(sys.stdout, index=False)
        if args.graph_output_format == "dot":
            G.to_dot("{}_corr.gv".format(dataset_name))
        elif args.graph_output_format == "gexf":
            G.to_gexf("{}_corr.gexf".format(dataset_name))
    except ValueError as e:
        logging.error(e)
        sys.exit(1)

    logging.info("… done!")


if __name__ == "__main__":
    _main()
