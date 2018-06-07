#!/usr/bin/env python3

import argparse
import logging
import re
import sys

from multiset import FrozenMultiset
import networkx as nx
import numpy as np
import pandas as pd

from glycan import PTMComposition
from glycoprotein import Glycoprotein


def calc_glycation_graph(glycan_library, glycoforms, glycation):
    """
    Assemble the glycoform graph from peptide mapping
    and glycation frequency data.

    :param str glycan_library: CSV file containing a glycan library
    :param pd.Series glycoforms: list of glycoforms with abundances
    :param pd.Series glycation: list of glycation counts with abundance
    :return: the glycation DAG
    :rtype: nx.DiGraph
    """

    # regex for extracting the first glycoform from a string like
    # "A2G0F/A2G1F or A2G1F/A2G0F"
    re_first_glycoform = re.compile("([^\s]*)")

    # series with a monosaccharide set as index
    # and experimental abundances as values
    exp_abundances = glycoforms.rename("abundance").reset_index()
    exp_abundances["sugar_set"] = exp_abundances.glycoform.apply(
        lambda v: FrozenMultiset(v.split("/")))
    exp_abundances = exp_abundances.set_index("sugar_set").abundance

    # dict mapping PTM differences to abundances
    delta_ptm = {PTMComposition({"Hex": count}): abundance
                 for count, abundance in glycation.iteritems()
                 if count > 0}

    # compare monosaccharide set of glycan library and glycoforms
    # add glycans that only appear in the list of glycoforms to the library
    # but this only works if they have a valid name
    gp = Glycoprotein(sites=2, library=glycan_library)
    library_glycans = set([n.name for n in gp.glycan_library])
    glycoform_glycans = set()
    for v in exp_abundances.index.values:
        glycoform_glycans |= set(v)

    glycans_only_in_library = library_glycans - glycoform_glycans
    if glycans_only_in_library:
        logging.warning(
            "The following glycans only appear in the glycan library, "
            + "but not in the list of glycoforms: "
            + str(glycans_only_in_library)
            + ".")

    glycans_only_in_glycoforms = glycoform_glycans - library_glycans
    if glycans_only_in_glycoforms:
        logging.warning(
            "The following glycans only appear in the list of glycoforms, "
            + "but not in the glycan library: "
            + str(glycans_only_in_glycoforms)
            + ". They will be added to the library.")
        for g in glycans_only_in_glycoforms:
            try:
                gp.add_glycan(g)
            except ValueError as e:
                logging.error(e)
                sys.exit(1)

    G = nx.DiGraph()
    for glycoform in gp.unique_glycoforms():
        # get the experimental abundance of a glycoform
        # use a default value of 0 if unavailable
        abundance = 0.0
        for name in glycoform.name.split(" or "):
            try:
                abundance = exp_abundances[FrozenMultiset(name.split("/"))]
                break
            except KeyError:
                continue
        glycoform.abundance = abundance

        # add the current glycoform as a node to the graph;
        # generate an edge to previous nodes if the difference is described in
        # the dict of PTM differences
        G.add_node(
            glycoform,
            abundance=float(abundance),
            label=re_first_glycoform.match(glycoform.name).group())
        for n in G:
            d = glycoform - n
            try:
                c = delta_ptm[d]
                source = n
                sink = glycoform
            except KeyError:
                try:
                    d = -d
                    c = delta_ptm[d]
                    source = glycoform
                    sink = n
                except KeyError:
                    continue
            G.add_edge(source, sink, label=d.composition_str(), c=float(c))
    return G


def correct_abundances(G):
    """
    Correct abundances in the glycoform graph.

    :param nx.DiGraph G: a glycoform graph
    :return: nothing, G is modified in place
    """

    total_abundance = 0.0

    # first pass: abundance correction,
    # calculated for each node from source to sink
    for n in nx.topological_sort(G):
        in_abundance = 0.0
        for pred in G.predecessors(n):
            in_abundance += G.nodes[pred]["corr_abundance"] * G[pred][n]["c"]
        out_c = 0.0
        for succ in G.successors(n):
            out_c += G[n][succ]["c"]
        corr_abundance = (n.abundance - in_abundance) / (1 - out_c)
        total_abundance += corr_abundance
        G.nodes[n]["corr_abundance"] = float(corr_abundance)

    # second pass: calculate fractional abundance
    for n in G:
        G.nodes[n]["corr_abundance"] = float(G.nodes[n]["corr_abundance"]
                                             / total_abundance * 100)


def save_glycoform_list(G, outfile):
    """
    Convert the glycoform graph to a list of glycoforms
    with corrected abundances.

    :param nx.DiGraph G: glycoform graph
    :param str outfile: output file name
    :return: nothing
    """

    glycoforms = []
    for n in G:
        glycoforms.append((n.name,
                           n.composition_str(),
                           G.nodes[n]["abundance"],
                           G.nodes[n]["corr_abundance"]))
    (pd.DataFrame(glycoforms, columns=["glycoform", "composition",
                                       "abundance", "corr_abundance"])
       .sort_values("corr_abundance", ascending=False)
       .to_csv(outfile, index=False))


if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s",
                        level=logging.INFO)

    # set up options for the argument parser
    parser = argparse.ArgumentParser(
        description="Correct glycation influence on glycoform abundances.")
    parser.add_argument("-v", "--version",
                        action="version",
                        help="print the version number",
                        version="%(prog)s 1.0")
    parser.add_argument("-o", "--output-format",
                        action="store",
                        help="graph output format, either 'dot' or 'gexf' "
                             "(default: 'dot')",
                        metavar="FORMAT",
                        choices=["dot", "gexf"])
    parser.add_argument("-l", "--glycan-library",
                        action="store",
                        help="CSV file containing a glycan library",
                        required=True)
    parser.add_argument("-f", "--glycoforms",
                        action="store",
                        help="CSV file containing glycoform abundances",
                        required=True)
    parser.add_argument("-g", "--glycation",
                        action="store",
                        help="CSV file containing glycation abundances",
                        required=True)
    args = parser.parse_args()

    glycoforms = pd.read_csv(args.glycoforms, index_col="glycoform")
    glycation = pd.read_csv(args.glycation, index_col="count")

    # check whether column names for glycoform and glycation data are equal
    columns_identical = (
        glycoforms.columns == glycation.columns)  # type: np.ndarray
    if not all(columns_identical):
        error_msg = ["The following column names differ:",
                     "{}\t{}".format(args.glycoforms, args.glycation)]
        for i in np.where(~columns_identical)[0]:
            error_msg.append(
                "{}\t{}".format(glycoforms.columns[i], glycation.columns[i]))
        logging.error("\n".join(error_msg))
        sys.exit(1)

    for col in glycoforms.columns:
        # assemble the glycation graph and correct abundances
        logging.info("Correcting dataset '{}' …".format(col))
        G = calc_glycation_graph(
            glycan_library=args.glycan_library,
            glycoforms=glycoforms[col],
            glycation=glycation[col])
        correct_abundances(G)
        save_glycoform_list(
            G,
            outfile="{}_corr.csv".format(col))

        # output graph data
        if args.output_format == "dot":
            # create more informative labels for the dot format
            for n in G:
                G.nodes[n]["label"] = "{}|{:.2f}|{:.2f}".format(
                    n.name, n.abundance, G.nodes[n]["corr_abundance"])
                G.nodes[n]["shape"] = "record"
            for source, sink in G.edges:
                new_label = "{}: {:.2%}".format(
                    G[source][sink]["label"], G[source][sink]["c"])
                G[source][sink]["label"] = new_label
            nx.nx_pydot.write_dot(G, "{}_corr.gv".format(col))

        elif args.output_format == "gexf":
            nx.write_gexf(G, "{}_corr.gexf".format(col))

    logging.info("… done!")
