import argparse
import re
import sys

from multiset import FrozenMultiset
import networkx as nx
import numpy as np
import pandas as pd

from formula import PTMComposition
from protein import Glycoprotein


def calc_glycation_graph(protein, glycoforms, glycation):
    """
    Assemble the glycoform graph from peptide mapping
    and glycation frequency data.

    :param str protein: JSON file describing protein parameters
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
    exp_abundances["ptm_set"] = exp_abundances.glycoform.apply(
        lambda v: FrozenMultiset(v.split("/")))
    exp_abundances = exp_abundances.set_index("ptm_set").abundance
    # exp_abundances = (exp_abundances
    #                   / exp_abundances.max() * 100)  # type: pd.Series TODO remove?

    # dict mapping PTM differences to abundances
    delta_ptm = {PTMComposition({"Hex": count}): abundance
                 for count, abundance in glycation.iteritems()
                 if count > 0}

    G = nx.DiGraph()
    for glycoform in Glycoprotein(protein).unique_glycoforms():
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

    # max_corr_abundance = 0.0
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
        # max_corr_abundance = max(max_corr_abundance, corr_abundance)
        total_abundance += corr_abundance
        G.nodes[n]["corr_abundance"] = float(corr_abundance)

    # second pass: normalize abundance to 100
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
    # set up options for the argument parser
    parser = argparse.ArgumentParser(description="")
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
    parser.add_argument("-p", "--protein",
                        action="store",
                        help="JSON file specifying protein parameters",
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
        print("Error: The following column names differ:\n{}\t{}"
              .format(args.glycoforms, args.glycation))
        for i in np.where(~columns_identical)[0]:
            print("{}\t{}"
                  .format(glycoforms.columns[i], glycation.columns[i]))
        sys.exit(1)

    for col in glycoforms.columns:
        # assemble the glycation graph and correct abundances
        print("Correcting dataset '{}' …".format(col))
        G = calc_glycation_graph(
            protein=args.protein,
            glycoforms=glycoforms[col],
            glycation=glycation[col])
        correct_abundances(G)
        save_glycoform_list(
            G,
            outfile="{}_corr.csv".format(col))

        # output graph data
        if args.output_format == "dot":
            # create more infrmative labels for the dot format
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

    print("… done!")
