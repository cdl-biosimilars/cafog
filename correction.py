import logging
import re

from multiset import FrozenMultiset
import networkx as nx
import pandas as pd
from uncertainties import ufloat

from glycan import PTMComposition
from glycoprotein import Glycoprotein


def calc_glycation_graph(glycan_library: pd.DataFrame,
                         glycoforms: pd.Series,
                         glycation: pd.Series) -> nx.DiGraph:
    """
    Assemble the glycoform graph from peptide mapping
    and glycation frequency data.

    :param pd.DataFrame glycan_library: dataframe containing a glycan library
    :param pd.Series glycoforms: list of glycoforms with abundances/errors
    :param pd.Series glycation: list of glycations with abundances/errors
    :raises ValueError: if a glycan with unknown monosaccharide composition
                        is added
    :return: the glycation DAG
    :rtype: nx.DiGraph
    """

    # regex for extracting the first glycoform from a string like
    # "A2G0F/A2G1F or A2G1F/A2G0F"
    re_first_glycoform = re.compile("([^\s]*)")

    # series with a monosaccharide set as index
    # and abundances as values
    exp_abundances = glycoforms.reset_index()
    exp_abundances["sugar_set"] = exp_abundances["index_col"].apply(
        lambda v: FrozenMultiset(v.split("/")))
    exp_abundances = exp_abundances.set_index("sugar_set")["abundance"]

    # dict mapping hexose differences to abundances
    delta_ptm = {PTMComposition({"Hex": count}): abundance / 100
                 for count, abundance in glycation.iteritems()
                 if count > 0}

    gp = Glycoprotein(sites=2, library=glycan_library)
    glycoform_glycans = set()
    for v in exp_abundances.index.values:
        glycoform_glycans |= set(v)

    if glycan_library is None:
        # fill the glycan library from glycans in glycoforms
        logging.info("No glycan library specified. "
                     "Extracting glycans from glycoforms …")
        for g in glycoform_glycans:
            try:
                gp.add_glycan(g)
            except ValueError as e:
                raise e
    else:
        # compare monosaccharide set of glycan library and glycoforms
        # add glycans that only appear in the list of glycoforms to the library
        # but this only works if they have a valid name
        library_glycans = set([n.name for n in gp.glycan_library])

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
                    raise e

    G = nx.DiGraph()
    for glycoform in gp.unique_glycoforms():
        # print("x")
        # get the experimental abundance of a glycoform
        # use a default value of 0±0 if unavailable
        abundance = ufloat(0, 0)
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
            abundance=abundance,
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
            G.add_edge(source, sink, label=d.composition_str(), c=c)
    return G


def correct_abundances(G: nx.DiGraph) -> None:
    """
    Correct abundances in the glycoform graph.

    :param nx.DiGraph G: a glycoform graph
    :return: nothing, G is modified in place
    :rtype: None
    """

    # calculate corrected abundance for each node from source to sink
    for n in nx.topological_sort(G):
        in_abundance = 0.0
        for pred in G.predecessors(n):
            in_abundance += G.nodes[pred]["corr_abundance"] * G[pred][n]["c"]
        out_c = 0.0
        for succ in G.successors(n):
            out_c += G[n][succ]["c"]
        corr_abundance = (n.abundance - in_abundance) / (1 - out_c)
        G.nodes[n]["corr_abundance"] = corr_abundance


def save_glycoform_list(G: nx.DiGraph,
                        filename: str) -> None:
    """
    Convert the glycoform graph to a list of glycoforms
    with corrected abundances.

    :param nx.DiGraph G: glycoform graph
    :param str filename: name of the output file
    :return: nothing
    :rtype: None
    """

    glycoforms = []
    composition = []
    for n in G:
        glycoforms.append((n.name,
                           G.nodes[n]["abundance"].nominal_value,
                           G.nodes[n]["abundance"].std_dev,
                           G.nodes[n]["corr_abundance"].nominal_value,
                           G.nodes[n]["corr_abundance"].std_dev))
        composition.append(n.composition)
    composition = pd.concat(composition, axis=1)  # type: pd.DataFrame
    (pd.DataFrame(glycoforms, columns=["glycoform", "abundance",
                                       "abundance_error", "corr_abundance",
                                       "corr_abundance_error"])
       .join(composition.T)
       .sort_values("corr_abundance", ascending=False)
       .to_csv("{}_corr.csv".format(filename), index=False))


def save_graph(G: nx.DiGraph,
               filename: str,
               output_format: str) -> None:
    """
    Export the glycation graph in a graph file format.

    :param nx.DiGraph G: glycation graph to be exported
    :param str filename: name of the output file
    :param str output_format: either ``"dot"`` or ``"gexf"``
    :return: nothing
    :rtype: None
    """

    if output_format == "dot":
        # create more informative labels for the dot format
        for n in G:
            G.nodes[n]["label"] = "{}|{:.2f}|{:.2f}".format(
                n.name, n.abundance, G.nodes[n]["corr_abundance"])
            G.nodes[n]["shape"] = "record"
        for source, sink in G.edges:
            new_label = "{}: {:.2%}".format(
                G[source][sink]["label"], G[source][sink]["c"])
            G[source][sink]["label"] = new_label
        nx.nx_pydot.write_dot(G, "{}_corr.gv".format(filename))

    elif output_format == "gexf":
        # split each ufloat attribute into two float attributes
        for n in G:
            G.nodes[n]["abundance_error"] = float(
                G.nodes[n]["abundance"].std_dev)
            G.nodes[n]["abundance"] = float(
                G.nodes[n]["abundance"].nominal_value)
            G.nodes[n]["corr_abundance_error"] = float(
                G.nodes[n]["corr_abundance"].std_dev)
            G.nodes[n]["corr_abundance"] = float(
                G.nodes[n]["corr_abundance"].nominal_value)
        for source, sink in G.edges:
            G[source][sink]["c_error"] = float(G[source][sink]["c"].std_dev)
            G[source][sink]["c"] = float(G[source][sink]["c"].nominal_value)
        nx.write_gexf(G, "{}_corr.gexf".format(filename))


def read_clean_datasets(filename: str) -> pd.Series:
    """
    Read input datasets (glycoforms, glycations) and prepare for analysis,
    i.e., generate a single column containing abundances with uncertainties.

    :param str filename: name of the file containing the dataset
    :return: a series called "abundance" containing a value with uncertainty
             and an index named "index_col"
    :rtype: pd.Series
    :raises ValueError: if the input dataset contains too few columns
    """

    df = pd.read_csv(filename, comment="#", index_col=0, header=None)
    col_count = df.shape[1]

    if col_count == 0:  # too few columns
        raise ValueError("{} contains too few columns.".format(filename))
    elif col_count == 1:  # add error column
        logging.warning(
            "{} lacks a column containing errors. Assuming errors of zero."
            .format(filename))
        df["auto_error_column"] = 0
    elif col_count > 2:  # remove surplus columns
        logging.warning(
            "{} contains {} additional columns, which will be ignored."
            .format(filename, col_count-2))
        df = df.iloc[:, :3]
    df.index.name = "index_col"
    df["abundance"] = df.apply(lambda r: ufloat(r.iloc[0], r.iloc[1]), axis=1)
    return df["abundance"]


def read_library(filename: str=None) -> pd.DataFrame:
    """
    Read glycan library and prepare for analysis.

    :param str filename: name of the file containing the glycan library
    :return: a dataframe twith wo columns
            containing glycan names and compositions, respectively
    :rtype: pd.DataFrame
    :raises ValueError: if the input dataset contains too few columns
    """

    df = pd.read_csv(filename, comment="#", header=None)
    col_count = df.shape[1]

    if col_count < 2:  # too few columns
        raise ValueError("{} contains too few columns.".format(filename))
    elif col_count > 2:  # remove surplus columns
        logging.warning(
            "{} contains {} additional columns, which will be ignored."
            .format(filename, col_count-2))
        df = df.iloc[:, :3]
    return df
