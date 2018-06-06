from collections import Counter
from itertools import product
import re

import numpy as np
import pandas as pd

from formula import Formula, Glycan, PTMComposition


amino_acid_names = {
    "A": ("Alanine", "Ala"),
    "C": ("Cysteine", "Cys"),
    "D": ("Aspartic Acid", "Asp"),
    "E": ("Glutamic Acid", "Glu"),
    "F": ("Phenylalanine", "Phe"),
    "G": ("Glycine", "Gly"),
    "H": ("Histidine", "His"),
    "I": ("Isoleucine", "Ile"),
    "K": ("Lysine", "Lys"),
    "L": ("Leucine", "Leu"),
    "M": ("Methionine", "Met"),
    "N": ("Asparagine", "Asn"),
    "P": ("Proline", "Pro"),
    "Q": ("Glutamine", "Gln"),
    "R": ("Arginine", "Arg"),
    "S": ("Serine", "Ser"),
    "T": ("Threonine", "Thr"),
    "V": ("Valine", "Val"),
    "W": ("Tryptophane", "Trp"),
    "Y": ("Tyrosine", "Tyr"),
    "6": ("Pyroglutamic Acid", "Pyr"),
    "7": ("Hydroxyproline", "Hyp"),
    "J": ("Phosphoserine", "Sep"),
    "Z": ("Phosphothreonine", "Tpo")}

amino_acid_compositions = {
    "A": {"H": 5, "C": 3, "O": 1, "N": 1},
    "C": {"H": 5, "C": 3, "S": 1, "O": 1, "N": 1},
    "D": {"H": 5, "C": 4, "O": 3, "N": 1},
    "E": {"H": 7, "C": 5, "O": 3, "N": 1},
    "F": {"H": 9, "C": 9, "O": 1, "N": 1},
    "G": {"H": 3, "C": 2, "O": 1, "N": 1},
    "H": {"H": 7, "C": 6, "N": 3, "O": 1},
    "I": {"H": 11, "C": 6, "O": 1, "N": 1},
    "K": {"H": 12, "C": 6, "N": 2, "O": 1},
    "L": {"H": 11, "C": 6, "O": 1, "N": 1},
    "M": {"H": 9, "C": 5, "S": 1, "O": 1, "N": 1},
    "N": {"H": 6, "C": 4, "O": 2, "N": 2},
    "P": {"H": 7, "C": 5, "O": 1, "N": 1},
    "Q": {"H": 8, "C": 5, "O": 2, "N": 2},
    "R": {"H": 12, "C": 6, "N": 4, "O": 1},
    "S": {"H": 5, "C": 3, "O": 2, "N": 1},
    "T": {"H": 7, "C": 4, "O": 2, "N": 1},
    "V": {"H": 9, "C": 5, "O": 1, "N": 1},
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},
    "Y": {"H": 9, "C": 9, "O": 2, "N": 1},
    "6": {"H": 5, "C": 5, "O": 2, "N": 1},          # pyroglutamate
    "7": {"H": 7, "C": 5, "O": 2, "N": 1},          # hydroxyproline
    "J": {"H": 6, "C": 3, "O": 5, "N": 1, "P": 1},  # phosphoserine
    "Z": {"H": 8, "C": 4, "O": 5, "N": 1, "P": 1}}  # phosphothreonine


def get_sequence_atoms(sequence, chains=1, disulfide_bonds=0):
    """
    Calculates the atoms in an amino acid sequence.

    :param str sequence: String of amino acids in one-letter format
    :param int chains: number of chains; for each chain, the weight of a
                       water molecule must be added to the total weight
    :param int disulfide_bonds: number of disulfide bonds
    :return: the elemental composition of the chains
             (example: ``{"C": 100; "H": 50; "N": 20}``)
    :rtype: dict
    """

    chains = max(chains, 1)
    composition = {a: 0 for a in "CHNOPS"}
    sequence_aa_composition = Counter(sequence)
    for aa, count in sequence_aa_composition.items():
        for atom, atom_count in amino_acid_compositions[aa].items():
            composition[atom] += count * atom_count
    composition["H"] += chains * 2
    composition["O"] += chains * 1
    composition["H"] -= disulfide_bonds * 2
    return composition


class Protein:
    """
    A class which represents a protein.

    :ivar str name: name of the protein
    :ivar dict amino_acid_composition: {amino acid: count}
    :ivar formula.Formula formula: atomic formula of the protein
    :ivar float mass: the protein's molecular mass

    .. automethod:: __init__
    """

    def __init__(self, name=None, sequence=None, chains=1, disulfides=0):
        """
        Create a new protein instance.

        :param str name: name of the protein
        :param str sequence: sequence of the protein
        :param int chains: number of chains
        :param int disulfides: number of disulfide bonds
        """

        if name is None:
            name = "unnamed protien"
        self.name = name

        if sequence is None:
            sequence = ""
        self.amino_acid_composition = {a: sequence.count(a)
                                       for a in amino_acid_names}
        self.formula = Formula(
            get_sequence_atoms(sequence, chains, disulfides))

    def mass(self, mass_set=None):
        """
        Return the mass of the protein.

        :param dict mass_set: set of atomic masses to use for the calculation
        :return: the mass of the protein
        :rtype float
        """

        return self.formula.mass(mass_set)


class Glycoprotein(Protein):
    """
    A protein with glycans.

    :ivar dict glycosylation_sites: glycosylation sites

    .. automethod:: __init__
    .. automethod:: __str__
    """

    def __init__(self, filename):
        """
        Create a new glycoprotein.

        :param str filename: CSV file containing a glycan library
        :return: nothing
        """

        super().__init__()
        self.glycosylation_sites = {}
        re_chains = re.compile("([^,\s]+)+")

        library = pd.read_csv(filename).fillna("")
        for _, row in library.iterrows():
            chains = re_chains.findall(row.chain)
            if row.composition == "":
                row.composition = None
            glycan = Glycan(name=row.glycan, composition=row.composition)
            for chain in chains:
                try:
                    self.glycosylation_sites[chain].append(glycan)
                except KeyError:
                    self.glycosylation_sites[chain] = [glycan]

    def __str__(self):
        """
        Convert the glycoprotein to a human-readable string.

        :return: a string representing the glycoprotein
        :rtype: str
        """

        result = ["Glycoprotein: " + self.name, "Glycosylation sites:"]
        for site_name, mods in self.glycosylation_sites.items():
            result.append("\tName: " + site_name)
            for mod in mods:
                result.append("\t\t" + str(mod))
        return "\n".join(result)

    def unique_glycoforms(self):
        """
        Calculate all glycoforms unique in terms
        of monosaccharide composition.

        :return: a generator that yields all unique monosaccharide compositions
        :rtype: generator(PTMComposition)
        """

        # determine all glycoforms by calculating the cartesian product
        glycoforms = []
        for combination in list(
                product(*[s for s in self.glycosylation_sites.values()])):
            composition = sum([PTMComposition(ptm.composition)
                               for ptm in combination],
                              PTMComposition())
            name = "/".join([ptm.name for ptm in combination])
            mass = sum([ptm.mass() for ptm in combination])
            abundance = np.prod([ptm.abundance for ptm in combination])
            glycoforms.append((composition, name, mass, abundance))

        # eliminate glycoforms with equal monosaccharide composition
        glycoforms = pd.DataFrame(
            glycoforms, columns=["composition", "name", "mass", "abundance"])
        glycoforms["group_key"] = glycoforms.apply(
            lambda p: hash(p.composition), axis=1)
        glycoforms_agg = (
            glycoforms
            .groupby("group_key")
            .agg({"composition": lambda x: x.iloc[0],
                  "name": lambda n: " or ".join(n),
                  "mass": lambda x: x.iloc[0],
                  "abundance": sum})
            .reset_index(drop=True)
        )
        # noinspection PyUnresolvedReferences
        glycoforms_agg.abundance = (glycoforms_agg.abundance
                                    / glycoforms_agg.abundance.max()
                                    * 100)

        # annotate each representative glycoform by name
        # and relative abundance and create the generator
        for _, row in glycoforms_agg.iterrows():
            row.composition.name = row["name"]
            row.composition.mass = row.mass
            row.composition.abundance = row.abundance
            yield row.composition
