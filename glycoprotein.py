from itertools import product

import numpy as np
import pandas as pd

from formula import Glycan, PTMComposition


class Glycoprotein:
    """
    A protein with glycans.

    :ivar dict glycosylation_sites: glycosylation sites

    .. automethod:: __init__
    .. automethod:: __str__
    """

    def __init__(self, sites, library):
        """
        Create a new glycoprotein.

        :param str library: CSV file containing a glycan library
        :return: nothing
        """

        self.sites = sites
        self.glycan_library = []

        glycan_library = pd.read_csv(library).fillna("")
        for _, row in glycan_library.iterrows():
            if row.composition == "":
                row.composition = None
            self.add_glycan(name=row.glycan, composition=row.composition)

    def __str__(self):
        """
        Convert the glycoprotein to a human-readable string.

        :return: a string representing the glycoprotein
        :rtype: str
        """

        result = ["Glycoprotein with {} sites and glycans".format(self.sites)]
        for glycan in self.glycan_library:
            result.append("\t{}".format(glycan))
        return "\n".join(result)

    def add_glycan(self, name, composition=None):
        """
        Add a glycan to the library.

        :param str name: name of the glycan
        :param str composition: monosaccharide composition
        :return: nothing
        """

        self.glycan_library.append(Glycan(name=name, composition=composition))

    def unique_glycoforms(self):
        """
        Calculate all glycoforms unique in terms
        of monosaccharide composition.

        :return: a generator that yields all unique monosaccharide compositions
        :rtype: generator(PTMComposition)
        """

        # determine all glycoforms by calculating the cartesian product
        glycoforms = []
        for combination in product(self.glycan_library, repeat=self.sites):
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

# gp = Glycoprotein(2, "data/glycans.csv")
# # print(gp)
# for c in gp.unique_glycoforms():
#     print(c) TODO remove
