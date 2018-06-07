import re

import pandas as pd
import pandas.core.series


class Glycan:
    """
    An N- or O-glycan, described by its monosaccharide composition.

    :ivar name
    :ivar composition
    :ivar abundance

    .. automethod:: __init__
    .. automethod:: __str__
    """

    def __init__(self, name=None, composition=None, abundance=0.0):
        """
        Create a new glycan.

        :param str name: name of the glycan; converted to a composition
                         if composition is None and name is a valid Zhang name
        :param str composition: comma-separated list of monosaccharides
        :param float abundance: relative abundance
        """

        self.name = name

        if composition is None:
            composition = Glycan.extract_composition(name)
        self.composition = composition

        self.abundance = abundance

    @staticmethod
    def extract_composition(glycan):
        """
        Convert a glycan abbreviation in Zhang nomenclature (e.g., "A2G1F")
        to a composition string (e.g., "4 Hex, 3 HexNAc, 1 Fuc").

        :param str glycan: a glycan abbreviation
        :return: a composition string
        :rtype: str
        :raises ValueError: if the conversion fails
        """

        re_zhang_glycan = re.compile(r"""
            ^
            (?:A(?P<A>\d)+)?    # antennas
            (?:Sg(?P<Sg>\d)+)?  # Neu5Gc
            (?:S(?P<S>\d)+)?    # Neu5Ac
            (?:Ga(?P<Ga>\d)+)?  # alpha-Gal
            (?:G(?P<G>\d)+)?    # Gal
            (?:M(?P<M>\d)+)?    # Man
            (?P<F>F)?           # core Fuc
            (?P<B>B)?           # bisecting GlcNAc
            $
            """, re.VERBOSE)

        try:
            g = re_zhang_glycan.match(glycan).groupdict()
        except (AttributeError, TypeError):
            if glycan in ["non-glycosylated", "unglycosylated", "null"]:
                g = {}
            else:
                raise ValueError("Invalid glycan name: '{}'".format(glycan))
        if all(v is None for v in g.values()):
            counts = [0, 0, 0, 0, 0]  # input was empty string
        else:
            if g["M"] is None:
                g["M"] = 3  # there are always three Man
                if g["A"] is None:
                    g["A"] = 2  # handle abbreviations 'Gn' and 'GnF'
            for k, v in g.items():
                try:
                    g[k] = int(v)
                except TypeError:
                    g[k] = 0  # for elements that are not found
                except ValueError:
                    g[k] = 1  # for F and B
            # noinspection PyTypeChecker
            counts = [
                g["Sg"] + g["S"] + 2 * g["Ga"] + g["G"] + g["M"],  # Hex
                g["A"] + 2 + (1 if g["B"] else 0),  # HexNAc
                g["S"],  # Neu5Ac
                g["Sg"],  # Neu5Gc
                1 if g["F"] else 0]  # Fuc
        monosaccharides = ["Hex", "HexNAc", "Neu5Ac", "Neu5Gc", "Fuc"]
        return ", ".join("{} {}".format(c, m)
                         for m, c in zip(monosaccharides, counts)
                         if c > 0)

    def __str__(self):
        """
        Convert the glycan to a human-readable string.

        :return: a string representing the glycan
        :rtype: str
        """
        return "{} ({}) = {}".format(
            self.name, self.composition, self.abundance)


class PTMComposition:
    """
    A composition of post-translational modifications (PTMs).

    :ivar str name: name of this composition
    :ivar float abundance: abundance of the respective proteoform
    :ivar pd.Series composition: composition of PTMs

    .. automethod:: __init__
    .. automethod:: __add__
    .. automethod:: __sub__
    .. automethod:: __str__
    .. automethod:: __repr__
    .. automethod:: __hash__
    .. automethod:: __eq__
    .. automethod:: __ne__
    .. automethod:: __neg__
    """

    def __init__(self, mods=None, name=None, abundance=None):
        """
        Creates a new PTM composition.

        :param mods: one of the following:

               * a pd.Series like ``pd.Series(dict(Hex=1, HexNAc=2))``
               * a dict like ``{"Hex": 1; "HexNAc": 2}``
               * a string like ``"1 Hex, 2 HexNAc"``
               * ``None``, which generates an empty composition
        :param str name: name of the composition (default: ``""``)
        :param float abundance: abundance of the respective proteoform
                                (default: 0.0)
        :raises TypeError: if composition is of an unsupported type
        """

        if name is None:
            name = ""
        self.name = name

        if abundance is None:
            abundance = 0.0
        self.abundance = abundance

        if mods is None:
            self.composition = pd.Series()
        elif isinstance(mods, pd.core.series.Series):
            self.composition = mods.copy()
        elif isinstance(mods, str):
            self.composition = PTMComposition.extract_composition(mods)
        elif isinstance(mods, dict):
            self.composition = pd.Series(mods)
        else:
            raise TypeError("Type {} not supported".format(type(mods)))
        self.composition = self.composition[self.composition != 0]

    @staticmethod
    def extract_composition(mods):
        """
        Extract the PTM composition from a string like "1 Hex, 2 HexNAc, Fuc".

        :param str mods: string describing a PTM composition
        :return: a series describing the PTM
        :rtype: pd.Series
        """

        re_ptm_list = re.compile(r"""
            (\d*)     # optional count
            \s*       # optional space
            ([\w-]+)  # monosaccharide name
            (?:,|$)   # comma separator
            """, re.VERBOSE)

        composition = {}
        for c, m in re_ptm_list.findall(mods):
            try:
                c = int(c)
            except ValueError:  # if count is absent
                c = 1
            composition[m] = c

        return pd.Series(composition)

    def composition_str(self):
        """
        Returns a string representing the composition of self
        (e.g., "1 Hex, 2 HexNAc"), or "[no PTMs]" for an empty composition.

        :return: composition string
        :rtype: str
        """

        if self.composition.empty:
            return "[no PTMs]"
        else:
            return ", ".join(["{:d} {:s}".format(c, m)
                              for m, c in self.composition.items()])

    def __add__(self, other):
        """
        Add two PTM compositions:

          * Counts will be added for each PTM.
          * Name will be a concatenation of the summands' names
            separated by "+".
          * Abundance values will be added.

        :param PTMComposition other: PTM composition to be added
        :return: the sum of self and other
        :rtype: PTMComposition
        """

        return PTMComposition(
            mods=self.composition
                     .add(other.composition, fill_value=0)
                     .astype(int),
            name=self.name + "+" + other.name,
            abundance=self.abundance + other.abundance)

    def __sub__(self, other):
        """
        Subtract two PTM compositions:

          * Counts will be subtracted for each PTM.
          * Name will be a concatenation of the summands' names
            separated by "-".
          * Abundance values will be subtracted.

        :param PTMComposition other: PTM composition to be subtracted
        :return: self minus other
        :rtype: PTMComposition
        """

        return PTMComposition(
            mods=self.composition
                     .sub(other.composition, fill_value=0)
                     .astype(int),
            name=self.name + "-" + other.name,
            abundance=self.abundance - other.abundance)

    def __str__(self):
        """
        Convert the PTM to a human-readable string.

        :return: a string representing the PTM, like
                 PTMComposition(1 Hex, 2 HexNAc, name=glycan1, abundance=10)
        :rtype: str
        """

        return "PTMComposition({}, name={}, abundance={:.2f})".format(
            self.composition_str(),
            self.name,
            self.abundance)

    def __repr__(self):
        """
        Representation of the PTM object.

        :return: the string representation of self
        :rtype: str
        """

        return self.__str__()

    def __hash__(self):
        """
        Calculates the hash value, which is the hash of the tuple of tuples
        describing self's composition.

        :return: hash value
        :rtype: int
        """

        return hash(tuple(self.composition.iteritems()))

    def __eq__(self, other):
        """
        Determine equality.

        :param PTMComposition other: object to which self is compared
        :return: True if self and other have equal compositions
        :rtype: bool
        """

        return dict(self.composition) == dict(other.composition)

    def __ne__(self, other):
        """
        Determine inequality.

        :param PTMComposition other: object to which self is compared
        :return: True if self and other have different compositions
        :rtype: bool
        """

        return not self.__eq__(other)

    def __neg__(self):
        """
        Unary negation:

          * Each PTM count changes sign.
          * Name will be prefixed by "-".
          * Abundance will change sign.

        :return: negated self
        :rtype: PTMComposition
        """

        return PTMComposition(
            mods=-self.composition,
            name="-" + self.name,
            abundance=-self.abundance)

if __name__ == "__main__":
    print(Glycan("A2S2G0F"))
