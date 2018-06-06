"""
Functions/classes for handling atomic masses.
"""

import re

import numpy as np
import pandas as pd
import pandas.core.series


class Formula:
    """
    A molecular formula comprising the elements C, H, N, O, P and S.

    :ivar pd.Series _composition: atomic composition
    :cvar default_mass_set: a default set of atomic masses (Zhang)

    .. automethod:: __init__
    .. automethod:: __add__
    .. automethod:: __sub__
    .. automethod:: __mul__
    .. automethod:: __rmul__
    .. automethod:: __bool__
    .. automethod:: __str__
    """

    @staticmethod
    def extract_composition(formula):
        """
        Converts a formula string to a Series.

        :param str formula: space-separated collection of elements,
                            each of which may include a count
                            (example: "C50 H100 N-3 Cl")
        :return: a Series with the elements as index and their counts as values
                 (example: C: 50, H: 100, N: -3, Cl: 1)
        :rtype: pd.Series
        :raises ValueError: if the composition string is invalid
        """

        pattern = re.compile(r"([A-Z][a-z]?)(-?\d*)$")
        composition = {}
        for element in formula.split():
            match = pattern.match(element)
            if match:
                atom, count = match.groups()
                if count:
                    composition[atom] = int(count)
                else:
                    composition[atom] = 1
            else:
                raise ValueError("Invalid formula part: {}".format(element))
        return pd.Series(composition)

    default_mass_set = {
        "C": 12.010790,
        "H": 1.007968,
        "N": 14.006690,
        "O": 15.999370,
        "P": 30.973763,
        "S": 32.063900,
        "Cl": 35.45,
        "Na": 22.99
    }

    def __init__(self, composition=None):
        """
        Create a new Formula.

        :param composition: one of the following:

               * a pd.Series like ``pd.Series(dict(C=6, H=12, O=6))``
               * a dict like ``{"C": 6; "H": 12; "O": 6}``
               * a string like ``"C6 H12 O6".``
        :raises ValueError: if composition is of a different type
        """

        if isinstance(composition, dict):
            self._composition = pd.Series(composition)
        elif isinstance(composition, pd.core.series.Series):
            self._composition = composition.copy()
        elif isinstance(composition, str):
            self._composition = Formula.extract_composition(composition)
        elif composition is None:
            self._composition = pd.Series()
        else:
            raise ValueError("Could not create formula from "
                             + str(composition))
        self._composition = self._composition.replace(np.nan, 0).astype(int)

    def mass(self, mass_set=None):
        """
        Calculate and return the mass of the formula.

        :param dict mass_set: {atom symbol: atomic mass} dict;
                              None if the default mass set should be used
        :return: the mass of the Formula
        :rtype: float
        """

        if mass_set is None:
            mass_set = Formula.default_mass_set
        mass_set = pd.Series(mass_set)

        return (self._composition
                .mul(mass_set.reindex(self._composition.index))
                .sum())

    def __add__(self, other):
        """
        Sum two formulas.

        :param Formula other: formula to be added
        :return: the sum of self and other
        :rtype: Formula
        """
        return Formula(self._composition
                       .add(other._composition, fill_value=0)
                       .astype(int))

    def __sub__(self, other):
        """
        Subtract two formulas.

        :param Formula other: formula to be subtracted
        :return: self minus other
        :rtype: Formula
        """
        return Formula(self._composition
                       .sub(other._composition, fill_value=0)
                       .astype(int))

    def __mul__(self, factor):
        """
        Multiply each atom by a factor.

        :param int factor: multiplicative factor
        :return: the multiplied formula
        :rtype: Formula
        """
        if type(factor) == int:
            return Formula(self._composition * factor)

    def __rmul__(self, factor):
        """
        Allow the formula to be the right operand in a multiplication.

        :param int factor: multiplicative factor
        :return: the multiplied formula
        :rtype: Formula
        """
        return self.__mul__(factor)

    def __bool__(self):
        """
        Implement truth value testing.

        :return: True if the formula contains at least one atom;
                 False otherwise
        :rtype: bool
        """

        return bool(self._composition.sum())

    def __str__(self):
        """
        Convert the formula to a human-readable string.

        :return: a string representing the formula
        :rtype: str
        """

        result = []
        for element, count in self._composition.items():
            if count == 1:
                result.append(element)
            elif count != 0:
                result.append("{:s}{:d}".format(element, count))
        return " ".join(result)


class Modification:
    """
    Class implementing a general post-translational modification.

    :ivar name
    :ivar composition
    :ivar formula
    :ivar abundance
    :ivar site

    .. automethod:: __init__
    .. automethod:: __str__
    """

    def __init__(self, name=None, composition=None, formula=None,
                 abundance=None, site=None):
        """
        Create a new modification.

        :param str name: name of the modification
        :param str composition: composition of the modification
        :param str formula: its formula
        :param abundance: its abundance; either a number
                          or a {count: abundance} dict
        :type abundance: float or dict
        :param str site: its site
        """

        self.name = name
        self.composition = composition
        self.formula = Formula(formula)
        self.site = site

        try:
            self.abundance = {}
            for count, abundance in abundance.items():
                # noinspection PyUnresolvedReferences
                self.abundance[int(count)] = float(abundance)
        except AttributeError:
            try:
                self.abundance = float(abundance)
            except TypeError:
                self.abundance = 0.0

    def __str__(self):
        """
        Convert the modification to a human-readable string.

        :return: a string representing the modification
        :rtype: str
        """
        return "{} ({}) = {}".format(self.name, self.formula, self.abundance)

    def mass(self, mass_set=None):
        """
        Return the mass of the modification.

        :param dict mass_set: set of atomic masses to use for the calculation
        :return: the mass of te modification
        :rtype float
        """

        return self.formula.mass(mass_set)


class Glycan(Modification):
    """
    Class implementing an N- or O-glycan.

    :cvar monosaccharide_formula: {monosaccharide: formula} dict
    :ivar name
    :ivar composition
    :ivar formula
    :ivar abundance
    :ivar site

    .. automethod:: __init__
    .. automethod:: __str__
    """

    monosaccharide_formula = dict(
        Hex=Formula("C6 H10 O5"),
        HexNAc=Formula("C8 H13 O5 N1"),
        Neu5Ac=Formula("C11 H17 O8 N1"),
        Neu5Gc=Formula("C11 H17 O9 N1"),
        Fuc=Formula("C6 H10 O4")
    )

    def __init__(self, name=None, composition=None,
                 formula=None, abundance=None, site=None):
        """
        Create a new glycan.

        :param str name: name of the glycan; converted to a composition
                         if composition is None and name is a valid Zhang name
        :param str composition: comma-separated list of monosaccharides
        :param str formula: formula of the glycan; if empty, automatically
                            deduced from the composition
        :param float abundance: relative abundance
        :param str site: glycosylation site
        """

        super().__init__(name=name, composition=composition, formula=formula,
                         abundance=abundance, site=site)
        if composition is None:
            try:
                self.composition = Glycan.extract_composition(name)
            except ValueError:
                pass

        if formula is None:
            self.formula = Glycan.extract_formula(self.composition)
        else:
            self.formula = Formula(formula)

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

    @staticmethod
    def extract_formula(glycan):
        """
        Convert a glycan described by its monosaccharide composition
        (e.g., "3 Hex, 4 HexNAc, 1 Fuc")
        to a Formula (e.g., "C56 H98 N4 O39").

        :param str glycan: a comma-separated list of monosaccharies
        :return: the corresponding formula
        :rtype: Formula
        :raises ValueError: if the conversion fails
        """

        re_monosaccharide_list = re.compile(r"""
            (\d*)     # optional count
            \s*       # optional space
            ([\w-]+)  # monosaccharide name
            (?:,|$)   # comma separator
            """, re.VERBOSE)

        formula = Formula()
        try:
            for (c, m) in re_monosaccharide_list.findall(glycan):
                try:
                    count = int(c)
                except ValueError:
                    count = 1
                formula += count * Glycan.monosaccharide_formula[m]
        except TypeError:
            pass

        return formula

    def __str__(self):
        """
        Convert the glycan to a human-readable string.

        :return: a string representing the glycan
        :rtype: str
        """
        return "{} ({} | {}) = {}".format(self.name, self.composition,
                                          self.formula, self.abundance)


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
    print(Glycan(name="unglycosylated", formula=None))
