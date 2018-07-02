import logging

import pandas as pd
from uncertainties import ufloat


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
