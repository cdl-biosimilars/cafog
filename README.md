# cafog - Unbiased relative quantification of protein N-glycosylation

## Requirements

General:

* Python 3.5
* numpy
* pandas
* matplotlib
* PyQt5
* PyQtChart
* Sphinx (for creating the documentation)
* sphinx_rtd_theme (for creating the documentation)
* sphinx-argparse
* uncertainties


For freezing:

* pyinstaller



## Installation

To run cafog directly in source:

* In the `docs` directory, execute `make html` (Linux) or `make.bat html` (Windows) to create the documentation.
* Run cafog from the command line (`python cafog.py`) or with a graphical user interface (`python cafog_gui.py`).

To run cafog from a frozen distribution:

* Extract the contents of the archive and start cafog by executing `cafog`.



## Creating frozen distributions

* Create the documentation as described above.
* Execute `pyinstaller cafog.spec`.
* The folder `dist/cafog` is now a self-contained cafog installation.
