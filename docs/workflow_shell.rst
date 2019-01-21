********************************
cafog - Workflow using the shell
********************************

.. argparse::
   :filename: ../cafog.py
   :func: setup_parser
   :prog: cafog.py

   -f --glycoforms
       Required columns:
       
       1. glycoform names (``glycan 1/glycan 2/…/glycan n``)
       2. abundances
       3. experimental errors


   -g --glycation
       The required columns for CSV files are described here.

       Required columns:
       
       1. glycation counts
       2. abundances
       3. experimental errors


   -l --glycan-library
       Required columns:

       1. glycan names
       2. monosaccharide compositions (e.g., ``1 Hex, 2 HexNAc, 3 Fuc``)



* A header is not allowed, but any line starting with ``#`` will be discarded.
* The monosaccharide composition for glycans whose name follows the Zhang nomenclature will be deduced automatically. The composition for other glycans must be provided in the glycan library.
