#!/bin/bash

rm -rf espresso.* qe.log tmp/
#
python ../../../ecut_convergence_check/cutoff_convergence.py --cif_file "ZnS.cif" --qe_path "/Users/mandal13/Downloads/softwares/q-e-qe-7.0/build/bin/pw.x" --np 2 --start_ecutwfc 20 --end_ecutwfc 50
