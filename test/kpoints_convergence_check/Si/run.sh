#!/bin/bash
#

rm espresso.* qe_kpts.log kpoint_convergence.png

python ../../../kpoints_convergence_check/kpoints_convergence.py --cif_file "Si.cif" --qe_path "/Users/mandal13/Downloads/softwares/q-e-qe-7.0/build/bin/pw.x" --max_k_multiplier 6 --json_pp_file "Si.json"
