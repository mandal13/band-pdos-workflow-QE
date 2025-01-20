#!/bin/bash

rm espresso.* qe.log cutoff_convergence.png

#python ../../../ecut_convergence_check/cutoff_convergence.py --cif_file "Si.cif" --qe_path "/Users/mandal13/Downloads/softwares/q-e-qe-7.0/build/bin/pw.x" 

python ../../../ecut_convergence_check/cutoff_convergence.py --cif_file "Si.cif" --qe_path "/Users/mandal13/Downloads/softwares/q-e-qe-7.0/build/bin/pw.x"  --kpts '(1,2,1)' --convergence_criteria 1e-4

