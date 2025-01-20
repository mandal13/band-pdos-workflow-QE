import os
import json
import argparse
from ase.io import read
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
from matplotlib import pyplot as plt
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 1.2

def scf_calculation(system: Atoms, profile: EspressoProfile, pseudopotentials: dict, input_data: dict, kpts):
    """
    Perform a self-consistent field (SCF) calculation.

    Parameters:
        system (Atoms): ASE Atoms object representing the system.
        profile (EspressoProfile): Quantum ESPRESSO execution profile.
        pseudopotentials (dict): Dictionary of pseudopotentials.
        input_data (dict): Input parameters for the calculation.
        kpts (tuple or None): K-point sampling.

    Returns:
        float: Total energy of the system.
    """
    calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=input_data, kpts=kpts)
    system.calc = calc
    return system.get_total_energy()

def parse_args():
    """
    Parse command-line arguments.

    Returns:
        Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Check convergence of total energy with respect to cutoff energy using Quantum ESPRESSO.")
    parser.add_argument("--cif_file", type=str, required=True, help="Input CIF file.")
    parser.add_argument("--qe_path", type=str, required=True, help="Path to the Quantum ESPRESSO executable.")
    parser.add_argument("--pseudo_dir", type=str, default="./", help="Directory containing pseudopotentials.")
    parser.add_argument("--np", type=int, default=1, help="Number of processors to use.")
    parser.add_argument("--start_ecutwfc", type=int, default=40, help="Starting value of ecutwfc (Ry).")
    parser.add_argument("--end_ecutwfc", type=int, default=150, help="Ending value of ecutwfc (Ry).")
    parser.add_argument("--step_ecutwfc", type=int, default=10, help="Step size for ecutwfc (Ry).")
    parser.add_argument("--ecutrho_factor", type=int, default=4, help="Factor to calculate ecutrho from ecutwfc.")
    parser.add_argument("--json_pp_file", type=str, default="pp.json", help="JSON file containing pseudopotential info.")
    parser.add_argument("--log_file", type=str, default="qe.log", help="Log file to save energies.")
    parser.add_argument("--prefix", type=str, default="pw", help="Prefix for output files.")
    parser.add_argument("--kpts", type=str, default=None, help="K-point sampling in the format '(kx, ky, kz)'.")
    parser.add_argument("--convergence_criteria", type=float, default=1e-3, help="Energy difference threshold for convergence (eV/atom)")
    return parser.parse_args()

def main():
    """
    Main function to perform convergence check of total energy with respect to cutoff energy.
    """
    # Parse arguments
    args = parse_args()

    # Read input parameters
    system = read(args.cif_file)

    # Number of atoms in the system
    num_atoms = len(system)
    
    # Parse k-point sampling
    kpts = tuple(map(int, args.kpts.strip("() ").split(","))) if args.kpts else None
    if kpts and len(kpts) != 3:
        raise ValueError("kpts must have exactly three values.")

    # Load pseudopotentials
    with open(args.json_pp_file, 'r') as f:
        pp_dict = json.load(f)

    # Prepare Quantum ESPRESSO profile
    command = f'mpiexec -n {args.np} {args.qe_path}'
    profile = EspressoProfile(command=command, pseudo_dir=args.pseudo_dir)

    # Set up input data template
    input_data = {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'etot_conv_thr': 1e-5,
        'forc_conv_thr': 1e-4,
        'ecutwfc': args.start_ecutwfc,
        'ecutrho': args.ecutrho_factor * args.start_ecutwfc,
        'input_dft': 'pbe',
        'occupations': 'smearing',
        'degauss': 0.01,
        'smearing': 'gaussian',
        'conv_thr': 1e-8,
        'mixing_beta': 0.35,
        'startingwfc': 'random',
        'outdir': './tmp',
        'prefix': args.prefix
    }

    # Perform calculations for a range of ecutwfc values
    ecuts = range(args.start_ecutwfc, args.end_ecutwfc + 1, args.step_ecutwfc)
    energies = []

    with open(args.log_file, 'w') as log_file:
        for ecut in ecuts:
            input_data['ecutwfc'] = ecut
            input_data['ecutrho'] = args.ecutrho_factor * ecut

            # Perform SCF calculation
            energy = scf_calculation(system, profile, pp_dict, input_data, kpts)/num_atoms
            energies.append(energy)

            # Log results
            log_file.write(f'{ecut} {energy}\n')

    # Clean up temporary files
    os.system(f'rm -rf {input_data["outdir"]}')

    # Determine convergence
    energy_diffs = [abs(energies[i+1] - energies[i]) for i in range(len(energies) - 1)]
    optimal_index = next((i for i, diff in enumerate(energy_diffs) if diff < args.convergence_criteria), len(energies) - 1)
    optimal_cutoff = ecuts[optimal_index]

    # Plot results
    plt.plot(ecuts, energies, 'o-', label='Total Energy')

    # Highlight optimal cutoff energy
    plt.axvline(optimal_cutoff, color='red', linestyle='--', label=f'Optimal Cutoff Energy: {optimal_cutoff} Ry')

    plt.xlabel('Cutoff Energy (Ry)')
    plt.ylabel('Total Energy (eV/atom)')
    plt.title('Total Energy vs Cutoff Energy')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig('cutoff_convergence.png', dpi=300, bbox_inches='tight')
    print(f"Convergence plot saved as 'cutoff_convergence.png'. Optimal cutoff energy: {optimal_cutoff} Ry")

if __name__ == '__main__':
    main()
