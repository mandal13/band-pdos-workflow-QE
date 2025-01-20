import os
import json
import argparse
from ase.io import read
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
from matplotlib import pyplot as plt
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 1.2

def scf_calculation(system: Atoms, profile: EspressoProfile, pseudopotentials: dict, input_data: dict, kpts: tuple):
    """
    Perform a single-point SCF calculation and return the total energy.

    Parameters:
        system (Atoms): ASE Atoms object representing the system.
        profile (EspressoProfile): Profile for Quantum Espresso execution.
        pseudopotentials (dict): Mapping of elements to pseudopotential files.
        input_data (dict): Input parameters for Quantum Espresso.
        kpts (tuple): K-point sampling.

    Returns:
        float: Total energy in eV.
    """
    calc = Espresso(profile=profile, pseudopotentials=pseudopotentials, input_data=input_data, kpts=kpts)
    system.calc = calc
    energy = system.get_total_energy()
    return energy

def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description='Perform k-point convergence test for a material using Quantum Espresso.')
    parser.add_argument("--cif_file", type=str, required=True, help="Input CIF file")
    parser.add_argument("--qe_path", type=str, required=True, help="Path to the Quantum Espresso executable")
    parser.add_argument("--pseudo_dir", default="./", type=str, help="Path to the directory containing the pseudopotentials")
    parser.add_argument("--np", default=1, type=int, help="Number of processors to use")
    parser.add_argument("--json_pp_file", default="pp.json", type=str, help="Path to the JSON file containing pseudopotential info.")
    parser.add_argument("--log_file", default="qe_kpts.log", type=str, help="Path to the log file.")
    parser.add_argument("--prefix", default="pw", type=str, help="Prefix for the output files.")
    parser.add_argument("--max_k_multiplier", type=int, default=5, help="Maximum multiplier for k-point generation")
    parser.add_argument("--ecutwfc", type=int, default=50, help="Cutoff energy for wavefunctions")
    parser.add_argument("--ecutrho", type=int, default=200, help="Cutoff energy for charge density")
    parser.add_argument("--convergence_criteria", type=float, default=1e-3, help="Energy difference threshold for convergence (eV/atom)")

    return parser.parse_args()

def generate_kpoints(box_lengths, max_multiplier):
    """
    Generate k-point sampling based on box lengths and a multiplier.

    Parameters:
        box_lengths (tuple): Lengths of the simulation cell along x, y, z directions.
        max_multiplier (int): Maximum multiplier for k-point generation. 

    Returns:
        list: List of k-point grids as tuples.
    """
    a, b, c = box_lengths
    
    inv_lengths = [1.0 / a, 1.0 / b, 1.0 / c]
    norm = max(inv_lengths)
    scaled_lengths = [l / norm for l in inv_lengths]
    return [
        (max(1, round(n * scaled_lengths[0])),
            max(1, round(n * scaled_lengths[1])),
            max(1, round(n * scaled_lengths[2])))
        for n in range(1, max_multiplier + 1)
    ]


if __name__ == '__main__':

    input_args = parse_args()

    # Default input parameters for Quantum Espresso
    input_data = {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'etot_conv_thr': 1e-5,
        'forc_conv_thr': 1e-4,
        'ecutwfc': input_args.ecutwfc,
        'ecutrho': input_args.ecutrho,
        'input_dft': 'pbe',
        'occupations': 'smearing',
        'degauss': 0.01,
        'smearing': 'gaussian',
        'conv_thr': 1e-8,
        'mixing_beta': 0.35,
        'startingwfc': 'random',
        'outdir': './tmp',
    }

    input_data.update({'prefix': input_args.prefix})

    # Read structure from CIF file
    system = read(input_args.cif_file)
    cell_lengths = system.cell.cellpar()[:3]  # Extract box lengths
    # Numer of atoms in the system
    num_atoms = len(system)

    # Load pseudopotential info from JSON file
    with open(input_args.json_pp_file, 'r') as f:
        pp_dict = json.load(f)

    # Prepare Quantum Espresso command
    command = f'mpiexec -n {input_args.np} {input_args.qe_path}'
    profile = EspressoProfile(command=command, pseudo_dir=input_args.pseudo_dir)

    # Parse k-points range
    kpoint_grids = generate_kpoints(cell_lengths, input_args.max_k_multiplier)

    energies = []
    log_file = open(input_args.log_file, 'w')

    # Perform calculations for different k-point samplings
    for kpts in kpoint_grids:
        energy = scf_calculation(system, profile, pp_dict, input_data, kpts)/num_atoms
        energies.append((kpts, energy))
        log_file.write(f'K-points: {kpts}, Energy: {energy:.6f} eV/atom\n')

    log_file.close()

    # Remove temporary files
    os.system(f'rm -rf {input_data["outdir"]}')

    # Determine optimal k-point grid
    energy_differences = [abs(energies[i + 1][1] - energies[i][1]) for i in range(len(energies) - 1)]
    optimal_index = next((i for i, diff in enumerate(energy_differences) if diff < input_args.convergence_criteria), len(energies) - 1)
    optimal_kpts, optimal_energy = energies[optimal_index]

    # Plot convergence results
    kpts_labels = ["x".join(map(str, k)) for k, _ in energies]
    energy_values = [e for _, e in energies]

    plt.figure(figsize=(8, 6))
    plt.plot(kpts_labels, energy_values, 'o-', label='Total Energy (eV/atom)')

    plt.axvline(optimal_index, color='r', linestyle='--', label=f'Optimal: {optimal_kpts}')
    plt.annotate(f'Optimal: {optimal_kpts}', xy=(optimal_index, optimal_energy), 
                 xytext=(optimal_index + 0.5, optimal_energy + 0.1),
                 arrowprops=dict(facecolor='red', arrowstyle='->'), fontsize=10)

    plt.xlabel('K-point Sampling', fontsize=12)
    plt.ylabel('Total Energy (eV/atom)', fontsize=12)
    plt.title('K-point Convergence', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig('kpoint_convergence.png', dpi=300)
    print(f"k-points Convergence plot saved as 'kpoint_convergence.png'. Optimal k-point grid: {optimal_kpts}.")

