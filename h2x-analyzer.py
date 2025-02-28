from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse

def read_input_files(molecule):
    if molecule.upper() == "H2O":
        folder_path = "H2Ooutfiles"
    elif molecule.upper() == "H2S":
        folder_path = "H2Soutfiles"
    else:
        print("Invalid molecule. Please choose either 'H2O' or 'H2S'.")
        exit(1)
    
    # Initialize lists to store values
    bond_energies = []
    bond_angles = []
    bond_lengths = []

    # Define the input folder path
    input_folder_path = Path(folder_path)

    # Iterate through all files in the folder
    for filename in input_folder_path.iterdir():
        with open(filename, "r") as file:
            for line in file:
                parts = line.split()

                # Extract bond length and bond angle by looking at lines starting with H and has 5 elements
                if len(parts) == 5 and parts[0] == "H":
                    try:
                        bond_lengths.append(float(parts[2]))  
                        bond_angles.append(float(parts[4]))   
                    except ValueError:
                        pass  

                # Extract SCF energy
                if "SCF Done:" in line:
                    try:
                        bond_energies.append(float(parts[4])) 
                    except ValueError:
                        pass  
                    break  # Stop internal looping once SCF energy is found


    # Convert to numpy arrays for better performance
    bond_lengths = np.array(bond_lengths)
    bond_angles = np.array(bond_angles)
    bond_energies = np.array(bond_energies)
    
    bond_angles_rad = np.radians(bond_angles)
    
    return bond_lengths, bond_angles_rad, bond_energies

def compute_vibfreq(bond_lengths, bond_angles_rad, bond_energies, threshold):
    # constant conversion factors
    hartree_to_J = 4.35974417e-18  
    amu_to_kg = 1.66053886e-27   
    # Speed of light in cm/s   
    c_cm_s = 2.99792458e10        
    angstrom_to_meter = 1.0e-10 

    # Find minimum energies, angle and length
    min_energy_index = np.argmin(bond_energies)
    eqm_energy = bond_energies[min_energy_index]
    eqm_bond_length = bond_lengths[min_energy_index]
    eqm_bond_angle = bond_angles_rad[min_energy_index]
    
    # Mask to set threshold of bond energies to fit subsection of array for more accurate results (harmonic assumption works best near eqm?)
    # This returns True or False values for each value in bond energies
    mask = np.abs(bond_energies - eqm_energy) <= threshold

    # Applying the mask to the existing arrays will filter the array elements based on the position i.e. True is included, False is excluded
    sub_bond_energies = bond_energies[mask]
    sub_bond_lengths = bond_lengths[mask]
    sub_bond_angles  = bond_angles_rad[mask]

    # Define function to fit
    def fit_function(data, E0, kr, ktheta):
        r, theta = data
        return E0 + 0.5 * kr * (r - eqm_bond_length) ** 2 + 0.5 * ktheta * (theta - eqm_bond_angle) ** 2

    # Use curve_fit to fit the data
    popt, pcov = curve_fit(fit_function, (sub_bond_lengths, sub_bond_angles), sub_bond_energies)

    E0_fit = popt[0]
    kr_fit = popt[1]
    ktheta_fit = popt[2]

    # print(f"E0 = {E0_fit}, kr = {kr_fit}, ktheta = {ktheta_fit}")

    # Hartree/angstrom^2 to J/m^2
    kr_SI = kr_fit * hartree_to_J / (angstrom_to_meter ** 2)

    # Convert Hartree to J
    ktheta_SI = ktheta_fit * hartree_to_J 

    # Convert reduced masses to kg
    mu1 = 2 * amu_to_kg
    mu2 = 0.5 * amu_to_kg

    v1_hz = (1 / (2 * np.pi)) * np.sqrt(kr_SI / mu1)
    v2_hz = (1 / (2 * np.pi)) * np.sqrt(ktheta_SI / ((eqm_bond_length * angstrom_to_meter) ** 2 * mu2))

    # Hz to cm^-1
    v1_cm1 = v1_hz / c_cm_s
    v2_cm1 = v2_hz / c_cm_s

    print(f"Equilibrium Energy: {eqm_energy:.3f} Hartrees")
    print(f"Equilibrium Bond Length: {eqm_bond_length:.3f} Å")
    print(f"Equilibrium Bond Angle: {np.degrees(eqm_bond_angle):.1f}°")
    print(f"Symmetric Stretch Frequency: {v1_cm1:.0f} cm\u207B\u00B9")
    print(f"Bending Mode Frequency: {v2_cm1:.0f} cm\u207B\u00B9")
    
    return sub_bond_lengths, sub_bond_angles, sub_bond_energies


def create_surface_plot(bond_lengths, bond_angles_rad, bond_energies, molecule_name, sub_bond_lengths, sub_bond_angles, sub_bond_energies):
    bond_angles = np.degrees(bond_angles_rad) 
    
    # Plotting all the values of h2s will result in a very steep curve which will make the energy minimum difficult to see
    # Excluding some part of the steep drop in energy makes the minimum more obvious
    if molecule_name.lower() == "h2s":
        mask = bond_lengths > 0.9
        
        bond_lengths = bond_lengths[mask]
        bond_angles  = bond_angles[mask]
        bond_energies = bond_energies[mask]
        
    # Create the 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Create the surface plot for raw data values
    surf = ax.plot_trisurf(bond_lengths, bond_angles, bond_energies, cmap='viridis', linewidth=0.1, antialiased=True, edgecolor='gray')

    # Add a color bar with a label
    cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
    cbar.set_label('Z values', rotation=270, labelpad=15)

    # Set labels and title
    ax.set_xlabel('r/Å')
    ax.set_ylabel('θ/degrees')
    ax.set_zlabel('Energy/Hartrees')
    ax.set_title(f'Energy Surface Plot for {molecule_name.upper()}')
    
    # Adjust the viewing angle
    ax.view_init(elev=35, azim=-70)
    
    # Overlapping plot    
    fig2 = plt.figure(figsize=(10,8))
    
    # Number of rows in the subplot grid, number of columns in subplot grid, index of subplot
    ax2 = fig2.add_subplot(111, projection='3d')
    
    # Create the surface plot for raw data values
    sub_bond_angles_deg = np.degrees(sub_bond_angles)
    surf2 = ax2.plot_trisurf(bond_lengths, bond_angles, bond_energies, cmap='viridis', linewidth=0.1, antialiased=True, edgecolor='gray', alpha=0.8)
    surf3 = ax2.plot_trisurf(sub_bond_lengths, sub_bond_angles_deg, sub_bond_energies, color="magenta", edgecolor='gray', linewidth=0.1, antialiased=True, alpha=1.0)

    # Add a color bar with a label
    cbar = fig2.colorbar(surf2, shrink=0.5, aspect=5)
    cbar.set_label('Z values', rotation=270, labelpad=15)

    # Set labels and title
    ax2.set_xlabel('r/Å')
    ax2.set_ylabel('θ/degrees')
    ax2.set_zlabel('Energy/Hartrees')
    ax2.set_title(f'Overlapping energy surface plot for values used to determine vibrational frequencies of {molecule_name.upper()}')
    

    # Adjust the viewing angle
    ax2.view_init(elev=35, azim=-70)

    # Show the plot
    plt.show(block=False)
    plt.show()
    
    
    

def analyze_h2x(molecule, energy_threshold):
    data = read_input_files(molecule)
    sub_data = compute_vibfreq(*data, threshold=energy_threshold)
    create_surface_plot(data[0], data[1], data[2], molecule, sub_data[0], sub_data[1], sub_data[2])
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and analyze H2X data")
    
    parser.add_argument("-molecule", type=str, help="Specify molecule to compute (H2O or H2S)")
    parser.add_argument("-threshold", type=float, default=0.015, help="Energy threshold for filtering data (default: 0.02)")
    

    args = parser.parse_args()

    # If molecule is not provided as an argument, ask the user interactively
    if not args.molecule:
        molecule = input("Enter molecule (H2O or H2S): ").strip().upper()
        # Check validity
        if molecule not in ("H2O", "H2S"):
            print("Invalid molecule. Please choose either 'H2O' or 'H2S'.")
            exit(1)
            
        print(molecule)
        print("Loading...\n")
    else:
        molecule = args.molecule
        print(molecule.upper())
        print("Loading...\n")

    analyze_h2x(molecule, args.threshold)
