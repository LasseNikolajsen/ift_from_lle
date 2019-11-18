from __future__ import print_function,division
import subprocess
import sys
import os
import numpy as np
import re
from functions_for_ift import *
from multiprocessing import Pool, cpu_count


# Run by: python "scriptname" "input_file_name"(without extensions) mix "user initials"(in caps)
# Water should be called "h2o" and vacuum should be called "vacuum"


def calculate_IFT_tot_and_coverage(input_file_name, phase_types, user, print_statements = True, debug = False, 
                                    multiprocess = False, delete_files = True, save_output_file = False):
    """ Calculate the total interfacial tension of the two input phases and 
        the surface coverage inbetween the phases.
    Args: 
        input_file_name: The name of the input file that the code should run either without extension, 
                         with extension or a path, as a string 
        phase_types: Input phase types of IFT calculation one char for each phase
                     as a string (L for liquid, G for gas, S for solid)
        user: The user initials in caps as a string
        print_statements: Print information from each iteration, boolean, optional, default = True
        debug: Print additional information from each COSMOtherm calculation, boolean, optional, default = False
        multiprocess: Run COSMOtherm simultaniously in the while loop, boolean, optional, default = True
        
    Return:
        coverage: The surface coverage between the two input phases as a numpy array
        IFT_tot: The total interfacial tension of the system as a float
    """
    # Correct to your own path to COSMOtherm
    if user == "LVN":
        COSMOtherm_path = "C:\Program Files\COSMOlogic\COSMOthermX19\COSMOtherm\BIN-WINDOWS\cosmotherm.exe"  # LVN
    elif user == "MPA":
        COSMOtherm_path = "/Applications/COSMOlogic/COSMOthermX18/COSMOtherm/BIN-LINUX/cosmotherm"  # MPA
    else:
        print("User not recognized, add new user and COSMOtherm path in the code.")
        quit()

    # Initial values

    start_ift = 20.
    IFT_write_length = 5
    scale_organic = 1  # /0.91/0.8
    R = 8.314*1e-3  # The gas constant in kJ/mol/K
    unit_converter = 1.66  # Converts to mN/m
    # Coverage
    max_CF = 2
    coverage_dampning = 0.5
    # IFT
    IFT_max_diff = 20
    IFT_dampning = 0.25
    # Convergence
    convergence_criteria = 3
    convergence_threshold = 1e-3

    # Output precision
    np.set_printoptions(precision = 6, suppress = True)
    float_precision = 10

    # Change input name
    input_file_name = change_input_name(input_file_name)

    # Check unit=si in input file
    check_units(input_file_name)
    
    # Check phase types
    check_phase_types(phase_types)
    
    # Check if this file water parameterization matches the input file
    scale_water, parameter = check_parameterization(input_file_name)

    # Run COSMOtherm
    subprocess.call([COSMOtherm_path, input_file_name+".inp"])

    # Read Number of compounds and Temperature from initial .inp file   
    N_compounds, T = get_N_compounds_and_T(input_file_name)
    
    if debug:    
        print("N_compounds:", N_compounds, "Temperature:", T, "[K]")

    # Get the composition of the two phases from the .tab file      
    compound_list, phase1, phase2 = get_comp_and_phases(input_file_name, N_compounds)
    phase1 = phase1/np.sum(phase1)
    phase2 = phase2/np.sum(phase2)
    # If there is a 0 in the coverage, convert it to 10^-16
    if 0 in phase1:
        for i in np.where(phase1==0)[0]:
            phase1[i] = 1e-16
    # If there is a 0 in the coverage, convert it to 10^-16
    if 0 in phase2:
        for i in np.where(phase2==0)[0]:
            phase2[i] = 1e-16
    
    
    if print_statements:
        print("Parameterization: {3} \nCompounds: {0} \nphase1: {1} \nphase2: {2}".format(compound_list, phase1, phase2, parameter))    

    # Create flatsurfAB file
    write_flatsurf_file(input_file_name, "flatsurfAB", phase1, phase2, T, start_ift, IFT_write_length)

    # Run COSMOtherm
    subprocess.call([COSMOtherm_path, "flatsurfAB.inp"])    

    # Extract Gtot and across,mean for each direction in the flatsurf file
    GtotAB, GtotBA, AreaAB, AreaBA = get_Gtot_and_Area("flatsurfAB", N_compounds)

    if debug:
        print("Gtot, AB:", GtotAB, "BA:", GtotBA)
        print("Area, AB:", AreaAB, "BA:", AreaBA)

    # Scale the calculated areas
    AreaAB, AreaBA = scale_area(compound_list, AreaAB, AreaBA, N_compounds, scale_water, scale_organic)

    # Calculate the coverage in the interface between A and B, using equation 1
    coverage = np.sqrt(calculate_coverage(phase1, GtotAB, R, T) * calculate_coverage(phase2, GtotBA, R, T))
    # Normalize coverage, using equation 2
    coverage /= np.sum(coverage)
    # If there is a 0 in the coverage, convert it to 10^-16
    if 0 in coverage:
        for i in np.where(coverage==0)[0]:
            coverage[i] = 1e-16

    IFT_A_value = start_ift
    IFT_B_value = start_ift
    IFT_tot = start_ift

    if debug:
        print("Coverage:", coverage, "IFT_A:", IFT_A_value, "IFT_B", IFT_B_value)

    N_cpu = cpu_count()        
    if N_cpu > 2:
        N_cpu = 2
    iterations = 0
    convergence_flag = 0
    if save_output_file:
        open("output.txt", "w").close()
    while convergence_flag < convergence_criteria:
        iterations += 1
            
        # Create flatsurf files
        write_flatsurf_file(input_file_name, "flatsurfAS", phase1, coverage, T, IFT_A_value, IFT_write_length)
        write_flatsurf_file(input_file_name, "flatsurfBS", phase2, coverage, T, IFT_B_value, IFT_write_length)
        
        # Run both COSMOtherm instances simultaniously 
        if multiprocess:
            pool = Pool(processes=N_cpu)
            pool.map(work, [[COSMOtherm_path, os.path.abspath("flatsurfAS.inp")], 
                            [COSMOtherm_path, os.path.abspath("flatsurfBS.inp")]])
            pool.close()
            pool.join()

        else:
            subprocess.call([COSMOtherm_path, "flatsurfAS.inp"]) 
            subprocess.call([COSMOtherm_path, "flatsurfBS.inp"])  
        
        # Extract Gtot and Area from the .tab file
        GtotAS, GtotSA, AreaAS, AreaSA = get_Gtot_and_Area("flatsurfAS", N_compounds)
        GtotBS, GtotSB, AreaBS, AreaSB = get_Gtot_and_Area("flatsurfBS", N_compounds)

        # Scale areas
        AreaAS, AreaSA = scale_area(compound_list, AreaAS, AreaSA, N_compounds, scale_water, scale_organic)
        AreaBS, AreaSB = scale_area(compound_list, AreaBS, AreaSB, N_compounds, scale_water, scale_organic)
        
        # Calculate coverages
        coverage_A = calculate_coverage(phase1, GtotAS, R, T)
        coverage_B = calculate_coverage(phase2, GtotBS, R, T)

        # Calculate coverage factor (CF) and replace the value if it is too high or too low
        CF = np.power((coverage_A*coverage_B/coverage**2), coverage_dampning)
        CF[CF>max_CF] = max_CF
        CF[CF<1/max_CF] = 1/max_CF
        # Calculate new coverage
        coverage = coverage*CF
        coverage /= np.sum(coverage)

        # Calculate IFT between phase and surface, using equation 3 for each direction
        IFT_A = calculate_IFT(phase1, GtotAS, AreaAS, coverage, R, T, unit_converter)
        IFT_B = calculate_IFT(phase2, GtotBS, AreaBS, coverage, R, T, unit_converter)
        
        # Dampning IFT
        IFT_A_value = calculate_IFT_dampning(IFT_A, IFT_A_value, IFT_max_diff, IFT_dampning)
        IFT_B_value = calculate_IFT_dampning(IFT_B, IFT_B_value, IFT_max_diff, IFT_dampning)
        
        # Calculate total system IFT
        IFT_tot_old = IFT_tot
        IFT_tot = IFT_A_value + IFT_B_value
        
        # Check convergence criteria
        if abs(IFT_tot_old-IFT_tot) < convergence_threshold:
            convergence_flag += 1
        else: 
            convergence_flag = 0
        if print_statements:
            print("Iterations: {0:>2} Coverage: {1} IFT_total: {2:.{3}}".format(str(iterations), coverage, IFT_tot, float_precision))
        
        if debug:
            print("Gtot, AS:", GtotAS, "SA:", GtotSA)
            print("Area, AS:", AreaAS, "SA:", AreaSA)
            print("Coverage_A:", coverage_A, "IFT_AS:", "IFT_A:", IFT_A, "IFT_A_value", IFT_A_value)
            print("Gtot, BS:", GtotBS, "SB:", GtotSB)
            print("Area, BS:", AreaBS, "SB:", AreaSB)
            print("Coverage_B:", coverage_B, "IFT_B:", IFT_B, "IFT_B_value", IFT_B_value)
            print("IFT_tot", IFT_tot)
            print("Coverage", coverage)
            print("IFT difference", IFT_tot-IFT_tot_old)
        if save_output_file:
            with open("output.txt", "a") as file:
                file.write(", ".join(map(str,coverage))+", {}\n".format(IFT_tot))
    
    if delete_files:
        files = ["flatsurfAB.inp", "flatsurfAB.out", "flatsurfAB.tab", "flatsurfAS.inp", "flatsurfAS.out", "flatsurfAS.tab", "flatsurfBS.inp", "flatsurfBS.out", "flatsurfBS.tab"]
        for i in range(len(files)):
            if os.path.exists(files[i]):
                os.remove(files[i])

    print("coverage:", coverage, "IFT_tot:", IFT_tot)        
    return coverage, IFT_tot


    
    
if __name__ == "__main__":
    # Inputs from the terminal
    try:
        input_file_name = sys.argv[1]
        mix = sys.argv[2]
        user = sys.argv[3]
    except:
        print("Incorrect inputs, run by: python \"script name\" \"input_file_name\"(without extension) mix \"user initials\"(in caps)")
        quit()
    
    coverage_final, IFT_final = calculate_IFT_tot_and_coverage(input_file_name, mix, user)