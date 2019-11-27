from __future__ import print_function,division
import subprocess
import sys
import os
import numpy as np
import re
from functions import *
from multiprocessing import Pool, cpu_count


# Run by: python "script name" "input_file_name"(without extensions) phase type (liquid (L), gas (G), solid (S)) "user initials"(in caps)
# Water should be called "h2o" and vacuum should be called "vacuum"


def calculate_IFT_tot_and_coverage(input_file_name, phase_types, user, print_statements = True, debug = False, 
                                    multiprocess = True, delete_files = True, save_output_file = True, forced_convergence = True):
    """ Calculate the total interfacial tension of the two input phases and 
        the surface coverage between the phases.
    Args: 
        input_file_name: The name of the input file that the code should run either without extension, 
                         with extension or a path, as a string 
        phase_types: Input phase types of IFT calculation one char for each phase
                     as a string (L for liquid, G for gas, S for solid)
        user: The user initials in caps as a string
        print_statements: Print information from each iteration, boolean, default = True
        debug: Print additional information from each COSMOtherm calculation, boolean, default = False
        multiprocess: Run COSMOtherm simultaneously in the while loop, boolean, default = True
        delete_files: Delete the intermediate files created during the calculation, boolean, default = True
        save_output_file: Save the direct output of the calculation, boolean, default = True
        forced_convergence: Force the iterative process to end prematurely, boolean, defalut = False
        
    Return:
        coverage: The surface coverage between the two input phases as a numpy array
        IFT_tot: The total interfacial tension of the system as a float
    """
    # Correct to your own path to COSMOtherm
    COSMOtherm_path = get_user_and_path(user)

    # Initial values

    start_ift = 20.
    IFT_write_length = 5
    scale_organic = 1  # /0.91/0.8
    R = 8.314*1e-3  # The gas constant in kJ/mol/K
    unit_converter = 1.66  # Converts to mN/m
    max_iterations = 2  # force converges the while loop after max_iterations
    # Coverage
    max_CF = 2
    coverage_dampning = 0.5
    # IFT
    IFT_max_diff = 20
    IFT_dampning = 0.25
    # Convergence
    convergence_criteria = 3
    convergence_threshold = 1e-3
    # Solids
    max_depth = 2.0

    # Output precision
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format}, suppress = True)
    float_precision = 6

    # Change input name
    input_file_name, output_path = change_input_name(input_file_name)

    # Check unit=si in input file
    liq_ex = check_units_get_liq_ex(input_file_name)
    
    # Check if the water parameterization matches the input file
    scale_water, parameter = check_parameterization(input_file_name)

    # Read Number of compounds and Temperature from initial .inp file   
    N_compounds, T = get_N_compounds_and_T(input_file_name)
    
    if debug:
        print("N_compounds:", N_compounds, "Temperature:", T, "[K]")
    
    # Check phase types
    phase_types = check_phase_types(phase_types, 2)

    # Get the composition of the two phases from the .tab file for LL after LLE or from the .inp file for everything els    
    if phase_types == "LL":
        subprocess.call([COSMOtherm_path, input_file_name+".inp"])
        compound_list, phase1, phase2 = get_comp_and_phases_for_LL(input_file_name, N_compounds)
    else:
        compound_list, phase1, phase2 = get_comp_and_phases(input_file_name, N_compounds)
        
    phase1 = phase1/np.sum(phase1)
    phase2 = phase2/np.sum(phase2)
    
    liquid_index, solid_index = get_liquid_index(phase1, phase2, phase_types)
    
    # If there is a 0 in phase1, convert it to 10^-16
    if 0 in phase1[liquid_index] and phase_types[0] == "L":
        for i in np.where(phase1[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a concentration in phase 1, which was 0.0")
            phase1[i] = 1e-16
    # If there is a 0 in phase2, convert it to 10^-16
    if 0 in phase2[liquid_index] and phase_types[1] == "L":
        for i in np.where(phase2[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a concentration in phase 2, which was 0.0")
            phase2[i] = 1e-16
   
    
    if print_statements:
        print_compound_list = "[ {}".format(compound_list[0])
        for i in compound_list[1:]:
            print_compound_list += "  " + i
        print_compound_list += "]"
        print("Parameterization: {0} \nCompounds: {1} \nPhase 1:   {2} {3} \nPhase 2:   {4} {5}".format(parameter, print_compound_list, phase1, phase_types[0], phase2, phase_types[1]))    

    # Create flatsurfAB file
    write_flatsurf_file(input_file_name, "flatsurfAB", phase1, phase2, T, start_ift, IFT_write_length, phase_types, max_depth)

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
    if phase_types == "LL":
        coverage = np.sqrt(calculate_coverage(phase1, GtotAB, R, T, liquid_index) * calculate_coverage(phase2, GtotBA, R, T, liquid_index))
    elif phase_types == "LS":
        coverage = calculate_coverage(phase1, GtotAB, R, T, liquid_index)
    elif phase_types == "SL":
        coverage = calculate_coverage(phase2, GtotBA, R, T, liquid_index)
    else:
        print("Coverage calculation is wrong")
        quit()
    
    # Normalize coverage, using equation 2
    coverage /= np.sum(coverage)
    # If there is a 0 in the coverage, convert it to 10^-16
    if 0 in coverage[liquid_index]:
        for i in np.where(coverage[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a value in coverage, which was 0.0")
            coverage[i] = 1e-16
    
    IFT_A_value = start_ift
    IFT_B_value = start_ift
    IFT_tot = start_ift

    if debug:
        print("Coverage:", coverage, "IFT_A:", IFT_A_value, "IFT_B", IFT_B_value)

    phase_types = phase_types[0]+"C"+phase_types[1]  # Add C (coverage) as the middle phase
    N_cpu = cpu_count()        
    if N_cpu > 2:
        N_cpu = 2
    iterations = 0
    convergence_flag = 0
    if save_output_file:
        open(output_path + "output.txt", "w").close()
    while convergence_flag < convergence_criteria:
        iterations += 1
        
        if iterations == max_iterations and forced_convergence:
            print("The script ended before convergence!\nPhase 1:  {} \nCoverage: {} \nPhase 2:  {} \nTotal IFT: {}".format(phase1, coverage, phase2, IFT_tot))
            break
        
        # Create flatsurf files
        write_flatsurf_file(input_file_name, "flatsurfAS", phase1, coverage, T, IFT_A_value, IFT_write_length, phase_types[:2], max_depth)
        write_flatsurf_file(input_file_name, "flatsurfBS", phase2, coverage, T, IFT_B_value, IFT_write_length, phase_types[1:], max_depth)

        # Run both COSMOtherm instances simultaneously 
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
        
        #print(GtotAS, GtotBS)
        # Calculate coverages
        if phase_types == "LCL":
            coverage_A = calculate_coverage(phase1, GtotAS, R, T, liquid_index)
            coverage_B = calculate_coverage(phase2, GtotBS, R, T, liquid_index)
            # Calculate coverage factor (CF) and replace the value if it is too high or too low
            CF = np.power((coverage_A*coverage_B/coverage**2), coverage_dampning)
            CF[CF>max_CF] = max_CF
            CF[CF<1/max_CF] = 1/max_CF
            # Calculate new coverage
            coverage = coverage*CF
            coverage /= np.sum(coverage)
        elif phase_types == "LCS":
            coverage_A = calculate_coverage(phase1, GtotAS, R, T, liquid_index)
            # Calculate coverage factor (CF) and replace the value if it is too high or too low
            CF = np.power((coverage_A[liquid_index]/coverage[liquid_index]), coverage_dampning)
            CF[CF>max_CF] = max_CF
            CF[CF<1/max_CF] = 1/max_CF
            # Calculate new coverage
            coverage[liquid_index] = coverage[liquid_index]*CF
            coverage /= np.sum(coverage)
        elif phase_types == "SCL":
            coverage_B = calculate_coverage(phase2, GtotBS, R, T, liquid_index)
            # Calculate coverage factor (CF) and replace the value if it is too high or too low
            CF = np.power((coverage_B[liquid_index]/coverage[liquid_index]), coverage_dampning)
            CF[CF>max_CF] = max_CF
            CF[CF<1/max_CF] = 1/max_CF
            # Calculate new coverage
            coverage[liquid_index] = coverage[liquid_index]*CF
            coverage /= np.sum(coverage)
        else:
            print("Something went wrong in the surface coverage calculation.")
            quit()

        # Calculate IFT between phase and surface, using equation 3 for each direction
        IFT_A = calculate_IFT(phase1, GtotAS, GtotSA, AreaAS, AreaSA, coverage, R, T, 
                              unit_converter, phase_types[:2], liquid_index)
        IFT_B = calculate_IFT(phase2, GtotBS, GtotSB, AreaBS, AreaSB, coverage, R, T, 
                              unit_converter, phase_types[1:], liquid_index)
        
        # Damping IFT
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
            with open(output_path + "output.txt", "a") as file:
                file.write(", ".join(map(str,coverage))+", {}\n".format(IFT_tot))
    
    if delete_files:
        files = ["flatsurfAB.inp", "flatsurfAB.out", "flatsurfAB.tab", "flatsurfAS.inp", "flatsurfAS.out", "flatsurfAS.tab", "flatsurfBS.inp", "flatsurfBS.out", "flatsurfBS.tab"]
        for i in range(len(files)):
            if os.path.exists(files[i]):
                os.remove(files[i])
    np.set_printoptions(suppress = True)
    if iterations < max_iterations or not forced_convergence:
        print("The script has converged!\nPhase 1:  {} \nCoverage: {} \nPhase 2:  {} \nTotal IFT: {}".format(phase1, coverage, phase2, IFT_tot))
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