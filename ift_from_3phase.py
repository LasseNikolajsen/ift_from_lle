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
                                    multiprocess = True, delete_files = True, save_output_file = True, forced_convergence = False, max_iterations = 3):
    """ Calculate the total interfacial tension of the two input phases and 
        the surface coverage between the phases.
    Args: 
        input_file_name: The name of the input file that the code should run either without extension, with extension or a path, as a string 
        phase_types: Input phase types of IFT calculation one char for each phase as a string (L for liquid, G for gas, S for solid)
        user: The user initials in caps as a string
        print_statements: Print information from each iteration, boolean, default = True
        debug: Print additional information from each COSMOtherm calculation, boolean, default = False
        multiprocess: Run COSMOtherm simultaneously in the while loop, boolean, default = True
        delete_files: Delete the intermediate files created during the calculation, boolean, default = True
        save_output_file: Save the direct output of the calculation, boolean, default = True
        forced_convergence: Force the iterative process to end prematurely, boolean, default = False
        max_iterations: The maximum amount of iterations for the iterative process if forced convergence is set to True, max_iterations should be an integer
        
    Return:
        coverage: The surface coverage between the two input phases as a numpy array
        IFT_tot: The total interfacial tension of the system as a float
    """
    # Add your own path to COSMOtherm and user name in the Users.txt file
    COSMOtherm_path = get_user_and_path(user)

    # Initial values
    start_ift = 20.  # Start_ift * 2 is the start position in the iterative process
    IFT_write_length = 5  # Decimals when writing IFT in the flatsurf files, more than 5 triggers an error
    scale_organic = 1.  # /0.91/0.8
    R = 8.314*1e-3  # The gas constant in kJ/mol/K
    unit_converter = 1.66  # Converts to mN/m
    # Coverage
    max_CF = 2.  # Limits the step length of coverage per iteration to twice the current concentration
    coverage_damping = 0.5  # Prevents oscillations
    # IFT
    IFT_max_diff = 40.  # Limits the step length of IFT per iteration
    IFT_damping = 0.25  # Prevents oscillations
    # Convergence
    convergence_criteria = 3  # Number of iterations with an IFT difference under convergence_threshold
    convergence_threshold = 1e-3
    inf_loop_precision = 3  # The precision for the infinite loop check, high number equals less likely to occur
    # Solids
    max_depth = 1.0  # max depth in Angstrom for the flatsurf calculations including a solid phase 

    # Output precision
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format}, suppress = True)
    float_precision = 4

    # Change input name
    input_file_name, output_path = change_input_name(input_file_name)

    # Check unit=si in input file
    liq_ex = check_units_get_liq_ex(input_file_name)
    
    # Check if the parameterization matches the input file
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
        compound_list, phases = get_comp_and_phases(input_file_name, N_compounds)
        phase1 = phases[0]
        phase2 = phases[1]
        
    # Get the indices of the liquid and solid compounds
    liquid_index, solid_index = get_liquid_index(phase1, phase2, phase_types)
    
    # If there is a 0 in the phase, convert it to 10^-16
    if 0 in phase1[liquid_index] and phase_types[0] == "L":
        for i in np.where(phase1[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a concentration in phase 1, which was 0.0")
            phase1[i] = 1e-16
    if 0 in phase2[liquid_index] and phase_types[1] == "L":
        for i in np.where(phase2[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a concentration in phase 2, which was 0.0")
            phase2[i] = 1e-16
   
    # Normalize the phases
    phase1 = phase1/np.sum(phase1)
    phase2 = phase2/np.sum(phase2)
   
    # Print initial values
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

    # Calculate the coverage in the interface between A and B, using equation 1 for LL and a reduced equation for LS and SL
    if phase_types == "LL":
        coverage = np.sqrt(calculate_coverage(phase1, GtotAB, R, T, liquid_index) * calculate_coverage(phase2, GtotBA, R, T, liquid_index))
    elif phase_types == "LS":
        coverage = calculate_coverage(phase1, GtotAB, R, T, liquid_index)
    elif phase_types == "SL":
        coverage = calculate_coverage(phase2, GtotBA, R, T, liquid_index)
    
    # If there is a 0 in the coverage, convert it to 10^-16
    if 0 in coverage[liquid_index]:
        for i in np.where(coverage[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a value in coverage, which was 0.0")
            coverage[i] = 1e-16
    
    # Normalize coverage
    coverage /= np.sum(coverage)
    
    # Initiate values for iterative process
    IFT_A_value = start_ift
    IFT_B_value = start_ift
    IFT_tot = start_ift
    phase_types = phase_types[0]+"C"+phase_types[1]  # Add C (coverage) as the middle phase
    iterations = 0
    convergence_flag = 0
    inf_loop_counter = 0
    IFT_tot_list = []
    # Multiprocessing
    N_cpu = cpu_count()        
    if N_cpu > 2:
        N_cpu = 2
    # Open output file
    if save_output_file:
        open(output_path + "output.txt", "w").close()
    while convergence_flag < convergence_criteria:
        iterations += 1
        
        # Check for forced convergence
        if iterations == max_iterations+1 and forced_convergence:
            print("The script ended before convergence!\nPhase 1:  {} \nCoverage: {} \nPhase 2:  {} \nTotal IFT: {}".format(phase1, coverage, phase2, IFT_tot))
            break
        
        # Create flatsurf files for phase1/coverage and coverage/phase2
        write_flatsurf_file(input_file_name, "flatsurfAS", phase1, coverage, T, IFT_A_value, IFT_write_length, phase_types[:2], max_depth)
        write_flatsurf_file(input_file_name, "flatsurfSB", coverage, phase2, T, IFT_B_value, IFT_write_length, phase_types[1:], max_depth)

        if multiprocess:  # Run both COSMOtherm instances simultaneously
            pool = Pool(processes=N_cpu)
            pool.map(work, [[COSMOtherm_path, os.path.abspath("flatsurfAS.inp")], [COSMOtherm_path, os.path.abspath("flatsurfSB.inp")]])
            pool.close()
            pool.join()

        else:  # One at a time
            subprocess.call([COSMOtherm_path, "flatsurfAS.inp"]) 
            subprocess.call([COSMOtherm_path, "flatsurfSB.inp"])  
        
        # Extract Gtot and Area from the .tab files
        GtotAS, GtotSA, AreaAS, AreaSA = get_Gtot_and_Area("flatsurfAS", N_compounds)
        GtotSB, GtotBS, AreaSB, AreaBS = get_Gtot_and_Area("flatsurfSB", N_compounds)

        # Scale areas
        AreaAS, AreaSA = scale_area(compound_list, AreaAS, AreaSA, N_compounds, scale_water, scale_organic)
        AreaBS, AreaSB = scale_area(compound_list, AreaBS, AreaSB, N_compounds, scale_water, scale_organic)
        
        # Calculate coverages
        if phase_types == "LCL":
            coverage_A = calculate_coverage(phase1, GtotAS, R, T, liquid_index)
            coverage_B = calculate_coverage(phase2, GtotBS, R, T, liquid_index)
            coverage = calculate_CF(coverage, [coverage_A, coverage_B], coverage_damping, max_CF, liquid_index)
        elif phase_types == "LCS":
            coverage_A = calculate_coverage(phase1, GtotAS, R, T, liquid_index)
            coverage = calculate_CF(coverage, coverage_A, coverage_damping, max_CF, liquid_index)
        elif phase_types == "SCL":
            coverage_B = calculate_coverage(phase2, GtotBS, R, T, liquid_index)
            coverage = calculate_CF(coverage, coverage_B, coverage_damping, max_CF, liquid_index)

        # Calculate IFT between phase and surface
        IFT_A = calculate_IFT(phase1, GtotAS, GtotSA, AreaAS, AreaSA, coverage, R, T, unit_converter, phase_types[:2], liquid_index)
        IFT_B = calculate_IFT(phase2, GtotBS, GtotSB, AreaBS, AreaSB, coverage, R, T, unit_converter, phase_types[1:], liquid_index)
        
        # Damping IFT
        IFT_A_value = calculate_IFT_damping(IFT_A, IFT_A_value, IFT_max_diff, IFT_damping)
        IFT_B_value = calculate_IFT_damping(IFT_B, IFT_B_value, IFT_max_diff, IFT_damping)
        
        # Calculate total system IFT
        IFT_tot_old = IFT_tot
        if phase_types[-1] == "S":
            IFT_tot = IFT_A_value + IFT_B_value
        else:
            IFT_tot = IFT_A_value + IFT_B_value
            
        # Check for out of bounds total IFT to prevent COSMOtherm error
        if IFT_tot < -95.0:
            print("Warning: The IFT is out of bounds, the iterative process will end without convergence")
            break
        
        # Check convergence criteria
        if abs(IFT_tot_old-IFT_tot) < convergence_threshold:
            convergence_flag += 1
        else: 
            convergence_flag = 0
        
        # Print current iteration results
        if print_statements:
            print("Iterations: {0:>2} Coverage: {1} IFT_total: {2:>8.{3}f}".format(str(iterations), coverage, IFT_tot, float_precision))

        # Check for infinite loop
        if "{:.{}f}".format(IFT_tot, inf_loop_precision) in IFT_tot_list:
            inf_loop_counter += 1
            if inf_loop_counter > 3:
                IFT_damping_new = IFT_damping * 0.5
                print("Infinite loop detected, changed the IFT_damping from {} to {}".format(IFT_damping, IFT_damping_new))
                IFT_damping = IFT_damping_new
                inf_loop_counter = 0
        IFT_tot_list.append("{:.{}f}".format(IFT_tot, inf_loop_precision))
        
        if debug:
            print("Gtot, AS:", GtotAS, "SA:", GtotSA)
            print("Area, AS:", AreaAS, "SA:", AreaSA)
            if phase_types[0] == "S":
                print("IFT_A:", IFT_A, "IFT_A_value", IFT_A_value)
            else:
                print("Coverage_A:", coverage_A, "IFT_AS:", "IFT_A:", IFT_A, "IFT_A_value", IFT_A_value)
            print("Gtot, SB:", GtotSB, "BS:", GtotBS)
            print("Area, SB:", AreaSB, "BS:", AreaBS)
            if phase_types[2] == "S":
                print("IFT_B:", IFT_B, "IFT_B_value", IFT_B_value)
            else:
                print("Coverage_B:", coverage_B, "IFT_B:", IFT_B, "IFT_B_value", IFT_B_value)
            print("\n")
        if save_output_file:
            with open(output_path + "output.txt", "a") as file:
                file.write(", ".join(map(str,coverage))+", {}\n".format(IFT_tot))
    
    # Delete files used in the calculation
    if delete_files:
        files = ["flatsurfAB.inp", "flatsurfAB.out", "flatsurfAB.tab", "flatsurfAS.inp", "flatsurfAS.out", "flatsurfAS.tab", "flatsurfSB.inp", "flatsurfSB.out", "flatsurfSB.tab"]
        for i in range(len(files)):
            if os.path.exists(files[i]):
                os.remove(files[i])
    np.set_printoptions(suppress = True)
    
    # Print final result
    if iterations < max_iterations or not forced_convergence:
        print("The script has converged!\nPhase 1:  {} \nCoverage: {} \nPhase 2:  {} \nTotal IFT: {}".format(phase1, coverage, phase2, IFT_tot))
    
    return coverage, IFT_tot

if __name__ == "__main__":
    # Inputs from the terminal
    try:
        input_file_name = sys.argv[1]
        phase_types = sys.argv[2]
        user = sys.argv[3]
    except:
        print("Incorrect inputs, run by: python \"script name\" \"input_file_name\"(without extension) \"phase_types(L, S or G)\" \"user_name\"")
        quit()
    
    coverage_final, IFT_final = calculate_IFT_tot_and_coverage(input_file_name, phase_types, user)