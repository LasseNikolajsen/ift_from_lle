from __future__ import print_function,division
import subprocess
import sys
import os
import numpy as np
import re
from functions import *

def calculate_IFT_tot(input_file_name, user):
    
    print_statements = True
    save_output_file = True
    forced_convergence = False
    max_iterations = 3
    
    phase_types = "CS"
    
    start_ift = 10.  # Start_ift * 2 is the start position in the iterative process
    IFT_write_length = 5  # Decimals when writing IFT in the flatsurf files, more than 5 triggers an error
    scale_organic = 1  # /0.91/0.8
    R = 8.314*1e-3  # The gas constant in kJ/mol/K
    unit_converter = 1.66  # Converts to mN/m
    # IFT
    IFT_max_diff = 40  # Limits the step length of IFT per iteration
    IFT_damping = 0.25  # Prevents oscillations
    # Convergence
    convergence_criteria = 3  # Number of iterations with an IFT difference under convergence_threshold
    convergence_threshold = 1e-3
    inf_loop_precision = 3  # The precision for the infinite loop check, high number equals less likely to occur
    # Solids
    max_depth = 2.0  # max depth in Angstrom for the flatsurf calculations including a solid phase 

    # Output precision
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format}, suppress = True)
    float_precision = 4
    
    COSMOtherm_path = get_user_and_path(user)
    
    input_file_name, output_path = change_input_name(input_file_name)
    
    liq_ex = check_units_get_liq_ex(input_file_name)
    
    scale_water, parameter = check_parameterization(input_file_name)
    
    N_compounds, T = get_N_compounds_and_T(input_file_name)

    compound_list, phases = get_comp_and_phases(input_file_name, N_compounds)
    coverage = phases[0]
    phase2 = phases[1]
    
    liquid_index, solid_index = get_liquid_index(coverage, phase2, phase_types)
    
    # If there is a 0 in the phase, convert it to 10^-16
    if 0 in coverage[liquid_index] and phase_types[0] == "C":
        for i in np.where(coverage[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a concentration in phase 1, which was 0.0")
            coverage[i] = 1e-16
    if 0 in phase2[liquid_index] and phase_types[1] == "C":
        for i in np.where(phase2[liquid_index]==0)[0]:
            print("Warning: Added 1e-16 to a concentration in phase 2, which was 0.0")
            phase2[i] = 1e-16
   
    # Normalize the phases
    coverage = coverage/np.sum(coverage)
    phase2 = phase2/np.sum(phase2)    
    
    if print_statements:
        print_compound_list = "[ {}".format(compound_list[0])
        for i in compound_list[1:]:
            print_compound_list += "  " + i
        print_compound_list += "]"
        print("Parameterization: {0} \nCompounds: {1} \nCoverage:  {2} {3} \nPhase 2:   {4} {5}".format(parameter, print_compound_list, coverage, phase_types[0], phase2, phase_types[1]))    

    IFT_tot = start_ift
    IFT_B_value = start_ift
    iterations = 0
    convergence_flag = 0
    inf_loop_counter = 0
    IFT_tot_list = []
    if save_output_file:
        open(output_path + "output.txt", "w").close()
    while convergence_flag < convergence_criteria:
        iterations += 1
    
        # Check for forced convergence
        if iterations == max_iterations+1 and forced_convergence:
            print("The script ended before convergence!\nCoverage: {} \nPhase 2:  {} \nTotal IFT: {}".format(coverage, phase2, IFT_tot))
            break
    
        write_flatsurf_file(input_file_name, "flatsurf_C_2", coverage, phase2, T, IFT_B_value, IFT_write_length, phase_types, max_depth)
        subprocess.call([COSMOtherm_path, "flatsurf_C_2.inp"]) 

        Gtot_C_2, Gtot_2_C, Area_C_2, Area_2_C = get_Gtot_and_Area("flatsurf_C_2", N_compounds)
        Area_C_2, Area_2_C = scale_area(compound_list, Area_C_2, Area_2_C, N_compounds, scale_water, scale_organic)
        
        print("Gtot C->2:", Gtot_C_2, "Gtot 2->C:", Gtot_2_C)
        
        IFT_B = calculate_IFT(phase2, Gtot_2_C, Gtot_C_2, Area_2_C, Area_C_2, coverage, R, T, unit_converter, phase_types, liquid_index)
        
        IFT_B_value = calculate_IFT_damping(IFT_B, IFT_B_value, IFT_max_diff, IFT_damping)
        
        IFT_tot_old = IFT_tot
        IFT_tot = IFT_B_value * 0.5
        
        # Check for out of bounds total IFT to prevent COSMOtherm error
        if IFT_tot < -95.0 or IFT_tot > 120.0:
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

        if save_output_file:
            with open(output_path + "output.txt", "a") as file:
                file.write(", ".join(map(str,coverage))+", {}\n".format(IFT_tot))
    
    
    print("The script has converged!\nCoverage: {} \nPhase 2:  {} \nTotal IFT: {}".format(coverage, phase2, IFT_tot))
    
    return coverage, IFT_tot
    
    
if __name__ == "__main__":
    # Inputs from the terminal
    try:
        input_file_name = sys.argv[1]
        user = sys.argv[2]
    except:
        print("Incorrect inputs, run by: python \"script name\" \"input_file_name\"(without extension) \"user_name\"")
        quit()
    
    coverage_final, IFT_final = calculate_IFT_tot(input_file_name, user)
    
    
    
    
    
    