from __future__ import print_function,division
import sys
from os import path, remove
import numpy as np
import pandas as pd
from ift_from_3phase import calculate_IFT_tot_and_coverage
from functions import change_input_name, get_comp_and_phases, get_N_compounds_and_T

def input_file_to_IFT(phase1, phase2, phase_types, types, input_file, output_path, user, N_comps, first_comp_line_index, N_lines_p_compound, 
                      phase1_compounds_index, phase2_compounds_index, phases):
    with open(input_file+".inp", "r") as file: 
        text = file.readlines()
        
        extra_lines = 0
        first_value = True
        with open(output_path+str(phase_types)+"_input.inp", "w") as WS_file:
            WS_file.writelines(text[:first_comp_line_index])
            for i in range(N_comps):
                if i in phase1_compounds_index+phase2_compounds_index:
                    for j in range(N_lines_p_compound[i]):
                        WS_file.write(text[first_comp_line_index+i+j+extra_lines])
                        
                extra_lines += N_lines_p_compound[i]-1
            
            last_line_modified = ""
            last_line = text[-1]
            last_line = last_line.split()
            last_line_modified += last_line[0] + " " + last_line[1][:-1]+"2 x1={"
            first_value = True
            for i in range(N_comps):
                if phase1[i] > 0.0 or phase2[i] > 0.0:
                    if first_value:
                        last_line_modified += str(phase1[i])
                        first_value = False
                    else:
                        last_line_modified += " " + str(phase1[i])
            last_line_modified += "} x2={"
            first_value = True
            for i in range(N_comps):
                if phase1[i] > 0.0 or phase2[i] > 0.0:
                    if first_value:
                        last_line_modified += str(phase2[i])
                        first_value = False
                    else:
                        last_line_modified += " " + str(phase2[i])
            last_line_modified += "}"
            for i in range(2+N_comps*len(phases), len(last_line)):
                last_line_modified += " " + last_line[i]
            WS_file.write(last_line_modified)

    IFT, coverage = calculate_IFT_tot_and_coverage(output_path+str(phase_types)+"_input.inp", types, user, save_output_file = False, forced_convergence = False, max_iterations = 2)
    return IFT, coverage

def main():
    
    # Use a 3 phase liquid extraction COSMOtherm input file and specify which phases are water, oil and sold.
    
    phase_types = "WOS"  # Water (O), Oil (O), Solid (S)
    
    WO_IFT = 0.0  # Water oil, if 0.0 run the calculation, else use specified value
    
    WS_IFT = 0.0  # Water solid, if 0.0 run the calculation, else use specified value
    
    OS_IFT = 0.0  # Oil solid, if 0.0 run the calculation, else use specified value
    
    if not(WS_IFT != 0.0 and OS_IFT != 0.0 and WO_IFT != 0.0):
        input_file = sys.argv[1]
        
        input, output = change_input_name(input_file)

        N_comps, T = get_N_compounds_and_T(input)
        
        comp_list, phases = get_comp_and_phases(input, N_comps)

        phase_types_input = ""
        for i in range(len(phase_types)):
            if phase_types[i] == "S":
                phase_types_input += "S"
            else:
                phase_types_input += "L"
        
        water_index = phase_types.index("W")
        oil_index = phase_types.index("O")
        solid_index = phase_types.index("S")

        water_phase = phases[water_index]
        water_compounds_index = []
        oil_phase = phases[oil_index]
        oil_compounds_index = []
        solid_phase = phases[solid_index]
        solid_compounds_index = []
        for i in range(len(comp_list)):
            if water_phase[i] != 0.0:
                water_compounds_index.append(i)
            if oil_phase[i] != 0.0:
                oil_compounds_index.append(i)
            if solid_phase[i] != 0.0:
                solid_compounds_index.append(i)
        
        
        N_lines_p_compound = []   
        count = 0
        first_value = True
        first_comp_line_index = 0
        
        
        with open(input+".inp", "r") as file: 
            text = file.readlines()
            
            
            for i in range(len(text)):
                if "VPfile" in text[i] or "liq_ex" in text[i]:
                    if first_comp_line_index == 0:
                        first_comp_line_index = i
                    if len(N_lines_p_compound) == 1 and first_value:
                        N_lines_p_compound[0] = (count)
                        first_value = False
                    else:
                        N_lines_p_compound.append(count)
                    count = 0
                count += 1

            if WO_IFT == 0.0:
                print("\nCalculating water/oil interface:\n")
                WO_coverage, WO_IFT = input_file_to_IFT(water_phase, oil_phase, phase_types[water_index]+phase_types[oil_index], "LL", input, output, "LVND", N_comps, 
                                                        first_comp_line_index, N_lines_p_compound, water_compounds_index, oil_compounds_index, phases)
            if WS_IFT == 0.0:
                print("\nCalculating water/solid interface:\n")
                WS_coverage, WS_IFT = input_file_to_IFT(water_phase, solid_phase, phase_types[water_index]+phase_types[solid_index], "LS", input, output, "LVND", N_comps, 
                                                        first_comp_line_index, N_lines_p_compound, water_compounds_index, solid_compounds_index, phases)
            if OS_IFT == 0.0:
                print("\nCalculating oil/solid interface:\n")
                OS_coverage, OS_IFT = input_file_to_IFT(oil_phase, solid_phase, phase_types[oil_index]+phase_types[solid_index], "LS", input, output, "LVND", N_comps,
                                                        first_comp_line_index, N_lines_p_compound, oil_compounds_index, solid_compounds_index, phases)

    # n-decane
    # IFT_W_S = 64.9449187143066
    # IFT_O_S = -12.500284722858611
    # IFT_O_W = 46.870301105580424
    
    # 1-octanol
    # IFT_W_S = -18.626208743461117
    # IFT_O_S = 10.935139924034504
    # IFT_O_W = 46.870301105580424

    # acetic acid
    # IFT_W_S = -14.774634533876746
    # IFT_O_S = 9.751316977037327
    # IFT_O_W = 46.870301105580424
    
    
    youngs_eq = (OS_IFT - WS_IFT) / WO_IFT
    if youngs_eq > 1:
        print("Error: Can not take arccos to a number ouside the range [-1,1]. Please check the calculated energies.")
        print("The calculated number is {}.".format(youngs_eq))
        quit()
    elif youngs_eq < -1:
        youngs_eq = -1
        
    contact_angle = np.arccos(youngs_eq) * (180/np.pi)
    print("\nContact angle [degrees]: {:.4}".format(contact_angle))


    # Create pandas data frame
    matrix = np.array((WO_IFT, WS_IFT, OS_IFT))
    df = pd.DataFrame(matrix, columns=["IFT"], index=["WO:", "WS:", "OS:"])

    # Print output in n_phase_output.txt
    with open(output+"IFT_contact_angle_output.txt", 'w') as file:
        sys.stdout = file
        print(df)
        line = "\nContact angle [degrees]: {}".format(contact_angle) + "\nPhase types: {}\n".format(phase_types)
        file.write(line)
        with open(input+".inp", "r") as input:
            text = input.read()
            file.write("\n\n\nInitial input:\n")
            file.write(text)
        if path.exists(output+"WO_input.inp"):
            with open(output+"WO_input.inp", "r") as input:
                text = input.read()
                file.write("\n\nWater/oil (LL) input:\n")
                file.write(text)
            remove(output+"WO_input.inp")
            remove(output+"WO_input.out")
            remove(output+"WO_input.tab")
        if path.exists(output+"WS_input.inp"):
            with open(output+"WS_input.inp", "r") as input:
                text = input.read()
                file.write("\n\nWater/solid (LS) input:\n")
                file.write(text)
            remove(output+"WS_input.inp")
        if path.exists(output+"OS_input.inp"):
            with open(output+"OS_input.inp", "r") as input:
                text = input.read()
                file.write("\n\nOil/solid (LS) input:\n")
                file.write(text)
            remove(output+"OS_input.inp")
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()

