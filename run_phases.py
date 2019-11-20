from __future__ import print_function,division
import sys
import os
import re
import traceback
import pandas as pd
import numpy as np
from ift_from_3phase import calculate_IFT_tot_and_coverage
from functions import change_input_name, get_comp_and_phases, check_phase_types

def run_IFT(input_file, error_attempts, phase_types, initials):
    """ Run IFT calculation again if a runtime error occurs

    Args:
        input_file: COSMOtherm input file for the IFT calculation
        error_attempts: The number of runtime errors the IFT script can encounter before terminating the calculation
        phase_types: The types of phases in the input file, liquid (L), gas (G) or solid (S)
		initials: Initials of the person running the script, so it can find the COSMOpath
        
    Return:
        IFT: The calculated IFT
        coverage: The calculated surface coverage
    """
    for k in range(1,error_attempts+1):
        try:
            coverage, IFT = calculate_IFT_tot_and_coverage(input_file, phase_types, initials, save_output_file = False)
            break
        except:
            print("An error occurred, trying again. Try number {}/{}.".format(k, error_attempts))
            traceback.print_exc()
            print(" \n")
            if k == error_attempts:
                quit()
            else:
                continue
    return IFT, coverage

def main():
    initials = "LVND"
    phase_types = "LLL"
    error_attempts = 2
    
    
    pd.options.display.float_format = '{:.10f}'.format
    pd.set_option('display.max_rows', 50)
    pd.set_option('display.max_columns', 50)
    pd.set_option('display.width', 200)
    path_to_COSMOfiles = r"C:\Users\lasse\OneDrive\KU\Kandidat\Projekt\COSMO\COSMOfiles"
    input_file = sys.argv[1]
    
    input_file, output_path = change_input_name(input_file)
    
    
    # Find number of liquid extractions
    with open(input_file+".inp","r") as file:
        text = file.read()
        # N compounds
        obj_comp = re.findall(r"\{(?:\d*\.\d*[\ ,\}]*)+", text)  # Find first { with numbers behind it
        N_compounds = len(obj_comp[0].split())  # Numbers in the {} in the file
        # Find number of liquid extractions
        obj_inp = re.findall(r"liq_ex=\d", text)
        liq_ex = int(obj_inp[0][-1])

    phase_types = check_phase_types(phase_types, liq_ex)
    # if liq_ex != len(phase_types):
        # print("Warning: Need to define the type of each phase in the calculation.")
        # quit()
    
	# Initialize calculated concentration lists
    conc = [[] for _ in range(liq_ex)]
    ift_list = []
    coverage_list = []
    # Initial IFT calculation
    ift, coverage = run_IFT(input_file, error_attempts, phase_types, initials)
    ift_list.append(ift)
    coverage_list.append(coverage)
    
    # Make the header for printout
    header = ["Phase 1 ({})".format(phase_types[0]), "Surface 1-2", "Phase 2 ({})".format(phase_types[1])]
    
    if phase_types == "LL":  # Read concentrations from liquid extraction in .tab file
        with open(input_file+".tab", "r") as file:
            tab_lines = file.readlines()
            for i in range(-N_compounds, 0, 1):
                for j in range(liq_ex):
                    conc[j].append(float(tab_lines[i].split()[j+2]))
    else:  # Read concentrations from .inp file
        phase = []
        with open(input_file+".inp","r") as file:
            lines = file.read()
            phase_object = re.findall(r"[^w]\d\ *=\ *\{[\d \. \ * e \-]*", lines)
            
            for i in range(len(phase_object)):
                phase = []
                for j in phase_object[i].split():
                    if re.findall(r"[^w]\d\ *=\ *\{", j) != []:
                        phase.append(float(j.split("{")[1]))
                    else:
                        phase.append(float(j))
                conc[i] = phase

	# If liq_ex >= 3, run the next interfaces
    if liq_ex >= 3:
        phase_counter = liq_ex
        save_last_line = 1
        # Revert the last line
        for k in range(1, liq_ex-1):  # The first 2 phases is done in the initial calculation 
            header.extend(["Surface {}-{}".format(k+1, k+2), "Phase {} ({})".format(k+2, phase_types[k+1])])
            with open(input_file+".inp", "r") as read:
                read_lines = read.readlines()
                if save_last_line == 1:
                    old_last_line = read_lines[-1]
                    save_last_line = 0
                last_line = ""
                phase_names = []
                phase_first_value = []
                for i in range(liq_ex):
                    phase_names.append(old_last_line.split()[2+i*N_compounds].split("=")[0]+"=")
                    phase_first_value.append(old_last_line.split()[2+i*N_compounds].split("=")[1])
                value_counter = 0
                phase_counter -= 1
                for j in old_last_line.split():
                    if "=" in j and "{" in j:
                        last_line += phase_names[phase_counter]
                        last_line += phase_first_value[value_counter]+" "
                        value_counter += 1
                        phase_counter += 1
                        if phase_counter == liq_ex:
                            phase_counter = 0
                    else:
                        last_line += j+" "
                print("\nNew input line:\n", last_line, "\n")
            # Write the reverted last line
            with open(input_file+".inp", "w") as write:
                write.writelines(read_lines[:-1])
                write.write(last_line)
                
            # Run the IFT on the reverted concentrations
            ift, coverage = run_IFT(input_file, error_attempts, phase_types[k:k+2], initials)
            ift_list.append(ift)
            coverage_list.append(coverage)
            
        # Set .inp file to initial conditions
        with open(input_file+".inp", "w") as write:
            write.writelines(read_lines[:-1])
            write.write(old_last_line)
        
    
	# Change lists to numpy arrays and reshape them
    for i in range(liq_ex):
        conc[i] = np.array(conc[i]).reshape((N_compounds,1))
        if i < liq_ex-1:
            coverage_list[i] = np.array(coverage_list[i].reshape((N_compounds,1)))
        
	# Concatenate arrays to a single matrix

    for i in range(liq_ex):
        if i == 0:  # First
            matrix = np.concatenate((conc[i], coverage_list[i]), axis=1)
        elif liq_ex == i+1:  # Last
            matrix = np.concatenate((matrix, conc[i]), axis=1)
        else:  # Rest
            matrix = np.concatenate((matrix, conc[i], coverage_list[i]), axis=1)


    
	# Create pandas data frame
    comp_list, _, _ = get_comp_and_phases(input_file, N_compounds)
    df = pd.DataFrame(matrix, columns=header, index=comp_list)

    
    # Print output in n_phase_output.txt
    with open(output_path+str(liq_ex)+"_phase_output.txt", 'w') as file:
        sys.stdout = file
        print(df)
        line = "\n"
        for i in range(len(ift_list)):
            line += "IFT_"+str(i+1)+"-"+str(i+2)+": "+str(ift_list[i])+"   "
        file.write(line)
        with open(input_file+".inp", "r") as input:
            text = input.read()
            file.write("\n\n\nInput:\n")
            file.write(text)
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()