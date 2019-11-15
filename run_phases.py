from __future__ import print_function,division
import sys
import os
import re
import pandas as pd
import numpy as np
from ift_from_lle_3phase_LVN import calculate_IFT_tot_and_coverage, get_comp_and_phases


def run_IFT(input_file, error_attempts):
    k = 1
    while k <= error_attempts:
        try:
            coverage, IFT = calculate_IFT_tot_and_coverage(input_file, "mix", "LVN", save_output_file = False)
            break
        except:
            print("An error occured, trying again. Try number {}/{}.".format(k, error_attempts))
            print(sys.stderr)
            k += 1
            continue
    return IFT, coverage

def main():
    pd.options.display.float_format = '{:.10f}'.format
    pd.set_option('display.max_rows', 50)
    pd.set_option('display.max_columns', 50)
    pd.set_option('display.width', 200)
    path_to_COSMOfiles = r"C:\Users\lasse\OneDrive\KU\Kandidat\Projekt\COSMO\COSMOfiles"
    input_file = sys.argv[1]
    output_path = ""
    error_attempts = 2
    
    # Changes .\input_file.inp -> input_file
    if input_file[:2] == ".\\" and input_file[len(input_file)-4:] == ".inp":
        input_file = input_file[2:len(input_file)-4]
        
    # If input file is a path on windows
    if re.findall("\w:", input_file) != []:
        input_file = input_file[:len(input_file)-4]
        output_path = os.path.split(input_file)[0]+"\\"
    
    # Find number of liquid extractions
    with open(input_file+".inp","r") as file:
        text = file.read()
        # N compounds
        obj_comp = re.findall(r"\{(?:\d*\.\d*[\ ,\}]*)+", text)  # Find first { with numbers behind it
        N_compounds = len(obj_comp[0].split())  # Numbers in the {} in the file
        # Find number of liquid extractions
        obj_inp = re.findall(r"liq_ex=\d", text)
        liq_ex = int(obj_inp[0][-1])
        
    
	# Initialize calculated concentration lists
    conc = [[] for _ in range(liq_ex)]
    ift_list = []
    coverage_list = []
    # Initial IFT calculation
    ift, coverage = run_IFT(input_file, error_attempts)
    ift_list.append(ift)
    coverage_list.append(coverage)
    
    # Make the header for printout
    header = ["phase_1", "surf_1_2", "phase_2"]
    
	# Read concentrations from liquid extraction
    with open(input_file+".tab", "r") as file:
        tab_lines = file.readlines()
        for i in range(-N_compounds, 0, 1):
            for j in range(liq_ex):
                conc[j].append(float(tab_lines[i].split()[j+2]))

	# If liq_ex = 3, run the second interface
    if liq_ex >= 3:
        phase_counter = liq_ex
        save_last_line = 1
        # Revert the last line
        for k in range(1, liq_ex-1):  # The first 2 phases is done in the initial calculation 
            header.extend(["surf_"+str(k+1)+"_"+str(k+2), "phase_"+str(k+2)])
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
                print("New input line:\n", last_line)
            # Write the revereted last line
            with open(input_file+".inp", "w") as write:
                write.writelines(read_lines[:-1])
                write.write(last_line)
                
            # Run the IFT on the reverted concentrations
            ift, coverage = run_IFT(input_file, error_attempts)
            ift_list.append(ift)
            coverage_list.append(coverage)
            
        # Set .inp file to initial conditions
        with open(input_file+".inp", "w") as write:
            write.writelines(read_lines[:-1])
            write.write(old_last_line)
        
    
	# Change lists to np arrays and reshape them
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

    
    # Print output in 3_phase_output.txt
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