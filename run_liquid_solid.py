from __future__ import print_function,division
import sys
import numpy as np
import pandas as pd
from ift_from_3phase import calculate_IFT_tot_and_coverage
from functions import change_input_name, get_comp_and_phases, get_N_compounds_and_T


def main():
    
    # Use a 3 phase liquid extraction COSMOtherm input file and specify which phases are water, oil and sold.
    
    phase_types = "SWO"  # Water (O), Oil (O), Solid (S)
    
    WS_IFT = 0.0  # Water solid, if 0.0 run the calculation, else use specified value
    
    OS_IFT = 0.0  # Oil solid, if 0.0 run the calculation, else use specified value
    
    WO_IFT = 46.870301105580424  # Water oil, if 0.0 run the calculation, else use specified value
    
    input_file = sys.argv[1]
    
    input, output = change_input_name(input_file)
    print(output)
    N_comps, T = get_N_compounds_and_T(input)
    
    comp_list, phases = get_comp_and_phases(input, N_comps)
    
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

    
    
    with open(input+".inp", "r") as file: 
        text = file.readlines()
        if WS_IFT == 0.0:
            with open(output+"WS_input.inp", "w") as WS_file:
                WS_file.writelines(text[:3])
                for i in range(N_comps):
                    if i in water_compounds_index+solid_compounds_index:  # It needs to be multiple lines for compounds with multiple conformers
                        WS_file.write(text[3+i])
                last_line = text[-1]
                #print(last_line)
                WS_file.write(last_line)
        #print(text)
    
    # Water phase
    
    # Oil phase
    
    
    
    
    # n-decane
    # IFT_W_S = 64.9449187143066
    # IFT_O_S = -12.500284722858611
    # IFT_O_W = 46.870301105580424
    
    # 1-octanol
    # IFT_W_S = -18.626208743461117
    # IFT_O_S = 10.935139924034504
    # IFT_O_W = 46.870301105580424

    # acetic acid
    IFT_W_S = -14.774634533876746
    IFT_O_S = 9.751316977037327
    IFT_O_W = 46.870301105580424
    
    
    youngs_eq = (IFT_O_S - IFT_W_S) / IFT_O_W
    if youngs_eq > 1:
        print("Error: Can not take arccos to a number ouside the range [-1,1]. Please check the calculated energies.")
        print("The calculated number is {}.".format(youngs_eq))
        quit()
    elif youngs_eq < -1:
        youngs_eq = -1
        
    contact_angle = np.arccos(youngs_eq) * (180/np.pi)
    #print("Contact angle:", contact_angle)



















if __name__ == "__main__":
    main()

