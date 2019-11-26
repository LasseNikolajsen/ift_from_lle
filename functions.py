import subprocess
import sys
import os
import numpy as np
import re

# This document includes all the functions called in the IFT calculation script and some called in the run_phases support script

def get_liquid_index(phase1, phase2, phase_types):
    """ Get the indecies where the compounds are above 0.0 in the liquid phase
    
    Args:
        phase1: The first phase
        phase2: The second phase
        phase_types: The types of phase 1 and phase 2
        
    Return:
        liquid_index: The index for compounds in the liquid phase above 0.0, as a list
        solid_index: The index for compounds in the solid phase above 0.0, as a list
    """
    liquid_index = []
    solid_index = []
    
    if phase_types == "SL":
        for i in range(len(phase1)):
            if phase1[i] > 0.0:
                solid_index.append(i)
            else:
                liquid_index.append(i)
    if phase_types == "LS":
        for i in range(len(phase2)):
            if phase2[i] > 0.0:
                solid_index.append(i)
            else:
                liquid_index.append(i)
    return liquid_index, solid_index


def get_user_and_path(user_name):
    """ Get COSMOtherm path from Users.txt file or add a new user based on user input
    
    Args:
        user_name: User name as a string
    
    Return:
        user: User nickname as a string
        COSMOtherm_path: Path to the program as a real string
    """
    user_list = []
    path_list = []
    
    with open("Users.txt", "r") as file:
        text = file.read()
        user_object = re.findall("[Nn]ame:\ *\w*", text)
        path_obejct = re.findall("[Pp]ath:\ *[\S\ ]*", text)
        for (u, p) in zip(user_object, path_obejct):
            full_path = ""
            user_list.append(u.split()[1])
            p_split = p.split()
            for i in p_split[1:]:
                if i == len(p_split)-1:
                    full_path += iI
                else:
                    full_path += i + " "
            path_list.append(full_path)
    try:   # Try and find the user name in the Users.txt file
        index = user_list.index(user_name)
        COSMOtherm_path = path_list[index]
    except:  # If not found prompt the user to input a new name or terminate the script
        print("User name not recognized")
        print("Do you want to add a new name and COSMOtherm path to the Users.txt file? [Yes(y)/No(n)]")
        agreement = input()

        if agreement == "y":
            name = input("Name:")
            path = input("COSMOtherm path:")
            print("Name: ", name)
            print("Path:", path)
            with open("Users.txt", "a") as file:
                if path[-14:] != "cosmotherm.exe":
                    path += "cosmotherm.exe"
                file.write("\n \n")
                file.write("Name: "+name)
                file.write("\n")
                file.write("Path: "+path)
            print("Were added to your local version of Users.txt")
            COSMOtherm_path = path
        else:
            print("The script will terminate now")
            quit()
    if os.path.isfile(COSMOtherm_path) == False:
        print("Error: Could not find cosmotherm.exe at the specified path. Is the path correct for this computer or did you misspell something in the path?")
        quit()
    return COSMOtherm_path


def work(cmd):
    """ Run the process for multiprocessing
    
    Args:
        cmd: Command to the process in the function
    
    Return:
        The process
    """
    return subprocess.call(cmd, shell=False)
    

def change_input_name(name):
    """ Change the input file name from a path or with extension to the name without extension

    Args:
        name: The input file name as a string
        
    Return:
        name: The input file name without extension as a string
        path: The path to the input file if it has one
    """
    path = ""
    # Changes .\input_file.inp -> input_file
    if name[:2] == ".\\" and name[len(name)-4:] == ".inp":
        name = name[2:len(name)-4]
        
    # If input file is a path on windows
    if re.findall("\w:", name) != [] and name[len(name)-4:] == ".inp":
        name = name[:len(name)-4]
        path = os.path.split(name)[0]+"\\"
        
    return name, path
    
    
def check_units_get_liq_ex(input_file_name):
    """ Check if unit=si is present in the .inp file
    
    Args:
        input_file_name: The input file name without extension as a string
    
    Return:
        liq_ex: The number of phases in the calculation as an integer
    """
    with open(input_file_name+".inp","r") as file:
        text = file.read()
        if re.findall(r"unit\ *=\ *[sS][iI]", text) == []:
            print("Warning: unit=si is missing from the input file.")
        if len(re.findall(r"unit", text)) > 1:
            print("Warning: Multiple instances of unit in the input file, COSMOtherm might not use the correct units.")
    # Find number of liquid extractions
    with open(input_file_name+".inp","r") as file:
        text = file.read()
        # N compounds
        obj_comp = re.findall(r"\{(?:\d*\.\d*[\ ,\}]*)+", text)  # Find first { with numbers behind it
        N_compounds = len(obj_comp[0].split())  # Numbers in the {} in the file
        # Find number of liquid extractions
        obj_inp = re.findall(r"liq_ex=\d", text)
        liq_ex = int(obj_inp[0][-1])
    return liq_ex
    
    
def check_parameterization(input_file_name):
    """ Check parameterization between the water parameterization and the .inp file
    
    Args:
        input_file_name: The input file name without extension as a string
    
    Return:
        scale_water: Water scaling parameter as a float
    """
    parameter = ["BP_TZVP_C30_1601", "BP_TZVP_C30_1501", "BP_TZVP_C30_1401", "BP_TZVP_C30_1301", 
        "BP_TZVP_C21_0111", "DMOL3_PBE_C30_1301", "BP_TZVPD_FINE_C30_1501", "BP_TZVPD_FINE_C30_1401", 
        "BP_TZVPD_FINE_C30_1301", "add_parameter_here"]
    parameterization = [1/0.625, 1/0.25475*0.43061, 1/0.26753*0.43061, 1/0.25*0.43061,
        1/0.26697*0.43061, 1/0.2641*0.43061, 1, 1/0.31733*0.43061, 1/0.28649*0.43061, 
        "add_parameterization_here"]
    with open(input_file_name+".inp","r") as file:
        lines = file.readlines()
        try:
            index = parameter.index(lines[0].split()[2].split(".")[0])
            scale_water = parameterization[index]
        except:
            print("Warning: No matching parameterization found.\
            \nGo to check_parameterization to add new parameterizations.")
            quit()
    return scale_water, parameter[index]

    
def check_phase_types(types, N_phases):
    """ Check the input phase types
    
    Args:
        types: The input types as a string, liquid (L), gas (G) or solid (S)
        
    Return:
        types: As a formated string
    """
    # Check the input and call for new input if the input does not match the input file
    correct_phase = (len(types) != N_phases or len(re.findall("[LlGgSs]", types)) != N_phases)
    while(correct_phase):
        if(len(types) != N_phases):
            print("Warning: Phase types did not match the correct length of {}.".format(N_phases))
        elif(len(re.findall("[LlGgSs]", types)) != N_phases):
            print("Warning: Input types did not match phase types of liquid (L), gas (G) or solid (S).")
        else:
            print("Unknown error when checking the phase types")
            quit()
        print("Do you want to change the phase types for this calculation? Write X to quit")
        types = input("Write new types:")
        if types == "X":
            quit()
        if (len(types) == N_phases and len(re.findall("[LlGgSs]", types)) == N_phases):
            break
        else:
            continue

    # Format the input for future implementation
    types_formated = ""
    for i in types:
        if re.findall("[Ll]", i) != []:
            types_formated += "L"
            
        if re.findall("[Gg]", i) != []:
            types_formated += "G"
            
        if re.findall("[Ss]", i) != []:
            types_formated += "S"
    
    
    return types_formated


def get_N_compounds_and_T(input_file_name):
    """ Extract data from the .inp file
    
    Args:
        input_file_name: Filename as a string
    
    Return:
        N_compounds: Number of compounds as a float
        T: Temperature as floats
    """
    with open(input_file_name+".inp","r") as file:
        text = file.read()
        
        obj_T = re.findall(r"t[ck]=[0-9]+\.*[0-9]*", text)  # Find tc= or tk=
        T_list = obj_T[0].split("=")
        if T_list[0] == 'tc':
            T = float(T_list[1])+273.15
        else:
            T = float(T_list[1]) 

        obj_comp = re.findall(r"[^w]\d\ *=\ *\{[\d \. \ * e \-]*", text)  # Find first { with numbers behind it
        N_compounds = len(obj_comp[0].split())  # Numbers in the {} in the file

    return N_compounds, T

    
def get_comp_and_phases_for_LL(input_file_name, N_compounds):
    """ Extract data from the .tab file
    
    Args:
        input_file_name: The input file name without extension as a string
        N_compounds: The number of compounds in the system as an integer
        
    Return:
        compound_list: Compound names as a list
        phase1: Phase 1 as a np.array of floats
        phase2: Phase 2 as a np.array of floats
    """
    compound_list = []
    phase1 = []
    phase2 = []
    
    with open(input_file_name+".tab","r") as file:
        lines = file.readlines()

        # Find the index of phase_1_x and phase_2_x and use that index to get the values in their columns
        index1 = lines[-1-N_compounds].split().index("phase_1_x")
        index2 = lines[-1-N_compounds].split().index("phase_2_x")
        for i in range(-N_compounds,0,1):
            compound_list.append(lines[i].split()[1])
            phase1.append(float(lines[i].split()[index1]))
            phase2.append(float(lines[i].split()[index2]))
        phase1 = np.array(phase1)
        phase2 = np.array(phase2)
    return compound_list, phase1, phase2
    
    
def get_comp_and_phases(input_file_name, N_compounds):
    """ Extract data from the .inp file
    
    Args:
        input_file_name: The input file name without extension as a string
        N_compounds: The number of compounds in the system as an integer
        
    Return:
        compound_list: Compound names as a list
        phase1: Phase 1 as a np.array of floats
        phase2: Phase 2 as a np.array of floats
    """
    compound_list = []
    phase1 = []
    phase2 = []
    
    with open(input_file_name+".inp","r") as file:
        lines = file.read()
        compound_object = re.findall(r"f\ *=\ *[\w\S\ ]* VPfile", lines)
        #print(compound_object)
        phase_object = re.findall(r"[^w]\d\ *=\ *\{[\d \. \ * e \-]*", lines)
        for i in range(N_compounds):
            compound_list.append(compound_object[i].split()[2].split("_")[0])
        for j in phase_object:
            if j[1] == "1":
                for k in j.split():
                    if re.findall(r"[^w]\d\ *=\ *\{", k) != []:
                        phase1.append(float(k.split("{")[1]))
                    else:
                        phase1.append(float(k))
            if j[1] == "2":
                for k in j.split():
                    if re.findall(r"[^w]\d\ *=\ *\{", k) != []:
                        phase2.append(float(k.split("{")[1]))
                    else:
                        phase2.append(float(k))
        phase1 = np.array(phase1)
        phase2 = np.array(phase2)
    return compound_list, phase1, phase2
    
    
def write_flatsurf_file(input_file_name, output_input_file_name, phase1, phase2, T, IFT, IFT_write_length, phase_types, max_depth):
    """ Create new .inp files for flatsurf calculations
    
    Args:
        input_file_name: Filename as a string
        output_input_file_name: Output filename as string
        phase1: First phase as a list
        phase2: Second phase as a list
        T: Temperature as a float
        IFT: IFT as a float
        phase_types: Type of phases (Liquid L, Gas, G, Solid S) as a string
    
    Return:
        None
    """
    max_depth_str = ""
    if phase_types[0] == "S" or phase_types[1] == "S":
        max_depth_str = "maxdepth="+str(max_depth)
    
    with open(input_file_name+".inp", "r") as file:  # Read the inital input file
        lines = file.readlines()
        # Create output_input_file_name.inp file and write all lines except the last from initial file and write new last line
        with open(output_input_file_name+".inp", "w") as output:    
            output.writelines(lines[:-1])  # All lines except the last
            # Last line 
            (output.write("tk={0} FLATSURF xf1={{{1}}} xf2={{{2}}} IGNORE_CHARGE {5} IFT={3:.{4}f} \n".
            format(T, "  ".join(map(str,phase1)), "  ".join(map(str,phase2)), IFT, IFT_write_length, max_depth_str)))
    return
    

def get_Gtot_and_Area(input_file_name, N_compounds): 
    """ Extract data from the .tab file
    
    Args:
        input_file_name: The input file name without extension as a string
        N_compounds: The number of compounds in the system as an integer
    
    Return:
        GtotAB: Gtot from one side as a list of floats 
        GtotBA: Gtot from the other side as a list of floats
        AreaAB: Area from one side as a list of floats
        AreaBA: Area from the other side as a list of floats
    """
    GtotAB = []
    GtotBA = []
    AreaAB = []
    AreaBA = []
    with open(input_file_name+".tab","r") as file:
        text = file.read()
        for i in range(N_compounds):
            # find a line of numbers with more than 4 numbers
            obj_ABtab = re.findall(r"(?:[-+]?\d*\.\d*\s*){4,}", text)  
            # There is 2 times the N_compounds lines in obj_ABtab.
            # The first N_compounds lines are from one side, the rest are
            # from the other. Gtot is the index 1 and across,mean is index 2.
            GtotAB.append(float(obj_ABtab[i].split()[1]))
            GtotBA.append(float(obj_ABtab[i+N_compounds].split()[1]))
            AreaAB.append(float(obj_ABtab[i].split()[2]))
            AreaBA.append(float(obj_ABtab[i+N_compounds].split()[2]))
    GtotAB = np.array(GtotAB)
    GtotBA = np.array(GtotBA)
    AreaAB = np.array(AreaAB)
    AreaBA = np.array(AreaBA)
    return GtotAB, GtotBA, AreaAB, AreaBA
  

def scale_area(compound_list, AreaAB, AreaBA, N_compounds, scale_water, scale_organic):
    """ Scaling areas
    
    Args:
        compound_list: Compound names as a list
        AreaAB: Areas from one side as a list of floats
        AreaBA: Areas from the other side as a list of floats
        N_compounds: The number of compounds in the system as an integer
        scale_water: The parameterization for water as a float or integer
        scale_organic: The parameterization for organic as a float or integer
        
    Return:
        AreaAB: Scaled area from one side as a list of floats
        AreaBA: Scaled area from the other side as a list of floats
    """
    for i in range(N_compounds):
        if (compound_list[i]=='h2o'):  # Scale water
            AreaAB[i]*=scale_water;
            AreaBA[i]*=scale_water;
        elif (compound_list[i]=='vacuum'):  # Scale vacuum
            AreaAB[i]=1e1000
            AreaBA[i]=1e1000
        else:
            AreaAB[i]*=scale_organic;
            AreaBA[i]*=scale_organic;
    return AreaAB, AreaBA  

    
def calculate_coverage(phase, Gtot, R, T, liquid_index):
    """ Calculate surface coverage between a phase and surface
    
    Args:
        phase: Phase as an array
        Gtot: The input phases Gtot as an array
        R: The gas constant in kj/mol/K as a float
        T: The temperature in Kelvin as a float
        liquid_index: The index for the liquid phase, if a solid phase is present
        
    Return:
        coverage: Surface coverage as an array
    """
    coverage = np.zeros(len(phase))
    for i in range(len(phase)):
        if i in liquid_index:
            coverage[i] = phase[i]*np.exp(-Gtot[i]/(R*T))
        else:
            coverage[i] = 0.0
    return coverage 
    
    
def calculate_IFT(bulk_phase, Gtot_phase_bulk_surface, Gtot_phase_surface_bulk, area_phase_bulk_surface, area_phase_surface_bulk, coverage, R, T, unit_converter, phase_types, liquid_index):
    """ Calculate IFT between two phases for either Liquid (L) - Surface Coverage (C)c Gas (G) - Surface Coverage (C) or Solid (S)
    
    Args:
        bulk_phase: Phase as a list
        Gtot_phase_bulk_surface: Gtot from the bulk phase to the surface as a list
        Gtot_phase_surface_bulk: Gtot from the surface to the bulk face as a list
        area_phase_bulk_surface: Area from the bulk phase to the surface as a list
        area_phase_surface_bulk: Area from the surface to the bulk face as a list
        coverage: Surface coverage from the phase as a list
        R: The gas constant in kj/mol/K as a float
        T: The temperature in Kelvin as a float
        unit_converter: Converts the output to mN/m as a float
        phase_types: Type of phases (Liquid L, Gas, G, Solid S, Coverge C) as a string
        
    Return:
        IFT: The sum of all interfacial tensions between phase and surface as a float
    """
    if phase_types[0] == "L" or phase_types[1] == "L":
        coverage_part = coverage[liquid_index]*(Gtot_phase_bulk_surface[liquid_index]-R*T*np.log(bulk_phase[liquid_index])+R*T*np.log(coverage[liquid_index]))
        phase_part = bulk_phase[liquid_index]*Gtot_phase_bulk_surface[liquid_index]
        IFT = np.sum((coverage_part + phase_part)/(2*area_phase_bulk_surface[liquid_index])*unit_converter)
    elif phase_types[0] == "G" or phase_types[1] == "G":
        coverage_part = coverage*(Gtot_phase_bulk_surface-R*T*np.log(bulk_phase)+R*T*np.log(coverage))
        phase_part = bulk_phase*Gtot_phase_bulk_surface
        IFT = np.sum((coverage_part + phase_part)/(2*area_phase_bulk_surface)*unit_converter)
    elif phase_types[0] == "S" or phase_types[1] == "S":
        phase_part = bulk_phase*Gtot_phase_surface_bulk
        IFT = np.sum(phase_part/(2*area_phase_surface_bulk)*unit_converter)
    return IFT

    
def calculate_IFT_dampning(IFT, IFT_value, IFT_max_diff, IFT_dampning):
    """ Calculate IFT direction and dampen the value
    
    Args:
        IFT: As calculated by calculate_IFT as a float
        IFT_value: IFT_value as a float
        IFT_max_diff: The maximum difference i.e. step size as a float or integer
        IFT_dampning: The dampning effect as a float
    
    Return:
        IFT_value with 6 decimals, truncated to prevent memory error in COSMOthermX18
    """
    if IFT-IFT_value > IFT_max_diff:
        difference = IFT_max_diff
    elif IFT-IFT_value < -IFT_max_diff:
        difference = -IFT_max_diff
    else:
        difference = IFT-IFT_value
    
    IFT_value = IFT_value+difference*IFT_dampning
    return IFT_value