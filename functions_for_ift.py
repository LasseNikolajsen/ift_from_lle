import subprocess
import sys
import os
import numpy as np
import re

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
    """
    # Changes .\input_file.inp -> input_file
    if name[:2] == ".\\" and name[len(name)-4:] == ".inp":
        name = name[2:len(name)-4]
        
    # If input file is a path on windows
    if re.findall("\w:", name) != [] and name[len(name)-4:] == ".inp":
        name = name[:len(name)-4]
    return name
    
def check_units(input_file_name):
    """ Check if unit=si is present in the .inp file
    
    Args:
        input_file_name: The input file name without extension as a string
    
    Return:
        None
    """
    with open(input_file_name+".inp","r") as file:
        text = file.read()
        if re.findall(r"unit\ *=\ *[sS][iI]", text) == []:
            print("Warning: unit=si is missing from the input file.")
    return
    
    
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

def check_phase_types(types):
    """ Check the input phase types
    
    Args:
        types: The input types as a string
        
    Return:
        types: As a formated string
    """
    
    # Check the input
    if len(types) != 2:
        print("Warning: Input types did not match the correct length of 2.")
        quit()
    if len(re.findall("[LlGgSs]", types)) != 2:
        print("Warning: Input types did not match phase types of liquid (L), gas (G) or solid (S).")
        quit()
    
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

    
def get_comp_and_phases(input_file_name, N_compounds):
    """ Extract data from the .tab file
    
    Args:
        input_file_name: The input file name without extension as a string
        N_compounds: The number of compounds in the system as an integer
        
    Return:
        compound_list: Compound names as a list
        phase1: Phase 1 as a list of floats
        phase2: Phase 2 as a list of floats
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
    
    
def write_flatsurf_file(input_file_name, output_input_file_name, phase1, phase2, T, IFT, IFT_write_length):
    """ Create new .inp files for flatsurf calculations
    
    Args:
        input_file_name: Filename as a string
        output_input_file_name: Output filename as string
        phase1: First phase as a list
        phase2: Second phase as a list
        T: Temperature as a float
        IFT: IFT as a float
    
    Return:
        None
    """
    with open(input_file_name+".inp", "r") as file:  # Read the inital input file
        lines = file.readlines()
        # Create output_input_file_name.inp file and write all lines except the last from initial file and write new last line
        with open(output_input_file_name+".inp", "w") as output:    
            output.writelines(lines[:-1])  # All lines except the last
            # Last line 
            (output.write("tk={0} FLATSURF xf1={{{1}}} xf2={{{2}}} IGNORE_CHARGE IFT={3:.{4}f}\n".
            format(T, "  ".join(map(str,phase1)), "  ".join(map(str,phase2)), IFT, IFT_write_length)))
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

    
def calculate_coverage(phase, Gtot, R, T):
    """ Calculate surface coverage between a phase and surface
    
    Args:
        phase: Phase as an array
        Gtot: The input phases Gtot as an array
        R: The gas constant in kj/mol/K as a float
        T: The temperature in Kelvin as a float
        
    Return:
        coverage: Surface coverage as an array
    """
    coverage = phase*np.exp(-Gtot/(R*T))
    return coverage 
    
    
def calculate_IFT(phase, Gtot_phase, area_phase, coverage, R, T, unit_converter):
    """ Calculate IFT between two phases
    
    Args:
        phase: Phase as a list
        Gtot_phase: Gtot from the phase as a list
        coverage: Surface coverage from the phase as a list
        R: The gas constant in kj/mol/K as a float
        T: The temperature in Kelvin as a float
        unit_converter: Converts the output to mN/m as a float
        
    Return:
        IFT: The sum of all interfacial tensions between phase and surface
    """
    coverage_part = coverage*(Gtot_phase-R*T*np.log(phase)+R*T*np.log(coverage))
    phase_part = phase*Gtot_phase
    IFT = np.sum((coverage_part + phase_part)/(2*area_phase)*unit_converter)
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