from __future__ import print_function,division
import subprocess
import sys
import os
import numpy as np
import re
from multiprocessing import Pool, cpu_count

def work(cmd):
    """ Run the process for multiprocessing
    
    Args:
        cmd: Command to the process in the function
    
    Return:
        The process
    """
    return subprocess.call(cmd, shell=False)
    

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

# Run by: python "scriptname" "input_file_name"(without extensions) mix "user initials"(in caps)
# Water should be called "h2o" and vacuum should be called "vacuum"


def calculate_IFT_tot_and_coverage(input_file_name, mix, user, print_statements = True, debug = False, multiprocess = True, delete_files = True, save_output_file = True):
    """ Calculate the total interfacial tension of the two input phases and 
        the surface coverage inbetween the phases.
    Args: 
        input_file_name: The name of the input file that the code should run without extension as a string 
        mix: Input type of IFT calculation as a string
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

    # Check unit=si in input file
    check_units(input_file_name)
    
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