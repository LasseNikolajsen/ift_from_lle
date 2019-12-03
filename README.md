# ift_from_lle
IFT and surface coverage calculations
COSMOtherm is required to calculate solvation parameters, use a COSMO parameterization listed in check_parameterization in the main script

Call the main script by:
python "script_name" "COSMO_input_file" "phase_types "user_name"

script_name: Is the name of the .py file
COSMO_input_file: Is the name of the LLE input file generated from COSMOtherm
phase_types: Is the types of phases in the calculation (L for liquid, S for solid and G for gas). It should be two letters, so a liquid liquid IFT calculation would be "LL".
initials: Is the name of the user, which should correspond to a COSMOtherm path in the Users.txt file

The script can be controlled by following statements:

print_statements: Printing initial information, every iteration, and final result, default = True
debug: Additional information in every iteration, default = False
multiprocess: Using 2 cores during while loop if possible, default = True
delete_files: Deleting files produced by COSMOtherm when the calculation is complete, default = True
save_output_file: Save the raw output file from the calculation, defalut = True
forced_convergence: Force the script to end the iterative process at max_iterations, default = False
max_iterations: Sets the maximum iterations in the iterative process if forced_convergence is set to True, as an integer

There are two support scripts, which can help in certain calculation situations.
First is the run_multi_L_phases.py, which runs n liquid phases and prints the results in an easy to overview output file, including the input file. The input file can still just be generated as a LLE input from COSMOtherm.
Run it by specifying the user and phase types inside the script and call: python run_multi_L_phases.py "input_file_name"

Second is the run_liquid_solid.py, which runs a complete contact angle experimental calculation, including water/oil, water/solid and oil/solid calculations and finally the estimated contact angle. The input file can still just be generated as a LLE input from COSMOtherm.
Run it by specifying the user and phase types inside the script and call: python run_multi_L_phases.py "input_file_name"
Here the phase types are: Water (W), oil (O) and solid (s).

