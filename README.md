# ift_from_lle
IFT and surface coverage calculations
COSMOtherm is required to calculate solvation parameters, use a COSMO parameterization listed in check_parameterization in the main script

Call the main script by (remember adding your initals and COSMOtherm path to the code):
python "script_name" "COSMO_input_file" mix "initials"

The script can be controlled by following statements:
print_statements (Printing initial information, every iteration, and final result) default = True

debug (Additional information in every iteration) default = False

multiprocess (Using 2 cores during while loop) default = True

delete_files (Deleting files produced by COSMOtherm when the calculation is complete) default = True

Call the support script by 
(remember adding your initals and COSMOtherm path to the main code, and in the support script where the main script is called):
python "script_name" "COSMO_input_file"
