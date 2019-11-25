from __future__ import print_function,division
import numpy as np
import pandas as pd
from ift_from_3phase import calculate_IFT_tot_and_coverage
from functions import change_input_name


def main():
    
    IFT_W_S = 19.57251427315493
    IFT_O_S = 6.882476613958883
    IFT_O_W = 46.87030111026727

    right_side = (IFT_O_S - IFT_W_S) / IFT_O_W
    
    if abs(right_side) > 1:
        print("Error: Can not take arccos to a number ouside the range [-1,1]. Please check the calculated energies.")
        print("The calculated number is {}.".format(right_side))
    else:
        contact_angle = np.arccos(-1) * (180/np.pi)
        print("Contact angle:", contact_angle)



















if __name__ == "__main__":
    main()

