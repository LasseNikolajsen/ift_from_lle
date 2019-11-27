from __future__ import print_function,division
import numpy as np
import pandas as pd
from ift_from_3phase import calculate_IFT_tot_and_coverage
from functions import change_input_name


def main():
    
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
    print("Contact angle:", contact_angle)



















if __name__ == "__main__":
    main()

