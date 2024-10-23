1.) To run the analysis, please run: "rspc_analyse.py" 
2.) You will be asked to choose the inserts for the analysis. You have "gammex", "new_cirs" and "old_cirs". ** The analysis can only work if you input the insert name in the correct format.

The script will perform three key analysis as follows:

1.) Calculated CT number
    Required input a.) The elemental mass attenuation coefficients from 1 keV to 120 keV (XCOM) 
                   b.) The measured insert Hounsfield Units (HUs) 
                   c.) The insert composition

    output -> An equation to estimate material CT numbers from the Siemens scanner based on known material compositions. (Reference: Schneider et al. (1996) paper ) 
           -> two figures. One to plot calc and meas ct numbers together. Another to plot the percent. diff. between calc and meas numbers.

The script will first calculate the effective energy of the scanner. We calculate the insert LACs at the effective energy. The K constants from the LAC equation can then be estimated.

*** This equation is scanner-specific.

2.) Calculated prsp >>> Required input a.) The insert composition data b.) The elemental i-value data from icru37 (compound state ) c.) The insert measured prsps.

3.) Stochiometric calibration using icru46 human tissue data >>> Required input 1.) The icru46 human composition data

1.) please use exact the same format to input/ edit data in the excel files.
2.) the insert name has to be in the correct format ( no extra space please). 
