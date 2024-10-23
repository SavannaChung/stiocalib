import numpy as np
import pandas as pd

import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import shutil

import easygui as eg

import rspc_read_file as rf
import rspc_calc as rc
# from rspc_calc import Material,calc_mat_ctno_from_lac # water
import rspc_plot as rp
# import rspc_report as rr
# Read the xlsx containing the insert measured CT number


#======================================================================#
#----------------------------- CONSTANTS ------------------------------#
#======================================================================#
# atomic number
A_NUM = {"h":1, "c":6, "n":7, "o":8,"f":9, "na":11, "mg":12, "p":15, "s":16,\
         "cl":17, "k":19, "ca":20, "sb":51, "sn":50, "fe":26, "i":53, "si": 14, "gd":64 }
# atomic mass
A_MASS =  {"h":1.0079, "c":12.011, "n":14.006, "o":15.999, "f":18.998, "na":22.989, "mg":24.312, "p":30.973, "s":32.064,\
           "cl":35.450, "k":39.102, "ca":40.080, "sb":121.76, "sn":118.71, "fe":55.847, "i":126.9, "si": 28.085, "gd":157.25}
# Avogadro constant
NA = 6.02214e23

water_den = 1.0
water_comp =  {"h":11.2, "o":88.8}

# proton energy in MeV
PRO_EN = 200

# mass of a proton/ an electron in kg
MASS_P = 1.67e-27
MASS_E = 9.11e-31

# speed of light ms-2
SPEED_OF_LIGHT = 299792458

# water i-value from icru 90 (eV)
WATER_IVAL = 78

# CONSTANTS for PRSP
BETA = np.sqrt(2*(PRO_EN*1e6*1.6e-19)/MASS_P )/SPEED_OF_LIGHT
MEC2 = (MASS_E*SPEED_OF_LIGHT**2 )/1.6e-19 # unit = eV -> matching the I_value unit

# I_value from Janni , unit =eV 
# I_VALUE = {"h":20.4, "c":73.8, "n":97.8, "o":115.7, "f":124.8, \
#            "na":143, "mg":151.1, "p":179.1, "s":183.6, "si":174.5 , "cl":182.6, \
#            "k":186.8, "ca":191.9, "sb":488, "sn":500.2, "fe":278.2, "i":515.2}

# i-value from ICRU37, unit = eV, i-values when elements are in compound state
# we took the i-value from 
# h = 19.2 eV (page 24, table 5.1 icru 37)
# c = 81 eV (page 24, table 5.1 icru 37)
# n = 82 eV (page 24, table 5.1 icru 37)
# o = 106 eV (page 24, table 5.1 icru 37)
# f = 112 eV (page 24, table 5.1 icru 37)
# others = 1.13 * I (eV from table 4.3)
# na = 168.37 eV (149*1.13, page 17, table 4.3, icru 37)
# mg = 176.38 eV (156*1.13, page 17, table 4.3, icru 37)
# p = 195.49 eV (173*1.13, page 17, table 4.3, icru 37)
# s = 203.4 eV (180*1.13, page 17, table 4.3, icru 37)
# si = 195.49 eV (173*1.13, page 17, table 4.3, icru 37)
# cl = 180 eV (page 24, table 5.1 icru 37)
# k = 214.7 eV (190*1.13, page 17, table 4.3, icru 37)
# ca = 215.83 eV (191*1.13, page 17, table 4.3, icru 37)
# fe = 323.18 eV (286*1.13, page 17, table 4.3, icru 37)
# i = 554.83 eV (491*1.13, page 17, table 4.3, icru 37)
# gd = 667.83 eV (591*1.13, page 17, table 4.3, icru 37)
# sb = 550.31 eV (487*1.13, page 17, table 4.3, icru 37)
# sn = 551.44 eV (488*1.13, page 17, table 4.3, icru 37)



I_VALUE=  {"h":19.2, "c":81, "n":82.0, "o":106.0, "f": 112, \
         "na":168.4 , "mg":176.3 , "p":195.5, "s": 203.4, "si":195.5, "cl":180,\
          "k":214.7, "ca":215.8, "sb":550.3 , "sn":551.4, "fe":323.2, "i":554.8, "gd":667.8 }

# i-value from ICRU37, unit = eV, i-values when elements are in gas/ condensed state (page 17, table 4.3, icru37)
# I_VALUE=  {"h":19.2, "c":78, "n":82, "o":95, "f": 115, \
#          "na":149 , "mg":156 , "p":173, "s": 180, "si":173, "cl":174,\
#           "k":190, "ca":191, "sb":487 , "sn":488, "fe":286, "i":491, "gd":591 }


# >>> define the tissue composition data tc_source
# "tc_eb_newborn_female"
# "tc_eb_1yo_female"
# "tc_eb_5yo_female"
# "tc_eb_10yo_female"
# "tc_eb_15yo_female"
# "tc_eb_adult_female"
# tc_source = "tc_icrp37" # icrp37 data
# tc_source = "tc_icru46" # icru46 data

# tc_source = "tc_eb_5yo_female" # tissue compositions of a 5 years old female
# tc_source = "tc_eb_5yo_male"
# tc_source = "tc_icrp110"
# tc_source = "tc_icrp110"
# tc_source = "tc_icru46"
tc_source = "tc_beam_model"



def main():
    # >>> define the phantom used
    insert_name = input("enter the insert name: gammex / old_cirs / new_cirs ")
    print(f"insert_name: {insert_name}")

    # navigate to the work directory
    file_input = "rsp_input.xlsx"
    lib_file_path = eg.fileopenbox("Please select your your input!  >>> hint: rsp_library.xlsx ")
    path, file_library = os.path.split(lib_file_path)
    os.chdir(path)

    # excel_library = pd.ExcelFile(lib_file_path)
    # ls_tc_libraries = excel_library.sheet_names
    # tc_source = eg.choicebox("select your tissue composition data interested.", "A list of sheet names.", ls_tc_libraries )



    #======================================================================#
    #------------------------------Read data ------------------------------#
    #======================================================================#

    # >>> load the measured ctnos for gammex and cirs inserts
    meas_ctno = rf.read_data_meas_ctno(file_input, "meas_ctno" ,insert_name)

    sn_insert_comp = "_".join(["insert_comp", insert_name])
    insert_dens, insert_comp = rf.read_data_mat_comp(file_library, sn_insert_comp)

    # >>> read all mass attenution coefficients in a nested dictionary
    df_mac = pd.read_excel(file_library,sheet_name = "xcom" )
    df_mac.dropna(inplace = True) # drop nan items
    df_mac = df_mac.set_index("pho_en_MeV") # set photon energy as the index of the dataFrame
    dt_mac = df_mac.T.to_dict()

    # >>> read the tissue composition data
    df_icrp = pd.read_excel(file_library,sheet_name = tc_source , header = 0, index_col =0)
    df_icrp.dropna( inplace = True) # drop nan items

    icrp_den = df_icrp['density'].to_dict()
    df_icrp = df_icrp.drop( columns = ['density'])
    icrp_comp = df_icrp.to_dict('index')

    # >>> read the measured prsps
    mp_sheet_name = "_".join(["mprsp", insert_name])
    mprsp = rf.read_mprsp(file_input, sheet_name = mp_sheet_name)

    #======================================================================#
    #----------------------- create directory for tc ------------------------#
    #======================================================================#
    if os.path.isdir(tc_source):
        os.chdir(os.path.join(path, tc_source))
    else:
        os.mkdir(os.path.join(path, tc_source))
        os.chdir(os.path.join(path, tc_source))

    shutil.copyfile(os.path.join(path, file_input), os.path.join(os.path.join(path, tc_source), file_input))
    shutil.copyfile(os.path.join(path, file_library), os.path.join(os.path.join(path, tc_source), file_library))

    # copy the input and library used for the RSP calibration



    #======================================================================#
    #------------------------------ Analysis ------------------------------#
    #======================================================================#

    # >>> calculate the effective energy (ee), ct no at ee, insert lac at ee
    ctno_lsm, calc_ctno_e, mat_lac_e, ee= rc.calc_mat_ctno_from_lac(insert_dens, insert_comp, dt_mac, meas_ctno)
    md, mc = rc.find_com_mat(mat_lac_e, insert_dens, insert_comp)
    materials = rc.make_materials(md, mc)
    ls_lac_en = list(mat_lac_e.values())

    popt, pcov = curve_fit(rc.calc_lacs, materials, ls_lac_en, bounds=([0, 0, 0], [1, 1, 1]))
    calc_ctno_inserts = rc.calc_ctno_k(materials, popt[0], popt[1], popt[2])

    # >>> calculate the insert prsps
    md_p, mc_p = rc.find_com_mat(mprsp, insert_dens, insert_comp)
    cprsp = rc.calc_prsp(rc.make_materials(md_p, mc_p))

    # >>> icrp human tissue data
    mat_icrp = rc.make_materials(icrp_den, icrp_comp)
    calc_hu =  rc.calc_ctno_k(mat_icrp, popt[0], popt[1], popt[2])
    calc_prsp = rc.calc_prsp(mat_icrp)


    #======================================================================#
    #------------------------------- plots --------------------------------#
    #======================================================================#

    # >>> plot the effective energy of the scanner
    rp.plot_lsm_pho_en(ctno_lsm, insert_name )

    # >>> plot the material ct no in an ascending order
    # rp.plot_sorted_calc_ctno(calc_ctno)

    # # >>> plot the calc ( using xcom database ) and meas CT numbers
    # rp.plot_xcalc_meas_ctno(calc_ctno_e, meas_ctno, insert_name, ee)ga

    # # >>> plot the calc ( using k-constants) and meas CT numbers
    rp.plot_kcalc_meas_ctno(calc_ctno_inserts, meas_ctno, insert_name, \
                                ee, "{:.2e}".format(popt[0]), "{:.2e}".format(popt[1]), "{:.2e}".format(popt[2]))
    #
    # # >>> plot the calc and meas prsp
    rp.plot_calc_meas_prsp(cprsp, mprsp, insert_name)
    # # print(rc.calc_prsp(rc.make_materials(md_p, mc_p)))
    #
    # >>> plot the final icrp result (calc prsp as a function of calc ct no)
    slope, const = rp.plot_prsp_ctno(calc_hu, calc_prsp, insert_name)

    # >>> plot the final rsp with the valiated hu number
    # rp.plot_rsp_w_mea_hu(calc_hu, calc_prsp, insert_name)


    #======================================================================#
    #----------------------- Export data to excel ------------------------#
    #======================================================================#

    EL_DATA = {"Atomic number":A_NUM, "Atomic mass":A_MASS, "i values (icru37)": I_VALUE}
    EL_DATA = pd.DataFrame.from_dict(EL_DATA).T

    EL_CONST = {"Avogadro constant": NA, "water density (g/cm3)": water_den, \
                "proton mass (kg)": MASS_P , "electron mass (kg)": MASS_E, "speed of light": SPEED_OF_LIGHT, \
                 "WATER_IVAL (eV)": WATER_IVAL}
    EL_CONST = pd.DataFrame(pd.Series(EL_CONST, name = "values"))

    # Get insert electron densities
    ACC_INSERTS = list(dict(sorted(md.items(), key = lambda item: item[1])).keys())

    insert_edens = {key:materials[key].edens for key in ACC_INSERTS}

    RSP_CONST = pd.DataFrame([slope, const], index = ["m", "c"], columns = ["lung", "soft tissues", "bones"] )

    # rsp_const = {"lung"}

    INSERT_DATA = pd.DataFrame.from_dict(mc).reindex(columns = ACC_INSERTS)
    INSERT_DATA = INSERT_DATA.append(pd.Series(md, name = "insert density(g/cm3)" ))
    INSERT_DATA = INSERT_DATA.append(pd.Series(insert_edens, name = "electron density (e-/cm3)"))
    INSERT_DATA = INSERT_DATA.append(pd.Series(meas_ctno, name = "meas. CTno"))
    INSERT_DATA = INSERT_DATA.append(pd.Series(calc_ctno_inserts, name = "calc ctno"))
    INSERT_DATA = INSERT_DATA.append(pd.Series(mprsp, name = "meas prsp"))
    INSERT_DATA = INSERT_DATA.append(pd.Series(cprsp, name = "calc prsp"))

    icrp_edens = {key:mat_icrp[key].edens for key in icrp_den.keys()}

    ICRP_DATA = pd.DataFrame.from_dict(icrp_comp)
    ICRP_DATA = ICRP_DATA.T
    ICRP_DATA["tissue density (g/cm3)"] = pd.Series(icrp_den)
    ICRP_DATA["electron density (e-/cm3)"] = pd.Series(icrp_edens)
    ICRP_DATA["calc_hu"] = pd.Series(calc_hu)
    ICRP_DATA["calc_prsp"] = pd.Series(calc_prsp)

    KCONST_DATA  = pd.DataFrame(pd.Series({"Kph":popt[0], "Kcoh":popt[1], "KKN":popt[2]}, name = "k constants"))

    # >>> writing out data
    excel_filename = "data_%s_%s.xlsx" % (insert_name, tc_source)
    with pd.ExcelWriter(excel_filename) as excel_writer:
        EL_DATA.to_excel(excel_writer, sheet_name = "element_details")
        EL_CONST.to_excel(excel_writer, sheet_name = "list_of_constants")
        KCONST_DATA.to_excel(excel_writer, sheet_name = "k_constants")
        INSERT_DATA.to_excel(excel_writer, sheet_name = "insert_details")
        ICRP_DATA.to_excel(excel_writer, sheet_name = "tc_data" )
        RSP_CONST.to_excel(excel_writer, sheet_name = "rsp_constants")

    print(f">>> The data are written out in the {excel_filename} file ")





if __name__=="__main__":
    main()
