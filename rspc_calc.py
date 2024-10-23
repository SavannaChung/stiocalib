# # rspc
import numpy as np
import pandas as pd

from scipy.optimize import curve_fit

import rspc_analyse as ra

# Read the atomic number, mass number and material composition in dictionary
class Material:
    """Class definining a material object.

    Give density (float) and elemental composition (dictionary).
    Init method automatically calculates electron density"""

    def __init__(self, density, composition):
        self.density = density
        self.composition = composition # where this composition is a dictionary of elemental fractions

        # # define constants used in this Class
        # self.A_NUM = {"h":1, "c":6, "n":7, "o":8,"f":9, "na":11, "mg":12, "p":15, "s":16, \
        #         "cl":17, "k":19, "ca":20, "sb":51, "sn":50, "fe":26, "i":53, "si": 14 }
        # self.A_MASS =  {"h":1.0079, "c":12.011, "n":14.006, "o":15.999, "f":18.998, "na":22.989, "mg":24.312, "p":30.973, "s":32.064,\
        #         "cl":35.450, "k":39.102, "ca":40.080, "sb":121.76, "sn":118.71, "fe":55.847, "i":126.9, "si": 28.085}
        # self.NA = 6.02214e23

        # Use class method to calculate electron density
        self.edens = self.calc_edens()
        self.lda = self.calc_lambda()
        self.lda_sum = sum(self.lda.values())
        self.ean_za, self.ean_zc = self.calc_ean()
        self.ival_mat = self.calc_ival_mat()

    def calc_edens(self):

        sum_a = 0
        for el in self.composition:
            sum_a += (self.composition[el]/100.0) * ra.A_NUM[el] / ra.A_MASS[el]
        electron_density = sum_a * self.density * ra.NA

        return electron_density

    # Calculate the effective atomic number Z_a (EAN_approx) and Z_c(EAN_circumflex)
    def calc_lambda(self):
        lda_el = {}
#         lda_sum = 0

        for el in self.composition:
            lda_el[el] = self.composition[el]/100.0 * ( ra.A_NUM[el] / ra.A_MASS[el])

        return lda_el


    def calc_ean(self):
        Z_a = 0
        Z_c = 0
        for el in self.lda:
            Z_a += (self.lda[el]/self.lda_sum)*(ra.A_NUM[el])**(3.62)
            Z_c += (self.lda[el]/self.lda_sum)*(ra.A_NUM[el])**(1.86)

        return Z_a**(1/3.62), Z_c**(1/1.86)

    # def the I_values of different materials
    def calc_ival_mat(self):
        """ to calculate the i-value of different materials,
        you need to input:
        mat_object in Material(density, {composititon})

        """
        ival_sum = 0
        for el in self.composition:
            ival_sum += (self.lda[el]*np.log(ra.I_VALUE[el]))/self.lda_sum
        ival_mat = np.exp(ival_sum)
        return ival_mat

def calc_prsp(mo):
    """ calculate material prsps.
    -------------------------------input-------------------------------------
    mo : material objects in a dictionary
    -------------------------------output-------------------------------------
    mprsp: material prsp in a dictionary
    """
    water = Material(ra.water_den, ra.water_comp)

    mprsp = {mat:0 for mat in mo.keys()}

    for mat in mo.keys():
        sp_mat = (np.log((2*ra.MEC2*ra.BETA**2)/(mo[mat].ival_mat *(1 - ra.BETA**2))) - ra.BETA**2)
        # sp_water =  (np.log((2*ra.MEC2*ra.BETA**2)/(water.ival_mat *(1 - ra.BETA**2))) - ra.BETA**2)
        sp_water =  (np.log((2*ra.MEC2*ra.BETA**2)/(ra.WATER_IVAL *(1 - ra.BETA**2))) - ra.BETA**2)
        mprsp[mat] = (mo[mat].edens/water.edens)*(sp_mat/sp_water)
        # print(f"mat: {mat} rel. edens: {mo[mat].edens/water.edens}")
    return mprsp

def calc_mat_ctno_from_lac(mat_dens, mat_comp, dt_mac, dt_meas_ctno):
    """ calculate the linear attenuation coefficeints for different inserts
        theoretical input:
                            mat_dens : insert densities in dictionary
                            mat_comp: insert compositions in dictionary
                            dt_mac: element mac in dictionary
        measured input:
                            dt_meas_ctno : measured insert ctno in dictionary
        -----------------------------------------------------------------
        output: a few nested dictionaries

        dt_ctno_lsm = {energy1: lsv1, energy2: lsv2 ...}
        dt_mat_calc_ctno_effen = {insert_1: calc ctno@ effective energy, ...}}
        dt_mat_lac_eff_en = {insert_1: calc_lac_1 @ eff_en, insert_2: calc lac2 @ eff_en, ... }
        eff_en = a single float number
    """
    # define water object
    water = Material(ra.water_den, ra.water_comp)


    # list the common materials between the meas_ctno and calc_ctno
    ls_com_mat = list(sorted(set(dt_meas_ctno).intersection(set(mat_dens))))
    # print(ls_com_mat)

    # create a np array to cover the entire photon energy
    ls_pho_en = [keys for keys in dt_mac.keys()]
    # print(f"list of photon energy: \n :{ls_pho_en} \n")

    # create an empty dictionary to store the lac corresponding to each material
    dt_mat_lac = {keys: {e:0} for keys in ls_com_mat for e in ls_pho_en}
    dt_water_lac = {}

    # print(f"empty material dictionary: \n {dt_mat_lac} tyep: {type(dt_mat_lac)} \n")

    ##  ======= >>>  calculate the linear attenuation coefficient of water as a function of photon energy <<< =======
    for e in ls_pho_en:
        sum_mac = 0
        for el in list(water.composition.keys()):
            sum_mac += (water.composition[el]/100)*dt_mac[e][el]
        dt_water_lac[e] =sum_mac*water.density

    # print(f"water lac at 0.02= {dt_water_lac[0.02]}")

    ##  ======= >>> calculate the linear attenuation coefficient of different materials <<< =======
    sum_mac = 0
    for mat in ls_com_mat:
        for e in ls_pho_en:
            # print(f"pho_en: {e}")
            sum_mac =0
            for el in mat_comp[mat]:
                sum_mac += (mat_comp[mat][el]/100)*dt_mac[e][el]
            dt_mat_lac[mat][e] =  sum_mac*mat_dens[mat]
        # print(f"calc lac for materials {mat}: \n\n {dt_mat_lac[mat]} \n")

    ##  ======= >>> Calculate the theoretical CT numbers  <<< =======
    dt_calc_ctno = {keys: {e:0} for keys in ls_com_mat for e in ls_pho_en}
    for mat in ls_com_mat:
        ctno = 0
        for e in ls_pho_en:
            ctno = 1000*(dt_mat_lac[mat][e]/ dt_water_lac[e]) # ct no
            ctno = 1000*(dt_mat_lac[mat][e]/ dt_water_lac[e]) -1000 # Hounsfield unit
            dt_calc_ctno[mat][e] = ctno
    df_mat_calc_ctno = pd.DataFrame.from_dict(dt_calc_ctno)


    ##  ======= >>> Find the effective photon energy of the CT scanner  <<< =======
    dt_ctno_lsm = {keys:0 for keys in ls_pho_en}
    for e in ls_pho_en:
        sum_lsm = 0
        for mat in ls_com_mat:
            # print(f"dt_calc_ctno[mat][e]:{dt_calc_ctno[mat][e]} || dt_meas_ctno[mat]:{dt_meas_ctno[mat]} \n")
            sum_lsm += (dt_calc_ctno[mat][e]- dt_meas_ctno[mat])**2
        dt_ctno_lsm[e] = sum_lsm

    # print(f"final dt_ctno_lsm : {dt_ctno_lsm}")

    ##  ======= >>> Find the effective photon energy of the CT scanner  <<< =======
    eff_en = min(dt_ctno_lsm, key = dt_ctno_lsm.get)
    print(f"effective energy = {eff_en}")

    ##  ======= >>> Find the material calc ctno at effective energy  <<< =======
    mat_calc_ctno_effen = df_mat_calc_ctno.loc[eff_en].to_dict()
    # print(f"mat_calc_ctno_effen at eff_en : {mat_calc_ctno_effen} \n type(mat_calc_ctno_effen) : {type(mat_calc_ctno_effen)}")

    ##  ======= >>> get the linear attenuation coefficients of different materials at the effective energy  <<< =======
    dt_mat_lac_effen = {keys:0 for keys in ls_com_mat}
    for mat in ls_com_mat:
        dt_mat_lac_effen[mat] = dt_mat_lac[mat][eff_en]

    # print(f"rspc_calc, dt_mat_lac : {dt_mat_lac_effen}")

    return  dt_ctno_lsm, mat_calc_ctno_effen, dt_mat_lac_effen, eff_en

#================================= Data fitting =================================#
def find_com_mat(lac_effen, mat_dens, mat_comp):
    """ create two dictionaries (one for density, another one for composition of
        the selected materials)
        -------------------------------input-------------------------------------
        lac_effen : material linear atten. coeffs. at eff. en.
        mat_dens, mat_comp : material densities, compositions in dictionaries (general)
        -------------------------------output-------------------------------------
        md, mc : material densities and compositions in dictionaries
    """
    # create a dictionary for common material deneities and compositions
    ls_mat = list(lac_effen.keys())

    md = {keys:0 for keys in ls_mat}
    mc = {keys:0 for keys in ls_mat}

    for mat in ls_mat:
        md[mat] = mat_dens.get(mat)
        mc[mat] = mat_comp.get(mat)

    return md, mc

def make_materials(mat_dens, mat_comp):
    """ create a list of materials having meas. ct no.
        -------------------------------input-------------------------------------
        mat_dens, mat_comp : material densities, compositions in dictionaries
        -------------------------------output-------------------------------------
        mat_obj : material objects (from Material class) in a dicitonary.
    """
    # create material object using Material class
    mat_obj = {keys:0 for keys in mat_dens.keys()}

    for mat in mat_dens.keys():
        print(f'material: {mat}')
        mat_obj[mat] = Material(mat_dens[mat], mat_comp[mat])

    return mat_obj

# define the fitting function to calculate the k-constants
def calc_lacs(materials, Kph, Kcoh, KN):
    """ using this function to estimate the k-constants
        -------------------------------input-------------------------------------
        materials : a dictionary of material_objects from the materials we measured ct numbers from make_materials()
                    ( we may have iserts in the composition tables but we did not measure their ct no.)
        -------------------------------output-------------------------------------
        lacs : a list of estimated lac values
    """
    lacs = []
    for mat in materials.keys():
        Za = materials[mat].ean_za
        Zc = materials[mat].ean_zc

        lac = (materials[mat].edens / 1e24) * (Kph*Za**3.62 + Kcoh*Zc**1.86 + KN)
        lacs.append(lac)

    return lacs

def calc_ctno_k(materials, Kph, Kcoh, KN):
    """ using this function to calculate HU from the k constants
        -------------------------------input-------------------------------------
        materials : a dictionary of material_objects from the materials we measured HU from make_materials()
                    ( we may have iserts in the composition tables but we did not measure their HUs.)
        -------------------------------output-------------------------------------
        ct_hu : a dictionary of estimated HUs
    """
    water = Material(ra.water_den, ra.water_comp)
    ct_hu = {keys:0 for keys in materials.keys()}
    for mat in materials.keys():
        Za = materials[mat].ean_za
        Zc = materials[mat].ean_zc

        ctno = 1000*(materials[mat].edens / water.edens) *((Kph*Za**3.62 + Kcoh*Zc**1.86 + KN) / (Kph*(water.ean_za)**3.62 + Kcoh*(water.ean_zc)**1.86 + KN))
        ct_hu[mat] = ctno - 1000

    return ct_hu
