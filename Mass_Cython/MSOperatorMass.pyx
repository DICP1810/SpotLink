import numpy as np
cimport numpy as np
from System.MSData import CDataPack
from Tool.MSTool import toolDecompositionFormular, toolInsertDictionary


def operatorCalculateAAMass(data_package=CDataPack()):

    for amino_acid in data_package.my_ini.DICT1_AA_COM:
        formula = data_package.my_ini.DICT1_AA_COM[amino_acid]
        decomposition_formula = toolDecompositionFormular(formula)
        mass = 0
        for element, number in zip(decomposition_formula[0], decomposition_formula[1]):
            mass += (data_package.my_ini.DICT0_ELEMENT_MASS[element][0] * number)
        data_package.my_ini.DICT1_AA_MASS[amino_acid] = mass


def operatorSumUpPeptide(peptide_sequence, aa_mass_dic):

    cdef float total_mass = 0
    for item in peptide_sequence:
        total_mass += aa_mass_dic[item]
    return total_mass


cpdef np.ndarray operatorCalculateMZFromZeroCharge(np.ndarray mass_list, list charge_list,float proton_mass):

    cdef int length_of_mass=len(mass_list)
    cdef int length_of_charge=len(charge_list)
    cdef list temp_results
    if length_of_charge:
        # charge_list不为空
        temp_results = [(mass_list + charge * proton_mass) / charge for charge in
                        charge_list if charge!=0]
        if 0 in charge_list:
            temp_results.append(mass_list)
        return np.concatenate(temp_results)
    else:
        return mass_list

cpdef np.ndarray operatorCalculateSearchResidueFromZeroCharge(np.ndarray mass_list, list charge_list):

    cdef int length_of_mass=len(mass_list)
    cdef int length_of_charge=len(charge_list)
    cdef list temp_results
    if length_of_charge:
        # charge_list不为空
        temp_results = [(mass_list + charge * 1000) / charge for charge in
                        charge_list if charge!=0]
        if 0 in charge_list:
            temp_results.append(mass_list)
        return np.concatenate(temp_results)
    else:
        return mass_list

def operatorCalculateSingleChargeMass(mz_list, charge_list, data_package=CDataPack()):

    return charge_list * (mz_list - data_package.my_ini.MASS_PROTON_MONO) + data_package.my_ini.MASS_PROTON_MONO


cpdef np.ndarray operatorOnePepOneCleavageHCD(str peptide, dict modification_dictionary, float h2o_mono_mass, dict aa_mass_dictionary):

    cdef list results = []
    cdef int length = len(peptide)
    cdef float b_ion_mass = 0.0
    cdef float y_ion_mass = h2o_mono_mass
    cdef float b_modi_mass = 0.0
    cdef int b_modi_order = 0
    cdef float y_modi_mass = 0.0
    cdef int y_modi_order = len(modification_dictionary) - 1
    cdef int order

    cdef list modi_position
    modi_position = list(modification_dictionary.keys())
    modi_position.sort()

    for order in range(length):
        if b_modi_order < len(modi_position) and order >= modi_position[b_modi_order]:
            b_modi_mass += modification_dictionary[modi_position[b_modi_order]]
            b_modi_order += 1
        b_ion_mass += aa_mass_dictionary[peptide[order]]
        results.append(b_ion_mass + b_modi_mass)
        if y_modi_order >= 0 and length - 1 - order <= modi_position[y_modi_order]:
            y_modi_mass += modification_dictionary[modi_position[y_modi_order]]
            y_modi_order -= 1
        y_ion_mass += aa_mass_dictionary[peptide[length - 1 - order]]
        results.append(y_ion_mass + y_modi_mass)
    return np.array(results)

cpdef np.ndarray operatorOnePepOneCleavageHCDaion(str peptide, dict modification_dictionary, float h2o_mono_mass, dict aa_mass_dictionary):

    cdef list results = []
    cdef int length = len(peptide)
    cdef float a_ion_mass = -27.994915
    cdef float a_modi_mass = 0.0
    cdef int a_modi_order = 0
    cdef int order

    cdef list modi_position
    modi_position = list(modification_dictionary.keys())
    modi_position.sort()

    for order in range(length):
        # 计算b离子的质量
        if a_modi_order < len(modi_position) and order >= modi_position[a_modi_order]:
            a_modi_mass += modification_dictionary[modi_position[a_modi_order]]
            a_modi_order += 1
        a_ion_mass += aa_mass_dictionary[peptide[order]]
        results.append(a_ion_mass + a_modi_mass)
    return np.array(results)

cpdef np.ndarray operatorOnePepOneCleavageResidue(str peptide, dict modification_dictionary, dict aa_mass_dictionary):

    cdef list results = []
    cdef int length = len(peptide)
    cdef float b_ion_mass = 0.0
    cdef float y_ion_mass = 0.0
    cdef float b_modi_mass = 0.0
    cdef int b_modi_order = 0
    cdef float y_modi_mass = 0.0
    cdef int y_modi_order = len(modification_dictionary) - 1
    cdef int order

    cdef list modi_position
    modi_position = list(modification_dictionary.keys())
    modi_position.sort()

    for order in range(length):
        if b_modi_order < len(modi_position) and order >= modi_position[b_modi_order]:
            b_modi_mass += modification_dictionary[modi_position[b_modi_order]]
            b_modi_order += 1
        b_ion_mass += aa_mass_dictionary[peptide[order]]
        results.append(b_ion_mass + b_modi_mass)
        if y_modi_order >= 0 and length - 1 - order <= modi_position[y_modi_order]:
            y_modi_mass += modification_dictionary[modi_position[y_modi_order]]
            y_modi_order -= 1
        y_ion_mass += aa_mass_dictionary[peptide[length - 1 - order]]
        results.append(y_ion_mass + y_modi_mass)
    return np.array(results)


cpdef np.ndarray operatorSinglePeptideCalculatorUnsort(str peptide, list charge_list, dict modification_dictionary,
                                                       dict aa_mass_dictionary, float h2o_mono_mass,
                                                       float proton_mono_mass):
    cdef np.ndarray result
    result = operatorOnePepOneCleavageHCD(peptide, modification_dictionary,h2o_mono_mass,aa_mass_dictionary)
    result = operatorCalculateMZFromZeroCharge(result, charge_list, proton_mono_mass)
    return result

cpdef np.ndarray operatorSinglePeptideResidueCalculatorUnsort(str peptide, list charge_list,
                                                              dict modification_dictionary, dict aa_mass_dictionary,
                                                              float proton_mono_mass):

    cdef np.ndarray result
    result = operatorOnePepOneCleavageResidue(peptide, modification_dictionary,aa_mass_dictionary)
    result = operatorCalculateSearchResidueFromZeroCharge(result, charge_list)
    return result

def operatorCrossLinkingPeptideCalculator(alpha_peptide, beta_peptide, alpha_site, beta_site, linker_mass=0.0,
                                          charge_list=None, alpha_modification_dictionary=None,
                                          beta_modification_dictionary=None, h20_mono_mass= 18.010565,proton_mono_mass=1.00727645224,aa_mass_dic={}):

    linking_modification_for_alpha = operatorSumUpPeptide(
        beta_peptide, aa_mass_dic) + h20_mono_mass + linker_mass + sum(
        beta_modification_dictionary.values())
    linking_modification_for_beta = operatorSumUpPeptide(
        alpha_peptide, aa_mass_dic) + h20_mono_mass + linker_mass + sum(
        alpha_modification_dictionary.values())
    toolInsertDictionary(alpha_modification_dictionary, alpha_site,
                                                         linking_modification_for_alpha)
    toolInsertDictionary(beta_modification_dictionary, beta_site,
                                                        linking_modification_for_beta)
    alpha_result_nparray = operatorOnePepOneCleavageHCD(alpha_peptide, alpha_modification_dictionary, h20_mono_mass,aa_mass_dic)
    beta_result_nparray = operatorOnePepOneCleavageHCD(beta_peptide, beta_modification_dictionary, h20_mono_mass,aa_mass_dic)
    results = np.concatenate([alpha_result_nparray, beta_result_nparray])
    results = operatorCalculateMZFromZeroCharge(results, charge_list, proton_mono_mass)
    return np.unique(results)

def operatorCrossLinkingPeptideCalculatorMatrix(alpha_peptide, beta_peptide, alpha_site, beta_site, linker_mass=0.0,
                                          charge_list=None, alpha_modification_dictionary=None,
                                          beta_modification_dictionary=None, aa_mass_dic={},h2o_mass_mono=0,mass_proton_mono=0):

    linking_modification_for_alpha = operatorSumUpPeptide(
        beta_peptide, aa_mass_dic) + h2o_mass_mono + linker_mass + sum(
        beta_modification_dictionary.values())
    linking_modification_for_beta = operatorSumUpPeptide(
        alpha_peptide, aa_mass_dic) + h2o_mass_mono + linker_mass + sum(
        alpha_modification_dictionary.values())
    toolInsertDictionary(alpha_modification_dictionary, alpha_site,
                                                         linking_modification_for_alpha)
    toolInsertDictionary(beta_modification_dictionary, beta_site,
                                                        linking_modification_for_beta)
    alpha_result_nparray = operatorOnePepOneCleavageHCD(alpha_peptide, alpha_modification_dictionary,h2o_mass_mono,aa_mass_dic)
    alpha_peptide_zero_charge_b=np.array([0]+[item for item in alpha_result_nparray[::2]])
    alpha_peptide_zero_charge_y=np.array([item for item in alpha_result_nparray[1::2][::-1]]+[0])
    #alpha_result_a_ion_nparray=operatorOnePepOneCleavageHCDaion(alpha_peptide, alpha_modification_dictionary, data_package.my_ini.MASS_H2O_MONO,data_package.my_ini.DICT1_AA_MASS)
    #alpha_result_a_ion_nparray=np.array([0]+[item for item in alpha_result_a_ion_nparray])
    beta_result_nparray = operatorOnePepOneCleavageHCD(beta_peptide, beta_modification_dictionary, h2o_mass_mono,aa_mass_dic)
    beta_peptide_zero_charge_b = np.array([0]+[item for item in beta_result_nparray[::2]])
    beta_peptide_zero_charge_y = np.array([item for item in beta_result_nparray[1::2][::-1]]+[0])
    #beta_result_a_ion_nparray=operatorOnePepOneCleavageHCDaion(beta_peptide, beta_modification_dictionary, data_package.my_ini.MASS_H2O_MONO,data_package.my_ini.DICT1_AA_MASS)
    #beta_result_a_ion_nparray=np.array([0]+[item for item in beta_result_a_ion_nparray])
    alpha_peptide_one_charge_b = alpha_peptide_zero_charge_b + mass_proton_mono
    alpha_peptide_two_charge_b = (alpha_peptide_zero_charge_b + 2 * mass_proton_mono) / 2
    #alpha_peptide_three_charge_b = (alpha_peptide_zero_charge_b + 3 * data_package.my_ini.MASS_PROTON_MONO) / 3
    alpha_peptide_one_charge_y = alpha_peptide_zero_charge_y + mass_proton_mono
    alpha_peptide_two_charge_y = (alpha_peptide_zero_charge_y + 2 * mass_proton_mono) / 2
    #alpha_peptide_three_charge_y = (alpha_peptide_zero_charge_y + 3 * data_package.my_ini.MASS_PROTON_MONO) / 3
    #alpha_peptide_one_charge_a=alpha_result_a_ion_nparray+data_package.my_ini.MASS_PROTON_MONO
    #alpha_peptide_two_charge_a=(alpha_result_a_ion_nparray+2*data_package.my_ini.MASS_PROTON_MONO)/2

    beta_peptide_one_charge_b = beta_peptide_zero_charge_b + mass_proton_mono
    beta_peptide_two_charge_b = (beta_peptide_zero_charge_b + 2 * mass_proton_mono) / 2
    #beta_peptide_three_charge_b = (beta_peptide_zero_charge_b + 3 * data_package.my_ini.MASS_PROTON_MONO) / 3
    beta_peptide_one_charge_y = beta_peptide_zero_charge_y + mass_proton_mono
    beta_peptide_two_charge_y = (beta_peptide_zero_charge_y + 2 * mass_proton_mono) / 2
    #beta_peptide_three_charge_y = (beta_peptide_zero_charge_y + 3 * data_package.my_ini.MASS_PROTON_MONO) / 3
    #beta_peptide_one_charge_a=beta_result_a_ion_nparray+data_package.my_ini.MASS_PROTON_MONO
    #beta_peptide_two_charge_a=(beta_result_a_ion_nparray+2*data_package.my_ini.MASS_PROTON_MONO)/2

    return np.array([alpha_peptide_one_charge_b,alpha_peptide_two_charge_b,alpha_peptide_one_charge_y,alpha_peptide_two_charge_y]),\
           np.array([beta_peptide_one_charge_b,beta_peptide_two_charge_b,beta_peptide_one_charge_y,beta_peptide_two_charge_y])
