import math
from concurrent import futures
import numpy as np
from System.MSData import CFileMS2, CDataPack, CLocalIndex
from System.MSLogging import log_get_error
from Operator.MSOperatorGetter import operatorGetSequence
from Operator.MSOperatorModi import operatorDFSVariableSites, operatorGenerateVariablePeptides
from Operator.MSOperatorSetter import operatorInitialB7UserModification
from Operator.MSOperatorMass import operatorCalculateSingleChargeMass
from Tool.MSTool import toolInsertDictionary, toolGoldelValue


class CFunctionMS2Processing:

    def convertAndChooseTop(self, ms2_object=CFileMS2()):
        self.__convertMozIntensityToNumpy(ms2_object)
        self.__chooseTop(ms2_object)

    def assumeCharge(self, ms2_object=CFileMS2(), data_package=CDataPack()):
        self.__assumeCharge(ms2_object, data_package)

    @staticmethod
    def assumeCXPeaks(ms2_object=CFileMS2(), data_package=CDataPack()):
        mass_h2o = data_package.my_ini.MASS_H2O_MONO
        mass_proton = data_package.my_ini.MASS_PROTON_MONO
        for x_index in ms2_object.INDEX_SCAN:
            for y_index, single_charge_p_mass in enumerate(
                    operatorCalculateSingleChargeMass(np.array(ms2_object.MATRIX_MZ[x_index]),
                                                      np.array(ms2_object.MATRIX_CHARGE[x_index]))):
                # 假定峰为ｂ离子
                assume_b = ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] - mass_proton
                # 假定峰为y离子
                assume_y = ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] - mass_proton - mass_h2o
                # 假定峰为单碎裂交联ｂ离子
                assume_cx_b = ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] * (
                    -1) + single_charge_p_mass - mass_h2o + 1000
                # 假定峰为单碎裂交联y离子
                assume_cx_y = ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] * (-1) + single_charge_p_mass + 1000

                ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] = np.array(
                    [assume_b, assume_y, assume_cx_b, assume_cx_y]).flatten()
                ms2_object.MATRIX_PEAK_INT[x_index][y_index] = np.array(
                    [ms2_object.MATRIX_PEAK_INT[x_index][y_index]] * 2 +
                    [data_package.my_config.C3H_PREPROCESSING_ASSUME_CX_ION_INTENSITY_RATE *
                     ms2_object.MATRIX_PEAK_INT[x_index][
                         y_index]] * 2).flatten()

    @staticmethod
    def sortByMz(ms2_object=CFileMS2(), data_package=CDataPack()):
        for x_index in ms2_object.INDEX_SCAN:
            for y_index in range(len(ms2_object.MATRIX_PEAK_MOZ[x_index])):
                mz_list = ms2_object.MATRIX_PEAK_MOZ[x_index][y_index]
                intensity_list = ms2_object.MATRIX_PEAK_INT[x_index][y_index]
                satisfied_index = np.where((mz_list > data_package.my_config.C34_MS2_MOZ_LOWER_LIMIT) & (
                        mz_list < data_package.my_config.C34_MS2_MOZ_UPPER_LIMIT))
                args = np.argsort(mz_list[satisfied_index])
                ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] = mz_list[satisfied_index][args]
                ms2_object.MATRIX_PEAK_INT[x_index][y_index] = intensity_list[satisfied_index][args]

    @staticmethod
    def convertInt(ms2_object=CFileMS2()):
        for x_index in ms2_object.INDEX_SCAN:
            for y_index in range(len(ms2_object.MATRIX_PEAK_MOZ[x_index])):
                ms2_object.MATRIX_PEAK_MOZ[x_index][y_index] = np.array(
                    1e4 * ms2_object.MATRIX_PEAK_MOZ[x_index][y_index], dtype='u4')

    @staticmethod
    def __chooseTop(ms2_object):
        """
        选择并保留CFileMS2对象中Top-200的峰
        :param ms2_object:
        :return:
        """
        for ms2_index in ms2_object.INDEX_SCAN:
            order = np.lexsort((ms2_object.MATRIX_PEAK_MOZ[ms2_index], ms2_object.MATRIX_PEAK_INT[ms2_index]))[::-1][
                    :200]
            ms2_object.MATRIX_PEAK_MOZ[ms2_index] = ms2_object.MATRIX_PEAK_MOZ[ms2_index][order]
            ms2_object.MATRIX_PEAK_INT[ms2_index] = ms2_object.MATRIX_PEAK_INT[ms2_index][order]

    @staticmethod
    def __convertMozIntensityToNumpy(ms2_object):
        """
        把CFileMS2对象中的质荷比信息和强度信息转化为numpy数组
        :param ms2_object:CFileMS2对象
        :return: None
        """
        for ms2_index in ms2_object.INDEX_SCAN:
            ms2_object.MATRIX_PEAK_MOZ[ms2_index] = np.array(ms2_object.MATRIX_PEAK_MOZ[ms2_index])
            ms2_object.MATRIX_PEAK_INT[ms2_index] = np.array(ms2_object.MATRIX_PEAK_INT[ms2_index])

    @staticmethod
    def __assumeCharge(ms2_object, data_package):
        """
        假定每个峰为不超过最大母离子电荷的可能值，计算该峰在这个电荷下的单电荷质荷比
        :param ms2_object:
        :param data_package:
        :return:
        """
        for ms2_index in ms2_object.INDEX_SCAN:
            mz_result = []
            intensity_result = []
            if data_package.my_config.C3H_PREPROCESSING_PEAK_HIGHER_CHARGE:
                for max_charge in ms2_object.MATRIX_CHARGE[ms2_index]:
                    mz_result_one_precursor = []
                    intensity_result_one_precursor = []
                    for charge in range(1, max_charge + 1):
                        mz_result_one_precursor.append(ms2_object.MATRIX_PEAK_MOZ[ms2_index] * charge - (
                                charge - 1) * data_package.my_ini.MASS_PROTON_MONO)
                        intensity_result_one_precursor.append(ms2_object.MATRIX_PEAK_INT[ms2_index])
                    mz_result.append(np.array(mz_result_one_precursor).flatten())
                    intensity_result.append(np.array(intensity_result_one_precursor).flatten())
            else:
                for _ in ms2_object.MATRIX_CHARGE[ms2_index]:
                    mz_result_one_precursor = [ms2_object.MATRIX_PEAK_MOZ[ms2_index]]
                    intensity_result_one_precursor = [ms2_object.MATRIX_PEAK_INT[ms2_index]]
                    mz_result.append(np.array(mz_result_one_precursor).flatten())
                    intensity_result.append(np.array(intensity_result_one_precursor).flatten())
            ms2_object.MATRIX_PEAK_MOZ[ms2_index] = mz_result
            ms2_object.MATRIX_PEAK_INT[ms2_index] = intensity_result
