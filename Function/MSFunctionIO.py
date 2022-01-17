import csv
import os.path
import pickle
import numpy as np
import random

from System.MSData import CDataPack, CFileMS2, CLocalIndex, CDataPSMCX
from System.MSLogging import log_get_error

from Operator.MSOperatorIO import operatorLoadMS2, operatorSavePSM, operatorLoadPSM,operatorProxlMxlGenerator
from Operator.MSOperatorSetter import operatorInitCfileMS2

from Tool.MSTool import toolGetWord, toolStr2List



class CFunctionIniParser:

    def file_to_ini(self, data_package=CDataPack()):
        # 读取元素的ini文件，氨基酸的ini文件，修饰的ini文件，交联剂的ini文件，酶的ini文件
        self.__captainFileToElement(data_package)
        self.__captainFileToAA(data_package)
        self.__captainFile2Mod(data_package)
        self.__captainFile2Linker(data_package)
        self.__captainFile2Enzyme(data_package)

    @staticmethod
    def __captainFileToElement(data_package):
        try:
            with open(data_package.my_config.I0_INI_PATH_ELEMENT) as element_file:
                for item in element_file.readlines():
                    if item.strip():
                        str_name = toolGetWord(item, 0, '|')  # 元素名称
                        str_mass = toolGetWord(item, 1, '|')  # 元素质量
                        str_abdc = toolGetWord(item, 2, '|')  # 元素丰度

                        list_mass = toolStr2List(str_mass, ',')
                        list_abdc = toolStr2List(str_abdc, ',')

                        data_package.my_ini.DICT0_ELEMENT_MASS[str_name] = list_mass
                        data_package.my_ini.DICT0_ELEMENT_ABDC[str_name] = list_abdc
        except FileNotFoundError:
            log_get_error(f'Cannot find file {data_package.my_config.I0_INI_PATH_ELEMENT}',
                          data_package.my_config.E5_LOG_FILE_PATH)

    @staticmethod
    def __captainFileToAA(data_package):

        try:
            with open(data_package.my_config.I1_INI_PATH_AA) as aa_file:
                for item in aa_file.readlines():
                    if len(item) > 1:
                        str_name = toolGetWord(item, 0, '|')
                        str_comp = toolGetWord(item, 1, '|')
                        data_package.my_ini.DICT1_AA_COM[str_name] = str_comp
        except FileNotFoundError:
            log_get_error(f'Cannot find file {data_package.my_config.I1_INI_PATH_AA}',
                          data_package.my_config.E5_LOG_FILE_PATH)

    @staticmethod
    def __captainFile2Mod(data_package):
        try:
            with open(data_package.my_config.I2_INI_PATH_MOD) as modi_file:
                for item in modi_file.readlines():
                    if item[:4] != 'name' and item.strip() and item[0] != '@' and '=' in item:
                        modification_name, modification_description = item.strip().split('=')
                        modification_description = modification_description.split(' ')
                        modification_aa = modification_description[0]
                        modification_type = modification_description[1]
                        modification_mass = float(modification_description[2])
                        data_package.my_ini.DICT2_MOD_COM[modification_name] = (
                            modification_type, modification_aa, modification_mass)
        except FileNotFoundError:
            log_get_error(f'Cannot find file {data_package.my_config.I2_INI_PATH_MOD}',
                          data_package.my_config.E5_LOG_FILE_PATH)

    @staticmethod
    def __captainFile2Enzyme(data_package):
        try:
            with open(data_package.my_config.I5_INI_PATH_ENZYME) as enzyme_file:
                for item in enzyme_file.readlines():
                    if item[0] == '@' or item.strip() == '':
                        continue
                    enzyme_name = item.split()[0].split('=')[1]
                    if enzyme_name == data_package.my_config.B1_NAME_ENZYME:
                        data_package.my_ini.ENZYME_CUT_LIST = [aa for aa in item.split()[1] if aa != '_']
                        data_package.my_ini.ENZYME_CUT_IGNORE_LIST = [aa for aa in item.split()[2] if aa != '_']
                        data_package.my_ini.ENZYME_TERMINAL = item.strip().split()[3]
        except FileNotFoundError:
            log_get_error(f'Cannot find file {data_package.my_config.I5_INI_PATH_ENZYME}',
                          data_package.my_config.E5_LOG_FILE_PATH)

    @staticmethod
    def __captainFile2Linker(data_package):
        with open(data_package.my_config.I4_INI_PATH_LINKER) as linker_file:
            for item in linker_file.readlines():
                if item.strip():
                    linker_name = toolGetWord(item, 0, '=')
                    other_info = toolGetWord(item, 1, '=')
                    alpha_sites, beta_sites = toolGetWord(other_info, -1, None), toolGetWord(other_info, 0, None)
                    cx_mass = toolGetWord(other_info, 1, None)
                    alpha_mass, beta_mass = toolGetWord(other_info, 2, None), toolGetWord(other_info, 3, None)
                    data_package.my_ini.DICT4_LINKER_COM[linker_name] = (
                        cx_mass, alpha_sites, beta_sites, alpha_mass, beta_mass)


class CFunctionConfigParser:
    """
    从文件中读取config文件所存储的参数
    """

    @staticmethod
    def file_to_config(file_path, data_package_my_config):
        with open(file_path, 'r', encoding='utf8') as f:
            for line in f:
                line = line.strip()
                if '=' not in line:
                    continue
                if '#' in line:
                    line = line[:line.find('#')].strip()
                tmp = line.split('=')
                if len(tmp) < 2:
                    continue
                name = tmp[0].strip()
                value = tmp[1]
                if name == 'C_TYPE_SEARCH':
                    data_package_my_config.C_TYPE_SEARCH = int(value.split(';')[0])
                elif name == 'C_REBUILD_INDEX':
                    data_package_my_config.C_REBUILD_INDEX = int(value.split(';')[0])
                elif name == 'F1_THRESHOLD_FOR_PSM':
                    data_package_my_config.F1_THRESHOLD_FOR_PSM = float(value.split(';')[0])*0.01
                elif name == 'F3_THRESHOLD_FOR_SITE':
                    data_package_my_config.F3_THRESHOLD_FOR_SITE = float(value.split(';')[0])*0.01
                elif name == 'A1_PATH_MS2':
                    data_package_my_config.A1_PATH_MS2 = [eval("r'%s'" % item.replace(' ', '')) for item in
                                                          value.split(';') if item.strip() != '']
                elif name == "I0_INI_PATH_ELEMENT":
                    data_package_my_config.I0_INI_PATH_ELEMENT = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == 'I4_INI_PATH_LINKER':
                    data_package_my_config.I4_INI_PATH_LINKER = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == 'I5_INI_PATH_ENZYME':
                    data_package_my_config.I5_INI_PATH_ENZYME = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "I1_INI_PATH_AA":
                    data_package_my_config.I1_INI_PATH_AA = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "I2_INI_PATH_MOD":
                    data_package_my_config.I2_INI_PATH_MOD = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "A3_PATH_FASTA":
                    data_package_my_config.A3_PATH_FASTA = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "E1_PATH_EXPORT":
                    data_package_my_config.E1_PATH_EXPORT = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                    if not os.path.exists(data_package_my_config.E1_PATH_EXPORT):
                        os.makedirs(data_package_my_config.E1_PATH_EXPORT)
                elif name == "B1_NAME_ENZYME":
                    data_package_my_config.B1_NAME_ENZYME = value.split(';')[0].replace(' ', '')
                elif name == "C31_LINKER_NAME":
                    data_package_my_config.C31_LINKER_NAME = value.split(';')[0].replace(' ', '')
                elif name == "B2_TYPE_DIGEST":
                    data_package_my_config.B2_TYPE_DIGEST = int(value.split(';')[0])
                elif name == "B3_NUMBER_MAX_MISS_CLV":
                    data_package_my_config.B3_NUMBER_MAX_MISS_CLV = int(value.split(';')[0])
                elif name == "P5_MASS_PEP_LOW":
                    data_package_my_config.P5_MASS_PEP_LOW = float(value.split(';')[0])
                elif name == "P6_MASS_PEP_UP":
                    data_package_my_config.P6_MASS_PEP_UP = float(value.split(';')[0])
                elif name == "P7_LEN_PEP_LOW":
                    data_package_my_config.P7_LEN_PEP_LOW = int(value.split(';')[0])
                elif name == "P8_LEN_PEP_UP":
                    data_package_my_config.P8_LEN_PEP_UP = int(value.split(';')[0])
                elif name == "P1_NUMBER_THREAD":
                    data_package_my_config.P1_NUMBER_THREAD = int(value.split(';')[0])
                elif name == "B4_NAME_MOD_FIX":
                    data_package_my_config.B4_NAME_MOD_FIX = [item.replace(' ', '') for item in value.split(';') if item.strip()]
                elif name=="B5_NAME_MOD_VAR":
                    data_package_my_config.B5_NAME_MOD_VAR = [item.replace(' ', '') for item in value.split(';') if item.strip()]
                elif name=="B6_NUMBER_MAX_MOD":
                    data_package_my_config.B6_NUMBER_MAX_MOD=int(value.split(';')[0])
                elif name == "C33_MS1_TOLERANCE_TYPE":
                    data_package_my_config.C33_MS1_TOLERANCE_TYPE = int(value.split(';')[0])
                elif name == "C33_MS1_ABSOLUTE_TOLERANCE":
                    if data_package_my_config.C33_MS1_TOLERANCE_TYPE == 1:
                        data_package_my_config.C33_MS1_ABSOLUTE_TOLERANCE = float(value.split(';')[0]) * 1e-6
                    elif data_package_my_config.C33_MS1_TOLERANCE_TYPE == 0:
                        data_package_my_config.C33_MS1_ABSOLUTE_TOLERANCE = float(value.split(';')[0])
                elif name == "C33_MS1_RELATIVE_TOLERANCE":
                    if data_package_my_config.C33_MS1_TOLERANCE_TYPE == 1:
                        data_package_my_config.C33_MS1_RELATIVE_TOLERANCE = float(value.split(';')[0]) * 1e-6
                    elif data_package_my_config.C33_MS1_TOLERANCE_TYPE == 0:
                        data_package_my_config.C33_MS1_RELATIVE_TOLERANCE = float(value.split(';')[0])
                elif name == "C32_MS2_TOLERANCE_TYPE":
                    data_package_my_config.C32_MS2_TOLERANCE_TYPE = int(value.split(';')[0])
                elif name == "C32_MS2_ABSOLUTE_TOLERANCE":
                    if data_package_my_config.C32_MS2_TOLERANCE_TYPE == 1:
                        data_package_my_config.C32_MS2_ABSOLUTE_TOLERANCE = float(value.split(';')[0]) * 1e-6
                    elif data_package_my_config.C32_MS2_TOLERANCE_TYPE == 0:
                        data_package_my_config.C32_MS2_ABSOLUTE_TOLERANCE = float(value.split(';')[0])
                elif name == "C32_MS2_RELATIVE_TOLERANCE":
                    if data_package_my_config.C32_MS2_TOLERANCE_TYPE == 1:
                        data_package_my_config.C32_MS2_RELATIVE_TOLERANCE = float(value.split(';')[0]) * 1e-6
                    elif data_package_my_config.C32_MS2_TOLERANCE_TYPE == 0:
                        data_package_my_config.C32_MS2_RELATIVE_TOLERANCE = float(value.split(';')[0])
                elif name=="A4_PATH_MS2_PKL":
                    data_package_my_config.A4_PATH_MS2_PKL = [eval("r'%s'" % value.replace(' ', '')) for item in
                                                          value.split(';') if item.strip() != '']
                elif name=="C36_RAW_PSM_PKL":
                    data_package_my_config.C36_RAW_PSM_PKL=[eval("r'%s'" % value.replace(' ', '')) for item in
                                                          value.split(';') if item.strip() != '']
                elif name=="C36_SITE_SCORING_PSM_PKL":
                    data_package_my_config.C36_SITE_SCORING_PSM_PKL=[eval("r'%s'" % value.replace(' ', '')) for item in
                                                          value.split(';') if item.strip() != '']
                elif name=="C36_FINE_SCORING_PSM_PKL":
                    data_package_my_config.C36_FINE_SCORING_PSM_PKL=[eval("r'%s'" % value.replace(' ', '')) for item in
                                                          value.split(';') if item.strip() != '']


'''
与序列，质谱文件相关的函数
'''


class CFunctionFastaFileParser:
    """
    读取fasta文件至内存以进行后续处理
    """

    def read_fasta(self, data_package=CDataPack(), indexes=CLocalIndex()):
        """
        读取fasta序列，并将其中信息保存到indexes中，同时生成反库信息
        :param data_package:
        :param indexes:
        :return:
        """
        protein_index = self.__proteinIndexWithDescription(data_package)
        if protein_index:
            length_of_index = len(protein_index)
            for index, protein_name in enumerate(protein_index):
                indexes.CFLOW3_PROTEIN_INDEX[index] = protein_name
                indexes.CFLOW3_SEQUENCE_INDEX[index] = protein_index[protein_name][0]
                indexes.CFLOW3_DESCRIPTION_INDEX[index] = protein_index[protein_name][1]
                if data_package.my_config.F2_ENABLE_TARGET_DECOY:
                    indexes.CFLOW3_SEQUENCE_INDEX[index + length_of_index] = protein_index[protein_name][0][::-1]
                    indexes.CFLOW3_PROTEIN_INDEX[index + length_of_index] = 'DECOY_' + protein_name
                    indexes.CFLOW3_DESCRIPTION_INDEX[index + length_of_index] = ''
        else:
            log_get_error('Cannot build protein index.',data_package.my_config.E5_LOG_FILE_PATH)

    @staticmethod
    def __proteinIndexWithDescription(data_package):
        """
        将fasta文件一次性的读入到内存当中
        :return: 一个字典，字典的键为蛋白的名称，字典的值为（蛋白所对应的序列，蛋白对应的描述）
        """
        with open(data_package.my_config.A3_PATH_FASTA, 'r') as file_to_read:
            protein_index_dictionary = {}
            temp_sequence = ''
            temp_protein_name = ''
            temp_protein_description = ''
            for item in file_to_read.readlines():
                if item.strip():
                    if '>' in item:
                        protein_index_dictionary[temp_protein_name] = (temp_sequence, temp_protein_description)
                        temp_protein_name = item.split('|')[1]
                        if item.count('|') >= 2:
                            temp_protein_description = item.split('|')[2]
                        else:
                            temp_protein_description = ''
                        temp_sequence = ''
                    else:
                        temp_sequence += item.strip()
        protein_index_dictionary[temp_protein_name] = (temp_sequence, temp_protein_description)
        del protein_index_dictionary['']
        return protein_index_dictionary


class CFunctionMS2FileParser:
    """
    对于ms2文件的操作
    """

    def batch_convert_to_ms2(self, file_path_list,data_package):
        """
        读取一个列表的二级谱文件并以CFileMS2格式保存在内存中
        :param file_path_list: 二级谱文件的路径列表
        :return: CFileMS2对象列表
        """
        result_list = []
        for path in file_path_list:
            extension_name = os.path.splitext(path)[1]
            if extension_name == '.ms2' or extension_name == '.MS2':
                result_list.append(self.__readMs2FromMs2(path))
            elif extension_name == '.mgf' or extension_name == '.MGF':
                result_list.append(self.__readMs2FromMgf(path))
            else:
                log_get_error('Unknown MS2 format',data_package.my_config.E5_LOG_FILE_PATH)
        return result_list

    def batch_save_to_drive(self, file_path_list, ms2_object_list, save_type):
        """
        将多个二级谱对象保存到硬盘中，与原始二级谱处于同一目录下，文件名相同
        :param file_path_list: 二级谱文件路径列表
        :param ms2_object_list: CFileMS2对象列表
        :param save_type: 保存类型，如果为0代表保存的是原始谱图信息，如果为1代表的是预处理后的谱图信息
        :return: pkl文件名列表
        """
        return [self.__saveDrive(path, ms2, save_type) for path, ms2 in zip(file_path_list, ms2_object_list)]

    @staticmethod
    def loadMS2Pickle(file_path):
        """
        将保存于硬盘中的ms2文件读取到内存中
        :param file_path: 二级谱文件
        :return: CFileMS2对象
        """
        return operatorLoadMS2(file_path)

    @staticmethod
    def __readMs2FromMs2(file_path):
        """
        读取MS2文件，本函数继承自示例代码
        :param file_path:文件路径
        :return:CFileMS2对象
        """
        data_ms2 = CFileMS2()
        operatorInitCfileMS2(data_ms2)
        try:
            with open(file_path, 'r') as f:
                for line in f.readlines():
                    len_line = len(line)
                    if len_line > 1:
                        if line.startswith('S'):
                            temp_scan = int(line.split()[1])
                            data_ms2.MATRIX_PEAK_MOZ[temp_scan] = []
                            data_ms2.MATRIX_PEAK_INT[temp_scan] = []
                            if temp_scan not in data_ms2.INDEX_SCAN:
                                data_ms2.MATRIX_MZ[temp_scan] = []
                                data_ms2.MATRIX_CHARGE[temp_scan] = []
                                data_ms2.INDEX_SCAN.append(temp_scan)
                        elif line.startswith('H'):
                            continue
                        elif line.startswith('I'):
                            if 'PrecursorScan' in line:
                                precursor_scan = int(line.split()[2])
                                data_ms2.LIST_PRECURSOR_SCAN[temp_scan] = precursor_scan
                            elif 'IonInjectionTime' in line:
                                time = line.split()[2]
                                data_ms2.LIST_ION_INJECTION_TIME[temp_scan] = float(time)
                            elif 'RetTime' in line:
                                r_time = line.split()[2]
                                data_ms2.LIST_RET_TIME[temp_scan] = float(r_time)
                            elif 'ActivationCenter' in line:
                                center = line.split()[2]
                                data_ms2.LIST_ACTIVATION_CENTER[temp_scan] = float(center)
                            elif 'MonoiosotopicMz' in line or 'MonoisotopicMz' in line:
                                mz = line.split()[2]
                                data_ms2.MATRIX_MZ[temp_scan].append(float(mz))
                        elif line.startswith('Z'):
                            charge = line.split()[1]
                            data_ms2.MATRIX_CHARGE[temp_scan].append(int(charge))
                        else:
                            data_ms2.MATRIX_PEAK_MOZ[temp_scan].append(float(toolGetWord(line, 0, ' ')))
                            data_ms2.MATRIX_PEAK_INT[temp_scan].append(float(toolGetWord(line, 1, ' ')))
        except:
            log_get_error('[SpotLink] Failed to parse MS2 file {}!'.format(file_path),data_package.my_config.E5_LOG_FILE_PATH)
        return data_ms2

    @staticmethod
    def __readMs2FromMgf(file_path):
        """
        读取mgf文件
        :param file_path: 文件路径
        :return: CFileMS2对象
        """
        data_ms2 = CFileMS2()
        operatorInitCfileMS2(data_ms2)
        spectra_order=0
        with open(file_path, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if line:
                    if line.startswith('TITLE'):
                        if '.dta' in line:
                            scan_number = int(line.split('.')[1])
                        else:
                            spectra_order+=1
                            scan_number=spectra_order
                        data_ms2.MATRIX_PEAK_MOZ[scan_number] = []
                        data_ms2.MATRIX_PEAK_INT[scan_number] = []
                        if scan_number not in data_ms2.INDEX_SCAN:
                            data_ms2.MATRIX_MZ[scan_number] = []
                            data_ms2.MATRIX_CHARGE[scan_number] = []
                            data_ms2.INDEX_SCAN.append(scan_number)
                    elif line.startswith('CHARGE'):
                        data_ms2.MATRIX_CHARGE[scan_number].append(int(line.split('=')[1][:-1]))
                    elif line.startswith('RTINSECONDS'):
                        if '-' in line:
                            data_ms2.LIST_RET_TIME[scan_number] = float(line.split('=')[1].split('-')[0])
                        else:
                            data_ms2.LIST_RET_TIME[scan_number] = float(line.split('=')[1])
                    elif line.startswith('PEPMASS'):
                        if ' ' in line:
                            data_ms2.MATRIX_MZ[scan_number].append(float(line.split('=')[1].split(' ')[0]))
                        else:
                            data_ms2.MATRIX_MZ[scan_number].append(float(line.split('=')[1]))
                    elif line[0].isdigit():
                        data_ms2.MATRIX_PEAK_MOZ[scan_number].append(float(line.split()[0]))
                        data_ms2.MATRIX_PEAK_INT[scan_number].append(float(line.split()[1]))
        return data_ms2

    @staticmethod
    def __saveDrive(file_path, ms2_object, save_type):
        """
        将二级谱对象保存至硬盘中
        :param file_path: 二级谱文件的路径
        :param ms2_object: CFileMS2对象
        :param save_type: 保存类型，如果为0代表保存的是原始谱图信息，如果为1代表的是预处理后的谱图信息
        :return: None
        """
        file_dir_path, file_name = os.path.split(file_path)
        pure_file_name = os.path.splitext(file_name)[0]
        if save_type:
            output_name = os.path.normpath(os.path.join(file_dir_path, pure_file_name + '.sms2'))
        else:
            output_name = os.path.normpath(os.path.join(file_dir_path, pure_file_name + '.rms2'))
        output_pickle_file = open(output_name, 'wb')
        pickle.dump(ms2_object, output_pickle_file)
        output_pickle_file.close()
        return output_name


class CFunctionPSMIO:

    @staticmethod
    def loadPSMPickle(psm_path):
        psm_object = operatorLoadPSM(psm_path)
        return psm_object

    @staticmethod
    def saveBinary(ms2_path, psm_data=CDataPSMCX(), data_package=CDataPack(), extension_name='rpsm'):
        if extension_name == 'rpsm':
            data_package.my_config.C36_RAW_PSM_PKL.append(
                operatorSavePSM(ms2_path, psm_data, data_package, extension_name))
        elif extension_name == 'fpsm':
            data_package.my_config.C36_FINE_SCORING_PSM_PKL.append(
                operatorSavePSM(ms2_path, psm_data, data_package, extension_name))
        elif extension_name == 'spsm':
            data_package.my_config.C36_SITE_SCORING_PSM_PKL.append(
                operatorSavePSM(ms2_path, psm_data, data_package, extension_name))
        elif extension_name == 'cxsvmpsm':
            data_package.my_config.C36_CX_SVM_PSM_PKL.append(
                operatorSavePSM(ms2_path, psm_data, data_package, extension_name))
        elif extension_name=='sfdrpsm':
            data_package.my_config.C36_SITE_FDR_PKL.append(
                operatorSavePSM(ms2_path, psm_data, data_package, extension_name))

    def writeCoarsePSMToCsv(self, psm_path, psm_data=CDataPSMCX()):
        filename = os.path.splitext(psm_path)[0] + '_result_coarse_score.csv'
        with open(filename, 'wt', newline='') as output_file:
            output_file_writer = csv.writer(output_file)
            output_file_writer.writerow(
                ('Order', 'Spectra Title', 'Precursor Charge', 'Precursor Mz',
                 'Peptide', 'Peptide Mz', 'Alpha Peptide Protein', 'Beta Peptide Protein',
                 'Alpha Score', 'Beta Score', 'Coarse Score'))
            for counter, i in enumerate(psm_data.PSM_RANK_INDEX):
                index, spec_title = psm_data.PSM_SPECTRA_INDEX[i], psm_data.PSM_SPECTRA_TITLE[i]
                mass, charge, score = psm_data.PSM_MASS[i], psm_data.PSM_CHARGE[i], psm_data.PSM_SCORE[i]
                alpha_score, beta_score = psm_data.PSM_PEP1_COARSE_SCORE[i], psm_data.PSM_PEP2_COARSE_SCORE[i]
                pep_mz = psm_data.PSM_OBSERVE_MOZ[i]
                # 整合肽段信息
                pep1, pep2 = psm_data.PSM_PEP1[i], psm_data.PSM_PEP2[i]
                peptide_info = '{}()-{}()'.format(pep1, pep2)
                pep1_pro, pep2_pro = psm_data.PSM_PEP1_PROTEIN[i], psm_data.PSM_PEP2_PROTEIN[i]
                output_file_writer.writerow(
                    (counter + 1, spec_title, charge, mass, peptide_info, pep_mz, pep1_pro, pep2_pro,
                     alpha_score, beta_score, score))

    def writeFinePSMToCsv(self,psm_path,psm_data=CDataPSMCX()):
        rank_by_score_dictionary={}
        for item in psm_data.PSM_RANK_INDEX:
            if psm_data.PSM_SPECTRA_TITLE[item] in rank_by_score_dictionary:
                rank_by_score_dictionary[psm_data.PSM_SPECTRA_TITLE[item]].append(item)
            else:
                rank_by_score_dictionary[psm_data.PSM_SPECTRA_TITLE[item]]=[item]
        for item in rank_by_score_dictionary:
            raw_order=np.array(rank_by_score_dictionary[item])
            rank_by_score_dictionary[item]=raw_order[np.argsort(np.array(psm_data.PSM_SCORE)[raw_order])[::-1]]
        filename = os.path.splitext(psm_path)[0] + '_result_fine_score.csv'
        counter=1
        with open(filename,'wt',newline='') as output_file:
            output_file_writer=csv.writer(output_file)
            output_file_writer.writerow(
                ('Order', 'Spectra Title', 'Precursor Charge', 'Precursor Mz',
                 'Peptide', 'Peptide Mz', 'Alpha Peptide Protein', 'Beta Peptide Protein', 'Score'))
            for title in rank_by_score_dictionary:
                for item in rank_by_score_dictionary[title]:
                    mass, charge, score = psm_data.PSM_MASS[item], psm_data.PSM_CHARGE[item], psm_data.PSM_SCORE[item]
                    pep_mz = psm_data.PSM_OBSERVE_MOZ[item]
                    pep1, pep2 = psm_data.PSM_PEP1[item], psm_data.PSM_PEP2[item]
                    peptide_info = '{}()-{}()'.format(pep1, pep2)
                    pep1_pro, pep2_pro = psm_data.PSM_PEP1_PROTEIN[item], psm_data.PSM_PEP2_PROTEIN[item]
                    output_file_writer.writerow((counter, title, charge, mass, peptide_info, pep_mz,pep1_pro, pep2_pro,score))
                    counter+=1

    def writePSMToCsv(self, psm_path, mass_name_dic, psm_data=CDataPSMCX(), filtered=0):
        if filtered == 0:
            filename = os.path.splitext(psm_path)[0] + '_result_notop1.csv'
        elif filtered == 1:
            filename = os.path.splitext(psm_path)[0] + '_result.csv'
        elif filtered == 2:
            filename = os.path.splitext(psm_path)[0] + '_result_filtered.csv'
        try:
            with open(filename,'wt') as f:
                pass
        except PermissionError:
            filename+=f'.{int(1e8*random.random())}.csv'
        with open(filename, 'wt', newline='') as output_file:
            output_file_writer = csv.writer(output_file)
            output_file_writer.writerow(
                ('Order', 'Spectra Title', 'Precursor Charge', 'Precursor Mz', 'Type',
                 'Peptide','Peptide Mz', 'Modifications', 'Alpha Peptide Protein', 'Beta Peptide Protein', 'Dot Score','Score',
                 'Coarse Score','SVM Score', 'q Value'))
            if type(psm_data.PSM_RANK_INDEX[0]) == list:
                flatten_index = [i for item_list in psm_data.PSM_RANK_INDEX for i in item_list]
            else:
                flatten_index = psm_data.PSM_RANK_INDEX
            counter=0
            for i in flatten_index:
                index, spec_title = psm_data.PSM_SPECTRA_INDEX[i], psm_data.PSM_SPECTRA_TITLE[i]
                mass, charge, score,dot_score = psm_data.PSM_MASS[i], psm_data.PSM_CHARGE[i], psm_data.PSM_SCORE[i],psm_data.PSM_DOT_SCORE[i][0]
                pep_mz = psm_data.PSM_OBSERVE_MOZ[i]
                coarse_score=psm_data.PSM_PEP1_COARSE_SCORE[i]+psm_data.PSM_PEP2_COARSE_SCORE[i]
                pep1, pep2 = psm_data.PSM_PEP1[i], psm_data.PSM_PEP2[i]
                if len(psm_data.PSM_CX_SITES):
                    link_sites = psm_data.PSM_CX_SITES[i]
                peptide_info = '{}({})-{}({})'.format(pep1, int(link_sites[0][0]) + 1, pep2, int(link_sites[0][1]) + 1)
                pep1_modi, pep2_modi = psm_data.PSM_PEP1_MODIFICATION[i], psm_data.PSM_PEP2_MODIFICATION[i]
                modi_info = self.__organizeModificationInfo(pep1, pep1_modi, pep2_modi, mass_name_dic)

                pep1_pro, pep2_pro = psm_data.PSM_PEP1_PROTEIN[i], psm_data.PSM_PEP2_PROTEIN[i]
                pep1_pro_all_decoy=True
                pep2_pro_all_decoy=True
                for item in pep1_pro:
                    if 'DECOY' not in item:
                        pep1_pro_all_decoy=False
                        break
                for item in pep2_pro:
                    if 'DECOY' not in item:
                        pep2_pro_all_decoy=False
                        break
                if pep1_pro_all_decoy or pep2_pro_all_decoy:
                    continue
                svm_score = psm_data.PSM_SVM_SCORE[i]
                if len(psm_data.PSM_q_VALUE):
                    q_value = psm_data.PSM_q_VALUE[i]
                else:
                    q_value = 1
                output_file_writer.writerow(
                    (counter + 1, spec_title, charge, mass, 3, peptide_info,
                     pep_mz, modi_info, pep1_pro, pep2_pro,dot_score,
                     score,coarse_score, svm_score, q_value))
                counter+=1

    def writeSiteToCsv(self, psm_path, mass_name_dic, psm_data=CDataPSMCX(),data_package=CDataPack()):
        cx_dictionary={}
        for counter, j in enumerate(psm_data.PSM_RANK_INDEX):
            i, site = psm_data.PSM_SITE_q_VALUE_INFO[j]
            if psm_data.PSM_q_VALUE[i]>data_package.my_config.F1_THRESHOLD_FOR_PSM:
                continue
            spec_title = psm_data.PSM_SPECTRA_TITLE[i]
            pep1, pep2 = psm_data.PSM_PEP1[i], psm_data.PSM_PEP2[i]
            site_list = psm_data.PSM_CX_SITE_RELIABILITY[i]
            pep1_modi, pep2_modi = psm_data.PSM_PEP1_MODIFICATION[i], psm_data.PSM_PEP2_MODIFICATION[i]
            modi_info = self.__organizeModificationInfo(pep1, pep1_modi, pep2_modi, mass_name_dic)
            if type(site_list[0]) == dict:
                peptide_info = '{}({})-{}({})'.format(pep1, int(site) + 1, pep2, int(site_list[1]) + 1)
            else:
                peptide_info = '{}({})-{}({})'.format(pep1, int(site_list[0]) + 1, pep2, int(site) + 1)
            if peptide_info in cx_dictionary:
                cx_dictionary[peptide_info].append(spec_title)
            else:
                cx_dictionary[peptide_info]=[spec_title]
        cx_list=[]
        filename = os.path.splitext(psm_path)[0] + '_site_filtered.csv'
        try:
            with open(filename,'wt') as f:
                pass
        except PermissionError:
            filename+=f'.{int(1e8*random.random())}.csv'
        with open(filename, 'wt', newline='') as output_file:
            output_file_writer = csv.writer(output_file)
            output_file_writer.writerow(
                ('Order', 'Type','Peptide',
                 'Peptide Mz', 'Modifications', 'Alpha Peptide Protein', 'Beta Peptide Protein', 'Score','Site q Value','Support spectra number', 'Spectra Title'))
            counter=0
            for j in psm_data.PSM_RANK_INDEX:
                i,site=psm_data.PSM_SITE_q_VALUE_INFO[j]
                if psm_data.PSM_q_VALUE[i] > data_package.my_config.F1_THRESHOLD_FOR_PSM:
                    continue
                score = psm_data.PSM_SCORE[i]
                pep_mz = psm_data.PSM_OBSERVE_MOZ[i]
                # 整合肽段信息
                pep1, pep2 = psm_data.PSM_PEP1[i], psm_data.PSM_PEP2[i]
                pep1_pro, pep2_pro = psm_data.PSM_PEP1_PROTEIN[i], psm_data.PSM_PEP2_PROTEIN[i]
                if 'DECOY' in str(pep1_pro) or 'DECOY' in str(pep2_pro):
                    continue
                site_list=psm_data.PSM_CX_SITE_RELIABILITY[i]
                if type(site_list[0])==dict:
                    peptide_info = '{}({})-{}({})'.format(pep1, int(site)+1,pep2, int(site_list[1])+1)
                else:
                    peptide_info = '{}({})-{}({})'.format(pep1,int(site_list[0]) + 1, pep2, int(site) + 1 )
                # 整合修饰信息
                pep1_modi, pep2_modi = psm_data.PSM_PEP1_MODIFICATION[i], psm_data.PSM_PEP2_MODIFICATION[i]
                modi_info = self.__organizeModificationInfo(pep1, pep1_modi, pep2_modi, mass_name_dic)
                if peptide_info in cx_list:
                    continue
                else:
                    spec_title=';'.join(cx_dictionary[peptide_info])
                    cx_list.append(peptide_info)
                q_value = psm_data.PSM_SITE_q_VALUE[j]
                output_file_writer.writerow(
                    (counter + 1,  3, peptide_info,
                     pep_mz, modi_info, pep1_pro, pep2_pro,
                     score, q_value, len(cx_dictionary[peptide_info]),spec_title))
                counter+=1

    def writeProxl(self,psm_path,psm_data=CDataPSMCX(),data_package=CDataPack()):
        filename = os.path.splitext(psm_path)[0] + '_proxl.proxl.xml'
        try:
            with open(filename,'wt') as f:
                pass
        except PermissionError:
            filename+=f'.{int(1e8*random.random())}.xml'
        peptide_dic=self.__organizeIntoPeptide(psm_data)
        xml_string=operatorProxlMxlGenerator(psm_data,peptide_dic,data_package)
        with open(filename, 'w') as f:
            f.write(xml_string)


    @staticmethod
    def __organizeModificationInfo(pep1, pep1_modi, pep2_modi, mass_name_dic):
        modification_infor = 'pep1:'
        for i in pep1_modi:
            modi_mass = '{:.5f}'.format(pep1_modi[i])
            if modi_mass in mass_name_dic:
                modi_name = mass_name_dic[modi_mass]
            else:
                modi_name = modi_mass
            modification_infor += '{}({}),'.format(modi_name, int(i) + 1)
        if modification_infor[-1]==':':
            modification_infor+=';'
        elif modification_infor[-1]==',':
            modification_infor=modification_infor[:-1]+';'
        modification_infor += 'pep2:'
        for i in pep2_modi:
            modi_mass = '{:.5f}'.format(pep2_modi[i])
            if modi_mass in mass_name_dic:
                modi_name = mass_name_dic[modi_mass]
            else:
                modi_name = modi_mass
            modification_infor += '{}({}),'.format(modi_name, int(i) + 1)
        if modification_infor[-1]==':':
            modification_infor+=';'
        elif modification_infor[-1]==',':
            modification_infor=modification_infor[:-1]+';'
        return modification_infor

    @staticmethod
    def __organizeIntoPeptide(psm_data):
        peptide_dictionary={}
        for order in range(len(psm_data.PSM_q_VALUE)):
            if psm_data.PSM_DECOY[order] in ['TD','DD']:
                continue
            pep1=psm_data.PSM_PEP1[order]
            pep2=psm_data.PSM_PEP2[order]
            modi1=psm_data.PSM_PEP1_MODIFICATION[order]
            modi2=psm_data.PSM_PEP2_MODIFICATION[order]
            if type(psm_data.PSM_CX_SITE_RELIABILITY[order][0])==dict:
                b_s,b_s_s=0,-9999
                for item in psm_data.PSM_CX_SITE_RELIABILITY[order][0]:
                    if psm_data.PSM_CX_SITE_RELIABILITY[order][0][item]>b_s_s:
                        b_s,b_s_s=item,psm_data.PSM_CX_SITE_RELIABILITY[order][0][item]
                site1,site2=b_s,psm_data.PSM_CX_SITE_RELIABILITY[order][1]
            else:
                b_s, b_s_s = 0, -9999
                for item in psm_data.PSM_CX_SITE_RELIABILITY[order][1]:
                    if psm_data.PSM_CX_SITE_RELIABILITY[order][1][item] > b_s_s:
                        b_s, b_s_s = item, psm_data.PSM_CX_SITE_RELIABILITY[order][1][item]
                site1, site2 = psm_data.PSM_CX_SITE_RELIABILITY[order][0],b_s
            site1,site2=str(int(site1)+1),str(int(site2)+1)
            pep_str1=''
            for i,item in enumerate(pep1):
                if i in modi1:
                    pep_str1+='{}[{:.2f}]'.format(item,modi1[i])
                else:
                    pep_str1 += '{}'.format(item)
            pep_str1+='({})-'.format(site1)
            for i,item in enumerate(pep2):
                if i in modi2:
                    pep_str1+='{}[{:.2f}]'.format(item,modi2[i])
                else:
                    pep_str1+='{}'.format(item)
            pep_str1+='({})'.format(site2)
            pep_str2 = ''
            for i, item in enumerate(pep2):
                if i in modi2:
                    pep_str2 += '{}[{:.2f}]'.format(item, modi2[i])
                else:
                    pep_str2+='{}'.format(item)
            pep_str2 += '({})-'.format(site2)
            for i, item in enumerate(pep1):
                if i in modi1:
                    pep_str2 += '{}[{:.2f}]'.format(item, modi1[i])
                else:
                    pep_str2+=''
            pep_str2 += '({})'.format(site1)
            if pep_str1 not in peptide_dictionary and pep_str2 not in peptide_dictionary:
                peptide_dictionary[pep_str1]=[[order],pep1,pep2,site1,site2,modi1,modi2]
            elif pep_str1 in peptide_dictionary:
                peptide_dictionary[pep_str1][0].append(order)
            elif pep_str2 in peptide_dictionary:
                peptide_dictionary[pep_str2][0].append(order)
        return peptide_dictionary
