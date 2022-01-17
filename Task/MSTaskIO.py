import os
from System.MSData import CDataPack, CLocalIndex
from System.MSLogging import create_logger
from System.MSSystem import MS2_OUTPUT_TYPE,PSM_OUTPUT_TYPE
from Function.MSFunctionIO import CFunctionFastaFileParser, CFunctionIniParser, CFunctionMS2FileParser, CFunctionPSMIO
from Function.MSFunctionPreprocessing import CFunctionMS2Processing, CFunctionModification
from Function.MSFunctionQC import CFunctionBeforeOutput
from Operator.MSOperatorMass import operatorCalculateAAMass


class CTaskPreparing:
    def work(self, data_package=CDataPack()):
        file_name = self.__clearLogFile(data_package)
        self.__createLocalLogger(file_name, data_package)

    def __clearLogFile(self, data_package=CDataPack()):
        file_name = os.path.join(data_package.my_config.E1_PATH_EXPORT, 'SpotLink.log')
        data_package.my_config.E6_LOCAL_LOG_FILE_PATH = file_name
        with open(file_name, 'w') as f:
            f.write('\n')
        return file_name

    def __createLocalLogger(self, file_name, data_package=CDataPack()):
        data_package.my_config.E5_LOCAL_LOGGER = create_logger('local_logger', file_name)


class CTaskPrepareIni:

    def work(self, data_package=CDataPack()):
        self.__readIniFile(data_package)
        self.__aaMassPrepare(data_package)

    @staticmethod
    def __readIniFile(data_package):

        ini_class = CFunctionIniParser()
        ini_class.file_to_ini(data_package)

    @staticmethod
    def __aaMassPrepare(data_package):
        operatorCalculateAAMass(data_package)


class CTaskPrePareMS2:

    def work(self, data_package=CDataPack()):
        work_object = CFunctionMS2FileParser()
        ms2_list = self.__convertAllMS2(work_object, data_package)
        data_package.my_config.A4_PATH_MS2_PKL = self.__saveAllMS2(work_object, ms2_list, data_package, save_type=MS2_OUTPUT_TYPE['RAW'])
        self.__preprocessMS2List(ms2_list, data_package)
        data_package.my_config.A5_PATH_PROCESSED_MS2_PKL = self.__saveAllMS2(work_object, ms2_list, data_package,
                                                                             save_type=MS2_OUTPUT_TYPE['MS2_PROCESSED'])

    @staticmethod
    def __convertAllMS2(work_object, data_package):
        ms2_list = work_object.batch_convert_to_ms2(data_package.my_config.A1_PATH_MS2,data_package)
        return ms2_list

    @staticmethod
    def __preprocessMS2List(ms2_list, data_package):
        process = CFunctionMS2Processing()
        for ms2_object in ms2_list:
            process.convertAndChooseTop(ms2_object)
            process.assumeCharge(ms2_object, data_package)
            process.assumeCXPeaks(ms2_object, data_package)
            process.sortByMz(ms2_object, data_package)
            process.convertInt(ms2_object)

    @staticmethod
    def __saveAllMS2(work_object, ms2_list, data_package, save_type):
        return work_object.batch_save_to_drive(data_package.my_config.A1_PATH_MS2, ms2_list, save_type)


class CTaskProteinIndex:

    @staticmethod
    def work(data_package=CDataPack(), indexes=CLocalIndex()):
        fasta_data = CFunctionFastaFileParser()
        fasta_data.read_fasta(data_package, indexes)


class CTaskOrganizeAndOutputCx:

    @staticmethod
    def psm_work(data_package=CDataPack()):
        psm_organize = CFunctionPSMIO()
        before_output = CFunctionBeforeOutput()
        modification_output = CFunctionModification()
        for raw_psm_path, raw_mgf_path in zip(data_package.my_config.C36_CX_SVM_PSM_PKL,
                                              data_package.my_config.A1_PATH_MS2):
            psm_object = psm_organize.loadPSMPickle(raw_psm_path)
            mass_name_dic = modification_output.generateMassNameIndex(data_package)
            psm_organize.writePSMToCsv(raw_psm_path, mass_name_dic, psm_object, filtered=PSM_OUTPUT_TYPE['TOP1_UNFILTERED'])
            before_output.filteredByqValueCx(psm_object, data_package)
            psm_organize.writePSMToCsv(raw_psm_path, mass_name_dic, psm_object, filtered=PSM_OUTPUT_TYPE['TOP1_FILTERED'])

    @staticmethod
    def site_work(data_package=CDataPack()):
        psm_organize = CFunctionPSMIO()
        before_output = CFunctionBeforeOutput()
        modification_output = CFunctionModification()
        for raw_psm_path, raw_mgf_path in zip(data_package.my_config.C36_SITE_FDR_PKL,
                                              data_package.my_config.A1_PATH_MS2):
            psm_object = psm_organize.loadPSMPickle(raw_psm_path)
            mass_name_dic = modification_output.generateMassNameIndex(data_package)
            before_output.filteredBySiteqValueCx(psm_object, data_package)
            psm_organize.writeSiteToCsv(raw_psm_path, mass_name_dic, psm_object,data_package)
            psm_organize.writeProxl(raw_psm_path, psm_object, data_package)
