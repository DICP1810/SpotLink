from collections import deque

import numpy as np

from System.MSData import CLocalIndex, CSearcherStatic,CSearcherDynamic, CPeptideInfo, CFileMS2, CDataPSMCX
from System.MSSystem import VALUE_ILLEGAL, VALUE_MAX_SCAN
from Tool.MSTool import toolGetSortIndex

"""""""""
CConfig的setter方法
"""""""""


def operatorInitialB7UserModification(data_package):
    """
    初始化B7UserModification，用于存储用户定义的修饰
    :param data_package: data package
    :return: None
    """
    data_package.my_config.B7_USER_MOD = {
        # 用户设定的蛋白N端固定修饰
        'fix_PN': {},
        # 用户设定的蛋白C端固定修饰
        'fix_PC': {},
        # 用户设定的肽段N端固定修饰
        'fix_N': {},
        # 用户设定的肽段N端固定修饰
        'fix_C': {},
        # 用户设定的固定修饰
        'fix': {},
        # 用户设定的蛋白N端可变修饰
        'var_PN': {},
        # 用户设定的蛋白C端可变修饰
        'var_PC': {},
        # 用户设定的肽段N端可变修饰
        'var_N': {},
        # 用户设定的肽段N端可变修饰
        'var_C': {},
        # 用户设定的可变修饰
        'var': {}
    }


"""""""""
CSearcher的setter方法
"""""""""


def operatorInitialSearcherStatic(searcher=CSearcherStatic(), local_index=CLocalIndex()):
    """
    初始化Searcher对象中的静态项目(一次搜索只需要初始化一次)
    :param searcher: searcher对象
    :param local_index: local index对象
    :param data_package: data package
    :return: None
    """
    searcher.CFLOW3_ORDERED_MASS, searcher.CFLOW3_ORDERED_PEP_ORDER = toolGetSortIndex(
        local_index.CFLOW3_PEPTIDE.PEPTIDE_MASS)
    searcher.CFLOW3_L_AVG = local_index.CFLOW3_THEORETICAL_SPECTRA.size / len(
        local_index.CFLOW3_THEORETICAL_SPECTRA_INDEX)
    searcher.CFLOW3_ORDERED_LINKABLILITY = [local_index.CFLOW3_PEPTIDE.PEPTIDE_LINKABILITY[item[0]] for item in
                                            searcher.CFLOW3_ORDERED_PEP_ORDER]

def operatorInitialSearcherDynamic(length_of_ordered_pep_order,searcher=CSearcherStatic()):
    """
    初始化Searcher对象(对每一张谱图都需要初始化的内容)
    :param searcher: searcher对象
    :return: None
    """
    searcher.CFLOW3_ALPHA = 0
    searcher.CFLOW3_BETA_FRONT = length_of_ordered_pep_order
    searcher.CFLOW3_BETA_END = searcher.CFLOW3_BETA_FRONT - 1
    searcher.CFLOW3_DEQUE_0b01 = deque()
    searcher.CFLOW3_DEQUE_0b10 = deque()
    searcher.CFLOW3_DEQUE_0b11 = deque()
    searcher.CFLOW3_SCORE_CACHE = {}


def operatorInitialSearcherAlphaLast(linker_mass,precursor_mass, shared_peptide_mass,shared_search_static_ordered_pep_order,
                                     searcher=CSearcherDynamic()):
    """
    按照一级谱参数，找到符合质量约束的最大的肽段索引值，即{i|M[i]<=(PM-XM+tolerance)/2}}
    M为肽段质量列表，PM为母离子质量，XM为交联剂质量
    使用的算法为线性查找算法
    :param precursor_mass: 母离子质量
    :param shared_peptide_mass:
    :param data_package: data package
    :param searcher: Searcher
    :return: None
    """
    cut_value = (precursor_mass - linker_mass + searcher.CFLOW3_SPECTRUM_TOLERANCE) / 2
    for order, item in enumerate(shared_search_static_ordered_pep_order):
        if shared_peptide_mass[item[0]][item[1]] > cut_value:
            searcher.CFLOW3_ALPHA_LAST = order - 1
            break


"""""""""
CPepInfo的setter方法
"""""""""


def operatorDeleteItemInModificationAndMass(unsatisfied_item, peptide_info_class=CPeptideInfo()):
    """
    按照unsatified_item列表删除对应的质量列表和修饰列表
    :param unsatisfied_item:
    :param peptide_info_class:
    :return:
    """
    for item in unsatisfied_item:
        peptide_info_class.PEPTIDE_MODIFICATIONS[item[0]][item[1]] = None
        peptide_info_class.PEPTIDE_MASS[item[0]][item[1]] = None
    peptide_info_class.PEPTIDE_MODIFICATIONS = [list(filter(lambda x: x is not None, itema)) for itema in
                                                peptide_info_class.PEPTIDE_MODIFICATIONS]
    peptide_info_class.PEPTIDE_MASS = [list(filter(lambda x: x is not None, itema)) for itema in
                                       peptide_info_class.PEPTIDE_MASS]


def operatorDeleteItemInCPeptideInfo(remove_list, peptide_info_class=CPeptideInfo()):
    """
    按照remove_list移除CpepInfo中的项目
    :param remove_list: 需要移除的项目列表
    :param peptide_info_class: 待处理的CPepInfo对象
    :return: None
    """
    for item in remove_list:
        peptide_info_class.PEPTIDE_RANGE[item] = None
        peptide_info_class.PEPTIDE_MASS[item] = None
        peptide_info_class.PEPTIDE_PROTEIN[item] = None
        peptide_info_class.PEPTIDE_MODIFICATIONS[item] = None
        peptide_info_class.PEPTIDE_RAW_MASS[item] = None
        peptide_info_class.PEPTIDE_LINKABILITY[item]=None
    peptide_info_class.PEPTIDE_RANGE = list(filter(lambda x: x is not None, peptide_info_class.PEPTIDE_RANGE))
    peptide_info_class.PEPTIDE_MASS = list(filter(lambda x: x is not None, peptide_info_class.PEPTIDE_MASS))
    peptide_info_class.PEPTIDE_PROTEIN = list(filter(lambda x: x is not None, peptide_info_class.PEPTIDE_PROTEIN))
    peptide_info_class.PEPTIDE_RAW_MASS = list(filter(lambda x: x is not None, peptide_info_class.PEPTIDE_RAW_MASS))
    peptide_info_class.PEPTIDE_MODIFICATIONS = list(
        filter(lambda x: x is not None, peptide_info_class.PEPTIDE_MODIFICATIONS))
    peptide_info_class.PEPTIDE_LINKABILITY=list(filter(lambda x: x is not None, peptide_info_class.PEPTIDE_LINKABILITY))


"""""""""
CFileMS2的setter方法
"""""""""


def operatorInitCfileMS2(input_ms2=CFileMS2()):
    """
    初始化FileMS2类以便存储
    :param input_ms2:
    :return:
    """
    input_ms2.INDEX_SCAN = []
    input_ms2.INDEX_RT = []

    input_ms2.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    input_ms2.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    input_ms2.LIST_ACTIVATION_CENTER = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    input_ms2.LIST_PRECURSOR_SCAN = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    input_ms2.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN  # 每一行是个list
    input_ms2.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN
    input_ms2.MATRIX_CHARGE = [[] * 1] * VALUE_MAX_SCAN  # 相同的scan，可能有多个母离子状态（质量+电荷）
    input_ms2.MATRIX_MZ = [[] * 1] * VALUE_MAX_SCAN  # 相同的scan，可能有多个母离子状态（质量+电荷）


"""""""""
CDataPSM的setter方法
"""""""""

def operatorCDataPSMInitial(psm_data=CDataPSMCX()):
    psm_data.PSM_FILE_NAME = []
    psm_data.PSM_SPECTRA_TITLE = []
    psm_data.PSM_SPECTRA_INDEX = []
    psm_data.PSM_PEP1 = []
    psm_data.PSM_PEP2 = []
    psm_data.PSM_PEP1_N_TERM = []
    psm_data.PSM_PEP2_N_TERM = []
    psm_data.PSM_PEP1_MISS_CLEAVAGE = []
    psm_data.PSM_PEP2_MISS_CLEAVAGE = []
    psm_data.PSM_PEP1_MODIFICATION = []
    psm_data.PSM_PEP2_MODIFICATION = []
    psm_data.PSM_PEP1_PROTEIN = []
    psm_data.PSM_PEP2_PROTEIN = []
    psm_data.PSM_MASS = []
    psm_data.PSM_CHARGE = []
    psm_data.PSM_PEP1_COARSE_SCORE = []
    psm_data.PSM_PEP2_COARSE_SCORE = []
    psm_data.PSM_SCORE = []
    psm_data.PSM_DELTA_FINE_SCORE=[]
    psm_data.PSM_DOT_SCORE=[]
    psm_data.PSM_DECOY = []
    psm_data.PSM_PRECURSOR_MOZ = []
    psm_data.PSM_OBSERVE_MOZ = []
    psm_data.PSM_q_VALUE = []
    psm_data.PSM_SITE_q_VALUE=[]
    psm_data.PSM_SITE_q_VALUE_INFO=[]
    psm_data.PSM_CX_SITES = []
    psm_data.PSM_PEP1_IF_LINKER_ALPHA = []
    psm_data.PSM_CX_SITES_INFO = []
    psm_data.PSM_PEP1_SITE_PROBABILITY = []
    psm_data.PSM_PEP2_SITE_PROBABILITY = []
    psm_data.PSM_CX_SITE_RELIABILITY = []

def operatorCDataPSMConvertToNumpy(psm_data=CDataPSMCX()):
    """
    将CDataPSM中的各成员转化为numpy.array的形式
    除了PSM_SVM_SCORE和PSM_RANK_INDEX
    :param psm_data: CDataPSM对象
    :return: None
    """
    psm_data.PSM_FILE_NAME = np.array(psm_data.PSM_FILE_NAME)
    psm_data.PSM_SPECTRA_TITLE=np.array(psm_data.PSM_SPECTRA_TITLE)
    psm_data.PSM_SPECTRA_INDEX = np.array(psm_data.PSM_SPECTRA_INDEX)
    psm_data.PSM_PEP1 = np.array(psm_data.PSM_PEP1)
    psm_data.PSM_PEP2 = np.array(psm_data.PSM_PEP2)
    psm_data.PSM_PEP1_N_TERM=np.array(psm_data.PSM_PEP1_N_TERM)
    psm_data.PSM_PEP2_N_TERM=np.array(psm_data.PSM_PEP2_N_TERM)
    psm_data.PSM_PEP1_MISS_CLEAVAGE=np.array(psm_data.PSM_PEP1_MISS_CLEAVAGE)
    psm_data.PSM_PEP2_MISS_CLEAVAGE = np.array(psm_data.PSM_PEP2_MISS_CLEAVAGE)
    psm_data.PSM_PEP1_MODIFICATION = np.array(psm_data.PSM_PEP1_MODIFICATION)
    psm_data.PSM_PEP2_MODIFICATION = np.array(psm_data.PSM_PEP2_MODIFICATION)
    psm_data.PSM_PEP1_PROTEIN = np.array(psm_data.PSM_PEP1_PROTEIN,dtype=object)
    psm_data.PSM_PEP2_PROTEIN = np.array(psm_data.PSM_PEP2_PROTEIN,dtype=object)
    psm_data.PSM_MASS = np.array(psm_data.PSM_MASS)
    psm_data.PSM_CHARGE = np.array(psm_data.PSM_CHARGE)
    psm_data.PSM_PEP1_COARSE_SCORE=np.array(psm_data.PSM_PEP1_COARSE_SCORE)
    psm_data.PSM_PEP2_COARSE_SCORE=np.array(psm_data.PSM_PEP2_COARSE_SCORE)
    psm_data.PSM_SCORE = np.array(psm_data.PSM_SCORE)
    psm_data.PSM_DELTA_FINE_SCORE=np.array(psm_data.PSM_DELTA_FINE_SCORE)
    psm_data.PSM_DOT_SCORE=np.array(psm_data.PSM_DOT_SCORE)
    psm_data.PSM_DECOY = np.array(psm_data.PSM_DECOY)
    psm_data.PSM_PRECURSOR_MOZ=np.array(psm_data.PSM_PRECURSOR_MOZ)
    psm_data.PSM_OBSERVE_MOZ=np.array(psm_data.PSM_OBSERVE_MOZ)
    # q值和交联位点在某些步骤之前，可能为空
    psm_data.PSM_q_VALUE = np.array(psm_data.PSM_q_VALUE)
    psm_data.PSM_SITE_q_VALUE=np.array(psm_data.PSM_SITE_q_VALUE)
    psm_data.PSM_SITE_q_VALUE_INFO=np.array(psm_data.PSM_SITE_q_VALUE_INFO,dtype=object)
    psm_data.PSM_CX_SITES = np.array(psm_data.PSM_CX_SITES,dtype=object)
    psm_data.PSM_PEP1_IF_LINKER_ALPHA=np.array(psm_data.PSM_PEP1_IF_LINKER_ALPHA)
    psm_data.PSM_CX_SITES_INFO=np.array(psm_data.PSM_CX_SITES_INFO,dtype=object)
    psm_data.PSM_PEP1_SITE_PROBABILITY=np.array(psm_data.PSM_PEP1_SITE_PROBABILITY)
    psm_data.PSM_PEP2_SITE_PROBABILITY=np.array(psm_data.PSM_PEP2_SITE_PROBABILITY)
    psm_data.PSM_CX_SITE_RELIABILITY=np.array(psm_data.PSM_CX_SITE_RELIABILITY,dtype=object)


def operatorCDataPSMExtractFromIndex(index_list, psm_data=CDataPSMCX()):
    """
    提取内容，除了PSM_SVM_SCORE和PSM_q_VALUE,PSM_SITE_q_VALUE, PSM_SITE_q_VALUE_INFO
    :param index_list:
    :param psm_data:
    :return:
    """
    if len(psm_data.PSM_RANK_INDEX):
        psm_data.PSM_RANK_INDEX=list(range(len(index_list)))
    if len(psm_data.PSM_FILE_NAME):
        psm_data.PSM_FILE_NAME = psm_data.PSM_FILE_NAME[index_list]
    if len(psm_data.PSM_SPECTRA_INDEX):
        psm_data.PSM_SPECTRA_INDEX = psm_data.PSM_SPECTRA_INDEX[index_list]
    if len(psm_data.PSM_SPECTRA_TITLE):
        psm_data.PSM_SPECTRA_TITLE = psm_data.PSM_SPECTRA_TITLE[index_list]
    if len(psm_data.PSM_PEP1):
        psm_data.PSM_PEP1 = psm_data.PSM_PEP1[index_list]
    if len(psm_data.PSM_PEP2):
        psm_data.PSM_PEP2 = psm_data.PSM_PEP2[index_list]
    if len(psm_data.PSM_PEP1_N_TERM):
        psm_data.PSM_PEP1_N_TERM=psm_data.PSM_PEP1_N_TERM[index_list]
    if len(psm_data.PSM_PEP2_N_TERM):
        psm_data.PSM_PEP2_N_TERM=psm_data.PSM_PEP2_N_TERM[index_list]
    if len(psm_data.PSM_PEP1_MISS_CLEAVAGE):
        psm_data.PSM_PEP1_MISS_CLEAVAGE = psm_data.PSM_PEP1_MISS_CLEAVAGE[index_list]
    if len(psm_data.PSM_PEP2_MISS_CLEAVAGE):
        psm_data.PSM_PEP2_MISS_CLEAVAGE = psm_data.PSM_PEP2_MISS_CLEAVAGE[index_list]
    if len(psm_data.PSM_PEP1_MODIFICATION):
        psm_data.PSM_PEP1_MODIFICATION = psm_data.PSM_PEP1_MODIFICATION[index_list]
    if len(psm_data.PSM_PEP2_MODIFICATION):
        psm_data.PSM_PEP2_MODIFICATION = psm_data.PSM_PEP2_MODIFICATION[index_list]
    if len(psm_data.PSM_PEP1_PROTEIN):
        psm_data.PSM_PEP1_PROTEIN = psm_data.PSM_PEP1_PROTEIN[index_list]
    if len(psm_data.PSM_PEP2_PROTEIN):
        psm_data.PSM_PEP2_PROTEIN = psm_data.PSM_PEP2_PROTEIN[index_list]
    if len(psm_data.PSM_MASS):
        psm_data.PSM_MASS = psm_data.PSM_MASS[index_list]
    if len(psm_data.PSM_CHARGE):
        psm_data.PSM_CHARGE = psm_data.PSM_CHARGE[index_list]
    if len(psm_data.PSM_PEP1_COARSE_SCORE):
        psm_data.PSM_PEP1_COARSE_SCORE=psm_data.PSM_PEP1_COARSE_SCORE[index_list]
    if len(psm_data.PSM_PEP2_COARSE_SCORE):
        psm_data.PSM_PEP2_COARSE_SCORE = psm_data.PSM_PEP2_COARSE_SCORE[index_list]
    if len(psm_data.PSM_SCORE):
        psm_data.PSM_SCORE = psm_data.PSM_SCORE[index_list]
    if len(psm_data.PSM_DELTA_FINE_SCORE):
        psm_data.PSM_DELTA_FINE_SCORE=psm_data.PSM_DELTA_FINE_SCORE[index_list]
    if len(psm_data.PSM_DOT_SCORE):
        psm_data.PSM_DOT_SCORE=psm_data.PSM_DOT_SCORE[index_list]
    if len(psm_data.PSM_DECOY):
        psm_data.PSM_DECOY = psm_data.PSM_DECOY[index_list]
    if len(psm_data.PSM_PRECURSOR_MOZ):
        psm_data.PSM_PRECURSOR_MOZ = psm_data.PSM_PRECURSOR_MOZ[index_list]
    if len(psm_data.PSM_OBSERVE_MOZ):
        psm_data.PSM_OBSERVE_MOZ=psm_data.PSM_OBSERVE_MOZ[index_list]
    if len(psm_data.PSM_CX_SITES):
        psm_data.PSM_CX_SITES = psm_data.PSM_CX_SITES[index_list]
    if len(psm_data.PSM_CX_SITES_INFO):
        psm_data.PSM_CX_SITES_INFO = psm_data.PSM_CX_SITES_INFO[index_list]
    if len(psm_data.PSM_PEP1_IF_LINKER_ALPHA):
        psm_data.PSM_PEP1_IF_LINKER_ALPHA=psm_data.PSM_PEP1_IF_LINKER_ALPHA[index_list]
    if len(psm_data.PSM_PEP1_SITE_PROBABILITY):
        psm_data.PSM_PEP1_SITE_PROBABILITY=psm_data.PSM_PEP1_SITE_PROBABILITY[index_list]
    if len(psm_data.PSM_PEP2_SITE_PROBABILITY):
        psm_data.PSM_PEP2_SITE_PROBABILITY=psm_data.PSM_PEP2_SITE_PROBABILITY[index_list]
    if len(psm_data.PSM_CX_SITE_RELIABILITY):
        psm_data.PSM_CX_SITE_RELIABILITY=psm_data.PSM_CX_SITE_RELIABILITY[index_list]



