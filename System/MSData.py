class CINI:
    # 与质量相关的常量
    MASS_ELECTRON = 0.0005485799
    MASS_PROTON_MONO = 1.00727645224
    MASS_PROTON_AVERAGE = 1.0025
    MASS_H2O_MONO = 18.010565
    # 与元素质量相关的常量
    DICT0_ELEMENT_MASS = {}
    DICT0_ELEMENT_ABDC = {}
    # 氨基酸的化学组成常量
    DICT1_AA_COM = {}
    # 氨基酸的质量常量
    DICT1_AA_MASS = {'I': 113.08406, 'E': 129.04259, 'A': 71.03711, 'K': 128.09496, 'G': 57.02146, 'R': 156.10111,
                     'D': 115.02694, 'T': 101.04768, 'Q': 128.05858, 'L': 113.08406, 'V': 99.06841, 'M': 131.004049,
                     'N': 114.04293, 'F': 147.06841}
    # 与修饰相关的常量
    DICT2_MOD_COM = {}
    # 与糖相关的常量
    DICT3_GLYCO_COM = {}
    # 与交联剂相关的参数
    DICT4_LINKER_COM = {}
    # 与酶相关的参数
    ENZYME_CUT_LIST = []
    ENZYME_CUT_IGNORE_LIST = []
    ENZYME_TERMINAL = 'N'


class CConfig:
    # 与质谱数据相关的参数
    A1_PATH_MS2 = [r'']
    A2_TYPE_MS2 = 0
    A3_PATH_FASTA = r''
    A4_PATH_MS2_PKL = []
    A5_PATH_PROCESSED_MS2_PKL = []

    # 与酶，修饰相关的参数
    B1_NAME_ENZYME = ''
    B2_TYPE_DIGEST = 0
    B3_NUMBER_MAX_MISS_CLV = 3
    B4_NAME_MOD_FIX = [""]
    B5_NAME_MOD_VAR = [""]
    B6_NUMBER_MAX_MOD = 5
    B7_USER_MOD = {}

    # 与检索类型相关的参数
    C_TYPE_SEARCH = 3
    C_REBUILD_INDEX = 0

    # 与Flow1相关的参数
    C11_UUA_SEQ = ""
    C12_UUA_AA = ""
    C13_UUA_LEN_LOW = 6
    C14_UUA_LEN_UP = 100
    C15_UUA_COM = "C(12)H(13)N(1)O(2)"
    C16_UUA_NAME_ENZYME = "trypsin KR C"
    C17_UUA_TYPE_DIGEST = 0
    C18_UUA_NUMBER_MAX_MISS_CLV = 0
    C19_UUA_LINKED_AA = "C"

    # 与Flow2相关的参数
    C21_NG_GLYCO = ""

    # 与Flow3相关的参数
    C31_LINKER_NAME = ''
    # 与Flow3相关的参数，打分函数计算时的误差控制
    C32_MS2_TOLERANCE_TYPE = -1
    C32_MS2_ABSOLUTE_TOLERANCE = 0
    C32_MS2_RELATIVE_TOLERANCE = 0
    C33_MS1_TOLERANCE_TYPE = -1
    C33_MS1_ABSOLUTE_TOLERANCE = 0
    C33_MS1_RELATIVE_TOLERANCE = 0
    # 与Flow3相关的参数，粗打分过程控制的参数（检索谱图的质荷比范围和保留粗打分的结果数目）
    C34_MS2_MOZ_UPPER_LIMIT = 60000
    C34_MS2_MOZ_LOWER_LIMIT = 0
    C34_RSV_K1 = 5.0
    C34_RSV_B = 1.0
    C35_ALPHA_PURNING_SCORE = 0
    C35_BETA_DEQUE_MIN_LEN = 0
    C35_TOP_N_FOR_ONE_SPECTRUM = 15
    C35_MIN_COARSE_FOR_ONE_SPECTRUM=0
    # 与Flow3相关的参数，PSM对象的保存路径
    C36_RAW_PSM_PKL = []
    C36_FINE_SCORING_PSM_PKL = []
    C36_CX_SVM_PSM_PKL = []
    C36_SITE_SCORING_PSM_PKL = []
    C36_SITE_FDR_PKL = []
    # 与Flow3相关的参数，位点打分过程的控制参数
    C37_MS2_MOZ_UPPER_LIMIT = 1500
    C37_MS2_MOZ_LOWER_LIMIT = 100
    C37_MS2_SPLIT_RANGE = 100
    C37_MS2_MAX_SPECTRA_DEPTH = 10
    # 与Flow3相关的隐藏参数
    C3H_WRITE_INDEX = 0
    C3H_PREPROCESSING_PEAK_HIGHER_CHARGE = 0
    C3H_PREPROCESSING_ASSUME_CX_ION_INTENSITY_RATE = 0
    C3H_WRITE_COARSE_SCORE_RESULT = 0

    # 与并行计算蛋白质长度相关的参数
    P1_NUMBER_THREAD = 10
    P2_TYPE_THREAD = 0
    P3_NUMBER_SELECT_PACK = 200
    P4_LEN_MAX_PROTEIN = 100000
    P5_MASS_PEP_LOW = 600
    P6_MASS_PEP_UP = 10000
    P7_LEN_PEP_LOW = 6
    P8_LEN_PEP_UP = 100
    P9_INDEX_SPLIT_MASS = 100
    P10_NUMBER_TOP_RESULT = 3
    P11_MULTI_MASS = 1000
    P12_TYPE_START = 0

    # 与质量控制相关的参数
    F1_THRESHOLD_FOR_PSM = 0.5
    F2_ENABLE_TARGET_DECOY = 1
    F3_THRESHOLD_FOR_SITE = 1

    # 与质量基本参数相关的参数
    I0_INI_PATH_ELEMENT = ''
    I1_INI_PATH_AA = ''
    I2_INI_PATH_MOD = ''
    I3_INI_PATH_GLYCO = r'glyco.ini'
    I4_INI_PATH_LINKER = r'.\ini\linker.ini'
    I5_INI_PATH_ENZYME = r'C:\pFindStudio\pFind3\bin\enzyme.ini'

    # 与输出相关的参数
    E1_PATH_EXPORT = r''
    E2_TYPE_EXPORT = 0
    E3_FLAG_CREATE_NEM_FOLDER = 0
    E4_FLAG_EXPORT_EVIDENCE = 0
    E5_LOCAL_LOGGER = r''
    E5_LOG_FILE_PATH=r''
    E6_LOCAL_LOG_FILE_PATH = r''


class CFileMS2:
    # 二级谱文件的类
    INDEX_SCAN = [0]
    INDEX_RT = []

    LIST_RET_TIME = []
    LIST_ION_INJECTION_TIME = []
    LIST_ACTIVATION_CENTER = []
    LIST_PRECURSOR_SCAN = []

    # 每一行是一个list
    MATRIX_PEAK_MOZ = []
    MATRIX_PEAK_INT = []
    MATRIX_CHARGE = []
    MATRIX_MZ = []


class CSearcherStatic:
    # 搜索核心模块所需要的数据结构
    CFLOW3_ORDERED_PEP_ORDER = []
    CFLOW3_ORDERED_MASS = []
    CFLOW3_ORDERED_LINKABLILITY = []
    CFLOW3_L_AVG = 0


class CSearcherDynamic:
    CFLOW3_ALPHA = 0
    CFLOW3_ALPHA_LAST = 0
    CFLOW3_BETA_FRONT = 0
    CFLOW3_BETA_END = 0

    CFLOW3_SPECTRUM_TOLERANCE = 0
    CFLOW3_DEQUE_0b01 = 0
    CFLOW3_DEQUE_0b10 = 0
    CFLOW3_DEQUE_0b11 = 0
    CFLOW3_SCORE_CACHE = {}


class CPeptideInfo:
    # 存储肽段所对应的蛋白质编号，列表中每一项是每个肽段的一个列表，列表中存储了该肽段对应的蛋白质编号
    PEPTIDE_PROTEIN = []
    # 存储肽段的序列编码，列表中每一项是每个肽段一个元组，元组中存储了其对应PEPTIDE_PROTEIN中第一个蛋白的（起，止）
    PEPTIDE_RANGE = []
    # 肽段的质量，列表中每一项为对应肽段的质量信息
    PEPTIDE_RAW_MASS = []
    # 含有修饰的肽段质量，列表中每一项为对应肽段的修饰字典列表的质量列表
    PEPTIDE_MASS = []
    # 存储肽段的修饰字典，列表中每一项是一个修饰字典列表
    PEPTIDE_MODIFICATIONS = []
    # 存储肽段的可交联性信息
    PEPTIDE_LINKABILITY = []


class CLocalIndex:
    # 蛋白编号索引
    CFLOW3_PROTEIN_INDEX = {}
    # 蛋白质描述索引
    CFLOW3_DESCRIPTION_INDEX = {}
    # 蛋白质序列索引
    CFLOW3_SEQUENCE_INDEX = {}

    # 肽段信息
    CFLOW3_PEPTIDE = CPeptideInfo()

    # 理论二级谱的质荷比信息
    CFLOW3_THEORETICAL_SPECTRA = []
    # 理论二级谱图中b,y粒子的原始索引
    CFLOW3_THEORETICAL_SPECTRA_SORT_BY_INDEX = []
    # 理论二级谱的长度信息
    CFLOW3_THEORETICAL_SPECTRA_INDEX = []


class CDataPSMCX:
    # 排序的索引
    PSM_RANK_INDEX = []
    # 谱图所属文件名
    PSM_FILE_NAME = []
    # 谱图扫描编号列表
    PSM_SPECTRA_INDEX = []
    # 谱图的真实标题
    PSM_SPECTRA_TITLE = []
    # 肽段1序列
    PSM_PEP1 = []
    # 肽段1是否是N端
    PSM_PEP1_N_TERM = []
    # 肽段2序列
    PSM_PEP2 = []
    # 肽段2是否是N端
    PSM_PEP2_N_TERM = []
    # 肽段1序列的漏切数目
    PSM_PEP1_MISS_CLEAVAGE = []
    # 肽段2序列的漏切数目
    PSM_PEP2_MISS_CLEAVAGE = []
    # 肽段的交联位点
    PSM_CX_SITES = []
    # 肽段1是否和交联剂的alpha端进行了反应
    PSM_PEP1_IF_LINKER_ALPHA = []
    # 肽段的交联位点的详细打分信息
    PSM_CX_SITES_INFO = []
    # 肽段1的蛋白归属
    PSM_PEP1_PROTEIN = []
    # 肽段2的蛋白归属
    PSM_PEP2_PROTEIN = []
    # 肽段1的修饰字典
    PSM_PEP1_MODIFICATION = []
    # 肽段2的修饰字典
    PSM_PEP2_MODIFICATION = []
    # 肽段1的粗打分值
    PSM_PEP1_COARSE_SCORE = []
    # 肽段2的粗打分值
    PSM_PEP2_COARSE_SCORE = []
    # 肽段1的位点概率分布
    PSM_PEP1_SITE_PROBABILITY = []
    # 肽段2的位点概率分布
    PSM_PEP2_SITE_PROBABILITY = []
    # 交联肽的位点可信度
    PSM_CX_SITE_RELIABILITY = []
    # 细打分打分值
    PSM_SCORE = []
    # 细打分差值
    PSM_DELTA_FINE_SCORE = []
    # 点积分数
    PSM_DOT_SCORE = []
    # SVM打分值
    PSM_SVM_SCORE = []
    # 母离子质量
    PSM_MASS = []
    # 母离子电荷
    PSM_CHARGE = []
    # 母离子质荷比
    PSM_PRECURSOR_MOZ = []
    # 是否存在反库中的蛋白
    PSM_DECOY = []
    # q-值计算
    PSM_q_VALUE = []
    # 位点q值计算
    PSM_SITE_q_VALUE = []
    # 位点q值对应的位点信息
    PSM_SITE_q_VALUE_INFO = []
    # 实际检索到的结果的母离子质荷比
    PSM_OBSERVE_MOZ = []


class CRescore:
    # 支持向量机的特征矩阵————交联两条肽段都有碎裂的情况
    SVM_RANK_FEATURE_CX_TWO_CLEAVAGE = []
    # 支持向量机的特征矩阵————交联只有一条肽段都有碎裂的情况
    SVM_RANK_FEATURE_CX_ONE_CLEAVAGE = []
    # 支持向量机中使用的CDataPSM样本序号————交联两条肽段都有碎裂的情况
    SVM_PSM_INDEX_CX_TWO_CLEAVAGE = []
    # 支持向量机中使用的CDataPSM样本序号————交联一条肽段都有碎裂的情况
    SVM_PSM_INDEX_CX_ONE_CLEAVAGE = []
    # 支持向量机中选取样本时，取top-N的进行分析(默认为top1)
    SVM_SAMPLE_TOP_N = 1
    # 支持向量机中选取正例样本时，取FDR小于此值的作为正例（默认为1%）
    SVM_SAMPLE_POSITIVE_FDR = 0.010
    # 支持向量机的迭代次数
    SVM_ITERATION_TIME = 20
    # 支持向量机的C值
    SVM_C_VALUE = 0.5


class CDataPack:
    # config文件和ini文件的存储列表
    my_config = CConfig()
    my_ini = CINI()

    # 保存索引内容以及检索核心部分数据
    my_index = CLocalIndex()

    # 保存重打分核心部分数据
    my_svm = CRescore()
