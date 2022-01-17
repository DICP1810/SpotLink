import os
import pickle
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

from System.MSData import CDataPack, CFileMS2


def operatorDumpIndexes(file_name, data_package=CDataPack()):
    """
    存储蛋白质肽段索引至硬盘
    :param file_name: 文件名称
    :param data_package: data package
    :return: None
    """
    with open(file_name, 'wb') as index_pickle:
        pickle.dump((data_package.my_index.CFLOW3_DESCRIPTION_INDEX,
                     data_package.my_index.CFLOW3_PEPTIDE,
                     data_package.my_index.CFLOW3_SEQUENCE_INDEX,
                     data_package.my_index.CFLOW3_PROTEIN_INDEX,
                     data_package.my_index.CFLOW3_THEORETICAL_SPECTRA,
                     data_package.my_index.CFLOW3_THEORETICAL_SPECTRA_SORT_BY_INDEX,
                     data_package.my_index.CFLOW3_THEORETICAL_SPECTRA_INDEX), index_pickle)


def operatorLoadIndexes(file_name, data_package=CDataPack()):
    """
    加载蛋白质肽段索引至硬盘
    :param file_name: 文件名称
    :param data_package: data package
    :return: None
    """
    with open(file_name, 'rb') as index_pickle:
        data_package.my_index.CFLOW3_DESCRIPTION_INDEX, \
        data_package.my_index.CFLOW3_PEPTIDE, \
        data_package.my_index.CFLOW3_SEQUENCE_INDEX, \
        data_package.my_index.CFLOW3_PROTEIN_INDEX, \
        data_package.my_index.CFLOW3_THEORETICAL_SPECTRA, \
        data_package.my_index.CFLOW3_THEORETICAL_SPECTRA_SORT_BY_INDEX, \
        data_package.my_index.CFLOW3_THEORETICAL_SPECTRA_INDEX = pickle.load(index_pickle)


def operatorSavePSM(path, ms2_object=CFileMS2(), data_package=CDataPack(), extension_name='rpsm'):
    """
    将一个二级谱文件对应的PSM对象保存到硬盘中，并将文件名记录在data package中
    :param path:二级谱文件名（可以含有路径）
    :param ms2_object: CFileMS2对象
    :param data_package: data package
    :param extension_name: 所需要保存的二进制文件名后缀 （rpsm:原始PSM数据，fpsm:经过细打分后的PSM数据)
    :return: 返回psm包含完整路径的文件名称
    """
    base_name = os.path.splitext(os.path.basename(path))[0]
    psm_name = os.path.normpath(os.path.join(data_package.my_config.E1_PATH_EXPORT, base_name + '.' + extension_name))
    with open(psm_name, 'wb') as psm_f:
        pickle.dump(ms2_object, psm_f, protocol=4)
    return psm_name


def operatorLoadPSM(psm_file_path):
    """
    加载PSM信息到内存中
    :param psm_file_path: psm对象存储文件的路径
    :return: None
    """
    with open(psm_file_path, 'rb') as psm_f:
        ms2_object = pickle.load(psm_f)
    return ms2_object


def operatorLoadMS2(ms2_file_path):
    """
    将保存于硬盘中的ms2文件读取到内存中
    :param ms2_file_path: 二级谱pkl文件
    :return: CFileMS2对象
    """
    with open(ms2_file_path, 'rb') as ms2_f:
        ms2_object = pickle.load(ms2_f)
    return ms2_object


def operatorGetTitleIndexFromMgf(file_path):
    """
    建立一个从（谱图扫描号，质荷比，电荷）到谱图真实Title的映射，这部分信息在将谱图转化为MS2对象的过程中丢失了
    :param file_path: mgf文件的路径
    :return: 一个索引
    """
    title_index = {}
    with open(file_path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line:
                if line.startswith('TITLE'):
                    title = line.split('=')[1]
                    scan_number = line.split('.')[1]
                elif line.startswith('CHARGE'):
                    charge = line.split('=')[1][:-1]
                elif line.startswith('PEPMASS'):
                    moz = float(line.split('=')[1])
                elif line.startswith('END'):
                    index_string = '{}-{:.5f}-{}'.format(scan_number, moz, charge)
                    title_index[index_string] = title
    return title_index


def operatorProxlMxlGenerator(psm_object, peptide_dic, data_package):
    fasta_name = os.path.basename(data_package.my_config.A3_PATH_FASTA)
    linker_name = data_package.my_config.C31_LINKER_NAME
    linker_mass = data_package.my_ini.DICT4_LINKER_COM[data_package.my_config.C31_LINKER_NAME][0]
    # 根节点
    level_1_node = Element('proxl_input')
    level_1_node.set('fasta_filename', fasta_name)
    # SpotLink相关参数
    level_2_1_node = SubElement(level_1_node, 'search_program_info')
    level_2_1_1_node = SubElement(level_2_1_node, 'search_programs')
    level_2_1_1_1_node = SubElement(level_2_1_1_node, 'search_program', name='SpotLink', display_name='SpotLink',
                                    version='1')
    level_2_1_1_1_1_node = SubElement(level_2_1_1_1_node, 'psm_annotation_types')
    level_2_1_1_1_1_1_node = SubElement(level_2_1_1_1_1_node, 'filterable_psm_annotation_types')
    level_2_1_1_1_1_1_1_node = SubElement(level_2_1_1_1_1_1_node, 'filterable_psm_annotation_type', name='FDR',
                                          description='False discovery rate', filter_direction='below',
                                          default_filter='true', default_filter_value='0.1')
    level_2_1_1_1_1_2_node = SubElement(level_2_1_1_1_1_node, 'descriptive_psm_annotation_types')
    level_2_1_1_1_1_2_1_node = SubElement(level_2_1_1_1_1_2_node, 'descriptive_psm_annotation_type', name='m/z',
                                          description='m/z')
    level_2_1_1_1_1_2_2_node = SubElement(level_2_1_1_1_1_2_node, 'descriptive_psm_annotation_type', name='scan num.',
                                          description='scan num.')
    level_2_1_2_node = SubElement(level_2_1_node, 'default_visible_annotations')
    level_2_1_2_1_node = SubElement(level_2_1_2_node, 'visible_psm_annotations')
    level_2_1_2_1_1_node = SubElement(level_2_1_2_1_node, 'search_annotation', search_program='SpotLink',
                                      annotation_name='scan num.')
    level_2_1_2_1_2_node = SubElement(level_2_1_2_1_node, 'search_annotation', search_program='SpotLink',
                                      annotation_name='score')
    level_2_1_2_1_3_node = SubElement(level_2_1_2_1_node, 'search_annotation', search_program='SpotLink',
                                      annotation_name='FDR')
    level_2_1_2_1_4_node = SubElement(level_2_1_2_1_node, 'search_annotation', search_program='SpotLink',
                                      annotation_name='m/z')
    # 交联剂参数
    level_3_1_node = SubElement(level_1_node, 'linkers')
    level_3_1_1_node = SubElement(level_3_1_node, 'linker', name=linker_name)
    level_3_1_1_1_node = SubElement(level_3_1_1_node, 'crosslink_masses')
    level_3_1_1_1_1_node = SubElement(level_3_1_1_1_node, 'crosslink_mass', mass=linker_mass)
    # PSM参数
    level_4_1_node = SubElement(level_1_node, 'reported_peptides')
    level_4_1_children_nodes = [operatorCreateSinglePeptideItem(item, peptide_dic[item], psm_object, data_package) for
                                item in peptide_dic]
    level_4_1_node.extend(level_4_1_children_nodes)
    # 蛋白质参数
    level_5_1_node = SubElement(level_1_node, 'matched_proteins')
    level_5_1_children_nodes = operatorCreateProteinChildrenNode(peptide_dic, psm_object, data_package)
    level_5_1_node.extend(level_5_1_children_nodes)
    # 修饰参数
    # level_6_1_node=SubElement(level_1_node,'static_modifications')
    # 组织与格式化
    rough_string = tostring(level_1_node, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    result_string = reparsed.toprettyxml(indent=" ")
    return result_string


def operatorCreateSinglePeptideItem(item, spectra_index_list, psm_object, data_package):
    node1 = Element('reported_peptide', reported_peptide_string=item, type='crosslink')
    node2 = SubElement(node1, 'peptides')
    pep1, pep2, site1, site2, modi1, modi2 = spectra_index_list[1:]
    # 肽段1
    node3 = SubElement(node2, 'peptide', sequence=pep1)
    if modi1:
        node6 = SubElement(node3, 'modifications')
        modi1_nodes = [Element('modification', mass='{:.2f}'.format(modi1[site]), position=str(int(site) + 1), isMonolink='false')
                       for site in
                       modi1]
        node6.extend(modi1_nodes)
    node4 = SubElement(node3, 'linked_positions')
    node5 = SubElement(node4, 'linked_position', position=site1)
    # 肽段2
    node7 = SubElement(node2, 'peptide', sequence=pep2)
    if modi2:
        node10 = SubElement(node7, 'modifications')
        modi2_nodes = [Element('modification',  mass='{:.2f}'.format(modi2[site]), position=str(int(site) + 1), isMonolink='false')
                       for site in
                       modi2]
        node10.extend(modi2_nodes)
    node8 = SubElement(node7, 'linked_positions')
    node9 = SubElement(node8, 'linked_position', position=site2)
    node11 = SubElement(node1, 'psms')
    psms_nodes = []
    for psm_index in spectra_index_list[0]:
        temp_psm = Element('psm', scan_file_name=psm_object.PSM_SPECTRA_TITLE[psm_index],
                           scan_number=str(psm_object.PSM_SPECTRA_INDEX[psm_index]),
                           precursor_charge=str(psm_object.PSM_CHARGE[psm_index]),
                           linker_mass=str(
                               data_package.my_ini.DICT4_LINKER_COM[data_package.my_config.C31_LINKER_NAME][0]))
        filter_node_1 = SubElement(temp_psm, 'filterable_psm_annotations')
        filter_node_2 = SubElement(filter_node_1, 'filterable_psm_annotation', search_program='SpotLink',
                                   annotation_name='FDR', value=str(psm_object.PSM_q_VALUE[psm_index]))
        describe_node_1 = SubElement(temp_psm, 'descriptive_psm_annotations')
        describe_node_2 = SubElement(describe_node_1, 'descriptive_psm_annotation', search_program='SpotLink',
                                     annotation_name='m/z', value=str(psm_object.PSM_OBSERVE_MOZ[psm_index]))
        psms_nodes.append(temp_psm)
    node11.extend(psms_nodes)
    return node1


def operatorCreateProteinChildrenNode(peptide_dic, psm_object, data_package):
    total_protein_sequence_dictionary = {}
    total_sequence_protein_dictionary = {}
    protein_children_node_list = []
    for order in data_package.my_index.CFLOW3_PROTEIN_INDEX:
        protein_name = data_package.my_index.CFLOW3_PROTEIN_INDEX[order]
        protein_sequence = data_package.my_index.CFLOW3_SEQUENCE_INDEX[order]
        protein_description = data_package.my_index.CFLOW3_DESCRIPTION_INDEX[order]
        total_protein_sequence_dictionary[protein_name] = protein_sequence
        if protein_sequence not in total_sequence_protein_dictionary:
            total_sequence_protein_dictionary[protein_sequence] = [(protein_name, protein_description)]
        else:
            total_sequence_protein_dictionary[protein_sequence].append((protein_name, protein_description))
    sequence_list = []
    for cx_pair in peptide_dic:
        psm_indexes = peptide_dic[cx_pair][0]
        for psm_index in psm_indexes:
            pep1_protein_list, pep2_protein_list = psm_object.PSM_PEP1_PROTEIN[psm_index], psm_object.PSM_PEP2_PROTEIN[
                psm_index]
            for pep_protein_name in pep1_protein_list:
                sequence_list.append(total_protein_sequence_dictionary[pep_protein_name])
            for pep_protein_name in pep2_protein_list:
                sequence_list.append(total_protein_sequence_dictionary[pep_protein_name])
    for sequence_item in set(sequence_list):
        temp_protein_node = Element('protein', sequence=sequence_item)
        protein_nodes = [Element('protein_annotation', name=protein_info_item[0], description=protein_info_item[1]) for
                         protein_info_item in total_sequence_protein_dictionary[sequence_item]]
        temp_protein_node.extend(protein_nodes)
        protein_children_node_list.append(temp_protein_node)
    return protein_children_node_list
