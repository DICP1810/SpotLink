from Tool.MSTool import toolInsertDictionary
from copy import deepcopy

def operatorDFSVariableSites(var_sites, max_var_num, recur_depth, recur_sites, result_site_combination):
    """
    使用深度优先遍历算法进行可能发生可变修饰位点的遍历
    :param var_sites: 一个列表，存储了肽段中可能发生修饰的位点
    :param max_var_num: 这个肽段中可能发生的可变修饰的最大数目
    :param recur_depth: 当前递归深度
    :param recur_sites: 当前递归层所对应位点列表
    :param result_site_combination: 最终返回的所有可能的可变修饰的位点组合
    :return: None
    """
    if recur_depth >= len(var_sites):
        result_site_combination.append(recur_sites)
        return
    operatorDFSVariableSites(var_sites, max_var_num, recur_depth + 1, recur_sites, result_site_combination)
    if len(recur_sites) < max_var_num:
        new_recur_sites = recur_sites.copy()
        new_recur_sites.append(var_sites[recur_depth])
        operatorDFSVariableSites(var_sites, max_var_num, recur_depth + 1, new_recur_sites, result_site_combination)

def operatorDFSVariablePeptides(site_combination,var_modi_dic,recur_modi_number,recur_modi, variable_peptide_list):
    """
    使用深度优先遍历算法进行生成所有可能的肽段修饰形式
    :param site_combination: 一个肽段的可变修饰位点列表，如[2, 3, 4]
    :param var_modi_dic: 一个肽段的可变修饰位点字典，字典的键是位点编号，字典的值是所有可能的修饰元组，如('NORMAL','C',57.02)
    :param recur_modi_number: 当前递归层所代表的site_combination序号
    :param recur_modi: 当前递归层所代表的修饰形式
    :param variable_peptide_list: 最终返回的所有可能的可变修饰的组合
    :return: None
    """
    if recur_modi_number==len(site_combination):
        variable_peptide_list.append(recur_modi)
        return
    modi_list=var_modi_dic[site_combination[recur_modi_number]]
    for single_modi in modi_list[site_combination[recur_modi_number]]:
        new_modi=deepcopy(recur_modi)
        toolInsertDictionary(new_modi,site_combination[recur_modi_number],single_modi[2])
        operatorDFSVariablePeptides(site_combination,var_modi_dic,recur_modi_number+1,new_modi,variable_peptide_list)


def operatorGenerateVariablePeptides(result_site_combination, var_modi_dic, fix_modification_dic):
    """
    生成可变修饰肽段
    :param result_site_combination:待计算肽段所有发生可变修饰的位点组合列表
    :param var_modi_dic: 可变修饰的位点所对应的修饰字典列表
    :param fix_modification_dic: 待计算肽段已发生的固定修饰字典
    :return: variable_final_list
    """
    variable_final_list=[]
    for site_combination in result_site_combination:
        if len(site_combination)==0:
            variable_final_list.append(fix_modification_dic)
            continue
        variable_modification_list=[]
        operatorDFSVariablePeptides(site_combination,var_modi_dic,0, {},variable_modification_list)
        for variable_dic in variable_modification_list:
            new_variable_dic=deepcopy(fix_modification_dic)
            new_variable_dic.update(variable_dic)
            variable_final_list.append(new_variable_dic)
    return variable_final_list
