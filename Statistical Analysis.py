import os
from statsmodels.api import OLS
import re
from statsmodels.discrete.discrete_model import Logit
from pySmartDL import SmartDL
import networkx as nx
import statistics
import random
import matplotlib.image as mpimg
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from collections import Counter
from sklearn import feature_selection, preprocessing, metrics, linear_model
from scipy.stats import kstest, ks_2samp, spearmanr, pearsonr, zscore, ttest_ind, normaltest, chi2, chi2_contingency, \
    f_oneway, linregress, gmean, mannwhitneyu, fisher_exact, sem, gmean, mstats, binomtest, linregress


def add_meta(region, trait, type, group):
    """
    Adds the metadata of the tRF data
    :param region: brain region
    :param trait: which trait should be analyzed
    :param type: sncRNA type
    :param group: Females / Males / All
    :return: Created new file with metadata
    """
    df = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\results\KStest\%s\%s\%s\kstest_%s_9_chol.csv' % (
    region, trait, type, group)).reset_index().iloc[:, 1:].set_index(type)
    path = r'C:\Users\danas\PycharmProjects\project\data\tRF\%s\tRF_meta_%s.csv' % (region, region)
    meta = pd.read_csv(path).set_index(type)
    meta['amino_acid'] = meta['Sequence locations in tRNA space (comma deliminated)'].str.split('_').str[1].str[:3]
    meta['codon'] = meta['Sequence locations in tRNA space (comma deliminated)'].str.split('_').str[1].str[3:6]
    meta['origin'] = meta['Sequence locations in tRNA space (comma deliminated)'].str.split('_').str[0].str[4:]
    meta.loc[meta.origin.str.startswith('lookalike'), 'origin'] = 'MT'
    meta.loc[meta.origin != 'MT', 'origin'] = 'Nuclear'
    meta = meta[meta.index.isin(df.index)]
    df = df[df.index.isin(meta.index)]
    comb = pd.concat([df, meta], axis=1)
    comb = comb[
        ['control_mean', 'ad_mean', 'log_fold_change', 'p_value', 'adj_p_value', 'significant', 'tRF sequence',
         'tRF type(s)', 'amino_acid', 'codon', 'origin', 'core_cholinergic_Production',
         'cholinergic_targets_Production', 'cholinergic_Production']]
    comb = comb.sort_values('adj_p_value')
    comb.to_csv(r'C:\Users\danas\PycharmProjects\project\results\KStest\%s\%s\%s\kstest_%s_9_chol.csv' % (
    region, trait, type, group))


def run_kstest(trait, type, cutoff, a, brain_region, group, sig, create_file, random_bool):
    """
    Runs the Kolmogorov-Smirnov test on RNA-Seq data from different brain regions
    :param trait: The trait the should be analyzed
    :param type: miRNA / tRF / mRNA
    :param cutoff: The median expression threshold
    :param a: Significance threshold
    :param brain_region: Brain tissue from which the data was taken
    :param group: Females / Males / All
    :param sig: Run the analysis only on features from the given list
    :param create_file: Boolean value
    :param random_bool: Whether random shuffling should be performed
    :return:
    """
    cond_dict = {'cognition': 'cog_cond', 'braak': 'braak_cond', 'condition': 'condition', 'cerad': 'cerad_cond',
                 'chol':'cog_cond'}
    chol_genes = list(pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\cholinergic_genes.csv')['Gene'])
    info = pd.read_csv(
        r'C:\Users\danas\PycharmProjects\project\info\%s\%s_%s_info.csv' % (brain_region, brain_region, group))
    if random_bool:
        info[cond_dict[trait]] = info[cond_dict[trait]].sample(frac=1).values
    control_info = info[info[cond_dict[trait]] == 'control'].reset_index()
    ad_info = info[info[cond_dict[trait]] == 'AD'].reset_index()
    control_buids = list(map(str, list(control_info['BUID'])))
    ad_buids = list(map(str, list(ad_info['BUID'])))
    buids = control_buids + ad_buids

    counts = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\%s\%s\norm_counts_%s_%s.csv' % (
    type, brain_region, brain_region, group)).set_index(type)[buids]
    if type == 'mRNA':
        counts = counts[counts.index.isin(chol_genes)]
    if sig:
        counts = counts[counts.index.isin(sig)]
    else:
        counts = counts[counts.median(axis=1) >= cutoff]
    control_counts = counts[control_buids].values.tolist()
    ad_counts = counts[ad_buids].values.tolist()

    p_values = []
    log_changes = []
    ad_means = []
    control_means = []
    for i in range(len(control_counts)):
        p_value = kstest(ad_counts[i], control_counts[i])[1]
        p_values.append(p_value)
        ad_mean = np.mean(ad_counts[i])
        control_mean = np.mean(control_counts[i])
        if ad_mean != 0 and control_mean != 0:
            log_changes.append(math.log2(ad_mean) - math.log2(control_mean))
        else:
            log_changes.append(0)
        control_means.append(control_mean)
        ad_means.append(ad_mean)

    counts['log_fold_change'] = log_changes
    counts['p_value'] = p_values
    counts['ad_mean'] = ad_means
    counts['control_mean'] = control_means
    results = counts.sort_values('p_value')
    sigs, adj_p_values, a, b = multipletests(list(results['p_value']), alpha=a, method='fdr_bh', is_sorted=True,
                                             returnsorted=True)
    results['adj_p_value'] = adj_p_values
    results['significant'] = sigs
    results = results[['control_mean', 'ad_mean', 'log_fold_change', 'p_value', 'adj_p_value', 'significant']]
    results1 = results[results['significant'] == True]
    if create_file:
        if sig:
            results.to_csv(r'C:\Users\danas\PycharmProjects\project\results\KStest\%s\%s\%s\kstest_%s_sig.csv' % (
            brain_region, trait, type, group))
        else:
            results.to_csv(
                r'C:\Users\danas\PycharmProjects\project\results\KStest\%s\%s\%s\kstest_%s_%d.csv' % (
                brain_region, trait, type, group, cutoff))

    return len(results1)


def create_chol_table(brain_region, cond, th, type, role):
    """
    Defining CholinotRFs and CholinomiRs based on the number and role of their cholinergic predicted targets
    :param brain_region: Analyzed brain region
    :param cond: By which trait the condition was set
    :param th: The median epression threshold
    :param type: The type of the data
    :param role: Production vs. Intake
    :return: Adds to the results files the cholinergic targets of each sncRNA and boolean value that represents
    their potential role in cholinergic rgulation
    """
    r = {'nuc': 'Nucleus accumbens', 'hyp': 'Hypothalamus', 'ITG': 'Cortex'}
    chol = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\cholinergic_genes.csv')
    chol = chol[chol['Role'] == role].reset_index().iloc[:, 1:]
    core_genes = list(chol[chol['Core'] == True]['Gene'])
    chol_genes = list(chol['Gene'])
    groups = ['females', 'males', 'all']
    dfs = []
    chol_dfs = []

    for group in groups:
        dfs.append(pd.read_csv(r'C:\Users\danas\PycharmProjects\project\results\KStest\%s\%s\%s\kstest_%s_%d.csv' %
                               (brain_region, cond, type, group, th)))

    if type == 'tRF':
        path = r'C:\Users\danas\PycharmProjects\tRF_targets\tRF_targets\%s' % r[brain_region]
        idx = -8
        suffix = '0.8.csv'
    else:
        path = r'C:\Users\danas\PycharmProjects\tRF_targets\miRNA_targets\%s' % r[brain_region]
        idx = -4
        suffix = 'p.csv'

    for i, df in enumerate(dfs):
        df_trfs = list(df[type])
        chol_df = pd.DataFrame(df_trfs, columns=[type])
        chol_df['Cholinergic_Score_%s' % role] = 0
        chol_df['Cholinergic_Targets_%s' % role] = ''
        chol_df['Core_Targets_%s' % role] = 0
        chol_dfs.append(chol_df)

    for file in os.listdir(path):
        if file.endswith(suffix):
            df1 = pd.read_csv(path + '\\%s' % file)
            df_chol = df1[df1['Gene'].isin(chol_genes)].reset_index().iloc[:, 1:]
            df_core = df1[df1['Gene'].isin(core_genes)].reset_index().iloc[:, 1:]
            for i, df in enumerate(chol_dfs):
                if file[:idx] in list(df.index):
                    df.loc[df[type] == file[:idx], 'Cholinergic_Targets_%s' % role] = ','.join(list(df_chol['Gene']))
                    df.loc[df[type] == file[:idx], 'Cholinergic_Score_%s' % role] = len(list(df_chol['Gene']))
                    df.loc[df[type] == file[:idx], 'Core_Targets_%s' % role] = len(list(df_core['Gene']))

    t = 1
    if role == 'Production':
        t = chol_a_df['Cholinergic_Score_%s' % role].quantile(0.8)
        print(t)

    for i, df in enumerate(chol_dfs):
        df['chol'] = df['Cholinergic_Score_%s' % role] >= t
        df['core'] = df['Core_Targets_%s' % role] >= 1
        dfs[i]['core_cholinergic_%s' % role] = df['core'].copy()
        dfs[i]['cholinergic_targets_%s' % role] = df['Cholinergic_Targets_%s' % role].copy()
        dfs[i]['cholinergic_%s' % role] = df['chol'] | df['core']
        dfs[i].to_csv(r'C:\Users\danas\PycharmProjects\project\results\KStest\%s\%s\%s\kstest_%s_%d_chol.csv' %
                      (brain_region, cond, type, groups[i], th), index=False)


def run_statistics_on_cells(time, sex):
    """
    Analysis of the cell files from two time points and females and males saparetly
    :param time: 2 or 4 days
    :param sex: 2 for females, 5 for males
    :return: Creates analysis files for the cell lines
    """
    d = {2: 'Female', 5: 'Male'}
    info = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\Cells\info_%d.csv' % time)
    counts = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\Cells\tRF_norm_counts_%d.csv' % time)\
        .set_index('tRF')
    counts = counts[counts.median(axis=1) >= 9]
    p_vals = []
    log_changes = []
    control_means = []
    cond_means = []

    controls = list(info[info['name'].str.startswith('LA%dCO%d' % (sex, time))]['name'])
    cond = list(info[info['name'].str.startswith('LA%dCN%d' % (sex, time))]['name'])
    for trf in list(counts.index):
        trf_counts = counts[counts.index == trf]
        control_counts = trf_counts[controls].values[0].astype(np.float64)
        chol_counts = trf_counts[cond].values[0].astype(np.float64)
        p_value = ttest_ind(chol_counts, control_counts)[1]
        p_vals.append(p_value)
        cond_mean = np.mean(chol_counts)
        control_mean = np.mean(control_counts)
        if cond_mean != 0 and control_mean != 0:
            log_changes.append(math.log2(cond_mean) - math.log2(control_mean))
        else:
            log_changes.append(0)
        control_means.append(control_mean)
        cond_means.append(cond_mean)
    df = pd.DataFrame(list(counts.index), columns=['tRF'])
    df['log_fold_change'] = log_changes
    df['p_value'] = p_vals
    df['cond_mean'] = cond_means
    df['control_mean'] = control_means
    results = df.sort_values('p_value')
    sigs, adj_p_values, a, b = multipletests(list(results['p_value']), alpha=0.1, method='fdr_bh',
                                             is_sorted=True,
                                             returnsorted=True)
    results['adj_p_value'] = adj_p_values
    results['significant'] = sigs
    results.to_csv(
        r'C:\Users\danas\PycharmProjects\project\results\cells\%s_time%d_kstest.csv' % (d[sex], time),
        index=False)


def main():
    run_statistics_on_cells(4, 2)


if __name__ == '__main__':
    main()
