import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import linregress, spearmanr, mannwhitneyu
import matplotlib.pyplot as plt
import importlib
import math
import os
from scipy.stats.stats import pearsonr


def expected_change(file_name):
    """
    Calculates the LFC of the genes expressed as linear combinations of MTG single cell data,
    and calculates the linear regression representing the consistency with STG data
    :param file_name: Group; Females / Males
    :return: Creates figures presenting the expected change in MTG single cell data
    """
    file_name = file_name.replace('/', '_').replace(' ', '_')
    d = {'AD': 'ad', 'Control': 'control'}
    with open(r'single cell/all_data.p', "rb") as f:
        df = pickle.load(f)
    cells = [x.split('_')[0] for x in df.keys() if x != 'var']
    cells = sorted(list(set(cells)))

    chol_df = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\cholinergic_genes.csv')
    chol_genes = list(chol_df[chol_df['Role'] == 'Production']['Gene'])

    per = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\single cell\percentages.csv').set_index('Group')
    per = per.reindex(sorted(per.columns), axis=1)

    for sex in ['Female', 'Male']:
        all_genes = []
        for cond in ['Control', 'AD']:
            numbers = per[per.index == '%s %s' % (cond, sex)].values[0]
            weights = list(numbers / np.sum(numbers))
            gene_vals = []
            for i, gene in enumerate(chol_genes):
                expression = []
                for cell in cells:
                    expression.append(
                        np.expm1(df[f'{cell}_{sex}_ad'][d[cond]]['X'].mean(axis=0, dtype=np.float64)[i]) + 1e-9)
                gene_vals.append(np.sum([x * y for x, y in zip(weights, expression)]))
            all_genes.append(gene_vals)
        changes = [np.log2(x / y) for x, y in zip(all_genes[1], all_genes[0])]
        df2 = pd.DataFrame(changes, index=chol_genes, columns=['MTG'])
        c = {'Female': 'deeppink', 'Male': 'darkblue'}
        df1 = pd.read_csv(
            r'C:\Users\danas\PycharmProjects\project\single cell\consistent\%s_STG_changes_100.csv' % sex).set_index(
            'mRNA')
        comb = pd.concat([df1, df2], axis=1)
        comb = comb[comb.index.isin(list(df1.index))]
        sns.regplot(x=list(comb['STG']), y=list(comb['MTG']), color=c[sex])
        plt.ylim(-0.45, 0.45)
        plt.xlim(-0.55, 0.35)
        p = linregress(x=list(comb['STG']), y=list(comb['MTG']), alternative='greater').pvalue
        s = linregress(x=list(comb['STG']), y=list(comb['MTG']), alternative='greater').slope
        plt.savefig(
            r'C:\Users\danas\PycharmProjects\project\single cell\consistent\Figures\
            %s_MTG_changes_%s_p=%.3f_s=%.3f.png' % (sex, file_name, p, s))
        plt.clf()
        df2.to_csv(r'C:\Users\danas\PycharmProjects\project\single cell\consistent\
        %s_MTG_changes_%s.csv' % (sex, file_name))


def pearson():
    """
    Creates the figures that represents the changes in specific cell populations as a function of
    changes in cell amount and mean LFC of cholinergic genes
    :return: Figures representing the changes in specific cell populations
    """
    color_dict = {
        'L2/3 IT': 'yellowgreen',
        'L4 IT': 'cyan',
        'L5 IT': 'lightseagreen',
        'L6 IT': 'darkgoldenrod',
        'L6 IT Car3': 'mediumblue',
        'L6b': 'rebeccapurple',
        'L6 CT': 'steelblue',
        'L5 ET': 'darkcyan',
        'L5/6 NP': 'forestgreen',
        'Lamp5 Lhx6': 'brown',
        'Lamp5': 'lightcoral',
        'Sncg': 'violet',
        'Pax6': 'purple',
        'Vip': 'darkorchid',
        'Pvalb': 'red',
        'Sst': 'orange',
        'Chandelier': 'deeppink',
        'OPC': 'darkslategray',
        'Oligodendrocyte': 'slategray',
        'Sst Chodl': 'goldenrod',
        'Astrocyte': 'sienna',
        'Microglia-PVM': 'darkseagreen',
        'Endothelial': 'saddlebrown',
        'VLMC': 'olive'
    }
    d = {'AD': 'ad', 'Control': 'control'}
    with open(r'single cell/all_data.p', "rb") as f:
        df = pickle.load(f)
    cells = [x.split('_')[0] for x in df.keys() if x != 'var']
    cells = sorted(list(set(cells)))

    chol_genes = list(
        pd.read_csv(r'C:\Users\danas\PycharmProjects\project\single cell\consistent\Female_STG_changes_100.csv')[
            'mRNA'])
    non_neruons = ['Astrocyte', 'Oligodendrocyte', 'Microglia-PVM', 'VLMC', 'OPC', 'Endothelial']

    number_df = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\single cell\percentages.csv').set_index('Group')
    number_df = number_df.reindex(sorted(number_df.columns), axis=1)
    number_df.columns = cells

    for sex in ['Female', 'Male']:
        cor_vals = []
        num_vals = []
        sub_vals = []
        for cell in cells:
            if cell in non_neruons:
                sub_vals.append('Non-Neuron')
            else:
                sub_vals.append('Neuron')
            cell_vals = []
            number_vals = []
            for cond in ['Control', 'AD']:
                cond_cell_vals = []
                for i, gene in enumerate(chol_genes):
                    cond_cell_vals.append(
                            np.expm1(df[f'{cell}_{sex}_ad'][d[cond]]['X'].mean(axis=0, dtype=np.float64)[i]) + 1e-9)
                cell_vals.append(cond_cell_vals)
                number_vals.append(number_df[number_df.index == '%s %s' % (cond, sex)][cell].values[0])
            cor_vals.append(pearsonr(cell_vals[1], cell_vals[0])[0])
            num_vals.append(number_vals[1] / number_vals[0])
        new_df = pd.DataFrame(cor_vals, columns=['Pearson Coefficient'], index=cells)
        new_df['Change in the relative %\nof the number of cells'] = num_vals
        new_df['Subclass'] = sub_vals
        new_df['Color'] = cells
        new_df['Color'] = new_df['Color'].map(color_dict)
        sns.scatterplot(x=new_df['Change in the relative %\nof the number of cells'], y=new_df['Pearson Coefficient'],
                    s=200, c=new_df['Color'], style=new_df['Subclass'], legend=None)
        plt.ylim(0.984, 1.0005)
        plt.xlim(0.45, 2.38)
        plt.title(sex)
        plt.show()


def run_analysis(ct_interest):
    """
    Runs U-test on all genes in each distinct cell population, and returns the changes in the
    desired genes after FDR correction
    :param ct_interest: The specific cell populations
    :return: Creates files with the U-test results in each cell population from ct_interest
    """
    df = pd.read_csv(r'C:\Users\danas\PycharmProjects\project\data\cholinergic_genes.csv')
    genes_of_interest = list(df[df['Role'] == 'Production']['Gene'])
    adata = sc.read_h5ad('single cell/SEAAD_MTG_RNAseq.h5ad')
    address = "single cell/Utest/"
    for i, ct in enumerate(ct_interest):
        data_sbst = adata[adata.obs['Subclass'].isin([ct])]
        for sex in ['Female', 'Male']:
            data_sbst_sex = data_sbst[data_sbst.obs['Sex'].isin([sex])]
            data_sbst_sex = data_sbst_sex[data_sbst_sex.obs['Cognitive Status'].isin(['Dementia', 'No dementia'])]
            sc.tl.rank_genes_groups(data_sbst_sex, groupby='Cognitive Status', reference='No dementia',
                                    method='wilcoxon')
            result = data_sbst_sex.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            de_df = pd.DataFrame(
                {group + '_' + key[:1]: result[key][group]
                 for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
            de_df_interest = de_df[de_df['Dementia_n'].isin(genes_of_interest)]
            de_df_interest.to_csv(address + 'DE_wilcoxon_%s_%s.csv' % (sex, ct.replace('/', '_').replace(' ', '_')))
        print('finished for: %s' % ct)


def n_samples(ct_interest):
    """
    Calculates the number of cells from each cell population in each group (AD / Control)
    :param ct_interest: list of desired cell populations
    :return: Creates files with the number of cells from each cell population from ct_interest
    """
    adata = sc.read_h5ad('single cell/SEAAD_MTG_RNAseq.h5ad')
    address = "single cell/"
    n_samples_dic = {}
    for ct in ct_interest:
        data_sbst = adata[adata.obs['Subclass'].isin([ct])]
        n_samples_dic[ct] = []
        for sex in ['Female', 'Male']:
            data_sbst_sex = data_sbst[data_sbst.obs['Sex'].isin([sex])]
            data_sbst_sex = data_sbst_sex[data_sbst_sex.obs['Cognitive Status'].isin(['Dementia', 'No dementia'])]
            dem = np.sum(data_sbst_sex.obs['Cognitive Status'] == 'Dementia')
            non_dem = np.sum(data_sbst_sex.obs['Cognitive Status'] == 'No dementia')
            n_samples_dic[ct].append((dem, non_dem))
    df=pd.DataFrame(n_samples_dic)
    df.index=['Female', 'Male']
    df.to_csv(address + 'n_samples.csv')


def main():
    pearson()


if __name__ == '__main__':
    main()
