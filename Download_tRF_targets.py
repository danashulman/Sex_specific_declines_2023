import pandas as pd
import urllib.request
import numpy as np
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import os
import warnings
import multiprocessing as mp

PROCESS_NUMBER = mp.cpu_count()


def download_targets_for_all(tissue, prediction_threshold):
    """
    Gets the arguments and runs the inner loop in parallel as the number of the cores.
    :param tissue: The wanted tissue out of the 'tissue_options.csv' file. Can also be 'All' for no tissue filtration.
    :param prediction_threshold: Float describing the wanted threshold for prediction by DIANA.
    :return: file of the predicted targets and file of the Cholinergic targets above the given threshold of
    the tRFs in the list.
    """
    if not os.path.exists('tRF_targets'):
        os.makedirs('tRF_targets')
    tissue_path = 'tRF_targets\\%s' % tissue
    if not os.path.exists(tissue_path):
        os.makedirs(tissue_path)
    trfs = list(pd.read_csv(r'C:\Users\danas\PycharmProjects\tRF_targets\tRFs_of_interest.csv')['tRF'])
    new_trfs = []
    for trf in trfs:
        if not os.path.exists(tissue_path + '\\%s_%s.csv' % (trf, prediction_threshold)):
            new_trfs.append(trf)
    split = np.array_split(new_trfs, PROCESS_NUMBER)

    procs = []
    for i in range(PROCESS_NUMBER):
        procs.append(mp.Process(target=download_targets_for_sublist,
                                args=(list(split[i]), tissue,
                                      prediction_threshold)))
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join()


def download_targets_for_sublist(trf_sublist, tissue, prediction_threshold):
    """
    Runs the loop for the single tRF on a sublist of tRFs, given for parallel running.
    :param trf_sublist: A list holding a subgroup of tRFs.
    :param tissue: The tissue that was given as an argument.
    :param prediction_threshold: Float describing the wanted threshold for prediction by DIANA.
    :return: file of the predicted targets and file of the Cholinergic targets above the given threshold of
    the tRFs in the list.
    """
    meta = pd.read_csv(r'C:\Users\danas\PycharmProjects\tRF_targets\tRF_meta.csv')
    meta = meta[meta['tRF'].isin(trf_sublist)].reset_index().iloc[:, 1:]
    transcript_coding = pd.read_csv('transcript_to_gene_ID.txt', sep='\t')
    gene_coding = pd.read_csv('gene_ID_to_symbol.csv')
    transcript_to_gene = pd.Series(transcript_coding['Gene stable ID'].values,
                                   index=transcript_coding['Transcript stable ID']).to_dict()
    id_to_symbol = pd.Series(gene_coding['external_gene_name'].values, index=gene_coding['X']).to_dict()

    for trf in trf_sublist:
        try:
            options = webdriver.ChromeOptions()
            options.add_argument("headless")
            driver = webdriver.Chrome(chrome_options=options)
            diana_url = 'http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=mrmicrot/index'
            download_targets_for_tRF(diana_url, driver, id_to_symbol, meta, prediction_threshold,
                                     transcript_to_gene, trf, tissue)
        except Exception as e:
            print(str(e))
            continue


def download_targets_for_tRF(diana_url, driver, id_to_symbol, meta, prediction_threshold,
                             transcript_to_gene, trf, tissue):
    """
    Download the targets for a single tRF.
    :param diana_url: The DIANA url for each tRFs.
    :param driver: The driver object for downloading the targets from DIANA.
    :param id_to_symbol: A dictionary allowing the transition from the gene ID to the gene symbol.
    :param meta: The updated meta for the tRFs.
    :param prediction_threshold: Float describing the wanted threshold for prediction by DIANA.
    :param transcript_to_gene: A dictionary allowing the transition from the transcript ID to the gene ID.
    :param trf: The current tRF.
    :param tissue: The tissue that was given as an argument.
    :return:
    """
    seq = list(meta[meta['tRF'] == trf]['tRF sequence'])[0]
    driver.get(diana_url)
    textarea = driver.find_element_by_tag_name("textarea")
    textarea.send_keys(seq)
    driver.find_element_by_name("yt0").click()
    timeout = 240
    try:
        element_present = EC.presence_of_element_located((By.ID, 'mrmicrot_result_download_button'))
        WebDriverWait(driver, timeout).until(element_present)
    except Exception:
        raise Exception("Time out for %s" % trf)
    try:
        trf_url = 'http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=mrmicrot/download&species=human&seq=%s' \
                  % seq
        response = urllib.request.urlopen(trf_url)
        web_content = response.read()
        web_string = web_content.decode('utf8')
        web_string = web_string.split('\n')
        new_lines = []
        for line in web_string:
            if line.startswith('>'):
                new_lines.append(line.split('|')[1:])
            else:
                continue
        df = pd.DataFrame(new_lines, columns=['Gene', 'Score'])

        df2 = df[df['Score'].astype('float') >= float(prediction_threshold)].reset_index().iloc[:, 1:]
        df2 = df2.replace({'Gene': transcript_to_gene})
        df2 = df2.replace({'Gene': id_to_symbol})
        df2.to_csv('tRF_targets\\%s\\%s_%s.csv' % (tissue, trf, prediction_threshold), index=False)
        print("Finished download for %s" % trf)

    except Exception:
        raise Exception("Couldn't reach the web for %s" % trf)


def main():
    tissue, prediction_threshold = "cells", "0.8"
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    download_targets_for_all(tissue, prediction_threshold)


if __name__ == '__main__':
    main()
