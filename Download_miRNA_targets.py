import pandas as pd
import urllib.request
import numpy as np
import os
import warnings
import multiprocessing as mp

# PROCESS_NUMBER = mp.cpu_count()
PROCESS_NUMBER = 6


def download_targets_for_all(tissue):
    """
    Gets the arguments and runs the inner loop in parallel as the number of the cores.
    :param tissue: The wanted tissue out of the 'tissue_options.csv' file. Can also be 'All' for no tissue filtration.
    :param prediction_threshold: Float describing the wanted threshold for prediction by DIANA.
    :return: file of the predicted targets and file of the Cholinergic targets above the given threshold of
    the tRFs in the list.
    """
    if not os.path.exists('miRNA_targets'):
        os.makedirs('miRNA_targets')
    tissue_path = 'miRNA_targets\\%s' % tissue
    if not os.path.exists(tissue_path):
        os.makedirs(tissue_path)
    mirs = list(pd.read_csv(r'C:\Users\danas\PycharmProjects\tRF_targets\miRNAs_of_interest.csv')['miRNA'])
    print(len(mirs))
    new_mirs = []
    for mir in mirs:
        if not os.path.exists(tissue_path + '\\%s.csv' % mir):
            new_mirs.append(mir)
    print(len(new_mirs))

    split = np.array_split(new_mirs, PROCESS_NUMBER)
    procs = []
    for i in range(PROCESS_NUMBER):
        procs.append(mp.Process(target=download_targets_for_sublist,
                                args=(list(split[i]), tissue)))
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join()


def download_targets_for_sublist(mir_sublist, tissue):
    """
    Runs the loop for the single tRF on a sublist of tRFs, given for parallel running.
    :param mir_sublist: A list holding a subgroup of tRFs.
    :param tissue: The tissue that was given as an argument.
    :param prediction_threshold: Float describing the wanted threshold for prediction by DIANA.
    :return: file of the predicted targets and file of the Cholinergic targets above the given threshold of
    the tRFs in the list.
    """
    for mir in mir_sublist:
        try:
            download_targets_for_miR(mir, tissue)
        except Exception as e:
            print(str(e))
            continue


def download_targets_for_miR(mir, tissue):
    try:
        mir_url = \
            'https://dianalab.e-ce.uth.gr/html/dianauniverse/index.php?r=download/microT_' \
            'CDS&keywords=%sp&genes=&mirnas=%s' % (mir, mir) + '%20&descr=&threshold=0.8'
        response = urllib.request.urlopen(mir_url)
        web_content = response.read()
        web_string = web_content.decode('utf8')
        web_string = web_string.split('\n')
        new_lines = []
        for line in web_string:
            if line.startswith('ENST'):
                split_list = line.split(',')
                new_lines.append([split_list[1][17:-1], split_list[3]])
            else:
                continue
        df = pd.DataFrame(new_lines, columns=['Gene', 'Score'])
        df.to_csv(r'miRNA_targets\\%s\\%s.csv' % (tissue, mir), index=False)
        print('Finished for %s' % mir)

    except Exception:
        raise Exception("Couldn't reach the web for %s" % mir)


def main():
    tissue = "cells"
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    download_targets_for_all(tissue)


if __name__ == '__main__':
    main()
