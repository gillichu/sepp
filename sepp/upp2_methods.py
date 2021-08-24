from datetime import timedelta
import glob
import os
import subprocess
from timeit import default_timer

import numpy as np

from sepp.filemgr import get_root_temp_dir
from sepp.hmm_concurrent import *
from sepp.scheduler import JobPool


_LOG = get_logger(__name__)

def create_dirs(path_name):
    if not os.path.exists(path_name):
        subprocess.call(['mkdir', path_name])


def parse_strats(stratsfile):
    stratlist = []
    with open(stratsfile, "r") as reader:
        lines = reader.readlines()
        for line in lines:
            stratlist.append(line.strip())
    return stratlist


def makedirstruct(dirpath):
    print("[makedirstruct]")
    # clear the folders
    dirpath = os.path.join(dirpath, "ensembleData")
    create_dirs(dirpath)
    for firstlvl in ["trueAlignment","temporaryFileSave","subsetTrueAln", "Searcher","queryToHmm","newHMM","ML","initialHMM","hmmSeqAlign","hmmScores","hmmQueryList","fullPrediction", "seqFileNames"]:
        firstdir = '%s/%s' % (dirpath, firstlvl)
        if os.path.isdir(firstdir):
            subprocess.call(['rm', '-rf', firstdir])

        subprocess.call(['mkdir', firstdir])
        if firstlvl == 'fullPrediction':
            subprocess.call(['mkdir', firstdir + '/sequences/'])
        elif firstlvl == 'hmmQueryList':
            for z in ['inputQuery', 'merged', 'predictedQuery']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'hmmScores':
            for z in ['newRawOld', 'processedOld', 'rawOld', 'scoresFull', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'ML':
            for z in ['hmm', 'models', 'queryHMM', 'scores', 'sequences', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
            subprocess.call(['mkdir', firstdir + '/sequences/np'])
        elif firstlvl == 'newHMM':
            for z in ['columnSets', 'hmm', 'newHMMseq']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'queryToHmm':
            for z in ['original', 'withUPP']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'Searcher':
            subprocess.call(['mkdir', '%s/scoreFiles' % firstdir])
        elif firstlvl == 'trueAlignment':
            for z in ['original', 'subset']:
                subprocess.call(['mkdir', firstdir + '/' + z])


def run_upp_strats(abstract_algorithm, dirname, hier_upp, adjusted_bitscore, early_stop, doResort=False):
    ''' Note: abstract_algorithm can be used to do commands like below.
        The moment you call addHMMBuildJob, jobs will get enqueued and run eventually
    addHMMBuildJob(abstract_algorithm, <hmmbuild profile output file path>, <fasta file path>)
    JobPool().wait_for_all_jobs()
    '''
    preprocessing_start_time = default_timer()
    _LOG.info("[run_upp_strats]")


    # globdir = glob.glob(os.path.join(dirname, "output*"))
    # assert len(globdir) == 1
    # outputdirname = globdir[0]
    root_temp_dir = get_root_temp_dir()
    hmm_folder_prefix = "%s/root/P_0/" % root_temp_dir
    # Nit: We should be able to retrieve the same information from options().fragment_file
    fragment_file_suffix = os.listdir(root_temp_dir + "/fragment_chunks/")[0]
    fragment_file = "%s/fragment_chunks/%s" % (root_temp_dir, fragment_file_suffix)
    # predictionName = './%s/UPPoutput/%s_output_alignment.fasta' % dirpath

    # TODO: Prediction name may get passed in
    prediction_name = ""
    # TODO: this is something we need to remove
    true_alignment = ""
    dataset_name = "default-value-not-empty"

    set_all_file_names(hmm_folder_prefix, fragment_file, true_alignment, prediction_name, root_temp_dir, dataset_name)
    save_initial_steps(abstract_algorithm)
    preprocessing_end_time = default_timer()
    time_output_suffix = "upp2.time"
    with open(abstract_algorithm.get_output_filename(time_output_suffix), "a+") as f:
        f.write(f"2-preprocessing: {preprocessing_end_time - preprocessing_start_time} seconds ({timedelta(seconds=preprocessing_end_time - preprocessing_start_time)})\n")

    search_and_align_start_time = default_timer()
    hierarchy_search_results = {
        "num_hmms_looked_at": None,
        "per_query_hmms_looked_at": None,
        "num_hmms": None,
        "per_query_hmm_chosen": None,
    }
    strat = None
    if hier_upp:
        strat = "hierarchical"
        hierarchy_search_results = hierchySearch(abstract_algorithm, adjusted_bitscore, early_stop)
    elif adjusted_bitscore:
        # this forces UPP2 to load fullAdjusted.npy instead of score.np
        strat = "adjusted_bitscore"
    else:
        strat = "upp"

    _LOG.info("[processing %s]" % strat)
    _LOG.info("[running scoresToHMMSeq]")
    scoresToHMMSeq(abstract_algorithm, strat)
    _LOG.info("[running buildAlignMerge, doResort is %s]" % doResort)
    buildAlignMerge(abstract_algorithm, strat, doResort=doResort)
    search_and_align_end_time = default_timer()

    num_hmms_output_suffix = "upp2.num_hmms"
    with open(abstract_algorithm.get_output_filename(time_output_suffix), "a+") as f:
        f.write(f"3-search and align: {search_and_align_end_time - search_and_align_start_time} seconds ({timedelta(seconds=search_and_align_end_time - search_and_align_start_time)})")
    if (hier_upp):
        with open(abstract_algorithm.get_output_filename(num_hmms_output_suffix), "w") as f:
            f.write(f"num hmms looked at: {hierarchy_search_results['num_hmms_looked_at']}\n")
            f.write(f"total num hmms: {hierarchy_search_results['num_hmms']}\n")
            f.write(f"per query index hmm chosen: {list(hierarchy_search_results['per_query_hmm_chosen'])}\n")
            f.write(f"per query num hmm looked at: {list(hierarchy_search_results['per_query_hmms_looked_at'])}\n")
            f.write(f"per query avg num hmm looked at: {np.average(hierarchy_search_results['per_query_hmms_looked_at'])}\n")

    _LOG.info(f"Time written to {abstract_algorithm.get_output_filename(time_output_suffix)}")
    _LOG.info(f"Num hmms written to {abstract_algorithm.get_output_filename(num_hmms_output_suffix)}")
