import copy
from datetime import timedelta
import random
import sys
from timeit import default_timer

import numpy as np
import os
import pathlib
from subprocess import run

from Bio import SeqIO

from sepp import get_logger
from sepp.filemgr import get_root_temp_dir
from sepp.scheduler import JobPool
from sepp.jobs import HMMAlignJob,HMMBuildJob,HMMSearchJob


hmmSeqFile = ''
queryName = ''
trueAlignment = ''
predictionName = ''
dirName = ''
dataSetName = ''

_LOG = get_logger(__name__)

def ensureFolder(fileName):
    pathlib.Path(fileName).parents[0].mkdir(parents=True, exist_ok=True)
    # def removeEnd(fileName):
    #     fileName = '/'.join(fileName.split('/')[:-1])
    #     return (fileName)
    # fileName = removeEnd(fileName)
    # if os.path.exists(fileName) == False:
    #     toMake = []
    #     while os.path.exists(fileName) == False:
    #         toMake.append(fileName)
    #         fileName = removeEnd(fileName)
    #     toMake = np.array(toMake)[-1::-1]
    #     for a in range(len(toMake)):
    #         os.mkdir(toMake[a])

def seqToArray(seqs):
    seqArray = []
    for seq in seqs:
        seqArray.append(list(seq))
    seqArray = np.array(seqArray)
    return seqArray

def arrayToSeq(array):
    seqs = []
    for a in range(len(array)):
        seq = ''.join(list(array[a]))
        seqs.append(seq)
    seqs = np.array(seqs)
    return seqs

def set_all_file_names(hmm_folder_prefix, fragment_file, true_alignment, prediction_name, root_temp_dir, dataset_name):
    global hmmSeqFile
    hmmSeqFile = hmm_folder_prefix
    global queryName
    queryName = fragment_file
    # TODO: this is something we need to remove
    global trueAlignment
    trueAlignment = true_alignment
    global predictionName
    predictionName = prediction_name
    global dirName
    dirName = root_temp_dir
    global dataSetName
    dataSetName = dataset_name

def giveAllFileNames():
    global hmmSeqFile
    global queryName
    global trueAlignment
    global predictionName
    global dataSetName

    fileNames = [hmmSeqFile, queryName, trueAlignment, predictionName, dataSetName]
    return fileNames

def giveQueries():
    ''' This file returns every other line without the first character
    starting with the first line. This will return the sequence names in some
    fasta files but will break for other fasta files. Some reasons for why I believe this is true below.
    1. Comments: semicolons, although less common, can appear in fasta files
    2. New lines: There can be any number of new lines to my knowledge inside a fasta files between
    the end of one sequnece and the beginning of next
    3. Interleaved vs sequential: A very common fasta file representation actually puts fasta sequneces
    on multiple lines limited to 80-120 characters. Programs commonly output their sequences in this format
    as well.
    The reason why this code correctly is because upp outputs fragments in a sequential format with
    no additional new lines or comments and seemingly this function is only called on UPP outputs.
    Not a huge concern but it's something we should consider changing since there are much more robust ways of
    doing this through libraries
    '''
    queryName = giveQueryFileName()
    file_obj = open(queryName)
    queries = []
    count1 = 0
    for line in file_obj:
        if line[0] == ">":
            queries.append(line[1:-1])
    queries = np.array(queries)
    return queries

def removeMultInvert(seqs):
    lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))
    seqArray = seqToArray(seqs)
    seqArray = seqArray.T
    seqArray2 = []
    for a in range(0, len(seqArray)):
        seq = seqArray[a]
        argsLower = np.argwhere(np.isin(seq, lowercase))[:, 0]
        if argsLower.shape[0] == 0:
            seqArray2.append(np.copy(seqArray[a]))
        else:
            for b in range(len(argsLower)):
                seq2 = np.zeros(len(seq)).astype(str)
                seq2[:] = '-'
                seq2[argsLower[b]] = seq[argsLower[b]]
                seqArray2.append(np.copy(seq2))
    seqArray2 = np.array(seqArray2).T
    seqs = arrayToSeq(seqArray2)
    return seqs

def removeEmptyColumns(seqs):
    seqArray = []
    for seq in seqs:
        seqArray.append(copy.copy(list(seq)))
    seqArray = np.array(seqArray).T
    seqArray2 = []
    for a in range(0, len(seqArray)):
        size1 = seqArray[a][seqArray[a] != '-'].shape[0]
        if size1 != 0:
            seqArray2.append(np.copy(seqArray[a]))
    seqArray2 = np.array(seqArray2).T
    seq2 = []
    for a in range(0, len(seqArray2)):
        str1 = ''.join(list(seqArray2[a]))
        seq2.append(str1)
    seq2 = np.array(seq2)
    return seq2

def saveFasta(name, formatedData):
    ensureFolder(name)
    with open(name, 'w') as f:
        f.write('\n'.join(formatedData))

def loadFastaFormat(name):
    file_obj = open(name)
    data = []
    name = ''
    for line_number, i in enumerate(file_obj):
        name = name + i
        if line_number % 2 == 1:
            name = name[:-1]
            data.append(name)
            name = ''
    return data

def loadFastaBasic(name):
    file_obj = open(name)
    keys = []
    seqs = []
    count1 = 0
    for line_number, i in enumerate(file_obj):
        if line_number % 2 == 0:
            keys.append(i[1:-1])
        else:
            seq1 = i.replace('\n', '')
            seqs.append(seq1)
    keys, seqs = np.array(keys), np.array(seqs)
    return keys, seqs

def saveFastaBasic(name, keys, seqs):
    data = []
    for a in range(0, len(keys)):
        data.append('>' + keys[a] + '\n' + seqs[a])
    data[-1] = data[-1] + '\n'
    saveFasta(name, data)

def save_sequence_filenames():
    fileNames = giveAllFileNames()
    hmm_folder_prefix = fileNames[0]
    sequenceFileNames = []
    Anums = []
    namesA = os.listdir(hmm_folder_prefix)
    for nameA in namesA:
        if nameA[:2] == 'A_':
            Anum = nameA.split('_')
            Anum = int(Anum[-1])
            src1 = hmm_folder_prefix + nameA + '/'
            names = os.listdir(src1)
            for name in names:
                if name[-5:] == 'fasta':
                    src2 = src1 + name
                    sequenceFileNames.append(src2)
                    Anums.append(Anum)
    # this contains all the .fasta files that were found by using
    # hmm_folder_prefix/A_*/*.fasta
    sequenceFileNames = np.array(sequenceFileNames)
    # from here to 28 lines below this line does deduplication
    # if the length of one fasta file equals the length of another file, then they are never going to
    # have intersections so np.intersect1d(keys1, keys2) == len(keys1) would always be false if upp decomposition
    # works properly. However, UPP decomposition, after looking at the files manually, is not exactly disjoint at the leaf level
    # why is upp decomposition like this?
    Anums = np.array(Anums)
    
    removalArgs = []

    numSeq = []
    firstKey = []
    for a in range(len(sequenceFileNames)):
        keys1, _ = loadFastaBasic(sequenceFileNames[a])
        keys1 = np.sort(np.array(keys1))
        numSeq.append(len(keys1))
        firstKey.append(keys1[0])


    for a in range(len(sequenceFileNames)):
        for b in range(len(sequenceFileNames)):
            if b > a:
                if numSeq[a] == numSeq[b]:
                    if firstKey[a] == firstKey[b]:
                        keys1, _ = loadFastaBasic(sequenceFileNames[a])
                        keys2, _ = loadFastaBasic(sequenceFileNames[b])
                        if len(keys1) == len(keys2):
                            if len(np.intersect1d(keys1, keys2)) == len(keys1):
                                removalArgs.append(b)
        
        
    keepArgs = np.arange(len(sequenceFileNames))
    # one effect i observed was that having removalArgs makes it so that
    # further down the edges and children no longer consider themselves their own parents
    # not sure why this happpens in the first place since i feel like it shouldn't
    removalArgs = np.array(removalArgs)
    keepArgs = keepArgs[np.isin(keepArgs, removalArgs) == False]

    Anums = Anums[keepArgs]
    sequenceFileNames = sequenceFileNames[keepArgs]
    sequenceFileNames = sequenceFileNames[np.argsort(Anums)]

    inputFiles = giveAllFileNames()
    dataFileName = inputFiles[4]

    sequenceFileNames_saveName = get_root_temp_dir() + '/data/internalData/' + dataFileName + '/seqFileNames/names.npy'
    ensureFolder(sequenceFileNames_saveName)
    np.save(sequenceFileNames_saveName, sequenceFileNames)

def giveSequenceFileNames():
    ''' This function returns all the sequence file names inside root_temp_dir/A_*/*.fasta
    '''
    inputFiles = giveAllFileNames()
    dataFileName = inputFiles[4]
    sequenceFileNames_saveName = get_root_temp_dir() + '/data/internalData/' + dataFileName + '/seqFileNames/names.npy'
    sequenceFileNames = np.load(sequenceFileNames_saveName)
    return sequenceFileNames

def giveQueryFileName():
    fileNames = giveAllFileNames()
    queryName = fileNames[1]
    return queryName

def save_initial_hmm(abstract_algorithm):
    ''' This function builds hmms on every single input inside sequenceFileNames
    which is equivalent to building hmms on every single fasta file in root_temp_dir/A_*/*.fasta
    '''
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        src = sequenceFileNames[a]
        inputFiles = giveAllFileNames()
        dataFileName = inputFiles[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFileName + "/initialHMM/test" + str(a) + ".hmm"
        ensureFolder(hmmName)
        addHMMBuildJob(abstract_algorithm, hmmName, src)
    JobPool().wait_for_all_jobs()

def giveHMMversion():
    hmmVersion = ''
    return hmmVersion

def addHMMSearchJob(abstract_algorithm, hmmName, queryName, scoreName):
    ensureFolder(scoreName)
    hmmsearch_job = HMMSearchJob(tblout=True, **vars(abstract_algorithm.options.hmmsearch))
    hmmsearch_job.setup(hmmName, queryName, scoreName, filters=False, elim=abstract_algorithm.elim) #, trim=False)
    hmmsearch_job.results_on_temp = False
    JobPool().enqueue_job(hmmsearch_job)

def addHMMBuildJob(abstract_algorithm, hmmName, seqName):
    ensureFolder(hmmName)
    hmmbuild_job = HMMBuildJob()
    if(abstract_algorithm.options.molecule is not None):
        hmmbuild_job.setup(seqName, hmmName, molecule=abstract_algorithm.options.molecule, **vars(abstract_algorithm.options.hmmbuild))
    else:
        hmmbuild_job.setup(seqName, hmmName, molecule=abstract_algorithm.options.molecule, **vars(abstract_algorithm.options.hmmbuild))
    JobPool().enqueue_job(hmmbuild_job)

def addHMMAlignJob(abstract_algorithm, hmmName, queryName, predictionName):
    ensureFolder(predictionName)
    hmmalign_job = HMMAlignJob(**vars(abstract_algorithm.options.hmmalign))
    # _LOG.debug(str(abstract_algorithm.options.hmmalign))
    if(abstract_algorithm.options.molecule is not None):
        hmmalign_job.setup(hmmName, queryName, predictionName, molecule=abstract_algorithm.options.molecule, trim=False)
    else:
        hmmalign_job.setup(hmmName, queryName, predictionName, trim=False)
    JobPool().enqueue_job(hmmalign_job)

def runHMMbuild(abstract_algorithm, hmmName, seqName):
    # ensureFolder(hmmName)
    hmmbuild_job = HMMBuildJob()
    hmmbuild_job.setup(seqName, hmmName, **vars(abstract_algorithm.options.hmmbuild))
    hmmbuild_job.run()

def runHMMalign(abstract_algorithm, hmmName, queryName, predictionName):
    # ensureFolder(predictionName)
    hmmalign_job = HMMAlignJob()
    hmmalign_job.setup(hmmName, queryName, predictionName, **vars(abstract_algorithm.options.hmmalign), trim=False)
    hmmalign_job.run()

def saveScoreSimple(abstract_algorithm, hmmNames, queryNames, scoreName):
    ''' This function saves to scoreName a 2d numpy array where the rows define the
    query sequences ordered by the order inside giveQueries which is the order inside the actual fasta file.
    The columns are the hmm files ordered by the array hmmNames which is the order of fasta files returned
    whnen calling giveSequenceFileNames. The values are the hmmsearch scores.
    '''
    queryDict = {}
    queries = giveQueries()
    for queryNum in range(0, len(queries)):
        query = queries[queryNum]
        queryDict[query] = queryNum
    sequenceFileNames = giveSequenceFileNames()
    hmmNum = len(hmmNames)
    data = np.zeros((len(queries), hmmNum))
    randIntArr = random.sample(range(100000000000), len(hmmNames))
    for a in range(0, len(hmmNames)):
        hmmName = hmmNames[a]
        scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + str(randIntArr[a]) + ".txt"
        queryName = queryNames[a]
        ensureFolder(scoreNameTemp)
        # BUG: this should be hmmmsearch I'm pretty sure
        addHMMSearchJob(abstract_algorithm, hmmName, queryName, scoreNameTemp)
    JobPool().wait_for_all_jobs()

    for a in range(0, len(hmmNames)):
        hmmName = hmmNames[a]
        scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + str(randIntArr[a]) + ".txt"
        queryName = queryNames[a]
        ensureFolder(scoreNameTemp)
        dataOld = []
        file_obj = open(scoreNameTemp)
        count1 = 0
        for line_number, i in enumerate(file_obj):
            ar = i.split('  ')
            listI = []
            bitscoreList = []
            if line_number == 1:
                posScore = i.find('score')
            if line_number == 2:
                argSpace = np.argwhere(np.array(list(i)) == ' ')[:, 0]
            if (len(ar) >= 30) and (line_number >= 3):
                # sequenceName = i[:20].replace(' ', '')
                sequenceName = ""
                for current_char in i:
                    if current_char == " " or current_char == "\t":
                        break
                    sequenceName += current_char

                bitScore = float(i[posScore-1:posScore+5].replace(' ', ''))
                dataOld.append([sequenceName, bitScore])

        dataOld = np.array(dataOld)
        for b in range(0, len(dataOld)):
            queryN = queryDict[dataOld[b, 0]]
            data[queryN, a] = float(dataOld[b, 1])
    ensureFolder(scoreName)
    np.save(scoreName, data)
    # os.remove(scoreNameTemp)

def saveNewScores(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    queryNameO = giveQueryFileName()
    queryData = loadFastaFormat(queryNameO)
    queryNames = []
    hmmNames = []

    queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    numHMM = int(np.max(queryToHmm[:, 1])+1)
    for b in range(0, numHMM):
        argsHMM = queryToHmm[:, 0][queryToHmm[:, 1] == b]
        dataNow = []
        keysNow = []
        for a in argsHMM:
             dataNow.append(queryData[a])
        dataFolderName = giveAllFileNames()[4]
        queryName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/inputQuery/'  + str(b) + '.fasta'
        ensureFolder(queryName)
        saveFasta(queryName, dataNow)
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(b) + ".hmm"
        ensureFolder(hmmName)
        queryNames.append(queryName)
        hmmNames.append(hmmName)
    dataFolderName = giveAllFileNames()[4]
    scoreName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/hmmScores/scoresFull/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(abstract_algorithm, hmmNames, queryNames, scoreName)

def saveScoreFromBool(abstract_algorithm, hmmNames, queryHMM, scoreName, noSave=False):
    queryDict = {}
    queries = giveQueries()
    for queryNum in range(0, len(queries)):
        query = queries[queryNum]
        queryDict[query] = queryNum
    sequenceFileNames = giveSequenceFileNames()
    hmmNum = len(hmmNames)
    data = np.zeros((len(queries), hmmNum))
    randInt1Arr = random.sample(range(100000000000), len(hmmNames))
    randInt2Arr = random.sample(range(100000000000), len(hmmNames))
    for a in range(0, len(hmmNames)):
        argsInclude = queryHMM[:, 0][queryHMM[:, 1] == a]
        if argsInclude.shape[0] > 0:
            hmmName = hmmNames[a]
            scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + str(randInt1Arr[a]) + ".txt"
            queryName = giveQueryFileName()
            queryData = loadFastaFormat(queryName)
            queryNameTemp = get_root_temp_dir() + '/data/temporaryStorage/temp_' + str(randInt2Arr[a]) + '.fasta'
            queryNow = []
            for b in argsInclude:
                queryNow.append(queryData[b])
            ensureFolder(queryNameTemp)
            saveFasta(queryNameTemp, queryNow)
            ensureFolder(scoreNameTemp)
            # BUG: this should be hmmmsearch I'm pretty sure
            addHMMSearchJob(abstract_algorithm, hmmName, queryNameTemp, scoreNameTemp)
    JobPool().wait_for_all_jobs()
    _LOG.debug(f"saveScoreFromBool: number of hmms - {len(hmmNames)}")
    for a in range(0, len(hmmNames)):
        argsInclude = queryHMM[:, 0][queryHMM[:, 1] == a]
        if argsInclude.shape[0] > 0:
            scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + str(randInt1Arr[a]) + ".txt"
            dataOld = []
            file_obj = open(scoreNameTemp)
            count1 = 0
            for line_number, i in enumerate(file_obj):
                # make sure this is literally two white spcae characters because otherwise it won't work
                # it also hinges on the assumption that the hmmsearch output is the --tblout format and not the default -o
                ar = i.split("  ")
                listI = []
                bitscoreList = []
                if line_number == 1:
                    posScore = i.find('score')
                if line_number == 2:
                    argSpace = np.argwhere(np.array(list(i)) == ' ')[:, 0]
                if (len(ar) >= 30) and (line_number >= 3):
                    sequenceName = i[:20].replace(' ', '')
                    # try:
                    bitScore = float(i[posScore-1:posScore+5].replace(' ', ''))
                    # except:
                    #     bitScore = 42112358
                    dataOld.append([sequenceName, bitScore])
            dataOld = np.array(dataOld)
            for b in range(0, len(dataOld)):
                queryN = queryDict[dataOld[b, 0]]
                data[queryN, a] = float(dataOld[b, 1])
    # this remove no longermakes sense due to the name changing in every iteration
    # os.remove(scoreNameTemp)
    # os.remove(queryNameTemp)
    if noSave:
        return data
    else:
        ensureFolder(scoreName)
        np.save(scoreName, data)

def save_scores_original(abstract_algorithm):
    ''' This function saves the results of hmmsearch bitscores in full.npy.
    see saveScoreSimple for more information (rows = query sequences, columns = hmm files)
    '''
    sequenceFileNames = giveSequenceFileNames()
    hmmNames = []
    queryNames = []

    for a in range(0, len(sequenceFileNames)):
        dataFolderName = giveAllFileNames()[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/initialHMM/test" + str(a) + ".hmm"
        hmmNames.append(hmmName)
        queryName = giveQueryFileName()
        queryNames.append(queryName)

    dataFolderName = giveAllFileNames()[4]
    scoreName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(abstract_algorithm, hmmNames, queryNames, scoreName)

def processScores():
    ''' Note: is this deprecated
    '''
    sequenceFileNames = giveSequenceFileNames()
    hmmNames = []
    queryNames = []

    for a in range(0, len(sequenceFileNames)):
        dataFolderName = giveAllFileNames()[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/initialHMM/test" + str(a) + ".hmm"
        hmmNames.append(hmmName)
        queryName = giveQueryFileName()
        queryNames.append(queryName)

    dataFolderName = giveAllFileNames()[4]
    scoreName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(hmmNames, queryNames, scoreName)

def processScoresOld():
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        data = []
        file_obj = open("%s/hmmScores/raw/" % dirName + str(a) + ".txt")
        count1 = 0
        for line_number, i in enumerate(file_obj):
            ar = i.split('  ')
            if (len(ar) >= 30) and (line_number >= 3):
                bitScore = float(i[76:83].replace(' ', ''))
                sequenceName = i[:20].replace(' ', '')
                data.append([sequenceName, bitScore])
        data = np.array(data)
        np.save("%s/hmmScores/processed/" %  dirName + str(a) + ".npy", data)

def saveScoresBySeq():
    queryDict = {}
    queries = giveQueries()
    for queryNum in range(0, len(queries)):
        query = queries[queryNum]
        queryDict[query] = queryNum
    sequenceFileNames = giveSequenceFileNames()
    hmmNum = len(sequenceFileNames)
    data = np.zeros((len(queries), hmmNum))
    for a in range(0, hmmNum):
        dataFolderName = giveAllFileNames()[4]
        dataOld = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/processed/" + str(a) + ".npy")
        for b in range(0, len(dataOld)):
            queryN = queryDict[dataOld[b, 0]]
            data[queryN, a] = float(dataOld[b, 1])

    dataFolderName = giveAllFileNames()[4]
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy", data)

def save_decomposition():
    ''' This functions stores into treeDecomp.npy which is a numpy array of tree decompositions where
    treeDecomp[a,b] = 1 if a is a parent of b in the rooted representation of the decomposition tree
    hmmSets.npy stores the sequence names for each hmm in the order defined by giveSequenceFileNames
    which is the same order as the HMMs in most functions.
    '''
    dataFolderName = giveAllFileNames()[4]
    sets = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)



    sizes = []
    for a in range(0, len(sets)):
        sizes.append(a)
    sizes = np.array(sizes).astype(int)

    sizes_argsort = np.argsort(sizes)
    argsort_inverse = np.argsort(sizes_argsort)

    parentNode = np.zeros(len(sets)).astype(int)


    for a in range(1, len(sets)):
        currentParent = 0
        for b in range(1, a):
            if parentNode[b] == currentParent:
                set1, set2 = np.array(sets[sizes_argsort[a]]), np.array(sets[sizes_argsort[b]])
                if np.intersect1d(set1, set2).shape[0] == set1.shape[0]:
                    currentParent = b
        parentNode[a] = currentParent

    treeData = np.zeros((len(sets), len(sets)))
    for a1 in range(0, len(sets)):
        a = int(len(sets) - 1 - a1)
        b = parentNode[a]
        treeData[a, a] = 1
        treeData[b] = treeData[b] + treeData[a]

    treeData[treeData > 1] = 1

    treeData = treeData[argsort_inverse, :]
    treeData = treeData[:, argsort_inverse]

    
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/treeDecomp.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/treeDecomp.npy", treeData)

def save_hmm_sets():
    ''' This function saves a 2d array whose rows are the sequence file names order which is
    also the same order as the HMMs in a lot of functions. The columns are variable but store
    the sequence names for each sequence file in a row. Each row probbaly contains different numbers of
    columns.
    '''
    sequenceFileNames = giveSequenceFileNames()
    sets = []
    for a in range(0, len(sequenceFileNames)):
        # print(f"{a}: {sequenceFileNames[a]}")
        src = sequenceFileNames[a]
        set1 = []
        file_obj = open(src)
        count1 = 0
        for line_number, i in enumerate(file_obj):
            if line_number % 2 == 0:
                set1.append(i[1:-1])

        sets.append(set1)


    #print ("Q")

    #It gives warning. Fine for now but maybe change at some point.

    dataFolderName = giveAllFileNames()[4]

    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", sets)

def findDecomposition():
    dataFolderName = giveAllFileNames()[4]
    treeData = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/treeDecomp.npy")
    return treeData

def findTreeEdges():
    ''' This function returns a 2darray where the array[a,b] = 1
    if there is an edge from a to b. From what I can tell,
    argMax is the node with the most children so right now
    and argMin is the node with the least number of children from the
    list of argMax's parents. There essentially is an edge iteratively constructed
    from selecting whoever has the most children right now so whoever is the highest up in the tree
    to whoever has the least children right now from the parents of the first selection. In the second search,
    we don't exclude nodes we have already visited.
    '''
    treeData = findDecomposition()
    treeData[np.arange(treeData.shape[0]), np.arange(treeData.shape[0])] = 0
    edges = np.zeros(treeData.shape)
    # this would be counting the number of children each HMM has
    treeSum = np.sum(treeData, axis=1)
    nodes = []
    for a in range(treeData.shape[0]):
        treeSumRemain = np.copy(treeSum)
        nodesArray = np.array(nodes).astype(int)
        treeSumRemain[nodesArray] = -1

        argMax = int(np.argmax(treeSumRemain))
        nodes.append(argMax)
        treeSumIn = np.copy(treeSum) + ((1-treeData[:, argMax]) * 1000000)
        argMin = int(np.argmin(treeSumIn))

        if treeData[argMin, argMax] == 1:
            edges[argMin, argMax] = 1
    return edges

def findFatherNode():
    edges = findTreeEdges()
    fatherNode = np.argmax(edges, axis=0)
    return fatherNode

def findBrotherNode():
    edges = findTreeEdges()
    fatherNode = findFatherNode()
    kidsOfFather = edges[fatherNode]
    kidsOfFather[np.arange(edges.shape[0]), np.arange(edges.shape[0])] = 0
    brotherNode = np.argmax(kidsOfFather, axis=1)
    return brotherNode

def interweive(data):
    M = data.shape[0]
    arange1 = np.arange(data.shape[1])
    arange2 = arange1 // M
    arange3 = arange1 % M
    data4 = data[arange3, arange2]
    return data4

def realign_hmm_seq():
    ''' This function takes individual sequences from sequenceFileNames
    and induces them on the root sequence which is the fasta file in the
    0th hmm folder(A_0_0). This gets saved to n.fasta where n is the index in the
    sequneceFileNames, which is the same as the suffix number in the folder structure A_0_n.
    '''
    sequenceFileNames = giveSequenceFileNames()
    rootName = sequenceFileNames[0]
    rootKeys, rootSeqs = loadFastaBasic(rootName)

    for a in range(0, len(sequenceFileNames)):
        name = sequenceFileNames[a]
        nowKeys, nowSeqs = loadFastaBasic(name)
        argsIn = np.argwhere(np.isin(rootKeys, nowKeys))[:, 0]
        nowKeys, nowSeqs = rootKeys[argsIn], rootSeqs[argsIn]
        dataFolderName = giveAllFileNames()[4]
        ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSeqAlign/" + str(a) + '.fasta')
        saveFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSeqAlign/" + str(a) + '.fasta', nowKeys, nowSeqs)

def compareHMM():
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        HMM1 = "%s/initialHMM/testStandardAlign" % dirName + str(a) + ".hmm"
        HMM2 = "%s/initialHMM/test" % dirName + str(a) + ".hmm"
        cmd = 'diff ' + HMM1 + ' ' + HMM2
        print (cmd)

def save_adjusted_score():
    ''' This function takes every single hmmsearch bitscore result thaht was stored in full.npy
    and multiplies them by the log base 2 of the size of the sequence set the hmm was built on.
    These adjusted scores get saved to fullAdjusted.npy
    '''
    dataFolderName = giveAllFileNames()[4]
    sets = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)
    sizes1 = []
    for a in range(len(sets)):
        sizes1.append(len(sets[a]))
    sizes1 = np.array(sizes1)
    scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    heightShift = np.log2(sizes1) * 1
    for a in range(len(scores[0])):
        scores[:, a] = scores[:, a] + heightShift[a]
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy", scores)

def save_initial_steps(abstract_algorithm):
    deduplication_start = default_timer()
    save_sequence_filenames()
    deduplication_end = default_timer()
    deduplication_time = deduplication_end - deduplication_start

    build_hmm_start = default_timer()
    save_initial_hmm(abstract_algorithm)
    build_hmm_end = default_timer()
    build_hmm_time = build_hmm_end - build_hmm_start

    build_array_start = default_timer()
    save_hmm_sets()
    build_array_end = default_timer()
    build_array_time = build_array_end - build_array_start

    build_ancestor_matrix_start = default_timer()
    save_decomposition()
    build_ancestor_matrix_end = default_timer()
    build_ancestor_matrix_time = build_ancestor_matrix_end - build_ancestor_matrix_start

    induce_start = default_timer()
    realign_hmm_seq()
    induce_end = default_timer()
    induce_time = induce_end - induce_start
    # save_scores_original(abstract_algorithm)
    # save_adjusted_score()

    time_output_suffix = "time.data"
    with open(abstract_algorithm.get_upp2_output_filename(time_output_suffix), "a+") as f:
        f.write(f"2-a-deduplication: {deduplication_time} seconds ({timedelta(seconds=deduplication_time)})\n")
        f.write(f"2-b-build_hmm: {build_hmm_time} seconds ({timedelta(seconds=build_hmm_time)})\n")
        f.write(f"2-c-build_array: {build_array_time} seconds ({timedelta(seconds=build_array_time)})\n")
        f.write(f"2-d-build_ancestor_matrix: {build_ancestor_matrix_time} seconds ({timedelta(seconds=build_ancestor_matrix_time)})\n")
        f.write(f"2-e-induce: {induce_time} seconds ({timedelta(seconds=induce_time)})\n")

def doPoisonRemoval(maxHMM, scores, queryNum):
    treeData = findDecomposition()
    treeSum = np.sum(treeData, axis=1)
    brotherNode = findBrotherNode()
    fatherNode = findFatherNode()
    treeDataMinusSame = np.copy(treeData)
    treeDataMinusSame[np.arange(treeData.shape[0]), np.arange(treeData.shape[0])] = 0
    nodeToCheck = treeDataMinusSame[maxHMM]
    argsCheck = np.argwhere(nodeToCheck == 1)

    queryCheck = queryNum[argsCheck[:, 0]]
    nodeCheck = argsCheck[:, 1]

    nodeFather = fatherNode[nodeCheck]
    nodeBrother = brotherNode[nodeCheck]

    scoreNode = scores[queryCheck, nodeCheck]
    scoreFather = scores[queryCheck, nodeFather]
    scoreBrother = scores[queryCheck, nodeBrother]

    ep1 = 1e-3
    diff1 = ((scoreFather + ep1) / (scoreNode + ep1)) - 1.1
    diff2 = ((scoreBrother + ep1) / (scoreFather + ep1)) - 1.1

    argsPoison = np.argwhere(np.logical_and(  diff1 > 0, diff2 > 0 ))[:, 0]

    queryPoison = queryCheck[argsPoison]
    nodePoison = nodeCheck[argsPoison]

    keepMask = treeData[maxHMM]
    keepMask[queryPoison, nodePoison] = 0

    return keepMask

def saveFromBooleanInclude(hmmLeafs, strategyName):
    dataFolderName = giveAllFileNames()[4]
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(hmmLeafs)):
        newFasta = []
        for b in range(0, len(hmmLeafs[0])):
            if hmmLeafs[a, b] == 1:
                name = get_root_temp_dir() + '/alignData/hmmSeqAlign/' + str(b) + '.fasta'
                data = loadFastaFormat(name)
                keysUsed, _ = loadFastaBasic(name)
                newFasta = newFasta + data
        saveFasta(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)

def scoresToHMMSeq(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/")
    if strategyName in ["adjusted_bitscore", "upp"]: # ['stefan_UPP', 'stefan_UPPadjusted']:
        save_scores_original(abstract_algorithm)
        save_adjusted_score()
        scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")
        if strategyName == "upp":
            scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy")
        scores = scores / np.max(scores)
        treeData = findDecomposition()
        treeSum = np.sum(treeData, axis=1)
        maxHMM = np.argmax(scores, axis=1).astype(int)
        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)
        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T
        _LOG.debug(f"saving to {get_root_temp_dir()}/data/internalData/{dataFolderName}/{strategyName}/queryToHmm/original.npy")
        np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)
        for a in range(0, len(HMMunique)):
            HMMnum = HMMunique[a]
            # name = get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(a) + '.fasta'
            name = get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            keys, seqs = loadFastaBasic(name)
            saveFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)
    elif strategyName in ["hierarchical"]: # ['stefan_fastUPP', 'stefan_fastUPPexception', 'stefan_fastUPPearly']:
        if strategyName == 'stefan_fastUPPearly':
            scoresFull = np.load('%s/ensembleData/Searcher/scoreFiles/Early_score.npy' % dirName)
        else:
            # normally score.npy stores the raw bitscores while fullAdjusted.npy stores the adjusted bitscores
            # however if the
            scoresFull = np.load('%s/ensembleData/Searcher/scoreFiles/score.npy' % dirName)
        maxHMM = np.argmax(scoresFull, axis=1).astype(int)
        if strategyName == 'stefan_fastUPPexception':
            scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")
            scores_original = np.load("%s/hmmScores/full.npy" % dirName)
            treeData = findDecomposition()
            treeSum = np.sum(treeData, axis=1)
            treeTop = np.unique(treeSum)[-3:]
            maxHMM2 = np.argmax(scores_original, axis=1).astype(int)
            argsReplace = np.argwhere(np.isin(treeSum[maxHMM], treeTop))[:, 0]
            maxHMM[argsReplace] = maxHMM2[argsReplace]
        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)
        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T
        np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)
        for a in range(0, len(HMMunique)):
            HMMnum = HMMunique[a]
            # this is the induced version of each backbone subalignment on the full backbone alignment
            name = get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            # this is probably space and time inefficient. We could most likely just move or create a symlink
            keys, seqs = loadFastaBasic(name)
            saveFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)

def generateNewHMM(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    numHMM = int(np.max(np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")[:, 1])+1)
    for a in range(0, numHMM):
        src = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta"
        keys, seqs = loadFastaBasic(src)
        seqMatrix = []
        for seq in seqs:
            seqMatrix.append(list(seq))
        seqMatrix = np.array(seqMatrix)
        seqMatrix2 = seqMatrix.T
        usedCols = []
        for b in range(seqMatrix2.shape[0]):
            list1 = seqMatrix2[b]
            sizeReal = list1[list1 != '-'].shape[0]
            if sizeReal != 0:
                usedCols.append(b)
        usedCols = np.array(usedCols)
        usedCols_filename = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/columnSets/" + str(a) + ".npy"
        ensureFolder(usedCols_filename)
        np.save(usedCols_filename, usedCols)
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(a) + ".hmm"
        ensureFolder(hmmName)
        addHMMBuildJob(abstract_algorithm, hmmName, src)
    _LOG.debug(f"JobPool actual cpus: {JobPool().cpus}")
    JobPool().wait_for_all_jobs()


def stockholmToFasta(stockholmName, fastaName):
    fakeFastaName = '%s/temporaryFileSave/fakeFasta_1.fasta' % dirName
    records = SeqIO.parse(stockholmName, "stockholm")
    count = SeqIO.write(records, fakeFastaName, "fasta")
    dataAlign = loadFastaFormat(fakeFastaName)

    newData = ''.join(dataAlign)
    newData = newData.replace('\n', '')
    newData = newData.split('>')[1:]

    keys = []
    seqs = []
    for a in range(0, len(newData)):
        keys.append(newData[a][:10])
        seqs.append(newData[a][10:])

    saveFastaBasic(fastaName, keys, seqs)

def loadStockholm(stockholmName):
    data = []
    started = False
    second = False
    startedLine = 0
    secondLine = 10000000000
    file_obj = open(stockholmName)
    count1 = 0
    for line_number, i in enumerate(file_obj):
        if started:
            position = (line_number - startedLine - 1) % (secondLine - startedLine)
            if position >= len(data):
                data.append([])
            data[position].append(i)

        if (i == '\n') and ((started == True) and (second == False)):
            second = True
            secondLine = line_number

        if (i == '\n') and (started == False):
            started = True
            startedLine = line_number

    data2 = []
    keys = []

    for a in range(0, len(data) - 1):
        str1 = ''
        key =  data[a][0]
        keyList = np.array(list(key))
        argFirst = np.argwhere(keyList == ' ')[-1, 0] + 1
        key = key[:argFirst]
        keys.append(key)

        for b in range(0, len(data[a])):
            str2 = data[a][b][argFirst:-1]
            str1 += str2
        data2.append(str1)

    data2 = np.array(data2)
    keys = np.array(keys)
    return keys, data2

def loadStockholmOnlySeqs(stockholmName):
    keys, seqs = loadStockholm(stockholmName)
    seqNum = (len(seqs) - 2) // 2
    seqsNew = []
    keysNew = []
    for a in range(0, seqNum):
        seqsNew.append(seqs[a * 2])
        key = keys[a*2]
        keyList = np.array(list(key))
        argFirst = np.argwhere(keyList == ' ')[0, 0]
        key = key[:argFirst]
        keysNew.append(key)

    seqsNew = np.array(seqsNew)
    keysNew = np.array(keysNew)

    return keysNew, seqsNew

def txtToFasta(txtName, fastaName):
    file_obj = open(txtName)
    dataAlign = []
    count1 = 0
    for line_number, i in enumerate(file_obj):
        dataAlign.append(i)
    newData = ''.join(dataAlign)
    newDataList = np.array(list(newData))
    allbreaks = np.argwhere(newDataList == '\n')[0::2, 0]
    allCarrot = np.argwhere(newDataList == '>')[:, 0]
    keyLengths = (allbreaks - allCarrot - 1).astype(int)

    newData = newData.replace('\n', '')
    newData = newData.split('>')[1:]

    keys = []
    seqs = []
    for a in range(0, len(newData)):
        keys.append(newData[a][:keyLengths[a]])
        seqs.append(newData[a][keyLengths[a]:])
    keys, seqs = np.array(keys), np.array(seqs)
    saveFastaBasic(fastaName, keys, seqs)

def alignQueries(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    queryName = giveQueryFileName()
    queryData = loadFastaFormat(queryName)
    queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")[:, 1]
    #queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy")
    #uppHMM = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/UPPused.npy")

    predictionDataFull = []
    for b in range(0, int(np.max(queryToHmm)) + 1):
        argsHMM = np.argwhere(queryToHmm == b)[:, 0]
        if argsHMM.shape[0] > 0:
            dataNow = []
            keysNow = []
            for a in range(0, len(queryData)):
                if queryToHmm[a] == b:
                    dataNow.append(queryData[a])
            assert len(dataNow) == argsHMM.shape[0]
            queryName1 = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/inputQuery/' + str(b) + '.fasta'
            ensureFolder(queryName1)
            saveFasta(queryName1, dataNow)
            predictionName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/predictedQuery/' + str(b) + '.sto'
            ensureFolder(predictionName)
            #if uppHMM[argsHMM[0]] == -1:
            hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(b) + ".hmm"
            #else:
            #    uppNum = int(uppHMM[argsHMM[0]])
            #    hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/initialHMM/test" + str(uppNum) + ".hmm"
            ensureFolder(hmmName)
            addHMMAlignJob(abstract_algorithm, hmmName, queryName1, predictionName)
    JobPool().wait_for_all_jobs()

def mergeAlignments(abstract_algorithm, strategyName, overlapLowercase=True):
    dataFolderName = giveAllFileNames()[4]
    #queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy")
    queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")[:, 1]
    predictionDataFull = []
    backboneKeys, backboneSeqs = loadFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(0) + '.fasta')
    backBoneChoice = np.zeros((queryToHmm.shape[0], len(backboneSeqs[0]))).astype(str)
    backBoneChoice[:] = '-'
    backBoneChoiceBool = np.zeros(backBoneChoice.shape)

    fullInsertions = []
    for a in range(len(backboneSeqs[0])+1):
        fullInsertions.append([])
    fullInsertionsIndex = []
    for a in range(len(backboneSeqs[0])+1):
        fullInsertionsIndex.append([])

    fullInsertionsNumber = np.zeros(len(backboneSeqs[0])+1).astype(int)
    queryNames = np.zeros(queryToHmm.shape[0]).astype(str)
    
    for a in range(0, int(np.max(queryToHmm)) + 1):
        argsHMM = np.argwhere(queryToHmm == a)[:, 0]
        if argsHMM.shape[0] > 0:
            #This loads in the insertions from the alignment stockholm, and converts it to a binary where 1 is an insertion and 0 is a match state. 
            predictionName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/predictedQuery/'  + str(a) + '.sto'
            _, insertions = loadStockholm(predictionName)
            insertions = insertions[-1]
            insertions = np.array(list(insertions))
            insertions[insertions == '.'] = 1
            insertions[insertions == 'x'] = 0
            insertions = insertions.astype(int)
        
            #This loads in the predicted alignments and converts it to a 2 dimensional array
            keys, seqs = loadStockholmOnlySeqs(predictionName)
            queryNames[argsHMM] = keys
            seqsArray = []
            for seq in seqs:
                seqsArray.append(list(seq))
            seqsArray = np.array(seqsArray)
            
            #This makes sure their is consistancy between the input queries and the prediction stockholm file. 
            queryFileName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/inputQuery/'  + str(a) + '.fasta'
            queryKey, querySeq = loadFastaBasic(queryFileName)
            assert len(keys) == len(queryKey)
            for b in range(len(queryKey)):
                assert queryKey[b] == keys[b]
            
            #This loads in the columns of the backbone which are used by this hmm. 
            usedCols = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/newHMM/columnSets/' + str(a) + ".npy")
            #This indicates where the alignment has match states rather than insertions.
            matchPosition = np.argwhere(insertions == 0)[:, 0]

            try:
                assert len(usedCols) == len(matchPosition)
            except:
                print ("This issue means that the number of match states from the HMM is")
                print ("not equal to the number of non-empty columns on the sequences")
                print ("the hmm was trained on.")
                assert len(usedCols) == len(matchPosition)
            
            #This sets the value of the alignment in the backbone positions. 
            backBoneChoice[np.ix_(argsHMM, usedCols)] = seqsArray[:, matchPosition]
            backBoneChoiceBool[np.ix_(argsHMM, usedCols)] = 1
            
            #This is used to determine insrtNumber which is the number of insertions between the match states.
            #Specifically, insertionNumber[i] is the number of insertions between the ith match state and i+1 match state.
            insertionsNew = np.concatenate((np.zeros(1), insertions))
            insertionsNew = np.concatenate((insertionsNew, np.zeros(1)))
            insertionsNew = insertionsNew.astype(int)

            insertions_cumsum = np.cumsum(insertionsNew)
            matchPosition = np.argwhere(insertionsNew == 0)[:, 0]
            insertNumber = insertions_cumsum[matchPosition[1:]] - insertions_cumsum[matchPosition[:-1]]
            usedColsExtra = np.concatenate((np.zeros(1), usedCols+1)).astype(int)

           
            usedColsExtra[-1] = len(backboneSeqs[0])

            for b in range(len(insertNumber)):
                if insertNumber[b] > 0:
                    #If there are insertions to make (insertNumber[b] > 0) the following code saves the insertions.
                    toInsert = seqsArray[:, matchPosition[b]:matchPosition[b]+insertNumber[b]]
                    insertPosition = usedColsExtra[b]
                    fullInsertions[insertPosition].append(np.copy(toInsert)) #Saves the values to insert
                    fullInsertionsIndex[insertPosition].append(np.copy(argsHMM)) #Saves the indices of the queries where the insertions are made.
                    if overlapLowercase:
                        #If a compressed format is used, insertions overlap with each other, so we only need as many insertions as the
                        #maximum number of insertions which occurs between these match states.
                        fullInsertionsNumber[insertPosition] = max(insertNumber[b], fullInsertionsNumber[insertPosition])
                    else:
                        #If a compressed format is not used, each new insertion requires a new column.
                        fullInsertionsNumber[insertPosition] += insertNumber[b]
    
    #This will be the new alignment with all insertions included. 
    #newAlignment = np.zeros((queryToHmm.shape[0], int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)
    newAlignment = np.zeros((queryToHmm.shape[0] + len(backboneKeys), int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)

    newAlignment[:] = '-'
    newAlignmentBool = np.zeros(newAlignment.shape)
    #This gives the index of the match states in the new alignment. 
    #If there are N insertions between match i and i+1, then the index of match state i+1 must be N+1 larger than
    #the index of match state i.
    colIndex = np.cumsum(fullInsertionsNumber+1).astype(int)
    colIndex = np.concatenate((np.zeros(1), colIndex[:-1])).astype(int)
    #This puts the back bone into the new alignment
    newAlignment[:queryToHmm.shape[0], colIndex[1:]] = backBoneChoice 
    newAlignmentBool[:queryToHmm.shape[0], colIndex[1:]] = backBoneChoiceBool * 2
    
    #This whole loop adds the insertions into the new alignment
    for a in range(0, len(fullInsertions)):
        start1 = colIndex[a] + 1
        start2 = start1
        for b in range(0, len(fullInsertions[a])):
            #This says the values to insert
            toInsert = fullInsertions[a][b]
            #This says the queries where the insertions are being made
            insertQuery = fullInsertionsIndex[a][b]
            #This says how many columns will be inserted
            size1 = toInsert.shape[1]
            #This says which columns to insert the insertions
            insertPosition = np.arange(size1) + start2 
            #This places the insertion into the new alignment
            newAlignment[np.ix_(insertQuery, insertPosition)] = np.copy(toInsert)
            newAlignmentBool[np.ix_(insertQuery, insertPosition)] = 1

            if not overlapLowercase:
                start2 += size1
                
    #Skip states are represented by '-' in our output format.
    newAlignment[newAlignment == '.'] = '-'

    #This adds the backbone sequences into the new alignment. 
    queryNames = list(queryNames)
    for a in range(len(backboneKeys)):
        #backboneKeys, backboneSeqs
        seq1 = backboneSeqs[a]
        key1 = backboneKeys[a]
        
        #This adds the backbone sequence to the backbone column (colIndex). 
        newAlignment[a + queryToHmm.shape[0], colIndex[1:]] = np.copy(np.array(list(seq1)))
        queryNames.append(key1)
    queryNames = np.array(queryNames)
    
    #This removes an immaginary match state before all match states, which is a useful conceptual
    #fiction which allows all insertions to be after some match state. 
    newAlignment = newAlignment[:, 1:]
    newAlignmentBool = newAlignmentBool[:, 1:]
    colIndexTrue = (colIndex[1:] - 1).astype(int)
    
    #This converts the 2 dimensional matrix alignment to a list of sequences. 
    newAlignmentString = []
    for a in range(newAlignment.shape[0]):
        str1 = ''.join(list(newAlignment[a]))
        newAlignmentString.append(copy.copy(str1))
    newAlignmentString = np.array(newAlignmentString)
    
    #This saves the new alginemtn
    fastaName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentFasta.fasta'
    ensureFolder(fastaName)
    saveFastaBasic(abstract_algorithm.get_output_filename("alignment.fasta"), queryNames, newAlignmentString)
    saveFastaBasic(fastaName, queryNames, newAlignmentString)
    
    #This saves only the query sequence part of the new alignment, which is useful for scoring the results. 
    fastaNameQuery = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/query_alignmentFasta.fasta'
    saveFastaBasic(fastaNameQuery, queryNames[:queryToHmm.shape[0]], newAlignmentString[:queryToHmm.shape[0]])
    
    #This saves extra information in addition the alignment fasta. 
    #This is the new alignment in a 2 dimensional array
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignment.npy', newAlignment)
    #This is the elements of the new alignment which correspond to the old subalignments
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentBool.npy', newAlignmentBool)
    #This is the columns in the new alignment which correspond to the backbone. 
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/columnIndex.npy', colIndexTrue)

def InputMergeAlignments(queryNamesFull, alignments, allInsertions, columnSets, Ncolumns, overlapLowercase=True):
    Nquery = 0
    for alignment in alignments:
        Nquery += len(alignment)

    backBoneChoice = np.zeros((Nquery, Ncolumns)).astype(str)
    backBoneChoice[:] = '-'
    backBoneChoiceBool = np.zeros(backBoneChoice.shape)

    fullInsertions = []
    for a in range(Ncolumns+1):
        fullInsertions.append([])

    fullInsertionsIndex = []
    for a in range(Ncolumns+1):
        fullInsertionsIndex.append([])

    fullInsertionsNumber = np.zeros(Ncolumns+1).astype(int)
    queryNames = np.zeros(Nquery).astype(str)
    queryCount = 0
    for a in range(0, len(alignments)):
        sizeAlign = len(alignments[a])
        argsHMM = np.arange(sizeAlign) + queryCount
        queryCount += sizeAlign
        if argsHMM.shape[0] > 0:
            insertions = allInsertions[a]
            insertions = np.array(list(insertions))
            insertions[insertions == '.'] = 1
            insertions[insertions == 'x'] = 0
            insertions = insertions.astype(int)
            seqs = alignments[a]
            queryNames[argsHMM] = queryNamesFull[argsHMM]
            seqsArray = []
            for seq in seqs:
                seqsArray.append(list(seq))
            seqsArray = np.array(seqsArray)
            usedCols = columnSets[a]
            matchPosition = np.argwhere(insertions == 0)[:, 0]
            backBoneChoice[np.ix_(argsHMM, usedCols)] = seqsArray[:, matchPosition]
            backBoneChoiceBool[np.ix_(argsHMM, usedCols)] = 1

            insertionsNew = np.concatenate((np.zeros(1), insertions))
            insertionsNew = np.concatenate((insertionsNew, np.zeros(1)))
            insertionsNew = insertionsNew.astype(int)

            insertions_cumsum = np.cumsum(insertionsNew)
            matchPosition = np.argwhere(insertionsNew == 0)[:, 0]
            insertNumber = insertions_cumsum[matchPosition[1:]] - insertions_cumsum[matchPosition[:-1]]

            usedColsExtra = np.concatenate((np.zeros(1), usedCols+1)).astype(int)

            if True:
                usedColsExtra[-1] = Ncolumns

            for b in range(len(insertNumber)):
                if insertNumber[b] > 0:
                    toInsert = seqsArray[:, matchPosition[b]:matchPosition[b]+insertNumber[b]]
                    insertPosition = usedColsExtra[b]
                    fullInsertions[insertPosition].append(np.copy(toInsert))
                    fullInsertionsIndex[insertPosition].append(np.copy(argsHMM))
                    if overlapLowercase:
                        fullInsertionsNumber[insertPosition] = max(insertNumber[b], fullInsertionsNumber[insertPosition])
                    else:
                        fullInsertionsNumber[insertPosition] += insertNumber[b]

    newAlignment = np.zeros((Nquery, int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)
    newAlignment[:] = '-'
    newAlignmentBool = np.zeros(newAlignment.shape)
    colIndex = np.cumsum(fullInsertionsNumber+1).astype(int)
    colIndex = np.concatenate((np.zeros(1), colIndex[:-1])).astype(int)

    newAlignment[:, colIndex[1:]] = backBoneChoice #[1:] removes fake "before everything" match state.
    newAlignmentBool[:, colIndex[1:]] = backBoneChoiceBool * 2

    for a in range(0, len(fullInsertions)):
        start1 = colIndex[a] + 1
        start2 = start1
        for b in range(0, len(fullInsertions[a])):
            toInsert = fullInsertions[a][b]
            insertQuery = fullInsertionsIndex[a][b]
            size1 = toInsert.shape[1]
            insertPosition = np.arange(size1) + start2
            newAlignment[np.ix_(insertQuery, insertPosition)] = np.copy(toInsert)
            newAlignmentBool[np.ix_(insertQuery, insertPosition)] = 1
            if not overlapLowercase:
                start2 += size1

    newAlignment[newAlignment == '.'] = '-'
    newAlignment = newAlignment[:, 1:]
    newAlignmentBool = newAlignmentBool[:, 1:]
    colIndexTrue = (colIndex[1:] - 1).astype(int)
    newAlignmentString = []
    for a in range(newAlignment.shape[0]):
        str1 = ''.join(list(newAlignment[a]))
        newAlignmentString.append(copy.copy(str1))
    newAlignmentString = np.array(newAlignmentString)
    return newAlignmentString

def easyInputMerge(alignments, columnSets, overlapLowercase=True):
    Ncolumns = 0
    for set1 in columnSets:
        Ncolumns = max(Ncolumns, np.max(set1))
    Ncolumns = int(Ncolumns) + 1
    Nquery = 0
    for align in alignments:
        for seq in align:
            Nquery += 1
    lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))
    queryNamesFull = np.arange(Nquery).astype(str)
    allInsertions = []
    for a in range(len(alignments)):
        seqArray = []
        for seq in alignments[a]:
            seqArray.append(list(seq))
        seqArray = np.array(seqArray).T
        insertions1 = ''
        for b in range(len(seqArray)):
            seq1 = seqArray[b]
            if np.intersect1d(seq1, lowercase).shape[0] == 0:
                insertions1 = insertions1 + 'x'
            else:
                insertions1 = insertions1 + '.'
        allInsertions.append(insertions1)
    seqsNew = InputMergeAlignments(queryNamesFull, alignments, allInsertions, columnSets, Ncolumns, overlapLowercase=overlapLowercase)
    return seqsNew

def checkAlignmentsValid(strategyName):
    newAlignment = np.load('%s/hmmQueryList/merged/' % dirName + strategyName + '_alignment.npy')
    newAlignmentBool = np.load('%s/hmmQueryList/merged/' % dirName + strategyName + '_alignmentBool.npy')
    queryToHmm = np.load("%s/queryToHmm/" % dirName + strategyName + ".npy")
    for a in range(0, int(np.max(queryToHmm)) + 1):
        argsHMM = np.argwhere(queryToHmm == a)[:, 0]
        if argsHMM.shape[0] > 0:
            predictionName = get_root_temp_dir() + 'alignData/hmmQueryList/predictedQuery/' + strategyName + '_'  + str(a) + '.sto'
            keys, seqs = loadStockholmOnlySeqs(predictionName)
            seqsArray = []
            for b in range(0, len(seqs)):
                seq = seqs[b]
                indexInNew = argsHMM[b]
                seq2 = np.array(list(seq))
                newSeq = newAlignment[indexInNew][newAlignmentBool[indexInNew] != 0]
                newSeqBool = newAlignmentBool[indexInNew][newAlignmentBool[indexInNew] != 0]
                for c in range(0, len(seq2)):
                    if seq2[c] != newSeq[c]:
                        print ("Issue")
                        quit()

def buildAlignMerge(abstract_algorithm, strategyName, doResort=False):
    generateNewHMM(abstract_algorithm, strategyName)

    saveNewScores(abstract_algorithm, strategyName)

    alignQueries(abstract_algorithm, strategyName)
    mergeAlignments(abstract_algorithm, strategyName)

def saveTrueUPPSubset():
    fileNames = giveAllFileNames()
    dataFolderName = giveAllFileNames()[4]
    predictionName = fileNames[3]

    predictionNameStefan = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/' + 'stefan_UPP' '/hmmQueryList/merged/alignmentFasta.fasta'
    predictionNameNew = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'

    ensureFolder(predictionNameStefan)
    ensureFolder(predictionNameNew)
    trueKey, trueSeq = loadFastaBasic(predictionName)
    predKey, predSeq = loadFastaBasic(predictionNameStefan)

    trueKey, trueSeq = trueKey[np.isin(trueKey, predKey)], trueSeq[np.isin(trueKey, predKey)]
    saveFastaBasic(predictionNameNew, trueKey, trueSeq)

def saveTrueAlignFasta():
    fileNames = giveAllFileNames()
    trueAlignment = fileNames[2]
    dataFolderName = giveAllFileNames()[4]
    trueFasta = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/trueAlignment/true_align_fragged_fasta.fasta'
    ensureFolder(trueFasta)
    txtToFasta(trueAlignment, trueFasta)

def doAllSteps(strategyName):
    buildAlignMerge(strategyName)

def scoreAlignment(strategyName):
    dataFolderName = giveAllFileNames()[4]
    saveTrueUPPSubset()
    saveTrueAlignFasta()
    trueFasta = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/trueAlignment/true_align_fragged_fasta.fasta'
    predictionName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/query_alignmentFasta.fasta'
    ensureFolder(predictionName)

    UPPName = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'
    ensureFolder(UPPName)
    trueFastaNew = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/trueAlignment/subset/true_fasta.fasta'
    ensureFolder(trueFastaNew)

    queryName = giveQueryFileName()
    queryKey, _ = loadFastaBasic(queryName)
    trueKey, trueSeq = loadFastaBasic(trueFasta)
    predKey, predSeq = loadFastaBasic(predictionName)
    uppKey, uppSeq = loadFastaBasic(UPPName)

    argsIn1 = []
    for a in range(len(predKey)):
        argUse1 = np.argwhere(trueKey == predKey[a])[0, 0]
        argsIn1.append(argUse1)
    argsIn1 = np.array(argsIn1)
    trueKey, trueSeq = trueKey[argsIn1], trueSeq[argsIn1]

    saveFastaBasic(trueFastaNew, trueKey, trueSeq)

    cmd = "java -jar ./FastSP.jar -ml -r " + trueFastaNew + " -e " + predictionName

    scoreInfo  = os.popen(cmd).readlines()

    scoreInfo2 = []
    for a in range(0, len(scoreInfo)):
        string1 = scoreInfo[a].split(' ')[1].replace('\n', '')
        scoreInfo2.append(float(string1))

    print (strategyName)
    print (scoreInfo)

def compareToUPP():
    keys1, seqs1 = loadFastaBasic(predictionStefan)
    keys2, seqs2 = loadFastaBasic(predictionUPP)

    keys2, seqs2 = keys2[np.isin(keys2, keys1)], seqs2[np.isin(keys2, keys1)]

    print (keys1[0])
    print (keys2[0])
    print (seqs1[0])
    print ('')
    print (seqs2[0])


def hierchySearch(abstract_algorithm, adjusted_bitscore, early_stop, fakeSimulate=False):
    observeNephew = False

    dataFolderName = giveAllFileNames()[4]
    # np.load("./data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    # scores_original is the raw bitscores without the adjusted weights
    if fakeSimulate:
        scores_original = np.load(get_root_temp_dir() + "/data/internalData/%s/hmmScores/full.npy" % (dataFolderName))
    # so here the size of each hmm is stored in hmm_sizes where hmm_sizes[i] is the size of the ith hmm
    hmm_sets = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)
    hmm_sizes = []
    for current_hmm_set in hmm_sets:
        hmm_sizes.append(len(current_hmm_set))
    hmm_sizes = np.array(hmm_sizes)

    hmmNames = []
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        inputFiles = giveAllFileNames()
        dataFileName = inputFiles[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFileName + "/initialHMM/test" + str(a) + ".hmm"
        hmmNames.append(hmmName)

    queryName = giveQueryFileName()
    queryData = loadFastaFormat(queryName)
    Nquery = int(len(queryData))
    Nhmm = int(len(hmmNames))

    # this returns a 2d numpy array where array[a,b] = 1 when a has an edge to b
    edges = findTreeEdges()
    # since edges is only the immediate children, argsorting and getting the last two
    # is then the two children of the node so this is a 2d array where
    # every row represents an HMM and the two array elements are the children HMMs,
    # or the indices/number representation/order inside edges/whatever of those children HMMs
    children = np.argsort(edges, axis=1)[:, -2:]
    # this is the number of children that each HMM has so it's a 1d array where
    # each childNum[2] returns the number of direct children for the 2nd HMM
    childNum = np.sum(edges, axis=1)
    brotherNode = findBrotherNode()

    # this returns a 2d numpy array where array[a,b] = 1 when a is a parent(not just the direct parent) of b
    treeData = findDecomposition()
    # this has the number of descendents for for each HMM in a 1d array
    treeSum = np.sum(treeData, axis=1)
    sortedSize = np.argsort(treeSum)
    scoresFull = np.zeros((Nquery, Nhmm)) - 1000
    scoreAdjustment = np.zeros(scoresFull.shape)
    initialSearchSet = np.array([0])

    queryPart = np.arange(Nquery * len(initialSearchSet)) % Nquery
    hmmPart = np.array(initialSearchSet).repeat(Nquery)
    # here we run hmmsearch on every singly query sequence againstn the 0 hmm, or the root hmm
    queryHMM = np.array([queryPart, hmmPart]).T
    scoreName = ''
    if fakeSimulate:
        scores = np.copy(scores_original)
    else:
        # I guess this method is called saveScoreFromBool since it runs hmmsearch just like
        # saveScore but only runs query sequences according to a 2d table where
        # the rows are (query sequence, hmm) pair and we only run those pairs using boolean logic?
        scores = saveScoreFromBool(abstract_algorithm, hmmNames, queryHMM, scoreName, noSave=True)
    scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = scores[queryHMM[:, 0], queryHMM[:, 1]]
    if(adjusted_bitscore):
        scoresFull[queryHMM[:,0], queryHMM[:,1]] = scoresFull[queryHMM[:,0], queryHMM[:,1]] + np.log2(hmm_sizes[queryHMM[:,1]])

    queryHMM_original = np.copy(queryHMM)
    # maxHMM is initally the best scoring hmm which would always be zero since queryHMM[:1] is just zeros initially
    # if we do an argmax along axis 1 then we are looking at the indices for each row the maximum scoring HMM, or
    # in this case the only hmm, hmm 0.
    # so this argmax can be skipped entirely, especially since bool_done is hardcoded to zero anyways
    maxHMM = np.argmax(scoresFull, axis=1)
    bool_done = np.zeros((Nquery, Nhmm))
    # for every query sequence, we mark the 0th root hmm as done, or 1,
    # since we already ran them
    bool_done[:, 0] = 1
    maxHMMlist = []
    scoresRecord = []

    done = False
    num_iteration = 0
    while not done:
        # so at every iteration, maxNotLeaf is going to be updated with the
        # children of maxHMM that have at least one child
        # so this is called maxNotLeaf since it's the children of maxHMM that are not leaves
        # friendly reminder that maxHMM is an array where the ith element is the ith queries best HMM
        maxNotLeaf = np.argwhere(childNum[maxHMM] != 0)[:, 0]

        if observeNephew:
            maxHMMChild_original = children[maxHMM[maxNotLeaf]]
            maxHMMChild_original = np.concatenate((maxHMMChild_original, children[brotherNode[maxHMM[maxNotLeaf]]]), axis=1)
        else:
            maxHMMChild_original = children[maxHMM[maxNotLeaf]]
        # we setup the queryHMM to be all the immediate children of the maxHMM that are not leaves
        # with all the query sequences so we do the arange with division so that for the two children
        # we can give one query sequence to both of them and the array just gets flattened like 00 11 22 etc
        maxHMMChild = np.copy(maxHMMChild_original).reshape((maxHMMChild_original.size,))
        queryCorrespond = maxNotLeaf[np.arange(maxHMMChild.size) // (maxHMMChild_original.shape[1]) ]
        queryToHMM_try = np.array([queryCorrespond, maxHMMChild]).T
        queryHMM = queryToHMM_try[bool_done[queryToHMM_try[:, 0], queryToHMM_try[:, 1] ] == 0]
        # and now we can mark these HMMs done so that we don't pick these again via a mask next time
        bool_done[queryHMM[:, 0], queryHMM[:, 1] ] = 1

        if queryHMM.shape[0] != 0:
            num_iteration += 1
            scoreName = ''
            if fakeSimulate:
                scores = np.copy(scores_original)
            else:
                scores = saveScoreFromBool(abstract_algorithm, hmmNames, queryHMM, scoreName, noSave=True)

            scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = np.copy(scores[queryHMM[:, 0], queryHMM[:, 1]])
            newScores1 = scoresFull[maxNotLeaf, maxHMMChild_original[:, 0]]
            newScores2 = scoresFull[maxNotLeaf, maxHMMChild_original[:, 1]]
            if(adjusted_bitscore):
                # I'm getting the sizes of the hmms and then adjusting the bitscores here
                score1_size_weight = np.log2(hmm_sizes[maxHMMChild_original[:,0]])
                score2_size_weight = np.log2(hmm_sizes[maxHMMChild_original[:,1]])
                # print(queryHMM[:,1])
                # print(hmm_sizes[queryHMM[:,1]])
                # print(queryHMM[:,1].size)
                # print(hmm_sizes[queryHMM[:,1]].size)
                # print(np.log2(hmm_sizes[queryHMM[:,1]])[0])
                scoresFull[queryHMM[:,0], queryHMM[:,1]] = scoresFull[queryHMM[:,0], queryHMM[:,1]] + np.log2(hmm_sizes[queryHMM[:,1]])
                newScores1 = newScores1 + score1_size_weight
                newScores2 = newScores2 + score2_size_weight
            # new scores 1 stores the scores from one of the children while new scores 2 stores the other child's score
            # so if we do the argmax along the 0 axis then we get the best HMM for each query sequence after comparing the two children
            newScoresMax = np.argmax(np.array([newScores1, newScores2]), axis=0)
            # this new HMM is now used to get a list of new children from the original children list for every query sequence
            newChild = maxHMMChild_original[np.arange(newScores1.shape[0]), newScoresMax]
            # now the entries of maxHMM can be updated with the new children info since again
            # friendly reminder as stated above that maxHMM is an array where the ith element is the ith queries best HMM
            # if it's early stop and the argmax happens to be ourselves, then the bool_done array will already have been populated
            # for every query sequence  so we won't enter here since queryHMM.shape[0] would be zero and likewise in saveScoreFromBool
            # even if it did run
            if early_stop:
                maxHMM = np.argmax(scoresFull, axis=1)
            else:
                maxHMM[maxNotLeaf] = newChild

            scoreStep = scoresFull[np.arange(maxHMM.shape[0]), maxHMM]
            scoresRecord.append(np.copy(scoreStep))
            # Not sure what this is but it's unused
            # it's possible that previously the hmm indices were being recorded but now we just get them using argmax
            maxHMMlist.append(np.copy(maxHMM[:10]))
        else:
            done = True
    per_query_hmms_looked_at = bool_done.sum(axis=1)
    num_hmms_looked_at = np.where(bool_done.sum(axis=0) > 0, 1, 0).sum()
    num_hmms = len(bool_done[0])
    _LOG.info(f"num_hmms_looked at: {num_hmms_looked_at} hmms")
    _LOG.info(f"UPP would have looked at: {num_hmms} hmms")
    _LOG.info(f"num_iteration in hierarchical search: {num_iteration}")

    per_query_hmm_chosen= np.argmax(scoresFull, axis=1)
    # maxHMM_original = np.argmax(scores_original, axis=1)

    # print(f"maxHMM_guess: {maxHMM_guess}")
    # print(f"maxHMM_original: {maxHMM_original}")
    # so score.npy then is a 2d array where the rows are query sequences and the columns are hmms and it stores the hmmsearch scores inside
    np.save('%s/ensembleData/Searcher/scoreFiles/score.npy' % dirName, scoresFull)
    # Bool.npy is going to be a 2d array where array[a,b] = 1 if hmmsearch was run for query sequence a using hmm b
    np.save('%s/ensembleData/Searcher/scoreFiles/Bool.npy' % dirName, bool_done)

    return {
        "num_hmms_looked_at": num_hmms_looked_at,
        "per_query_hmms_looked_at": per_query_hmms_looked_at,
        "num_hmms": num_hmms,
        "per_query_hmm_chosen": per_query_hmm_chosen,
    }
    # this is also unused
    scoresRecord = np.array(scoresRecord)
    # this line and below are unused

    argDiff = np.argwhere((maxHMM_original-maxHMM_guess) != 0)[:, 0]
    argDiffTop = argDiff[treeSum[maxHMM_guess[argDiff]] == 999]
    print (np.mean(scores_original[argDiff, maxHMM_original[argDiff]] - scores_original[argDiff, maxHMM_guess[argDiff]]))
    print (argDiff.shape)
    print (np.unique(treeSum[maxHMM_guess[argDiff]], return_counts=True)[1][-3:])
    print (np.unique(treeSum[maxHMM_guess], return_counts=True)[1][-3:])
    bitscore_guess = scores_original[np.arange(Nquery), maxHMM_guess]
    bitscore_original = scores_original[np.arange(Nquery), maxHMM_original]
