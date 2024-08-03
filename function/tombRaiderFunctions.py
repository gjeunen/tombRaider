#! /usr/bin/env python3

##################
# IMPORT MODULES #
##################
import collections
import os
import numpy as np
import pandas as pd

########################
# tombRaider FUNCTIONS #
########################
def checkTaxonomyFiles(taxonomy_input_, blast_input_, bold_input_, sintax_input_, idtaxa_input_):
    taxonomy_file_type_mapping = {
    'taxonomy' : taxonomy_input_,
    'blast': blast_input_,
    'bold': bold_input_,
    'sintax': sintax_input_,
    'idtaxa': idtaxa_input_,
    }
    taxonomyFileType = None
    taxonomyInputFile = None
    for file_type, input_var in taxonomy_file_type_mapping.items():
        if input_var is not None and taxonomyFileType is None:
            taxonomyFileType = file_type
            taxonomyInputFile = input_var
        elif input_var is not None and taxonomyFileType is not None:
            taxonomyFileType = 'too-many-tax-files'
    return taxonomyInputFile, taxonomyFileType

def freqToMemory(frequency_input_, pbar, progress_bar, console, transpose_, omit_rows_, omit_columns_, sort_):
    '''
    Function parsing the frequency table input file into a dictionary
    '''
    frequencyTable = pd.read_csv(frequency_input_, sep='\t', index_col=0)
    progress_bar.update(pbar, advance=os.path.getsize(frequency_input_))
    if omit_rows_ != None:
        rowList = omit_rows_.split(',')
        try:
            frequencyTable = frequencyTable.drop(rowList)
        except KeyError as k:
            console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]{k}, aborting analysis...[/]\n")
            exit()
    if omit_columns_ != None:
        colList = omit_columns_.split(',')
        try:
            frequencyTable = frequencyTable.drop(colList, axis = 1)
        except KeyError as k:
            console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]{k}, aborting analysis...[/]\n")
            exit()
    if transpose_:
        frequencyTable = frequencyTable.transpose()
    sortOptions = {
        'total read count': frequencyTable.sum(axis = 1),
        'average read count': frequencyTable.mean(axis = 1),
        'detections': (frequencyTable > 0).sum(axis = 1)
    }
    if sort_ in sortOptions:
        frequencyTable = frequencyTable.loc[sortOptions[sort_].sort_values(ascending = False).index]
    elif sort_ == None:
        sort_ = 'None'
    else:
        console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]option for '--sort' not identified, aborting analysis...[/]\n")
        exit()
    return frequencyTable, pbar, progress_bar

def zotuToMemory(sequence_input_, frequencyTable, pbar, progress_bar):
    '''
    Function parsing the ZOTU sequence input file into a dictionary
    '''
    seqInputDict = {}
    count = 0
    seqName = ''
    sequence = ''
    with open(sequence_input_, 'r') as seqFile:
        for line in seqFile:
            progress_bar.update(pbar, advance=len(line))
            line = line.rstrip('\n')
            count += 1
            if line.startswith('>'):
                if count > 1:
                    seqInputDict[seqName] = sequence
                    seqName = ''
                    sequence = ''
                seqName = line.lstrip('>')
            else:
                sequence += line
    seqInputDict[seqName] = sequence
    return seqInputDict, pbar, progress_bar

def blastToMemory(taxonomyInputFile, blast_format_, use_accession_id_, seqInputDict, pbar, progress_bar, console):
    '''
    Function parsing the BLAST taxonomy file
    For now it only takes in the specific outfmt "6" structure
    '''
    taxonomyForMissingSeqs = []
    if use_accession_id_:
        neededBlastInfo = {'qaccver': None,
                        'qcovs': None,
                        'pident': None,
                        'saccver': None,
                        'evalue': None,
                        'length': None,
                        'gapopen': None,
                        'mismatch': None,
                        'staxid': None,
        }
    else:
        neededBlastInfo = {'qaccver': None,
                        'qcovs': None,
                        'pident': None,
                        'saccver': None,
                        'evalue': None,
                        'length': None,
                        'gapopen': None,
                        'mismatch': None,
                        'staxid': None,
        }
    blastFormattingList = blast_format_.split(' ')
    if blastFormattingList[0] != '6':
        console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]blast format identified as '{blastFormattingList[0]}', only format '6' is supported, aborting analysis...[/]\n")
        exit()
    for item in neededBlastInfo:
        if item in blastFormattingList:
            neededBlastInfo[item] = blastFormattingList.index(item) - 1
        else:
            console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]'{item}' not found in BLAST input file, aborting analysis...[/]\n")
            exit()
    taxIdInputDict = collections.defaultdict(list)
    taxPidentInputDict = collections.defaultdict(list)
    rawTaxDict = collections.defaultdict(list)
    taxEvalInputDict = collections.defaultdict(list)
    with open(taxonomyInputFile, 'r') as taxFile:
        for line in taxFile:
            progress_bar.update(pbar, advance=len(line))
            line = line.rstrip('\n')
            seqName = line.split('\t')[neededBlastInfo['qaccver']]
            taxQcov = int(line.split('\t')[neededBlastInfo['qcovs']])
            if taxQcov == 100:
                taxPident = float(line.split('\t')[neededBlastInfo['pident']])
            else:
                try:
                    taxPident = float(100 * ((int(line.split('\t')[neededBlastInfo['length']]) - int(line.split('\t')[neededBlastInfo['mismatch']]) - int(line.split('\t')[neededBlastInfo['gapopen']])) / len(seqInputDict[seqName])))
                except KeyError:
                    taxonomyForMissingSeqs.append(seqName)
            if use_accession_id_:
                taxID = line.split('\t')[neededBlastInfo['saccver']]
            else:
                taxID = line.split('\t')[neededBlastInfo['staxid']]
            evalNumber = float(line.split('\t')[neededBlastInfo['evalue']])
            rawTaxDict[seqName].append(line)
            if all(taxPident >= item for item in taxPidentInputDict[seqName]) and taxID not in taxIdInputDict[seqName]:
                taxIdInputDict[seqName].append(taxID)
                taxPidentInputDict[seqName].append(taxPident)
                taxEvalInputDict[seqName].append(evalNumber)
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, taxonomyForMissingSeqs, pbar, progress_bar

def boldToMemory(taxonomyInputFile, bold_format_, pbar, progress_bar, console):
    '''
    Function parsing the BOLD taxonomy file
    Either the summary report (website) or the complete BOLD IDE results (BOLDIGGER)
    '''
    taxIdInputDict = collections.defaultdict(list)
    taxPidentInputDict = collections.defaultdict(list)
    rawTaxDict = collections.defaultdict(list)
    if bold_format_ == 'summary':
        with open(taxonomyInputFile, 'r') as taxFile:
            for line in taxFile:
                progress_bar.update(pbar, advance = len(line))
                seqName = line.split('\t')[0]
                taxID = line.split('\t')[1]
                try:
                    taxPident = float(line.split('\t')[3])
                except IndexError:
                    taxPident = 0.0
                taxIdInputDict[seqName].append(taxID)
                taxPidentInputDict[seqName].append(taxPident)
                rawTaxDict[seqName].append(line)
    elif bold_format_ == 'complete':
        with open(taxonomyInputFile, 'r') as taxFile:
            next(taxFile)
            for line in taxFile:
                progress_bar.update(pbar, advance = len(line))
                if line.split('\t')[0] != '':
                    seqName = line.split('\t')[0].lstrip('>')
                try:
                    taxPident = float(line.split('\t')[8])
                    taxID = ','.join(line.split('\t')[1:8])
                except ValueError:
                    taxPident = 0.0
                    taxID = 'not assigned'
                rawTaxDict[seqName].append(line)
                if all(taxPident >= item for item in taxPidentInputDict[seqName]) and taxID not in taxIdInputDict[seqName]:
                    taxIdInputDict[seqName].append(taxID)
                    taxPidentInputDict[seqName].append(taxPident)
    else:
        console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]'{bold_format_}' not identified, aborting analysis...[/]\n")
        exit()
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar


def sintaxToMemory(taxonomyInputFile, sintax_threshold_, pbar, progress_bar):
    '''
    Function parsing the SINTAX taxonomy file
    For now, the file needs to be tab-delimited
    '''
    taxIdInputDict = collections.defaultdict(list)
    taxPidentInputDict = collections.defaultdict(list)
    rawTaxDict = collections.defaultdict(list)
    thresholdSplit = 1
    if sintax_threshold_:
        thresholdSplit = 3
    with open(taxonomyInputFile, 'r') as taxFile:
        for line in taxFile:
            progress_bar.update(pbar, advance=len(line))
            seqName = line.split('\t')[0]
            taxID = line.split('\t')[thresholdSplit].split(',')[-1].split(':')[1].split('(')[0]
            taxPident = float(line.split('\t')[thresholdSplit].split(',')[-1].split(':')[1].split('(')[1].rstrip(')'))
            taxIdInputDict[seqName].append(taxID)
            taxPidentInputDict[seqName].append(taxPident)
            rawTaxDict[seqName].append(line)
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar

def idtaxaToMemory(taxonomyInputFile, pbar, progress_bar):
    '''
    Function parsing the IDTAXA taxonomy file
    '''
    taxIdInputDict = collections.defaultdict(list)
    taxPidentInputDict = collections.defaultdict(list)
    rawTaxDict = collections.defaultdict(list)
    with open(taxonomyInputFile, 'r') as taxFile:
        for line in taxFile:
            progress_bar.update(pbar, advance=len(line))
            seqName = line.split('\t')[0]
            taxID = line.split('\t')[1].split('; ')[-1].split(' (')[0]
            taxPident = float(line.split('\t')[1].split('; ')[-1].split(' (')[1].split('%)')[0])
            taxIdInputDict[seqName].append(taxID)
            taxPidentInputDict[seqName].append(taxPident)
            rawTaxDict[seqName].append(line)

    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar

def identifyTaxonomyInput(taxonomyInputFile):
    '''
    identify the taxonomy file type when "--taxonomy-input" was given as parameter
    '''
    bold_format_ = None
    with open(taxonomyInputFile, 'r') as infile:
        firstLine = infile.readline().rstrip('\n')
        if firstLine.startswith('Query ID'):
            taxonomyFileType = 'bold'
            bold_format_ = 'summary'
        elif firstLine.startswith('You'):
            taxonomyFileType = 'bold'
            bold_format_ = 'complete'
        elif '%' in firstLine:
            taxonomyFileType = 'idtaxa'
        elif firstLine.split('\t')[2] == '+' or firstLine.split('\t')[2] == '-':
            taxonomyFileType = 'sintax'
        else:
            taxonomyFileType = 'blast'
    return taxonomyFileType, bold_format_

def taxonomyToMemory(taxonomyInputFile, taxonomyFileType, blast_format_, use_accession_id_, bold_format_, sintax_threshold_, seqInputDict, pbar, progress_bar, console):
    '''
    '''
    taxonomyForMissingSeqs = []
    if taxonomyFileType == 'taxonomy':
        taxonomyFileType, bold_format_ = identifyTaxonomyInput(taxonomyInputFile)
    if taxonomyFileType == 'blast':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, taxonomyForMissingSeqs, pbar, progress_bar = blastToMemory(taxonomyInputFile, blast_format_, use_accession_id_, seqInputDict, pbar, progress_bar, console)
    elif taxonomyFileType == 'sintax':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = sintaxToMemory(taxonomyInputFile, sintax_threshold_, pbar, progress_bar)
    elif taxonomyFileType == 'bold':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = boldToMemory(taxonomyInputFile, bold_format_, pbar, progress_bar, console)
    elif taxonomyFileType == 'idtaxa':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = idtaxaToMemory(taxonomyInputFile, pbar, progress_bar)
    else:
        console.print("[cyan]|               ERROR[/] | [bold yellow]taxonomy input file type not recognised, aborting analysis...[/]\n")
        exit()

    return taxIdInputDict, taxPidentInputDict, taxTotalDict, taxonomyFileType, taxonomyForMissingSeqs, pbar, progress_bar

def fillOutTaxonomyFiles(taxIdInputDict, taxPidentInputDict, taxTotalDict, frequencyTable):
    '''
    fill out missing taxonomic IDs
    '''
    for item in frequencyTable.index.tolist():
        if item not in taxIdInputDict:
            taxIdInputDict[item].append('not assigned')
            taxPidentInputDict[item].append(0.0)
            taxTotalDict[item].append('not assigned')
    return taxIdInputDict, taxPidentInputDict, taxTotalDict
    
def alignmentToMemory(alignment_input_, pbar, progress_bar):
    """
    Parses a Nexus alignment file and returns a dictionary with sequence IDs as keys and sequences as values.
    """
    sequences = {}
    in_matrix = False
    with open(alignment_input_, 'r') as file:
        for line in file:
            progress_bar.update(pbar, advance=len(line))
            line = line.strip()
            if line.startswith('MATRIX'):
                in_matrix = True
                continue
            if in_matrix:
                if line == ';':
                    break
                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    seq = ''.join(parts[1:])
                    if seq_id in sequences:
                        sequences[seq_id] += seq
                    else:
                        sequences[seq_id] = seq
    return sequences, pbar, progress_bar

def verifySequences(seqInputDict, frequencyTable):
    '''
    function to verify all sequences are present in seqInputDict
    '''
    seqVerification = []
    for item in frequencyTable.index.tolist():
        if item not in seqInputDict:
            seqVerification.append(item)
    return seqVerification

def verifyAlignment(alignmentInputDict, seqInputDict):
    '''
    function to verify if alignment includes:
        1. all sequences that are in --sequence-input
        2. sequences with equal length
        3. sequences without gaps are identical to --sequence-input
    '''
    alignmentVerification = {}
    missingSeqs = list(set(seqInputDict.keys()) - set(alignmentInputDict.keys()))
    equalLength = all(len(v) == len(next(iter(alignmentInputDict.values()))) for v in alignmentInputDict.values())
    identicalSeqs = all(k in alignmentInputDict and v == alignmentInputDict[k].replace('-', '') for k, v in seqInputDict.items()) and len(seqInputDict) == len(alignmentInputDict)
    if len(missingSeqs) > 0:
        alignmentVerification['not all sequences in alignment'] = missingSeqs
    if equalLength != True:
        alignmentVerification["sequences not identical to '--sequence-input'"] = equalLength
    if identicalSeqs != True:
        alignmentVerification["sequences not identical to '--sequence-input'"] = identicalSeqs
    return alignmentVerification

def removeNegativeSamples(negative, frequencyTableSubset):
    '''
    function to remove negative samples before the algorithm
    '''
    negativeList = negative.split('+')
    for item in negativeList:
        if '*' not in item:
            frequencyTableSubset = frequencyTableSubset.drop(item, axis = 1)
        elif item.startswith('*') and item.endswith('*'):
            itemMatch = item.rstrip('*').lstrip('*')
            droppedColumns = frequencyTableSubset.filter(like = itemMatch).columns
            frequencyTableSubset = frequencyTableSubset.drop(columns = droppedColumns)
        elif item.startswith('*'):
            itemMatch = item.lstrip('*')
            droppedColumns = frequencyTableSubset.filter(regex = f'{itemMatch}$', axis = 1).columns
            frequencyTableSubset = frequencyTableSubset.drop(columns = droppedColumns)
        elif item.endswith('*'):
            itemMatch = item.rstrip('*')
            droppedColumns = frequencyTableSubset.filter(regex = f'^{itemMatch}', axis = 1).columns
            frequencyTableSubset = frequencyTableSubset.drop(columns = droppedColumns)
    return frequencyTableSubset

def passingFunction(*args, **kwargs):
    '''
    function to skip a step in the analysis
    '''
    pass

def pseudogeneIdentificationFunction(alignmentInputDict, orf_, parentID, pseudogeneDict):
    '''
    function to read in a sequence and identify a stop codon based on the starting position
    '''
    stop_codons = {"TAA", "TAG", "TGA"}
    for i in range(orf_ - 1, len(alignmentInputDict[parentID]), 3):
        codon = alignmentInputDict[parentID][i:i+3]
        if codon in stop_codons:
            pseudogeneDict[parentID] = 1
            return pseudogeneDict
    return pseudogeneDict

def taxidIdentificationFunction(logDict, childID, parentID, taxIdInputDict):
    '''
    '''
    taxIDparent = taxIdInputDict[parentID]
    taxIDchild = taxIdInputDict[childID]
    if set(taxIDparent) & set(taxIDchild):
        logDict[childID][parentID].append(f'taxonomic IDs are matching between {parentID} and {childID} ({list(set(taxIdInputDict[parentID]) & set(taxIdInputDict[childID]))[0]}), continuing inspection')
        return True
    return False

def taxqualIdentificationFunction(logDict, childID, parentID, taxPidentInputDict):
    '''
    '''
    taxPidentChild = taxPidentInputDict[childID][0]
    taxPidentParent = taxPidentInputDict[parentID][0]
    if taxPidentChild > taxPidentParent:
        logDict[childID][parentID].append(f'taxonomic similarity score higher for {childID} ({taxPidentInputDict[childID][0]}) than {parentID} ({taxPidentInputDict[parentID][0]}), aborting inspection...')
        return False
    logDict[childID][parentID].append(f'taxonomic similarity score not higher for {childID} ({taxPidentInputDict[childID][0]}) than {parentID} ({taxPidentInputDict[parentID][0]}), continuing inspection')
    return True

def cooccurIdentificationFunction(console, frequencyTableSubset, child, detection_threshold_, parent, occurrence_type_, occurrence_ratio_, logDict, childID, parentID):
    '''
    '''
    positiveDetectionsChild = frequencyTableSubset.iloc[child][frequencyTableSubset.iloc[child] >= detection_threshold_].index.tolist()
    positiveDetectionsParent = frequencyTableSubset.iloc[parent][frequencyTableSubset.iloc[parent] >= detection_threshold_].index.tolist()
    if occurrence_type_ == 'presence-absence':
        missingCount = len(set(positiveDetectionsChild) - set(positiveDetectionsParent))
    elif occurrence_type_ == 'abundance':
        missingCount = ((frequencyTableSubset.iloc[child] >= detection_threshold_) & (frequencyTableSubset.iloc[child] > frequencyTableSubset.iloc[parent])).sum()
    else:
        console.print(f"[cyan]\n|               ERROR[/] | [bold yellow]'--occurrence-type' not specified as 'presence-absence' or 'abundance', aborting analysis...[/]\n")
        exit()
    if occurrence_ratio_.split(';')[0].upper() == 'COUNT':
        if int(occurrence_ratio_.split(';')[1]) < missingCount:
            logDict[childID][parentID].append(f'co-occurrence ratio (method: {occurrence_ratio_.split(";")[0]}) not met: {childID} found {missingCount} times without {parentID} (threshold: {int(occurrence_ratio_.split(";")[1])}), aborting inspection...')
            return False
        logDict[childID][parentID].append(f'co-occurrence ratio (method: {occurrence_ratio_.split(";")[0]}) met: {childID} found {missingCount} times without {parentID} (threshold: {int(occurrence_ratio_.split(";")[1])}), continuing inspection')
        return True
    elif occurrence_ratio_.split(';')[0].upper() == 'GLOBAL':
        if 1 - (missingCount / len(frequencyTableSubset.index.tolist())) < float(occurrence_ratio_.split(';')[1]):
            logDict[childID][parentID].append(f'co-occurrence ratio (method: {occurrence_ratio_.split(";")[0]}) not met: {childID} observed in {float("{:.2f}".format(1 - (missingCount / len(frequencyTableSubset.index.tolist()))))}% of samples without a positive detection of {parentID} (threshold: {int(occurrence_ratio_.split(";")[1])}), aborting inspection...')
            return False
        logDict[childID][parentID].append(f'co-occurrence ratio (method: {occurrence_ratio_.split(";")[0]}) met: {childID} observed in {float("{:.2f}".format(1 - (missingCount / len(frequencyTableSubset.index.tolist()))))}% of samples without a positive detection of {parentID} (threshold: {int(occurrence_ratio_.split(";")[1])}), continuing inspection')
        return True
    elif occurrence_ratio_.split(';')[0].upper() == 'LOCAL':
        if 1 - (missingCount / (len(positiveDetectionsParent) + missingCount)) < float(occurrence_ratio_.split(';')[1]):
            logDict[childID][parentID].append(f'co-occurrence ratio (method: {occurrence_ratio_.split(";")[0]}) not met: {childID} observed in {float("{:.2f}".format(1 - (missingCount / len(frequencyTableSubset.index.tolist()))))}% of samples without a positive detection of {parentID} (threshold: {int(occurrence_ratio_.split(";")[1])}), aborting inspection...')
            return False
        logDict[childID][parentID].append(f'co-occurrence ratio (method: {occurrence_ratio_.split(";")[0]}) met: {childID} observed in {float("{:.2f}".format(1 - (missingCount / len(frequencyTableSubset.index.tolist()))))}% of samples without a positive detection of {parentID} (threshold: {int(occurrence_ratio_.split(";")[1])}), continuing inspection')
        return True
    else:
        console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]--occurrence-ratio parameter not identified, aborting analysis...[/]\n")
        exit()
    
def seqSimIdentificationFunction(console, seqParent, seqChild, alignmentInputDict, childID, parentID, calculate_pairwise_, pairwise_alignment_, similarity_, logDict):
    '''
    '''
    if len(alignmentInputDict) == 0 or calculate_pairwise_:
        if pairwise_alignment_ == 'global':
            alignmentSeq1, alignmentSeq2 = needleman_wunsch(seqParent, seqChild)
        elif pairwise_alignment_ == 'local':
            alignmentSeq1, alignmentSeq2, max_score = smith_waterman(seqParent, seqChild)
        else:
            console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]--occurrence-ratio parameter not identified, aborting analysis...[/]\n")
            exit()
    else:
        alignmentSeq1 = alignmentInputDict[parentID]
        alignmentSeq2 = alignmentInputDict[childID]
    distanceCalculation = sum(1 for a, b in zip(alignmentSeq1, alignmentSeq2) if a != b)
    if 100 - (distanceCalculation/ max(len(seqParent), len(seqChild)) * 100) <= int(similarity_):
        logDict[childID][parentID].append(f'sequence similarity ratio not met: {float("{:.2f}".format(100 - (distanceCalculation/ max(len(seqParent), len(seqChild)) * 100)))}% (threshold: {similarity_}), aborting inspection...')
        return False
    logDict[childID][parentID].append(f'sequence similarity ratio met: {float("{:.2f}".format(100 - (distanceCalculation/ max(len(seqParent), len(seqChild)) * 100)))}% (threshold: {similarity_}), continuing inspection')
    return True

def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-5, gap_penalty=-5):
    '''
    local alignment function in base python (except Numpy) based on the Smith-Waterman algorithm
    '''
    len_seq1, len_seq2 = len(seq1), len(seq2)
    # Initialize the scoring matrix with zeros
    score_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1), dtype=int)
    # Fill in the scoring matrix
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            match = score_matrix[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1, j] + gap_penalty
            insert = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(0, match, delete, insert)
    # Find the maximum score in the matrix
    max_i, max_j = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
    max_score = score_matrix[max_i, max_j]
    # Trace back to find the alignment
    alignment_seq1, alignment_seq2 = "", ""
    while score_matrix[max_i, max_j] != 0:
        if max_i > 0 and max_j > 0 and score_matrix[max_i, max_j] == score_matrix[max_i - 1, max_j - 1] + (match_score if seq1[max_i - 1] == seq2[max_j - 1] else mismatch_penalty):
            alignment_seq1 = seq1[max_i - 1] + alignment_seq1
            alignment_seq2 = seq2[max_j - 1] + alignment_seq2
            max_i -= 1
            max_j -= 1
        elif max_i > 0 and score_matrix[max_i, max_j] == score_matrix[max_i - 1, max_j] + gap_penalty:
            alignment_seq1 = seq1[max_i - 1] + alignment_seq1
            alignment_seq2 = "-" + alignment_seq2
            max_i -= 1
        elif max_j > 0 and score_matrix[max_i, max_j] == score_matrix[max_i, max_j - 1] + gap_penalty:
            alignment_seq1 = "-" + alignment_seq1
            alignment_seq2 = seq2[max_j - 1] + alignment_seq2
            max_j -= 1
    return alignment_seq1, alignment_seq2, max_score

def needleman_wunsch(seq1, seq2, gap_penalty=-1, match_score=2, mismatch_penalty=-1):
    '''
    global alignment function in base python (except Numpy) based on the Needleman-Wunsch algorithm
    '''
    m, n = len(seq1), len(seq2)
    # Initialize the dynamic programming table
    F = np.zeros((m + 1, n + 1), dtype=np.int32)
    F[:, 0] = gap_penalty * np.arange(m + 1)
    F[0, :] = gap_penalty * np.arange(n + 1)
    # Fill in the dynamic programming table using vectorized operations
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score1 = F[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            score2 = F[i-1, j] + gap_penalty
            score3 = F[i, j-1] + gap_penalty
            F[i, j] = np.max([score1, score2, score3])
    # Trace back through the dynamic programming table to find the optimal alignment
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i, j] == F[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and F[i, j] == F[i-1, j] + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2