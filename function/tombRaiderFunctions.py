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
def checkTaxonomyFiles(blast_input_, bold_input_, sintax_input_, idtaxa_input_):
    taxonomy_file_type_mapping = {
    'blast': blast_input_,
    'bold': bold_input_,
    'sintax': sintax_input_,
    'idtaxa': idtaxa_input_,
    }
    taxonomyFileType = None
    taxonomyInputFile = None
    for file_type, input_var in taxonomy_file_type_mapping.items():
        if input_var is not None:
            taxonomyFileType = file_type
            taxonomyInputFile = input_var
            break
    return taxonomyInputFile, taxonomyFileType

def freqToMemory(frequency_input_, pbar, progress_bar, console, taxa_are_rows_, omit_rows_, omit_columns_, sort_):
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
    if taxa_are_rows_ != True:
        frequencyTable = frequencyTable.transpose()
    sortOptions = {
        'total read count': frequencyTable.sum(axis = 1),
        'average read count': frequencyTable.mean(axis = 1),
        'detections': (frequencyTable > 0).sum(axis = 1),
    }
    if sort_ in sortOptions:
        #freqTotalCountSortedDict = sortOptions[sort_].sort_values(ascending = False).to_dict()
        frequencyTable = frequencyTable.loc[sortOptions[sort_].sort_values(ascending = False).index]
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
    for item in frequencyTable.index.tolist():
        if item not in seqInputDict:
            seqInputDict[item] = ''
    return seqInputDict, pbar, progress_bar

def blastToMemory(taxonomyInputFile, frequencyTable, blast_format_, use_accession_id_, seqInputDict, pbar, progress_bar, console):
    '''
    Function parsing the BLAST taxonomy file
    For now it only takes in the specific outfmt "6" structure
    '''
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
                taxPident = float(100 * ((int(line.split('\t')[neededBlastInfo['length']]) - int(line.split('\t')[neededBlastInfo['mismatch']]) - int(line.split('\t')[neededBlastInfo['gapopen']])) / len(seqInputDict[seqName])))
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
    for item in frequencyTable.index.tolist():
        if item not in taxIdInputDict:
            taxIdInputDict[item].append('not assigned')
            taxPidentInputDict[item].append(0.0)
            rawTaxDict[item].append('not assigned')
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar

def boldToMemory(taxonomyInputFile, frequencyTable, bold_format_, pbar, progress_bar, console):
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
        for item in frequencyTable.index.tolist():
            if item not in taxIdInputDict:
                taxIdInputDict[item].append('not assigned')
                taxPidentInputDict[item].append(0.0)
                rawTaxDict[item].append('not assigned')
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
        for item in frequencyTable.index.tolist():
            if item not in taxIdInputDict:
                taxIdInputDict[item].append('not assigned')
                taxPidentInputDict[item].append(0.0)
                rawTaxDict[item].append('not assigned')
    else:
        console.print(f"\n[cyan]|               ERROR[/] | [bold yellow]'{bold_format_}' not identified, aborting analysis...[/]\n")
        exit()
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar


def sintaxToMemory(taxonomyInputFile, frequencyTable, sintax_threshold_, pbar, progress_bar):
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
    for item in frequencyTable.index.tolist():
        if item not in taxIdInputDict:
            taxIdInputDict[item].append('not assigned')
            taxPidentInputDict[item].append(0.0)
            rawTaxDict[item].append('not assigned')
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar

def idtaxaToMemory(taxonomyInputFile, frequencyTable, pbar, progress_bar):
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
    for item in frequencyTable.index.tolist():
        if item not in taxIdInputDict:
            taxIdInputDict[item].append('not assigned')
            taxPidentInputDict[item].append(0.0)
            rawTaxDict[item].append('not assigned')
    return taxIdInputDict, taxPidentInputDict, rawTaxDict, pbar, progress_bar

def taxonomyToMemory(taxonomyInputFile, taxonomyFileType, frequencyTable, blast_format_, use_accession_id_, bold_format_, sintax_threshold_, seqInputDict, pbar, progress_bar, console):
    '''
    '''
    if taxonomyFileType == 'blast':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = blastToMemory(taxonomyInputFile, frequencyTable, blast_format_, use_accession_id_, seqInputDict, pbar, progress_bar, console)
    elif taxonomyFileType == 'sintax':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = sintaxToMemory(taxonomyInputFile, frequencyTable, sintax_threshold_, pbar, progress_bar)
    elif taxonomyFileType == 'bold':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = boldToMemory(taxonomyInputFile, frequencyTable, bold_format_, pbar, progress_bar, console)
    elif taxonomyFileType == 'idtaxa':
        taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = idtaxaToMemory(taxonomyInputFile, frequencyTable, pbar, progress_bar)
    else:
        console.print("[cyan]|               ERROR[/] | [bold yellow]option for '--method' not identified, aborting analysis...[/]\n")
        exit()

    return taxIdInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar

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