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
def checkParamsNotNone(**kwargs):
    none_args = [arg_name for arg_name, arg_value in kwargs.items() if arg_value is None]
    return none_args

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

def blastToMemory(taxonomyInputFile, frequencyTable, seqname, taxid, pident, qcov, eval_, pbar, progress_bar):
    '''
    Function parsing the BLAST taxonomy file
    For now it only takes in the specific outfmt "6" structure
    '''
    taxIdInputDict = collections.defaultdict(list)
    taxQcovInputDict = collections.defaultdict(list)
    taxPidentInputDict = collections.defaultdict(list)
    taxTotalDict = collections.defaultdict(list)
    taxEvalInputDict = collections.defaultdict(list)
    with open(taxonomyInputFile, 'r') as taxFile:
        for line in taxFile:
            progress_bar.update(pbar, advance=len(line))
            line = line.rstrip('\n')
            seqName = line.split('\t')[seqname]
            taxQcov = int(line.split('\t')[qcov])
            if taxQcov == 100:
                taxPident = float(line.split('\t')[pident])
            else:
                print()
            taxID = line.split('\t')[taxid]
            evalNumber = float(line.split('\t')[eval_])
            taxTotalDict[seqName].append(line)
            print(taxPident, taxQcov, taxPident/taxQcov*100)
            if all(evalNumber <= item for item in taxEvalInputDict[seqName]) and taxID not in taxIdInputDict[seqName]:
                taxIdInputDict[seqName].append(taxID)
                taxQcovInputDict[seqName].append(taxQcov)
                taxPidentInputDict[seqName].append(taxPident)
                taxEvalInputDict[seqName].append(evalNumber)
    for item in frequencyTable.index.tolist():
        if item not in taxIdInputDict:
            taxIdInputDict[item] = ''
            taxQcovInputDict[item] = ''
            taxPidentInputDict[item] = ''
            taxTotalDict[item] = ''
    return taxIdInputDict, taxQcovInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar

def boldToMemory():
    '''
    '''

def sintaxToMemory():
    '''
    '''

def idtaxaToMemory():
    '''
    '''

def taxonomyProcessing(taxonomyFileType, taxonomyInputFile, frequencyTable, seqname, taxid, pident, qcov, eval_, pbar, progress_bar, console):
    '''
    Function parsing taxonomy input file 100*(186/198)
    '''
    if taxonomyFileType == 'blast':
        taxIdInputDict, taxQcovInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = blastToMemory(taxonomyInputFile, frequencyTable, seqname, taxid, pident, qcov, eval_, pbar, progress_bar)
    elif taxonomyFileType == 'bold':
        test = boldToMemory(frequencyTable)
    elif taxonomyFileType == 'sintax':
        test = sintaxToMemory(seqname)
    elif taxonomyFileType == 'idtaxa':
        test = idtaxaToMemory(eval_)

def taxToMemory(TAX, freqTotalCountDict, seqname, taxid, pident, qcov, eval_, pbar, progress_bar):
    '''
    Function parsing the BLAST taxonomy file
    For now it only takes in the specific outfmt "6" structure
    '''
    taxIdInputDict = collections.defaultdict(list)
    taxQcovInputDict = collections.defaultdict(list)
    taxPidentInputDict = collections.defaultdict(list)
    taxTotalDict = collections.defaultdict(list)
    taxEvalInputDict = collections.defaultdict(list)
    with open(TAX, 'r') as taxFile:
        for line in taxFile:
            progress_bar.update(pbar, advance=len(line))
            line = line.rstrip('\n')
            seqName = line.split('\t')[seqname]
            taxQcov = int(line.split('\t')[qcov])
            taxPident = float(line.split('\t')[pident])
            taxID = line.split('\t')[taxid]
            evalNumber = float(line.split('\t')[eval_])
            taxTotalDict[seqName].append(line)
            if all(evalNumber <= item for item in taxEvalInputDict[seqName]) and taxID not in taxIdInputDict[seqName]:
                taxIdInputDict[seqName].append(taxID)
                taxQcovInputDict[seqName].append(taxQcov)
                taxPidentInputDict[seqName].append(taxPident)
                taxEvalInputDict[seqName].append(evalNumber)
    for item in freqTotalCountDict:
        if item not in taxIdInputDict:
            taxIdInputDict[item] = ''
            taxQcovInputDict[item] = ''
            taxPidentInputDict[item] = ''
            taxTotalDict[item] = ''
    return taxIdInputDict, taxQcovInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar

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