#! /usr/bin/env python3

##################
# IMPORT MODULES #
##################
import collections


########################
# tombRaider FUNCTIONS #
########################
def freqToMemory(FREQ, pbar):
    '''
    Function parsing the frequency table input file into a dictionary
    '''
    freqInputDict = collections.defaultdict(dict)
    freqTotalCountDict = {}
    sampleNameList = []
    with open(FREQ, 'r') as freqFile:
        for line in freqFile:
            pbar.update(len(line))
            line = line.rstrip('\n')
            if line.startswith('#'):
                sampleNameList = line.split('\t')[1:]
            else:
                zotuName = line.split('\t')[0]
                sumCount = 0
                for i in range(len(sampleNameList)):
                    freqInputDict[zotuName][sampleNameList[i]] = int(line.split('\t')[i + 1])
                    sumCount += int(line.split('\t')[i + 1])
                freqTotalCountDict[zotuName] = sumCount
    freqTotalCountSortedDict = dict(sorted(freqTotalCountDict.items(), key = lambda x:x[1], reverse = True))
    return freqInputDict, freqTotalCountSortedDict, pbar

def zotuToMemory(ZOTU, pbar):
    '''
    Function parsing the ZOTU sequence input file into a dictionary
    '''
    seqInputDict = {}
    count = 0
    seqName = ''
    sequence = ''
    with open(ZOTU, 'r') as seqFile:
        for line in seqFile:
            pbar.update(len(line))
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
    return seqInputDict, pbar

def taxToMemory(TAX, pbar):
    '''
    Function parsing the BLAST taxonomy file
    For now it only takes in the top-BLAST-hit with the specific outfmt "6" structure
    '''
    taxIdInputDict = {}
    taxQcovInputDict = {}
    taxPidentInputDict = {}
    with open(TAX, 'r') as taxFile:
        for line in taxFile:
            pbar.update(len(line))
            line = line.rstrip('\n')
            seqName = line.split('\t')[0]
            taxAccession = line.split('\t')[1]
            taxQcov = int(line.split('\t')[5])
            taxPident = float(line.split('\t')[3])
            taxIdInputDict[seqName] = taxAccession
            taxQcovInputDict[seqName] = taxQcov
            taxPidentInputDict[seqName] = taxPident
    return taxIdInputDict, taxQcovInputDict, taxPidentInputDict, pbar


    
