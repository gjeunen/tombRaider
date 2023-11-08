#! /usr/bin/env python3

##################
# IMPORT MODULES #
##################
import collections
import rich
import rich.progress
import os
from Bio import pairwise2


########################
# tombRaider FUNCTIONS #
########################
def checkParamsNotNone(**kwargs):
    none_args = [arg_name for arg_name, arg_value in kwargs.items() if arg_value is None]
    return none_args

def freqToMemory(FREQ, pbar, progress_bar):
    '''
    Function parsing the frequency table input file into a dictionary
    '''
    freqInputDict = collections.defaultdict(dict)
    freqTotalCountDict = {}
    sampleNameList = []
    with open(FREQ, 'r') as freqFile:
        for line in freqFile:
            progress_bar.update(pbar, advance=len(line))
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
    return freqInputDict, freqTotalCountSortedDict, pbar, progress_bar

def zotuToMemory(ZOTU, freqTotalCountDict, pbar, progress_bar):
    '''
    Function parsing the ZOTU sequence input file into a dictionary
    '''
    seqInputDict = {}
    count = 0
    seqName = ''
    sequence = ''
    with open(ZOTU, 'r') as seqFile:
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
    for item in freqTotalCountDict:
        if item not in seqInputDict:
            seqInputDict[item] = ''
    return seqInputDict, pbar, progress_bar

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
    
########################
# tombRaider ALGORITHM #
########################
def taxonDependentCoOccurrenceAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, occurrence_type_, abundance, similarity, seqname, taxid, pident, qcov, eval_):
    '''
    The main function to identify and merge parent-child sequences using the taxon-dependent co-occurrence algorithm
    '''
    console = rich.console.Console(stderr=True, highlight=False)
    columns = [*rich.progress.Progress.get_default_columns(), rich.progress.TimeElapsedColumn()]

    # check if all parameters are provided
    missingArguments = checkParamsNotNone(frequency_input = frequency_input_, sequence_input = sequence_input_, taxonomy_input = taxonomy_input_, frequency_output = frequency_output_, sequence_output = sequence_output_, taxonomy_output = taxonomy_output_)
    if len(missingArguments) == 1:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]'--{''.join(missingArguments).replace('_', '-')}' parameter not specified, aborting analysis...[/]\n")
        exit()
    elif len(missingArguments) > 1:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]--{' and --'.join(missingArguments).replace('_', '-')} parameters not specified, aborting analysis...[/]\n")
        exit()

    # try reading in the files
    try:
        inputFilePaths = [frequency_input_, sequence_input_, taxonomy_input_]
        inputTotalFileSize = sum(os.path.getsize(inputFilePath) for inputFilePath in inputFilePaths)
        with rich.progress.Progress(*columns) as progress_bar:
            pbar = progress_bar.add_task(console = console, description="[cyan]|       Reading Files[/] |", total=inputTotalFileSize)
            freqInputDict, freqTotalCountDict, pbar, progress_bar = freqToMemory(frequency_input_, pbar, progress_bar)
            seqInputDict, pbar, progress_bar = zotuToMemory(sequence_input_, freqTotalCountDict, pbar, progress_bar)
            taxIdInputDict, taxQcovInputDict, taxPidentInputDict, taxTotalDict, pbar, progress_bar = taxToMemory(taxonomy_input_, freqTotalCountDict, seqname, taxid, pident, qcov, eval_, pbar, progress_bar)
    except TypeError as e:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]{e}, aborting analysis...[/]\n")
        exit()
    except FileNotFoundError as f:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]{f}, aborting analysis...[/]\n")
        exit()

    ## calculate number of unique combinations for progress bar (x = ((n*n)-n)/2)
    uniqueCombinations = int(((len(freqTotalCountDict) * len(freqTotalCountDict)) - len(freqTotalCountDict)) / 2)

    ## determine parent and child sequences
    newlyUpdatedCountDict = collections.defaultdict(list)
    combinedDict = collections.defaultdict(list)
    childParentComboDict = {}
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description="[cyan]|  Identify artefacts[/] |", total=uniqueCombinations)

        # 0. child needs to have lower abundance then parent
        for parent in range(len(freqTotalCountDict.keys())):
            parentName = list(freqTotalCountDict.keys())[parent]
            newParentDict = freqInputDict[parentName]
            for child in range(parent + 1, len(freqTotalCountDict.keys())):
                progress_bar.update(pbar, advance=1)
                childName = list(freqTotalCountDict.keys())[child]
                try:

                    # 1. check if childName already in childParentComboDict, skip if yes
                    if childName in childParentComboDict:
                        continue

                    # 2. check if BLAST taxonomic ID is matching between parent and child
                    taxIdParentSet = set(taxIdInputDict[parentName])
                    taxIdChildSet = set(taxIdInputDict[childName])
                    if not taxIdParentSet.intersection(taxIdChildSet):
                        continue

                    # 3. check BLAST quality on percent identity and query coverage are lower for child than parent
                    # only check top BLAST hit for now. Probably not accurate, so will need to be altered in future
                    if taxPidentInputDict[childName][0] > taxPidentInputDict[parentName][0] and taxQcovInputDict[childName][0] > taxPidentInputDict[parentName][0]:
                        continue

                    # 4. check co-occurrence pattern
                    # 4.1 check if child only appears in samples where parent is present
                    if occurrence_type_ == 'presence-absence':
                        positiveDetectionsChild = [k for k, v in freqInputDict[childName].items() if v >= int(abundance)]
                        positiveDetectionsParent = [k for k, v in freqInputDict[parentName].items() if v >= int(abundance)]
                        if not all(item in positiveDetectionsParent for item in positiveDetectionsChild):
                            continue

                    # 4.2 check if child only has lower abundance in samples compared to parent
                    elif occurrence_type_ == 'abundance':
                        count = 0
                        for item in freqInputDict[parentName]:
                            parentValue = freqInputDict[parentName][item]
                            childValue = freqInputDict[childName][item]
                            if parentValue < int(abundance):
                                parentValue = 0
                            if childValue < int(abundance):
                                childValue = 0
                            if parentValue < childValue:
                                count += 1
                        if count > 0:
                            continue
                    else:
                        console.print(f"[cyan]\n|               ERROR[/] | [bold yellow]'--occurrence-type' no specified as 'presence-absence' or 'abundance', aborting analysis...[/]\n")
                        exit()

                    # 5. check sequence similarity
                    alignment = pairwise2.align.globalxx(seqInputDict[parentName], seqInputDict[childName])
                    distanceCalculation = sum(1 for a, b in zip(alignment[0][0], alignment[0][1]) if a != b)
                    if 100 - (distanceCalculation/ max(len(seqInputDict[parentName]), len(seqInputDict[childName])) * 100) <= int(similarity):
                        continue

                    # if it passes all the checks, we need to determine how it can be combined --> several options
                    # first: if parent not identified as a child previously, we can combine child and parent data
                    if parentName not in childParentComboDict:
                        childParentComboDict[childName] = parentName
                        combinedDict[parentName].append(childName)
                        for item in newParentDict:
                            newValue = int(newParentDict[item]) + int(freqInputDict[childName][item])
                            newParentDict[item] = newValue

                    # second: childName not already identified as a child previously and parent identified as a child previously, add childName data to parent of parentName data
                    elif parentName in childParentComboDict:
                        combinedDict[childParentComboDict[parentName]].append(childName)
                        childParentComboDict[childName] = childParentComboDict[parentName]
                        for item in newlyUpdatedCountDict[childParentComboDict[parentName]]:
                            newValueGrandParent = int(newlyUpdatedCountDict[childParentComboDict[parentName]][item]) + int(freqInputDict[childName][item])
                            newlyUpdatedCountDict[childParentComboDict[parentName]][item] = newValueGrandParent
                except KeyError as k:
                    console.print(f"[cyan]|               ERROR[/] | [bold yellow]{k}, aborting analysis...[/]\n")
                    exit()
            if parentName not in childParentComboDict:
                newlyUpdatedCountDict[parentName] = newParentDict
                
    ## write updated frequency table to output
    count = 0
    with open(frequency_output_, 'w') as outfile:
        for item in newlyUpdatedCountDict:
            count += 1
            if count == 1:
                title = "\t".join(newlyUpdatedCountDict[item].keys())
                outfile.write(f'ID\t{title}\n')
            test = "\t".join(str(value) for value in newlyUpdatedCountDict[item].values())
            outfile.write(f'{item}\t{test}\n')

    ## write updated sequence file to output
    with open(sequence_output_, 'w') as seqoutfile:
        for item in newlyUpdatedCountDict:
            seqoutfile.write(f'>{item}\n{seqInputDict[item]}\n')

    ## write updated taxonomy file to output
    with open(taxonomy_output_, 'w') as taxoutfile:
        for item in newlyUpdatedCountDict:
            for subitem in taxTotalDict[item]:
                taxoutfile.write(f'{subitem}\n')

    ## write log
    console.print(f"[cyan]|                    [/] | [bold yellow][/]")
    console.print(f"[cyan]|  Summary Statistics[/] | [bold yellow][/]")
    console.print(f"[cyan]|     Total # of ASVs[/] | [bold yellow]{len(seqInputDict)}[/]")
    console.print(f"[cyan]|Total # of Artefacts[/] | [bold yellow]{len(childParentComboDict)} ({float('{:.2f}'.format(len(childParentComboDict) / len(seqInputDict) * 100))}%)[/]")
    console.print(f"[cyan]|   Parent-Child List[/] | [bold yellow][/]")
    for item in combinedDict:
        spaces = ' ' * (9 - len(item))
        if len(combinedDict[item]) == 1:
            console.print(f"[cyan]|    parent:{spaces}{item}[/] | [bold yellow]child:      {', '.join(combinedDict[item])}[/]")
        else:
            console.print(f"[cyan]|    parent:{spaces}{item}[/] | [bold yellow]children:   {', '.join(combinedDict[item])}[/]")

####################
# LULU ALTERNATIVE #
####################
def taxonIndependentCoOccurrenceAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, occurrence_type_, abundance, similarity):
    '''
    The function to identify and merge parent-child sequences based on the LULU algorithm (taxon-independent co-occurrence patterns)
    '''
    console = rich.console.Console(stderr=True, highlight=False)

####################
# MRCA ALTERNATIVE #
####################
def mostRecentCommonAncestorAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, occurrence_type_, abundance, similarity, seqname, taxid, pident, qcov, eval_):
    '''
    Thefunction to identify and merge parent-child sequences based on the taxonomic ID of sequences
    '''
    console = rich.console.Console(stderr=True, highlight=False)

###################
# MRCA CALCULATOR #
###################    
def mrcaCalculatorAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, occurrence_type_, abundance, similarity, seqname, taxid, pident, qcov, eval_):
    '''
    The function to calculate the Most Recent Common Ancestor from a standard BLAST output
    '''
    console = rich.console.Console(stderr=True, highlight=False)