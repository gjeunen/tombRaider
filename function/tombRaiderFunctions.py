#! /usr/bin/env python3

##################
# IMPORT MODULES #
##################
import collections
import rich
import rich.progress
import os
from Bio import pairwise2
import sys
import copy
import datetime


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
    return freqInputDict, freqTotalCountSortedDict, sampleNameList, pbar, progress_bar

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
def taxonDependentCoOccurrenceAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, condensed_log_, detailed_log_, occurrence_type_, detection_threshold_, similarity, negative, ratio, seqname, taxid, pident, qcov, eval_):
    '''
    The main function to identify and merge parent-child sequences using the taxon-dependent co-occurrence algorithm
    '''
    console = rich.console.Console(stderr=True, highlight=False)
    columns = [*rich.progress.Progress.get_default_columns(), rich.progress.TimeElapsedColumn()]
    startTime = datetime.datetime.now()
    formattedTime = startTime.strftime("%Y-%m-%d %H:%M:%S")
    commandLineInput = ' '.join(sys.argv[1:])

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
            freqInputDict, freqTotalCountDict, sampleNameList, pbar, progress_bar = freqToMemory(frequency_input_, pbar, progress_bar)
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

    ## get list of samples and exclude if {negative} != None
    fullNegativeList = []
    if negative == None:
        freqInputDictSubset = copy.deepcopy(freqInputDict)
    else:
        negativeList = negative.split('+')
        for item in negativeList:
            if item.startswith('*') and item.endswith('*'):
                itemMatch = item.rstrip('*').lstrip('*')
                for sampleName in sampleNameList:
                    if itemMatch in sampleName:
                        fullNegativeList.append(sampleName)
            elif item.startswith('*'):
                itemMatch = item.lstrip('*')
                for sampleName in sampleNameList:
                    if sampleName.endswith(itemMatch):
                        fullNegativeList.append(sampleName)
            elif item.endswith('*'):
                itemMatch = item.rstrip('*')
                for sampleName in sampleNameList:
                    if sampleName.startswith(itemMatch):
                        fullNegativeList.append(sampleName)
            else:
                for sampleName in sampleNameList:
                    if sampleName == item:
                        fullNegativeList.append(sampleName)
            freqInputDictSubset = copy.deepcopy(freqInputDict)
            for item in freqInputDictSubset:
                for sampleToRemove in fullNegativeList:
                    if sampleToRemove in freqInputDictSubset[item]:
                        del freqInputDictSubset[item][sampleToRemove]
            
    ## determine parent and child sequences
    newlyUpdatedCountDict = collections.defaultdict(list)
    combinedDict = collections.defaultdict(list)
    childParentComboDict = {}
    logDict = collections.defaultdict(lambda: collections.defaultdict(list))
    condensedLogDict = collections.defaultdict(lambda: collections.defaultdict(list))
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
                        try:
                            logDict[childName][parentName].append(f'non-matching tax IDs ({list(taxIdChildSet)[0]}; {list(taxIdParentSet)[0]})')
                        except IndexError:
                            logDict[childName][parentName].append(f'non-matching tax IDs (NA; NA)')
                        #condensedLogDict[childName]['non-matching tax IDs'].append(parentName)
                        continue

                    # 3. check BLAST quality on percent identity and query coverage are lower for child than parent
                    # only check top BLAST hit for now. Probably not accurate, so will need to be altered in future
                    logDict[childName][parentName].append(f'matching tax IDs ({list(taxIdParentSet.intersection(taxIdChildSet))[0]})')
                    condensedLogDict[childName]['matching tax IDs'].append(parentName)
                    if taxPidentInputDict[childName][0] > taxPidentInputDict[parentName][0] and taxQcovInputDict[childName][0] > taxQcovInputDict[parentName][0]:
                        logDict[childName][parentName].append(f'BLAST score threshold not met ({taxPidentInputDict[childName][0]}, {taxPidentInputDict[parentName][0]}; {taxQcovInputDict[childName][0]}, {taxQcovInputDict[parentName][0]})')
                        #condensedLogDict[childName]['BLAST score threshold not met'].append(parentName)
                        continue

                    # 4. check co-occurrence pattern
                    # 4.1 check if child only appears in samples where parent is present
                    logDict[childName][parentName].append(f'BLAST score threshold met ({taxPidentInputDict[childName][0]}, {taxPidentInputDict[parentName][0]}; {taxQcovInputDict[childName][0]}, {taxQcovInputDict[parentName][0]})')
                    condensedLogDict[childName]['BLAST score threshold met'].append(parentName)
                    if occurrence_type_ == 'presence-absence':
                        positiveDetectionsChild = [k for k, v in freqInputDictSubset[childName].items() if v >= int(detection_threshold_)]
                        positiveDetectionsParent = [k for k, v in freqInputDictSubset[parentName].items() if v >= int(detection_threshold_)]
                        missingCount = 0
                        for item in positiveDetectionsChild:
                            if item not in positiveDetectionsParent:
                                missingCount += 1
                        totalCount = len(positiveDetectionsParent) + missingCount
                        totalRatio = 1 - (missingCount / totalCount)
                        if totalRatio < ratio:
                            logDict[childName][parentName].append(f'co-occurrence ratio not met ({float("{:.2f}".format(totalRatio))}%)')
                            #condensedLogDict[childName]['co-occurrence ratio not met'].append(parentName)
                            continue

                    # 4.2 check if child only has lower abundance in samples compared to parent
                    elif occurrence_type_ == 'abundance':
                        count = 0
                        totalCount = 0
                        for item in freqInputDictSubset[parentName]:
                            parentValue = freqInputDictSubset[parentName][item]
                            childValue = freqInputDictSubset[childName][item]
                            if parentValue < int(detection_threshold_):
                                parentValue = 0
                            if childValue < int(detection_threshold_):
                                childValue = 0
                            if parentValue < childValue:
                                count += 1
                            if parentValue > 0:
                                totalCount += 1
                            if childValue > 0 and parentValue == 0:
                                totalCount += 1
                        totalRatio = 1 - (count / totalCount)
                        if totalRatio < ratio:
                            logDict[childName][parentName].append(f'co-occurrence ratio not met ({float("{:.2f}".format(totalRatio))}%)')
                            #condensedLogDict[childName]['co-occurrence ratio not met'].append(parentName)
                            continue
                    else:
                        console.print(f"[cyan]\n|               ERROR[/] | [bold yellow]'--occurrence-type' not specified as 'presence-absence' or 'abundance', aborting analysis...[/]\n")
                        exit()

                    # 5. check sequence similarity
                    logDict[childName][parentName].append(f'co-occurrence ratio met ({float("{:.2f}".format(totalRatio))}%)')
                    condensedLogDict[childName]['co-occurrence rate met'].append(parentName)
                    alignment = pairwise2.align.globalxx(seqInputDict[parentName], seqInputDict[childName])
                    distanceCalculation = sum(1 for a, b in zip(alignment[0][0], alignment[0][1]) if a != b)
                    if 100 - (distanceCalculation/ max(len(seqInputDict[parentName]), len(seqInputDict[childName])) * 100) <= int(similarity):
                        logDict[childName][parentName].append(f'sequence similarity threshold not met ({float("{:.2f}".format(100 - (distanceCalculation/ max(len(seqInputDict[parentName]), len(seqInputDict[childName])) * 100)))}%)')
                        #condensedLogDict[childName]['sequence similarity threshold not met'].append(parentName)
                        continue

                    # if it passes all the checks, we need to determine how it can be combined --> several options
                    # first: if parent not identified as a child previously, we can combine child and parent data
                    logDict[childName][parentName].append(f'sequence similarity threshold met ({float("{:.2f}".format(100 - (distanceCalculation/ max(len(seqInputDict[parentName]), len(seqInputDict[childName])) * 100)))}%)')
                    condensedLogDict[childName]['sequence similarity threshold met'].append(parentName)
                    if parentName not in childParentComboDict:
                        childParentComboDict[childName] = parentName
                        logDict[childName][parentName].append(f'parent identified!')
                        condensedLogDict[childName]['parent identified!'].append(parentName)
                        combinedDict[parentName].append(childName)
                        for item in newParentDict:
                            newValue = int(newParentDict[item]) + int(freqInputDict[childName][item])
                            newParentDict[item] = newValue

                    # second: childName not already identified as a child previously and parent identified as a child previously, add childName data to parent of parentName data
                    elif parentName in childParentComboDict:
                        combinedDict[childParentComboDict[parentName]].append(childName)
                        childParentComboDict[childName] = childParentComboDict[parentName]
                        logDict[childName][parentName].append(f'grandparent identified ({childParentComboDict[parentName]})!')
                        condensedLogDict[childName]['grandparent identified!'].append(childParentComboDict[parentName])
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
    
    ## write detailed log file
    try:
        with open(detailed_log_, 'w') as logOutFile:
            logOutFile.write('#################\n#### SUMMARY ####\n#################\n\n')
            logOutFile.write(f'date-time: {formattedTime}\n\n')
            logOutFile.write(f'parameters:\n')
            logOutFile.write(f'--method: taxon-dependent co-occurrence (default)\n')
            logOutFile.write(f'--occurrence type: {occurrence_type_}\n')
            logOutFile.write(f'--detection threshold: {detection_threshold_}\n')
            logOutFile.write(f'--similarity threshold: {similarity}\n')
            logOutFile.write(f'--co-occurrence ratio: {ratio}\n')
            logOutFile.write(f'--sample exclusion list: {", ".join(fullNegativeList)}\n\n')
            logOutFile.write(f'results:\n')
            logOutFile.write(f'--total seqs: {len(seqInputDict)}\n')
            logOutFile.write(f'--total artefacts: {len(childParentComboDict)} ({float("{:.2f}".format(len(childParentComboDict) / len(seqInputDict) * 100))}%)\n')
            for item in combinedDict:
                logOutFile.write(f'--parent {item}: {", ".join(combinedDict[item])}\n')
            logOutFile.write(f'\ncode: tombRaider {commandLineInput}\n\n\n')
            logOutFile.write('###########################\n#### DETAILED ANALYSIS ####\n###########################\n\n')
            for item in logDict:
                logOutFile.write(f'### analysing: {item} ###\n')
                for subitem in logDict[item]:
                    logOutFile.write(f'{subitem}:\t')
                    outputString = "\t".join(logDict[item][subitem])
                    logOutFile.write(f'{outputString}\n')
                logOutFile.write('\n')
    except TypeError:
        pass

    ## write condensed log file
    try:
        with open(condensed_log_, 'w') as logOut:
            logOut.write('#################\n#### SUMMARY ####\n#################\n\n')
            logOut.write(f'date-time: {formattedTime}\n\n')
            logOut.write(f'parameters:\n')
            logOut.write(f'--method: taxon-dependent co-occurrence (default)\n')
            logOut.write(f'--occurrence type: {occurrence_type_}\n')
            logOut.write(f'--detection threshold: {detection_threshold_}\n')
            logOut.write(f'--similarity threshold: {similarity}\n')
            logOut.write(f'--co-occurrence ratio: {ratio}\n')
            logOut.write(f'--sample exclusion list: {", ".join(fullNegativeList)}\n\n')
            logOut.write(f'results:\n')
            logOut.write(f'--total seqs: {len(seqInputDict)}\n')
            logOut.write(f'--total artefacts: {len(childParentComboDict)} ({float("{:.2f}".format(len(childParentComboDict) / len(seqInputDict) * 100))}%)\n')
            for item in combinedDict:
                logOut.write(f'--parent {item}: {", ".join(combinedDict[item])}\n')
            logOut.write(f'\ncode: tombRaider {commandLineInput}\n\n\n')
            logOut.write('############################\n#### CONDENSED ANALYSIS ####\n############################\n\n')
            for seq in freqTotalCountDict:
                if seq in condensedLogDict:
                    logOut.write(f'### analysing: {seq} ###\n')
                    for subitem in condensedLogDict[seq]:
                        logOut.write(f'{subitem}:\t')
                        outputString = ", ".join(condensedLogDict[seq][subitem])
                        logOut.write(f'{outputString}\n')
                    logOut.write('\n')
    except TypeError:
        pass

####################
# LULU ALTERNATIVE #
####################
def taxonIndependentCoOccurrenceAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, condensed_log_, detailed_log_, occurrence_type_, detection_threshold_, similarity, negative, ratio, seqname, taxid, pident, qcov, eval_):
    '''
    The function to identify and merge parent-child sequences based on the LULU algorithm (taxon-independent co-occurrence patterns)
    '''
    console = rich.console.Console(stderr=True, highlight=False)

####################
# MRCA ALTERNATIVE #
####################
def taxonDependentMergingAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, condensed_log_, detailed_log_, occurrence_type_, detection_threshold_, similarity, negative, ratio, seqname, taxid, pident, qcov, eval_):
    '''
    Thefunction to identify and merge parent-child sequences based on the taxonomic ID of sequences
    '''
    console = rich.console.Console(stderr=True, highlight=False)

###################
# MRCA CALCULATOR #
###################    
def mrcaCalculatorAlgorithm(frequency_input_, sequence_input_, taxonomy_input_, frequency_output_, sequence_output_, taxonomy_output_, condensed_log_, detailed_log_, occurrence_type_, detection_threshold_, similarity, negative, ratio, seqname, taxid, pident, qcov, eval_):
    '''
    The function to calculate the Most Recent Common Ancestor from a standard BLAST output
    '''
    console = rich.console.Console(stderr=True, highlight=False)