#!/usr/bin/python
from __future__ import print_function
import sys
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import re
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import statsmodels.stats.multitest as smm

pd.set_option('display.max_columns', 500)
    
def main(argv):
    timestamp = datetime.now().strftime('%Y%m%H%M')
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metafile", dest="metafile", help="Metadata file for BarSeq runs. A tab-delimited file with columns titled Sample, Group, Reference, OutputDir, PoolCountFile.  Entires in the Sample and Reference columns must correspond to column heading in the poolCountFile. Samples from biological replicates ",default="metafile.txt")
    parser.add_argument("-l", "--logFile", dest="logFile", help="File to write run log to. Default is Fitness_TIMESTAMP.log",default="Fitness_"+timestamp+".log")
    parser.add_argument("-L", "--normLocal", dest="normLocal", help="-L [int]: Normalize fitness scores such that local windows of [int] contiguous insertions have a median fitness score of zero.  Default behavior is to normalize fitness scores accross full contigs such that the median fitness score is zero", default=0, type=int)
    parser.add_argument("-s", "--fitnessStartBoundary", dest="fitnessStartBoundary", help="Integer 0-100. Fraction of the region between coding start and coding stop to set as start boundary for calculating gene fitness.  Default is 5: i.e. insertions in the first 5 percent of the region between start and stop will not be included in calculations for gene fitness scores or statistical tests", default=5, type=int)
    parser.add_argument("-e", "--fitnessEndBoundary", dest="fitnessEndBoundary", help="Integer 0-100. Fraction of the region between coding start and coding stop to set as end boundary for calculating gene fitness.  Default is 95: i.e. insertions in the last 5 percent of the region between start and stop will not be included in calculations for gene fitness scores or statistical tests", default=95, type=int)
    parser.add_argument("-i", "--minInsertionCounts", dest="minInsertionCounts", help="Integer. Minimum counts each insertion must have between both the test condition and the reference condition to be included in calculations for gene fitness scores or statistical tests. Default is 3", default=3, type=int)
    parser.add_argument("-g", "--minGeneCounts", dest="minGeneCounts", help="Integer. Minimum total counts across all insertions in a given gene before gene fitness scores are calculated. Default is 20", default=20, type=int)
    parser.add_argument("-I", "--minInsertions", dest="minInsertions", help="Integer. Minimum number of insertions with sufficient counts in a given gene before gene fitness scores are calculated. Default is 3", default=3, type=int)
    parser.add_argument("-W", "--maxWeightCounts", dest="maxWeightCounts", help="Integer. Maximum number of counts to be used for gene fitness score calculation.  E.g. if set at 50 then two insertions with 50 and 100 counts would be weighted equally in computing gene fitness, but an insertion with 10 counts would have a smaller weight ", default=50, type=int)
    parser.add_argument("-P", "--noPseudoCounts", dest="smartPseudoCounts", action='store_false', help="If this flag is passed, fitness scores will NOT be computed with 'smart' pseudocounts as in Wetmore et al 2015.", default=True)
    parser.add_argument("-B", "--fitnessBrowserOutput", dest="fitnessBrowserOutput", action='store_true', help="If this flag is passed, output files will be saved in formats compatible with the Arkin Lab Fitness Browser", default=False)
    parser.add_argument("-C", "--centerOnMedian", dest="centerOnMedian", action='store_true', help="By default strain fitness scores are noramlized to a mean of zero in a given sample.  If this flag is passed they will be normalized to the median", default=True)

    options = parser.parse_args()


    statusUpdate = 'RBseq_Calculate_Fitness.py  Samuel Coradetti 2019.'
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = 'Version 1.0.7'
    printUpdate(options.logFile,statusUpdate)

    optionDict = options.__dict__
    statusUpdate = 'Options passed: '
    for option in optionDict:
        statusUpdate+=" "+option+":"+str(optionDict[option])+" "
    printUpdate(options.logFile,statusUpdate)
    
    statusUpdate = 'Logging status updates in '+options.logFile
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = 'Loading BarSeq library metadata from '+options.metafile
    printUpdate(options.logFile,statusUpdate)
     
    #Load up experiment metadata
    fileToOpen = options.metafile
    if (not fileToOpen):
        statusUpdate = options.metafile + ' could not be loaded... exiting.'
        printUpdate(options.logFile,statusUpdate)
        sys.exit()
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            metaFrame = pd.read_csv(FileHandle,low_memory=False,sep='\t')
            metaFrame = metaFrame[~metaFrame[metaFrame.columns[0]].isnull()]
            FileHandle.close()
    except IOError:
        statusUpdate = "Could not read file:"+fileToOpen+" ...exiting."
        printUpdate(options.logFile,statusUpdate)
        sys.exit()

    requiredColumns=['Sample','Reference','OutputDir','Condition','PoolCountFile','Paired']
    for requiredColumn in requiredColumns:
        if requiredColumn not in metaFrame.columns:
            statusUpdate = "Metadata file is missing column "+requiredColumn+" ...exiting."
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

    if options.fitnessBrowserOutput:
        requiredColumns=['SetName','Group','Index','short','Condition_1','Date']
        for requiredColumn in requiredColumns:
            if requiredColumn not in metaFrame.columns:
                statusUpdate = "-B option passed but metadata file is missing column "+requiredColumn+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()
            

    #Get a list of fastqs and matching model files from the metadata
    sampleNames=metaFrame['Sample'].values
    sampleNumbers = list(range(0,len(sampleNames)))
    references=metaFrame['Reference'].values
    outputDirs=metaFrame['OutputDir'].values
    groupNames=metaFrame['Condition'].values
    poolCountFiles=metaFrame['PoolCountFile'].values
    pairedFlags=metaFrame['Paired'].values
    groupSet = []

    groupSampNumbers = {}
    for sampNum,groupName in enumerate(groupNames):
        if groupName in groupSampNumbers:
            groupSampNumbers[groupName].append(sampNum)
        else:
            groupSampNumbers[groupName] = [sampNum]
            groupSet.append(groupName)
        

    for outputN,dir in enumerate(outputDirs):
        if not dir[-1] == "/":
            outputDirs[outputN] = dir+"/"
        if not os.path.exists(outputDirs[outputN]):
            os.makedirs(outputDirs[outputN])

    if not os.path.exists(outputDirs[outputN]+'QCplots/'):
            os.makedirs(outputDirs[outputN]+'QCplots/')

    numSample = len(sampleNames)
    numGroups = len(set(groupNames))
    
    statusUpdate = 'Loaded '+str(numSample)+' samples from ' +str(numGroups)+ ' groups of biological replicates from '+ options.metafile
    printUpdate(options.logFile,statusUpdate)


    #Load and parse poolCount file
    fileToOpen = poolCountFiles[0]
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            poolCounts = pd.read_csv(fileToOpen, low_memory=False, dtype={'NearestGene':str},sep='\t')
            FileHandle.close()
    except IOError:
        statusUpdate = " Could not read file:"+fileToOpen+" ...exiting."
        printUpdate(options.logFile,statusUpdate)
        sys.exit()

    #extract approximate gene positions from poolcount file
    geneInfo = poolCounts[['NearestGene','scaffold','pos']]
    genePositions = geneInfo.groupby('NearestGene')['scaffold','pos'].min()


    #Make new averaged counts columns for sets of reference conditions used in non-paired comparisons.
    for replicateGroup in groupSet:
        if False in pairedFlags[groupSampNumbers[replicateGroup]]:
            newColumn = replicateGroup + '_averaged_references'
            poolCounts[newColumn] = poolCounts[references[groupSampNumbers[replicateGroup]]].mean(axis=1)
            for sampNum in groupSampNumbers[replicateGroup]:
                references[sampNum] = newColumn
    

    statusUpdate = 'Loaded counts for ' + str(len(poolCounts)-1) + ' barcodes from ' + poolCountFiles[0]
    printUpdate(options.logFile,statusUpdate)

    statusUpdate = 'Finding barcodes in coding regions with sufficient counts for fitness analysis'
    printUpdate(options.logFile,statusUpdate)

    infoColumns = poolCounts.columns.values[1:5]
    dataColumns = poolCounts.columns.values[5:]
    barcodesDict = {'barcode':poolCounts['barcode'].values}

    totalCounts = pd.DataFrame.from_dict(barcodesDict)
    strainFit = pd.DataFrame.from_dict(barcodesDict)
    sampleCounts = pd.DataFrame.from_dict(barcodesDict)
    refCounts = pd.DataFrame.from_dict(barcodesDict)
    pseudoCounts = pd.DataFrame.from_dict(barcodesDict)
    
    
    for column in infoColumns:
        totalCounts[column] = poolCounts[column].values
        strainFit[column] = poolCounts[column].values
        sampleCounts[column] = poolCounts[column].values
        refCounts[column] = poolCounts[column].values
        pseudoCounts[column] = poolCounts[column].values
    
    #identify insertions within acceptable range between start and stop codons
    eligibleCounts = totalCounts.copy()
    eligibleCounts1 = totalCounts.copy()
    eligibleCounts2 = totalCounts.copy()
    centralCodingMask = np.ma.masked_inside(poolCounts['CodingFraction'].values, options.fitnessStartBoundary,options.fitnessEndBoundary).mask
    firstHalfMask = np.ma.masked_inside(poolCounts['CodingFraction'].values, 2,49).mask
    secondHalfMask = np.ma.masked_inside(poolCounts['CodingFraction'].values, 50,98).mask

    weightingCounts = totalCounts.copy()
    weightingCounts1 = totalCounts.copy()
    weightingCounts2 = totalCounts.copy()
    includeMask = totalCounts.copy()
    includeMask1 = totalCounts.copy()
    includeMask2 = totalCounts.copy()
    weightingFactors = totalCounts.copy()
    weightingFactors1 = totalCounts.copy()
    weightingFactors2 = totalCounts.copy()

    
    for sampNum,sampleName in enumerate(sampleNames):
        #total up all counts between experimental and reference conditions for each comparison
        totalCounts[sampNum] = poolCounts[sampleNames[sampNum]] + poolCounts[references[sampNum]]

        #mask out insertions with low counts in each experiment
        eligibleCounts[sampNum] = np.ma.masked_inside(totalCounts[sampNum],0,options.minInsertionCounts).filled(0)

        #mask out insertions outside acceptable range in coding regions
        eligibleCounts[sampNum] = eligibleCounts[sampNum] * centralCodingMask.astype('int')
        eligibleCounts1[sampNum] = eligibleCounts[sampNum] * firstHalfMask.astype('int')
        eligibleCounts2[sampNum] = eligibleCounts[sampNum] * secondHalfMask.astype('int')
        
        includeMask[sampNum] = eligibleCounts[sampNum].astype('bool')
        includeMask1[sampNum] = eligibleCounts1[sampNum].astype('bool')
        includeMask2[sampNum] = eligibleCounts2[sampNum].astype('bool')

        #calculate weighting factors for each strain fitness score, will need to mask out unused strains later
        weightingFactors[sampNum] = (np.log(2)**2) / (1.0 / (1 + poolCounts[sampleNames[sampNum]]) + 1.0 / (1 + poolCounts[references[sampNum]]))

    #maximum weights for any strain is that of a strain with maximumWeightCounts (e.g. 20) in both sample and reference condition
    maxWeight = (np.log(2)**2) / (1.0 / (1 + options.maxWeightCounts) + 1.0 / (1 + options.maxWeightCounts))
    weightingFactors[sampleNumbers] = np.ma.masked_outside(weightingFactors[sampleNumbers],0,maxWeight).filled(maxWeight)
    

    #mask out insertions in genes without sufficient total counts per gene
    #for comparing 1st and 2nd half of genes, both halves must meet minimum counts per gene requirement
    summedByGene = eligibleCounts.groupby('NearestGene')[sampleNumbers].transform(np.sum)
    summedByGene1 = eligibleCounts1.groupby('NearestGene')[sampleNumbers].transform(np.sum)
    summedByGene2 = eligibleCounts2.groupby('NearestGene')[sampleNumbers].transform(np.sum)
    summedByGeneHalves = pd.concat([summedByGene1, summedByGene2]).min(level=0)
    includeMask[sampleNumbers] = includeMask[sampleNumbers] & np.ma.masked_outside(summedByGene,0,options.minGeneCounts).mask
    includeMask1[sampleNumbers] = includeMask1[sampleNumbers] & np.ma.masked_outside(summedByGeneHalves,0,options.minGeneCounts).mask
    includeMask2[sampleNumbers] = includeMask2[sampleNumbers] & np.ma.masked_outside(summedByGeneHalves,0,options.minGeneCounts).mask

    #if inserts meet minimum criteria in one replicate of a set, compute fitness across all replicates
    for group,sampList in groupSampNumbers.items(): 
        for sampNum in sampList:
             includeMask[sampNum] = (includeMask[sampList].sum(axis=1)>0)
             includeMask1[sampNum] = (includeMask1[sampList].sum(axis=1)>0)
             includeMask2[sampNum] = (includeMask2[sampList].sum(axis=1)>0)

    #mask out insertions for genes without enough insertions that meet other critera for fitness calculations
    #for comparing 1st and 2nd half of genes, both halves must meet minimum insertions per gene requirement
    enoughInserts = np.ma.masked_outside(includeMask.groupby('NearestGene')[sampleNumbers].transform(np.sum),0,options.minInsertions-1).mask
    enoughInserts1 = np.ma.masked_outside(includeMask1.groupby('NearestGene')[sampleNumbers].transform(np.sum),0,options.minInsertions-1).mask
    enoughInserts2 = np.ma.masked_outside(includeMask2.groupby('NearestGene')[sampleNumbers].transform(np.sum),0,options.minInsertions-1).mask
    enoughInsertsHalves = enoughInserts1 & enoughInserts2
    
    includeMask[sampleNumbers] = includeMask[sampleNumbers] & enoughInserts
    includeMask1[sampleNumbers] = includeMask1[sampleNumbers] & enoughInsertsHalves
    includeMask2[sampleNumbers] = includeMask2[sampleNumbers] & enoughInsertsHalves
    
    statusUpdate = 'Genes with sufficient counts for fitness analysis in each experiment:'
    printUpdate(options.logFile,statusUpdate)
    genesPassed = includeMask.groupby('NearestGene')[sampleNumbers].max().astype('int').sum()
    for experiment,sampNums in groupSampNumbers.items():
        statusUpdate = '  ' + experiment + ': ' + str(genesPassed[sampNums[0]])
        printUpdate(options.logFile,statusUpdate)

    statusUpdate = 'Calculating strain fitness'
    printUpdate(options.logFile,statusUpdate)
    
    #Normalize counts accross samples
    analyzedSamples = list(sampleNames)
    analyzedSamples.extend(list(references))
    analyzedSamples = list(set(analyzedSamples))

    medianCounts = poolCounts[analyzedSamples].sum().median()
    normFactors = 100000000/poolCounts[analyzedSamples].sum()
    poolCountsNotNormalized = poolCounts.copy()
    poolCounts[analyzedSamples] = poolCounts[analyzedSamples] * normFactors
    
    #calculate intial strain fitness. 
    comparisons = []
    for sampNum,sampleName in enumerate(sampleNames):
        comparisonText = sampleName + "_vs_" + references[sampNum]
        comparisons.append(comparisonText)
        #Replace zero counts with 0.1 to prevent infinite results
        strainFit[sampNum] = np.log2((poolCounts[sampleName] + 0.1) / (poolCounts[references[sampNum]] + 0.1))
        if options.centerOnMedian:
            strainFit[sampNum] = strainFit[sampNum] - strainFit[sampNum].median()
        else:
            strainFit[sampNum] = strainFit[sampNum] - strainFit[sampNum].mean()


    statusUpdate = 'Calculating gene fitness'
    printUpdate(options.logFile,statusUpdate)

    #filter weighting factors based on includeMask
    weightingFactors1[sampleNumbers] = weightingFactors[sampleNumbers] * includeMask1[sampleNumbers].astype('int')
    weightingFactors2[sampleNumbers] = weightingFactors[sampleNumbers] * includeMask2[sampleNumbers].astype('int')
    weightingFactors[sampleNumbers] = weightingFactors[sampleNumbers] * includeMask[sampleNumbers].astype('int')
    
    weights = totalCounts.copy()
    weights1 = totalCounts.copy()
    weights2 = totalCounts.copy()
    weightedFit = totalCounts.copy()
    weightedFit1 = totalCounts.copy()
    weightedFit2 = totalCounts.copy()
    geneGroups = weightingFactors.groupby('NearestGene')
    geneGroups1 = weightingFactors1.groupby('NearestGene')
    geneGroups2 = weightingFactors2.groupby('NearestGene')
    weights[sampleNumbers] = weightingFactors[sampleNumbers]/geneGroups[sampleNumbers].transform(np.sum)
    weights1[sampleNumbers] = weightingFactors1[sampleNumbers]/geneGroups1[sampleNumbers].transform(np.sum)
    weights2[sampleNumbers] = weightingFactors2[sampleNumbers]/geneGroups2[sampleNumbers].transform(np.sum)
    weightedFit[sampleNumbers] = strainFit[sampleNumbers] * weights[sampleNumbers]
    weightedFit1[sampleNumbers] = strainFit[sampleNumbers] * weights1[sampleNumbers]
    weightedFit2[sampleNumbers] = strainFit[sampleNumbers] * weights2[sampleNumbers]
    geneFit = weightedFit.groupby('NearestGene')[sampleNumbers].sum()
    geneFit1 = weightedFit1.groupby('NearestGene')[sampleNumbers].sum()
    geneFit2 = weightedFit2.groupby('NearestGene')[sampleNumbers].sum()

    #Compute 'smart' pseudocounts a la Wetmore 2015.
    #Need sample counts and reference counts for later even if not using pseudocounts
    for sampNum,sampleName in enumerate(sampleNames):
            sampleCounts[sampNum] = includeMask[sampNum].astype('int').values*poolCounts[sampleName]
            refCounts[sampNum] = includeMask[sampNum].astype('int').values*poolCounts[references[sampNum]]
    #only need this code if using pseudocounts
    if options.smartPseudoCounts:
        statusUpdate = 'Recomputing fitness scores with smart pseudocounts'
        printUpdate(options.logFile,statusUpdate)
        readRatios = (sampleCounts.groupby('NearestGene')[sampleNumbers].sum() + 0.1)/(refCounts.groupby('NearestGene')[sampleNumbers].sum() + 0.1)
        pseudoCountsByGene = (2**geneFit)*readRatios
        pseudoCounts = pseudoCounts.merge(pseudoCountsByGene, on='NearestGene').fillna(0.1)
        pseudoCounts[sampleNumbers] = pseudoCounts[sampleNumbers].replace(0.0, 0.1)
        
        for sampNum,sampleName in enumerate(sampleNames):
            strainFit[sampNum] = np.log2((poolCounts[sampleName] + pseudoCounts[sampNum]**0.5) / (poolCounts[references[sampNum]] + pseudoCounts[sampNum]**-0.5))
            if options.centerOnMedian:
                strainFit[sampNum] = strainFit[sampNum] - strainFit[sampNum].median()
            else:
                strainFit[sampNum] = strainFit[sampNum] - strainFit[sampNum].mean()

    #Option to normalize strain fitness over local neighborhood
    if options.normLocal > 0:
        statusUpdate = 'Normalizing over rolling window'
        printUpdate(options.logFile,statusUpdate)
        window = options.normLocal
        for sampNum in sampleNumbers:
            rollingMean = []
            for scaffold,scaffoldGroup in strainFit.groupby('scaffold'):
                if options.centerOnMedian:
                    rollingMean.extend(scaffoldGroup[sampNum].rolling(window,center=True).median().fillna(method='ffill').fillna(method='bfill').values)
                else:
                    rollingMean.extend(scaffoldGroup[sampNum].rolling(window,center=True).mean().fillna(method='ffill').fillna(method='bfill').values)
            strainFit[sampNum] = strainFit[sampNum] - rollingMean

    
    #Recalculate gene fitness based on adjusted strain fitness        
    weightedFit[sampleNumbers] = strainFit[sampleNumbers] * weights[sampleNumbers]
    weightedFit1[sampleNumbers] = strainFit[sampleNumbers] * weights1[sampleNumbers]
    weightedFit2[sampleNumbers] = strainFit[sampleNumbers] * weights2[sampleNumbers]
    geneFit = weightedFit.groupby('NearestGene')[sampleNumbers].sum()
    geneFit1 = weightedFit1.groupby('NearestGene')[sampleNumbers].sum()
    geneFit2 = weightedFit2.groupby('NearestGene')[sampleNumbers].sum()

    #filter
    geneWeights = weights.groupby('NearestGene')[sampleNumbers].sum()
    geneWeightNonZero = (geneWeights[sampleNumbers] != 0).astype('int').replace(0.0, np.nan)
    geneFit[sampleNumbers] = geneFit[sampleNumbers]*geneWeightNonZero

    uniqGeneAnnotations = poolCounts[['NearestGene','AlternateID','Annotation']].copy().drop_duplicates()
    uniqGeneAnnotations = uniqGeneAnnotations.set_index('NearestGene')

    geneFitAnnotated = geneFit.merge(uniqGeneAnnotations, left_index=True, right_on='NearestGene',how='left')
    geneFitAnnotated = geneFitAnnotated[['AlternateID','Annotation']+sampleNumbers]
        
    #Save fitness data
    headerText = list(strainFit.columns[:-len(sampleNumbers)].values)
    headerText.extend(comparisons)
    statusUpdate = 'Saving all strain fitness scores in ' + outputDirs[0]+'strainFit.txt'
    printUpdate(options.logFile,statusUpdate)
    strainFit.to_csv(outputDirs[0]+'strainFit.txt',sep="\t",index=False,header=headerText)

    strainFitUsed = strainFit.copy()
    strainFitUsed[sampleNumbers] = strainFit[sampleNumbers] * includeMask[sampleNumbers].astype('int')

    statusUpdate = 'Saving strain fitness scores used in fitness calculations in ' + outputDirs[0]+'strainFitUsed.txt'
    printUpdate(options.logFile,statusUpdate)
    strainFitUsed.to_csv(outputDirs[0]+'strainFitUsed.txt',sep="\t",index=False,header=headerText)
    
    statusUpdate = 'Saving gene fitness scores in ' + outputDirs[0]+'geneFit_individual_replicates.txt'
    printUpdate(options.logFile,statusUpdate)
    geneFitAnnotated.to_csv(outputDirs[0]+'geneFit_individual_replicates.txt',sep="\t",header=['AlternateID','Annotation']+comparisons)

    #Aggregate accross replicates
    geneFitAverages = geneFitAnnotated[['AlternateID','Annotation']].copy()
    geneFitAverages1 = pd.DataFrame(index=geneFit.index)
    geneFitAverages2 = pd.DataFrame(index=geneFit.index)

    
    for replicateGroup in groupSet:
        geneFitAverages[replicateGroup] = geneFit[groupSampNumbers[replicateGroup]].mean(axis=1)
        geneFitAverages1[replicateGroup] = geneFit1[groupSampNumbers[replicateGroup]].mean(axis=1)
        geneFitAverages2[replicateGroup] = geneFit2[groupSampNumbers[replicateGroup]].mean(axis=1)

     
    statusUpdate = 'Saving average gene fitness scores for conditions with replicates in ' + outputDirs[0]+'geneFit.txt'
    printUpdate(options.logFile,statusUpdate)
    geneFitAverages.to_csv(outputDirs[0]+'geneFit.txt',sep="\t")
    
    #Calculate T-like statistics as in wetmore et al 2015
    #Naive variance per gene Vn based on poisson noise
    Nafter = sampleCounts.groupby('NearestGene')[sampleNumbers].sum()
    Nbefore = refCounts.groupby('NearestGene')[sampleNumbers].sum()
    
    NafterAverages = pd.DataFrame(index=Nafter.index)
    NbeforeAverages = pd.DataFrame(index=Nbefore.index)
    for replicateGroup in groupSet:
        NafterAverages[replicateGroup] = Nafter[groupSampNumbers[replicateGroup]].sum(axis=1)
        NbeforeAverages[replicateGroup] = Nbefore[groupSampNumbers[replicateGroup]].sum(axis=1)
    Vn = (1/(1 + Nafter) + 1/(1 + Nbefore))/(np.log(2)**2)
    VnAverages = (1/(1 + NafterAverages) + 1/(1 + NbeforeAverages))/(np.log(2)**2)

    

    #MAD12 Mean absolute difference in fitness scores between first and second half of genes.
    MAD12 = (geneFit1[sampleNumbers] - geneFit2[sampleNumbers]).abs().mean()
    
    MAD12Averages = pd.Series()
    for replicateGroup in groupSet:
        MAD12Averages[replicateGroup] = MAD12[groupSampNumbers[replicateGroup]].mean()


    #Typical variance in each experiment
    Vt = MAD12**2/(2*0.672)**2
    VtAverages = MAD12Averages[groupSampNumbers.keys()]**2/(2*0.672)**2

    #a priori estimate of gene fitness
    Vg = Vn.copy()
    VgAverages = VnAverages.copy()
    for sampNum,sampleName in enumerate(sampleNames):
        Vg[sampNum] = Vt[sampNum] * Vn[sampNum]/Vn[sampNum].median()
    for replicateGroup in groupSet:
        VgAverages[replicateGroup] = VtAverages[replicateGroup] * VnAverages[replicateGroup]/VnAverages[replicateGroup].median()

    #Observed variance
    weightedResidules = weightedFit.copy()
    residules = strainFit[sampleNumbers] - weightedFit.groupby('NearestGene')[sampleNumbers].transform(np.sum)
    weightedResidules[sampleNumbers] = residules**2 * weights[sampleNumbers]
    geneResidules = weightedResidules.groupby('NearestGene')[sampleNumbers].sum()
    Nused = includeMask.groupby('NearestGene')[sampleNumbers].sum()
    Ve = (geneResidules + Vg)/Nused
        
    geneResidulesAverages = pd.DataFrame(index=geneResidules.index)
    NusedAverages = pd.DataFrame(index=geneResidules.index)
    NusedTotals = pd.DataFrame(index=geneResidules.index)
    VeAverages = VgAverages.copy
    for replicateGroup in groupSet:
        geneResidulesAverages[replicateGroup] = geneResidules[groupSampNumbers[replicateGroup]].mean(axis=1)
        NusedTotals[replicateGroup] = Nused[groupSampNumbers[replicateGroup]].sum(axis=1)
        NusedAverages[replicateGroup] = Nused[groupSampNumbers[replicateGroup]].mean(axis=1)
    VeAverages = (geneResidulesAverages + VgAverages)/NusedTotals

    #calculate T-stats
    maxVar = pd.concat([Vn, Ve]).max(level=0) + 0.01
    Tstats = geneFit.copy()
    Tstats[sampleNumbers] = geneFit[sampleNumbers]/(maxVar**(0.5))

    maxVarAverages = pd.concat([VnAverages, VeAverages]).max(level=0) + 0.01
    TstatsAverages = geneFitAverages.copy()
    TstatsAverages[groupSet] = geneFitAverages[groupSet]/(maxVarAverages**(0.5))

    statusUpdate = 'Saving T-like statistics (Wetmore et al 2015) in ' + outputDirs[0]+'Tstats.txt and '  + outputDirs[0]+'Tstats_individual_replicates.txt'
    printUpdate(options.logFile,statusUpdate)
    Tstats.to_csv(outputDirs[0]+'Tstats_individual_replicates.txt',sep="\t",header=comparisons)
    TstatsAverages.to_csv(outputDirs[0]+'Tstats.txt',sep="\t")

    pVals = TstatsAverages.copy()
    
    for replicateGroup in groupSet:
        analyzedGenes = pd.notnull(TstatsAverages[replicateGroup])
        analyzedGenes = analyzedGenes & (abs(TstatsAverages[replicateGroup]) > 0)
        result = stats.t.sf(TstatsAverages[analyzedGenes][replicateGroup].abs(), NusedTotals[analyzedGenes][replicateGroup]-1)*2
        adjusted = smm.multipletests(result, alpha=0.1, method='fdr_bh')[1]
        pVals[replicateGroup] = 1
        pVals.loc[pVals[analyzedGenes].index,replicateGroup] = adjusted

    statusUpdate = 'Saving pValues from T-like statistics in ' + outputDirs[0]+'pVals.txt'
    printUpdate(options.logFile,statusUpdate)
    pVals.to_csv(outputDirs[0]+'pVals.txt',sep="\t")

    statusUpdate = 'Generating quality control plots and saving in ' + outputDirs[0]+'QCplots/'
    printUpdate(options.logFile,statusUpdate)
    
    GCcor = []
    for replicateGroup in groupSet:
        #Generate plots of fitness in first and second half of the gene

        bothHalvesFilter = ((geneFitAverages1[replicateGroup] != 0) & (geneFitAverages2[replicateGroup] != 0))
        sigFilter01 = pVals[replicateGroup] < 0.05
        sigFilter001 = pVals[replicateGroup] < 0.001
        sig1 = geneFitAverages1[bothHalvesFilter & sigFilter01]
        sig2 = geneFitAverages2[bothHalvesFilter & sigFilter01]
        sig1b = geneFitAverages1[bothHalvesFilter & sigFilter001]
        sig2b = geneFitAverages2[bothHalvesFilter & sigFilter001]
        
        pp = PdfPages(outputDirs[0] + 'QCplots/' + replicateGroup + '_cor12.pdf')

        fig, ax = plt.subplots()
        fig.set_size_inches(6,6.6)
        ax.axhline(y=0, color='grey',linewidth=0.5,dashes=[1,1],zorder=-1)
        ax.axvline(x=0, color='grey',linewidth=0.5,dashes=[1,1],zorder=-1)
        
        plt.scatter(geneFitAverages1[replicateGroup],geneFitAverages2[replicateGroup],color='C0',alpha=0.5,s=2, edgecolor='none',label='All Genes')
        plt.scatter(sig1[replicateGroup],sig2[replicateGroup],color='C1',alpha=0.6,s=8, marker='o',facecolor='none',label='pVal < 0.05')
        plt.scatter(sig1b[replicateGroup],sig2b[replicateGroup],color='C2',alpha=0.8,s=20, marker='s',facecolor='none',label='pVal < 0.001')
        plt.xticks(fontsize=8)
        plt.xlabel('First Half of Gene',fontsize=8,labelpad=1)
        plt.yticks(fontsize=8)
        plt.ylabel('Seond Half of Gene',fontsize=8,labelpad=1)
        plt.title('Fitness Scores: First vs Second Half of Gene',fontsize=10)
        legend = plt.legend(loc='upper left',edgecolor=None,frameon=False, fontsize=8)
        plt.xticks([-6,-4,-2,0,2],fontsize=8)
        plt.yticks([-6,-4,-2,0,2],fontsize=8)
        plt.xlim(-8,3)
        plt.ylim(-8,3)
        
        plt.gcf().subplots_adjust(bottom=0.12, left=0.12, right=0.98, top=0.91)

        pp.savefig()
        plt.clf()
        plt.close()
        pp.close()

        #Generate QQ plot of T-statistics vs normal distribution
        pp = PdfPages(outputDirs[0] + 'QCplots/' + replicateGroup + '_qq.pdf')
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6.6)
        
        observed = np.sort(TstatsAverages[replicateGroup].values)
        observed = observed[~np.isnan(observed)]
        n = len(observed)
        intervals = np.arange(0.0,1.0,1.0/(n-0.00001))
        expected = stats.norm.ppf(intervals)
        expected = expected[:n]
        plt.plot(expected,expected,label='Normal Distribution',color='grey',dashes=[2,2])
        plt.plot(expected,observed,label=replicateGroup,linewidth=3,alpha=0.75)

        for sampNum in groupSampNumbers[replicateGroup]:
            observed = np.sort(Tstats[sampNum].values)
            observed = observed[~np.isnan(observed)]
            plt.plot(expected,observed,label=sampleNames[sampNum]+" vs "+references[sampNum],linewidth=1,alpha=0.75)
            

        plt.title('QQ plot for T-stats:' + replicateGroup,fontsize=10)
        legend = plt.legend(loc='upper left',edgecolor=None,frameon=False, fontsize=8,labelspacing=0.1,handletextpad=1)
        plt.xlim(-6,6)
        plt.ylim(-6,6)
        plt.xticks([-4,-2,0,2,4],fontsize=8)
        plt.yticks([-4,-2,0,2,4],fontsize=8)
        plt.xlabel('Standard Normal',fontsize=8,labelpad=1)
        plt.ylabel('T-statistic',fontsize=8,labelpad=1)
        plt.gcf().subplots_adjust(bottom=0.12, left=0.12, right=0.98, top=0.91)
        pp.savefig()
        plt.clf()
        plt.close()
        pp.close()

        #Generate plot of fitness vs position on largest chcromosome
        pp = PdfPages(outputDirs[0] + 'QCplots/' + replicateGroup + '_fit_v_position.pdf')
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6.6)

        longestScaffold = strainFit.loc[strainFit['pos'].idxmax(),'scaffold']
        scaffData = strainFit[strainFit['scaffold'] == longestScaffold]
        geneList = scaffData['NearestGene']
        posList = scaffData['pos']

        plotGenes = []
        plotFit = []
        plotPos = []
        for idx,gene in enumerate(geneList):
            if gene not in plotGenes:
                fitEntry = geneFitAverages.loc[gene,replicateGroup]
                if not np.isnan(fitEntry):
                    plotGenes.append(gene)
                    plotFit.append(geneFitAverages.loc[gene,replicateGroup])
                    plotPos.append(posList[idx])

        plt.scatter(plotPos,plotFit,s=2,alpha=0.5)

        window = 100
        cumSumFit = np.cumsum(plotFit, dtype=float)
        cumSumFit[window:] = cumSumFit[window:] - cumSumFit[:-window]
        runningAve = cumSumFit[window - 1:] / window
        plt.plot(plotPos[49:-50], runningAve, label='100 gene rolling average')
        
            

        plt.title('Fitness scores on largest scaffold:' + replicateGroup,fontsize=10)
        legend = plt.legend(loc='upper left',edgecolor=None,frameon=False, fontsize=8,labelspacing=0.1,handletextpad=1)
        plt.ylim(-2,2)
        plt.xticks(fontsize=8)
        plt.yticks([-2,-1,0,1,2],fontsize=8)
        plt.xlabel('Position on scaffold',fontsize=8,labelpad=1)
        plt.ylabel('Fitness',fontsize=8,labelpad=1)
        plt.gcf().subplots_adjust(bottom=0.12, left=0.12, right=0.98, top=0.91)
        pp.savefig()
        plt.clf()
        plt.close()
        pp.close()

        #Generate plots of fitness versus GC content around insertions

        pp = PdfPages(outputDirs[0] + 'QCplots/' + replicateGroup + '_GCcor.pdf')
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6.6)

        fitScores = strainFit[groupSampNumbers[replicateGroup]].mean(axis=1)

        GCcor_sample = np.corrcoef(fitScores,poolCounts['LocalGCpercent'])[0,1]
        GCcor.append(GCcor_sample)

        means = []
        meansPlus = []
        meansMinus = []
        GCpctRange = []
        for GCpct in list(range(35,75)):
            if len(fitScores[poolCounts['LocalGCpercent'] == GCpct]) > 50:
                GCpctRange.append(GCpct)
                mean = fitScores[poolCounts['LocalGCpercent'] == GCpct].mean(axis=0)
                upperQuartile = fitScores[poolCounts['LocalGCpercent'] == GCpct].quantile(0.75)
                lowerQuartile = fitScores[poolCounts['LocalGCpercent'] == GCpct].quantile(0.25)
                means.append(mean)
                meansPlus.append(upperQuartile)
                meansMinus.append(lowerQuartile)
        plt.plot(GCpctRange,means,color='C0',label='Average')
        plt.plot(GCpctRange,meansPlus,color='grey',dashes=[5,5],label='Upper Quartile')
        plt.plot(GCpctRange,meansMinus,color='grey',dashes=[3,3],label='Lower Quartile')
        midpoint = int(len(GCpctRange)/2)
        plt.annotate("R = {0:.2g}".format(GCcor_sample),(GCpctRange[midpoint],means[midpoint]-0.2),horizontalalignment='left', verticalalignment='bottom',fontsize=8)

        plt.title('Fitness Versus Local GC percent: ' + replicateGroup+'\n',fontsize=10)
        plt.xticks(fontsize=8)
        plt.yticks([-2,-1,0,1,2],fontsize=8)
        plt.xlabel('GC percent',fontsize=8,labelpad=1)
        plt.ylabel('Fitness',fontsize=8,labelpad=1)
        plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.98, top=0.90)
        plt.legend(bbox_to_anchor=(0.9, 1.045),ncol=3,fontsize=7,frameon=False)
        pp.savefig()
        plt.clf()
        plt.close()
        pp.close()

    genicCounts = totalCounts.copy()
    usedCounts = totalCounts.copy()
    usedCountsRefs = totalCounts.copy()
    for sampNum,sampleName in enumerate(sampleNames):
        genicCounts[sampNum] = poolCountsNotNormalized[sampleNames[sampNum]].values * centralCodingMask.astype('int')
        usedCounts[sampNum] = poolCountsNotNormalized[sampleNames[sampNum]].values * includeMask[sampNum].astype('int')
        usedCountsRefs[sampNum] = poolCountsNotNormalized[references[sampNum]].values * includeMask[sampNum].astype('int')
    
    Cor12Averages = pd.Series()
    MaxFit = pd.Series()
    groupTotals = pd.Series()
    groupGenic = pd.Series()
    groupUsed = pd.Series()
    groupUsedCounts = genicCounts[['barcode','NearestGene']].copy()
    groupUsedRefCounts = genicCounts[['barcode','NearestGene']].copy()
    groupMedian = pd.Series()
    groupMean = pd.Series()
    groupRefMedian = pd.Series()
    groupAdjCor = pd.Series()
    positionalFit = geneFitAverages.merge(genePositions, left_index=True, right_index=True)
    positionalFit = positionalFit.sort_values(['scaffold','pos'])
    for replicateGroup in groupSet:
        newFrame = pd.DataFrame(index=geneFit.index)
        newFrame['A'] = geneFitAverages1[replicateGroup]
        newFrame['B'] = geneFitAverages2[replicateGroup]
        Cor12Averages[replicateGroup] = newFrame.corr().iloc[0,1]
        MaxFit[replicateGroup] = strainFit[groupSampNumbers[replicateGroup]].mean(axis=1).max()
        groupTotals[replicateGroup] = totalCounts[groupSampNumbers[replicateGroup]].sum(axis=1).astype('bool').astype('int').sum()
        groupGenic[replicateGroup] = genicCounts[groupSampNumbers[replicateGroup]].sum(axis=1).astype('bool').astype('int').sum()
        groupUsed[replicateGroup] = usedCounts[groupSampNumbers[replicateGroup]].sum(axis=1).astype('bool').astype('int').sum()
        groupUsedCounts[replicateGroup] = usedCounts[groupSampNumbers[replicateGroup]].sum(axis=1)
        groupMean[replicateGroup] = groupUsedCounts.groupby('NearestGene')[replicateGroup].sum().replace(0, np.NaN).mean()
        groupMedian[replicateGroup] = groupUsedCounts.groupby('NearestGene')[replicateGroup].sum().replace(0, np.NaN).median()
        groupUsedRefCounts[replicateGroup] = usedCountsRefs[groupSampNumbers[replicateGroup]].sum(axis=1)
        groupRefMedian[replicateGroup] = groupUsedRefCounts.groupby('NearestGene')[replicateGroup].sum().replace(0, np.NaN).median()
        adjacentFit = positionalFit[replicateGroup].dropna()
        groupAdjCor[replicateGroup] = np.corrcoef(adjacentFit,np.roll(adjacentFit,1))[0,1]

    qualityStats = pd.DataFrame({'Condition':groupSet})
    qualityStats['nMapped'] = groupTotals.values
    qualityStats['nGenic'] = groupGenic.values
    qualityStats['nUsed'] = groupUsed.values
    qualityStats['gMed'] = groupMedian.values
    qualityStats['gMedt0'] = groupRefMedian.values
    qualityStats['gMean'] = groupMean.values
    qualityStats['nMapped'] = groupTotals.values
    qualityStats['cor12'] = Cor12Averages.values
    qualityStats['mad12'] = MAD12Averages.values
    qualityStats['gccor'] = GCcor
    qualityStats['adjcor'] = groupAdjCor.values
    qualityStats['maxFit'] = MaxFit.values
        

    statusUpdate = 'Saving quality statistics in ' + outputDirs[0]+'qualityStats.txt'
    printUpdate(options.logFile,statusUpdate)
    qualityStats.to_csv(outputDirs[0]+'qualityStats.txt',sep="\t",index=False)

    #find genes with specific phenotypes
    absFit = geneFitAverages.copy().fillna(value=0)
    absFit[groupSet] = absFit[groupSet].abs()
    percentileFit = absFit.copy()
    percentile95 = np.percentile(absFit[groupSet], 95, axis=1)

    for replicateGroup in groupSet:
        percentileFit[replicateGroup] = (percentileFit[replicateGroup] > (percentile95 + 0.5))
    
    magnitudeFit = absFit.copy()
    magnitudeFit[groupSet] = magnitudeFit[groupSet] > 1
    magnitudeT = TstatsAverages.copy()
    magnitudeT[groupSet] = magnitudeT[groupSet].abs()
    magnitudeT[groupSet] = magnitudeT[groupSet] > 3
    specificFit = geneFitAverages.copy()
    specificFit[groupSet] = percentileFit[groupSet] & magnitudeFit[groupSet] & magnitudeT[groupSet] 
    specificPass = specificFit[groupSet].max(axis=1).astype('bool')

    specificFitGenes = []
    specificFitConds = []
    for replicateGroup in groupSet:
        geneList = specificFit.index[specificFit[replicateGroup]].values
        specificFitConds.extend([replicateGroup]*len(geneList))
        specificFitGenes.extend(geneList)

    specificPhenotypes = pd.DataFrame({'NearestGene':specificFitGenes})
    specificPhenotypes = specificPhenotypes.merge(geneFitAverages[['AlternateID','Annotation']], right_index=True, left_on='NearestGene',how='left')
    specificPhenotypes['Condition'] = specificFitConds

    fitscores = []
    for idx,gene in enumerate(specificFitGenes):
        fitscores.append(geneFitAverages.loc[gene,specificFitConds[idx]])
    specificPhenotypes['Fitness'] = fitscores

    Tstats = []
    for idx,gene in enumerate(specificFitGenes):
        Tstats.append(TstatsAverages.loc[gene,specificFitConds[idx]])
    specificPhenotypes['Tstat'] = Tstats

    specificPhenotypes = specificPhenotypes.merge(metaFrame[['Condition','Experiment']].drop_duplicates(), left_on='Condition', right_on='Condition', how='left')

    statusUpdate = 'Saving specific phenotypes in ' + outputDirs[0]+'specificPhenotypes.txt'
    printUpdate(options.logFile,statusUpdate)
    specificPhenotypes.to_csv(outputDirs[0]+'specificPhenotypes.txt',sep="\t",index=False)

    
    #Create inputfiles for LBNL fitness browswer
    if options.fitnessBrowserOutput:
        if not os.path.exists(outputDirs[outputN]+'FitnessBrowser/'):
            os.makedirs(outputDirs[outputN]+'FitnessBrowser/')

        statusUpdate = 'Saving files for the fitness browser in:' + outputDirs[outputN]+'FitnessBrowser/'
        printUpdate(options.logFile,statusUpdate)

        #expsUsed
        firstReplicateNames = []
        for replicateGroup in groupSet:
            firstReplicateNames.append(sampleNames[groupSampNumbers[replicateGroup][0]])
        expsUsed = metaFrame[metaFrame['Sample'].isin(firstReplicateNames)].copy()
        expsUsed['name'] = expsUsed['SetName'].map(str) + expsUsed['Index']
        expsUsed['t0set'] = expsUsed['SetName'].map(str) + " " + expsUsed['Date']
        if 'Drop' not in expsUsed.columns:
            expsUsed['Drop'] = 'FALSE'
        if 'Description' not in expsUsed.columns:
            expsUsed['Description'] = expsUsed['Condition']
        if 'Date_pool_expt_started' not in expsUsed.columns:
            expsUsed['Date_pool_expt_started'] = expsUsed['Date']
        if 'Mutant.Library' not in expsUsed.columns:
            expsUsed['Mutant.Library'] = expsUsed['Mutant Library']
        if 'num' not in expsUsed.columns:
            expsUsed['num'] = list(range(1,len(firstReplicateNames)+1))

        expectedColumns = ['Description','SetName','Date_pool_expt_started','Person','Mutant.Library','Drop','gDNA.plate','gDNA.well',
                           'Index','Sequenced.At','Media','Growth.Method','Group','Temperature','pH','Liquid.v..solid','Aerobic_v_Anaerobic',
                           'Shaking','Condition_1','Concentration_1','Units_1','Condition_2','Concentration_2','Units_2','Timecourse','Timecourse.Sample',
                           'Growth.Plate.ID','Growth.Plate.wells','StartOD','EndOD','Total.Generations','num','name','short']

        for expectedColumn in expectedColumns:
            if expectedColumn not in expsUsed.columns:
                expsUsed[expectedColumn] = 'Not Passed'
        expsUsed[expectedColumns].to_csv(outputDirs[outputN]+'FitnessBrowser/expsUsed',sep="\t",index=False)

        #fit_quality.tab
        qualityStats = qualityStats.merge(expsUsed,left_on='Condition',right_on='Condition',how='left')
        qualityStats['mad12c'] = 'NA'
        qualityStats['mad12c_t0'] = 'NA'
        qualityStats['opcor'] = 'NA'
        qualityStats['nPastEnd'] = 'NA'
        qualityStats['u'] = "TRUE"

        expectedColumns = ['name','short','t0set','num','nMapped','nPastEnd','nGenic','nUsed','gMed','gMedt0','gMean','cor12','mad12','mad12c','mad12c_t0','gccor','opcor','adjcor','maxFit','u']
        qualityStats[expectedColumns].to_csv(outputDirs[outputN]+'FitnessBrowser/fit_quality.tab',sep="\t",index=False)
        
        
        #specific_phenotypes
        browserSpecific = specificPhenotypes.copy()
        browserSpecific.columns = ['locusId','sysName','desc','Condition','lrn','t','Experiment']
        browserSpecific = browserSpecific.merge(expsUsed, left_on='Condition', right_on='Condition', how='left')

        expectedColumns = ['locusId','sysName','desc','name','lrn','t','Group','short','Condition_1','Concentration_1','Units_1','Condition_2','Concentration_2','Units_2']
        browserSpecific[expectedColumns].to_csv(outputDirs[outputN]+'FitnessBrowser/specific_phenotypes',sep="\t",index=False)

        #fit_logratios_good.tab
        newNames = dict(zip(expsUsed['Condition'],expsUsed['name']))
        columnLabels = geneFitAverages.columns.tolist()
        columnLabels = ['locusId','sysName','desc'] + columnLabels[2:]

        newColumnLabels = []
        for idx,columnLabel in enumerate(columnLabels):
            if columnLabel in newNames:
                newColumnLabels.append(newNames[columnLabel])
            else:
                newColumnLabels.append(columnLabel)
                
        #Replace empty cells with zeros for cofitness analysis
        geneFitNoBlank = geneFitAverages.copy()
        noBlanks = geneFitNoBlank[groupSet].mean(axis=1).notna()
        geneFitNoBlank = geneFitNoBlank[noBlanks]
        geneFitNoBlank[groupSet] = geneFitNoBlank[groupSet].fillna(value=0)        
        geneFitNoBlank.to_csv(outputDirs[outputN]+'FitnessBrowser/fit_logratios_good.tab',sep="\t",index=True,index_label=newColumnLabels[0],header=newColumnLabels[1:])

        #fit_t.tab
        TstatsNoBlank = TstatsAverages.copy()
        TstatsNoBlank = TstatsNoBlank[noBlanks]
        TstatsNoBlank[groupSet] =TstatsAverages[groupSet].fillna(value=0)
        TstatsNoBlank.to_csv(outputDirs[outputN]+'FitnessBrowser/fit_t.tab',sep="\t",index=True,index_label=newColumnLabels[0],header=newColumnLabels[1:])

        #.FEBA.success
        successFile = outputDirs[outputN]+'FitnessBrowser/.FEBA.success'
        try:
            with open(successFile, 'w') as FileHandle:
                FileHandle.close()
        except IOError:
            statusUpdate = "Could not write "+successFile+" required by fitness browswer"
            printUpdate(options.logFile,statusUpdate)


    statusUpdate = 'Calculating alternative pValues with Wilcoxon Signed Rank Test'
    printUpdate(options.logFile,statusUpdate)


    #make dataFrame of just rows with counts used in fitness calculations
    wilcoxonTuples = usedCounts.copy()
    wilcoxonPvals = pVals.copy()
    for sampNum in sampleNumbers:
        wilcoxonTuples[sampNum] = list(zip(usedCounts[sampNum], usedCountsRefs[sampNum]))

    for replicateGroup in groupSet:
        analyzedGenes = []
        analyzedPvals = []
        analyzedGroupLabel = []
        for gene, group in wilcoxonTuples.groupby('NearestGene')[sampleNumbers]:
            tupleList = group[groupSampNumbers[replicateGroup]].values.flatten()
            #remove data for unused insertions in this condition
            nonZeroTupleList = [x for x in tupleList if x != (0,0)]
            if len(nonZeroTupleList) > 17:
                analyzedGenes.append(gene)
                condCounts,refCounts = zip(*nonZeroTupleList)
                wilcoxResult = stats.wilcoxon(condCounts,refCounts)[1]
                analyzedPvals.append(wilcoxResult)
                analyzedGroupLabel.append(replicateGroup)
            elif len(nonZeroTupleList) > 9:
                analyzedGenes.append(gene)
                condCounts,refCounts = zip(*nonZeroTupleList)
                wilcoxResult = wilcoxonExact(condCounts,refCounts)
                analyzedPvals.append(wilcoxResult)
                analyzedGroupLabel.append(replicateGroup)

        adjusted = smm.multipletests(analyzedPvals, alpha=0.1, method='fdr_bh')[1]
        
        wilcoxonPvals[replicateGroup] = 1
        for idx,gene in enumerate(analyzedGenes):
            wilcoxonPvals.loc[gene,replicateGroup] = adjusted[idx]

    statusUpdate = 'Saving pValues from Wilcoxon Signed Rank Tests in ' + outputDirs[0]+'pVals_Wilcoxon.txt'
    printUpdate(options.logFile,statusUpdate)
    wilcoxonPvals.to_csv(outputDirs[0]+'pVals_Wilcoxon.txt',sep="\t")


        
def printUpdate(logfile,update):
    timestamp = datetime.now().strftime('%Y-%m-%d %X')
    update = timestamp + ' ' + update
    file = open(logfile, "a")
    file.write(update+'\n')
    file.close()
    print(update)
    return update

def wilcoxonExact(cond1,cond2):
    pcrits = [0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001,0.0005,0.0001]

    Tcrits = {5:[2, 0],6:[3, 2, 0],7:[5, 3, 2, 0],8:[8, 5, 3, 1, 0],
              9:[10, 8, 5, 3, 1, 0],10:[14, 10, 8, 5, 3, 1, 0],
              11:[17, 13, 10, 7, 5, 3, 1, 0],
              12:[21, 17, 13, 9, 7, 5, 2, 1, 0],
              13:[26, 21, 17, 12, 9, 7, 4, 2, 0],
              14:[31, 25, 21, 15, 12, 9, 6, 4, 2, 1],
              15:[36, 30, 25, 19, 15, 12, 8, 6, 4, 1],
              16:[42, 35, 29, 23, 19, 15, 11, 8, 6, 2],
              17:[48, 41, 34, 27, 23, 19, 14, 11, 8, 3],
              18:[55, 47, 40, 32, 27, 23, 18, 14, 10, 4],
              19:[62, 53, 46, 37, 32, 27, 21, 18, 13, 6]}

    diffs = np.array(cond1) - np.array(cond2)

    diffs = [x for x in diffs if x != 0]
    N = len(diffs)

    if N < 5:
        return 1

    else:
    
        diffFrame = pd.DataFrame(diffs,columns=['diffs'])
        diffFrame['abs'] = np.abs(diffs)
        diffFrame['sign'] = diffs/np.abs(diffs)
        signSum = diffFrame['sign'].sum()

        if signSum == 0:
            mindir = 1
        else:
            mindir = -signSum/abs(signSum)

        if mindir == 0:
            return 1
        else:
            sortedFrame = diffFrame.sort_values('abs')
            absList = sortedFrame['abs']
            ranks = list(range(1,N+1))
            sortedFrame['rank'] = ranks

            for absVal,group in sortedFrame.groupby('abs'):
                if len(group) > 1:
                    newRank = group['rank'].mean()
                    sortedFrame.loc[group.index,'rank'] = newRank

            Tw = sortedFrame[sortedFrame['sign'] == mindir]['rank'].sum()

            Pval = 1
            for idx,Tcrit in enumerate(Tcrits[N]):
                if Tw <= Tcrit:
                    Pval = pcrits[idx]

            return Pval 

if __name__ == "__main__":
    main(sys.argv[1:])
exit()

