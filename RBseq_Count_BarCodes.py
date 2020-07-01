#!/usr/bin/python
import sys
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import re
import os
import json

Version = '1.1.3'
ReleaseDate = 'July 1, 2020'


def OffByOneList(seq):
    if seq[0] in ("A","T","G","C"):
        char_set = ("A","T","G","C")
    elif seq[0] in ("a","t","g","c"):
        char_set = ("a","t","g","c")
    else:
        return False

    variants = []
    for chari in range(len(seq)):
        if chari == 0:
            preseq = ""
        else:
            preseq = seq[0:chari]
        if chari == len(seq)-1:
            postseq = ""
        else:
            postseq = seq[chari+1:]
        for char in char_set:
            if seq[chari] != char:
                variants.append(preseq+char+postseq)

    return(variants)

def printUpdate(logfile,update):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    update = timestamp + ' ' + update
    file = open(logfile, "a")
    file.write(update+'\n')
    file.close()
    print(update)
    return update


def main(argv):
    timestamp = datetime.now().strftime('%Y%m%H%M')
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metafile", dest="metafile", help="Metadata file for BarSeq runs. A tab-delimited file with columns titled Fastq,SampleName,BarSeqModel,Poolfile,OutputDir.",default="metafile.txt")
    parser.add_argument("-l", "--logFile", dest="logFile", help="File to write run log to. Default is Count_TIMESTAMP.log",default="Count_"+timestamp+".log")
    parser.add_argument("-q", "--qual", dest="minQual", help="Minimum quality score for the barcode region for a read to counted. Default is 10",default=10, type=int)
    parser.add_argument("-b", "--matchBefore", dest="matchBefore", help="Number of bases before the barcode to match. Default is 6", default=6, type=int)
    parser.add_argument("-a", "--matchAfter", dest="matchAfter", help="Number of bases after the barcode to match. Default is 6", default=6, type=int)
    parser.add_argument("-Q", "--quietMode", action='store_true', dest="quietMode", help="Give fewer details about barcodes found during the run and in the log file.  Summary stats will still be reported in summaryStats.txt", default=False)
    options = parser.parse_args()

    statusUpdate = 'RBseq_Count_BarCodes.py'
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = 'Version: ' + Version
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = 'Release Date: ' + ReleaseDate
    printUpdate(options.logFile,statusUpdate)

    optionDict = options.__dict__
    statusUpdate = 'Options passed: '
    for option in optionDict:
        statusUpdate+=" "+option+":"+str(optionDict[option])+" "
    printUpdate(options.logFile,statusUpdate)

    if (not options.quietMode):
        statusUpdate = 'Logging status updates in '+options.logFile
        printUpdate(options.logFile,statusUpdate)
        statusUpdate = 'Loading TnSeq library metadata from '+options.metafile
        printUpdate(options.logFile,statusUpdate)
     
    #Load up experiment metadata
    fileToOpen = options.metafile
    if (not fileToOpen):
        statusUpdate = options.metafile + ' could not be loaded... exiting.'
        printUpdate(options.logFile,statusUpdate)
        sys.exit()
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            metaFrame = pd.read_csv(FileHandle,sep='\t')
            FileHandle.close()
            metaFrame = metaFrame[metaFrame[metaFrame.columns[0]].notna()]
            metaFrame = metaFrame.fillna("NO INPUT")
            
    except IOError:
        statusUpdate = "Could not read file:"+fileToOpen+" ...exiting."
        printUpdate(options.logFile,statusUpdate)
        sys.exit()
    
    requiredColumns=['FileIndex','Fastq','SampleName','UsePrecounted','Poolfile','OutputDir','minRandom','maxRandom','DualIndex','BeforeBarcode','BarcodeLengths','AfterBarcode']
    for requiredColumn in requiredColumns:
        if requiredColumn not in metaFrame.columns:
            statusUpdate = "Metadata file is missing column "+requiredColumn+" ...exiting."
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

    #Get a list of fastqs and matching model files from the metadata
    fileIndexes=metaFrame['FileIndex'].values
    metaFrame = metaFrame.set_index('FileIndex')
    

    outputDir = metaFrame.loc[fileIndexes[0]]['OutputDir']
    if (not options.quietMode):
        statusUpdate = "Setting output directory as: " + outputDir + " (from first line of metadata file)"
        printUpdate(options.logFile,statusUpdate)
    
    if not outputDir[-1] == "/":
        outputDir = outputDir+"/"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    if not os.path.exists(outputDir+"countsFiles/"):
        os.makedirs(outputDir+"countsFiles/")

    statList = ['FileIndex','TotalReads','BarcodeFoundReads','BarcodeInPoolfileReads','BarcodeNotFoundReads','NoDualIndexReads',
                'NoPreSeqReads','NoPostSeqReads','NonCompliantBarcodeReads','PoorQualityBarcodeReads','TotalBarcodes','TotalBarcodesInPoolfile',
                'MostAbundantBarcode','CountsForMostAbudantBarcode','ErrorRatePercent','OneRead','TwoReads','ThreeReadsOrMore','EstimatedPopulationSize']
               
        
    sumStats = {}
    for stat in statList:
        sumStats[stat] = []
        

    poolFileName = str(metaFrame.loc[fileIndexes[0]]['Poolfile'])
    if (not options.quietMode):
        statusUpdate = "Loading mapped barcodes in mutant pool from: " + poolFileName + " (from first line of metadata file)"
        printUpdate(options.logFile,statusUpdate)
    poolFileFound = False
    if not (poolFileName == "NO INPUT"):
        try:
            with open(poolFileName, 'rb') as poolFileHandle:
                poolFrame = pd.read_csv(poolFileHandle,low_memory=False,dtype={'NearestGene':str,'CodingFraction':str},sep='\t')
                poolFileHandle.close()
                poolFrame.dropna(how='all')
                statusUpdate =  "Read "+ str(len(poolFrame)) + " barcodes from "+poolFileName
                printUpdate(options.logFile,statusUpdate)
                poolFileFound = True
        except IOError:
            statusUpdate =  "Could not read file: "+poolFileName
            printUpdate(options.logFile,statusUpdate)
            sys.exit()
  
    if poolFileFound:
        poolCounts = poolFrame[['rcbarcode','scaffold','pos','NearestGene','CodingFraction','LocalGCpercent','AlternateID','Annotation']].copy()
        poolCounts = poolCounts.rename(index=str, columns={'rcbarcode':'barcode'})
        poolCounts_byIndex = poolCounts.copy()
    else:
        statusUpdate =  "No valid poolfile found on the first line of metadata.  Barcodes will be extracted from sequence data, but no poolcounts file will be generated."
        printUpdate(options.logFile,statusUpdate)

    if (not options.quietMode):
        statusUpdate = 'Finding barcodes in fastqs and counting occurances'
        printUpdate(options.logFile,statusUpdate)
        statusUpdate = "---------------------"
        printUpdate(options.logFile,statusUpdate)

    for fileIndex in fileIndexes:
        sumStats['FileIndex'].append(fileIndex)
        indexRow = metaFrame.loc[fileIndex]
        noPreSeq = 0
        noPostSeq = 0
        noDualIndex = 0
        nonCompliantBarcode = 0
        poorQualityBarcode = 0
        
        barcodeCounts={}

        if  indexRow['UsePrecounted']:
            statusUpdate = '  Loading previously counted reads for '+ indexRow['OutputDir']
            printUpdate(options.logFile,statusUpdate)
            fileToOpen = outputDir+"countsFiles/"+str(fileIndex)+'.counts'
            try:
                with open(fileToOpen, 'r') as file:
                    barcodeCounts = json.load(file)
                    file.close()
            except IOError:
                statusUpdate = " Could not read file:"+fileToOpen+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            for stat in ['TotalReads','BarcodeFoundReads','BarcodeNotFoundReads','NoDualIndexReads','NoPreSeqReads','NoPostSeqReads',
                         'NonCompliantBarcodeReads','PoorQualityBarcodeReads']:
                sumStats[stat].append(np.nan)


        else:
            statusUpdate = '  Mapping reads from ' + indexRow['Fastq']
            printUpdate(options.logFile,statusUpdate)

            skipBases = indexRow['minRandom']
            phasingBases = indexRow['maxRandom'] - skipBases
            
            dualIndex = indexRow['DualIndex']
            if not bool(re.match('^[AGCT]+$', str(dualIndex))):
                dualIndex = ""
            else:
                if (not options.quietMode):
                    statusUpdate = "  Filtering reads on second index: "+str(dualIndex)
                    printUpdate(options.logFile,statusUpdate)

            beforeBarcode = indexRow['BeforeBarcode']

            if len(beforeBarcode) < options.matchBefore:
                statusUpdate = " Metafile defines a shorter sequence before the barcode than required by --matchBefore... exiting"
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            #Define the short sequences to search for flanking barcodes
            searchBefore = beforeBarcode[-options.matchBefore:]
            

            #Find the expected positions for search sequences
            offsetBefore = skipBases + len(dualIndex) + len(beforeBarcode) - len(searchBefore)
            offsetBarcode = skipBases + len(dualIndex) + len(beforeBarcode)

            endSearch = indexRow['AfterBarcode'][:options.matchAfter]
            matchAfterWarning = 0

            allowedLengths = [ int(x) for x in indexRow['BarcodeLengths'].split(",") ]

            if (not options.quietMode):
                statusUpdate = "  Looking for sequence barcode preceeded by "+searchBefore
                printUpdate(options.logFile,statusUpdate)

            fileToOpen = indexRow['Fastq']
            try:
                with open(fileToOpen, 'r') as FileHandle:
                    readName=""
                    readCount=0
                    readSeq=""
                    readQual=""
                    NbarcodeFound = 0
                    NbarcodeNotFound = 0
                    lowQualScores = {}

                    for readLine in FileHandle:
                        compliantBarcode=False
                        if readLine[0] == "@":
                            readName = readLine.strip()
                            readCount = 0
                        elif readCount == 1:
                            readSeq = readLine.strip()

                        #Third line after a readname is the quality score, now weve loaded info for read and can filter for barcode presence before moving to next read.
                        elif readCount == 3:
                            readQual = readLine.strip()

                            passIndex = True
                            if not dualIndex == "":
                                if readSeq[skipBases:skipBases + len(dualIndex) + phasingBases].find(dualIndex) == -1:                           
                                    passIndex = False
                                    noDualIndex+=1

                            compliantBarcode = False
                            if passIndex:
                                matchedBefore = readSeq[offsetBefore:offsetBarcode + phasingBases].find(searchBefore)
                                if matchedBefore > -1:
                                    prebarcode = readSeq[offsetBarcode + matchedBefore:]
                                    n = 0
                                    while n < len(allowedLengths) and not compliantBarcode:
                                        allowedlength = allowedLengths[n]
                                        afterbarcode = prebarcode[allowedlength:]
                                        searchableBases = len(afterbarcode)

                                        if searchableBases < options.matchAfter:
                                            matchAfterWarning += 1
                                            searchAfter = searchableBases
                                        else:
                                            searchAfter = options.matchAfter

                                        if searchAfter > 1 and afterbarcode[:searchAfter] == endSearch[:searchAfter]:
                                            barcode = prebarcode[:allowedlength]
                                            if bool(re.match('^[AGCT]+$', barcode)):
                                                compliantBarcode = True
                                            else:
                                                nonCompliantBarcode+=1
                                        elif n == len(allowedLengths) - 1:
                                            noPostSeq+=1
                                        n += 1
                                    barQual = readQual[offsetBarcode + matchedBefore:offsetBarcode + matchedBefore + allowedlength]
                                    if min(map(ord,list(barQual)))-33 < options.minQual:
                                        compliantBarcode = False
                                        poorQualityBarcode+=1
                                else:
                                    noPreSeq+=1


                            if compliantBarcode:
                                NbarcodeFound+=1
                                if barcode in barcodeCounts:
                                    barcodeCounts[barcode]+=1
                                else:
                                    barcodeCounts[barcode]=1
                            else:
                                NbarcodeNotFound+=1

                        readCount+=1

                    FileHandle.close()


            except IOError:
                statusUpdate = "  Could not read file:"+fileToOpen+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()


            sumStats['TotalReads'].append(NbarcodeNotFound + NbarcodeFound)
            sumStats['BarcodeNotFoundReads'].append(NbarcodeNotFound)
            sumStats['BarcodeFoundReads'].append(NbarcodeFound)
            sumStats['NoDualIndexReads'].append(noDualIndex)
            sumStats['NoPreSeqReads'].append(noPreSeq)
            sumStats['NoPostSeqReads'].append(noPostSeq)
            sumStats['NonCompliantBarcodeReads'].append(nonCompliantBarcode)
            sumStats['PoorQualityBarcodeReads'].append(poorQualityBarcode)

            statusUpdate = "    " + str(NbarcodeNotFound + NbarcodeFound) + " reads processed."
            printUpdate(options.logFile,statusUpdate)

            if matchAfterWarning:
                statusUpdate = "    " + str(matchAfterWarning) + " reads were too short to check " + str(options.matchAfter) + " bases after barcodes as specified with --matchAfter.  For these reads all remaining bases will be checked against the expected sequence. Reads with less than two remaining bases after the barcode will be counted as non-complaint." 
                printUpdate(options.logFile,statusUpdate)

            statusUpdate = "    " + str(NbarcodeFound) + " reads with compliant barcodes."
            printUpdate(options.logFile,statusUpdate)

            if (not options.quietMode):
                statusUpdate = "    " + str(NbarcodeNotFound) + " reads without recognizable, compliant barcodes. Of those:"
                printUpdate(options.logFile,statusUpdate)
                if noDualIndex:
                    statusUpdate = "      " + str(noDualIndex) + " reads without expected Dual Index."
                    printUpdate(options.logFile,statusUpdate)
                statusUpdate = "      " + str(noPreSeq) + " reads without expected sequence before the barcode region."
                printUpdate(options.logFile,statusUpdate)
                statusUpdate = "      " + str(noPostSeq) + " reads without expeced sequence after the barcode region."
                printUpdate(options.logFile,statusUpdate)
                statusUpdate = "      " + str(nonCompliantBarcode) + " reads with noncompliant barcdes. (Contains Ns, etc)."
                printUpdate(options.logFile,statusUpdate)
                statusUpdate = "      " + str(poorQualityBarcode) + " reads with quality scores less than " + str(options.minQual)
                printUpdate(options.logFile,statusUpdate)

            fileToSave = outputDir+"countsFiles/"+str(fileIndex)+'.counts'
            if (not options.quietMode):
                statusUpdate = "  Saving barcode counts to " + fileToSave
                printUpdate(options.logFile,statusUpdate)

            try:
                with open(fileToSave, 'w') as file:
                    file.write(json.dumps(barcodeCounts))
            except IOError:
                statusUpdate = " Could not read file:"+fileToSave+" ...exiting."
                printUpdate(options.logFile,statusUpdate)

        sumStats['TotalBarcodes'].append(len(barcodeCounts))
        statusUpdate = "  Total barcodes seen (incudes sequencing errors): " + str(len(barcodeCounts))
        printUpdate(options.logFile,statusUpdate)

        if poolFileFound:
            if (not options.quietMode):
                statusUpdate = "  Matching barcodes to poolfile"
                printUpdate(options.logFile,statusUpdate)

            barcodeCountsFrame = pd.DataFrame.from_dict(barcodeCounts, orient='index',columns=[fileIndex],dtype=int)
            poolCounts_byIndex[fileIndex] = barcodeCountsFrame.reindex(poolCounts_byIndex['barcode'],fill_value=0)[fileIndex].values

            totalReadsInPool = poolCounts_byIndex[fileIndex].sum()
            NbarcodesInPool = len(poolCounts_byIndex[poolCounts_byIndex[fileIndex] > 0])

            statusUpdate = "    Number of barcodes from poolfile seen: " + str(NbarcodesInPool)
            printUpdate(options.logFile,statusUpdate)
            if (not options.quietMode):
                statusUpdate = "    Reads with barcodes from poolfile: " + str(totalReadsInPool)
                printUpdate(options.logFile,statusUpdate)

            sumStats['TotalBarcodesInPoolfile'].append(NbarcodesInPool)
            sumStats['BarcodeInPoolfileReads'].append(totalReadsInPool)
        else:
            sumStats['TotalBarcodesInPoolfile'].append(np.nan)
            sumStats['BarcodeInPoolfileReads'].append(np.nan)

        mostAbundantBarcode = barcodeCountsFrame[fileIndex].idxmax()
        mostAbundantCounts = barcodeCountsFrame.loc[mostAbundantBarcode][fileIndex]

        sumStats['MostAbundantBarcode'].append(mostAbundantBarcode)
        sumStats['CountsForMostAbudantBarcode'].append(mostAbundantCounts)

        if not options.quietMode:
            statusUpdate = "  Most abundant barcode: " + mostAbundantBarcode + " seen " + str(mostAbundantCounts) + " times."
            printUpdate(options.logFile,statusUpdate)
            

        if mostAbundantCounts > 1000:
            errorList = OffByOneList(mostAbundantBarcode)
            errorsSeen = []
            for barcode in errorList:
                if barcode in  barcodeCountsFrame.index:
                    errorsSeen.append(barcode)
            errorCounts = barcodeCountsFrame.loc[errorsSeen][fileIndex].sum()
            errorRate = max(.005,float(errorCounts)/(mostAbundantCounts+errorCounts))*100

            if not options.quietMode:
                statusUpdate = "  Number of reads that differ from this barcode by one base pair (likely sequencing errors):" + str(errorCounts)
                printUpdate(options.logFile,statusUpdate)

        else:
            if not options.quietMode:
                statusUpdate = "  Not enough reads for most common barcode to estimate sequencing error rate in barcodes. Assuming 1%."
                printUpdate(options.logFile,statusUpdate)
            errorRate = 1.0

        sumStats['ErrorRatePercent'].append(errorRate)

        ones = (barcodeCountsFrame[fileIndex] == 1).sum()
        twos = (barcodeCountsFrame[fileIndex] == 2).sum()
        threes = len(barcodeCountsFrame) - ones - twos

        sumStats['OneRead'].append(ones)
        sumStats['TwoReads'].append(twos)
        sumStats['ThreeReadsOrMore'].append(threes)       

        #Chao population size estimate, conservative
        Nch = ones**2/(2*twos)

        #Round to two significant figures
        errorRate = round(errorRate, -(int(np.floor(np.log10(abs(errorRate)))))+1)
        Nch = int(round(Nch, -(int(np.floor(np.log10(abs(Nch)))))+1))
        sumStats['EstimatedPopulationSize'].append(Nch)

        if (not options.quietMode):
            statusUpdate = "  Estimated sequencing error rate for barcodes: " + str(errorRate) + "%"
            printUpdate(options.logFile,statusUpdate)
            statusUpdate = "  Barcodes seen once (highly inflated by sequencing errors): " + str(ones)
            printUpdate(options.logFile,statusUpdate)
            statusUpdate = "  Barcodes seen twice (slightly inflated by sequencing errors): " + str(twos)
            printUpdate(options.logFile,statusUpdate)
            statusUpdate = "  Barcodes seen three times or more: " + str(threes)
            printUpdate(options.logFile,statusUpdate)
            statusUpdate = "  Chao estimate of population size (ones^2/2*twos): " + str(Nch)
            printUpdate(options.logFile,statusUpdate)


        if (not options.quietMode):
            statusUpdate = "---------------------"
            printUpdate(options.logFile,statusUpdate)
            
    sumStatFrame = pd.DataFrame.from_dict(sumStats)
    fileToSave = outputDir+"fastqSummaryStats.txt"
    statusUpdate = "Saving summary statistics for fastqs to: " + fileToSave
    printUpdate(options.logFile,statusUpdate)
    sumStatFrame.to_csv(fileToSave,sep='\t',index=None)
            
    sampleGroups = metaFrame.groupby("SampleName")
    sampleList = list(sampleGroups.groups.keys())

    statusUpdate = "Processed " + str(len(fileIndexes)) + " sequence files corresponding to " + str(len(sampleList)) + " biological samples."
    printUpdate(options.logFile,statusUpdate)
            
    if poolFileFound:
        for SampleName,group in sampleGroups:
            poolCounts[SampleName] = poolCounts_byIndex[group.index].sum(axis=1)

        totalCounts = poolCounts[sampleList].sum(axis=1)
        nSeenInPool = totalCounts.astype(bool).sum()
        nCounts = totalCounts.sum()

        statusUpdate = "Counted " + str(nCounts) + " total reads for " + str(nSeenInPool) + " mapped barcodes from the mutant pool."
        printUpdate(options.logFile,statusUpdate)

        fileToSave = outputDir+"poolCount.txt"
        statusUpdate = "Saving poolcount file to: " + fileToSave
        printUpdate(options.logFile,statusUpdate)
        poolCounts.to_csv(fileToSave,sep='\t',index=None)
    


if __name__ == "__main__":
    main(sys.argv[1:])
exit()
