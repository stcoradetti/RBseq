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
import json
from collections import Counter


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
    options = parser.parse_args()

    statusUpdate = 'RBseq_Count_BarCodes.py  Samuel Coradetti 2019.'
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = 'Version 1.0.8'
    printUpdate(options.logFile,statusUpdate)

    optionDict = options.__dict__
    statusUpdate = 'Options passed: '
    for option in optionDict:
        statusUpdate+=" "+option+":"+str(optionDict[option])+" "
    printUpdate(options.logFile,statusUpdate)
    
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
            metaFrame = metaFrame[metaFrame[metaFrame.columns[0]].notna()]
            FileHandle.close()
    except IOError:
        statusUpdate = "Could not read file:"+fileToOpen+" ...exiting."
        printUpdate(options.logFile,statusUpdate)
        sys.exit()
    
    requiredColumns=['Fastq','SampleName','UsePrecounted','Poolfile','OutputDir','minRandom','maxRandom','DualIndex','BeforeBarcode','BarcodeLengths','AfterBarcode']
    for requiredColumn in requiredColumns:
        if requiredColumn not in metaFrame.columns:
            statusUpdate = "Metadata file is missing column "+requiredColumn+" ...exiting."
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

    #Get a list of fastqs and matching model files from the metadata
    fastqFiles=metaFrame['Fastq'].values
    useCounted=metaFrame['UsePrecounted'].values
    outputDirs=metaFrame['OutputDir'].values
    sampleNames=metaFrame['SampleName'].values
    poolNames=metaFrame['Poolfile'].values
    minRandoms=metaFrame['minRandom'].values
    DualIndexes=metaFrame['DualIndex'].values
    BeforeBarcodes=metaFrame['BeforeBarcode'].values
    BarcodeLenths=metaFrame['BarcodeLengths'].values
    AfterBarcodes=metaFrame['AfterBarcode'].values
    maxRandoms=metaFrame['maxRandom'].values

    for outputN,dir in enumerate(outputDirs):
        if not dir[-1] == "/":
            outputDirs[outputN] = dir+"/"
        if not os.path.exists(outputDirs[outputN]):
            os.makedirs(outputDirs[outputN])
        if not os.path.exists(outputDirs[outputN]+"/countsFiles/"):
            os.makedirs(outputDirs[outputN]+"/countsFiles/")


    #Loop through fastq files and map TnSeq reads to the genome
    statusUpdate = 'Finding barcodes in fastqs and counting occurances'
    printUpdate(options.logFile,statusUpdate)

    allCounts = {}


    for fastqNum,fastqFile in enumerate(fastqFiles):
        sampleName = sampleNames[fastqNum]
        noPreSeq = 0
        noPostSeq = 0
        noDualIndex = 0
        nonCompliantBarcode = 0
        poorQualityBarcode = 0
        if useCounted[fastqNum]:
            statusUpdate = '  Loading previously counted reads for '+fastqFile
            printUpdate(options.logFile,statusUpdate)
            fileToOpen = outputDirs[fastqNum]+"countsFiles/"+sampleName+'.counts'
            try:
                with open(fileToOpen, 'r') as file:
                    barcodesCountsLoad = Counter(json.load(file))
                    file.close()
            except IOError:
                statusUpdate = " Could not read file:"+fileToOpen+" ...exiting."
                printUpdate(options.logFile,statusUpdate)



            if sampleName in allCounts:
                barcodeCounts = allCounts[sampleName] + barcodesCountsLoad
            else:
                barcodeCounts = Counter(barcodesCountsLoad)


        else:
            statusUpdate = '  Mapping reads from '+fastqFile
            printUpdate(options.logFile,statusUpdate)

            #first check if this fastq is for a sample with counts already loaded from another sequence file
            if sampleNames[fastqNum] in allCounts:
                barcodeCounts = allCounts[sampleNames[fastqNum]]
            else:
                barcodeCounts=Counter({})

            skipBases = minRandoms[fastqNum]
            dualIndex = DualIndexes[fastqNum]

            if not bool(re.match('^[AGCT]+$', str(dualIndex))):
                dualIndex = ""
            else:
                statusUpdate = "  Filtering reads on second index: "+str(dualIndex)
                printUpdate(options.logFile,statusUpdate)

            beforeBarcode = BeforeBarcodes[fastqNum]
            phasingBases = maxRandoms[fastqNum] - skipBases

            if len(beforeBarcode) < options.matchBefore:
                statusUpdate = " Metafile defines a shorter sequence before the barcode than required by --matchBefore... exiting"
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            #Define the short sequences to search for flanking barcodes
            searchBefore = beforeBarcode[-options.matchBefore:]
            

            #Find the expected positions for search sequences
            offsetBefore = skipBases + len(dualIndex) + len(beforeBarcode) - len(searchBefore)
            offsetBarcode = skipBases + len(dualIndex) + len(beforeBarcode)

            endSearch = AfterBarcodes[fastqNum][:options.matchAfter]
            matchAfterWarning = 0

            allowedLengths = [ int(x) for x in BarcodeLenths[fastqNum].split(",") ]

            statusUpdate = "  Looking for sequence barcode preceeded by "+searchBefore
            printUpdate(options.logFile,statusUpdate)

            fileToOpen = fastqFile
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


            if matchAfterWarning:
                statusUpdate = "    Warning: " + str(matchAfterWarning) + " reads were too short to check " + str(options.matchAfter) + " bases after barcodes as specified with --matchAfter." 
                printUpdate(options.logFile,statusUpdate)
                statusUpdate = "      At least two bases following these barcodes matched the expected sequence." 
                printUpdate(options.logFile,statusUpdate)


            statusUpdate = "    " + str(NbarcodeNotFound) + " reads without recognizable, compliant barcodes."
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
            statusUpdate = "      " + str(poorQualityBarcode) + " reads with quality scores less than " + str(options.minQual) + "."
            printUpdate(options.logFile,statusUpdate)

            
            statusUpdate = "    " + str(NbarcodeFound) + " reads with compliant barcodes"
            printUpdate(options.logFile,statusUpdate)

        allCounts[sampleNames[fastqNum]] = barcodeCounts

    statusUpdate = "  Saving barcode counts to " + outputDirs[0]+"countsFiles/"
    printUpdate(options.logFile,statusUpdate)

    for sampleName in allCounts:
        fileToOpen = outputDirs[fastqNum]+"countsFiles/"+sampleName+'.counts'
        try:
            with open(fileToOpen, 'w') as file:
                file.write(json.dumps(allCounts[sampleName]))
        except IOError:
            statusUpdate = " Could not read file:"+fileToOpen+" ...exiting."
            printUpdate(options.logFile,statusUpdate)

    statusUpdate = "  Combining counts over all samples"
    printUpdate(options.logFile,statusUpdate)

    barcodeTotals = {}
    for sample in allCounts:
        for barcode in allCounts[sample]:
            if barcode in barcodeTotals:
                barcodeTotals[barcode] = barcodeTotals[barcode] + allCounts[sample][barcode]
            else:
                barcodeTotals[barcode] = allCounts[sample][barcode]

    allBarcodes = list(barcodeTotals.keys())

    combinedCounts = pd.DataFrame(allBarcodes,columns=['barcode'])
    for sample in allCounts:
        sampleCounts = []
        for barcode in allBarcodes:
            if barcode in allCounts[sample]:
                sampleCounts.append(allCounts[sample][barcode])
            else:
                sampleCounts.append(0)
        combinedCounts[sample] = sampleCounts

    

    #estimate error rate in barcodes and number of true barcodes with one or two counts
 
    totalCounts = combinedCounts.sum(axis=1).values
    combinedCounts['totalCounts'] = totalCounts
    totalCountSum = np.sum(totalCounts)

    statusUpdate = "Sequence Files processed: " + str(len(fastqFiles)) 
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = "Total reads with compliant barcodes: " + str(totalCountSum)
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = "Total barcodes seen (incudes sequencing errors): " + str(len(combinedCounts))
    printUpdate(options.logFile,statusUpdate)


    allBarcodes = combinedCounts['barcode'].values
    maxIndex = np.argmax(totalCounts)
    maxBarcode = allBarcodes[maxIndex]
    maxCounts = totalCounts[maxIndex]

    errorList = OffByOneList(maxBarcode)
    errorCounts = 0
    for index in range(0,len(totalCounts)):
        if allBarcodes[index] in errorList:
            errorCounts += totalCounts[index]
    errorRate = max(.01,float(errorCounts)/(maxCounts+errorCounts))

    statusUpdate = "Most abundant barcode: " + maxBarcode + " seen " + str(maxCounts) + " times."
    printUpdate(options.logFile,statusUpdate)

    if maxCounts > 10000:
        statusUpdate = "Number of reads that differ from this barcode by one base pair (likely sequencing errors):" + str(errorCounts)
        printUpdate(options.logFile,statusUpdate)
        statusUpdate = "Estimated sequencing error rate for barcodes: " + str(errorRate*100) + "%"
        printUpdate(options.logFile,statusUpdate)

    else:
        statusUpdate = "Not enough reads for most common barcode to estimate sequencing error rate in barcodes. Assuming 1%."
        printUpdate(options.logFile,statusUpdate)
        erroRate = 0.01
        

    bins = range(1,maxCounts+2)
    countHistogram, division = np.histogram(totalCounts, bins = bins)
    old_countHistogram = countHistogram.copy()

    #For a barcode with given number of counts, remove expected number of erroneous counts from the barcode distribution, starting with most abundant barcodes
    for counts in range(maxCounts,0,-1):
        if countHistogram[counts-1] > 0:
            lambdaCounts = counts*errorRate/60
            if lambdaCounts >= 1:
                nSeqErrBarcodes = 60 * countHistogram[counts-1]
                SeqErrBarcodsCounts = np.random.poisson(lambdaCounts,nSeqErrBarcodes)
                for ErrBarcodeCount in SeqErrBarcodsCounts:
                    countHistogram[ErrBarcodeCount-1] -= 1
            else:
                countHistogram[0] -= int(countHistogram[counts-1] * errorRate)

    twos = countHistogram[1]
    ones = countHistogram[0]

    #Chao population size estimate, conservative
    Nch = ones**2/(2*twos)
    Nch = 10**int(np.log10(Nch))

    statusUpdate = "  Estimate of real barcodes seen once: " + str(ones)
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = "  Estimate of real barcodes seen twice: " + str(twos)
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = "  Chao estimate of true population size (VERY approximate): " + str(Nch)
    printUpdate(options.logFile,statusUpdate)

    poolFileName = poolNames[0]
    statusUpdate = "Compiling counts for barcodes in " + poolFileName
    printUpdate(options.logFile,statusUpdate)

    try:
        with open(poolFileName, 'rb') as poolFileHandle:
            poolFrame = pd.read_csv(poolFileHandle,low_memory=False,dtype={'NearestGene':str,'CodingFraction':str},sep='\t')
            poolFileHandle.close()
            poolFrame.dropna(how='all')
            statusUpdate =  "Read "+ str(len(poolFrame)) + " barcodes from "+poolFileName
            printUpdate(options.logFile,statusUpdate)
    except IOError:
        statusUpdate =  "Could not read file: "+poolFileName
        printUpdate(options.logFile,statusUpdate)
        sys.exit()

    poolCounts = poolFrame[['rcbarcode','scaffold','pos','NearestGene','CodingFraction','LocalGCpercent','AlternateID','Annotation']].copy()
    poolCounts = poolCounts.rename(index=str, columns={'rcbarcode':'barcode'})
    poolCounts = poolCounts.merge(combinedCounts,how='left',on='barcode').fillna(0)

    for sample in sampleNames:
        poolCounts[sample] = poolCounts[sample].astype(int)

    totalCounts = poolCounts[sampleNames].sum(axis=1)

    nSeenInPool = totalCounts.astype(bool).sum()
    nCounts = totalCounts.sum()


    statusUpdate = "Counted " + str(nCounts) + " reads across all BarSeq samples for " + str(nSeenInPool) + " mapped barcodes from the mutant pool."
    printUpdate(options.logFile,statusUpdate)


    poolCounts.to_csv(outputDirs[0]+'/poolCount.txt',sep="\t",index=False)

    statusUpdate = "Wrote results to: " + outputDirs[0] + "poolCount.txt"
    printUpdate(options.logFile,statusUpdate)

if __name__ == "__main__":
    main(sys.argv[1:])
exit()
