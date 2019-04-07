#!/usr/bin/python
from __future__ import print_function
import sys
import argparse
import numpy as np
import operator
import pandas as pd
from datetime import datetime
import re
from Bio.Blast.Applications import NcbiblastnCommandline
import os

#http://code.activestate.com/recipes/576874-levenshtein-distance/
def levenshtein(s1, s2):
    l1 = len(s1)
    l2 = len(s2)
    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in list(range(l2 + 1)):
        matrix[zz] = list(range(zz,zz + l1 + 1))
    for zz in list(range(0,l2)):
        for sz in list(range(0,l1)):
            if s1[sz] == s2[zz]:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
            else:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]

def ReverseComplement(seq):
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = {}
    for i in list(range(16)):
        if i < 4 or 8<=i<12:
            seq_dict[seq1[i]] = seq1[i+4]
    return "".join([seq_dict[base] for base in reversed(seq)])


def OffByOneList(seq):
    if seq[0] in ("A","T","G","C"):
        char_set = ("A","T","G","C")
    elif seq[0] in ("a","t","g","c"):
        char_set = ("a","t","g","c")
    else:
        return False

    variants = []
    for chari in list(range(len(seq))):
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
    timestamp = datetime.now().strftime('%Y-%m-%d %X')
    update = timestamp + ' ' + update
    file = open(logfile, "a")
    file.write(update+'\n')
    file.close()
    print(update)
    return update

def matchExpected(searchSeq,readSeq,offset,wobbleAllowed):
    matched=False
    tryPositions = []
    for tryPosition in list(range(0,wobbleAllowed)):
        tryPositions.append(tryPosition)
        if tryPosition > 0:
            if offset >= tryPosition:
                tryPositions.append(-tryPosition)
    tryInd = 0
    while (not matched) and tryInd < len(tryPositions):
        tryPosition = tryPositions[tryInd]+offset
        matched = (readSeq[tryPosition:tryPosition+len(searchSeq)] == searchSeq)
        tryInd+=1
    return [matched,tryPosition-offset]
    

def main(argv):
    timestamp = datetime.now().strftime('%Y%m%H%M')
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--metafile", dest="metafile", help="Metadata file of TnSeq runs. A tab-delimited file with columns titled Pool, ShortName, Fastq, ReadModel, InsertionSequence, GenomeSequence, OutputDir, GeneLocations, GeneAnnotations, GeneIdentifier.  See README for more detail",default="metafile.txt")
    parser.add_argument("-l", "--logFile", dest="logFile", help="File to write run log to. Default is Map_TIMESTAMP.log",default="Map_"+timestamp+".log")
    parser.add_argument("-q", "--qual", dest="minQual", help="Minimum quality score for the barcode region for a read to be mapped. Default is 10",default=10)
    parser.add_argument("-b", "--matchBefore", dest="matchBefore", help="Number of bases before the barcode to match. Default is 6", default=6, type=int)
    parser.add_argument("-a", "--matchAfter", dest="matchAfter", help="Number of bases in the constant region after the barcode to match. Default is 6. Note that if read is too short to search 6 bases after the barcode, but at least 3 bases remain to be matched, then those bases will be used.", default=6, type=int)
    parser.add_argument("-j", "--matchJunction", dest="matchJunction", help="Number of bases to match at the junction with the genome. Default is 5. Because under many library construction strategies there is often loss of sequence at the end of the inserted fragment, this script attempts to find a short segment of sequence matching the end of the insertion, scanning for the closest matching segment to the end of the sequence.  Sequence after that furthest match is used as the junction with the native genome sequence.", default=5, type=int)
    parser.add_argument("-f", "--maxFillerSeq", dest="maxFillerSeq", help="NOT IMPLEMENTED YET. The maximum number of bases allowed between the detected end of the insertion and the start of a match to the genome. This could be used in future versions to report insertions with long pieces of filler DNA added at insertion sites.",default=100, type=int)
    parser.add_argument("-v", "--barcodeVariation", dest="barcodeVariation", help="Allowed deviation from expected barcode length. Default is 2. Depending on the method of library construction, significant numbers of barcoded insertions may have slightly shorter or longer random barcodes.  Allowing flexibility in barcode length can recover these insertions, but slows down processing.", default=2, type=int)
    parser.add_argument("-w", "--wobbleAllowed", dest="wobbleAllowed", help="Allowed deviation from expected position in the read for the matching sequence around the barcode and the junction with the genome. Default is 1.  Indels introduced near the barcode can leave some insertion unmapped without some allowance for small deviations from expected positions of search sequences.", default=1, type=int)
    parser.add_argument("-s", "--scoreDiff", dest="scoreDiff", help="For a read mapping to the genome to be counted as unique, the best-hit bitscore must be at least this much greater than the second-best hit bitscore. Default is 10", default=10, type=int)
    parser.add_argument("-F", "--minFrac", dest="minFraction", help="Minimum fraction of reads mapping to a single location before barcode is classified as multilocus. Default is 0.9", default=0.9, type=float)
    parser.add_argument("-N", "--filterNeighborhood", dest="filterNeighborhood", help="Number of basepairs around each insertion to filter out less abundant insertions with similair barcodes. Default is 10", default=10, type=int)
    parser.add_argument("-D", "--filterEditDistance", dest="filterEditDistance", help="Nearby barcodes within this edit distance will be filtered. Default is 5. This argument is ignored if --noBarcodes is passed (all insertions within --filterNeighborhood will be filtered)", default=5, type=int)
    parser.add_argument("-p", "--minPercentID", dest="minPercentID", help="Minimum percent identiy required when BLASTing sequences to the genome or insert. Default is 95", default=95, type=int)
    parser.add_argument("-e", "--maxEvalue", dest="maxEvalue", help="Maximum Evalue allowed when BLASTing sequences to the genome or insert. Default is 0.1", default=0.1, type=float)
    parser.add_argument("-u", "--useMappedFiles", dest="useMappedFiles",  action='store_true', help="Bypass the mapping step and work from previously mapped files in the output directory", default=False)
    parser.add_argument("-I", "--noInsertHits", dest="noInsertHits",action='store_true', help="Skip search for additional insertion sequence after the expected junction with the genome. Default behavior is to look for these sequences, typically indicating concatameric insertions. Turning it off will save time if concatamers aren't common in your mutant pool.", default=False)
    parser.add_argument("-n", "--noBarcodes", dest="noBarcodes", action='store_true', help="This flag is used for optimizing library construction protocols with non-barcoded sequences.  As generating a diverse pool of insertion constructs with unique barcodes can be time consuming, a test insertion library can be created with constructs containing just one constant barcode. If TnSeq is performed on the resulting library, the distinct locations of insertions can be mapped, but they will all have the same sequence barcode associated with the insertion.  To allow characterization of the library with this script, this flag signals to use a snippet of the associated genomic sequence as a substitute for a sequence barcode on the insertion.  This allows for distinct insertions to be mapped and estimation of library diversity and identification of biases in insertion position, etc.", default=False)

    options = parser.parse_args()

    statusUpdate = 'RBseq_Map_Insertions.py  Samuel Coradetti 2019.'
    printUpdate(options.logFile,statusUpdate)
    statusUpdate = 'Version 1.0.3'
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
        with open(fileToOpen, 'r') as FileHandle:
            metaFrame = pd.read_table(FileHandle,low_memory=False)
            metaFrame = metaFrame[~metaFrame[metaFrame.columns[0]].isnull()]
            FileHandle.close()
    except IOError:
        statusUpdate = "Could not read file:"+fileToOpen+" ...exiting."
        printUpdate(options.logFile,statusUpdate)
        sys.exit()

    requiredColumns=['Pool','ShortName','Fastq','ReadModel','InsertionSequence','GenomeSequence','OutputDir']
    for requiredColumn in requiredColumns:
        if requiredColumn not in metaFrame.columns:
            statusUpdate = "Metadata file is missing column "+requiredColumn+" ...exiting."
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

    #Get a list of fastqs and matching model files from the metadata
    fastqFiles=metaFrame['Fastq'].values
    modelFiles=metaFrame['ReadModel'].values
    insertionFiles=metaFrame['InsertionSequence'].values
    genomeFiles=metaFrame['GenomeSequence'].values
    outputDirs=metaFrame['OutputDir'].values
    shortNames=metaFrame['ShortName'].values
    poolNames=metaFrame['Pool'].values

    for outputN,dir in enumerate(outputDirs):
        if not dir[-1] == "/":
            outputDirs[outputN] = dir+"/"
        if not os.path.exists(outputDirs[outputN]):
            os.makedirs(outputDirs[outputN])

    #Loop through fastq files and map TnSeq reads to the genome
    if options.useMappedFiles:
        statusUpdate = '--useMappedFiles option passed. Skipping mapping step and proceeding to pool analysis'
        printUpdate(options.logFile,statusUpdate)
    else:
        statusUpdate = 'Finding barcodes in fastqs and mapping insertion locations'
        printUpdate(options.logFile,statusUpdate)
        if options.noBarcodes:
            statusUpdate = '  --noBarcodes option passed.  Using pseudo-barcodes from expected genomic sequence region of reads.'
            printUpdate(options.logFile,statusUpdate)
        
    
        for fastqNum,fastqFile in enumerate(fastqFiles):
            statusUpdate = '  Mapping reads from '+fastqFile+' using insertion model '+modelFiles[fastqNum]
            printUpdate(options.logFile,statusUpdate)

            #Load and parse modelfile
            fileToOpen = modelFiles[fastqNum]
            try:
                with open(fileToOpen, 'r') as FileHandle:
                    insertionModel = FileHandle.readlines()
                    FileHandle.close()
            except IOError:
                statusUpdate = " Could not read file:"+fileToOpen+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            for jdx,line in enumerate(insertionModel):
                insertionModel[jdx] = line.upper()
                if not bool(re.match('^[AGCTN]+$', insertionModel[jdx])):
                    statusUpdate = " Model file "+modelFiles[jdx]+" contains illegal characters (only AGCTN allowed) on line "+str(fastqNum+1)+", exiting..."
                    printUpdate(options.logFile,statusUpdate)
                    sys.exit()

            if len(insertionModel) < 4:
                statusUpdate = " Model file "+modelFiles[fastqNum]+" is too short, exiting..."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            skipBases = len(insertionModel[0].strip())
            beforeBarcode = insertionModel[1].strip()
            expectedBarcodeLength = len(insertionModel[2].strip())
            afterBarcode = insertionModel[3].strip()

            if len(beforeBarcode) < options.matchBefore:
                statusUpdate = " Model file defines a shorter sequence before the barcode than required by --matchBefore, exiting"
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            if len(afterBarcode) < max(options.matchAfter,options.matchJunction):
                statusUpdate = " Model file defines a shorter sequence after the barcode than required by --matchAfter or --matchJunction, exiting"
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            #Define the short sequences to search for flanking barcodes
            searchBefore = beforeBarcode[-options.matchBefore:]
            searchAfter = afterBarcode[:options.matchAfter]

            #Find the expected positions for search sequences and genome junction
            offsetBefore = skipBases + len(beforeBarcode) - len(searchBefore)
            offsetBarcode = skipBases + len(beforeBarcode)
            offsetAfter = offsetBarcode + expectedBarcodeLength
            offsetGenome = offsetAfter + len(afterBarcode)

            statusUpdate = "  Looking for "+str(expectedBarcodeLength)+"bp sequence barcode flanked by "+searchBefore+" and "+searchAfter
            printUpdate(options.logFile,statusUpdate)

            #Create output files for mapped reads and reads without barcodes
            try:
                blastQueryFileName = outputDirs[fastqNum] + shortNames[fastqNum] + '_blastquery.fasta'
                with open(blastQueryFileName, 'w') as blastQueryFile:
                    blastQueryFile.write("")
                    blastQueryFile.close()
            except IOError:
                statusUpdate = "  Could not read file:"+blastQueryFileName+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()
            try:
                noBarcodeFoundFileName = outputDirs[fastqNum] + shortNames[fastqNum] + '_noBarcodeFound.txt'
                with open(noBarcodeFoundFileName, 'w') as noBarcodeFile:
                    noBarcodeFile.write("ReadId\tReadSeq\n")
                    noBarcodeFile.close()
            except IOError:
                statusUpdate = "  Could not read file:"+noBarcodeFoundFileName+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            #Load fastq and filter for reads with the expected barcode sequence
            fileToOpen = fastqFile
            try:
                with open(fileToOpen, 'r') as FileHandle:
                    readName=""
                    readCount=0
                    readSeq=""
                    readQuql=""
                    NbarcodeFound = 0
                    NbarcodeNotFound = 0
                    lowQualScores = {}
                    for readLine in FileHandle:
                        if readLine[0] == "@":
                            readName = readLine.strip()
                            readCount = 0
                        elif readCount == 1:
                            readSeq = readLine.strip()

                        #Third line after a readname is the quality score, now weve loaded info for read and can filter for barcode presence before moving to next read.
                        elif readCount == 3:
                            readQual = readLine.strip()
                            #First look for barcode with flanking sequences
                            passedBefore, deviationBefore = matchExpected(searchBefore,readSeq,offsetBefore,options.wobbleAllowed)
                            passedAfter = False
                            foundEnd = False
                            endOfInsertion = 0
                            barcode=""
                            compliantBarcode = False
                            if passedBefore:
                                passedAfter, deviationAfter = matchExpected(searchAfter,readSeq,offsetAfter+deviationBefore,options.wobbleAllowed)
                                if passedAfter:
                                    barcode = readSeq[offsetBarcode+deviationBefore:offsetBarcode+expectedBarcodeLength+deviationBefore+deviationAfter]
                                    compliantBarcode = bool(re.match('^[AGCT]+$', barcode)) and (abs(len(barcode)-expectedBarcodeLength) <= options.barcodeVariation)
                                    barQual = readQual[offsetBarcode+deviationBefore:offsetBarcode+expectedBarcodeLength+deviationBefore+deviationAfter]
                                    if min(map(ord,list(barQual)))-33 < options.minQual:
                                        compliantBarcode = False

                            #If the barcode is found, try to find the end of the inserted sequences
                            if compliantBarcode:
                                for scanPosition in list(range(len(afterBarcode)-options.matchJunction,-1,-1)):
                                    if not foundEnd:
                                        searchJunction = afterBarcode[scanPosition:scanPosition+options.matchJunction]
                                        foundEnd, deviation = matchExpected(searchJunction,readSeq,offsetAfter+deviationBefore+deviationAfter+scanPosition,options.wobbleAllowed)
                                        if foundEnd:
                                            endOfInsertion = offsetAfter + deviationBefore + deviationAfter + scanPosition + options.matchJunction

                            #Use the end of sequence matching the insertion model as the junction point with the genome.
                            genomeQuery = ""
                            if foundEnd:
                                NbarcodeFound+=1
                                genomeQuery = readSeq[endOfInsertion:]
                                if options.noBarcodes:
                                    barcode = genomeQuery[0:20]
                                try:
                                    with open(blastQueryFileName, 'a') as blastQueryFile:
                                        blastQueryFile.write(">" + barcode + "__"  + readName + "\n")
                                        blastQueryFile.write(genomeQuery+"\n")
                                        blastQueryFile.close()
                                except IOError:
                                    statusUpdate = "  Could not read file:"+blastQueryFileName+" ...exiting."
                                    printUpdate(options.logFile,statusUpdate)
                                    sys.exit()



                            #Save reads without barcodes for troubleshooting TnSeq
                            else:
                                NbarcodeNotFound+=1
                                try:
                                    with open(noBarcodeFoundFileName, 'a') as noBarcodeFile:
                                        noBarcodeFile.write(readName+"\t"+readSeq+"\n")
                                        noBarcodeFile.close()
                                except IOError:
                                    statusUpdate = "  Could not read file:"+noBarcodeFoundFileName+" ...exiting."
                                    printUpdate(options.logFile,statusUpdate)
                                    sys.exit()

                        readCount+=1

                    FileHandle.close()

            except IOError:
                statusUpdate = "  Could not read file:"+fileToOpen+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

            statusUpdate = "    " + str(NbarcodeNotFound) + " reads without recognizable, compliant barcodes written to " + noBarcodeFoundFileName
            printUpdate(options.logFile,statusUpdate)
            statusUpdate = "    " + str(NbarcodeFound) + " query sequences with compliant barcodes written to " + blastQueryFileName
            printUpdate(options.logFile,statusUpdate)

            #Blast query sequences against genome
            statusUpdate = "  BLASTING reads against "+genomeFiles[fastqNum]
            printUpdate(options.logFile,statusUpdate)
            blastGenome = NcbiblastnCommandline(query=blastQueryFileName, db=genomeFiles[fastqNum], perc_identity=options.minPercentID, evalue=options.maxEvalue, outfmt=6, out=outputDirs[fastqNum]+shortNames[fastqNum]+"_blastGenome.txt")
            BlastOutputFile = outputDirs[fastqNum]+shortNames[fastqNum]+"_blastGenome.txt"
            stdout, stderr = blastGenome()

            mappedReads = {}
            try:
                BlastOutputFile = outputDirs[fastqNum]+shortNames[fastqNum]+"_blastGenome.txt"
                with open(BlastOutputFile, 'r') as BlastOutput:
                    for outputLine in BlastOutput:
                        blastArray = outputLine.strip().split("\t")
                        barcode,readName = blastArray[0].split("__")
                        scaffold = blastArray[1]
                        percentID = blastArray[2]
                        matchLength = int(blastArray[3])
                        qstart = int(blastArray[6])
                        qend = int(blastArray[7])
                        sstart = int(blastArray[8])
                        send = int(blastArray[9])
                        bitscore = float(blastArray[11])
                        if sstart < send:
                            sign = '+'
                        else:
                            sign = '-'

                        if (readName in mappedReads):
                            if mappedReads[readName]['loc2'] == 'None':
                                #possibly unnecessary filter for repeated lines in BLAST ouput
                                if  not (mappedReads[readName]['scaffold'] == scaffold and abs(mappedReads[readName]['loc'] - sstart) < 10):
                                    mappedReads[readName]['scaffold2']=scaffold
                                    mappedReads[readName]['loc2']=sstart
                                    mappedReads[readName]['strand2']=sign
                                    mappedReads[readName]['bitscore2']=bitscore
                                    if bitscore > (mappedReads[readName]['bitscore'] - options.scoreDiff):
                                        mappedReads[readName]['unique'] = 0
                        else:
                            mappedReads[readName] = {'barcode':barcode,'scaffold':scaffold,'loc':sstart,'strand':sign,'unique':1,
                                                    'qstart':qstart,'qend':qend,'bitscore':bitscore,'percentID':percentID,
                                                    'scaffold2':'None','loc2':'None','strand2':'None','bitscore2':'None',
                                                    'InsertPlasmid':'None','InsLoc':'None','InsStrand':'None','InsScore':'None'}
                    BlastOutput.close()

                    statusUpdate = "    " + str(len(mappedReads)) + " Reads mapped to genome."
                    printUpdate(options.logFile,statusUpdate)

            except IOError:
                statusUpdate = " Could not read file:"+BlastOutputFile+"...exiting."
                printUpdate(options.logFile,statusUpdate)
                
                sys.exit()


            #Blast query sequences against rest of insert
            #All insert hits scored as 'unique'. Only care about reads with ambiguous mapping to different places in the genome.
            if not options.noInsertHits:
                statusUpdate = "  BLASTING reads against full insertion sequence in "+insertionFiles[fastqNum]
                printUpdate(options.logFile,statusUpdate)
                blastInsert = NcbiblastnCommandline(query=blastQueryFileName, db=insertionFiles[fastqNum], perc_identity=options.minPercentID, evalue=options.maxEvalue, outfmt=6, out=outputDirs[fastqNum]+shortNames[fastqNum]+"_blastInsert.txt")
                stdout, stderr = blastInsert()

                newReadsMapped = 0
                try:
                    BlastOutputFile = outputDirs[fastqNum]+shortNames[fastqNum]+"_blastInsert.txt"
                    with open(BlastOutputFile, 'r') as BlastOutput:
                        for outputLine in BlastOutput:
                            blastArray = outputLine.strip().split("\t")
                            barcode,readName = blastArray[0].split("__")
                            scaffold = blastArray[1]
                            percentID = blastArray[2]
                            matchLength = int(blastArray[3])
                            qstart = int(blastArray[6])
                            qend = int(blastArray[7])
                            sstart = int(blastArray[8])
                            send = int(blastArray[9])
                            bitscore = float(blastArray[11])
                            if sstart < send:
                                sign = '+'
                            else:
                                sign = '-'

                            if (readName in mappedReads):
                                if mappedReads[readName]['InsLoc'] == 'None':
                                    mappedReads[readName]['InsertPlasmid']=scaffold
                                    mappedReads[readName]['InsLoc']=sstart
                                    mappedReads[readName]['InsStrand']=sign
                                    mappedReads[readName]['InsScore']=bitscore

                            else:
                                newReadsMapped+=1
                                mappedReads[readName] = {'barcode':barcode,'scaffold':'None','loc':'None','strand':'None','unique':1,
                                                    'qstart':'None','qend':'None','bitscore':'None','percentID':'None',
                                                    'scaffold2':'None','loc2':'None','strand2':'None','bitscore2':'None',
                                                    'InsertPlasmid':scaffold,'InsLoc':sstart,'InsStrand':sign,'InsScore':bitscore}
                        BlastOutput.close()

                        statusUpdate = "    "+str(newReadsMapped) + " Reads mapped to insertion/plasmid sequence."
                        printUpdate(options.logFile,statusUpdate)

                except IOError:
                    statusUpdate = " Could not read file:"+BlastOutputFile+"...exiting."
                    printUpdate(options.logFile,statusUpdate)
                    sys.exit()


            try:
                mappedReadsFileName = outputDirs[fastqNum] + shortNames[fastqNum] + '_mapped.txt'
                with open(mappedReadsFileName, 'w') as mappedReadsFile:
                    mappedReadsFile.write("ReadName\tBarcode\tScaffold\tLocation\tStrand\tUniquelyMapping\tQueryStart\tQueryEnd\tBitscore\tPercentID")
                    mappedReadsFile.write("\tScaffold2ndHit\tLocation2ndHit\tStrand2ndHit\tBitscore2ndHit\tInsertOrPlasmidHit\tInsertHitLocation\tInsertHitStrand\tInsertHitScore\n")
                    for readName in mappedReads:
                        linetowrite = [readName]
                        for item in ['barcode','scaffold','loc','strand','unique','qstart','qend','bitscore','percentID',
                                     'scaffold2','loc2','strand2','bitscore2','InsertPlasmid','InsLoc','InsStrand','InsScore']:
                            linetowrite.append(str(mappedReads[readName][item]))
                        mappedReadsFile.write('\t'.join(linetowrite)+'\n')
                    mappedReadsFile.close()

            except IOError:
                statusUpdate = "Could not read file:"+mappedReadsFileName+"...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()


    statusUpdate = "Processing mapped reads to characterize mutant pool(s)"
    printUpdate(options.logFile,statusUpdate)

    mappedFilesPerPool={}
    poolOutPutDirs={}
    for mappedNum,pool in enumerate(poolNames):
        mappedFileName = outputDirs[mappedNum]+shortNames[mappedNum]+"_mapped.txt"
        if not pool in mappedFilesPerPool:
            mappedFilesPerPool[pool] = [mappedFileName]
            poolOutPutDirs[pool] = outputDirs[mappedNum]
        else:
            mappedFilesPerPool[pool].append(mappedFileName)
        
    for pool in mappedFilesPerPool:
        statusUpdate = "Analyzing pool: " + pool
        printUpdate(options.logFile,statusUpdate)

        totals={}
        barcode_to_genome = {}
        barcode_to_insert = {}
        barcode_flags = {}
        barcode_primary_location = {}
        barcode_summaries = {}
        offbytwo = []

        #Load up mapped inserts from one or more files
        minFraction = float(options.minFraction)
        maxQbeg = int(options.maxFillerSeq)

        totalLines = 0
        nFiles = 0
        for tabFileName in mappedFilesPerPool[pool]:
            nLines = 0
            statusUpdate = "  Reading entries from: " + tabFileName
            printUpdate(options.logFile,statusUpdate)
            try:
                with open(tabFileName, 'r') as tabFileHandle:
                    first_line = tabFileHandle.readline()
                    for line in tabFileHandle:
                        nLines += 1
                        if (nLines % 1000000 == 0 ):
                            statusUpdate = "      " + str(nLines/1000000) + " million lines read... "
                            printUpdate(options.logFile,statusUpdate)
 
                        readMapped = line.rstrip().split()
                        location = readMapped[2]+":"+readMapped[3]+":"+readMapped[4]

                        #First check if barcode is new or has been seen before
                        newMappedBarcode = readMapped[1]
                        if "N" not in newMappedBarcode:
                            if newMappedBarcode in totals:
                                    totals[newMappedBarcode] += 1
                                    if readMapped[2] == "None":
                                        location = readMapped[14]+":"+readMapped[15]+":"+readMapped[16]
                                        barcode_flags[newMappedBarcode]['insTot'] += 1
                                        if location in barcode_to_insert[newMappedBarcode]:
                                            barcode_to_insert[newMappedBarcode][location] += 1
                                        else:
                                            barcode_to_insert[newMappedBarcode][location] = 1
                                    else:
                                        if location in barcode_to_genome[newMappedBarcode]:
                                            barcode_to_genome[newMappedBarcode][location] += 1
                                        else:
                                            barcode_to_genome[newMappedBarcode][location] = 1
                            else:
                                    totals[newMappedBarcode] = 1
                                    if readMapped[2] == "None":
                                        location = readMapped[14]+":"+readMapped[15]+":"+readMapped[16]
                                        barcode_to_insert[newMappedBarcode] = {location:1}
                                        barcode_to_genome[newMappedBarcode] = {'Null:0:+':0}
                                        barcode_flags[newMappedBarcode] = {'ambReads':0,'highQbeg':0,'insTot':1}
                                    else:
                                        barcode_to_genome[newMappedBarcode] = {location:1}
                                        barcode_to_insert[newMappedBarcode] = {'Null:0:+':0}
                                        barcode_flags[newMappedBarcode] = {'ambReads':0,'highQbeg':0,'insTot':0}

                            if not readMapped[2] == "None" and int(readMapped[6]) > maxQbeg:
                                barcode_flags[newMappedBarcode]['highQbeg'] += 1

                            if (int(readMapped[5]) != 1):
                                barcode_flags[newMappedBarcode]['ambReads'] += 1

                    tabFileHandle.close()
                    totalLines = totalLines + nLines
                    nFiles += 1

            except IOError:
                statusUpdate = "Could not read file:"+tabFileName+"...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

        
        statusUpdate = "    "+str(totalLines)+" lines read from "+str(nFiles)+" files"
        printUpdate(options.logFile,statusUpdate)


        # Mask out off-by-one variants
        statusUpdate = "  Masking out off-by-one barcodes. (Most likely sequencing errors)"
        printUpdate(options.logFile,statusUpdate)

        offByOneList = []
        for barcode in totals:
            variants = OffByOneList(barcode)
            offByOne = False
            for variantBarcode in variants:
                if (not offByOne) & (variantBarcode in totals):
                    if (totals[variantBarcode] > totals[barcode]*100):
                        offByOne = True
                        offByOneList.append(barcode)

        outPutFileName = poolOutPutDirs[pool] + pool + "_poolfile"
        offByOneFileHandle = open(outPutFileName+"_offByOne", 'w')
        barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n"])
        offByOneFileHandle.write(barcodeSummary+"\n")
        nOffByOne = 0
        OffByOneReads = 0
        for barcode in offByOneList:
            OffByOneReads += totals[barcode]
            barcodeSummary = "\t".join([barcode,ReverseComplement(barcode),str(totals[barcode]), str(totals[barcode]-barcode_flags[barcode]['highQbeg'])])
            offByOneFileHandle.write(barcodeSummary+"\n")
            del barcode_to_insert[barcode]
            del barcode_to_genome[barcode]
            del totals[barcode]
            del barcode_flags[barcode]
            nOffByOne +=1
        offByOneFileHandle.close()

        statusUpdate = "    Wrote " + str(nOffByOne) + " off-by-one barcodes ("+str(OffByOneReads)+" reads) to "+outPutFileName+"_offByOne"
        printUpdate(options.logFile,statusUpdate)

        #Set up pool file
        
        outPutFileHandle = open(outPutFileName, 'w')
        ambFileHandle = open(outPutFileName+"_ambiguous", 'w')
        multiFileHandle = open(outPutFileName+"_multiLocus", 'w')
        barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n", "scaffold", "strand", "pos","type","nMainLocation", "nInsert", "All genomic mappings", "All insert mappings"])
        outPutFileHandle.write(barcodeSummary+"\n")
        ambFileHandle.write(barcodeSummary+"\n")
        multiFileHandle.write(barcodeSummary+"\n")

        statusUpdate = "  Processing barcodes with mapped genomic locations"
        printUpdate(options.logFile,statusUpdate)

        nMulti = 0
        nAmb = 0
        nInsert = 0
        nMapped = 0
        MultiReads = 0
        AmbReads = 0
        allInsertReads = 0
        MappedReads = 0

        for barcode, barcode_locations in barcode_to_genome.items():

            sorted_locations = sorted(barcode_to_genome[barcode].items(), key=operator.itemgetter(1), reverse=True)
            firstLocation = sorted_locations[0][0].split(":")

            ambiguousReads = int(barcode_flags[barcode]['ambReads'])
            insertReads = int(barcode_flags[barcode]['insTot'])
            highQbegreads = int(barcode_flags[barcode]['highQbeg'])
            insertionType = ""
            primeLocationAbundance = 0

            #Ignore barcodes with vast majority of reads in insert
            if (insertReads > (sorted_locations[0][1] * 7)):
                insertionType = "insert"
                primeLocationAbundance = insertReads

            #If ambiguous reads are more common than top unambiguous location tag, barcode with ambiguous mapping
            elif (sorted_locations[0][1] < 2 * ambiguousReads):
                insertionType = "Amb"
                primeLocationAbundance = ambiguousReads

            #If not insert, ambiguous, or multilocus hit, classify the type of insertion
            else:
                #If second most common location is on opposite strand nearby, call as inverted concatamer
                if (len(sorted_locations) > 1) and (sorted_locations[1][0] != 'Null'):
                    secondLocation = sorted_locations[1][0].split(":")
                    if (firstLocation[0] == secondLocation[0]) & (firstLocation[2] != secondLocation[2]) & (abs(int(firstLocation[1])-int(secondLocation[1])) < 10):
                        insertionType = "InCat"
                        primeLocationAbundance = sorted_locations[0][1]+sorted_locations[1][1]
                    else:
                        insertionType = "Single"
                        primeLocationAbundance = sorted_locations[0][1]
                else:
                    insertionType = "Single"
                    primeLocationAbundance = sorted_locations[0][1]

                #If not inverted concatamer, but both insert hits and genome hits found, tag as concatamer
                if (insertReads > (sorted_locations[0][1])/10) & (insertionType != "InCat"):
                    insertionType = "Concat"
                    primeLocationAbundance = sorted_locations[0][1]

                #Most common location must accout for at least minFrac of reads (not counting insert hits)
                if ((float(primeLocationAbundance)/float(totals[barcode]-insertReads)) < minFraction):
                    insertionType = "Multi"
                    primeLocationAbundance = sorted_locations[0][1]

            #Print out summary
            barcodeSummary = "\t".join([barcode,ReverseComplement(barcode),str(totals[barcode]), str(totals[barcode]-highQbegreads), firstLocation[0], firstLocation[2], firstLocation[1],insertionType,str(primeLocationAbundance), str(insertReads), str(sorted_locations), str(barcode_to_insert[barcode])])

            if insertionType in ["Single","Concat","InCat"]:
                outPutFileHandle.write(barcodeSummary+"\n")
                nMapped += 1
                MappedReads += totals[barcode]
                barcode_primary_location[barcode] = {'scaffold':firstLocation[0],'strand':firstLocation[2],'position':int(firstLocation[1]),'counts':int(totals[barcode])}
                barcode_summaries[barcode] = barcodeSummary
            elif insertionType == "Multi":
                multiFileHandle.write(barcodeSummary+"\n")
                nMulti += 1
                MultiReads += totals[barcode]
            elif insertionType == "Amb":
                ambFileHandle.write(barcodeSummary+"\n")
                nAmb += 1
                AmbReads += totals[barcode]

            #Remove barcode from insert list to prevent redundant reporting
            if insertionType != "insert":
                del barcode_to_insert[barcode]

        statusUpdate = "    Wrote "+str(nMapped)+" barcodes ("+str(MappedReads)+" reads) to "+outPutFileName
        printUpdate(options.logFile,statusUpdate)
        outPutFileHandle.close()

        statusUpdate = "    Wrote "+str(nMulti)+" barcodes ("+str(MultiReads)+" reads) to "+outPutFileName+"_multiLocus"
        printUpdate(options.logFile,statusUpdate)
        multiFileHandle.close()

        statusUpdate = "    Wrote "+str(nAmb)+" barcodes ("+str(AmbReads)+" reads) to "+outPutFileName+"_ambiguous"
        printUpdate(options.logFile,statusUpdate)
        ambFileHandle.close()

        #Set up pool file for insert hits
        outPutFileHandle = open(outPutFileName+"_insertReadsOnly", 'w')
        barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n", "scaffold", "strand", "pos","type","nHeadtoTail", "nHeadtoHead", "nUnprocessed", "Misprocessed","All Insert Mappings","All Genomic Mappings"])
        outPutFileHandle.write(barcodeSummary+"\n")

        statusUpdate = "  Processing hits to insert only"
        printUpdate(options.logFile,statusUpdate)

        #Report barcodes with only insert reads
        for barcode in barcode_to_insert:
            HeadtoTail = 0
            HeadtoHead = 0
            Unprocessed = 0
            Misprocessed = 0

            sorted_locations = sorted(barcode_to_insert[barcode].items(), key=operator.itemgetter(1), reverse=True)
            firstLocation = sorted_locations[0][0].split(":")
            insertReads = barcode_flags[barcode]['insTot']

            #Classify the type of insertion
            for insertLocation in sorted_locations:
                insertLocationDetails = insertLocation[0].split(":")

                if (insertLocationDetails[2] == '-') & (insertLocationDetails[0].upper() == 'INSERT') :
                    HeadtoHead += int(insertLocation[1])
                elif (insertLocationDetails[2] == '+') & (insertLocationDetails[0].upper() == 'INSERT'):
                    HeadtoTail += int(insertLocation[1])
                elif (insertLocationDetails[2] == '+') &  (insertLocationDetails[0].upper() == 'PLASMID') & (int(insertLocationDetails[1]) < 10):
                    Unprocessed += int(insertLocation[1])
                else:
                    Misprocessed += int(insertLocation[1])


            if HeadtoHead > insertReads * 0.9:
                insertionType = "HeadtoHead"
            elif Unprocessed > insertReads * 0.9 :
                insertionType = "Unprocessed"
            elif HeadtoTail > insertReads * 0.9:
                insertionType = "HeadtoTail"
            elif Misprocessed > insertReads * 0.9:
                insertionType = "Misprocessed"
            else:
                insertionType = "Mixed"


            #Print out summary
            barcodeSummary = "\t".join([barcode,ReverseComplement(barcode),str(insertReads), str(insertReads), "insert", firstLocation[2], firstLocation[1],insertionType,str(HeadtoTail), str(HeadtoHead), str(Unprocessed),str(Misprocessed), str(sorted_locations), str(barcode_to_genome[barcode])])
            outPutFileHandle.write(barcodeSummary+"\n")
            nInsert += 1
            allInsertReads += totals[barcode]
        outPutFileHandle.close()

        statusUpdate = "    Wrote "+str(nInsert)+" barcodes ("+str(allInsertReads)+" reads) to "+outPutFileName+"_insertHitsOnly"
        printUpdate(options.logFile,statusUpdate)


        if options.filterNeighborhood > 0 and len(barcode_primary_location) > 0:
            #Find similair sequences (off-by-two, short indels) at same position (computationaly intensive).
            statusUpdate = "  Finding similair barcodes that mapped to the same location (off-by-two, indels, etc)."
            printUpdate(options.logFile,statusUpdate)

            #Loop through scaffolds
            locationFrame = pd.DataFrame.from_dict(barcode_primary_location,orient='index')
            scaffoldGroups = locationFrame.groupby('scaffold')

            for scaffold, barcodes in scaffoldGroups:
                statusUpdate = "    Checking scaffold "+scaffold
                printUpdate(options.logFile,statusUpdate)

                
                #For each scaffold get a list of barcodes sorted by position
                barcodes = barcodes.sort_values(by='position')
                scaffold_barcode_list = barcodes.index.values
                scaffold_barcode_positions = barcodes['position'].values
                scaffold_barcode_counts = barcodes['counts'].values

                #Scroll a moving window through each scaffold checking local neighborhood for similair barcodes
                for idx, ref_position in enumerate(scaffold_barcode_positions[0:-1]):
                    window_end = min(len(scaffold_barcode_positions),idx+20)
                    check_idx = idx+1
                    #only need to compare barcodes that map very close to each other
                    while check_idx < len(scaffold_barcode_positions) and scaffold_barcode_positions[check_idx] - ref_position < int(options.filterNeighborhood):
                        
                        filterTest = False

                        #If the library isn't barcoded, no way to distinguish truly different insertions at the same locus
                        if options.noBarcodes:
                            filterTest = True

                        #Calculate minimum edit distance for nearby by barcodes
                        else:
                            filterTest = (levenshtein(scaffold_barcode_list[idx], scaffold_barcode_list[check_idx]) < int(options.filterEditDistance))
                    
                        if filterTest:
                            #Keep the most abundant barcode/insertion location
                            if scaffold_barcode_counts[idx] > scaffold_barcode_counts[check_idx]:
                                if not scaffold_barcode_list[check_idx] in offbytwo:
                                    offbytwo.append(scaffold_barcode_list[check_idx])
                            else:
                                if not scaffold_barcode_list[idx] in offbytwo:
                                    offbytwo.append(scaffold_barcode_list[idx])
                        check_idx += 1

            outPutFileHandle = open(outPutFileName, 'w')
            obtFileHandle = open(outPutFileName+"_filteredLocal", 'w')
            barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n", "scaffold", "strand", "pos","type","nMainLocation", "nInsert", "All genomic mappings", "All insert mappings"])
            obtFileHandle.write(barcodeSummary+"\n")
            outPutFileHandle.write(barcodeSummary+"\n")
            nGood = 0
            nBad = 0
            for barcode in barcode_summaries:
                summary_arr = barcode_summaries[barcode].split("\t")
                if barcode in offbytwo:
                    summary_arr[7] = "offByTwo"
                    obtFileHandle.write("\t".join(summary_arr)+"\n")
                    nBad+=1
                else:
                    outPutFileHandle.write("\t".join(summary_arr)+"\n")
                    nGood+=1
            outPutFileHandle.close()
            obtFileHandle.close()
            statusUpdate = "    Moved " + str(nBad) + " barcodes with likely sequencing errors to " + outPutFileName + "_filteredLocal"
            printUpdate(options.logFile,statusUpdate)

            statusUpdate = "    Wrote "+str(nGood)+" barcodes ("+str(MappedReads)+" reads) to "+outPutFileName
            printUpdate(options.logFile,statusUpdate)


if __name__ == "__main__":
    main(sys.argv[1:])
exit()
