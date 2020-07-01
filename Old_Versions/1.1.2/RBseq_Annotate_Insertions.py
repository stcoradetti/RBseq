#!/usr/bin/python
import sys
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import operator
from datetime import datetime
import collections
import matplotlib
import matplotlib.pyplot as plt
import argparse

Version = '1.1.2'
ReleaseDate = 'Apr 13, 2020'


plt.style.use('seaborn-colorblind')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42


def printUpdate(logfile,update):
    timestamp = datetime.now().strftime('%Y-%m-%d %X')
    update = timestamp + ' ' + update
    file = open(logfile, "a")
    file.write(update+'\n')
    file.close()
    print(update)
    return update

def main(argv):
    timestamp = datetime.now().strftime('%Y%m%H%M')
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metafile", dest="metafile", help="Metadata file of TnSeq runs.",default="metafile.txt")
    parser.add_argument("-l", "--logFile", dest="logFile", help="File to write run log to. Default is Annotate_TIMESTAMP.log",default="Annotate_"+timestamp+".log")
    parser.add_argument("-P", "--promoterLength", dest="promoterLength", help="Number of bases before transcriptional start counted as promoter sequence",default=500)
    parser.add_argument("-T", "--terminatorLength", dest="terminatorLength", help="Number of bases after transcriptional stop counted as terminator sequence",default=100)

    options  = parser.parse_args()

    statusUpdate = 'RBseq_Annotate_Insertions.py'
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
            metaFrame = pd.read_csv(FileHandle,low_memory=False,sep='\t')
            metaFrame = metaFrame[metaFrame[metaFrame.columns[0]].notna()]
            FileHandle.close()
    except IOError:
        statusUpdate = "Could not read file:"+fileToOpen+" ...exiting."
        printUpdate(options.logFile,statusUpdate)
        sys.exit()

    requiredColumns=['Pool','GenomeSequence','OutputDir','GeneLocations','GeneAnnotations','GeneIdentifier']
    for requiredColumn in requiredColumns:
        if requiredColumn not in metaFrame.columns:
            statusUpdate = "Metadata file is missing column "+requiredColumn+" ...exiting."
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

    poolNames=metaFrame['Pool'].values
    gffNames=metaFrame['GeneLocations'].values
    outputDirs=metaFrame['OutputDir'].values
    genomeFiles=metaFrame['GenomeSequence'].values
    annotationFiles=metaFrame['GeneAnnotations'].values
    GeneIdentifiers=metaFrame['GeneIdentifier'].values

    for outputN,dir in enumerate(outputDirs):
        if not dir[-1] == "/":
            outputDirs[outputN] = dir+"/"
        if not os.path.exists(outputDirs[outputN]):
            os.makedirs(outputDirs[outputN])

    poolFileNames=[]
    gffFileNames=[]
    annotationFileNames=[]
    poolGeneIdentifiers=[]
    
    for lineNum,name in enumerate(poolNames):
        poolFileName = outputDirs[lineNum]+name+"_poolfile"
        if not poolFileName in poolFileNames:
            poolFileNames.append(poolFileName)
            gffFileNames.append(gffNames[lineNum])
            annotationFileNames.append(annotationFiles[lineNum])
            poolGeneIdentifiers.append(GeneIdentifiers[lineNum])



    for poolNum, poolFileName in enumerate(poolFileNames):
        gffFileName = gffFileNames[poolNum]
        annotationFile = annotationFileNames[poolNum]
        GeneIdentifier = poolGeneIdentifiers[poolNum]
        if GeneIdentifier in ['',' ','\t','NA',np.nan]:
            GeneIdentier = 'Parent'
        scaffoldSequences = SeqIO.to_dict(SeqIO.parse(genomeFiles[poolNum]+".fasta", "fasta"))
        
        #Load gene locations
        try:
            with open(gffFileName, 'rb') as FileHandle:
                genestext = FileHandle.readlines()
                FileHandle.close()
                statusUpdate = "Read "+str(len(genestext))+" features from "+gffFileName
                printUpdate(options.logFile,statusUpdate)
        except IOError:
            statusUpdate = "Could not read file: "+gffFileName
            printUpdate(options.logFile,statusUpdate)
            sys.exit()
            
        #Load pool information
        try:
            with open(poolFileName, 'rb') as poolFileHandle:
                poolFrame = pd.read_csv(poolFileHandle,sep='\t')
                poolFileHandle.close()
                poolFrame.dropna(how='all')
                statusUpdate =  "Read "+ str(len(poolFrame)) + " barcodes from "+poolFileName
                printUpdate(options.logFile,statusUpdate)
        except IOError:
            statusUpdate =  "Could not read file: "+poolFileName
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

        #Load gene annotations
        try:
            with open(annotationFile, 'rb') as FileHandle:
                annotFrame = pd.read_csv(FileHandle,sep='\t')
                annotFrame = annotFrame[annotFrame[annotFrame.columns[0]].notna()]
                FileHandle.close()
                statusUpdate = "Read "+str(len(annotFrame))+" annotations from "+annotationFile
                printUpdate(options.logFile,statusUpdate)
        except IOError:
            statusUpdate = "Could not read file: "+annotationFile
            printUpdate(options.logFile,statusUpdate)
            sys.exit()

        requiredColumns=['ID','AlternateID','Annotation']
        for requiredColumn in requiredColumns:
            if requiredColumn not in annotFrame.columns:
                statusUpdate = "Gene annotation file is missing column "+requiredColumn+" ...exiting."
                printUpdate(options.logFile,statusUpdate)
                sys.exit()

        


        #Sort pool by location
        statusUpdate =  "Sorting pool entries by location."
        printUpdate(options.logFile,statusUpdate)
        poolFrame.sort_values(['scaffold','pos'],ascending=[1,1],inplace=True)
        

        scaffoldLengths = {}
        totalGenome = 0
        for scaffold in scaffoldSequences:
            scaffoldLengths[scaffold] = len(scaffoldSequences[scaffold].seq)
            totalGenome+=len(scaffoldSequences[scaffold].seq)

        sorted_scaffoldLengths = sorted(scaffoldLengths.items(), key=operator.itemgetter(1), reverse=True)
        scaffoldOrder = []
        for scaff in sorted_scaffoldLengths:
            scaffoldOrder.append(scaff[0])

        statusUpdate =  "Parsing GFF for gene locations and attributes"
        printUpdate(options.logFile,statusUpdate)
        #Build list of RNAs
        Genes = collections.OrderedDict()
        Counter = 0
        for genesline in genestext:
            genesfeilds = genesline.decode().split('\t')
            if len(genesfeilds) == 9:
                scaffold,evidence,featureType,start,end,flag,strand,frame,comments = genesfeilds
                #translate start/end to zero based indexing with non-inclusive end
                start = int(start)-1
                end = int(end)
                if featureType == 'mRNA':
                    notes = dict((k.strip(), v.strip()) for k,v in
                                  (item.split('=') for item in comments.split(';')))

                    if strand == '+':
                        promoter = [start-options.promoterLength,start-1]
                        terminator = [end+1,end+options.terminatorLength]
                    else:
                        promoter = [end+1,end+options.promoterLength]
                        terminator = [start-options.terminatorLength,start-1]

                    terminator[0] = max(0,terminator[0])
                    promoter[0] = max(0,promoter[0])
                    terminator[1] = min(terminator[1],scaffoldLengths[scaffold]-1)
                    promoter[1] = min(promoter[1],scaffoldLengths[scaffold]-1)
                    
                    Genes[notes['ID']] = {'scaffold':scaffold,
                                                 'start':int(start),
                                                 'end':int(end),
                                                 'strand':strand,
                                                 'exon':[],
                                                 'intron':[],
                                                 'promoter':promoter,
                                                 'terminator':terminator,
                                                 '5pu':[],
                                                 '3pu':[],
                                                 'Parent':notes['Parent'],
                                                 'GeneIdentifier':notes[GeneIdentifier]}
                elif featureType == 'CDS':
                    notes = dict((k.strip(), v.strip()) for k,v in
                                  (item.split('=') for item in comments.split(';')))
                    if len(Genes[notes['Parent']]['exon']) > 0:
                        Genes[notes['Parent']]['intron'].append([Genes[notes['Parent']]['exon'][-1][1]+1,start-1])
                    Genes[notes['Parent']]['exon'].append([start,end])
                elif featureType == 'five_prime_UTR':
                    notes = dict((k.strip(), v.strip()) for k,v in
                                  (item.split('=') for item in comments.split(';')))
                    Genes[notes['Parent']]['5pu'] = [start,end]
                elif featureType == 'three_prime_UTR':
                    notes = dict((k.strip(), v.strip()) for k,v in
                                  (item.split('=') for item in comments.split(';')))
                    Genes[notes['Parent']]['3pu'] = [start,end]


        #Build list of locations in the genome and in genic and non-genic regions
        locationTypes = {}
        nearestGeneDict = {}
        for scaffold in scaffoldOrder:
            locationTypes[scaffold] = ['intergenic']*scaffoldLengths[scaffold]
            nearestGeneDict[scaffold] = ['intergenic']*scaffoldLengths[scaffold]

        for locationType in ['terminator','promoter','3pu','5pu']:
            for Gene in Genes:
                entry = Genes[Gene][locationType]
                if len(entry) == 2:
                    for position in range(entry[0],entry[1]):
                        locationTypes[Genes[Gene]['scaffold']][position] = locationType
                        nearestGeneDict[Genes[Gene]['scaffold']][position] = Gene

        for locationType in ['intron','exon']:
            for Gene in Genes:
                for entry in Genes[Gene][locationType]:
                    for position in range(entry[0],entry[1]):
                        locationTypes[Genes[Gene]['scaffold']][position] = locationType
                        nearestGeneDict[Genes[Gene]['scaffold']][position] = Gene

        #Count up locations in promoters, exons, etc
        locationsByType = {}
        locationTypeTotals = {}
        insertionTypeTotals = {}
        locationTypeClasses = ['terminator','promoter','3pu','5pu','intron','exon','intergenic']
        for locationType in locationTypeClasses:
            locationTypeTotals[locationType] = 0
            insertionTypeTotals[locationType] = 0
            locationsByType[locationType] = {}
            for scaffold in scaffoldOrder:
                locationsByType[locationType][scaffold] = []
                for location in range(0,scaffoldLengths[scaffold]):
                    if locationTypes[scaffold][location] == locationType:
                        locationsByType[locationType][scaffold].append(location)
                        locationTypeTotals[locationType]+=1

        #Count up insertions in promoters, exons, etc.
        codingFractions = []
        nearestGeneIDs = []
        insertionTypesOrdered = []
        statusUpdate =  "Classifying insertion locations"
        printUpdate(options.logFile,statusUpdate)
        for index in poolFrame.index:
            scaffold = poolFrame.at[index,'scaffold']
            position = int(poolFrame.at[index,'pos'])-1
            insertionType = locationTypes[scaffold][position]
            insertionTypeTotals[insertionType] += 1
            insertionTypesOrdered.append(insertionType)

            nearestGene = nearestGeneDict[scaffold][position]
            if nearestGene == 'intergenic':
                nearestGeneIDs.append('intergenic')
            else:
                nearestGeneIDs.append(Genes[nearestGene]['GeneIdentifier'])
            if insertionType in ['intron','exon']:
                flattened_exons = [item for sublist in Genes[nearestGene]['exon'] for item in sublist]
                codingStart = min(flattened_exons)
                codingStop = max(flattened_exons)
                codingFractions.append(int((float(position - codingStart) / (codingStop - codingStart)) * 100))

                
            else:
                codingFractions.append("NA")


        poolFrame['InsertionType'] = insertionTypesOrdered
        poolFrame['NearestGene'] = nearestGeneIDs
        poolFrame['CodingFraction'] = codingFractions

        annotFrame['ID'] = annotFrame['ID'].astype(str)

        #print(annotFrame,poolFrame)

        poolFrame = poolFrame.merge(annotFrame, left_on='NearestGene', right_on='ID',how='left')

        poolFrame.to_csv(poolFileName+'_annotated',sep="\t",index=False)
        statusUpdate =  "Wrote poolfile with gene information to "+ poolFileName+'_annotated'
        printUpdate(options.logFile,statusUpdate)

        file = open(poolFileName+"_insertionDistribution.txt", "w")
        file.write("Type\tInsertions\tGenome\n")
        for locationType in locationTypeClasses:
            file.write(locationType+"\t"+str(insertionTypeTotals[locationType])+"\t"+str(locationTypeTotals[locationType])+"\n")
        file.write("\nTotal Genome: "+str(totalGenome))
        file.close()

        statusUpdate =  "Wrote summary of insertion rates across gene features to "+poolFileName+"_insertionDistribution.txt"
        printUpdate(options.logFile,statusUpdate)


        #Make chart of insertion rates in different features
        plotLabels =[]
        genomeCounts = []
        genomePercents = []
        insertionCounts = []
        insertionPercents = []
        for locationType,total in sorted(locationTypeTotals.items(), key=operator.itemgetter(1)):
            plotLabels.append(locationType)
            genomeCounts.append(total)
            insertionCounts.append(insertionTypeTotals[locationType])

        genomePercents = np.array(genomeCounts)/float(np.sum(genomeCounts))*100
        insertionPercents = np.array(insertionCounts)/float(np.sum(insertionCounts))*100

        fig, ax = plt.subplots()
        fig.set_size_inches(3.42,2)
        barWidth=0.4
        plt.bar(np.arange(len(plotLabels))-barWidth/2,genomePercents,barWidth,label='In Genome')
        plt.bar(np.arange(len(plotLabels))+barWidth/2,insertionPercents,barWidth,color='C2',label='Insertions')
        plt.xticks(np.arange(len(plotLabels)),plotLabels,fontsize=8,rotation=20)
        plt.yticks(fontsize=8)
        plt.ylabel('Fraction (%)',fontsize=8,labelpad=1)
        plt.title('Relative Insertion Frequency by Gene Feature',fontsize=10)
        plt.legend(fontsize=8)
        plt.gcf().subplots_adjust(bottom=0.2, left=0.11, right=0.98, top=0.90)
        plt.savefig(poolFileName+'_insertionDistribution.pdf')
        plt.close()

        statusUpdate =  "Saved summary chart of insertion rates across gene features to "+poolFileName+"_insertionDistribution.pdf"
        printUpdate(options.logFile,statusUpdate)
        

        #Check insertion rates accross each scaffold
        insertionDict = {}
        Ninserts = []
        Lengths = []
        for scaffold in scaffoldOrder:
            insertionDict[scaffold] = poolFrame[poolFrame['scaffold'] == scaffold]['pos'].values
            Lengths.append(scaffoldLengths[scaffold])
            Ninserts.append(len(insertionDict[scaffold]))

        fig, ax = plt.subplots()
        fig.set_size_inches(3.42,2)

        lengthsMb = np.array(Lengths)/1000000.0

        plt.scatter(lengthsMb,Ninserts,s=10)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.xlabel('Scaffold length (Mb)',fontsize=8,labelpad=1)
        plt.ylabel('Insertions',fontsize=8,labelpad=1)
        plt.title('Inserts per scaffold',fontsize=10)
        plt.gcf().subplots_adjust(bottom=0.16, left=0.13, right=0.98, top=0.90)
        plt.savefig(poolFileName+'_InsertsPerScaffold.pdf')
        plt.close()

        statusUpdate =  "Saved summary chart of insertion rates across scaffolds "+poolFileName+"_InsertsPerScaffold.pdf"
        printUpdate(options.logFile,statusUpdate)

        #Check GC content around insertion sites
        scaffoldArr = poolFrame['scaffold'].values
        posArr = poolFrame['pos'].values
        flankingWindow = 50
        GCpcts = []
        scaffGC = []
        for idx, scaffold in enumerate(scaffoldArr):
            pos = int(posArr[idx])
            limits = [max(0,pos-flankingWindow),min(pos+flankingWindow,scaffoldLengths[scaffold])]
            localSeq = str(scaffoldSequences[scaffold].seq[limits[0]:limits[1]])
            Ccount = localSeq.count('C')
            Gcount = localSeq.count('G')
            GCpcts.append(100*(Gcount + Ccount)/len(localSeq))
        poolFrame['LocalGCpercent'] = GCpcts

        poolFrame.to_csv(poolFileName+'_annotated',sep="\t",index=False)
        statusUpdate =  "Wrote poolfile with gene information to "+ poolFileName+'_annotated'
        printUpdate(options.logFile,statusUpdate)

        for scaffold in scaffoldOrder:
            randomLocations = np.random.randint(0, high=scaffoldLengths[scaffold], size=len(insertionDict[scaffold]))
            for pos in randomLocations:
                limits = [max(0,pos-flankingWindow),min(pos+flankingWindow,scaffoldLengths[scaffold])]
                localSeq = str(scaffoldSequences[scaffold].seq[limits[0]:limits[1]])
                Ccount = localSeq.count('C')
                Gcount = localSeq.count('G')
                scaffGC.append(100*(Gcount + Ccount)/len(localSeq))


        fig, ax = plt.subplots()
        fig.set_size_inches(3.42,2)

        bins = np.arange(25,75,1)
        plt.hist(scaffGC,bins, density=True,label='Random sites',color='C0',alpha=0.5)
        plt.hist(GCpcts,bins, density=True,label='Insertion sites',color='C2',alpha=0.5)

        plt.xticks([30,40,50,60,70],fontsize=8)
        plt.xlabel('GC content (%)',fontsize=8,labelpad=1)
        plt.yticks([0,0.05,0.1],[0,5,10],fontsize=8)
        plt.ylabel('Fraction (%)',fontsize=8,labelpad=1)
        plt.title('GC content around insertion sites',fontsize=10)
        legend = plt.legend(loc='upper left',edgecolor=None,frameon=False, fontsize=8)
        plt.gcf().subplots_adjust(bottom=0.16, left=0.11, right=0.98, top=0.90)
        plt.savefig(poolFileName+'_GChistogram.pdf')
        plt.close()

        statusUpdate =  "Saved histogram of insertion location GC content to "+poolFileName+"_GChistogram.pdf"
        printUpdate(options.logFile,statusUpdate)

        #Simulate random insertins with observed biases towardsd promoters/exons etc
        statusUpdate =  "Simulating random insertions across the largest scaffold with observed biases towards promoters,exons, etc"
        printUpdate(options.logFile,statusUpdate)
        scaffold = scaffoldOrder[0]
        biasedInsertions = []
        NscaffoldInserts = len(insertionDict[scaffold])
        TotalInserts = len(poolFrame)
        
        randomInsertions = np.random.randint(1, high=scaffoldLengths[scaffold]-1, size=NscaffoldInserts)
        for locationType in locationTypeClasses:
            classIndexes = np.random.randint(0, high=len(locationsByType[locationType][scaffold])-1, size=int(insertionTypeTotals[locationType] * NscaffoldInserts/float(TotalInserts)))
            for i in classIndexes:
                biasedInsertions.append(locationsByType[locationType][scaffold][i])

        #Taking a rolling window along first major chromosome and find the number of insertions
        window = 500

        
        
        hotspotlocations = {}
        hotspotlocations = []

        scaffold_length = scaffoldLengths[scaffold]
        positions = insertionDict[scaffold]

        biasedPositions = np.array(biasedInsertions)
        locations = np.array(range(scaffold_length+window))
        scaffoldDensity = np.zeros(scaffold_length+window)
        randomDensity = np.zeros(scaffold_length+window)
        biasedDensity = np.zeros(scaffold_length+window)

        statusUpdate =  "Summing insertions in the largest scaffold in a rolling window of "+str(2*window)+" base pairs."
        printUpdate(options.logFile,statusUpdate)
        for i in locations:
            limits = [max(0,i-window),min(scaffold_length,i+window)]
            scaffoldDensity[i] = np.log10(max(((limits[0] < positions) & (positions < limits[1])).sum(),0.3))
            biasedDensity[i] = np.log10(max(((limits[0] < biasedPositions) & (biasedPositions < limits[1])).sum(),0.3))
            if scaffoldDensity[i] > 50:
                hotspotlocations[scaffold].append(i)

        fig, ax = plt.subplots()
        fig.set_size_inches(3.42,2)

        MBlocations = locations/float(1000000)
        
        plt.plot(locations,scaffoldDensity,label='Observed',lw=0.25)
        plt.plot(locations,biasedDensity, '-C2', label='Simulated',lw=0.25)

        plt.xticks(fontsize=8)
        plt.xlabel('Position (Mb)',fontsize=8,labelpad=1)
        plt.ylabel('Inserts/Kb',fontsize=8,labelpad=1)
        plt.yticks([-0.5,0,1,2,3],[0,1,10,100],fontsize=8)
        plt.title('Insert density across scaffold '+scaffold.split(" ")[0],fontsize=10)
        plt.gcf().subplots_adjust(bottom=0.16, left=0.14, right=0.98, top=0.90)
        legend = plt.legend(loc='upper center',edgecolor=None,frameon=False,fontsize=8,ncol=2)
        plt.savefig(poolFileName+"_LargestScaffold.pdf")
        plt.close()

        statusUpdate =  "Saved plot of insertion density across the largest scaffold in "+poolFileName+"_LargestScaffold.pdf"
        printUpdate(options.logFile,statusUpdate)
 

if __name__ == "__main__":
    main(sys.argv[1:])
exit()




