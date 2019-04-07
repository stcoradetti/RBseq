# RBseq
Track fitness of deletion mutants with a randomly barcoded random insertion library

Please note, this is scientific software for research purposes only.  It is offered freely with no guarantees whatsoever.

Updated 5 April 2019

--OVERVIEW--
RBseq is a python-based pipeline for performing genome-wide experiments to quantify the relative fitness of random insertion mutants with unique sequence barcodes.  This software was developed specifically with fungal genomes in mind, but theoretically should work with any eukaryotic genome.

The RBseq package consists of four programs which each carry out a distinct step required for fitness analysis. The process is (1) map insertion locations and associate barcodes, (2) annotate the genomic context of those insertions, (3) count barcodes in different conditions, and then (4) calculate fitness scores for each gene, aggregating data from several barcoded insertions per gene.  Once a given 'pool' of mutant strains has been characterized by (1) and (2), (3) and (4) can be repeated numerous times on diverse experiments starting with the same mutant pool.  For a detailed example of this process see "Rapid Quantification of Mutant Fitness in Diverse Bacteria by Sequencing Randomly Bar-Coded Transposons" by Wetmore et al 2015 (DOI: 10.1128/mBio.00306-15) and "Functional genomics of lipid metabolism in the oleaginous yeast Rhodosporidium toruloides" by Coradetti et al 2018 (DOI: 10.7554/eLife.32110).

RBseq_Map_Insertions.py
This program takes Illumina sequence reads as input and associates unique sequence barcodes with unique genomic insertion locations.  Typically these are sequence reads from a TnSeq experiment or similar variation in which a known sequence in the inserted DNA is used to amplify the insertion/genome junction.  While it is agnostic as to the method of insertion library construction, it assumes some rate of concatamer formation.  That is, it does not assume inserted DNA is always incorporated into the genome in single copy at any given location.  The type of insertion is tracked (concatameric vs single copy) and barcode that cannot be unambiguously mapped to single locus are filtered out. It also does not assume perfect editing of the insertion, some loss of sequence at the insertion/genome junction is tolerated.  These features were included specifically for analysis of libraries constructed with Agrobacterium tumefaciens mediated transformation, but the software should be compatible with any random insertion library in which a defined length of random sequence barcodes are placed within a short sequence read of identifiable genomic sequence.

RBseq_Annotate_Insertions.py
This program uses a gff formatted file of gene locations to assign the mapped insertion locations as intergenic or as interrupting a protein coding gene.  It also computes some summary information about the insertion pool for method development, such as the rate insertion in various genomic features, biases for GC content or across different scaffolds, etc.

RBseq_Count_BarCodes.py
This program takes as input Illumina sequence reads of barcode sequences (BarSeq) from a collection of experimental conditions.  It computes some summary statistics on the observed barcode sequences and collects counts for barcodes mapped and annotated by RBseq_Map_Insertions.py and RBseq_Annotate_Insertions.py

RBseq_Calculate_Fitness.py
This program takes as input the barcode counts from RBseq_Count.py and calculates a fitness score for each gene [log2(condition abundance/reference abundance)] with a sufficient number of barcoded insertions with sufficient sequence depth, as well as statistical significance of those fitness scores in the form of T-like statistics (as described in Wetmore et al 2015, DOI: 10.1128/mBio.00306-15) P-values from those T-like statistics and a more conservative P-value from a Wilcoxon Signed Rank Test.  


--PREREQUISITES--
-Python 3.7
-Biopython 1.72 or later
-Anaconda scientific computing environment or the following python libraries: numpy, pandas, matplotlib, json, scipy, and stats models
-NCBI BLAST+ 2.2.30 or later

--USAGE--
Each script requires a metadata file to list input sequencing files, required reference files and some common options, including output file locations.  These text files are passed to the script with the -m or --metafile flag.  The scripts provide a series of status updates to STDOUT and simultaneously save them in log file for easy reference later. By default, the log files are given a unique name based on the date, but you can change this with the optional -l or --logFile flag.

--
RBseq_Map_Insertions.py Inputs
python RBseq_Map_Insertions.py --metafile METAFILE [--logFile LOGFILE --qual 10 --matchBefore 6 --matchAfter 6 --matchJunction 6 --maxFillerSeq 10 --barcodeVariation 2 --wobbleAllowed 1 --scoreDiff 10 --minFrac 0.8  --filterNeighborhood 10 --filterEditDistance 5 --minPercentID 95 --maxEvalue 0.1 --noInsertHits --noBarcodes --useMappedFiles]

  -m/--metafile METAFILE
    Metadata file for TnSeq runs. A tab-delimited file with columns as follows:

      Pool
      A short name for the mutant pool being analyzed. Should represent a physical sample of thousands to millions of unique mutation events in one mixed culture.  Multiple mutant pools can be analyzed with one metafile.

      ShortName	
      A short, unique name for each fastq file.
        
      Fastq
      Full or relative path to fastq file, e.g sequenceFolder/experiment/run.fastq

      ReadModel	
      Full or relative path to a text file describing the expected structure of TnSeq reads.  The file consists of 4 lines.  Line 1 denotes the number of initial variable bases to ignore as 'N's.  Line 2 is the sequence preceding the barcode.  Line 3 is the expected length of sequence barcode denoted as 'N's.  Line 4 is the sequence between the barcode and the junction with the genome.  Example:

nnnnnn
GATGTCCACGAGGTCTCT
NNNNNNNNNNNNNNNNNNNN
CGTACGCTGCAGGTCGACAATGATCCAAACTATCAGTGTTTGA	

      InsertionSequence
      Full or relative path to a text file describing giving the full sequence of each insertion.  This should be a fasta format file with two entries labeled 'insert' and 'plasmid' giving the full sequence of the expected insertion into the genome and and the remaining source plasmid sequence, if appropriate. Example:

>insert
agtcac...
>plasmid
actgact... 
      
      GenomeSequence
      Full or relative path to the basename of a BLAST database for the relevant genome.  This should be the same basename of the fasta source file.  Example: 'path_to/dir/genome' should point to a folder 'dir' containing 'genome.fasta', 'genome.nhr', 'genome.nin', and 'genome.nsq'

      OutputDir	
      Full or relative path to the desired directory to save output files.

      GeneLocations
      Full or relative path to a GFF file describing gene locations

      GeneIdentifier
      Name of unique gene identifiers from the GFF file to use as gene names.  This should be one of the semicolon-delimited fields in mRNA entries in the GFF.  If nothing is entered in this column 'Parent' will be used as default.

      GeneAnnotations
      Full or relative path to a tab-delimited file with columns labeled 'ID', 'AlternateID', and 'Annotation'.  ID should be the unique gene identifier you are using as gene names and should match the field you specified in the metafile column 'GeneIdentifier'.  AlternateID is an alternate, perhaps more convenient gene name.  ID will likely be a cryptic systematic name such as 'RTO4_01234', but Alternate ID might be 'Arg9' or something more descriptive.  'Annotation' is a brief text describing the predicted gene function.  This file allows you to set your preferred gene name and annotations for final output with fitness data later.
	
  -l/--logFile LOGFILE
    File to write run log to. Default is Map_TIMESTAMP.log

  -q/--qual MINQUAL
    Minimum quality score for the barcode region for a read to be mapped. Default is 10

  -b/--matchBefore MATCHBEFORE
    Number of bases in the constant region before the barcode to match. Default is 6

  -a/--matchAfter MATCHAFTER
    Number of bases in the constant region after the barcode to match. Default is 6. 

  -j/--matchJunction MATCHJUNCTION
    Number of bases to match at the junction with the genome. Default is 5. Because under many library construction strategies there is often loss of sequence at the end of the inserted fragment, this script attempts to find a short segment of sequence matching the end of the insertion, scanning for the closest matching segment to the end of the sequence.  Sequence after that furthest match is used as the junction with the native genome sequence.

  -f/--maxFillerSeq MAXFILLERSEQ
    NOT IMPLEMENTED YET. The maximum number of bases allowed between the detected end of the insertion and the start of a match to the genome. This could be used in future versions to report insertions with long pieces of filler DNA added at insertion sites.

  -v/--barcodeVariation BARCODEVARIATION
    Allowed deviation from expected barcode length. Default is 2. Depending on the method of library construction, significant numbers of barcoded insertions may have slightly shorter or longer random barcodes.  Allowing flexibility in barcode length can recover these insertions, but slows down processing.

  -w/--wobbleAllowed WOBBLEALLOWED
    Allowed deviation from expected position in the read for the matching sequence around the barcode and the junction with the genome. Default is 1.  Indels introduced near the barcode can leave some insertion unmapped without some allowance for small deviations from expected positions of search sequences.

  -s/--scoreDiff SCOREDIFF
    For a read mapping to the genome to be counted as unique, the best-hit bitscore must be at least this much greater than the second-best hit bitscore. Default is 10

  -F/--minFrac MINFRACTION
    Minimum fraction of reads mapping to a single location before barcode is classified as multilocus. Default is 0.9

  -N/--filterNeighborhood FILTERNEIGHBORHOOD
    Number of basepairs around each insertion to filter out less abundant insertions with similar barcodes. Default is 10

  -D/--filterEditDistance FILTEREDITDISTANCE
    Nearby barcodes within this Levenshtein edit distance will be filtered. Default is 5. This argument is ignored if --noBarcodes is passed (all insertions within --filterNeighborhood will be filtered)

  -p/--minPercentID MINPERCENTID
    Minimum percent identity required when BLASTing sequences to the genome or insert. Default is 95

  -e/--maxEvalue MAXEVALUE
    Maximum E-value allowed when BLASTing sequences to the genome or insert. Default is 0.1

  -u/--useMappedFiles
    Bypass the mapping step and work from previously mapped files in the output directory. Useful for testing the effects of some parameters without repeating time consuming mapping of reads.

  -I/--noInsertHits
    Skip search for additional insertion sequence after the expected junction with the genome. Default behavior is to look for these sequences, typically indicating concatameric insertions. Turning it off will save time if concatamers aren't common in your mutant pool.

  -n/--noBarcodes
    This flag is used for optimizing library construction protocols with non-barcoded sequences.  As generating a diverse pool of insertion constructs with unique barcodes can be time consuming, a test insertion library can be created with constructs containing just one constant barcode. If TnSeq is performed on the resulting library, the distinct locations of insertions can be mapped, but they will all have the same sequence barcode associated with the insertion.  To allow characterization of the library with this script, this flag signals to use a snippet of the associated genomic sequence as a substitute for a sequence barcode on the insertion.  This allows for distinct insertions to be mapped and estimation of library diversity and identification of biases in insertion position, etc.

--
RBseq_Map_Insertions.py Outputs


  POOL_poolfile
    Tab-delimited file describing the pool of insertion mutants for which unique barcodes mapped to a single unambiguous genomic location.  This file is used for downstream analysis. Other outputs are provided for troubleshooting analysis and library construction.
  POOL_poolfile_ambiguous
    Barcodes ambiguously mapping to two or more locations.  These are barcodes inserted into sequence with one or more very similar sequences in the genome, or with poor quality reads for mapping the junction to the genome, creating uncertainty as to where exactly the barcode is inserted.
  POOL_poolfile_filteredLocal
    Barcodes that were mapped close to a very similar, more abundant barcode
  POOL_poolfile_insertReadsOnly
    Barcodes that only mapped to insert sequences.  These could be barcodes in concatameric insertions for which the junction with the genome is difficult to amplify or in multi-barcode concatamers
  POOL_poolfile_multiLocus
    Barcodes with reads mapping to two or more distinct locations.  Most likely the same sequence barcode inserted in two or more different genomic locations
  POOL_poolfile_offByOne    Barcodes that differed by just one base from much more abundant barcodes.  Almost certainly sequencing errors.  
  TNSEQ_FILE_blastGenome.txt
    Results from the BLAST search of sequences from TNSEQ_FILE against the genome.
  TNSEQ_FILE_blastInsert.txt
    Results from the BLAST search of sequences from TNSEQ_FILE against the genome.
  TNSEQ_FILE_blastquery.fasta
    Sequences blasted agains the genome/insert from TNSEQ_FILE
  TNSEQ_FILE_mapped.txt
    Summary of each read blasted against the genome/insert, barcode sequence found and information about the blast results
  TNSEQ_FILE_noBarcodeFound.txt
    Reads for which the compliant barcode sequences were not found.

--
RBseq_Annotate_Insertions.py Inputs
python RBseq_Annotate_insertions.py -m metadatafile [-l logfile -P 500 -T 100]

  -m/--metafile METAFILE
    Metadata file for TnSeq runs. Same as for RBseq_Map_Insertions.py.  For each unique entry in the Pool column of the metafile, this script will look for files named POOL_poolfile in the matching 'OutputDir' column.  Gene annotations are taken from the GFF file specified in the column 'GeneLocations'.  This GFF should contain mRNA features and CDS features to identify transcript and exon boundaries.  The Parent field in the mRNA features is the default source of gene ids to add to the pool file, but any field specified by Name=value in the GFF's mRNA features can be used by setting the 'GeneIdentifier' field in the metafile to Name.  The 'GeneAnnotations' field is the path and name of tab-delimited text file containing custom annotations and alternate names for your genes. This file should contain the headers 'ID','AlternateID', and 'Annotation' where 'ID' contains the same identifiers as present in the GFF file.

  -l/--logFile LOGFILE
    File to write run log to. Default is Annotate_TIMESTAMP.log

  -P/--promoterLength PROMOTERLENGTH
    Number of bases before transcriptional start to count as promoter sequences. Default is 500

  -T/--terminatorLength TERMINATORLENGTH
    Number of bases after transcriptional stop counted as terminator sequence. Default is 100
    
--
RBseq_Annotate_Insertions.py Outputs
  
  POOL_poolfile_annotated
    The input poolfile with information on the nearest gene to each insertion and the type of genomic feature the insertion is in (promoter, exon, intergenic, etc)

  POOL_poolfile_GChistogram.pdf
    A histogram of GC content around insertion sites versus GC content in the genome.  Included to detect GC bias in insertion rates.
  POOL_poolfile_insertionDistribution.pdf
    A graph comparing frequency of insertions in promoters, exons, etc versus percent of the genome those features comprise. Included to detect feature bias in insertion rates.  With many insertional strategies, higher insertion rates are common in more open regions such as promoters. In some species or with some methods this bias may be extreme enough as to make generation of a large number of insertions in coding regions much more difficult.
  POOL_poolfile_insertionDistribution.txt
    A text file giving raw data for the insertion distribution graph
  POOL_poolfile_InsertsPerScaffold.pdf
    Graphs of insertion number versus scaffold length. Useful for detecting bias in insertion rates or problems with genome assembly or insertion mapping.
  POOL_poolfile_LargestScaffold.pdf
    Graph of insertion density in a 1000 base-pair rolling window across the largest scaffold, and predicted insertion density if insertions occurred randomly at the same relative rate in promoters, eons, etc as measured in this library.  Useful for understanding how well distributed insertions occur.  Empirically most insertion libraries are not quite random, that is less evenly distributed than a random model would predict, but useful libraries don't have any large scale biases in insertion rates or particularly dominant insertion hot-spots.

--
RBseq_Count_Barcodes.py Inputs
python RBseq_Count_Barcodes.py -m metadatafile [-l logfile -b 6 -a 6 -q 10]

  -m/--metafile METAFILE
    Metadata file for BarSeq runs. Contains the following columns:

      Fastq
      Full or relative path to barseq fastq files.

      SampleName
      Short unique name for each physical sample.  One physical sample can have multiple fastq files with their own parameters and dual indexes. All counts across the different fastq files will be aggregated for fitness analysis.  So biological replicates should have different sample names such as treatment_a, treatment_b, treatment_c, control_a, etc.  Whereas technical replicates should have the same sample names.

      UsePrecounted
      TRUE or FALSE.  If TRUE then computed barcode counts from a previous run will be loaded from the output directory.  Useful when re-running analysis on an updated pool mapping or annotation.

      Poolfile
      Full or relative path to annotated poolfile from RBseq_Annotate_Insertions.py

      OutputDir
      Full or relative path to output directory

      minRandom
      The minimum number of expected random bases in the barseq reads.  Often barseq primers contain a number of leading random bases to assist in clustering on the illumina sequencer

      maxRandom
      The maximum number of expected random bases.  In some implementations barseq primers are mixtures with different lengths of leading random bases.  In these cases minRandom and maxRandom should cover the expected range of leading random bases, otherwise maxRandom should be equal to minRandom.

      DualIndex
      In some applications BarSeq primers will contain a second index sequence to guard against 'index hopping' on the illumine sequencer.  In cases for which the dual indexes have not be de-multiplexed by the sequencer, this field can be used to specify the expected dual index sequence.  This is assumed to appear in the BarSeq read after the leading random bases and before the sequence specified in 'BeforeBarcode'
      
      BeforeBarcode
      The constant sequence before the barcode sequence.  Not all of this sequence need be matched to identify the barcode, on the range specified with the -b flag.

      BarcodeLengths
      The expected length of the barcode.  Can be a single number e.g. 20, or a comma-separated list with no spaces of allowed lengths, e.g. 20,19,21 This list should be in order of expected abundance for fastest processing.

      AfterBarcode
      The constant sequence after the barcode sequence.  Not all of this sequence need be matched to identify the barcode, on the range specified with the -a flag.


  -l/--logFile LOGFILE
    File to write run log to. Default is Count_TIMESTAMP.log

  -q/--qual MINQUAL
    Minimum quality score for the barcode region for a read to counted. Default is 10

  -b/--matchBefore
    Number of bases before the barcode to match. Default is 6

  -a/--matchAfter
    Number of bases to match after the barcode. Default is 6.  Note that in some cases this search region may be approaching the end of 50bp BarSeq reads.  If matchAfter happens to span the end of the read, and thus there aren't as many bases to match as set here, but there are 3 or more bases to mach, as many bases are available up to the number specified will be matched.  If there aren't 3 or more bases to match, then the read will be discarded.

RBseq_Count_Barcodes.py Outputs
  
  RBseq_Count_Barcodes.py will summarize the number of barcodes seen for each sample, and use the most abundant barcode to estimate overall error rate.  It will also offer an extremely approximate guess at to the real diversity of barcodes present in the physical sample using the formula postulated by (Chao 1987): N seen once squared divided by two times N seen twice.  This number is at best a guess within a few orders of magnitude in this context and is really just confirmation that there is likely a large diverse population of unique barcodes in your sample or there isn't.

  poolCount.txt
    Tab-delimited file with uniquely mapped barcodes in the insertion pool, some information about their location and genomic context, and the number of times each was seen in each sample.

  countsFiles
    A directory of files with raw barcode counts seen in each sample.  These include barcodes not in the mapped mutant pool. These files can be used to repeat counts with an updated poolfile without the time consuming step of finding and counting barcodes in all reads.

--
RBseq_Calculate_Fitness.py Inputs
python RBseq_Count_Barcodes.py -m metadatafile [-l logfile -L 0 -s 5 -e 95 -i 3 -g 20 -I 3 -W 50 -P -B]

  -m/--metafile METAFILE
    Metadata file for fitness calculation. Contains the following columns:

    parser.add_argument("-m", "--metafile", dest="metafile", help="Metadata file for BarSeq runs. A tab-delimited file with columns titled Sample, Group, Reference, OutputDir, PoolCountFile.  Entires in the Sample and Reference columns must correspond to column heading in the poolCountFile. Samples from biological replicates ",default="metafile.txt")
    parser.add_argument("-l", "--logFile", dest="logFile", help="File to write run log to. Default is FitnessLog_TIMESTAMP.log",default="FitnessLog_"+timestamp+".log")
    parser.add_argument("-L", "--normLocal", dest="normLocal", help="-L [int]: Normalize fitness scores such that local windows of [int] contiguous insertions have a median fitness score of zero.  Default behavior is to normalize fitness scores accross full contigs such that the median fitness score is zero", default=0, type=int)
    parser.add_argument("-s", "--fitnessStartBoundary", dest="fitnessStartBoundary", help="Integer 0-100. Fraction of the region between coding start and coding stop to set as start boundary for calculating gene fitness.  Default is 5: i.e. insertions in the first 5 percent of the region between start and stop will not be included in calculations for gene fitness scores or statistical tests", default=5, type=int)
    parser.add_argument("-e", "--fitnessEndBoundary", dest="fitnessEndBoundary", help="Integer 0-100. Fraction of the region between coding start and coding stop to set as end boundary for calculating gene fitness.  Default is 95: i.e. insertions in the last 5 percent of the region between start and stop will not be included in calculations for gene fitness scores or statistical tests", default=95, type=int)
    parser.add_argument("-i", "--minInsertionCounts", dest="minInsertionCounts", help="Integer. Minimum counts each insertion must have between both the test condition and the reference condition to be included in calculations for gene fitness scores or statistical tests. Default is 3", default=3, type=int)
    parser.add_argument("-g", "--minGeneCounts", dest="minGeneCounts", help="Integer. Minimum total counts across all insertions in a given gene before gene fitness scores are calculated. Default is 20", default=20, type=int)
    parser.add_argument("-I", "--minInsertions", dest="minInsertions", help="Integer. Minimum number of insertions with sufficient counts in a given gene before gene fitness scores are calculated. Default is 3", default=3, type=int)
    parser.add_argument("-W", "--maxWeightCounts", dest="maxWeightCounts", help="Integer. Maximum number of counts to be used for gene fitness score calculation.  E.g. if set at 50 then two insertions with 50 and 100 counts would be wieghted equally in computing gene fitness, but an insertion with 10 counts would have a smaller weight ", default=50, type=int)
    parser.add_argument("-P", "--smartPseudoCounts", dest="smartPseudoCounts", action='store_false', help="If this flag is passed, fitness scores will be computed without 'smart' pseudocounts as in Wetmore et al 2015.", default=True)
    parser.add_argument("-B", "--fitnessBrowserOutput", dest="fitnessBrowserOutput", action='store_true', help="If this flag is passed, output files will be saved in formats compatible with the Arkin Lab Fitness Browser", default=False)

RBseq_Calculate_Fitness.py Outputs


--LICENSE--
Copyright (C) 2019 Samuel Coradetti and the United States Department of Energy. All rights reserved.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

--DISCLAIMER--
NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY NOR THE REGENTS OF THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.