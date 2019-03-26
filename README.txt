# RBseq
Track fitness of deletion mutants with a randomly barcoded random insertion library

Please note, this is scientific software for research purposes only.  It is offered freely with no guarantees whatsoever.

Updated 25 March 2019

--OVERVIEW--
RBseq is a python-based pipeline for performing genome-wide experiments to quantify the relative fitness of random insertion mutants with unique sequence barcodes.  This software was developed specifically with fungal genomes in mind, but theoretically should work with any eukaryotic genome.

The RBseq package consists of four programs which each carry out a distinct step required for fitness analysis.

RBseq_Map_Insertions.py
This program takes Illumina sequence reads as input and associates unique sequence barcodes with unique genomic insertion locations.  Typically these are sequence reads from a TnSeq experiment or similar variation in which a known sequence in the inserted DNA is used to amplify the insertion/genome junction.  While it is agnostic as to the method of insertion library construction, it assumes some rate of concatamer formation.  That is, it does not assume inserted DNA is always incorporated into the genome in single copy at any given location.  The type of insertion is tracked (concatameric vs single copy) and barcode that cannot be unambiguously mapped to single locus are filtered out. It also does not assume perfect editing of the insertion, some loss of sequence at the insertion/genome junction is tolerated.  These features were included specifically for analysis of libraries constructed with Agrobacterium tumefaciens mediated transformation, but the software should be compatible with any random insertion library in which a defined length of random sequence barcodes are placed within a short sequence read of identifiable genomic sequence.

RBseq_Annotate_Insertions.py
This program uses a gff formatted file of gene locations to assign the mapped insertion locations as intergenic or as interrupting a protein coding gene.  It also computes some summary information about the insertion pool for method development, such as the rate insertion in various genomic features, biases for GC content or across different scaffolds, etc.

RBseq_Count_BarCodes.py
This program takes as input Illumina sequence reads of barcode sequences (BarSeq) from a collection of experimental conditions.  It computes some summary statistics on the observed barcode sequences and collects counts for barcodes mapped and annotated by RBseq_Map_Insertions.py and RBseq_Annotate_Insertions.py

RBseq_Calculate_Fitness.py
This program takes as input the barcode counts from RBseq_Count.py and calculates a fitness score for each gene [log2(condition abundance/reference abundance)] with a sufficient number of barcoded insertions with sufficient sequence depth, as well as statistical significance of those fitness scores in the form of T-like statistics (as described in Wetmore et al 2015, DOI: 10.1128/mBio.00306-15) P-values from those T-like statistics and a more conservative P-value from a Wilcoxon Signed Rank Test.  


--PREREQUISITES--
Python 3.7
Anaconda scientific computing environment
Biopython 1.72 or later
NCBI BLAST+ 2.2.30 or later


--LICENSE--
Copyright (C) 2019 Samuel Coradetti and the United States Department of Energy. All rights reserved.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

--DISCLAIMER--
NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY NOR THE REGENTS OF THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.