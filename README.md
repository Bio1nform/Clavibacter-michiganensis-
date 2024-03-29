# Genomic Analysis of Clavibacter michiganensis Reveals Insight Into Virulence Strategies and Genetic Diversity of a Gram-Positive Bacterial Pathogen

Shree P. Thapa, Sivakumar Pattathil, Michael G. Hahn, Marie-Agnès Jacques, Robert L. Gilbertson, Gitta Coaker


## Comparative genomics

#### Genome sequencing and annotation
Sequenced genomes were de novo assembled with the SPAdes v. 3.10 (Bankevich et al., 2012). Draft genomes were annotated with Prokka v. 1.11 and the NCBI Prokaryotic Genome Annotation Pipeline (Seemann, 2014; Tatusova et al., 2016).

#### Command used for SPAdes: 
    spades.py --careful -k 21,33,55,77,99,127 -s sample.fastq -o sample.out
    
#### Command used for PROKKA: 
    prokka --compliant --prefix sample  --centre CoUCD --force --addgenes --genus Clavibacter --species michiganensis --strain sample --locustag $sample $sample.fasta --outdir  samplePRoKKA

Bioinformatics analysis

Sec-dependent effectors, protein sequences were screened for Sec signal peptides. Proteins possessing signal pep-tides for the Sec-dependent pathway were identified using Signalv. 3.0, SignalP v. 4.0, and Phobius v. 1.01 (Bendsten et al., 2004; Kall et al., 2004; 2007; Petersen et al., 2011). 

#### Command used for SignalP3:  
    signalp -t gram- -f short -u 0.44 /Whole_genome_analysis/Clavibacter/$sample.fasta > / Whole_genome_analysis/Clavibacter/INTERPROSCAN/Signalp4.1/$sample_signP4_OPR.out

#### Command used for SignalP4:
    interproscan.sh -f TSV -appl SignalP-GRAM_NEGATIVE -i /Whole_genome_analysis/Clavibacter/Proteins/$sample.fasta -b /Whole_genome_analysis/Clavibacter/INTERPROSCAN/Signalp4.1/$sample.iprscan.signalp_4

#### Command used for Phobius:
    interproscan.sh -f TSV -appl Phobius -i /Whole_genome_analysis/Clavibacter/Proteins/$sample.fasta -b /Whole_genome_analysis/Clavibacter/INTERPROSCAN/Phobius/$sample.phobius

Predicted lipoproteins and transmembrane proteins were filtered from the Sec secretomes.

#### Command used for lipoP: 
    perl LipoP -short -html /Whole_genome_analysis/Proteins/$sample.faa > $sample_.lipoP

Transmembrane topology was predicted using TMHMM v. 2.0 (Sonnhammer et al., 1998). 
#### Command used for TMHMM:
     interproscan.sh -f TSV -appl TMHMM -i /Whole_genome_analysis/Clavibacter/Proteins/compliantFasta/$sample.fasta -b /Whole_genome_analysis/Clavibacter/INTERPROSCAN/TMHMM/$sample.TMHMM


## Phylogenetic analyses
Orthologous genes of Las isolates were predicted using the OrthoMCL v. 2.0 pipeline (Li et al., 2003). 

#### orthomcl clustering
    blastall -d goodNucleotides.fasta -p blastn -i goodNucleotides.fasta -m 8 -e 1e-5 -o goodNucleotides_blast_result.out
    orthomclBlastParser /Clavibacter/OrthoMCL/Genes/goodprotein_blast_result.out /Clavibacter/OrthoMCL/Genes/compliantFasta > similarSequences.txt
    orthomclInstallSchema orthomcl.config mysql.log
    orthomclLoadBlast orthomcl.config similarSequences.txt
    orthomclPairs orthomcl.config orthomcl_pairs.log cleanup=no
    orthomclDumpPairsFiles orthomcl.config > dump.log
    mcl mclInput --abc -t 10 -I 1.5 -o mcl_Output
    orthomclMclToGroups CLas 00001 < mcl_Output > groups.txt

Multiple alignments of gene sequences were done with PRANK v. 170,427 (Löytynoja, 2014). 

#### Command used for PRANK: 
      prank -d=CLas.fasta -o=CLas_prank.fasta -DNA -F

#### All the alignments were concatenated by FASconCAT v. 1.1, yielding a gene supermatrix (Kuck and Meusemann, 2010). 
    perl FASconCATv.1.1.pl

A maximum-likelihood approach was used to reconstruct the phylogenetic tree using RAxML v. 8.2 software (Stamatakis, 2006). 
#### Command used for RAxML:
      raxmlHPC-PTHREADS-SSE3 -k -f a -m GTRGAMMA -d -p 12345 -x 12345 -s FcC_smatrix.fasta -n FcC_GAMMA_raxML_N.out -N 1000 -T 12

## Citation

You can cite our paper as:

Thapa SP, Pattathil S, Hahn MG, Jacques MA, Gilbertson RL, Coaker G, 2017. Genomic analysis of Clavibacter michiganensis reveals insight into virulence strategies and genetic diversity of a Gram‐positive bacterial pathogen. Molecular Plant‐Microbe Interactions 30, 786–802. https://doi.org/10.1094/MPMI-06-17-0146-R
