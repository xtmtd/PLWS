# PLWS
Phylogenomics from Low-coverage Whole-genome Sequencing
We present a novel WGS-based pipeline for extracting essential phylogenomic markers through rapid genome assembling from low-coverage genome data, employing a series of computationally efficient bioinformatic tools. We tested the pipeline on a Hexapoda dataset and a more focused Phthiraptera dataset (genome sizes 0.1‒2 Gbp), and further investigated the effects of sequencing depth on target assembly success rate based on raw data of six insect genomes (0.1‒1 Gbp). Each genome assembly was completed in 2‒24 hours on desktop PCs. We extracted 872‒1,615 near-universal single-copy orthologs (BUSCOs) per species. This method also enables development of ultraconserved element (UCE) probe sets; we generated probes for Phthiraptera based on our WGS assemblies, containing 55,030 baits targeting 2,832 loci, from which we extracted 2,125‒2,272 UCEs. We also showed that only 10‒20× sequencing coverage was sufficient to produce hundreds to thousands of targeted loci from BUSCO sets, and even lower coverage (5×) was required for UCEs. Our study demonstrates the feasibility of conducting phylogenomics from low-coverage WGS for a wide range of organisms. This new approach has major advantages in data collection, particularly in reducing sequencing cost and computing consumption, while expanding loci choices.

![image](https://raw.githubusercontent.com/xtmtd/image/master/WGS_pipeline.png)
----------------------------------------------------------------------------------------------------------------------------------------
All four bash scripts are only tested in Centos operating system, and may be available for other Linux systems. Their functions are described below.
  1. script1_Genome_assembly.sh: Script 1 performs main rapid genome assembly steps, i.e. read quality trimming and normalization (BBtools),  error correction (Lighter), multi-k-mer assembly (Minia3), reduction of heterozygous contigs (Redundans), scaffolding (BESST), gap filling (GapCloser) and BUSCO assessment (BUSCO).
  2. script2_BUSCO_extraction.sh: Script 2 extracts single-copy orthologs (BUSCOs) from previous BUSCO assessments and generates nuclotide/protein alignment matrices of 50%-100% completeness (MAFFT, trimAl and FASconCAT-G) and partitioning schemes for phylogenetic analyses.
  3. script3_UCE_probe_design.sh: Script 3 can design UCE probe with tools samtools, art, BBtools, faToTwoBit and stampy.
  4. script4_UCE_extraction.sh: Script 4 extracts UCE loci from genome assemblies and generates alignment matrices of 50%-100% completeness and partitioning schemes for phylogenetic analyses. Tools MAFFT, seqkit, trimAl, faToTwoBit and FASconCAT are required.

----------------------------------------------------------------------------------------------------------------------------------------
Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source ate listed below.

  SRA-tools v2.9.0 (https://github.com/ncbi/sratoolkit), parallel-fastq-dump v0.6.3 (https://github.com/rvalieris/parallel-fastq-dump)
  BBTools v37.93 (https://sourceforge.net/projects/bbmap/), Lighter v1.1.1 (https://github.com/mourisl/Lighter), Minia v3.00-alpha1	(https://github.com/GATB/minia), Redundans v0.13c	(https://github.com/lpryszcz/redundans), Minimap2 v2.9	(https://github.com/lh3/minimap2), Samtools v1.7	(http://www.htslib.org/), BESST v2.2.8	(https://github.com/ksahlin/BESST)
  GapCloser v1.12	(http://soap.genomics.org.cn/), BUSCO v3.0.2	(http://busco.ezlab.org/), MAFFT v7.394	(https://mafft.cbrc.jp/alignment/software/), trimAl v1.4.1	(http://trimal.cgenomics.org/), FASconCAT-G v1.04	(https://github.com/PatrickKueck/FASconCAT-G), PHYLUCE v1.5.0	(http://phyluce.readthedocs.io/en/latest/index.html), faToTwoBit	(http://hgdownload.soe.ucsc.edu/admin/exe/), ART-20160605	(https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm), Stampy v1.0.32	(http://www.well.ox.ac.uk/project-stampy), SeqKit v0.8.0	(https://github.com/shenwei356/seqkit)

----------------------------------------------------------------------------------------------------------------------------------------
User manual

The requirements for each script have been described in the beginning of the bash script text. Values of variables/parameters, such as read length, tool path, number of threads etc., can be modified prior to analyses according to the status of sequencing data, computers etc. When everything is ready, just type 'sh ***.sh'.

----------------------------------------------------------------------------------------------------------------------------------------
Contact

Please send emails to Dr. Feng Zhang (xtmtd.zf at gmail.com).
