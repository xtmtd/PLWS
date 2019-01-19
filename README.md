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

  SRA-tools v2.9.0 (https://github.com/ncbi/sratoolkit)  
  parallel-fastq-dump v0.6.3 (https://github.com/rvalieris/parallel-fastq-dump)  
  BBTools v38.32 (https://sourceforge.net/projects/bbmap/)  
  Lighter v1.1.2 (https://github.com/mourisl/Lighter)  
  Minia v3.2	(https://github.com/GATB/minia)  
  Redundans v0.14c	(https://github.com/lpryszcz/redundans)  
  Minimap2 v2.12	(https://github.com/lh3/minimap2)  
  Samtools v1.9	(http://www.htslib.org/)  
  BESST v2.2.8	(https://github.com/ksahlin/BESST)  
  GapCloser v1.12	(http://soap.genomics.org.cn/)  
  BUSCO v3.0.2	(http://busco.ezlab.org/)  
  MAFFT v7.407	(https://mafft.cbrc.jp/alignment/software/)  
  trimAl v1.4.1	(http://trimal.cgenomics.org/)  
  FASconCAT-G v1.04	(https://github.com/PatrickKueck/FASconCAT-G)  
  PHYLUCE v1.6.6	(http://phyluce.readthedocs.io/en/latest/index.html)  
  faToTwoBit	(http://hgdownload.soe.ucsc.edu/admin/exe/)  
  ART-20160605	(https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)  
  Stampy v1.0.32	(http://www.well.ox.ac.uk/project-stampy)  
  SeqKit v0.9.3	(https://github.com/shenwei356/seqkit)
  
 To faciliate the use of our pipeline and simplify the installation of various packages, we also prepared the virtual mirror (.vmdk), including all the bioinformatic tools and custom scripts ready in the CENTOS 7.3 system. It can be directly opened by VMware or VirtualBox, which are often installed on Windows systems. Please download it here: https://1drv.ms/f/s!Ak7sQGBGqlguq1gdk1KKIR7uwAfj. More than 100 Gb (at least 50) disk space is recommented for its use. Below is related information on CENTOS system:  

  operating system: centos 7.3  
  passwords of accounts, i.e. root and zf: 1  
  Most packages are installed under /usr/local/bin, /home/zf/install/ etc. Environmental variables and paths can be checked by "vi ~/.bashrc" and "printenv". Most executables can be directly performed without typing the installation paths. Some should be executed by initiating conda environments (/home/zf/install/miniconda/envs/), such as phyluce, just type 'source activate phyluce'.

----------------------------------------------------------------------------------------------------------------------------------------
User manual

The requirements for each script have been described in the beginning of the bash script text. Values of variables/parameters, such as read length, tool path, number of threads etc., can be modified prior to analyses according to the status of sequencing data, computers etc. When everything is ready, just type 'sh ***.sh'.

Below are some notes for each script although more details have been described in the script.

Script1: 
1) Modify the names of input raw paired read files as 1.raw.fq.gz and 2.raw.fq.gz in the working folder.  
2) The default starting kmer value is 21 and thus kmer values are 21, 21+20, 21+2*20....  
3) All the statistics of input reads and assemblies generated in the whole asembly progress are summerized respectively /bbtool/bbtool.log and /minia/assembly.log.  
4) Tools pigz, BBTools, Lighter, Minia, redundans, Minimap2, samtools, BESST, GapCloser and BUSCO may be used in this script.

Script2:  
1) Copy all the BUSCO results (run_* folders) into the same subfolder ./raw_busco.  
2) Tools MAFFT, trimAl and FASconCAT-G are used in this script and the former two are recommended to be installed in environmental paths.

Script3:   
1) Copy all the genome assemblies to working_folder/DIR_ASSEMBLY/ and replace assembly name using "species_name.fa". 
2) Tools samtools, art, BBtools, faToTwoBit and stampy are used in this script and the former two are recommended to be installed in environmental paths.

Script4:  
1) Prepare a probe set file (***.fasta) in the initial working folder.
2) Modify all genome assemblies endding with .fa ("species_names.fa") and copy them to working_folder/genomes/ 
3) Tools MAFFT, seqkit, trimAl, faToTwoBit and FASconCAT are used in this script and the former three were recommended to be installed in environmental paths.

----------------------------------------------------------------------------------------------------------------------------------------
Citation
Feng Zhang, Yinhuan Ding, Chaodong Zhu, Xin Zhou, Michael C. Orr, Stefan Scheu, Yun-Xia Luan. 2018. Phylogenomics from Low-coverage Whole-genome Sequencing. Methods in Ecology and Evolution. (accepted)

----------------------------------------------------------------------------------------------------------------------------------------
Contact

Please send emails to Dr. Feng Zhang (xtmtd.zf at gmail.com).
