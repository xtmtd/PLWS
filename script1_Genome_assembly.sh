#!/bin/bash
#2020.07.08 by ZF

#Type 'sh script1_Genome_assembly.sh forward_reads_file reverse_reads_file', e.g. sh script1_Genome_assembly.sh 1.raw.fq.gz 2.raw.fq.gz
#Most executables are recommended to be added into the environmental paths
#Tools pigz, BBTools, Minia, redundans, Minimap2, samtools, BESST, GapCloser and BUSCO may be used and will be automatically checked prior to formal analyses in this script
#the default starting kmer value is 21 and thus kmer values are 21, 21+20, 21+2*20....
#the statistics of assemblies generated in the asembly procedure are summerized in assembly.statistics



##Checking the package dependency
echo "Checking the package dependency......" | tee -a log.txt
echo -e "\n" >> log.txt

#check pigz
if [ $(which pigz) ]
    then
      echo "pigz ...... OK" | tee -a log.txt
      EXE_PIGZ=$(which pigz)
      DIR_PIGZ=${EXE_PIGZ%/*}
    else
      until [ -x $DIR_PIGZ/pigz ]
        do
          read -p "Pigz is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_PIGZ
        done
      echo "pigz ...... OK" | tee -a log.txt
fi

#check BBtools
if [ $(which bbduk.sh) ]
    then
      echo "BBtools ...... OK" | tee -a log.txt
      EXE_BBTOOLS=$(which bbduk.sh)
      DIR_BBTOOLS=${EXE_BBTOOLS%/*}
    else
      until [ -x $DIR_BBTOOLS/bbduk.sh ]
        do
          read -p "BBtools is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/bbtools):      " DIR_BBTOOLS
        done
      echo "BBtools ...... OK" | tee -a log.txt
fi

#check Minia
if [ $(which minia) ]
    then
      echo "Minia ...... OK" | tee -a log.txt
      EXE_MINIA=$(which minia)
      DIR_MINIA=${EXE_MINIA%/*}
    else
      until [ -x $DIR_MINIA/minia ]
        do
          read -p "Minia is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_MINIA
        done
      echo "Minia ...... OK" | tee -a log.txt
fi

#check the step of normalization
read -p "Execute the step of normalization of sequencing data? It can be helpful for reducing the data size and often improve the assembly. It can be skipped but the assembly time would be longer:  y or n      " NORM
echo "NORMALIZATION=$NORM" >> parameters.cfg

#check Redundans
read -p "Execute the step of reduction of heterozygous contigs? It can be skipped in most cases. If heterozygosity of the final assembly is high (i.e. value of D% in BUSCO assessments is high), you can use some tools, such as Redundans, to reducen them: y or n      " HETER
echo "REDUNDANS=$HETER" >> parameters.cfg

if [ "$HETER" == "y" ]
  then
    if [ $(which redundans.py) ]
      then
        echo "Redundans ...... OK" | tee -a log.txt
        EXE_REDUNDANS=$(which redundans.py)
        DIR_REDUNDANS=${EXE_REDUNDANS%/*}
      else
        until [ -x $DIR_REDUNDANS/redundans.py ]
          do
            read -p "Redundans is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/redundans):      " DIR_REDUNDANS
          done
        echo "Redundans ...... OK" | tee -a log.txt
    fi
  else
    echo "Reduction of heterozygous contigs will be skipped."
fi

#check Minimap2
if [ $(which minimap2) ]
    then
      echo "Minimap2 ...... OK" | tee -a log.txt
      EXE_MINIMAP2=$(which minimap2)
      DIR_MINIMAP2=${EXE_MINIMAP2%/*}
    else
      until [ -x $DIR_MINIMAP2/minimap2 ]
        do
          read -p "Minimap2 is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_MINIMAP2
        done
      echo "Minimap2 ...... OK" | tee -a log.txt
fi

#check Samtools
if [ $(which samtools) ]
    then
      echo "Samtools ...... OK" | tee -a log.txt
      EXE_SAMTOOLS=$(which samtools)
      DIR_SAMTOOLS=${EXE_SAMTOOLS%/*}
    else
      until [ -x $DIR_SAMTOOLS/samtools ]
        do
          read -p "Samtools is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_SAMTOOLS
        done
      echo "Samtools ...... OK" | tee -a log.txt
fi

#check BESST
if [ $(which runBESST) ]
    then
      echo "BESST ...... OK" | tee -a log.txt
      EXE_BESST=$(which runBESST)
      DIR_BESST=${EXE_BESST%/*}
    else
      until [ -x $DIR_BESST/runBESST ]
        do
          read -p "BESST is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_BESST
        done
      echo "BESST ...... OK" | tee -a log.txt
fi

#check GapCloser
if [ $(which GapCloser) ]
    then
      echo "GapCloser ...... OK" | tee -a log.txt
      EXE_GAPCLOSER=$(which GapCloser)
      DIR_GAPCLOSER=${EXE_GAPCLOSER%/*}
    else
      until [ -x $DIR_GAPCLOSER/GapCloser ]
        do
          read -p "GapCloser is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_GAPCLOSER
        done
      echo "GapCloser ...... OK" | tee -a log.txt
fi

#check BUSCO
if [ $(which run_BUSCO.py) ]
    then
      echo "BUSCO ...... OK" | tee -a log.txt
      EXE_BUSCO=$(which run_BUSCO.py)
      DIR_BUSCO=${EXE_BUSCO%/*}
    else
      until [ -x $DIR_BUSCO/run_BUSCO.py ]
        do
          read -p "BUSCO is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/busco-3.0.2/scripts):      " DIR_BUSCO
        done
      echo "BUSCO ...... OK" | tee -a log.txt
fi
echo -e "\n" >> log.txt

#Check species name
read -p "Please input the species name as the prefix of resulting files (e.g. Escherichia_coli):      " SPECIES
echo "SPECIES=$SPECIES" >> parameters.cfg

#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done
echo "THREADS=$THREADS" >> parameters.cfg

#Check the read length of sequencing data
read -p "Please input the read length of sequencing data (e.g. 100, 150):      " READ_LENGTH
until [ $READ_LENGTH -gt 0 ]
    do
      read -p "Please input the correct integer for the read length of sequencing data (e.g. 100, 150 ):      " READ_LENGTH
    done
echo "READ_LENGTH=$READ_LENGTH" >> parameters.cfg

#Check the memory
read -p "Please input the memory size (MB) available for assembly (e.g. 20000):      " MEMORY
until [ $MEMORY -gt 0 ]
    do
      read -p "Please input the correct integer for the memory size (MB) available for assembly (e.g. 20000):      " MEMORY
    done
echo "MEMORY=$MEMORY" >> parameters.cfg

#Define the target normalization depth
read -p "Please input the target normalization depth, which may be ranged from 10 to 20 or higher for low-coverage reads (e.g. 10):      " NORMALIZATION_TARGET
until [ $NORMALIZATION_TARGET -gt 0 ]
    do
      read -p "Please input the correct integer for the target normalization depth (e.g. 10):      " NORMALIZATION_TARGET
    done
echo "NORMALIZATION_TARGET=$NORMALIZATION_TARGET" >> parameters.cfg

#Check AUGUSTUS species for BUSCO analyses
until [ -s $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES/"$AUGUSTUS_SPECIES"_parameters.cfg ]
    do
      read -p "Please input the AUGUSTUS species (they can be checked in "$AUGUSTUS_CONFIG_PATH"/species/) for BUSCO analyses (e.g. fly, honeybee1 ...):      " AUGUSTUS_SPECIES
    done
echo "AUGUSTUS_SPECIES=$AUGUSTUS_SPECIES" >> parameters.cfg

#Check Predefined BUSCO set for BUSCO analyses
until [ -s $DIR_BUSCO_LINEAGE/ancestral ]
    do
      read -p "Please input the directory of Predefined BUSCO lineage set for BUSCO analyses (e.g. /home/zf/install/busco-3.0.2/lineage/arthropoda_odb9):      " DIR_BUSCO_LINEAGE
    done
echo "DIR_BUSCO_LINEAGE=$DIR_BUSCO_LINEAGE" >> parameters.cfg


##prepare data for assembly
mkdir 0-data

#group overlapping reads into clumps and remove duplicates
echo "Clumpify and remove duplicates......" | tee -a log.txt
$DIR_BBTOOLS/clumpify.sh in1=$1 in2=$2 out1=0-data/1.clumped.fq.gz out2=0-data/2.clumped.fq.gz pigz dedupe 1>>log.txt 2>&1
echo -e "\n"  

#Quality trimming
cd 0-data/
echo "Quality trimming......" | tee -a ../log.txt
$DIR_BBTOOLS/bbduk.sh in1=1.clumped.fq.gz in2=2.clumped.fq.gz out1=1.trim.fq.gz out2=2.trim.fq.gz ziplevel=5 pigz ordered qtrim=rl trimq=15 minlen=15 ecco=t maxns=5 trimpolya=10 trimpolyg=10 trimpolyc=10 1>>../log.txt 2>&1
echo -e "\n" >> ../log.txt

#Normalize coverage by down-sampling reads over high-depth areas of a genome
if [ "$NORM" == "y" ]
  then
    echo "Normalizing coverage of raw data......" | tee -a ../log.txt
    $DIR_BBTOOLS/bbnorm.sh in1=1.trim.fq.gz in2=2.trim.fq.gz out1=1.nor.fq.gz out2=2.nor.fq.gz target="$NORMALIZATION_TARGET" min=2 histcol=2 khist=khist.txt peaks=peaks.txt 1>>../log.txt 2>&1
    echo -e "\n" >> ../log.txt
  else
    mv 1.trim.fq.gz 1.nor.fq.gz
    mv 2.trim.fq.gz 2.nor.fq.gz
fi

mv 1.nor.fq.gz 1.fq.gz && mv 2.nor.fq.gz 2.fq.gz
rm -rf *clump* *trim* *raw* *.fastq.gz
cd ..


##multi-kmer genome assembly

#Prepare the input read list
touch reads.list
echo "../0-data/1.fq.gz" >> reads.list
echo "../0-data/2.fq.gz" >> reads.list

#Prepare the kmer values used for assembling
touch kmer.list
if [ $READ_LENGTH -lt 120 ]
    then
      for n in 0 1 2 3; do echo "21+20*$n"|bc; done >> kmer.list
    else
      for n in 0 1 2 3 4 5; do echo "21+20*$n"|bc; done >> kmer.list
fi

#assembly with minia
#KMER values can be adjusted according to the genome size
mkdir 1-assembly && cd 1-assembly/
cp ../reads.list reads.list
for KMER in $(cat ../kmer.list)
do
      echo "Minia assembling at k=$KMER ......" | tee -a ../log.txt
      $DIR_MINIA/minia -in reads.list -kmer-size $KMER -abundance-min 2 -out k$KMER -max-memory $MEMORY
      test -s k$KMER.contigs.fa && echo || ($DIR_MINIA/minia -in reads.list -kmer-size $KMER -abundance-min 2 -storage-type file -out k$KMER -max-memory $MEMORY)
      test -s k$KMER.contigs.fa && echo || ($DIR_MINIA/minia -in reads.list -kmer-size `expr $KMER - 2` -abundance-min 2 -out k$KMER -max-memory $MEMORY)
      test -s k$KMER.contigs.fa && echo || ($DIR_MINIA/minia -in reads.list -kmer-size `expr $KMER - 2` -abundance-min 2 -storage-type file -out k$KMER -max-memory $MEMORY)
      test -s k$KMER.contigs.fa && (echo "Finished Minia k=$KMER" | tee -a ../log.txt) || (echo "Minia assembly with k-mer value of $KMER failed......" && exit)
      echo -e "\n"  | tee -a ../log.txt
      rm -rf dummy* *unitigs* *h5 trashme*
      cp ../reads.list reads.list
      for n in 1 2 3; do echo "./k$KMER.contigs.fa" >> reads.list; done
done
echo "Finished Multi-kmer Minia assembly" | tee -a ../log.txt
echo -e "\n"  | tee -a ../log.txt

#Generate basic assembly statistics for all contig assemblies
for kmer in $(cat ../kmer.list); do $DIR_BBTOOLS/statswrapper.sh in="k$kmer".contigs.fa format=6 >> ../assembly.statistics; done

#Reduction of heterozygous contigs
if [ "$HETER" == "y" ]
  then
    echo "Reduce heterozygous sequences with Redundans......" | tee -a ../log.txt
    KMER_MAX=$(sed -n '$p' ../kmer.list)
    python2 $DIR_REDUNDANS/redundans.py -v -f k$KMER_MAX.contigs.fa -o ./reduced -t $THREADS --log log_redundans.txt --noscaffolding --nogapclosing --identity 0.7
    cat log* >> ../log.txt
    mv reduced/scaffolds.reduced.fa k$KMER_MAX.contigs.reduced.fa && rm -rf reduced *fa.fai log* reads.list
    echo -e "\n"  | tee -a ../log.txt
    CONTIGS="k$KMER_MAX.contigs.reduced.fa"
  else
    KMER_MAX=$(sed -n '$p' ../kmer.list)
    CONTIGS="k$KMER_MAX.contigs.fa"
    echo -e "\n"  | tee -a ../log.txt
fi

##Scaffolding assembled contigs
echo "Scaffolding......" | tee -a ../log.txt
#Generate the mapping file in sorted, indexed bam format
cd .. && mkdir 2-scaffolding && cd 2-scaffolding/
$DIR_MINIMAP2/minimap2 -ax sr ../1-assembly/$CONTIGS ../0-data/1.fq.gz ../0-data/2.fq.gz -t $THREADS | $DIR_SAMTOOLS/samtools view -bS -@ $THREADS | $DIR_SAMTOOLS/samtools sort -@ $THREADS > map.bam
$DIR_SAMTOOLS/samtools index map.bam -@ $THREADS

#Scaffolding with BESST
$DIR_BESST/runBESST -c ../1-assembly/$CONTIGS -f map.bam -o ./ -orientation fr --iter 10000
test -s BESST_output/pass1/Scaffolds_pass1.fa && echo || (read -p "Library parameters cannot be assessed by BESST, please input the mean insert size (bp) and standard deviation (bp) of libraries, for example 300 50:      " m s && echo "m=$m" "s=$s" >> parameters.cfg && $DIR_BESST/runBESST -c ../1-assembly/$CONTIGS -f map.bam -o ./ -orientation fr --iter 10000 -m $m -s $s)
test -s BESST_output/pass1/Scaffolds_pass1.fa && echo || (read -p "Input library parameters may be wrong, please input the mean insert size (bp) and standard deviation (bp) of libraries again, for example 300 50:      " m s && echo "m=$m" "s=$s" >> parameters.cfg && $DIR_BESST/runBESST -c ../1-assembly/$CONTIGS -f map.bam -o ./ -orientation fr --iter 10000 -m $m -s $s)
test -s BESST_output/pass1/Scaffolds_pass1.fa && (echo -e "Scaffolding has been finished\n" | tee -a ../log.txt) || (echo "Scaffolding with BESST failed. Please try other tools for scaffolding or use the assembly in directory 1-assembly for BUSCO analyses" && exit)
echo -e "\n"  | tee -a ../log.txt
cd ..


##Gap filling
echo "Gap filling......" | tee -a log.txt
mkdir 3-gapclosing && cd 3-gapclosing/
#Prepare config file forGapCloser
touch gapcloser.config
echo "[LIB]" >> gapcloser.config
echo "q1=../0-data/1.fq.gz" >> gapcloser.config
echo "q2=../0-data/2.fq.gz" >> gapcloser.config

#Gap filling with GapCloser
$DIR_GAPCLOSER/GapCloser -a ../2-scaffolding/BESST_output/pass1/Scaffolds_pass1.fa -b gapcloser.config -o scaffolds.gapcloser.fa -l $READ_LENGTH -t $THREADS
test -s scaffolds.gapcloser.fa && (echo -e "Gap filling has been finished\n" | tee -a ../log.txt) || (echo "Gap filling with GapCloser failed. Please try other tools for gap filling or use the assemblies generated in previous steps for BUSCO analyses" && exit)


#filter the scaffolds shorter than 1000 bp
$DIR_BBTOOLS/reformat.sh in=scaffolds.gapcloser.fa out=../scaffolds_"$SPECIES".fa minlength=1000
cd ..

#Generate basic assembly statistics
$DIR_BBTOOLS/statswrapper.sh in=1-assembly/$CONTIGS,2-scaffolding/BESST_output/pass1/Scaffolds_pass1.fa,scaffolds_$SPECIES.fa >> assembly.statistics

if [ $READ_LENGTH -lt 120 ]
    then
      for n in 9 7 5 3; do sed -i ""$n"d" assembly.statistics; done
    else
      for n in 13 11 9 7 5 3; do sed -i ""$n"d" assembly.statistics; done
fi

rm -rf 2-scaffolding/map* *list


##BUSCO assessment
mkdir 4-busco && cd 4-busco/
python $DIR_BUSCO/run_BUSCO.py -i ../scaffolds_$SPECIES.fa -c $THREADS -l $DIR_BUSCO_LINEAGE -m geno -o $SPECIES -sp $AUGUSTUS_SPECIES
