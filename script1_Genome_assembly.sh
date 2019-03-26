#!/bin/bash

#Most executables are recommended to be added into the environmental paths
#Tools pigz, BBTools, Lighter, Minia, redundans, Minimap2, samtools, BESST, GapCloser and BUSCO may be used in this script

#modify the names of input raw paired read files as 1.raw.fq.gz and 2.raw.fq.gz in the working folder
#the default starting kmer value is 21 and thus kmer values are 21, 21+20, 21+2*20....
#all the statistics of input reads and assemblies generated in the whole asembly progress are summerized respectively /bbtool/bbtool.log and /minia/assembly.log


#Define variables
SPECIES="Escherichia_coli"
THREADS="8"
READ_LENGTH="100" #bp
GENOME_SIZE="4000000"  #value usually slightly smaller than the real genome size
MEMORY="10000" #Mb
NORMALIZATION_TARGET="10"  #may be ranged from 10 to 20 or higher for low-coverage reads
DIR_REDUNDANS="/home/zf/install/redundans"
DIR_BUSCO="/home/zf/install/busco-3.0.2"
DIR_BUSCO_LINEAGE="/home/zf/install/busco-3.0.2/lineage/arthropoda_odb9"  #Predefined BUSCO set can be downloaded from BUSCO website
AUGUSTUS_SPECIES="fly"  #species parameter can be checked in /home/zf/install/Augustus-3.3.2/config/species

mkdir bbtool && cd bbtool/

#group overlapping reads into clumps and remove duplicates
clumpify.sh in1=../1.raw.fq.gz in2=../2.raw.fq.gz out1=1.clumped.fq.gz out2=2.clumped.fq.gz pigz dedupe

#Quality trimming
bbduk.sh in1=1.clumped.fq.gz in2=2.clumped.fq.gz out1=1.trim.fq.gz out2=2.trim.fq.gz ziplevel=5 pigz ordered qtrim=rl trimq=15 minlen=15 ecco=t maxns=5 trimpolya=10 1>>bbtool.log 2>&1

#Normalize coverage by down-sampling reads over high-depth areas of a genome
bbnorm.sh in1=1.trim.fq.gz in2=2.trim.fq.gz out1=1.nor.fq.gz out2=2.nor.fq.gz target="$NORMALIZATION_TARGET" min=2 histcol=2 khist=khist.txt peaks=peaks.txt 1>>bbtool.log 2>&1

#Error correction
lighter -r 1.nor.fq.gz -r 2.nor.fq.gz -K 27 $GENOME_SIZE -t $THREADS -zlib 5

mv 1.nor.cor.fq.gz 1.fq.gz && mv 2.nor.cor.fq.gz 2.fq.gz
rm -rf *clump* *trim* *nor* *raw* *.fastq.gz && cd ..

#Prepare the input read list
touch reads.list
echo "../bbtool/1.fq.gz" >> reads.list
echo "../bbtool/2.fq.gz" >> reads.list

#Prepare the kmer values used for assembling
touch kmer.list
if [ $READ_LENGTH -lt 120 ]
    then
      for n in 0 1 2 3; do echo "21+20*$n"|bc; done >> kmer.list
    else
      for n in 0 1 2 3 4 5; do echo "21+20*$n"|bc; done >> kmer.list
fi

#multi-kmer genome assembly
#KMER values can be adjusted according to the genome size
mkdir minia && cd minia/
cp ../reads.list reads.list
for KMER in $(cat ../kmer.list)
do
      echo "Minia assembling at k=$KMER"
      minia -in reads.list -kmer-size $KMER -abundance-min 2 -out k$KMER -max-memory $MEMORY
      echo "Finished Minia k=$KMER"
      rm -rf dummy* *unitigs* *h5 trashme*
      cp ../reads.list reads.list
      for n in 1 2 3; do echo "./k$KMER.contigs.fa" >> reads.list; done
done
echo "Finished Multi-kmer Minia assembly"

#Generate basic assembly statistics for all contig assemblies
for kmer in $(cat ../kmer.list); do statswrapper.sh in="k$kmer".contigs.fa format=6 >> assembly.log; done

#Reduction of heterozygous contigs
KMER_MAX=$(sed -n '$p' ../kmer.list)
python2 $DIR_REDUNDANS/redundans.py -v -f k$KMER_MAX.contigs.fa -o ./reduced -t $THREADS --log log_redundans.txt --noscaffolding --nogapclosing --identity 0.7

mv reduced/scaffolds.reduced.fa k$KMER_MAX.contigs.reduced.fa && rm -rf reduced *fa.fai

#Generate the mapping file in sorted, indexed bam format
cd .. && mkdir mapping && cd mapping/
minimap2 -ax sr ../minia/k$KMER_MAX.contigs.reduced.fa ../bbtool/1.fq.gz ../bbtool/2.fq.gz -t $THREADS | samtools view -bS -@ $THREADS | samtools sort -@ $THREADS > map.bam
samtools index map.bam -@ $THREADS
cd ..

#Scaffolding assembled contigs
runBESST -c ./minia/k$KMER_MAX.contigs.reduced.fa -f ./mapping/map.bam -o BESST -orientation fr --iter 10000
test -s BESST/BESST_output/pass1/Scaffolds_pass1.fa && echo -e "Scaffolding has been finished\n" || (read -p "Library parameters cannot be assessed by BESST, please input the mean insert size (bp) and standard deviation (bp) of libraries, for example 300 50:      " m s && runBESST -c ./minia/k$KMER_MAX.contigs.reduced.fa -f ./mapping/map.bam -o BESST -orientation fr --iter 10000 -m $m -s $s)

#Prepare config file forGapCloser
touch gapcloser.config
echo "[LIB]" >> gapcloser.config
echo "q1=./bbtool/1.fq.gz" >> gapcloser.config
echo "q2=./bbtool/2.fq.gz" >> gapcloser.config

#Gap filling
GapCloser -a BESST/BESST_output/pass1/Scaffolds_pass1.fa -b gapcloser.config -o scaffolds.gapcloser.fa -l $READ_LENGTH -t $THREADS

#filter the scaffolds shorter than 1000 bp
reformat.sh in=scaffolds.gapcloser.fa out=scaffolds_"$SPECIES".fa minlength=1000

#Generate basic assembly statistics
statswrapper.sh in=./minia/k$KMER_MAX.contigs.reduced.fa,./BESST/BESST_output/pass1/Scaffolds_pass1.fa,scaffolds_$SPECIES.fa >> ./minia/assembly.log

if [ $READ_LENGTH -lt 120 ]
    then
      for n in 9 7 5 3; do sed -i ""$n"d" minia/assembly.log; done
    else
      for n in 13 11 9 7 5 3; do sed -i ""$n"d" minia/assembly.log; done
fi

rm -rf mapping BESST gapcloser.config reads.list scaffolds.gapcloser*
#rm -rf bbtool/*fq.gz minia/*fa

#BUSCO assessment
mkdir busco && cd busco/
python $DIR_BUSCO/scripts/run_BUSCO.py -i ../scaffolds_$SPECIES.fa -c $THREADS -l $DIR_BUSCO_LINEAGE -m geno -o $SPECIES -sp $AUGUSTUS_SPECIES
