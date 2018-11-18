#!/bin/bash

#Prepare a species list file (species.list with base species Pediculus_humanus excluded) in the initial working folder
#Tools samtools, art, BBtools, faToTwoBit and stampy are used in this script and the former two have been installed in environmental paths

#Define variables
DIR_ASSEMBLY="/home/zf/test/phylogenomics/assembly"
DIR_fa2bit="/home/zf/test/phylogenomics/uce/probe"
DIR_STAMPY="/home/zf/stampy-1.0.32"
SPECIES_NAME=$(cat species.list)
base="Pediculus_humanus"
GROUP="Phthiraptera"
THREADS="8"

#Initiate the phyluce environment
source activate phyluce

mkdir genomes && cd genomes/

#Copy all the genome assemblies to genomes/ and replace assembly name using "species_name.fa"
for SPECIES in $SPECIES_NAME; do cp $DIR_ASSEMBLY/$SPECIES/scaffolds_$SPECIES.fa ./; mv scaffolds_$SPECIES.fa $SPECIES.fa; done

#Generate a python script simplify_headers.py to simplify the sequence head
echo '#!/usr/bin/python' >> ../simplify_headers.py
echo 'from Bio import SeqIO' >> ../simplify_headers.py
for SPECIES in $base $SPECIES_NAME
  do
    echo -e '\n' >> ../simplify_headers.py
    echo "# $SPECIES" >> ../simplify_headers.py
    echo "with open(\"$SPECIES.fa\", "'"rU"'") as infile:" >> ../simplify_headers.py
    echo "  with open(\"$SPECIES.fasta\", \"w\") as outf:" >> ../simplify_headers.py
    echo "    for seq in SeqIO.parse(infile, 'fasta'):" >> ../simplify_headers.py
    echo "      seq.name = \"\"" >> ../simplify_headers.py
    echo "      seq.description = \"\"" >> ../simplify_headers.py
    echo "      outf.write(seq.format('fasta'))" >> ../simplify_headers.py
  done

#Simplify the FASTA heads
python ../simplify_headers.py
rm *.fa

#Move each genome assembly to its own directory
for SPECIES in *; do mkdir ${SPECIES%.*}; mv $SPECIES ${SPECIES%.*}; done

#Convert genomes to 2bit format
for SPECIES in *; do $DIR_fa2bit/faToTwoBit $SPECIES/$SPECIES.fasta $SPECIES/${SPECIES%.*}.2bit; done

#Simulate interleaved reads from genome assemblies (base genome "Pediculus_humanus" excluded)
cd .. && mkdir reads && cd reads/
for SPECIES in $SPECIES_NAME
  do
    art_illumina --paired --in ../genomes/$SPECIES/$SPECIES.fasta --out $SPECIES-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 50 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
    reformat.sh in1=$SPECIES-pe100-reads1.fq in2=$SPECIES-pe100-reads2.fq out=$SPECIES-pe100-reads.fq.gz 1>>simulation.log 2>&1
    rm $SPECIES-pe100-reads1.fq $SPECIES-pe100-reads2.fq
  done

#Prepare base genome prior to alignment
cd .. && mkdir base && cd base/
cp ../genomes/$base/$base.fasta ./
$DIR_STAMPY/stampy.py --species="$base" --assembly="$base" -G $base $base.fasta
$DIR_STAMPY/stampy.py -g $base -H $base

#Align simulated reads to the base genome
cd .. && mkdir alignments && cd alignments/
for SPECIES in $SPECIES_NAME
  do
    echo "aligning sequencing reads of $SPECIES to the base genome"
    export reads=$SPECIES-pe100-reads.fq.gz
    mkdir $SPECIES
    cd $SPECIES
    $DIR_STAMPY/stampy.py --maxbasequal 93 -g ../../base/$base -h ../../base/$base --substitutionrate=0.05 -t $THREADS --insertsize=200 -M ../../reads/$reads | samtools view -Sb -@ $THREADS - > $SPECIES-to-$base.bam
    echo "mapping status of $SPECIES to the base genome" >> ../mapping.log
    samtools flagstat $SPECIES-to-$base.bam -@ $THREADS >> ../mapping.log
    cd ..
  done

#Remove those unmapped reads and symlink all of the reduced BAM files
mkdir all
for SPECIES in $SPECIES_NAME
  do
    samtools view -h -F 4 -b $SPECIES/$SPECIES-to-$base.bam -@ $THREADS > $SPECIES/$SPECIES-to-$base-MAPPING.bam
    rm $SPECIES/$SPECIES-to-$base.bam
    ln -s ../$SPECIES/$SPECIES-to-$base-MAPPING.bam all/$SPECIES-to-$base-MAPPING.bam
  done

#Convert bam to sorted bed
cd .. && mkdir bed && cd bed/
for i in ../alignments/all/*.bam; do echo $i; bedtools bamtobed -i $i -bed12 > `basename $i`.bed; done
for i in *.bed; do echo $i; bedtools sort -i $i > ${i%.*}.sort.bed; done

#Merge overlapping or nearly-overlapping intervals
for i in *.bam.sort.bed; do echo $i; bedtools merge -i $i > ${i%.*}.merge.bed; done
for i in *.bam.sort.merge.bed; do wc -l $i; done

#Remove repetitive intervals
for i in *.sort.merge.bed; do phyluce_probe_strip_masked_loci_from_set --bed $i --twobit ../genomes/$base/$base.2bit --output ${i%.*}.strip.bed --filter-mask 0.25 --min-length 80; done

#Generate a configure file for bed files
echo '[beds]' >> ../bed-files.conf
for SPECIES in $SPECIES_NAME; do echo "$SPECIES:$SPECIES-to-$base-MAPPING.bam.sort.merge.strip.bed" >> ../bed-files.conf; done

#Determining shared conserved loci
phyluce_probe_get_multi_merge_table --conf ../bed-files.conf --base-taxon $base --output $GROUP-to-$base.sqlite
phyluce_probe_query_multi_merge_table --db $GROUP-to-$base.sqlite --base-taxon $base
phyluce_probe_query_multi_merge_table --db $GROUP-to-$base.sqlite --base-taxon $base --output $base+9.bed --specific-counts 9

#Extract FASTA sequence from base genome for temp bait design
phyluce_probe_get_genome_sequences_from_bed --bed $base+9.bed --twobit ../genomes/$base/$base.2bit --buffer-to 160 --output $base+9.fasta

#Design a temporary bait set from the base taxon
phyluce_probe_get_tiled_probes --input $base+9.fasta --probe-prefix "uce-" --design $GROUP-v1 --designer zf --tiling-density 3 --two-probes --overlap middle --masking 0.25 --remove-gc --output $base+9.temp.probes

#Remove duplicates from our temporary bait set
phyluce_probe_easy_lastz --target $base+9.temp.probes --query $base+9.temp.probes --identity 50 --coverage 50 --output $base+9.temp.probes-TO-SELF-PROBES.lastz
phyluce_probe_remove_duplicate_hits_from_probes_using_lastz --fasta $base+9.temp.probes --lastz $base+9.temp.probes-TO-SELF-PROBES.lastz --probe-prefix=uce-

#Align duplicate-free temporary bait against exemplar genomes
cd .. && mkdir -p probe_design/$GROUP-genome-lastz && cd probe_design/
echo "$base" >> ../species.list
phyluce_probe_run_multiple_lastzs_sqlite --probefile ../bed/$base+9.temp-DUPE-SCREENED.probes --scaffoldlist $(cat ../species.list) --genome-base-path ../genomes --identity 50 --cores $THREADS --db $base+9+Osborniella.sqlite --output $GROUP-genome-lastz

#Generate the configure file of the genome data
echo '[scaffolds]' >> ../$GROUP-genome.conf
for SPECIES in $SPECIES_NAME; do echo "$SPECIES:../genomes/$SPECIES/$SPECIES.2bit" >> ../$GROUP-genome.conf; done
echo "$base:../genomes/$base/$base.2bit" >> ../$GROUP-genome.conf

#Extract sequence around conserved loci from exemplar genomes
phyluce_probe_slice_sequence_from_genomes --conf ../$GROUP-genome.conf --lastz $GROUP-genome-lastz --probes 180 --name-pattern "$base+9.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" --output $GROUP-genome-fasta

#Detect loci consistently across all of the exemplar taxa
phyluce_probe_get_multi_fasta_table --fastas $GROUP-genome-fasta --output multifastas.sqlite --base-taxon $base

#Get the distribution of hits among exemplar taxa
phyluce_probe_query_multi_fasta_table --db multifastas.sqlite --base-taxon $base

#Determine how many loci will be kept in final probe
phyluce_probe_query_multi_fasta_table --db multifastas.sqlite --base-taxon $base --output $base+9-back-to-10.conf --specific-counts 10

#Design the final bait set
phyluce_probe_get_tiled_probe_from_multiple_inputs --fastas $GROUP-genome-fasta --multi-fasta-output $base+9-back-to-10.conf --probe-prefix "uce-" --designer zf --design $GROUP-v1 --tiling-density 3 --overlap middle --masking 0.25 --remove-gc --two-probes --output $GROUP-v1-master-probe-list.fasta

#Remove duplicates from bait set
phyluce_probe_easy_lastz --target $GROUP-v1-master-probe-list.fasta --query $GROUP-v1-master-probe-list.fasta --identity 50 --coverage 50 --output $GROUP-v1-master-probe-list-TO-SELF-PROBES.lastz
phyluce_probe_remove_duplicate_hits_from_probes_using_lastz --fasta $GROUP-v1-master-probe-list.fasta --lastz $GROUP-v1-master-probe-list-TO-SELF-PROBES.lastz --probe-prefix=uce-

source deactivate

