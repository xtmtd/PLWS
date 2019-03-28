#!/bin/bash
#2019.03.28 by ZF

#Type 'sh script3_UCE_probe_design assemblies_folder', e.g. sh script3_UCE_probe_design.sh assemblies
#Copy all the genome assemblies to a folder (e.g. assemblies/) and replace assembly name using "SPECIES_NAME.fa"
#Tools art, BBtools, faToTwoBit and stampy, as well as PHYLUCE conda enviorment, are used and will be automatically checked prior to formal analyses in this script


##check the packages
echo "Checking the package dependency......" | tee -a log.txt
echo -e "\n" >> log.txt

#Check the phyluce conda environment
read -p "Please input the name of PHYLUCE conda environment ( e.g. phyluce):      " ENV
source activate $ENV
if [ $(which phyluce_probe_strip_masked_loci_from_set) ]
    then
      echo "PHYLUCE conda environment ...... OK" | tee -a log.txt
    else
      until [ $(which phyluce_probe_strip_masked_loci_from_set) ]
        do
          read -p "PHYLUCE conda environment is not found. Please input the correct conda environment name ( e.g. phyluce):      " ENV
          source activate $ENV
        done
      echo "PHYLUCE conda environment ...... OK" | tee -a log.txt
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

#check ART
if [ $(which art_illumina) ]
    then
      echo "ART ...... OK" | tee -a log.txt
      EXE_ART=$(which art_illumina)
      DIR_ART=${EXE_ART%/*}
    else
      until [ -x $DIR_ART/art_illumina ]
        do
          read -p "ART is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/art):      " DIR_ART
        done
      echo "ART ...... OK" | tee -a log.txt
fi

#check faToTwoBit
if [ $(which faToTwoBit) ]
    then
      echo "faToTwoBit ...... OK" | tee -a log.txt
      EXE_FATOTWOBIT=$(which faToTwoBit)
      DIR_FATOTWOBIT=${EXE_FATOTWOBIT%/*}
    else
      until [ -x $DIR_FATOTWOBIT/faToTwoBit ]
        do
          read -p "faToTwoBit is not found. Please input its installation directory (absolute path, e.g. /usr/local/bin):      " DIR_FATOTWOBIT
        done
      echo "faToTwoBit ...... OK" | tee -a log.txt
fi

#check stampy
if [ $(which stampy.py) ]
    then
      echo "stampy ...... OK" | tee -a log.txt
      EXE_STAMPY=$(which stampy.py)
      DIR_STAMPY=${EXE_STAMPY%/*}
    else
      until [ -x $DIR_STAMPY/stampy.py ]
        do
          read -p "stampy is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/stampy-1.0.32):      " DIR_STAMPY
        done
      echo "stampy ...... OK" | tee -a log.txt
fi

#Check the name of base species
read -p "Please input the name of the base species for UCE probe design (e.g. Mesaphorura_yosii):      " base
echo "base=$base" >> parameters.cfg

#Check the name of the lineage
read -p "Please input the name of the lineage for the prefix of UCE probe (e.g. Collembola, Diptera):      " GROUP
echo "GROUP=$GROUP" >> parameters.cfg

#Check the threads to be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done
echo "THREADS=$THREADS" >> parameters.cfg

#Check the simulated coverage
read -p "Please input the coverage to be simulated for each species, usually 2-5. It can be increased when more loci are required or the species are very divergent (e.g. 2):      " SIMULATE_COV
until [ $SIMULATE_COV -gt 0 ]
    do
      read -p "Please input the the correct integer for the coverage to be simulated, usually 2-5 (e.g. 2):      " SIMULATE_COV
    done
echo "SIMULATE_COV=$SIMULATE_COV" >> parameters.cfg

#Check the substitution rate for sequence aligning
read -p "Please input the substitution rate for sequence stampy aligning. Value of 0.05 may be OK for most cases, but for the species very divergent, rate value can be increased to 0.1 (e.g. 0.05):      " SUBSTITUTION_RATE
echo "SUBSTITUTION_RATE=$SUBSTITUTION_RATE" >> parameters.cfg

echo -e "\n" | tee -a log.txt

#copy the assemblies to 0-assemblies/
echo "copy the assemblies to 0-assemblies/" | tee -a log.txt
cp -r $1 0-genomes


#Generate a species list file (species.list with base species excluded) in the initial working folder
cd 0-genomes
for file in * ; do echo ${file%.*} >> ../species.list; done
sed -i "/"$base"/d" ../species.list
cd ..
SPECIES_NAME=$(cat species.list)


#Generate a python script simplify_headers.py to simplify the sequence head
echo "simplify the sequence head..." | tee -a log.txt
echo '#!/usr/bin/python' >> simplify_headers.py
echo 'from Bio import SeqIO' >> simplify_headers.py
for SPECIES in $base $SPECIES_NAME
  do
    echo -e '\n' >> simplify_headers.py
    echo "# $SPECIES" >> simplify_headers.py
    echo "with open(\"$SPECIES.fa\", "'"rU"'") as infile:" >> simplify_headers.py
    echo "  with open(\"$SPECIES.fasta\", \"w\") as outf:" >> simplify_headers.py
    echo "    for seq in SeqIO.parse(infile, 'fasta'):" >> simplify_headers.py
    echo "      seq.name = \"\"" >> simplify_headers.py
    echo "      seq.description = \"\"" >> simplify_headers.py
    echo "      outf.write(seq.format('fasta'))" >> simplify_headers.py
  done

#Simplify the FASTA heads
cd 0-genomes/
python ../simplify_headers.py
rm *.fa

#Move each genome assembly to its own directory
for SPECIES in *; do mkdir ${SPECIES%.*}; mv $SPECIES ${SPECIES%.*}; done

#Convert genomes to 2bit format
echo "Convert genomes to 2bit format..." | tee -a log.txt
for SPECIES in *; do $DIR_FATOTWOBIT/faToTwoBit $SPECIES/$SPECIES.fasta $SPECIES/${SPECIES%.*}.2bit; done

#Simulate interleaved reads from genome assemblies (base genome "Pediculus_humanus" excluded)
cd .. && mkdir 1-reads && cd 1-reads/
for SPECIES in $SPECIES_NAME
  do
    echo "Simulate interleaved reads from genome assembly of $SPECIES..." | tee -a ../log.txt
    $DIR_ART/art_illumina --paired --in ../0-genomes/$SPECIES/$SPECIES.fasta --out $SPECIES-pe100-reads --len 100 --fcov $SIMULATE_COV --mflen 200 --sdev 50 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
    reformat.sh in1=$SPECIES-pe100-reads1.fq in2=$SPECIES-pe100-reads2.fq out=$SPECIES-pe100-reads.fq.gz 1>>../log.txt 2>&1
    echo -e "\n" | tee -a ../log.txt
    rm $SPECIES-pe100-reads1.fq $SPECIES-pe100-reads2.fq
  done

#Prepare base genome prior to alignment
echo "Prepare base genome..." | tee -a ../log.txt
cd .. && mkdir 2-base && cd 2-base/
cp ../0-genomes/$base/$base.fasta ./
python2 $DIR_STAMPY/stampy.py --species="$base" --assembly="$base" -G $base $base.fasta
python2 $DIR_STAMPY/stampy.py -g $base -H $base --logfile log

cat log >> ../log.txt && rm log
echo -e "\n" | tee -a ../log.txt

#Align simulated reads to the base genome
cd .. && mkdir 3-alignments && cd 3-alignments/
for SPECIES in $SPECIES_NAME
  do
    echo "Aligning sequencing reads of $SPECIES to the base genome..." | tee -a ../log.txt
    mkdir $SPECIES
    cd $SPECIES
    python2 $DIR_STAMPY/stampy.py --maxbasequal 93 -g ../../2-base/$base -h ../../2-base/$base --substitutionrate="$SUBSTITUTION_RATE" -t $THREADS --insertsize=200 -M ../../1-reads/$SPECIES-pe100-reads.fq.gz | samtools view -Sb -@ $THREADS - > $SPECIES-to-$base.bam
    echo "mapping results of $SPECIES to the base genome" >> ../../log.txt
    samtools flagstat $SPECIES-to-$base.bam -@ $THREADS >> ../../log.txt
    cd ..
    echo -e "\n" | tee -a ../log.txt
  done

#Remove those unmapped reads and symlink all of the reduced BAM files
echo "Remove those unmapped reads and symlink all of the reduced BAM files..." | tee -a ../log.txt
mkdir all
for SPECIES in $SPECIES_NAME
  do
    samtools view -h -F 4 -b $SPECIES/$SPECIES-to-$base.bam -@ $THREADS > $SPECIES/$SPECIES-to-$base-MAPPING.bam
    rm $SPECIES/$SPECIES-to-$base.bam
    ln -s ../$SPECIES/$SPECIES-to-$base-MAPPING.bam all/$SPECIES-to-$base-MAPPING.bam
  done

#Convert bam to sorted bed
echo "Convert bam to sorted bed..." | tee -a ../log.txt
cd .. && mkdir 4-bed && cd 4-bed/
for i in ../3-alignments/all/*.bam; do echo $i; bedtools bamtobed -i $i -bed12 > `basename $i`.bed; done
for i in *.bed; do echo $i; bedtools sort -i $i > ${i%.*}.sort.bed; done

#Merge overlapping or nearly-overlapping intervals
echo "Merge overlapping intervals..." | tee -a ../log.txt
for i in *.bam.sort.bed; do echo $i; bedtools merge -i $i > ${i%.*}.merge.bed; done
for i in *.bam.sort.merge.bed; do wc -l $i; done

#Remove repetitive intervals
echo "Remove repetitive intervals..." | tee -a ../log.txt
for i in *.sort.merge.bed; do phyluce_probe_strip_masked_loci_from_set --bed $i --twobit ../0-genomes/$base/$base.2bit --output ${i%.*}.strip.bed --filter-mask 0.25 --min-length 80; done

#Generate a configure file for bed files
echo "Generate a configure file for bed files..." | tee -a ../log.txt
echo '[beds]' >> ../bed-files.conf
for SPECIES in $SPECIES_NAME; do echo "$SPECIES:$SPECIES-to-$base-MAPPING.bam.sort.merge.strip.bed" >> ../bed-files.conf; done

#Determining shared conserved loci
echo -e "\n" | tee -a ../log.txt
echo "Determining shared conserved loci..." | tee -a ../log.txt
phyluce_probe_get_multi_merge_table --conf ../bed-files.conf --base-taxon $base --output $GROUP-to-$base.sqlite
phyluce_probe_query_multi_merge_table --db $GROUP-to-$base.sqlite --base-taxon $base | tee -a ../log.txt

#Check the species count to be used
read -p "Please input the number of taxa to determine shared conserved loci. The maximum value is the number of species (base species excluded) (e.g. 10):      " SPECIFIC_COUNT
until [ $SPECIFIC_COUNT -gt 0 ]
    do
      read -p "Please input the correct integer for the number of taxa (e.g. 10):      " SPECIFIC_COUNT
    done
echo "SPECIFIC_COUNT=$SPECIFIC_COUNT" >> ../parameters.cfg

echo -e "\n" | tee -a ../log.txt
phyluce_probe_query_multi_merge_table --db $GROUP-to-$base.sqlite --base-taxon $base --output "$base"+"$SPECIFIC_COUNT".bed --specific-counts $SPECIFIC_COUNT | tee -a ../log.txt

#Extract FASTA sequence from base genome for temp bait design
echo -e "\n" | tee -a ../log.txt
echo "Extract FASTA sequence from base genome for temp bait design..." | tee -a ../log.txt
phyluce_probe_get_genome_sequences_from_bed --bed "$base"+"$SPECIFIC_COUNT".bed --twobit ../0-genomes/$base/$base.2bit --buffer-to 160 --output "$base"+"$SPECIFIC_COUNT".fasta | tee -a ../log.txt

#Design a temporary bait set from the base taxon
echo -e "\n" | tee -a ../log.txt
echo "Design a temporary bait set from the base taxon..." | tee -a ../log.txt
phyluce_probe_get_tiled_probes --input "$base"+"$SPECIFIC_COUNT".fasta --probe-prefix "uce-" --design $GROUP-v1 --designer zf --tiling-density 3 --two-probes --overlap middle --masking 0.25 --remove-gc --output "$base"+"$SPECIFIC_COUNT".temp.probes | tee -a ../log.txt

#Remove duplicates from our temporary bait set
echo -e "\n" | tee -a ../log.txt
echo "Remove duplicates from the temporary bait set..." | tee -a ../log.txt
phyluce_probe_easy_lastz --target "$base"+"$SPECIFIC_COUNT".temp.probes --query "$base"+"$SPECIFIC_COUNT".temp.probes --identity 50 --coverage 50 --output "$base"+"$SPECIFIC_COUNT".temp.probes-TO-SELF-PROBES.lastz | tee -a ../log.txt
echo -e "\n" | tee -a ../log.txt
phyluce_probe_remove_duplicate_hits_from_probes_using_lastz --fasta "$base"+"$SPECIFIC_COUNT".temp.probes --lastz "$base"+"$SPECIFIC_COUNT".temp.probes-TO-SELF-PROBES.lastz --probe-prefix=uce- | tee -a ../log.txt

#Align duplicate-free temporary bait against exemplar genomes
echo -e "\n" | tee -a ../log.txt
echo "Align duplicate-free temporary bait against exemplar genomes..." | tee -a ../log.txt
cd .. && mkdir -p 5-probe_design/$GROUP-genome-lastz && cd 5-probe_design/
echo "$base" >> ../species.list
phyluce_probe_run_multiple_lastzs_sqlite --probefile ../4-bed/"$base"+"$SPECIFIC_COUNT".temp-DUPE-SCREENED.probes --scaffoldlist $(cat ../species.list) --genome-base-path ../0-genomes --identity 50 --cores $THREADS --db "$base"+"$SPECIFIC_COUNT".sqlite --output $GROUP-genome-lastz

#Generate the configure file of the genome data
echo '[scaffolds]' >> ../$GROUP-genome.conf
for SPECIES in $SPECIES_NAME; do echo "$SPECIES:../0-genomes/$SPECIES/$SPECIES.2bit" >> ../$GROUP-genome.conf; done
echo "$base:../0-genomes/$base/$base.2bit" >> ../$GROUP-genome.conf

#Extract sequence around conserved loci from exemplar genomes
echo -e "\n" | tee -a ../log.txt
echo "Extract sequence around conserved loci from exemplar genomes..." | tee -a ../log.txt
phyluce_probe_slice_sequence_from_genomes --conf ../$GROUP-genome.conf --lastz $GROUP-genome-lastz --probes 180 --name-pattern ""$base"+"$SPECIFIC_COUNT".temp-DUPE-SCREENED.probes_v_{}.lastz.clean" --output $GROUP-genome-fasta | tee -a ../log.txt

#Detect loci consistently across all of the exemplar taxa
echo -e "\n" | tee -a ../log.txt
echo "Detect loci consistently across all of the exemplar taxa..." | tee -a ../log.txt
phyluce_probe_get_multi_fasta_table --fastas $GROUP-genome-fasta --output multifastas.sqlite --base-taxon $base | tee -a ../log.txt

#Get the distribution of hits among exemplar taxa
echo -e "\n" | tee -a ../log.txt
echo "Get the distribution of hits among exemplar taxa..." | tee -a ../log.txt
phyluce_probe_query_multi_fasta_table --db multifastas.sqlite --base-taxon $base | tee -a ../log.txt

#Determine how many loci will be kept in final probe
echo -e "\n" | tee -a ../log.txt
echo "Determine the number of loci kept in final probe..." | tee -a ../log.txt
#Check the final species count to be used
read -p "Please input the number of taxa to determine how many loci will be kept in final probe. The value is not less than previous "SPECIFIC_COUNT" (e.g. 10):      " FINAL_SPECIFIC_COUNT
until [ $FINAL_SPECIFIC_COUNT -gt 0 ]
    do
      read -p "Please input the correct integer for the number of taxa (e.g. 10):      " FINAL_SPECIFIC_COUNT
    done
echo "FINAL_SPECIFIC_COUNT=$FINAL_SPECIFIC_COUNT" >> ../parameters.cfg

phyluce_probe_query_multi_fasta_table --db multifastas.sqlite --base-taxon $base --output "$base"+"$SPECIFIC_COUNT"-back-to-"$FINAL_SPECIFIC_COUNT".conf --specific-counts "$FINAL_SPECIFIC_COUNT" | tee -a ../log.txt

#Design the final bait set
echo -e "\n" | tee -a ../log.txt
echo "Design the final bait set..." | tee -a ../log.txt
phyluce_probe_get_tiled_probe_from_multiple_inputs --fastas $GROUP-genome-fasta --multi-fasta-output "$base"+"$SPECIFIC_COUNT"-back-to-"$FINAL_SPECIFIC_COUNT".conf --probe-prefix "uce-" --designer zf --design $GROUP-v1 --tiling-density 3 --overlap middle --masking 0.25 --remove-gc --two-probes --output $GROUP-v1-master-probe-list.fasta | tee -a ../log.txt

#Remove duplicates from bait set
echo -e "\n" | tee -a ../log.txt
echo "Remove duplicates from bait set..." | tee -a ../log.txt
phyluce_probe_easy_lastz --target $GROUP-v1-master-probe-list.fasta --query $GROUP-v1-master-probe-list.fasta --identity 50 --coverage 50 --output $GROUP-v1-master-probe-list-TO-SELF-PROBES.lastz | tee -a ../log.txt
phyluce_probe_remove_duplicate_hits_from_probes_using_lastz --fasta $GROUP-v1-master-probe-list.fasta --lastz $GROUP-v1-master-probe-list-TO-SELF-PROBES.lastz --probe-prefix=uce- | tee -a ../log.txt

echo "The final bait set is /probe_design/"$GROUP"-v1-master-probe-list-DUPE-SCREENED.fasta" | tee -a ../log.txt

source deactivate
