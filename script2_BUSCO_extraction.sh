#!/bin/bash
#2019.07.26 by ZF

#Type 'sh script2_BUSCO_extraction.sh BUSCO_folder', e.g. sh script2_BUSCO_extraction.sh BUSCOs
#Tools TransDecoder, MAFFT, trimAl and FASconCAT-G are used in this script and will be automatically checked prior to formal analyses
#All the BUSCO results (run_* folders) are deposited in the same folder, e.g. BUSCOs/
#The final matrices and partition files are placed in 5-loci_concat/



##Checking the package dependency
echo "Checking the package dependency......"

#check TransDecoder
if [ $(which TransDecoder.LongOrfs) ]
    then
      echo "TransDecoder ...... OK"
      EXE_TRANSDECODER=$(which TransDecoder.LongOrfs)
      DIR_TRANSDECODER=${EXE_TRANSDECODER%/*}
    else
      until [ -x $DIR_TRANSDECODER/TransDecoder.LongOrfs ]
        do
          read -p "TransDecoder is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TransDecoder-v5.5.0):      " DIR_TRANSDECODER
        done
      echo "TransDecoder ...... OK"
fi

#check MAFFT
if [ $(which mafft) ]
    then
      echo "MAFFT ...... OK"
      EXE_MAFFT=$(which mafft)
      DIR_MAFFT=${EXE_MAFFT%/*}
    else
      until [ -x $DIR_MAFFT/mafft ]
        do
          read -p "MAFFT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_MAFFT
        done
      echo "MAFFT ...... OK"
fi

#check trimAL
if [ $(which trimal) ]
    then
      echo "trimal ...... OK"
      EXE_TRIMAL=$(which trimal)
      DIR_TRIMAL=${EXE_TRIMAL%/*}
    else
      until [ -x $DIR_TRIMAL/trimal ]
        do
          read -p "trimal is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/trimal-1.4.1/source):      " DIR_TRIMAL
        done
      echo "trimal ...... OK"
fi

#check FASconCAT
until [ -s $DIR_FASconCAT/FASconCAT*.pl ]
    do
      read -p "DIR_FASconCAT is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/FASconCAT-G):      " DIR_FASconCAT
    done
echo "FASconCAT ...... OK"

#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done


#copy the busco results to 0-raw_busco/
cp -r $1 0-raw_busco

#Generate a file containing species list
ls 0-raw_busco/ | sed 's/run_//g' > species.list

#Generate a file of loci list
SPECIES_NAME=$(cat species.list)
for SPECIES in $SPECIES_NAME
do
cd 0-raw_busco/run_$SPECIES/single_copy_busco_sequences
ls | sed 's/.fna\|.faa//g' >> ../../../temp.list; cd ../../..
done
sort -n temp.list | uniq > loci.list
rm temp.list

#Modify the head name of the fasta files for each locus
for SPECIES in $SPECIES_NAME
	do
	  cd 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/
	  sed -i -n '1,2p' EOG*
          for fna in *.fna
            do
              $DIR_TRANSDECODER/TransDecoder.LongOrfs -t $fna -m 20
              cp *dir/longest_orfs.cds $fna
              cp *dir/longest_orfs.pep $(ls $fna | cut -d . -f1).faa
              rm -rf *dir* *cmds
            done
	  sed -i -n '1,2p' EOG*
	  sed -i "1c >$SPECIES" EOG*
	  sed -i "s/*//g" EOG*.faa
	  cd ../../..
	done

#Merge sequences of the same locus into the fasta files
mkdir -p 1-raw_loci/fna 1-raw_loci/faa
LOCI_NAME=$(cat loci.list)
for LOCI in $LOCI_NAME; do touch 1-raw_loci/fna/$LOCI.fna && touch 1-raw_loci/faa/$LOCI.faa; done
for LOCI in $LOCI_NAME
do
    for SPECIES in $SPECIES_NAME
    do
	  if [ -f 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.fna ]
          then
            cat 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.fna >> 1-raw_loci/fna/$LOCI.fna 
            cat 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.faa >> 1-raw_loci/faa/$LOCI.faa
          else
            echo "$LOCI in $SPECIES does not exist" | tee -a 1-raw_loci/log.txt
       fi
    done
done

#Filter loci having too few taxa (less than three)
mkdir -p 2-loci_filter/fna 2-loci_filter/faa
TOTAL_TAXA=$(cat species.list | wc -l)
for LOCI in $LOCI_NAME
do
    echo -e $LOCI'\t'$(grep -o $LOCI 1-raw_loci/log.txt | wc -l) | tee -a 2-loci_filter/sequence_number.log
    if [ $(grep -o $LOCI 1-raw_loci/log.txt | wc -l) -lt `expr $TOTAL_TAXA - 2` ]
      then
        cp 1-raw_loci/fna/$LOCI.fna 2-loci_filter/fna
        cp 1-raw_loci/faa/$LOCI.faa 2-loci_filter/faa
        echo "$LOCI" >> 2-loci_filter/loci_name_filter.log
    fi
done

#Align all the faa and fna files
mkdir 3-faa_align
LOCI_NAME_FILTER=$(cat 2-loci_filter/loci_name_filter.log)
for LOCI_FILTER in $LOCI_NAME_FILTER
do
  $DIR_MAFFT/linsi --thread $THREADS 2-loci_filter/faa/$LOCI_FILTER.faa > 3-faa_align/$LOCI_FILTER.faa
  test -s 3-faa_align/$LOCI_FILTER.faa && echo "loci $LOCI_FILTER has been aligned" || $DIR_MAFFT/mafft --thread $THREADS 2-loci_filter/faa/$LOCI_FILTER.faa > 3-faa_align/$LOCI_FILTER.faa
done

#Trim the alignments
mkdir -p 4-loci_trim/fna 4-loci_trim/faa
for LOCI_FILTER in $LOCI_NAME_FILTER
do
  echo "trim protein sequence of loci $LOCI_FILTER ......"
  $DIR_TRIMAL/trimal -in 3-faa_align/$LOCI_FILTER.faa -out 4-loci_trim/faa/$LOCI_FILTER.aa.fas -automated1
  echo -e '\n'
  echo "trim nucleotide sequence of loci $LOCI_FILTER ......"
  $DIR_TRIMAL/trimal -in 3-faa_align/$LOCI_FILTER.faa -out 4-loci_trim/fna/$LOCI_FILTER.nuc.fas -automated1 -backtrans 2-loci_filter/fna/$LOCI_FILTER.fna
  echo -e '\n'
done

#Concatenate all the nucleotide/protein alignments in phylip format and generate partition files
mkdir -p 5-loci_concat/all/fna 5-loci_concat/all/faa
cd 4-loci_trim/fna
perl $DIR_FASconCAT/FASconCAT-G*pl -a -p -p -s -l
mv FcC* ../../5-loci_concat/all/fna/
cd ../faa
perl $DIR_FASconCAT/FASconCAT-G*pl -a -p -p -s -l
mv FcC* ../../5-loci_concat/all/faa/
cd ../../5-loci_concat/

#Generate the final matrices
for percent in 50 60 70 80 90 100
do
    mkdir -p matrix$percent/fna matrix$percent/faa
  for LOCI in $LOCI_NAME_FILTER
  do
    num1=$(echo "scale=1;$TOTAL_TAXA*(1-$percent/100)"|bc)
    num2=$(grep -o $LOCI ../1-raw_loci/log.txt | wc -l)
    c=$(echo "$num1 >= $num2" | bc)
    if [ $c -eq 1 ]
      then
        cp ../4-loci_trim/fna/$LOCI.nuc.fas matrix$percent/fna
        cp ../4-loci_trim/faa/$LOCI.aa.fas matrix$percent/faa
        echo $LOCI >> matrix$percent/busco$percent.loci.list
    fi 
  done
  cd matrix$percent/fna/
  perl $DIR_FASconCAT/FASconCAT-G*.pl -a -p -p -s -l
  cd ../faa
  perl $DIR_FASconCAT/FASconCAT-G*.pl -a -p -p -s -l
  cd ../../
done

#Summerize the matrices
echo -e '\n'
echo "protein BUSCO matrices" | tee -a summary.extraction
echo "primary matrix has $(awk 'NF{a=$0}END{print a}' all/faa/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat all/faa/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a summary.extraction
for percent in 50 60 70 80 90 100
do
  echo "$percent% matrix has $(awk 'NF{a=$0}END{print a}' matrix$percent/faa/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat matrix$percent/faa/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a summary.extraction
done
echo -e '\n'
echo "nucleotide BUSCO matrices" | tee -a summary.extraction
echo "primary matrix has $(awk 'NF{a=$0}END{print a}' all/fna/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat all/fna/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a summary.extraction
for percent in 50 60 70 80 90 100
do
  echo "$percent% matrix has $(awk 'NF{a=$0}END{print a}' matrix$percent/fna/FcC_supermatrix_partition.txt | sed -r 's/.*-(.*).*/\1/') sites and $(cat matrix$percent/fna/FcC_supermatrix_partition.txt | wc -l) loci" | tee -a summary.extraction
done
echo -e '\n'
echo "All the supermatrices and partition files are kept in 5-loci_concat/"
echo "Infmormation of matrices is saved in 5-loci_concat/summary.extraction"
