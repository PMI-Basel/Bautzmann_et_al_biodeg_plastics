#!/bin/bash
#SBATCH --qos=1day
#SBATCH --time=24:00:00
#SBATCH --mem=80g
#SBATCH --output=run.out
#SBATCH --error=run.error
#SBATCH --job-name=demultiplexing
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load modules
module load foss/2018b #interpreters
module load FastQC/0.11.8-Java-1.8
module load cutadapt/2.10-foss-2018b-Python-3.6.6
module load xlsx2csv/0.7.4-foss-2018b-Python-3.6.6

## --------------------------------------------------------------------
## Jan Wälchli | 15.02.2022 | Version 1.0
## --------------------------------------------------------------------

#running time notification
echo 'Start script'

#convert the design file from xlsx to tab
xlsx2csv ../../1_start/design.xlsx design.csv
cat design.csv | tr ',' '\t' | tail -n +2 > design.tab
rm design.csv

#get the runs
 runs=$(awk '{print $3}' design.tab | sort | uniq)

## --------------------------------------------------------------------
## A | Quality Control - FastQC
## --------------------------------------------------------------------

#create output folder
mkdir ../../4_output 2> /dev/null #suppress error message
rm -r  ../../4_output/qc 2> /dev/null
mkdir ../../4_output/qc

#quality control
fastqc -t 20 -k 0 -q ../../2_data/* -o ../../4_output/qc 2> /dev/null

#remove no longer needed files
rm ../../4_output/qc/*.zip

#running time notification
echo 'A - Quality Control done'

# --------------------------------------------------------------------
# B | Primer Files
# --------------------------------------------------------------------

#create folder
rm -r primer_cutted 2> /dev/null
mkdir primer_cutted
mkdir primer_cutted/primers

#create files
for run in ${runs}; do
	 grep ${run} design.tab | awk '{if ($2 == "b") print $0;}' | \
	 awk '{print $6, $7, $8, $9}' | sort | uniq \
	 > primer_cutted/primers/${run}_primers.txt
done

#running time notification
echo 'E - Primer Files done'

# ---------------------------------------------------------------------
# C | Primer cutting
# ---------------------------------------------------------------------

#repeat for each species in each run
for run in ${runs}; do

		while read p; do #loop over each primer combination
			#get sequence
			f_seq=$(echo ${p} | cut -f2 -d " ")
			r_seq=$(echo ${p} | cut -f4 -d " ")
			#input
			path_in=../../2_data/${run}*
			ls ${path_in} | cut -f4 -d "/" > filenames.txt
			#output
			path_out=primer_cutted/${run}/
			mkdir ${path_out}
			#cut primers for each file
			while read infile; do
				#part of the name to keep
				outname=$(echo ${infile} | cut -f1 -d ".")
				direction=$(echo ${infile} | cut -f2 -d "_" | cut -f1 -d ".")
				if [[ ${direction} = "r1" ]]; then
					cutadapt -g ${f_seq} -o ${path_out}${outname}"_cutted".fastq.gz ../../2_data/${infile}
				else cutadapt -g ${r_seq} -o ${path_out}${outname}"_cutted".fastq.gz ../../2_data/${infile}
				fi
			done < filenames.txt
			rm filenames.txt
		done < primer_cutted/primers/${run}_primers.txt

  done

  #running time notification
  echo 'F - Primer cutting done'

# --------------------------------------------------------------------
# D | Barcode Files
# --------------------------------------------------------------------

#create folder
rm -r  demultiplexed 2> /dev/null
mkdir demultiplexed
mkdir demultiplexed/barcodes

#create files
for run in ${runs}; do

	#barcodes
	grep ${run} design.tab | awk '{print $4, $5}' | \
	sort | uniq | tr ' ' '\n' | \
	sed 's'/'^F'/'>'${run}'-F'/'g' > demultiplexed/barcodes/${run}_barcodes.fasta

done

#running time notification
echo 'B - Barcode Files done'

## ---------------------------------------------------------------------
## E | Demultiplexing
## ---------------------------------------------------------------------

for run in ${runs}; do

	mkdir demultiplexed/${run}

  #read name and barcode
  while read name;
   do read bc;

   #remove '>'' from name
   outname=$(echo ${name} | tr -d '>')

   ###forward reads
   #open compressed file
   #grep reads with the barcode in the header
   #save in a new file
   gunzip -c primer_cutted/${run}/${run}_r1_cutted.fastq.gz | \
   grep -A3 -e ^@\.\*${bc} --no-group-separator\
   > demultiplexed/${run}/${outname}_F.fastq

   ###reverse reads
   gunzip -c primer_cutted/${run}/${run}_r2_cutted.fastq.gz | \
   grep -A3 -e ^@\.\*${bc} --no-group-separator\
   > demultiplexed/${run}/${outname}_R.fastq

 done < demultiplexed/barcodes/${run}_barcodes.fasta

done

#running time notification
echo 'C - Demultiplexing done'

# ---------------------------------------------------------------------
# F | Clean up
# ---------------------------------------------------------------------

#sort files by taxa

for run in ${runs}; do

	taxa=$(grep ${run} design.tab | awk '{print $2}' | sort | uniq)

	#bacteria
	 if [[ ${taxa} == *'b'* ]]; then
		mkdir demultiplexed/bacteria 2> /dev/null #may already exist
		mkdir demultiplexed/bacteria/${run};
		bac_primers=$(grep ${run} design.tab | awk '{if ($2 == "b") print $0;}' | awk '{print $4}' | sort | uniq)
		for b in ${bac_primers}; do
			mv demultiplexed/${run}/*${b}*.fastq.gz demultiplexed/bacteria/${run}
		done
	fi

	#fungi
	 if [[ ${taxa} == *'f'* ]]; then
		mkdir demultiplexed/fungi 2> /dev/null #may already exist
		mkdir demultiplexed/fungi/${run};
		fun_primers=$(grep ${run} design.tab | awk '{if ($2 == "f") print $0;}' | awk '{print $4}' | sort | uniq)
		for f in ${fun_primers}; do
			mv demultiplexed/${run}/*${f}*.fastq.gz demultiplexed/fungi/${run}
		done
	fi

done

#running time notification
echo 'D - Clean up done'

# ---------------------------------------------------------------------
# G | sequences tracking
# ---------------------------------------------------------------------

rm ../../2_data/seqs.txt 2> /dev/null
touch ../../2_data/seqs.txt

#number of raw sequences per run
for run in ${runs}; do

	#raw
	raw_lines=$(gunzip -c ../../2_data/${run}_r1.fastq.gz | wc -l)
	raw_seqs=$(echo ${raw_lines} / 4 | bc) #each seq has 4 lines

	for species in bacteria fungi; do

		#input
		path_in=demultiplexed/${species}/${run}/

		#count number of demultiplexed sequences
		if [ -d "${path_in}" ]; then #check if path exist, otherwise set to NA
			dem_lines=$(gunzip -c ${path_in}*r1* | wc -l)
			dem_seqs=$(echo ${dem_lines} / 4 | bc) #each seq has 4 lines
		else dem_seqs="NA"
		fi

		#save to corresponding taxa
		if [ ${species} == "bacteria" ]; then
			dem_seqs_bac=$(echo ${dem_seqs})
		else dem_seqs_fun=$(echo ${dem_seqs})
		fi

	done

	#output
	echo ${run} "raw" ${raw_seqs} "bac" ${dem_seqs_bac} "fun" ${dem_seqs_fun} >> ../../2_data/seqs.txt

done

echo 'F - sequences tracking done'
echo 'End script'
