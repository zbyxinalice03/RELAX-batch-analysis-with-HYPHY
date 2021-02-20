#!/bin/bash
# written by Zhuby
# This analysis did not detect any evidence of relaxed selection. 
# However, if it had, a significant K>1 would indicate intensified selection on test lineages
# and significant K<1 would indicate relaxed selection on test lineages
# conda activate hyphy # conda create -n hyphy, then, conda install hyphy 

# first step: rm Stop codon for single sequence
#for i in `ls *.fas`; do 
#	name=${i%.fas}; 
#	`less ${name}.fas | awk '{ if (NR%2==0) {print substr($0,1,length($0)-3)} else{print$0} }' > ${name}.fasta` 
#done

# next step: to cal Selection pressure, prepare fasta_file(no matter num) & this sh_file & tree_file, if have errors, rm all dir 
for fasta_file in `ls *.fasta`; do
	name_variable=$(basename $fasta_file | sed 's/\.fasta//');
	mkdir -p ${name_variable};
	wait;
	path=$(readlink -fz ${name_variable}); # to get full {path} to new dir
	cp $fasta_file ${path}/$fasta_file; # one phy_file for single dir 
	cp test.tree ${path}/test.tree;
	cat ${path}/$fasta_file ${path}/test.tree > ${path}/${name_variable}.fna;
	nohup hyphy relax  CPU=20 --alignment ${path}/${name_variable}.fna\
	--test test --output ${path}/${name_variable}.relax.json --grid-size 250 --starting-points 1 \
	--rates 3 --models All --srv No  --reference-group Unlabeled branches --syn-rates 3 --save-fit /dev/null > ${path}/${name_variable}.output 2>&1 &;
	# extract the K value & P value
	grep -h '>Loaded a multiple' ${path}/${name_variable}.output | sed 's/>.*\///g' | sed 's/\.fna`//g' > ${path}/genename;
	grep -h 'Relaxation/intensification parameter' ${path}/${name_variable}.output | awk '{print$6}'  > ${path}/Kvalue;
	grep -h 'Likelihood ratio' ${path}/${name_variable}.output | awk '{print$6}'| sed 's/\*\*\.//g' > ${path}/Pvalue;
	paste ${path}/genename ${path}/Kvalue ${path}/Pvalue > ${path}/${name_variable}.K_P.list;
	cp ${path}/${name_variable}.K_P.list ../${name_variable}.K_P.list;
	rm ${path}/genename ${path}/Kvalue ${path}/Pvalue;
	cat *.K_P.list > relax.genes.txt
done
	
ls ./*/*.output | wc -l 
	
	
