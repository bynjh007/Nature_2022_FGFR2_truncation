#!/bin/sh

# -p : path where star-fusion output is located
# -s : file including list of samples to be analyised (*** this should be unix format :dos2unix ***)
# -t : target gene
# -g : gtf file to extract genomic locus of target gene

while getopts "p:s:t:g:" opt; do
	case $opt in
		p) path="$OPTARG"
		;;
		s) samp="$OPTARG"
		;;
		t) target="$OPTARG"
		;;
		g) gtf="$OPTARG"
		;;
		# if argument is not one of -f, -t, -g, break the code
		*)
			help
			exit 0
		;;
	esac
done

# target gene location
target_chr=$(grep \"${target}\" ${gtf} | awk '$3 == "gene" {print $1}');
target_start=$(grep \"${target}\" ${gtf} | awk '$3 == "gene" {print $4}');
target_end=$(grep \"${target}\" ${gtf} | awk '$3 == "gene" {print $5}');

# read files for each sample (read from the second line as first line is "SAMPLE_ID")
for i in $(tail -n +2 ${samp}); do
	
	path_file="${path}/${i}/star-fusion.preliminary";
	
	# all the reads with brkpt in target region from the gene-mapped junction file
	chim_junc="${path_file}/star-fusion.junction_breakpts_to_genes.txt";
	awk -v chr="${target_chr}" -v start="${target_start}" -v end="${target_end}" '($1==chr && $2>=start && $2<=end) || ($4==chr && $5>=start && $5<=end)' ${chim_junc} > "${chim_junc}.${target}";

	# failed reads that contain selfie, homology, or repeat match
	fail_file="${path_file}/star-fusion.junction_breakpts_to_genes.txt.fail";
	sed -n '/Contains/,/^$/p' ${fail_file} | awk -v chr="${target_chr}" -v start="${target_start}" -v end="${target_end}" '($1==chr && $2>=start && $2<=end) || ($4==chr && $5>=start && $5<=end) {print $10}' | grep ":" | sort | uniq > "${fail_file}.list";

	# excluding the failed reads from the all the reads with brkpt in target region
	grep -v -f ${fail_file}.list ${chim_junc}.${target} > "${chim_junc}.${target}.filtered";

	# all the qualified reads for the target to predict fusion in Star-Fusion
	pass_file="${path_file}/star-fusion.junction_breakpts_to_genes.txt.pass";
	grep \"$target\^ ${pass_file} > "${pass_file}.${target}";

done;




	

