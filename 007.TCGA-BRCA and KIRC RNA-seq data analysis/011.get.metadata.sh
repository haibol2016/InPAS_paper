
 ls results/010.bedtools.bedgraph.out/BRCA/*bedgraph | \
      perl -p -e 's{(.+/((.+?)_([TN]))_.+)}{$3\t$2\t$4\t$1}' | \
	  awk 'BEGIN{FS=OFS="\t"} NR % 2 == 1 {print "tag\tcondition\tbedgraph_file" >> "docs/BRCA/"$1".metadata.txt"} {print $2, $3, $4 >> "docs/BRCA/"$1".metadata.txt"}'


ls results/010.bedtools.bedgraph.out/KIRC/*bedgraph | \
        perl -p -e 's{(.+/((.+?)_([TN]))_.+)}{$3\t$2\t$4\t$1}' | \
		awk 'BEGIN{FS=OFS="\t"} NR % 2 == 1 {print "tag\tcondition\tbedgraph_file" >> "docs/metadata/KIRC/"$1".metadata.txt"} {print $2, $3, $4 >> "docs/metadata/KIRC/"$1".metadata.txt"}'

