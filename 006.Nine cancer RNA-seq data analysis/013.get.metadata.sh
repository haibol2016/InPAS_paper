

ls ~/marcus_ruscetti/haibo/results/010.bedtools.bedgraph.out/* | \
      perl -p -e 's{(.+/((.+?)_([TN]))-.+)}{$3\t$2\t$4\t$1}' |\
      awk 'BEGIN{FS=OFS="\t"} NR % 2 == 1 {print "tag\tcondition\tbedgraph_file" >> "docs_2/007.nine.cancer/"$1".metadata.txt"} {print $2, $3, $4 >> "docs_2/007.nine.cancer/"$1".metadata.txt"}'
