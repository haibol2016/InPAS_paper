ls results/010.bedtools.bedgraph.out/01.polyester.*/*unique.map.bedgraph | \
   perl -p -e 'BEGIN{print "tag\tcondition\tbedgraph_file\n"} s{(.+?/(sample_\d+).+)}{$2\t$2\t$1}' > docs_2/polyester.simulation.metadata.full.txt
