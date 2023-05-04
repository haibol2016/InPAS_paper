

ls * | grep -F -w -f - ~/arthur_mercurio/Haibo/InPAS/TCGA_KIRC_BRCA/TCGA_BRCA.manifest.metadata.txt | \
             awk 'BEGIN{FS=OFS="\t"} {type = $13~/Normal/? "N" : "T";  print "mv", $2, $11 "_"  type "_"$1".bam"}' |sh

ls * | grep -F -w -f - ~/arthur_mercurio/Haibo/InPAS/TCGA_KIRC_BRCA/TCGA_KIRC.manifest.metadata.txt | \
             awk 'BEGIN{FS=OFS="\t"} {type = $13~/Normal/? "N" : "T";  print "mv", $2, $11 "_"  type "_"$1".bam"}' |sh