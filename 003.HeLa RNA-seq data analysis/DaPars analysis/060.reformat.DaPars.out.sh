awk 'BEGIN{FS=OFS="\t"} $22 == "Y" '  results/015.Dapars.out/4.HeLa.PE.RNA-seq/01.HeLa.Github.DaPars_All_Prediction_Results.txt  | \
   perl -p -e 's{\|}{\t}g' | \
   awk 'BEGIN{FS=OFS="\t"} NR == FNR {a[$2] = $3"\t"$4; next} $1 in a {print $0, a[$1]}' docs_2/005.gencode.Tx.gene.name.type.txt  - |\
   cat <(head -1 results/015.Dapars.out/4.HeLa.PE.RNA-seq/01.HeLa.Github.DaPars_All_Prediction_Results.txt | \
   perl -pe 's{^Gene}{Transcript\tGene_id\tChr\tStrand}; s{$}{\tGene_symbol\tTranscript_type}') -



awk 'BEGIN{FS=OFS="\t"}NR >1 { start = gensub(/chr[^:]+:(.+?)-.+/,  "\\1", "g", $9); end =  gensub(/chr[^:]+:.+?-(.+)/,  "\\1", "g", $9); if ($6 == "+") {print $5, start -1, end, $9"_long", 1000, $6; print $5, start -1, $8, $9"_short", 200, $6} else {print $5, start -1, end, $9"_long", 1000, $6; print $5, $8 -1, end, $9"_short", 200, $6}}' 00.GitHub.DaPars.Final.sig.APA.txt  |sort -k1,1V -k2,2 |cat <(echo "track useScore=1") - > 00.DaPars.sig.APA.bed

awk 'BEGIN{FS=OFS="\t"} $42 == "TRUE" { start = $2; end = $3 ; if ($4 == "+") {print $1, start -1, end, $1":"$2"-"$3"_long", 1000, $4; print $1, start -1, $18, $1":"$2"-"$3"_short", 200, $4} else {print $1, start -1, end, $1":"$2"-"$3"_long", 1000, $4; print $1, $18 -1, end, $1":"$2"-"$3"_short", 200, $4}}' 010.InPAS.MAQC.Fisher.sig.APA.txt  |sort -k1,1V -k2,2 |cat <(echo "track useScore=1") - > 00.InPAS.fisher.sig.APA.bed


