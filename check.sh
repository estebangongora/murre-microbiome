#!/bin/bash
awk 'BEGIN{OFS="\n"}
     NR>1{for (i=2;i<=NF;i++) a[i]+=$i}
     END{for (i=2;i<=NF;i++) printf a[i] OFS}' check.txt > sums_2.txt
awk 'NR==2' check.txt | awk 'BEGIN{RS="\t"} {$1 = $1} NR > 1 { print "" } 1' | sed '/^\s*$/d' | tail -n+2 > neg1_2.txt
awk 'NR==3' check.txt | awk 'BEGIN{RS="\t"} {$1 = $1} NR > 1 { print "" } 1' | sed '/^\s*$/d' | tail -n+2 > neg2_2.txt
awk '{print $1}' check.tsv | tail -n+3 > header_2.txt
paste header_2.txt sums_2.txt neg1_2.txt neg2_2.txt detail.txt > all.txt
awk '{$6 = $2 - $3 - $4} {$7 = $6 / $5} {$8 = ($3 / $2) * 100} {$9 = ($3 / $7) * 100} {$10 = ($4 / $2) * 100} {$11 = ($4 / $7) * 100}1' all.txt > prop.txt

awk '{
if ($8 > 10)
	print $1"\t"$2"\t"$3"\tCheck_Neg1_T";
else if($9 > 10)
	print $1"\t"$2"\t"$7"\t"$3"\tCheck_Neg1_P";
else print $1,"\tOK_Neg1";
}' prop.txt > checked_1.txt
awk '{
if ($10 > 10)
	print $1"\t"$2"\t"$4"\tCheck_Neg2_T";
else if($11 > 10)
	print $1"\t"$2"\t"$7"\t"$4"\tCheck_Neg2_P";
else print $1,"\tOK_Neg2";
}' prop.txt > checked_2.txt
paste checked_1.txt checked_2.txt > checked_3.txt

grep 'Check' checked_3.txt > checked.txt
