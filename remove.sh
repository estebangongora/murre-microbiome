#!/bin/bash
#For samples with Neg 2
awk 'BEGIN{OFS="\n"}
     NR>1{for (i=2;i<=NF;i++) a[i]+=$i}
     END{for (i=2;i<=NF;i++) printf a[i] OFS}' neg2_and_contams.txt > sums.txt
awk 'NR==2' neg2_and_contams.txt | awk 'BEGIN{RS="\t"} {$1 = $1} NR > 1 { print "" } 1' | sed '/^\s*$/d' | tail -n+2 > neg2.txt
awk '{print $1}' neg2_and_contams.tsv | tail -n+3 > header.txt
paste header.txt sums.txt neg2.txt > both.txt

awk '{
if ($3 < $2)
	print $1"\tCheck";
else if($3 = $2)
	print $1"\tDelete";
else print $1,"\tOK";
}' both.txt > comp_1.txt
awk '{
if ($3==0)
	print "OK";
else print "Check";
}' both.txt > comp_2.txt
paste comp_1.txt comp_2.txt > comp_3.txt
awk '{
if ($3=="OK")
	print $1"\t"$3;
else print $1"\t"$2;
}' comp_3.txt > comp.txt

grep 'Delete' comp.txt | grep -o 'ASV_[0-9]*' > remove.txt
sed -i '1s/^/FeatureID\n/' remove.txt
grep 'Check' comp.txt | grep -o 'ASV_[0-9]*' > check.txt
sed -i '1s/^/FeatureID\n/' check.txt

#####################################################################################

#For samples with Neg 1
awk 'BEGIN{OFS="\n"}
     NR>1{for (i=2;i<=NF;i++) a[i]+=$i}
     END{for (i=2;i<=NF;i++) printf a[i] OFS}' neg1_and_other-samples.txt > sums_N1.txt
awk 'NR==2' neg1_and_other-samples.txt | awk 'BEGIN{RS="\t"} {$1 = $1} NR > 1 { print "" } 1' | sed '/^\s*$/d' | tail -n+2 > neg1.txt
awk '{print $1}' neg1_and_other-samples.tsv | tail -n+3 > header_N1.txt
paste header_N1.txt sums_N1.txt neg1.txt > both_N1.txt

awk '{
if ($3 < $2)
	print $1"\tCheck";
else if($3 = $2)
	print $1"\tDelete";
else print $1,"\tOK";
}' both_N1.txt > comp_N1_1.txt
awk '{
if ($3==0)
	print "OK";
else print "Check";
}' both_N1.txt > comp_N1_2.txt
paste comp_N1_1.txt comp_N1_2.txt > comp_N1_3.txt
awk '{
if ($3=="OK")
	print $1"\t"$3;
else print $1"\t"$2;
}' comp_N1_3.txt > comp_N1.txt

grep 'Delete' comp_N1.txt | grep -o 'ASV_[0-9]*' >> remove.txt
awk '!seen[$0]++' remove.txt > asvs_to_remove.txt
grep 'Check' comp_N1.txt | grep -o 'ASV_[0-9]*' >> check.txt
awk '!seen[$0]++' check.txt > asvs_to_check.txt
