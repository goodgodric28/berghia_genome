cat protein_result_allRefSeq.txt | sed -e "s/^Record removed.*$//g" -e "s/^This sequence.*$//g" -e "/^$/d" | tr "\n" "\t" | sed -e "s/[0-9]*\. //g" -e "s/ \[/\t/g" -re "s/ (GI\:[0-9]*)\t/\t\1\n/g" -e "s/ aa protein//g" -e "s/\]//g" > protein_result_allRefSeq_edited.txt

awk -F '\t' '{ print $4, $1, $2, $3 }' OFS='\t' protein_result_allRefSeq_edited.txt | sort | uniq > protein_result_allRefSeq_edited_organized.txt

grep -v "uncharacterized protein"  protein_result_allRefSeq_edited_organized.txt > protein_result_allRefSeq_edited_organized_nouncharacterized.txt
