sed 's/\t\t/ /g' Orthogroups.tsv | sed 's/\t/ /g' | sed 's/,//g' | sed 's/^M$//g' | sed -r 's/(OG[0-9][0-9][0-9][0-9][0-9][0-9][0-9])\ +/\1\ /g' > Orthogroups_edited.txt

dos2unix Orthogroups_edited.txt

tail -n +2 > Orthogroups_edited.txt.temp
mv Orthogroups_edited.txt.temp Orthogroups_edited.txt
