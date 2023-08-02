#pull out the transcripts from the filtered gtf file
awk '{print $10}' augustus.hints.tsebra.gtf > tsebrasupport.genes

#filter so that you get a file containing unique names
cat tsebrasupport.genes | uniq > tsebrasupport.genes.uniq

#remove the " and ; from each name
sed -i 's/"//g' tsebrasupport.genes.uniq
sed -i 's/;//g' tsebrasupport.genes.uniq
sed -i 's/anno1.//g' tsebrasupport.genes.uniq

#perl one liner that uses the id file you just created to pull out corresponding sequences in the fasta file produced by braker
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' tsebrasupport.genes.uniq augustus.hints.aa > augustus.hints.tsebra.aa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' tsebrasupport.genes.uniq augustus.hints.codingseq > augustus.hints.tsebra.codingseq
