# use ChiaSig to identify significant TLD-TLD interaction 
# https://github.com/cheehongsg/ChiaSigScaled

TAD_files=`ls "/home/wanghy/deDoc2_clique/N_results/ted_interaction1"

for filename in $TAD_files
do

./ChiaSig -n 300 -t 8 "/home/wanghy/deDoc2_clique/N_results/ted_interaction1/$filename" > "/home/wanghy/deDoc2_clique/N_results/ted_sig_interaction/${filename%%T*d}sig.TAD" 

done

