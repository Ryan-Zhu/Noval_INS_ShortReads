wd=/tscc/lustre/ddn/scratch/yzhu/gencDNA_0924/analysis_1111

cd $wd

mkdir -p split_bed_genes
mkdir -p collapsed_res


awk '{print > "split_bed_genes/"$5".txt"}' test_genes_ins.bed

parallel -j4 "sort -k 1,1 -k2,2n {} | bedtools merge -i - -d 10 -c 4,5 -o count_distinct,distinct > collapsed_res/{/}"  ::: split_bed_genes/*.txt

