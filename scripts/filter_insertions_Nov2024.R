library(data.table)
library(stringr)
library(valr)
library(ggplot2)


setwd("/tscc/lustre/ddn/scratch/yzhu/gencDNA_0924/analysis_1111/")

check_intersect <- function(x, y){
  pattern = "^(chr.{1,2}):([0-9]+)-([0-9]+)$"
  chr.x = gsub(pattern, "\\1", x) 
  start.x = as.numeric(gsub(pattern, "\\2", x))
  end.x = as.numeric(gsub(pattern, "\\3", x))
  chr.y = gsub(pattern, "\\1", y) 
  start.y = as.numeric(gsub(pattern, "\\2", y))
  end.y = as.numeric(gsub(pattern, "\\3", y))
  out = F
  
  if (chr.x == chr.y){
    X = tibble::tribble(
      ~chrom, ~start, ~end,
      chr.x, start.x, end.x)
    
    Y = tibble::tribble(
      ~chrom, ~start, ~end,
      chr.y, start.y, end.y)
    
    
    n = nrow(bed_intersect(X, Y))
    
    if (n!=0){out=T}else{out=F}
    
  }
  
  return(out)
  
}

raw <- fread("/tscc/lustre/ddn/scratch/yzhu/gencDNA_0924/analysis_1111/combined_res/combined_res_Nov2024.txt")
raw[, Cate:="test_genes"]


# since overlap indicates TSD, we need to make this more stringent
# if overlap>0, then both blast and mate alignments must be precise (prct_exact_match=100,	query_cov=1, CIGAR="^[0-9]+M$")

# discard if overlap >= 40 because TSD is usually 20bp

filtered <- raw[Overlap>0 & Overlap<40 & prct_exact_match==100 & query_cov==1][grep("^[0-9]+M$", CIGAR)]

filtered <- rbind(filtered, raw[Overlap==0 & prct_exact_match>= 95 & query_cov >= 0.95])

# inser size should be 300~500 so discard anything with genomic_dist> 2000
filtered <- filtered[genomic_dist < 2000]

# clipped seq cannot be within mate seq

filtered[, overlap_start_on_mate:="NA"][, overlap_end_on_mate:="NA"]
for(i in 1:nrow(filtered)){
  if (filtered[i, query_length==Overlap]) {
    match_loc = str_locate(filtered[i, BAM_query_seq], pattern=filtered[i, blast_query_seq])
    filtered[i, overlap_start_on_mate:=match_loc[1,1]][i, overlap_end_on_mate:=match_loc[1,2]]
  }
}

filtered[, bam_len := nchar(BAM_query_seq)]
filtered <- rbind(filtered[query_length==Overlap & (overlap_start_on_mate == 1 | overlap_end_on_mate == bam_len)], filtered[overlap_start_on_mate=="NA" ])


# annotate side info of the insertion site
filtered[, Note:="2_sides"]
filtered[side_clip_adj=="right" & MATE_POS_end > ref_end, Note:="1_side"]
filtered[side_clip_adj=="left" & MATE_POS_end < ref_end, Note:="1_side"]

# annotate potential inversion; we make sure the clipped alignment has no mismatches to increase confidence

filtered[, keep:=1]
filtered[blast_strand==BAM_strand, Note:="potential inversion"] 

filtered[blast_strand==BAM_strand & (prct_exact_match!=100 | query_cov!=1), keep:=0]

filtered <- filtered[keep==1]
filtered[,keep:=NULL]

# exclude insertion sites within the gene

if_intersect <- c()

for (i in 1:nrow(filtered)){
  if_intersect[i] = (check_intersect(filtered$UTR_clip_mapping_POS[i],
                                     filtered$gene_POS[i]))
}


filtered <- filtered[!if_intersect]

anno <- fread("/tscc/nfs/home/yzhu/projects/gencDNA/2024june_TSCC_bwaPaths_V2.txt", head=F)
anno[, sample_name:=basename(V1)]
anno[, study:=basename(dirname(dirname(V1)))]
#anno[, sample_frac_pct:=round(.N/nrow(anno),3)*100, by=.(study)]
anno[, num_samples:=.N, by=.(study)]
setnames(anno, "V1", "path")

filtered <- merge(filtered, anno, by="sample_name")

#fwrite(filtered, "filtered_combined_ins_Nov2024.txt", sep = "\t", quote=F, na="NA")

#fwrite(filtered, "filtered_combined_ins_10152024.txt", sep = "\t", quote=F, na="NA")

# read distribution by study
filtered_uniq_reads <- unique(filtered[,.(sample_name, read_id, study, num_samples)])
study_dist <- filtered_uniq_reads[,.N, by=.(study,num_samples)]

setnames(study_dist, "N","num_ins_support_reads")
study_dist[, avg_ins_read_per_sam:=round(num_ins_support_reads/num_samples,2)]

#fwrite(study_dist, "Ins_support_reads_by_study.csv")

# merge neighbouring insertion sites (distance<=10bp) to aggreate read count at insertion site level by gene

bed <- filtered[Cate=="test_genes", list(ref_id, ref_start, ref_end, read_id, gene_name)][, ref_start:= ref_start-1]
fwrite(bed, "test_genes_ins.bed", sep = "\t", col.names = F)

# next run the bash commands in Step2_bash_component_bedtools.analysis.sh 

# after obtaining collapsed_res we continue to visualization
res <- rbindlist( lapply(list.files("collapsed_res/", full.names = T), fread)) 

# insertion-site level result

names(res) <- c("chrom", "start", "end", "num_uniq_reads", "gene_name")
res[, Ins_name := paste0(gene_name,"_",chrom,":",start,"-",end)]

res <- res[order(-num_uniq_reads)]

dir.create("plots")

ggplot(res[1:20], aes(x=reorder(Ins_name, num_uniq_reads), y=num_uniq_reads)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  theme_light() +
  coord_flip() +
  geom_text(aes(label=num_uniq_reads), hjust=-0.3, size=3) +
  labs(x="Insertion site", y="Num. supporting reads", title="Top 20 insertion sites by #supporting reads")
ggsave("plots/ins_by_read.jpg", dpi=300, width = 10, height = 5, units = "in")

fwrite(res, "test_genes_INS_Level_unique_reads_Nov2024.csv")


# gene-level result: num. unique reads by gene
tmp <- res[, .(num_reads=sum(num_uniq_reads)), by=gene_name][order(-num_reads)]


ggplot(tmp[num_reads>=8], aes(x=reorder(gene_name, -num_reads), y=num_reads)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  theme_light() +
  geom_text(aes(label=num_reads), vjust=-0.3, size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7)) +
  labs(x="Gene", y="Num. supporting reads", title="Genes with insertion-site-supporting-reads>=8")

ggsave("plots/gene_by_read.jpg", dpi=300, width = 12, height = 5, units = "in")



# distribution of # reads per ins

dist <- as.data.table(table(res$num_uniq_reads))
names(dist) <- c("Num_uniq_reads", "Num_uniq_ins")
dist[, freq:=round(Num_uniq_ins/sum(Num_uniq_ins),3)]

# > dist
# Num_uniq_reads Num_uniq_ins  freq
# 1:              1         1417 0.860
# 2:              2          168 0.102
# 3:              3           35 0.021
# 4:              4           11 0.007
# 5:              5            8 0.005
# 6:              6            1 0.001
# 7:              7            1 0.001
# 8:              8            2 0.001
# 9:             10            1 0.001
# 10:             12            1 0.001
# 11:             38            1 0.001
# 12:            501            1 0.001



# gene-level data: num. unique insertion sites by gene

# an insertion site must have >2 supporting reads
res.gene <- res[num_uniq_reads>2, .(num_uniq_ins_sites=.N), by=gene_name][order(-num_uniq_ins_sites)]

# fwrite(res.gene, "test_genes_num_ins_above2reads.txt", sep = "\t")
# fwrite(res.gene, "test_genes_num_ins_above2reads_PP.txt", sep = "\t")
