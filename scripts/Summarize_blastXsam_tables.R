# combine blast_X_sam outputs for all samples

local_lib <- file.path("/tscc/nfs/home/yzhu/miniconda3/lib/R/library")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse, lib.loc = local_lib))
suppressPackageStartupMessages(library(dplyr, lib.loc = local_lib))
suppressPackageStartupMessages(library(GenomicAlignments, lib.loc = local_lib))
suppressPackageStartupMessages(library(DescTools, lib.loc = local_lib))

option_list = list(
    make_option(c('-i', '--indir'),
                help = 'Input directory.'),
    make_option(c('-o', '--out'),
                help = 'Output file of final result.'),
    make_option(c('-d', '--distance'), default = 5000,
                help = 'Largest genomic distance cutoff.'),
    make_option(c('-p', '--pattern'), default = "highQual.raw",
                help = 'Postfix pattern to search for input files')
)

opt = parse_args(OptionParser(option_list=option_list))

#dir="/oasis/tscc/scratch/yzhu/3-11-inser_test_2_samples_3_genes"

files <- list.files(path=opt$indir, pattern=paste0(".*", opt$pattern, "$"), full.names=T)

combined <- rbindlist(lapply(files, fread, header=T))

combined1 <- combined[!is.na(BAM_strand), ] # pair with low quality mate reads
combined2 <- combined[is.na(BAM_strand), ] # pair with unmapped mate reads

combined1$read_len_ref <- cigarWidthAlongReferenceSpace(combined1$CIGAR)
setnames(combined1, "mapping_POS", "MATE_POS_start")
combined1[,("MATE_POS_end") := MATE_POS_start+read_len_ref-1]

# calculate overlap and span of  two reads
sub_comb <- combined1[, list(ref_start, ref_end, MATE_POS_start, MATE_POS_end)]
combined1$Overlap <- apply(sub_comb, 1, function(x) Overlap(x[1:2], x[3:4]))

# the Overlap function returns (overlap-1) so we need to adjust for it
combined1[Overlap != 0, Overlap := Overlap+1]
combined1[ref_end == MATE_POS_start | ref_start == MATE_POS_end, Overlap:=Overlap+1]

combined1$genomic_dist <- apply(sub_comb, 1, function(x) max(x)-min(x)+1)

combined1 <- combined1[genomic_dist <= opt$distance]

combined1$UTR_clip_mapping_POS <- paste0(combined1$ref_id,":", combined1$ref_start, "-", combined1$ref_end)
combined1$Mate_mapping_POS <- paste0(combined1$ref_id,":", combined1$MATE_POS_start, "-", combined1$MATE_POS_end)
combined1 <- combined1[order(read_id, ref_id)]
combined1 <- unique(combined1)


combined1 <- combined1[,c("read_id","ref_id","query_length","alignment_length","prct_exact_match","query_cov","query_start","query_end","ref_start","ref_end","Evalue","Bitscore","blast_strand","BAM_strand","utr_read_strand","CIGAR","FLAG_BAM","MAPQ","read_len_ref","ref_id_next","POS_next","MATE_POS_start","MATE_POS_end","UTR_clip_mapping_POS","Mate_mapping_POS","UTR_map_pos","gene_POS","gene_strand","sample_name","genomic_dist","Overlap","UTR_read_seq", "blast_query_seq","BAM_query_seq","side_clip", "side_clip_adj", "gene_name","UTR")]

fwrite(combined1, opt$out, quote=F, sep="\t", col.names = T, na="NA")

if(nrow(combined2)!=0){
	combined2$read_len_ref <- NA
	setnames(combined2, "mapping_POS", "MATE_POS_start")
	combined2$MATE_POS_end <- NA
	combined2$Overlap <- NA
	combined2$genomic_dist <- NA
	combined2$UTR_clip_mapping_POS <- paste0(combined2$ref_id,":", combined2$ref_start, "-", combined2$ref_end)
	combined2$Mate_mapping_POS <- NA
	combined2 <- combined2[order(read_id, ref_id)]
	combined2 <- unique(combined2)

	fwrite(combined2, paste0(opt$out, ".unmapped"), quote=F, sep="\t", col.names = T, na="NA")

}