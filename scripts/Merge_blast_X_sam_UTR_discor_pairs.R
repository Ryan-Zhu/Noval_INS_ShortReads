local_lib <- file.path("/tscc/nfs/home/yzhu/miniconda3/lib/R/library")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse, lib.loc = local_lib))

option_list = list(
    make_option(c('-i', '--blast'),
                help = 'Input file from blast results.'),
    make_option(c('-s', '--sam'),
                help = 'Input headless SAM file.'),
    make_option(c('-a', '--anno'),
                help = 'File with read-sample annotations'),
    make_option(c('-o', '--out'),
                help = 'Output matched insertion sites.'),
    make_option(c('-c', '--qcov'), default = 0.8,
                help = 'Cutoff for query coverage.'),
    make_option(c('-p', '--pem'), default = 90.0,
                help = 'Cutoff for percent exact match.'),
    make_option(c('-b', '--min_bp'), default = 20,
                help = 'min alignment length.')

)
opt = parse_args(OptionParser(option_list=option_list))

blast <- fread(opt$blast, header=T)
if(nrow(blast) == 0){
    print("No blast hits. Exiting.")
    quit()
}

sam <- fread(opt$sam, sep="\t", header=F)
anno <- fread(opt$anno, header=F)

blast <- unique(blast)
sam <- unique(sam)
anno <- unique(anno)

blast <- blast[prct_exact_match >= opt$pem & alignment_length >= opt$min_bp,]

blast[blast_strand == "plus", blast_strand := "+"]
blast[blast_strand == "minus", blast_strand := "-"]
#blast[, blast_strand := mapvalues(x=blast[, blast_strand], from=c("plus", "minus"), to=c("+", "-"))]

# convert NCBI chromosome names to that of UCSC 
chr <- gsub(".+NC_0+([1-2]{1}[0-4]{1})\\..*", "\\1", blast$ref_id)
chr <- gsub("NC_012920\\..*", "chrM", chr)
chr <- paste0("chr", chr)
blast$ref_id <- chr
blast[ref_id == "chr23", ref_id := "chrX"]
blast[ref_id == "chr24", ref_id := "chrY"]

names(sam) <- c("read_id","FLAG_BAM", "ref_id", "mapping_POS", "MAPQ", "CIGAR", "ref_id_next", "POS_next", "genomic_dist", "BAM_query_seq", "BAM_strand")
names(anno) <- c("read_id", "sample_name", "gene_name", "gene_POS", "UTR", "side_clip", "utr_read_strand", "UTR_map_pos", "UTR_read_seq", "blast_query_seq")

merge <- merge(blast, sam, by=c("read_id","ref_id"), all.y=T, allow.cartesian=T) # right join, dropping alternate chromosomes from the blast result
merge <- merge[!is.na(query_length),]
#merge$clip_mapping_POS <- paste0(merge$ref_id,":", merge$ref_start, "-", merge$ref_end)

merge_anno <- merge(merge, anno, by="read_id", all.x=T, allow.cartesian=T) # left join
#merge_anno <- merge_anno[!is.na(UTR),]
merge_anno <- unique(merge_anno)

if(nrow(merge_anno) > 0){

    # reformat the blast ref_start and ref_end to follow the + strand coordinates
    if(nrow(merge_anno[blast_strand == "-"]) > 0){
        orig.start <- merge_anno[blast_strand == "-", ref_start]
        orig.end <- merge_anno[blast_strand == "-", ref_end]
        merge_anno[blast_strand == "-", ref_start := orig.end]  
        merge_anno[blast_strand == "-", ref_end := orig.start]

        orig.seq <- merge_anno[blast_strand == "-", blast_query_seq]
        rev.seq <- reverseComplement(DNAStringSet(orig.seq))
        merge_anno[blast_strand == "-", blast_query_seq := as.character(rev.seq)]
        }
        
        # add source gene strand to output
        merge_anno$gene_strand <- "-"
        merge_anno[(UTR=="5utr" & side_clip=="left") | (UTR=="3utr" & side_clip=="right"), gene_strand:="+"]

        # add side_clip adjusted by source and insertion side strands
        merge_anno[, side_clip_adj := ifelse(side_clip == "left", -1, 1)]
        merge_anno[utr_read_strand != blast_strand, side_clip_adj := side_clip_adj * -1]
        merge_anno[, side_clip_adj := ifelse(side_clip_adj == -1, "left", "right")]

        merge_anno$query_cov <- round(merge_anno[,alignment_length/query_length],digits=2)

        fwrite(merge_anno[query_cov >= opt$qcov,], opt$out, quote=F, sep="\t", col.names = T, na="NA")
        fwrite(merge_anno[query_cov < opt$qcov,], paste0(opt$out,".polyA"), quote=F, sep="\t", col.names = T, na="NA")

}else{
    print(paste0(basename(opt$blast), " and ", basename(opt$sam), " have no overlap at qcov = ", opt$qcov, " and pem = ", opt$pem))
}
