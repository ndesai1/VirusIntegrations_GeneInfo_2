#first, need to get chromosome location from table with int locations
args <- commandArgs(TRUE)
vir = toString(args[1])

filedest = getwd()

filename=paste(c(filedest, "/", vir, "_virusint/", vir, "_samples_with_locus.csv"), collapse="")

chrm_data=read.csv(filename)
chrm_data["integration_location_type"] <-NA
chrm_data["gene_id"] <- NA
chrm_data["nearest_gene"] <-NA
chrm_data["dist_to_nearest_gene"] <- NA
chrm_data["int_point_upstream_or_downstream_of_nearest_gene"] <- NA
chrm_data["gene_symbol"] <- NA
chrm_data["gene_description"] <- NA
chrm_data["number_of_cancer_related_publications"] <- NA
chrm_data["cancer_related_pub_proportion"] <- NA
chrm_data["number_of_cancer_related_rifs"] <- NA
chrm_data["cancer_related_rif_proportion"] <- NA

#chrm_data["int_site_upstream_of"] <- NA
#chrm_data["distance_to_next_gene"] <- NA
#chrm_data["int_site_downstream_of"] <- NA
#chrm_data["distance_to_prev_gene"] <- NA

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(rentrez)
library(org.Hs.eg.db)


genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

for (entry in 1:NROW(chrm_data)) {
    tryCatch({
    if (chrm_data[entry,]$Chromosome.1 != 'chrVirus')
    {
       chr_num=toString(chrm_data[entry,]$Chromosome.1)
       start = as.numeric(as.character(chrm_data[entry,]$Position.1))
       strandch = toString(chrm_data[entry,]$Strand.1)
    }else if (chrm_data[entry,]$Chromosome.2 != 'chrVirus')
    {
	chr_num=toString(chrm_data[entry,]$Chromosome.2)
	start = as.numeric(as.character(chrm_data[entry,]$Position.2))
	strandch = toString(chrm_data[entry,]$Strand.2)
    }
    
    print (chr_num)
    print(start)
    print(strandch)
    
    chrmloc<-GRanges(chr_num, IRanges(start, start), strand = strandch)
    loc_all<- locateVariants(chrmloc, txdb, AllVariants())
    int_site_type<-paste(unique(loc_all$LOCATION), sep=",", collapse=",")
    

#    print(all_genes)
    nearest_gene_id=genes[nearest(chrmloc,genes, ignore.strand=FALSE)]$gene_id
#    print (nearest_gene_id)
    preceding_gene=genes[precede(chrmloc, genes, ignore.strand=FALSE)]$gene_id
    preceding_gene_info=genes[precede(chrmloc, genes, ignore.strand=FALSE)]
#    print ("preceding gene is:")
#    print (preceding_gene)
#    print (preceding_gene_info)
    following_gene=genes[follow(chrmloc, genes, ignore.strand=FALSE)]$gene_id
    following_gene_info=genes[follow(chrmloc,genes, ignore.strand=FALSE)]
#    print("following gene is:")
#    print(following_gene)
#    print (following_gene_info)
    nearest_gene_info=genes[nearest(chrmloc,genes)]
    distance=distance(chrmloc, nearest_gene_info)
    distance_to_preceding=distance(chrmloc, preceding_gene_info)
    distance_to_following=distance(chrmloc, following_gene_info)
    if (nearest_gene_id == preceding_gene)
    {
	pre_follow = "upstream"
    }else if (nearest_gene_id == following_gene)
    {
	pre_follow = "downstream"
    }else if (nearest_gene_id != preceding_gene && nearest_gene_id != following_gene)
    {
	pre_follow = "overlapping"		
    }
    
    chrm_data[entry,]$integration_location_type = int_site_type
    
    if (toString(int_site_type) != "intergenic")
    {
       print (toString(int_site_type))
       chrm_data[entry,]$gene_id <- paste(unique(loc_all$GENEID), sep=",", collapse=",")
    }
    
    chrm_data[entry,]$nearest_gene = nearest_gene_id
    chrm_data[entry,]$dist_to_nearest_gene = distance
    chrm_data[entry,]$int_point_upstream_or_downstream_of_nearest_gene=pre_follow

    }

    if (distance <= 50000)
    {
    entrez = toString(nearest_gene_id)
#    geneinfo = entrez_summary(db="gene", id = entrez)
    chrm_data[entry,]$gene_symbol = select(org.Hs.eg.db, keys = entrez, columns=c("SYMBOL"), keytype="ENTREZID")
#    chrm_data[entry,]$gene_description = geneinfo$summary

    #get information on number of publications with cancer related mentions
#    record <- entrez_link(dbfrom='gene', id = entrez, db='all')
#    toMatch <- c("cancer", "tumor", "tumour", "carcinogenic", "carcinoma", "oncogene", "oncogenic", "sarcoma", "tumorigenesis", "metastasis")

#    if (is.null(record$links$gene_pubmed) == FALSE){
#    publication_summary <- entrez_summary(db="pubmed", id = c(record$links$gene_pubmed))
    
#    publication_titles = extract_from_esummary(publication_summary, "title")
    
#    publication_matches <- unique(grep(paste(toMatch,collapse="|"), publication_titles, value=TRUE))

    #number of cancer-related publications and rifs, and proportion of publications and rifs related to cancer
#    chrm_data[entry,]$number_of_cancer_related_publications = length(publication_matches)
#    chrm_data[entry,]$cancer_related_pub_proportion = as.numeric(length(publication_matches)/length(publication_titles))
#    }
    
#    if (is.null(record$link$gene_pubmed_rif) == FALSE){
#    rif_summary <- entrez_summary(db="pubmed", id = c(record$links$gene_pubmed_rif))
#    rif_titles= extract_from_esummary(rif_summary, "title")
#    rif_matches = unique(grep(paste(toMatch,collapse="|"), rif_titles, value=TRUE))
#    chrm_data[entry,]$number_of_cancer_related_rifs = length(rif_matches)
#    chrm_data[entry,]$cancer_related_rif_proportion = as.numeric(length(rif_matches)/length(rif_titles))
    }   
    
    }
    })






     print (chrm_data[entry,])

}

write.csv(chrm_data, paste(c(filedest,"/", vir, "_virusint/", vir, "_samples_with_locus_and_genes_closestgenes.csv"), collapse=""))