library(rentrez)

entrez = c("1510","7031","4588","51208","56667","340485","23584","4586","7032","10551","135656","7033","11199", "340547", "389015", "200504", "56287")
entrez_all = c("1510","7031","8970","5644","760","5284","90993","5950","3006","1755","150","25878","1357","7503","10418","84699","9271","29850","4588","8618","5176","5328","7169", "1075")
gene_pub_info <- NULL

for (entry in entrez_all){
    record <- entrez_link(dbfrom='gene', id = entry, db='all')
    gene_info <- entrez_summary(db='gene', id=entry)
    gene_name <- gene_info$name
    #print(gene_name)
    publication_summary <- entrez_summary(db="pubmed", id = c(record$links$gene_pubmed), always_return_list=TRUE, web_history=NULL)
    rif_summary <- entrez_summary(db="pubmed", id = c(record$links$gene_pubmed_rif), always_return_list=TRUE)
    #print(entry)
    publication_titles = extract_from_esummary(publication_summary, "title")
rif_titles= extract_from_esummary(rif_summary, "title")
toMatch <- c("cancer", "tumor", "tumour", "carcinogenic", "carcinoma", "oncogene", "oncogenic", "sarcoma", "tumorigenesis", "metastasis")

publication_matches <- unique(grep(paste(toMatch,collapse="|"), publication_titles, value=TRUE))
rif_matches = unique(grep(paste(toMatch,collapse="|"), rif_titles, value=TRUE))

pub_number<- length(publication_matches)
rif_number<- length(rif_matches)

#get proportion of publications that match to cancer terms
pub_match_prop = as.numeric(length(publication_matches)/length(publication_titles))
rif_match_prop = as.numeric(length(rif_matches)/length(rif_titles))

row <- cbind(gene_name, entry, pub_number, pub_match_prop, rif_number, rif_match_prop)
gene_pub_info <- rbind(gene_pub_info, row)

print(row)

}

write.csv(gene_pub_info,"/u/ndesai/gene_pub_all26_info.csv", sep=",")









