genome_ann <- function(ann.gr, TE_file.f) {
  
  ##################################################
  
  ######################## annotate genes body and promoter regions from 'TAIR10' or 'TAIR10.1'
  ann.gr$type = as.character(ann.gr$type)
  
  ### annotate genes, CDS and promoters
  Genes = ann.gr[which(ann.gr$type == "gene")]
  CDS = ann.gr[which(ann.gr$type == "CDS")]
  
  Promoters = promoters(Genes, upstream=2000, downstream=0, use.names=TRUE)
  Promoters$type = "promoter"
  
  ### make 'TxDb' object and 'overlap' function - for intron and UTRs annotations (if needed)
  ###   ###    ###   ###    ###   ###    ###   ###    ###   ###  
  ###   ###  make just if its proper 'gff' file  ###   ### 
  if (F) {#any(!grepl("intron|_prime_UTR$", ann.gr$type))) {
    gff3_TxDb = makeTxDbFromGRanges(ann.gr)
    TxDb_overlap <- function(ann_TxDb, type_name) {
      # overlap to get 'intron' ranges
      m <- findOverlaps(ann_TxDb, ann.gr)
      x <- ann_TxDb[queryHits(m)]
      mcols(x) <- cbind.data.frame(
        mcols(x),
        mcols(ann.gr[subjectHits(m)]))
      mcols(x)$type = type_name
      
      return(unique(x))
    }
    
  }
  
  ### annotate 'introns'
  if (any(grepl("intron", ann.gr$type))) {
    Introns = ann.gr[which(ann.gr$type == "intron")]
  } else {
    Introns_0 = unlist(intronsByTranscript(gff3_TxDb))[,!1:3]
    Introns = TxDb_overlap(Introns_0, "intron")
  }
  
  ### annotate 'UTRs'
  if (any(grepl("_prime_UTR$", ann.gr$type))) {
    fiveUTRs = ann.gr[which(ann.gr$type == "five_prime_UTR")]
    threeUTRs = ann.gr[which(ann.gr$type == "three_prime_UTR")]
  } else {
    fiveUTRs_0 = unlist(fiveUTRsByTranscript(gff3_TxDb))[,!1:3]
    threeUTRs_0 = unlist(threeUTRsByTranscript(gff3_TxDb))[,!1:3]
    fiveUTRs = TxDb_overlap(fiveUTRs_0, "five_prime_UTR")
    threeUTRs = TxDb_overlap(threeUTRs_0, "three_prime_UTR")
  }
  
  ### annotate 'TEG' and 'pseudogene'
  contains_TEG_n_pseudogene <<- any(grepl("^transposable_element_gene$|^pseudogene$", ann.gr$type))
  if (contains_TEG_n_pseudogene) {
    TEG = ann.gr[which(ann.gr$type == "transposable_element_gene")]
    pseudogene = ann.gr[which(ann.gr$type == "pseudogene")]
  } 
  ##################################################
  annotation_vec = list(Genes = Genes,
                        Promoters = Promoters,
                        CDS = CDS,
                        Introns = Introns,
                        fiveUTRs = fiveUTRs,
                        threeUTRs = threeUTRs,
                        Transposable_Elements = TE_file.f)
  
  if (contains_TEG_n_pseudogene) {
    annotation_vec$TEG <- TEG
    annotation_vec$pseudogene <- pseudogene
  }
  
  return(annotation_vec)
}