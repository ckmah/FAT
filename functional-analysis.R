library(cummeRbund)
library(gage)
library(pathview)
library(preprocessCore)
library(clusterProfiler)

keggAnalysis <- function(fpkmData, cuff, kegg.gs, ref) {

  ids <- rownames(fpkmData)
  geneSet <- getGenes(cuff, ids)
  geneIDs <- featureNames(geneSet)
  
  # select genes that map and are not duplicates
  isUniqueMap <- !is.na(geneIDs[, 2]) & !duplicated(geneIDs[, 2])
  mapped <- fpkmData[isUniqueMap,]
  # get entrez ids
  rownames(mapped) <- geneIDs[isUniqueMap, 2]
  
  # reference and sample column indices
  fpkmData.ref <- ref
  # fpkmData.samp <- 2
  
  # convert from gene ID to entrez id
  entrezIDs <- pathview::id2eg(rownames(mapped), category = "symbol")
  entrezIDs.sel <- !is.na(entrezIDs[, 2])
  # select converted genes
  mapped <- mapped[entrezIDs.sel,]
  rownames(mapped) <- entrezIDs[entrezIDs.sel, 2]
  
  return(list("mapped" = mapped, "kegg.gage" = gage(mapped, kegg.gs, sam.dir = FALSE)))
  
  # gs=unique(unlist(kegg.gs[rownames(keggs.p$greater)[1:3]]))
  
  # cnts.d = cnts.norm[, samp.idx] - rowMeans(cnts.norm[, ref.idx])
  # sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$greater[, "q.val"])
  # path.ids <- rownames(cnts.kegg.p$greater)[sel]
  # sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$less[, "q.val"])
  # path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
  # path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
  # library(pathview)
  # pv.out.list <- sapply(path.ids2, function(pid)
  #   pathview(
  #     gene.data = cnts.d,
  #     pathway.id = pid,
  #     species = "hsa"
  #   ))
  #
}
#
# # extract gene names
# gnames = cuff.res$gene
#
# # select genes present in cuffdiff output
# sel = gnames != "-"
# gnames = as.character(gnames[sel])
# cuff.fc = cuff.fc[sel]
# names(cuff.fc) = gnames
#
# # convert to entrez gene id
# gnames.eg = pathview::id2eg(gnames, category = "symbol")
#
# # filter for genes with > 0 fold change
# sel2 = gnames > ""
# cuff.fc = cuff.fc[sel2]
#
# names(cuff.fc) = gnames[sel2]
# range(cuff.fc)
#
# # max of 10 fold change
# cuff.fc[cuff.fc > 10] = 10
# cuff.fc[cuff.fc < -10] = -10
# exp.fc = cuff.fc
# out.suffix = "cuff"
#
