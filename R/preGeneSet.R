#' Conver gmt format file to a list
#' 
#' Convert a gmt file (Could be downloaded from MSigDB) to a list of gene sets
#' information.
#' 
#' 
#' @param file The directory and file name of the gmt file.
#' @return This function returns an object of list containing gene set
#' information, which could be further processed by function "prepSupMoa" to
#' convert to the object that can be used as input of "sup.moa" or "mogsa".
#' @author Chen Meng
#' @seealso See Also as \code{\link{prepGraphite}} and
#' \code{\link{prepSupMoa}}.
#' @export
#' @examples
#' 
#' 	# not run
#' 	dir <- system.file(package = "mogsa")
#' 	preGS <- prepMsigDB(file=paste(dir, 
#' 		"/extdata/example_msigdb_data.gmt.gz", sep = ""))
#' 
prepMsigDB <- function(file) {
  gmt <- readLines(file)
  gmt <- lapply(gmt, strsplit, split="\t")
  gmt <- lapply(gmt, "[[", 1)
  gsn <- sapply(gmt, "[", 1)
  names(gmt) <- gsn
  gs <- lapply(gmt, "[", -(1:2))
  return(gs)
}



#' Prepare pathway gene sets from graphite package
#' 
#' Prepare pathway gene sets from "graphite" package, which could be passed to
#' "prepSupMoa" function.
#' 
#' Only support "entrez" or "symbol" output currently.
#' 
#' @param db The database to be used, an object of class either 'PathwayList'
#' create by "pathways" function.
#' @param id Which identifier for output, either "entrez" or "symbol".
#' @return This function returns an object of list containing gene set
#' information, which could be further processed by function "prepSupMoa" to
#' convert to the object that can be used as input of "sup.moa" or "mogsa".
#' @author Chen Meng
#' @seealso See Also as \code{\link{prepMsigDB}} and \code{\link{prepSupMoa}}.
#' @references Sales G, Calura E and Romualdi C (2014). graphite: GRAPH
#' Interaction from pathway Topological Environment. R package version 1.10.1.
#' @keywords pahtways graphite
#' @export
#' @examples
#' 
#'   library(graphite)
#'   keggdb <- prepGraphite(db = pathways("hsapiens", "kegg")[1:3], id = "entrez")
#' 
prepGraphite <- function(db, id = c("entrez", "symbol")) {

  if (!inherits(db, c("PathwayList")))
    stop("db should be an object of class either 'PathwayList'.")

  if (id %in% "symbol") {

    cat("converting identifiers!\n")
    suppressMessages(
      db <- lapply(db, convertIdentifiers, to="symbol")
      )
    cat("converting identifiers done!\n")
    gs <- lapply(db, nodes)
  
  } else if (id %in% "entrez") {

    database <- slot(db, name="name")
    if (tolower(database) %in% c("biocarta", "kegg")) {
    
      gs <- lapply(db, nodes)
      gs <- lapply(gs, function(x) gsub("EntrezGene:", "", x))
    
    } else if (tolower(database) %in% c("humancyc", "panther", "reactome", "nci")) {
    
      cat("converting identifiers!")
      suppressMessages(
        db <- lapply(db, convertIdentifiers, to="entrez")
        )
      cat("converting identifiers done!\n")
      gs <- lapply(db, nodes)
    
    } 

  } else
    stop("unknow identifiers selected")
  return (gs)
}




#' Prepare sumpplementary tables for projection by sup.moa or mogsa.
#' 
#' Convert a list of gene set information to a set of sumpplementary tables
#' that can be used as input of function "sup.moa" or "mogsa".
#' 
#' Details here
#' 
#' @param X A matrix/data.frame or a list of matrix/data.frame or a list of
#' character vector. If it is a list of matrix/data.frame, row names of
#' matrix/data.frame will be used to create the projection matrix. Otherwise
#' the charater vectors will used to create the supplementary matirx.
#' @param geneSets Gene sets list or an object of class "GeneSet" or
#' "GeneSetCollection".  A gene set list could be returned by prepGraphite or
#' prepMolsigDB.
#' @param minMatch The minimum match of geneset.
#' @param maxMatch The maximum match genesets.
#' @return A list of matrix could used as supplementary tables by "sup.moa" or
#' "mogsa".
#' @author Chen Meng
#' @seealso See Also as \code{\link{prepGraphite}} and
#' \code{\link{prepMsigDB}}.
#' @export
#' @examples
#' 
#'   library(graphite)
#'   data(NCI60_4arrays)
#'   kegg <- pathways(species = "hsapiens", "kegg")
#'   pw <- c("Purine metabolism", "MAPK signaling pathway")
#'   gss <- prepGraphite(db = kegg[pw], id="symbol")
#'   gss <- lapply(gss, function(x) sub("SYMBOL:", "", x))
#'   sup_data1 <- prepSupMoa(NCI60_4arrays, geneSets=gss)
#'   gene_list <- lapply(NCI60_4arrays, rownames)
#'   sup_data2 <- prepSupMoa(gene_list, geneSets=gss)
#' 
prepSupMoa <- function (X, geneSets, minMatch=10, maxMatch=500) {
  
  if (inherits(geneSets, "GeneSet")) {
      nm <- geneSets@setName
      geneSets <- list(geneSets@geneIds)
      names(geneSets) <- nm
    } else if (inherits(geneSets, "GeneSetCollection")) {
        nm <- sapply(geneSets, function(x) x@setName)
        geneSets <- lapply(geneSets, function(x) x@geneIds)
        names(geneSets) <- nm
      } else if (inherits(geneSets, "list")) {
          geneSets <- geneSets
        } else
          stop ("unknown class of input for geneSets")

  matchgs <- function(genes, geneSets) {
    mat <- sapply(geneSets, function(x) match(genes, x, nomatch=0))
    mat[mat > 0] <- 1  
    return(mat)
  }  

  if (is.character(X)) {

    r <- matchgs(X, geneSets)
    i <- colSums(r) >= minMatch & colSums(r) <= maxMatch
    r <- r[, ]

  } else if (inherits(X, c("data.frame", "matrix"))) {
  
    genes <- rownames(X)
    r <- matchgs(genes, geneSets)
    i <- colSums(r) >= minMatch & colSums(r) <= maxMatch
    r <- r[, i]
  
  } else if  (inherits(X, c("list"))) {

    if (all(sapply(X, is.character)))
      genes <- X else
        genes <- lapply(X, rownames)

    r <- lapply(genes, function(x) {
      r <- matchgs(x, geneSets)
      return (r)} )
    ii <- rowSums(sapply(r, colSums))
    i <- ii >= minMatch & ii <= maxMatch
    r <- lapply(r, function(x) x[, i, drop = FALSE])
  
  } else
    stop("unknown type of X")

  if (sum(i) < 1)
    stop("No geneset annotates the data, did you use the correct identifier?")


  return (r)
}







