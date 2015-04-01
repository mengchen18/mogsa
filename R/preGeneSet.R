prepMsigDB <- function(file) {
  gmt <- readLines(file)
  gmt <- lapply(gmt, strsplit, split="\t")
  gmt <- lapply(gmt, "[[", 1)
  gsn <- sapply(gmt, "[", 1)
  names(gmt) <- gsn
  gs <- lapply(gmt, "[", -(1:2))
  return(gs)
}

prepGraphite <- function(db, id = c("entrez", "symbol")) {
  if (id %in% "symbol") {
    cat("converting identifiers!")
    db <- lapply(db, convertIdentifiers, type="symbol")
    cat("converting identifiers done!")
    gs <- lapply(db, nodes)
  } else if (id %in% "entrez") {
    database <- slot(db[[1]], name="database")
    if (database %in% c("biocarta", "kegg")) {
      gs <- lapply(db, nodes)
      gs <- lapply(gs, function(x) gsub("EntrezGene:", "", x))
    } else if (database %in% c("reactome", "nci")) {
      cat("converting identifiers!")
      db <- lapply(db, convertIdentifiers, type="entrez")
      cat("converting identifiers done!")
      gs <- lapply(db, nodes)
    } else 
      stop("unknown database")
  } else
    stop("unknow identifiers selected")
  return (gs)
}


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
    r <- lapply(r, function(x) x[, i])
  
  } else
    stop("unknown type of X")

  if (sum(i) < 1)
    stop("No geneset annotates the data, did you use the correct identifier?")


  return (r)
}







