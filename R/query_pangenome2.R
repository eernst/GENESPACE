query_pangenome2 <- function(gsParam = NULL,
                             refGenome = NULL,
                             path2pgTxt = NULL,
                             pullEntriesWithTheseIDs = NULL,
                             onlyAnchorIDs = FALSE,
                             pullEntriesInTheseRegions = NULL,
                             mustHaveRefAnchor = FALSE,
                             directSynOnly = TRUE){

  # -- parameter checking
  if(is.null(gsParam) & is.null(path2pgTxt))
    stop("either gsParam or path2pgTxt file needs to be specified\n")

  pgColNames <- c("pgChr", "pgOrd", "pgID", "og", "ofID", "isRep", "genome",
                  "isDirectSyn", "isSynOgOnly", "isArrayRep", "isNSOrtho", "id")
  if(!is.null(path2pgTxt)){
    if(!file.exists(path2pgTxt)){
      if(!is.null(gsParam)){
        warning(sprintf("path2pgTxt %s does not exist, using path in gsParam\n", path2pgTxt))
      }else{
        stop(sprintf("path2pgTxt %s does not exist\n", path2pgTxt))
      }
    }
  }else{
    pgfs <- list.files(
      file.path(gsParam$paths$results),
      pattern = "pangenomeDB.txt.gz$")

    if(length(pgfs) == 0)
      stop(sprintf(strwrap(
        "can't find any files ending with pangenomeDB.txt.gz in %s. Have you run pangenome() yet?\n",
        gsParam$paths$results)))

    rgs <- gsub("_pangenomeDB.txt.gz", "", basename(pgfs), fixed = T)

    if(is.null(refGenome)){
      refGenome <- gsParam$genomes$genomeIDs
      refGenome <- refGenome[refGenome %in% rgs][1]
    }

    if(!refGenome %in% rgs)
      stop(sprintf(strwrap(
        "can't find files ending with %s_pangenomeDB.txt.gz in %s. Have you run pangenome(refGenome = %s) yet?\n",
        refGenome, gsParam$paths$results, refGenome)))

    path2pgTxt <- file.path(
      gsParam$paths$results,
      sprintf("%s_pangenomeDB.txt.gz", refGenome))
  }

  pg <- fread(path2pgTxt, na.strings = c("", "NA"), showProgress = F)
  if(!all(pgColNames %in% colnames(pg)) || nrow(pg) < 10)
    stop(sprintf(strwrap(
      "path2pgTxt %s contains a file that isn't a properly formatted pangenome source data file\n",
      path2pgTxt)))

  if(!is.null(pullEntriesWithTheseIDs)){
    pgo <- subset(pg, id %in% pullEntriesWithTheseIDs)
    if(nrow(pgo) == 0){
      warning("could not find any entries with ids in pullEntriesWithTheseIDs\n")
      return(NULL)
    }else{
      # -- pull entries with these genes
      if(directSynOnly)
        pgo <- subset(pgo, isDirectSyn)
      if(!includeNSOrthos)
        pgo <- subset(pgo, !isNSOrtho)

      pgtmp <- subset(pg, pgID %in% pgo$pgID)
      pgl <- list(pgtmp)
    }
  }else{
    if(is.null(pullEntriesInTheseRegions)){
      warning("neither regions nor geneIDs are specified - returning the full pan-genome annotation\n")
      # -- pull full pg
      #
      #
      pgl <- list(pgtmp)
    }else{
      bed <- data.table(pullEntriesInTheseRegions)
      if(nrow(bed) == 0)
        stop("could not coerce pullEntriesInTheseRegions into a data.table\n")
      if(!all(c("chr", "start", "end") %in% colnames(bed)))
        stop("pullEntriesInTheseRegions must be a data.frame/table with chr, start and end columns\n")
      gffp <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")

      if(!file.exists(gffp))
        stop("cant find the annotated gff file - has synteny been run with this gsParam object?\n")
      gff <- fread(gffp, na.strings = c("", "NA"), showProgress = F)

      gff <- subset(gff, genome == refGenome)
      if(nrow(gff) == 0)
        stop("something is wrong with the annotated gff file - has synteny been run with this gsParam object?\n")

      bed <- subset(bed, chr %in% gff$chr)
      gff <- subset(gff, ofID %in% pg$ofID)
      if(nrow(bed) == 0)
        stop(strwrap(sprintf(
          "could not find any specified chromosomes (%s) in %s refgenome. Chromsomes available for this refgenome are %s\n"
          paste(bed$chr, collapse = ","),
          refGenome,
          paste(unique(gff$chr), collapse = ", "))))

      pgl <- lapply(1:nrow(bed), function(i){
        tmpGff <- subset(
          gff, chr == bed$chr[i] & end >= bed$start[i] & start <= bed$end[i])
        pgRep <- subset(pg, pgChr == bed$chr[i] & isRep & ofID %in% tmpGff$ofID)
        pr <- range(pgRep$pgOrd)
        out <- subset(pg, pgChr == bed$chr[i] & pgOrd >= pr[1] & pgOrd <= pr[2])
        return(out)
      })
    }
  }

  pgoutl <- lapply(pgl, function(x){
    # -- function to reformat the data
  })

}
