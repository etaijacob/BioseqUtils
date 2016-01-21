#BioseqUtils.R
#Copyright (C) 2015 Etai Jacob, etai.jacob@gmail.com

clusterMSABySVD <- function(inputfile, energy.th = 0.7, maxNumOfClusters=12, write2file=F) {
  S <- return_alignment_as_matrix(inputfile = inputfile, letter_type = "asstring")
  M <- return_alignment_as_matrix(inputfile = inputfile, letter_type = "binary")

  mysvd = svd(scale(M, center = TRUE, scale = F))
  energy <- sapply(1:length(mysvd$d), function(x) sum(mysvd$d[1:x]^2)/sum(mysvd$d[1:length(mysvd$d)]^2))
  nu = max(which(energy <= energy.th))
  cat(sprintf("Taking the first %d components out of %d.\n", nu, length(energy)))
  u <- mysvd$u[,1:nu]
  rownames(u) <- rownames(M)

  dists <- dist(u)*dim(M)[1]

  fit <- hclust(dists, method="ward.D")


  df <- data.frame(ID = rownames(M), stringsAsFactors = F)
  for(k in 2:maxNumOfClusters) {
    groups <- cutree(fit, k=k)
    df[, sprintf("groups.%d", k)] <- groups[df$ID]
  }
  df <- cbind(df, sequence = S[df$ID], stringsAsFactors = F)
  if(write2file)
    write.table(df, file=sprintf("%s.clusters.e%f.csv", inputfile, energy.th), row.names = F, sep=",")

  return(df)
}

#maualFilterFile first line is a header. This header is ignored. First column is seq name and second is zero or one.
filter_MSA <- function(inputMSAfile, inputMSAfileCommentSep = "|",
                       manualFilterFile = NA,
                       filteringMethod = c("removeDuplicates", "manualFilterFile"),
                       fileType = "fasta") {
  if(filteringMethod == "removeDuplicates")
    return(return_alignment(inputfile = inputMSAfile, doUnique = T, fileType = fileType))
  if(filteringMethod == "manualFilterFile" & !is.na(manualFilterFile)) {
    flt <- read.table(file = manualFilterFile, sep = "\t", stringsAsFactors = F, header = T)
    msa <- return_alignment(inputfile = inputMSAfile, fileType = fileType)
    if(!is.na(inputMSAfileCommentSep)) {
      mynames <- as.vector(sapply(msa$nam,
                                  function(s) gsub(sprintf("^\\s+%s\\s+$", inputMSAfileCommentSep), "",
                                                   strsplit(x = s, sprintf("\\%s", inputMSAfileCommentSep))[[1]][1])))
      msa$nam <- mynames
    }
    inidxs <- which((msa$nam %in%  flt[flt[,2] == 1, 1]))
    msa$nb  <- length(msa$nb[inidxs])
    msa$nam <- msa$nam[inidxs]
    msa$seq <- msa$seq[inidxs]
    msa$com <- msa$com[inidxs]

    return(msa)
  }
}

write2Fasta <- function(msa, foutname) {
  write.fasta(sequences = strsplit(msa$seq, split = ""), names = msa$nam, file.out = foutname)

}

return_alignment_as_matrix <- function(inputfile,
                                       letter_type = c("letter", "number", "binary", "asstring")) {
  letter_type <- match.arg(letter_type)
  msa <- return_alignment(inputfile = inputfile)
  mynames <- as.vector(sapply(msa$nam,
                       function(s) gsub("^\\s+|\\s+$", "", strsplit(x = s, "\\|")[[1]][1])))
  if(letter_type == "number") {
    M <- do.call(rbind, lapply(strsplit(toupper(msa$seq), ""), function(x) factor(x, levels = aas)))
  } else if(letter_type == "letter") {
    M <- do.call(rbind, lapply(msa$seq, function(x) strsplit(toupper(x), "")[[1]]))
  } else if(letter_type =="binary") {
    M <- do.call(rbind, lapply(strsplit(toupper(msa$seq), ""), function(x) factor(x, levels = aas)))
    M <- aa_msa_to_binary_matrix(M)
  } else {
    M <- msa$seq
    names(M) <- mynames
    return(M)
  }

  rownames(M) <- mynames
  return(M)

}

return_alignment <- function(inputfile,
                             nuc = F, fileType = "fasta", doUnique=F) {

  if(fileType == "fasta") {
    msa <- seqinr::read.alignment(inputfile, format = "fasta", forceToLower = F)
    cat(sprintf("Original MSA length is: %d.\n", length(msa$seq)))
    if(doUnique) {
      inidxs <- which(!duplicated(msa$seq))
      msa$nb  <- length(msa$nb[inidxs])
      msa$nam <- msa$nam[inidxs]
      msa$seq <- msa$seq[inidxs]
      msa$com <- msa$com[inidxs]
    }
  }
  return(msa)
}

#reads a stockholm file format
# only lines which do not start with a # are considered
# returns a dataframe with name and seq
read.stockholm.alignment <- function(file) {
  mm <- readLines(file)
  mm <- do.call(rbind, strsplit(mm[-grep("^#", mm)], "\\s+", perl = T))
  mm <- data.frame(name=mm[,1], seq=mm[,2], row.names = mm[,1], stringsAsFactors = F)
  return(mm)
}



aa_msa_to_binary_matrix <- function(pmat) {
  tt <- lapply(1:dim(pmat)[1], function(x) aa_seq_to_binary_vec(pmat[x,]))
  tt <- do.call(rbind, tt)
  rownames(tt) <- rownames(pmat)
  return(tt)
}

aa_seq_to_binary_vec <- function(sequence) {
  unlist(lapply(sequence, aa_letter_to_binary_vec))
}

aa_letter_to_binary_vec <- function(aa) {
  (1:21==aa)+0
}

aas <- toupper(c(".","A", "C", "D", "E", "F", "G", "H", "I",
                 "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
                 "Y", "-", "X", "B"))

