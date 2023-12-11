############### WILDERLAB STANDARD BIOINFORMATIC PROCESSING SCRIPT ################
## V2.1.0 compiled by: Shaun Wilkinson 2023-12-11
############### SET RUN NUMBER ##############
runnumber <- "XXX" ## insert sequencing run number here

###############  LOAD PACKAGES AND SYSTEM ENVIRONMENT VARIABLES ########################
library(wildertools)
library(insect)
library(dada2)

############### GET CURRENT RUN TAGS AND PRIMER SEQUENCES ######################
con = connect()
qry <- paste0("SELECT * FROM indexing WHERE Run = '", runnumber, "';")
mtags <- DBI::dbGetQuery(con, qry)
qry <- paste0("SELECT * FROM tags;")
tags <- DBI::dbGetQuery(con, qry)
disconnect(con)

primers <- list(WV = c("GACGAGAAGACCCTWTGGAGC", "CCRYGGTCGCCCCAAC", "35", "200"),
                RV = c("TTAGATACCCCACTATGC", "TAGAACAGGCTCCTCTAG", "35", "200"),
                CI = c("DACWGGWTGAACWGTWTAYCCHCC", "GTTGTAATAAAATTAAYDGCYCCTARAATDGA", "70", "76"),
                LV = c("TCGTGCCAGCCRCCGC", "CATAGTGGGGTATCTAATCCCAGTTTG", "100", "200"),
                ZV = c("GACGAGAAGACCCTWTGGAGC", "GGATTGCGCTGTTATCCCT", "70", "300"),
                TP = c("GGGCAATCCTGAGCCAA", "CCATTGAGTCTCTGCACCTATC", "20", "200"),
                MZ = c("CTTCTTCAGGTGGAACTCCAG", "GTCACCACAAACAGAGACTAAAGCAAGT", "100", "112"),
                BE = c("CCCTGCCHTTTGTACACAC", "CCTTCYGCAGGTTCACCTAC", "50", "200"),
                BU = c("TTGTACACACCGCCC", "CCTTCYGCAGGTTCACCTAC", "50", "200"),
                ZC = c("GGACGATAAGACCCTATAAADCTT", "CGCTGTTATCCCTAAAGTAAYTT", "70", "250"),
                UM = c("GGATTAGATACCCTGGTA", "CCGTCAATTCMTTTRAGTTT", "70", "250"),
                DG = c("TCTTCGGTTGGGGCGAC", "GGATTGCGCTGTTATCCCT", "35", "200"),
                ZC = c("GGACGATAAGACCCTATAAADCTT", "CGCTGTTATCCCTAAAGTAAYTT", "35", "200"),
                EA = c("TATATAATGTTATTGTAACRGCGC", "CCCARCATCAAAGGAATCAAYCA", "50", "250"))


primers <- lapply(primers, function(v) gsub(" ", "", v))
mtags$Ftag_seq <- tags$Sequence[match(mtags$Ftag, tags$TID)]
mtags$Rtag_seq <- tags$Sequence[match(mtags$Rtag, tags$TID)]
mtags$idxs <- paste0(mtags$PrimerSet, "_",mtags$Ftag_seq, "_", mtags$Rtag_seq)



############### DOWNLOAD FASTQ FILE #####
tmpf <- tempfile()
system(paste0("/home/bin/bs list projects > ", tmpf))
tmp <- scan(tmpf, what = "", sep = "\n")
projID <- tmp[grepl(paste0("^\\| ", "WL0", runnumber, " "), tmp)]
projID <- sub(paste0("^\\| ","WL0", runnumber, "[ \\|]+([[:digit:]]+).+"), "\\1", projID)
destdir <- "~/Documents/Runs/"
if(!dir.exists(paste0("~/Documents/Runs/", runnumber, "/fastq"))){
  syscall <- paste0("/home/bin/bs download project -i ",projID," -o ", destdir, runnumber,"/fastq --extension=fastq.gz")
  system(syscall)
}



############### READ IN FASTQ FILE  #####
p2fq <- paste0("~/Documents/Runs/", runnumber, "/fastq/")
p2fq <- paste0(p2fq, dir(p2fq)[nchar(dir(p2fq)) > 30], "/")
p2fq <- p2fq[!grepl("\\.json", p2fq)]
p2fq <- paste0(p2fq, dir(p2fq)[1])
fq <- readFASTQ(p2fq, bin = F)


############### DEMULTIPLEX READS ####################
primersets <- unique(mtags$PrimerSet)
if(!dir.exists(paste0("~/Documents/Runs/", runnumber))) {
  dir.create(paste0("~/Documents/Runs/", runnumber))
}
setwd(paste0("~/Documents/Runs/", runnumber))

quals <- attr(fq, "quality")
attr(fq, "quality") <- NULL

if(!dir.exists(paste0("~/Documents/Runs/", runnumber, "/trimmed"))){
  dir.create(paste0("~/Documents/Runs/", runnumber, "/trimmed"))
}
for(h in seq_along(primersets)){
  dir.create(paste0("~/Documents/Runs/", runnumber, "/trimmed/", primersets[h]))
  primerset <- primersets[h]
  primerinfo <- primers[[match(primerset, names(primers))]]
  up <- primerinfo[1]
  down <- primerinfo[2]
  up <- gsub("I", "N", up)
  down <- gsub("I", "N", down)
  minlength <- as.integer(primerinfo[3])
  maxlength <- as.integer(primerinfo[4])
  ncup <- nchar(up)
  ncdn <- nchar(down)
  up <- disambiguate(up)
  rcdown <- disambiguate(rc(down))
  down <- disambiguate(down)
  hitsh <- grepl(paste0("^.{8}", up, ".+", rcdown, ".+"), fq)
  fp_only <- grepl(paste0("^.{8}", up, ".+"), fq)
  seqsh <- fq[hitsh]
  qualsh <- quals[hitsh]
  myftags <- sort(table(sub("^(.{8}).+", "\\1", seqsh)), decreasing = T)
  myrtags <- sort(table(rc(sub(paste0(".+", rcdown, "(.{8}).*"), "\\1", seqsh))), decreasing = T)
  mtagsh <- mtags[mtags$PrimerSet == primersets[h],]
  mtagsh$allreads <- 0
  mtagsh$passfilter <- 0
  ftags <- unique(mtagsh$Ftag_seq)
  rtags <- unique(mtagsh$Rtag_seq)
  ftaglengths <- nchar(ftags)
  rtaglengths <- nchar(rtags)
  fregexs <- paste0(ftags, up)
  rregexs <- paste0(rcdown, rc(rtags))
  rregexs <- paste0(rregexs, "ATCTCG")
  seqlengths <- integer(0)
  passfilter <- integer(0)
  for(i in seq_along(ftags)){
    hitsi <- grepl(fregexs[i], seqsh)
    if(sum(hitsi) > 0){
      seqstart <- ftaglengths[i] + ncup + 1
      seqsi <- substr(seqsh[hitsi], seqstart, 10000)
      qualsi <- substr(qualsh[hitsi], seqstart, 10000)
      if(TRUE){
        for(j in seq_along(rtags)){
          whichsample <- which(mtagsh$Ftag_seq == ftags[i] & mtagsh$Rtag_seq == rtags[j])
          if(length(whichsample) > 1) {
            whichsample <- whichsample[1]
          }
          if(length(whichsample) > 0 & length(seqsi) > 0){
            gregs <- gregexpr(rregexs[j], seqsi)
            gregs <- vapply(gregs, head, 0L, 1L)
            hitsij <- gregs > 0
            nhits <- sum(hitsij)
            mtagsh$allreads[whichsample] <- nhits
            if(nhits > 0){
              if(!all(gregs == 1L))
                seqsij <- substr(seqsi[hitsij], rep(1, nhits), gregs[hitsij] - 1)
                qualsij <- substr(qualsi[hitsij], rep(1, nhits), gregs[hitsij] - 1)
              }else{
                seqsij <- seqsi[hitsij]
                qualsij <- qualsi[hitsij]
              }
              oks <- nchar(seqsij) >= minlength & nchar(seqsij) <= maxlength
              mtagsh$passfilter[whichsample] <- sum(oks)
              seqlengths <- c(seqlengths, nchar(seqsij))
              passfilter <- c(passfilter, sum(oks))
              seqsij <- seqsij[oks]
              qualsij <- qualsij[oks]
              if(length(seqsij) > 0L){
                reslen <- length(seqsij) * 4
                res <- character(reslen)
                res[seq(1, reslen, by = 4)] <- paste0("@", names(seqsij))
                res[seq(2, reslen, by = 4)] <- seqsij
                res[seq(3, reslen, by = 4)] <- rep("+", length(seqsij))
                res[seq(4, reslen, by = 4)] <- qualsij
                if(!dir.exists(paste0("~/Documents/Runs/", runnumber, "/trimmed/", primersets[h], "/demux"))){
                  dir.create(paste0("~/Documents/Runs/", runnumber, "/trimmed/", primersets[h], "/demux"))
                }
                fpath <- paste0(paste0("~/Documents/Runs/", runnumber, "/trimmed/", primersets[h], "/demux/"),
                                runnumber, "_",
                                primersets[h], "_",
                                format(mtagsh$UID[whichsample], scientific = F), "_",
                                mtagsh$Ftag_seq[whichsample], "_",
                                mtagsh$Rtag_seq[whichsample],".fastq")
                message("writing ", length(seqsij), " sequences to ", fpath)
                cat(res, file = fpath, sep = "\n", append = FALSE)
              }
              seqsi <- seqsi[!hitsij]
              qualsi <- qualsi[!hitsij]
            }
          }
        }
      }
      seqsh <- seqsh[!hitsi]
      qualsh <- qualsh[!hitsi]
  }
}

############### DADA2 FILTERING ###########
primersets <- unique(mtags$PrimerSet)
for(h in seq_along(primersets)){
  setwd(paste0("~/Documents/Runs/", runnumber, "/trimmed/", primersets[h]))
  path <- "demux"
  fnFs <- sort(list.files(path, pattern="", full.names = TRUE))
  continue = T
  if(length(fnFs) > 0){
    sample.names <- basename(sub("\\.fastq", "", fnFs))
    filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    out <- filterAndTrim(fnFs, filtFs,  truncLen=0,
                         maxN=0, maxEE = 2, truncQ = 2, rm.phix=TRUE,
                         compress = TRUE, multithread = 12)
    some_missing <- !all(file.exists(filtFs))
    if(some_missing){
      if(any(file.exists(filtFs))){
        noread.names <- sample.names[!file.exists(filtFs)]
        nfails <- length(noread.names)
        noread.track <- data.frame(input = rep(0, nfails), filtered  = rep(0, nfails), denoisedF  = rep(0, nfails), nonchim  = rep(0, nfails))
        rownames(noread.track) <- noread.names
        sample.names <- sample.names[file.exists(filtFs)]
        fnFs <- fnFs[file.exists(filtFs)]
        filtFs <- filtFs[file.exists(filtFs)]
        out <- filterAndTrim(fnFs, filtFs,  truncLen=0,
                             maxN=0, maxEE = 2, truncQ = 2, rm.phix=TRUE,
                             compress = TRUE, multithread = 12)
      }else{
        continue = FALSE
      }
    }
    errF <- tryCatch(learnErrors(filtFs, multithread=12), error = function(e) return(NULL))
    if(sum(file.size(filtFs)) < 1000 | !continue | is.null(errF)){
      message("Too few sequences for primerset ", primersets[h])
    }else{
      saveRDS(errF, file = "errF.rds")
      derepFs <- derepFastq(filtFs, verbose=F)
      if(names(derepFs[1]) == "uniques") derepFs <- list(tmp = derepFs)
      names(derepFs) <- sample.names
      dadaFs <- dada(derepFs, err=errF, multithread=12, pool = F)
      if(names(dadaFs)[1] == "denoised"){
        dadaFs <- list(tmp = dadaFs)
        names(dadaFs) <- sample.names
      }
      seqtab <- makeSequenceTable(dadaFs)
      saveRDS(seqtab, file = "seqtab.rds")
      seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=12, verbose=TRUE)
      dim(seqtab.nochim)
      sum(seqtab.nochim)/sum(seqtab)
      getN <- function(x) sum(getUniques(x))
      if(length(derepFs) > 1){
        track <- cbind(out, sapply(dadaFs, getN),rowSums(seqtab.nochim))
        colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
        rownames(track) <- sample.names
        write.csv(track, file = "trackprocess.csv", row.names = TRUE)
      }
      alluids <- paste0(runnumber, "_", primersets[h], "_",
                        mtags$UID[mtags$PrimerSet == primersets[h]], "_",
                        mtags$Ftag_seq[mtags$PrimerSet == primersets[h]], "_",
                        mtags$Rtag_seq[mtags$PrimerSet == primersets[h]])
      alluids <- unique(alluids)
      missinguids <- alluids[!alluids %in% rownames(seqtab.nochim)]
      saveRDS(seqtab.nochim, file = "seqtab_nochim.rds")
      myrecords <- vector(mode = "list", length = nrow(seqtab.nochim))
      for(i in seq_len(nrow(seqtab.nochim))){
        tmp <- seqtab.nochim[i, ]
        if(is.null(names(tmp))) names(tmp) <- colnames(seqtab.nochim)
        tmp <- tmp[tmp > 0]
        if(length(tmp) > 0L){
          tmp <- sort(tmp, decreasing = TRUE)
          out <- names(tmp)
          UID <- as.integer(sub(".+_([[:digit:]]{6})_.+", "\\1", rownames(seqtab.nochim)[i]))
          newrecords <- data.frame(UID = rep(UID, length(tmp)))
          newrecords$Run <- runnumber
          newrecords$PrimerSet <- primersets[h]
          newrecords$Ftag <-  sub(".+_([ACGT]{8})_.+", "\\1", rownames(seqtab.nochim)[i])
          newrecords$Rtag <- sub(".+_([ACGT]{8})$", "\\1", rownames(seqtab.nochim)[i])
          newrecords$Sequence <- out
          newrecords$SequenceHash <- paste(openssl::md5(out))
          newrecords$Count <- unname(tmp)
          newrecords$Verified <- 1L
          newrecords$Notes = ""
          myrecords[[i]] <- newrecords
        }
      }
      newrecords <- do.call("rbind", myrecords)
      con = connect()
      test <- DBI::dbWriteTable(con, "records", newrecords, append=T, row.names=F)
      RMySQL::dbDisconnect(con)
    }
    setwd(paste0("~/Documents/Runs/", runnumber))
  }
}



############### TAXON ASSIGNMENT ##########################
con = connect()
qry <- "SELECT TaxID,Sequence,SequenceHash,PrimerSet FROM eDNA;"
eDNA <- DBI::dbGetQuery(con, qry)
eDNA_default <- eDNA[grepl("^00", eDNA$SequenceHash),]
qry <- paste0("SELECT DISTINCT PrimerSet FROM records WHERE Run = '", runnumber, "';")
primersets <- DBI::dbGetQuery(con, qry)
primersets <- primersets$PrimerSet
primersets
qry <- "SELECT * FROM taxa;"
tx <- DBI::dbGetQuery(con, qry)
colnames(tx)[1:4] <- c("taxID", "parent_taxID", "rank", "name")
tdb <- tx[,1:4]
colnames(tdb) <- c("taxID", "parent_taxID", "rank", "name")

## NCBI common names
# tdbcnf <- NCBI_commonnames() ## last update 17 May 2021
# saveRDS(tdbcnf, file = "C:/Users/shaun/Dropbox/Runs/NCBI_commonnames.rds")
tdbcnf <- readRDS("~/Documents/Taxonomy/NCBI_commonnames.rds")

runuids <- paste0(unique(mtags$UID), collapse = "','")
qry <- paste0("SELECT * FROM records WHERE UID IN ('", runuids,"');" )
runrecs <- DBI::dbGetQuery(con, qry)
runrecs <- runrecs[runrecs$Run == runnumber, ]
disconnect(con)

## Assign taxonomy and update eDNA table
for(h in seq_along(primersets)){
  seqstoadd <- runrecs[runrecs$PrimerSet == primersets[h], ]
  seqstoadd <- unique(seqstoadd$Sequence[!seqstoadd$Sequence %in% eDNA$Sequence])
  refseqs_available <- TRUE
  if(!file.exists(paste0("~/Documents/Refseqs/", primersets[h], ".rds"))){
    con = connect()
    qry <- paste0("SELECT * FROM refseqs WHERE PrimerSet ='", primersets[h],"';")
    rs <- DBI::dbGetQuery(con, qry)
    RMySQL::dbDisconnect(con)
    if(nrow(rs) > 0){
      rs <- rs[rs$Use == 1L, ]
      z <- char2dna(rs$Sequence)
      names(z) <- paste0(rs$Accession, "|", rs$TaxID)
      message("Found ",length(z)," reference sequences for primer set ", primersets[h])
      taxids <- as.integer(sub(".+\\|", "", names(z)))
      lineages <- insect::get_lineage(taxids, db = tdb, numbers = TRUE, simplify=F)
      lineages <- vapply(lineages, paste0, "", collapse = "; ")
      names(lineages) <- sub("\\|.+", "", names(z))
      attr(z, "lineages") <- lineages
      message("Saving classifier trainingset for assay ", primersets[h])
      saveRDS(z, file = paste0("~/Documents/Refseqs/", primersets[h], ".rds"))
      message("writing sintax file for primerset ",  primersets[h])
      sintax_headers <- format_sintax(lineages)
      tsx <- z
      names(tsx) <- sintax_headers
      writeFASTA(tsx, file = paste0("~/Documents/Refseqs/",  primersets[h], ".fa"))
    }else{
      refseqs_available <- FALSE
      warning("No refseqs available in database for primerset ", primersets[h])
    }
  }
  if(length(seqstoadd) > 0L){
    if(!refseqs_available){
      message("Adding new sequences to eDNA table with root TaxID")
      new_eDNA <- data.frame(TaxID = 1,
                             Sequence = seqstoadd,
                             SequenceHash = sub("^..", "00", paste(openssl::md5(seqstoadd))),
                             PrimerSet = primersets[h],
                             Curated = 0,
                             Deposited = 0,
                             Predicted = 0,
                             Notes = paste0("WL0", runnumber)
      )
      message("pushing eDNA to Aurora")
      con = connect()
      test <- DBI::dbWriteTable(con, "eDNA", new_eDNA, append = TRUE, row.names=FALSE)
      disconnect(con)
    }else{
      z <- readRDS(file = paste0("~/Documents/Refseqs/", primersets[h], ".rds"))
      message("Found ",length(z)," reference sequences for primer set ", primersets[h])
      if(is.null(attr(z, "lineages"))){
        taxids <- as.integer(sub(".+\\|", "", names(z)))
        if(!all(taxids %in% tdb$taxID)){
          warning("some taxids are missing from taxon database - removing affected reference seqs")
          z <- z[taxids %in% tdb$taxID]
          taxids <- taxids[taxids %in% tdb$taxID]
        }
        lineages <- insect::get_lineage(taxids, db = tdb, numbers = TRUE, simplify=F)
        lineages <- vapply(lineages, paste0, "", collapse = "; ")
        names(lineages) <- sub("\\|.+", "", names(z))
        attr(z, "lineages") <- lineages
        message("Overwriting classifier trainingset for assay ", primersets[h])
        saveRDS(z, file = paste0("~/Documents/Refseqs/", primersets[h], ".rds"))
      }
      newtaxids <- classify2(char2dna(seqstoadd), z, tdb, oneMM = F)
      if(any(nchar(newtaxids$sequence) > 50 & newtaxids$taxID == 1L)){
        sxf <-paste0("~/Documents/Refseqs/", primersets[h], ".fa")
        if(file.exists(sxf)){
          thirdpass <- sintax(char2dna(newtaxids$sequence[nchar(newtaxids$sequence) > 50 & newtaxids$taxID == 1L]),sxf, tdb, threshold = 0.99)
        }else{
          thirdpass <- sintax(char2dna(newtaxids$sequence[nchar(newtaxids$sequence) > 50 & newtaxids$taxID == 1L]),z, tdb, threshold = 0.99)
        }
        ranks <- tx$rank[match(thirdpass, tx$taxID)]
        parent_taxids <- tx$parent_taxID[match(thirdpass, tx$taxID)]
        # Don't classify to species level with sintax
        thirdpass[ranks == "species"] <- parent_taxids[ranks == "species"]
        newtaxids$predicted[nchar(newtaxids$sequence) > 50 & newtaxids$taxID == 1L] <- TRUE
        newtaxids$taxID[nchar(newtaxids$sequence) > 50 & newtaxids$taxID == 1L] <- thirdpass
      }
      newtaxids$curated <- FALSE
      newtaxids$deposited <- newtaxids$sequence %in% dna2char(z)
      newtaxids <- newtaxids[c("sequence", "taxID", "curated", "deposited", "predicted")]
      newlineages <- insect::get_lineage(newtaxids$taxID, tdb)
      newlins <- sapply(newlineages, paste0, collapse = ";")
      taxnames <- tx$name[match(newtaxids$taxID, tx$taxID)]
      if(nrow(newtaxids) > 0L){
        colnames(newtaxids) <- c("Sequence", "TaxID", "Curated", "Deposited", "Predicted")
        seqhashes <- data.frame(SequenceHash = sub("^..", "00",paste(openssl::md5(newtaxids$Sequence))),
                                PrimerSet = rep(primersets[h], nrow(newtaxids)))
        new_eDNA <- cbind(seqhashes, newtaxids)
        new_eDNA$Notes <- paste0("WL0", runnumber)
        new_eDNA$Curated <- as.integer(new_eDNA$Curated)
        new_eDNA$Deposited <- as.integer(new_eDNA$Deposited)
        new_eDNA$Predicted <- as.integer(new_eDNA$Predicted)
        con = connect()
        test <- DBI::dbWriteTable(con, "eDNA", new_eDNA, append = TRUE, row.names=FALSE)
        RMySQL::dbDisconnect(con)
      }
    }
  }
}


############### RESET TO CLEAR MEMORY ################
rm(list = ls())
.rs.restartR()
############### DONE #################

