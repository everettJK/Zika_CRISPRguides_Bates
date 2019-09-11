library(ShortRead)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(tidyverse)
library(ggrepel)

options(stringsAsFactors = FALSE)
CPUs <- 30
report <- list()

# Capture session info.
write(capture.output(sessionInfo()), file = 'sessionInfo.txt')

#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

# [ Manual step ]
# Join together R1 and R2 reads with PEAR software.
# source ~/miniconda3/bin/activate sequenceTools
# pear -f Undetermined_S0_L001_R1_001.fastq -r Undetermined_S0_L001_R2_001.fastq -o pearOutput -j 20 -v 40 -m 160 -n 100

# Read in sequencing data.
# (!) This data needs to be exported to a public resource and pull in when this analysis in published.

reads      <- readFastq('data.large/180301_M03249_0321_000000000-BNDFG/pearOutput.assembled.fastq')
indexReads <- readFastq('data.large/180301_M03249_0321_000000000-BNDFG/Undetermined_S0_L001_I1_001.fastq')
report$totalReads <- length(reads)


# Match ids between the PEAR concatenated reads and the index reads reorder/subset the index reads.
readsIDs      <- as.character(ShortRead::id(reads))
indexReadsIDs <- as.character(ShortRead::id(indexReads))

i <- match(readsIDs, indexReadsIDs)
if(any(is.na(i))) stop('The index reads could not be matched to the target reads.')
indexReads <- indexReads[i]

# Sanity check
table(as.character(ShortRead::id(reads)) == as.character(ShortRead::id(indexReads)))


#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

# Now the reads have been read in and ordered so that the index reads match the target reads.

# Based on manually aligning reads to one another, these sequence boundaries demark the start
# and end of of the CRSPR guide sequences: CGAAACACC -> NNNNNN <- GTTTTAGAG 

# The positions of the anchor sequences are not fixed and we need to search for them.
# Determine the boundaries of the anchoring sequences allowing for 1 mismatch.

boundary1 <- vmatchPattern('CGAAACACC', sread(reads), max.mismatch=0)
boundary2 <- vmatchPattern('GTTTTAGAG', sread(reads), max.mismatch=0)

# For each indentified boundary range, return the sequence positions for the extraction 
# of the CRISPR sequence. NA will be returned for sequences where the anchor positions 
# were not found so that the returned vectors will be the same length as the read object.

boundary1.position <- suppressWarnings(unlist(lapply(endIndex(boundary1), min))) + 1
boundary2.position <- suppressWarnings(unlist(lapply(startIndex(boundary2), max))) - 1

# Correct for max and min acting on NULL
boundary1.position[boundary1.position == Inf]  <- NA
boundary2.position[boundary2.position == -Inf] <- NA

# Create a save point.
save.image(file='savePoints/SP1.RData')


#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~



# Create an index vector to subset reads where both boundaries were found.
indexBothBoundariesFound <- base::intersect(which(! is.na(boundary1.position)), which(! is.na(boundary2.position)))

# Now use this index vector to subset the two boundary1, boundary2, reads, and index read vectors.

# Subset the boundary vectors so that only instances where both boundaries remain.
boundary1.position.both <- boundary1.position[indexBothBoundariesFound]
boundary2.position.both <- boundary2.position[indexBothBoundariesFound]

# Subset the reads so that only reads with both boundaries remain.
reads.withBothBoundaries      <- reads[indexBothBoundariesFound]
indexReads.withBothBoundaries <- indexReads[indexBothBoundariesFound]

# Remove instances where boundary 2 comes before boundary 1. (just a handful of cases)
i <- which(boundary1.position.both >= boundary2.position.both)

boundary1.position.both       <- boundary1.position.both[-i]
boundary2.position.both       <- boundary2.position.both[-i]

reads.withBothBoundaries      <- reads.withBothBoundaries[-i]
indexReads.withBothBoundaries <- indexReads.withBothBoundaries[-i]


# Determine the extent to which the boundaries where found in the reads.
# These values will be used in the report.
report$percentBoundary1Found      <- sprintf("%0.2f", (length(which(! is.na(boundary1.position))) / length(boundary1.position)*100))
report$percentBoundary2Found      <- sprintf("%0.2f", (length(which(! is.na(boundary2.position))) / length(boundary2.position)*100))
report$percentBothBoundariesFound <- sprintf("%0.2f", (length(indexBothBoundariesFound) / length(boundary1.position))*100)


# Read in demultiplexing information.
sampleInfo <- read.table(file = 'data/sampleInfo.tsv', header = TRUE, sep = '\t')


# Read in expected guideSeq sequences.
# All of the guide sequences in this file are unique.
expectedGuides <- read.csv('data/librarySequences.csv', header = TRUE)


# Use biomart to create a data.frame of the transcript ids (NM_xxx) that were provided with the 
# expected guide sequences and associated gene symbols then use this table to add gene symbols 
# to the expectedGuides data frame.

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
o <- getBM(filters="refseq_mrna", 
           attributes=c("refseq_mrna", "refseq_mrna_predicted", "hgnc_symbol"), 
           values=unique(expectedGuides$gene), mart=mart, uniqueRows=FALSE)

expectedGuides$symbol <- o[match(expectedGuides$gene, o$refseq_mrna),]$hgnc_symbol

expectedGuides[which(is.na(expectedGuides$symbol)),]$symbol <- expectedGuides[which(is.na(expectedGuides$symbol)),]$gene
expectedGuides[which(expectedGuides$symbol == ''),]$symbol <- expectedGuides[which(expectedGuides$symbol == ''),]$gene

# Create a data frame where each row represents a single read.
# Barcodes and read heads are added in order to check demultiplexing.

s <- subseq(sread(reads.withBothBoundaries), boundary1.position.both, boundary2.position.both)

guideSeqs <- data.frame(conditionBarCode  = as.character(sread(indexReads.withBothBoundaries)),
                        readHead          = subseq(sread(reads.withBothBoundaries), 1, 20),
                        guideSeq          = s,
                        guideSeqLength    = width(s),
                        condition         = NA,
                        replicate         = NA)


# Test each read for each of the 10 replicate bar codes.
# Each test will return the number of matches and the results are bound into a table 
# as shown below. Rows that sum to 1 have a single match to the bar code table.
# Row sums of more than 1 indicate an ambigous match and row sums of 0 indicate no matches.
# This proceadure is then repeated for the condition bar codes.
#
# read  bc1  bc2  bc3  bc4  bc5  bc6  bc7  bc8  bc9  bc10
# 1     0    0    0    1    0    0    0    0    0    0
# 2     0    0    0    1    0    0    0    0    0    0
# 3     0    0    0    0    1    0    0    0    0    0
# 4     1    0    0    0    0    0    0    0    0    0
# 5     0    0    0    0    0    1    0    0    0    0

demultiplex <- function(s, barcodes){
  m <- lapply(barcodes, function(x){ vcountPattern(x, s, max.mismatch=1) })
  a <- do.call(cbind, m)
  apply(a, 1, function(x){ ifelse(sum(x) == 1, which(x == 1), 0) })
}

replicates <- dplyr::select(sampleInfo, replicate, barcode1) %>% dplyr::distinct()
conditions <- dplyr::select(sampleInfo, condition, barcode2) %>% dplyr::distinct()

guideSeqs$replicate <- demultiplex(subseq(sread(reads.withBothBoundaries), 1, 25), replicates$barcode1)
guideSeqs$condition <- demultiplex(sread(indexReads.withBothBoundaries), conditions$barcode2)

guideSeqs$conditionName <- 'Unknown'
guideSeqs[guideSeqs$condition != 0,]$conditionName <- conditions[guideSeqs[guideSeqs$condition != 0,]$condition,]$condition

# Create a save point
save.image(file='savePoints/SP2.RData')


# Sanity check
#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

# Should resolve to a single (correct) condition.
table(subset(guideSeqs, conditionBarCode == 'GTGCGTAA')$conditionName)
table(subset(guideSeqs, conditionBarCode == 'CTATTCAA')$conditionName)

# Should all appear to be unknown
head(sort(table(subset(guideSeqs, conditionName == 'Unknown')$conditionBarCode), decreasing = TRUE), n = 10)

# Should all resemble Zika condition barcode with 1 mismatch
head(sort(table(subset(guideSeqs, conditionName == 'Zika' & conditionBarCode != 'CTATTCAA')$conditionBarCode), decreasing = TRUE))

#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~






#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~


# Statistics on how many reads map to the expected sequences.
ppPercent <- function(n) paste0(sprintf("%.2f", n*100), '%')

stats <- summarize(guideSeqs, 
                   readsDemultiplexed       = ppPercent(sum(conditionName != 'Unknown') / n()),
                   expected.percentReads    = ppPercent(sum(guideSeq %in% expectedGuides$guideSeq) / n()),
                   expected.percentReads.20 = ppPercent(sum(guideSeq[guideSeqLength == 20] %in% expectedGuides$guideSeq) / length(guideSeq[guideSeqLength==20])),
                   expected.percentSeqs     = ppPercent(sum(unique(guideSeq) %in% expectedGuides$guideSeq) / n_distinct(guideSeq)),
                   expected.percentSeqs.20  = ppPercent(sum(unique(guideSeq[guideSeqLength == 20]) %in% expectedGuides$guideSeq) / n_distinct(guideSeq[guideSeqLength == 20])))


# Now that we have demultiplexed the reads, we take a look at the distirubtions of guide sequence
# lenghts and their experimental conditions.

d <- group_by(guideSeqs, guideSeqLength, conditionName) %>%
     summarise(n=n()) %>%
     ungroup() %>% 
     mutate(conditionName = factor(conditionName, levels = c("Unknown", "ParLib", "UnInf", "Zika"))) %>%
     filter(guideSeqLength <= 125) %>% data.frame()
             
guideSeqLengthPlot1 <- ggplot(d, aes(guideSeqLength, n, fill=conditionName)) + 
  theme_bw() +
  scale_fill_manual(values=c('gray50', 'dodgerblue2', 'green3', 'red1')) + 
  geom_bar(stat='identity') + 
  xlim(c(0, 40)) +
  labs(x='Length of excised guide sequence', y='Number of reads') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Set all guide sequencs with length 16 -> 21 to zero in order to show the distribution tails.
d[which(d$guideSeqLength %in% c(16:21)),]$n <- 0

guideSeqLengthPlot2 <- ggplot(d, aes(guideSeqLength, n, fill=conditionName)) + 
  theme_bw() +
  scale_fill_manual(values=c('gray50', 'dodgerblue2', 'green3', 'red1')) + 
  geom_bar(stat='identity') + 
  xlim(c(0, 40)) +
  labs(x='Length of excised guide sequence', y='Number of reads') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

#
#  From this point onward, we are only focusing on the guide sequences that are 20 NT long.
#

# Subset the guideSeq data frame such that it only contains expected guide sequences.
guideSeqs <- subset(guideSeqs, guideSeq %in% expectedGuides$guideSeq)


# Create a rarefaction plot data.
s <- split(guideSeqs,  guideSeqs$conditionName)
rarificationData <- bind_rows(lapply(s, function(x){
  reads <- table(x$guideSeq)
  rarificationSteps <-seq(from=1, to=sum(reads), by=10000)
  y <- sapply(rarificationSteps, function(x){ as.numeric(vegan::rarefy(reads, x)) })
  data.frame(sample=x$conditionName[1], x=rarificationSteps, y=y)
}))

rarificationData$sample <- factor(as.character(rarificationData$sample), levels = c('UnInf', 'Zika', 'Unknown'))

# Create a rarefaction plot.
rarificationPlot <- ggplot(rarificationData , aes(x, y, group=sample, color=sample)) +
  geom_line(aes(color=sample), size=0.75) +
  scale_color_manual(values=c('dodgerblue2', 'green3', 'gray50')) + 
  scale_x_continuous(labels=scales::comma) +
  scale_y_continuous(labels=scales::comma) +
  labs(x='Reads', y='Unique guide sequences') +
  theme_bw() +
  guides(color=guide_legend(override.aes=list(size=5))) + theme(legend.title=element_blank())



#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~


guideSeqs$gene <- expectedGuides[match(guideSeqs$guideSeq, expectedGuides$guideSeq),]$symbol


# Transform the reads data frame where each row is a single read to a read tally data frame
# where each row is a unique guide sequences and the number of reads for each condition is summed.

guideSeqs.seqAgg            <- dcast(guideSeqs, guideSeq~conditionName, fun.aggregate = length, value.var = 'conditionName', fill = 0)
guideSeqs.seqAgg$gene       <- expectedGuides[match(guideSeqs.seqAgg$guideSeq, expectedGuides$guideSeq),]$symbol
guideSeqs.seqAgg$totalUnInf <- sum(guideSeqs.seqAgg$UnInf)
guideSeqs.seqAgg$totalZika  <- sum(guideSeqs.seqAgg$Zika)
guideSeqs.seqAgg$UnInfFreq  <- guideSeqs.seqAgg$UnInf / sum(guideSeqs.seqAgg$UnInf)
guideSeqs.seqAgg$ZikaFreq   <- guideSeqs.seqAgg$Zika / sum(guideSeqs.seqAgg$Zika)
guideSeqs.seqAgg$label1     <- ''
guideSeqs.seqAgg$label2     <- ''




guideSeqs.geneAgg            <- dcast(guideSeqs, gene~conditionName, fun.aggregate = length, value.var = 'conditionName', fill = 0)
guideSeqs.geneAgg$totalUnInf <- sum(guideSeqs.geneAgg$UnInf)
guideSeqs.geneAgg$totalZika  <- sum(guideSeqs.geneAgg$Zika)
guideSeqs.geneAgg$UnInfFreq  <- guideSeqs.geneAgg$UnInf / sum(guideSeqs.geneAgg$UnInf)
guideSeqs.geneAgg$ZikaFreq   <- guideSeqs.geneAgg$Zika / sum(guideSeqs.geneAgg$Zika)
guideSeqs.geneAgg$label1     <- ''
guideSeqs.geneAgg$label2     <- ''


#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~



# For each aggregated row (counts per guideSeq), perform one sided Fisher's Exact test. 
guideSeqs.seqAgg.tests <- bind_rows(lapply(1:nrow(guideSeqs.seqAgg), function(x){
     d <- guideSeqs.seqAgg[x,]
     m <- matrix(c(d$UnInf, d$Zika, d$totalUnInf-d$UnInf,  d$totalZika-d$Zika),  byrow = TRUE, nrow = 2)
     f <- fisher.test(m, alternative = 'less')
     bind_cols(d, data.frame(pval           = f$p.value,
                             oddsRatio      = f$estimate,
                             confInt95start = f$conf.int[1],
                             confInt95End   = f$conf.int[2])) }))

guideSeqs.seqAgg.tests <- arrange(guideSeqs.seqAgg.tests, confInt95End)
write.csv(guideSeqs.seqAgg.tests, file = 'data/guideSeqs.seqAgg.tests.csv')


guideSeqs.geneAgg.tests <- bind_rows(lapply(1:nrow(guideSeqs.geneAgg), function(x){
  d <- guideSeqs.geneAgg[x,]
  m <- matrix(c(d$UnInf, d$Zika, d$totalUnInf-d$UnInf,  d$totalZika-d$Zika),  byrow = TRUE, nrow = 2)
  f <- fisher.test(m, alternative = 'less')
  bind_cols(d, data.frame(pval           = f$p.value,
                          oddsRatio      = f$estimate,
                          confInt95start = f$conf.int[1],
                          confInt95End   = f$conf.int[2])) }))

guideSeqs.geneAgg.tests <- arrange(guideSeqs.geneAgg.tests, confInt95End)
write.csv(guideSeqs.geneAgg.tests, file = 'data/guideSeqs.geneAgg.tests.csv')


#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~


# Select data points to label.
i <- which(guideSeqs.seqAgg.tests$UnInfFreq > 0.0005)
guideSeqs.seqAgg.tests[i,]$label1 <- guideSeqs.seqAgg.tests[i,]$gene

guideSeqs.seqAgg.tests.plot1 <- ggplot(guideSeqs.seqAgg.tests, aes(UnInfFreq, ZikaFreq)) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  labs(x = 'Uninfected read frequency', y = 'Zika read frequency') +
  geom_text_repel(aes(label = label1)) 


guideSeqs.seqAgg.tests$label1 <- 'p > 1e-4'

i <- which(guideSeqs.seqAgg.tests$pval <= 1e-4)
guideSeqs.seqAgg.tests[i,]$label1 <- 'p < 1e-4'

i <- which(guideSeqs.seqAgg.tests$oddsRatio <= 1/2)
guideSeqs.seqAgg.tests[i,]$label1 <- 'oddsRatio < 1/2'

i <- which(guideSeqs.seqAgg.tests$confInt95End <= 1/2)
guideSeqs.seqAgg.tests[i,]$label1 <- 'oddsRatio 95 conf margin <= 1/2'

guideSeqs.seqAgg.tests <- bind_rows(subset(guideSeqs.seqAgg.tests, label1 == 'p > 1e-4'),
                                    subset(guideSeqs.seqAgg.tests, label1 == 'p < 1e-4'),
                                    subset(guideSeqs.seqAgg.tests, label1 == 'oddsRatio < 1/2'),
                                    subset(guideSeqs.seqAgg.tests, label1 == 'oddsRatio 95 conf margin <= 1/2'))

guideSeqs.seqAgg.tests$label1 <- factor(guideSeqs.seqAgg.tests$label1, levels = unique(guideSeqs.seqAgg.tests$label1))

guideSeqs.seqAgg.tests.plot2 <- ggplot(guideSeqs.seqAgg.tests, aes(UnInfFreq, ZikaFreq)) +
  theme_bw() +
  geom_point(aes(color=label1)) +
  scale_color_manual(name = 'Significance', values=c('gray70', 'dodgerblue', 'green3', 'red')) +
  labs(x = 'Uninfected read frequency', y = 'Zika read frequency') +
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 0.0002)) + ylim(c(0,0.0002))


write(subset(guideSeqs.seqAgg.tests, confInt95End <= 1/2)$gene, file='data/seqAgg.OR_95CI_UL_0.5.genes.txt')
write(subset(guideSeqs.seqAgg.tests, oddsRatio <= 1/4)$gene, file='data/seqAgg.OR_0.25.genes.txt')

#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~



# Select data points to label.
i <- which(guideSeqs.geneAgg.tests$gene %in% c('GAGE1', 'ZC3HAV1', 'FBXW2'))
guideSeqs.geneAgg.tests[i,]$label1 <- guideSeqs.geneAgg.tests[i,]$gene

guideSeqs.geneAgg.tests.plot1 <- ggplot(guideSeqs.geneAgg.tests, aes(UnInfFreq, ZikaFreq)) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  labs(x = 'Uninfected read frequency', y = 'Zika read frequency') +
  geom_text_repel(aes(label = label1)) 


guideSeqs.geneAgg.tests$label1 <- 'p > 1e-3'

i <- which(guideSeqs.geneAgg.tests$pval <= 1e-3)
guideSeqs.geneAgg.tests[i,]$label1 <- 'p < 1e-3'

i <- which(guideSeqs.geneAgg.tests$oddsRatio <= 1/2)
guideSeqs.geneAgg.tests[i,]$label1 <- 'oddsRatio < 1/2'

i <- which(guideSeqs.geneAgg.tests$confInt95End <= 1/2)
guideSeqs.geneAgg.tests[i,]$label1 <- 'oddsRatio 95 conf margin <= 1/2'

guideSeqs.geneAgg.tests <- bind_rows(subset(guideSeqs.geneAgg.tests, label1 == 'p > 1e-3'),
                                     subset(guideSeqs.geneAgg.tests, label1 == 'p < 1e-3'),
                                     subset(guideSeqs.geneAgg.tests, label1 == 'oddsRatio < 1/2'),
                                     subset(guideSeqs.geneAgg.tests, label1 == 'oddsRatio 95 conf margin <= 1/2'))

guideSeqs.geneAgg.tests$label1 <- factor(guideSeqs.geneAgg.tests$label1, levels = unique(guideSeqs.geneAgg.tests$label1))

guideSeqs.geneAgg.tests.plot2 <- ggplot(guideSeqs.geneAgg.tests, aes(UnInfFreq, ZikaFreq)) +
  theme_bw() +
  geom_point(aes(color=label1)) +
  scale_color_manual(name = 'Significance', values=c('gray70', 'dodgerblue', 'green3', 'red')) +
  labs(x = 'Uninfected read frequency', y = 'Zika read frequency') +
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 0.0005)) + ylim(c(0,0.0005))



# Create a save point
save.image(file='savePoints/SP3.RData')


# Save the data needed to generate the report to a separate file.
save(stats, report, guideSeqLengthPlot1, guideSeqLengthPlot2, rarificationPlot, 
     expectedGuides, guideSeqs.seqAgg.tests, guideSeqs.seqAgg.tests.plot1, 
     guideSeqs.seqAgg.tests.plot2, guideSeqs.geneAgg.tests.plot1, guideSeqs.geneAgg.tests.plot2, 
     guideSeqs, guideSeqs.geneAgg.tests, file='savePoints/project.RData')


guideSeqs.seqAgg.tests$label1  <- NULL
guideSeqs.seqAgg.tests$label2  <- NULL
guideSeqs.geneAgg.tests$label1 <- NULL
guideSeqs.geneAgg.tests$label2 <- NULL

write.csv(guideSeqs.seqAgg.tests, file='data/guideSeqs.seqAgg.tests.csv')
write.csv(guideSeqs.geneAgg.tests, file='data/guideSeqs.geneAgg.tests.csv')


