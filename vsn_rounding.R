# test case for rounding issue
# Aurelien Ginolhac LSRU
suppressMessages(library("reshape2"))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("argparse"))
# from bioconductor
suppressMessages(library("vsn"))
suppressMessages(library("limma"))
suppressMessages(library("qvalue"))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-n", "--normalisation", action="store_true", default=FALSE,
                    help="perform VSN normalisation [%(default)s]")
parser$add_argument("-c", "--count", type="integer", default=900, 
                    help="Number of proteins to use [default %(default)s]",
                    metavar="number")

args <- parser$parse_args()

# FUNCTIONS

n_replicate2 <- function(df, cond){
  table(df) %>% 
    as.data.frame %>%
    filter(df == cond) %>% 
    use_series(Freq) %>%
    return()
}
# code from Enrico Glaab
contrasts_limma <- function(dat, groups, ID, control, treatment){
  
  design <- model.matrix(~ -1 + 
                           factor(letters[ifelse(groups[c(which(groups == control),
                                                          which(groups == treatment))] == control, 0, 1) + 1]))
  
  outnum <- as.vector(ifelse(groups[c(which(groups == control),
                                      which(groups == treatment))] == control,
                             0, 1))
  colnames(design) <- unique(letters[outnum + 1])
  
  corfit <- duplicateCorrelation(dat[,c(which(groups == control),
                                        which(groups == treatment))],
                                 design,
                                 block = ID)
  #corfit$consensus
  fit <- lmFit(dat[,c(which(groups == control),
                      which(groups == treatment))],
               design,
               block = ID,
               correlation = corfit$consensus)
  
  cm <- makeContrasts(comp = a - b, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  return(eBayes(fit2))
}

# MAIN
set.seed(1234)

mass.cast <- readRDS("24h_test.rds")

# convert to matrix
mass.mat <- mass.cast %>%
  select(-c(Accession)) %>%
  as.matrix  
if (args$normalisation){
  # do variance stabilization 
  vsndat <- mass.mat %>% 
    vsn2 %>% 
    exprs
} else {
  vsndat <- mass.mat
}
cat("matrix:\n")
print(head(vsndat))

groups <- lapply(strsplit(colnames(vsndat), "_"), function(x) x[1]) %>% as.character
sampleID <- c(1:n_replicate2(groups, "control"),
              1:n_replicate2(groups, "IL27calprotectin"))
# constrat using limma
fit <- contrasts_limma(vsndat, groups, sampleID, "IL27calprotectin", "control")
res <- topTable(fit, n = nrow(fit)) %>%
  mutate(qvalue = qvalue(P.Value)$qvalues)

cat("top 10\n")
print(head(res, 10))
cat(paste(res %>% filter(qvalue <= 0.1) %>% nrow, "proteins qvalue <= 0.1\n"))
