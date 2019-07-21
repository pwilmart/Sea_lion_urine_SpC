
library("tidyverse")
library("edgeR")
library("limma")
library("psych")
library("gridExtra")
library("scales")

# read in the prepped data
paw_spc <- read_tsv("QC_check_669x19.txt")

# save accessions in vector and remove from data table
accession <- paw_spc$Accession
paw_spc <- select(paw_spc, -Accession)

# replace zeros with 0.25
paw_spc[paw_spc == 0] <- 0.25
head(paw_spc)
nrow(paw_spc)

# function for simple normalization
SL_Norm <- function(df, color = NULL, plot = TRUE) {
    # This makes each channel sum to the average grand total
        # df - data frame of TMT intensities
        # returns a new data frame with normalized values
    
    # compute scaling factors to make colsums match the average sum
    norm_facs <- mean(c(colSums(df))) / colSums(df)
    cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
    df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(df_sl), col = color, notch = TRUE, main = "SL Normalized data")
    }
    df_sl
}

# normalize the data before plotting
color <- c(rep('blue', 8), rep('red', 11))
paw_sl <- SL_Norm(paw_spc, color)

# load the data into edgeR data structures
# group labels need to be factors
group <- c(rep("Control", 8), rep("lepto", 11))
y <- DGEList(counts = paw_spc, group = group, genes = accession)

# run the TMM normalization (and library size corrections)
y <- calcNormFactors(y)

# check normalizations
y$samples

apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
    # computes the tmm normalized data from the DGEList object
        # y - DGEList object
        # returns a dataframe with normalized intensities
    
    # compute grand total (library size) scalings
    lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size

    # the TMM factors are library adjustment factors (so divide by them)
    norm_facs <- lib_facs / y$samples$norm.factors
    cat("Overall Factors (lib.size+TMM):\n", sprintf("%-5s -> %f\n", 
                                                     colnames(y$counts), norm_facs))

    # compute the normalized data as a new data frame
    df_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
    colnames(df_tmm) <- str_c(colnames(y$counts), "_tmm")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(df_tmm), col = color, notch = TRUE, main = "TMM Normalized data")
    }
    df_tmm
}

# get normalized data and check boxplots
paw_spc_tmm <- apply_tmm_factors(y, color)

# check clustering
plotMDS(y, col = color, main = "Control vs Lepto")

CV <- function(df) {
    # Computes CVs of data frame rows
        # df - data frame, 
        # returns vector of CVs (%)
    ave <- rowMeans(df)    # compute averages
    sd <- apply(df, 1, sd) # compute standard deviations
    cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}

labeled_boxplot <- function(df, ylim, title) {
    # Makes a box plot with the median value labeled
        # df - data frame with data to compute CVs of
        # ylim - upper limit for y-axis
        # title - plot title
    cv = CV(df)
    boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
    text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
         labels = round(boxplot.stats(cv)$stats[3], 1))
}

# define the columns for each condition
C <- 1:8
L <- 9:19
CD <- c(1, 2, 5, 6, 7, 8)
LD <- c(9, 11, 12, 13, 15, 16, 17, 19)


# see what effect TMM had on CV distributions
par(mfrow = c(2, 2))
labeled_boxplot(paw_spc[C], 200, "All Control before")
labeled_boxplot(paw_spc[L], 200, "All Lepto before")
labeled_boxplot(paw_spc_tmm[C], 200, "All Control after")
labeled_boxplot(paw_spc_tmm[L], 200, "All Lepto after")

# dropping 2 controls and 3 lepto samples
par(mfrow = c(2, 2))
labeled_boxplot(paw_spc[CD], 200, "Control before")
labeled_boxplot(paw_spc[LD], 200, "Lepto before")
labeled_boxplot(paw_spc_tmm[CD], 200, "Control after")
labeled_boxplot(paw_spc_tmm[LD], 200, "Lepto after")
par(mfrow = c(1, 1))

# load the data into edgeR data structures
# group labels need to be factors
countsd <- paw_spc[c(C, LD)]
colord <- c(rep('blue', 8), rep('red', 8))
groupd <- c(rep("Control", 8), rep("lepto", 8))
yd <- DGEList(counts = countsd, group = groupd, genes = accession)

# run the TMM normalization (and library size corrections)
yd <- calcNormFactors(yd)

# check normalizations
yd$samples

# get normalized data and make boxplots
paw_spc_tmm <- apply_tmm_factors(yd, colord)

# check clustering
plotMDS(yd, col = colord, main = "ControlD vs LeptoD")

# log the session
sessionInfo()


