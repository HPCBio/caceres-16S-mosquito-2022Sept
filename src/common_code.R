# Common code 
# This should be moved into the HPCBio microbiome folder

# Note that not all libraries will be needed.  Most phyloseq code uses ggplot and tidyverse internally, therefore we explicitly load here
library(tidyverse)
library(phyloseq)
library(plotly)
library(scales)
library(knitr)
library(gridExtra)

# this seems to have issues with caching and phyloseq
# library(ggtree) 

# For normalization
# library(metagenomeSeq)

# phylogenetic tree input
library(ape)

# read/modify BIOM 
# library(biomformat)

# ggplot functions for trees and dendrograms
# library(ggdendro)

# distance measures, PERMANOVA, ANOSIM
library(vegan)

# generation of stats values for graphs
library(ggpubr)

# normalization (CLR)
library(mixOmics)

# to get labels2color
library(WGCNA)

# mixed models (needs to be updated)
# library(lme4)
# library(lmerTest)
# library(nlme)
# to get post-hoc tests for mixed-model tests 
# library(lsmeans)

# sample decontamination 
# library(decontam)

# library(devtools)

# needed in case we want to use ANCOM
#library(exactRankTests)

#Other libraries I added later
##library(BiocManager)
##BiocManager::install("microbiome")
##library(devtools)
##devtools::install_github("gauravsk/ranacapa")
##devtools::install_github("hpcbio/plotly_microbiome")
# library(plotly.microbiome)
# library(microbiome)
library(ranacapa)

# this is to load some extension helper code, see: https://github.com/HPCBio/phyloseq-extended
devtools::load_all('~/src/phyloseq-extended/')

# Remove the tags on the taxonomic ranks, which are redundant with the column headers.
stripTaxaTags <- function(physeq) {
  oldMA <- as(tax_table(physeq), "matrix")
  newMA <- apply(oldMA, 2, function(x) {sub('\\w__','', x)})
  if (inherits(physeq, "taxonomyTable")) {
    return(tax_table(newMA))
  }
  else {
    tax_table(physeq) <- tax_table(newMA)
    return(physeq)
  }
}

# Convert sequences to names (culled from https://github.com/LangilleLab/microbiome_helper/blob/master/convert_dada2_out.R) 

renameTaxIds <- function(physeq, file.name="seqs.fasta") {
  suppressMessages(require("ShortRead"))
  seqtab.physeq <- otu_table(physeq)
  seqs <- colnames(seqtab.physeq)
  ids_study <- paste("seq", 1:ncol(seqtab.physeq), sep = "_")
  seqs.dna <- ShortRead(sread = DNAStringSet(seqs), id = BStringSet(ids_study))
  # Write out fasta file.
  writeFasta(seqs.dna, file = file.name)
  taxa_names(physeq) <- ids_study
  # TODO: add the sequences back to the phyloseq instance
  # physeq <- merge_phyloseq(physeq)
  return(physeq)
}

# original code: https://github.com/twbattaglia/btools/blob/master/R/estimate_pd.R
estimate_pd <- function(phylo) {
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  if (!require('picante')) stop("Function requires the picante library.")
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from phyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  #
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

# CLR normalization 
# (from McMurdie (Meth Mol Bio 2018) supplemental package)
zero_comp = function(x){
  if(taxa_are_rows(x)){x <- t(x)}
  matx = otu_table(x)
  # `zCompositions::cmultRepl` expects the samples to be in rows and OTUs to be in columns
  matxzc = zCompositions::cmultRepl(matx, method="CZM", output="p-counts")
  otu_table(x) <- otu_table(matxzc, taxa_are_rows = FALSE)
  return(x)
}
# CLR definition
geometric_mean = function(x){
  exp(mean(log(x)))
}
clr = function(x, base=2){
  x <- log((x / geometric_mean(x)), base)
}
phyloseq_CLR = function(physeq){
  suppressMessages({physeq <- zero_comp(physeq)})
  return(transform_sample_counts(physeq, fun = clr))
}

# this is a modification of the plot_richness function from phyloseq, but takes as input a pre-generated matrix of estimates from `estimate_richness` or any other function, plus the phyloseq instance.  
plot_richness_estimates = function(physeq, 
                                   erDF, 
                                   x="samples", 
                                   color=NULL, 
                                   shape=NULL, 
                                   title=NULL,
                                   scales="free_y", 
                                   nrow=1, 
                                   sortby=NULL) {
  # TODO: add sanity check on matrix (e.g. rows == sample IDs, sample names, and column names)
  
  # Measures may have been renamed in `erDF`. Replace it with the name from erDF
  measures = colnames(erDF)
  # Define "measure" variables and s.e. labels, for melting.
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  # Remove any S.E. from `measures`
  measures = measures[!measures %in% ses]
  # Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if( !is.null(sample_data(physeq, errorIfNULL=FALSE)) ){
    # Include the sample data, if it is there.
    DF <- data.frame(erDF, sample_data(physeq))
  } else {
    # If no sample data, leave it out.
    DF <- data.frame(erDF)
  }
  if( !"samples" %in% colnames(DF) ){
    # If there is no "samples" variable in DF, add it
    DF$samples <- sample_names(physeq)
  }
  # sample_names used to be default, and should also work.
  # #backwardcompatibility
  if( !is.null(x) ){
    if( x %in% c("sample", "samples", "sample_names", "sample.names") ){
      x <- "samples"
    }
  } else {
    # If x was NULL for some reason, set it to "samples"
    x <- "samples"
  }
  # melt to display different alpha-measures separately
  mdf = reshape2::melt(DF, measure.vars=measures)
  # Initialize the se column. Helpful even if not used.
  mdf$se <- NA_integer_
  if( length(ses) > 0 ){
    ## Merge s.e. into one "se" column
    # Define conversion vector, `selabs`
    selabs = ses
    # Trim the "se." from the names
    names(selabs) <- substr(selabs, 4, 100)
    # Make first letter of selabs' names uppercase
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 1, 1))
    # use selabs conversion vector to process `mdf`
    mdf$wse <- sapply(as.character(mdf$variable), function(i, selabs){selabs[i]}, selabs)
    for( i in 1:nrow(mdf) ){
      if( !is.na(mdf[i, "wse"]) ){
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    # prune the redundant columns
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  ## Interpret measures
  # If not provided (default), keep all 
  if( !is.null(measures) ){
    if( any(measures %in% as.character(mdf$variable)) ){
      # If any measures were in mdf, then subset to just those.
      mdf <- mdf[as.character(mdf$variable) %in% measures, ]
    } else {
      # Else, print warning about bad option choice for measures, keeping all.
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  # Address `sortby` argument
  if(!is.null(sortby)){
    if(!all(sortby %in% levels(mdf$variable))){
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if(!is.discrete(mdf[, x])){
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if(all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, x])){
      # Replace x-factor with same factor that has levels re-ordered according to `sortby`
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x],
                         levels = names(sort(tapply(X = mdf[wh.sortby, "value"],
                                                    INDEX = mdf[wh.sortby, x],
                                                    mean,
                                                    na.rm=TRUE, simplify = TRUE))))
    }
  }
  # Define variable mapping
  richness_map = aes_string(x=x, y="value", colour=color, shape=shape)
  # Make the ggplot.
  p = ggplot(mdf, richness_map) + geom_point(na.rm=TRUE)  
  # Add error bars if mdf$se is not all NA
  if( any(!is.na(mdf[, "se"])) ){
    p = p + geom_errorbar(aes(ymax=value + se, ymin=value - se), width=0.1) 
  }
  # Rotate horizontal axis labels, and adjust
  p = p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0))
  # Add y-label 
  p = p + ylab('Alpha Diversity Measure') 
  # Facet wrap using user-options
  p = p + facet_wrap(~variable, nrow=nrow, scales=scales)
  # Optionally add a title to the plot
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }
  return(p)
}

options(stringsAsFactors = FALSE)
theme_set(theme_bw())