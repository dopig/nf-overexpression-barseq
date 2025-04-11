#!/usr/bin/env Rscript

# This is an adaptation of fit2.R by Morgan Price
# This link should download the original R script directly
# https://figshare.com/ndownloader/files/42456840
# Or you can find details/its directory here:
# https://figshare.com/articles/dataset/Barcoded_overexpression_screens_in_gut_Bacteroidales_identify_genes_with_new_roles_in_carbon_utilization_and_stress_resistance_/24195054
#
# Modifications included:
# - Swapping out hardcoded configurations for values in a JSON file
# - (Temporarily) blocking use of protein similarity, with intention of allowing a flag for it later.
#   The reasoning here is that it won't be useful when working with simulated data.
# - Removing checks for presence in E. coli donor strain. Longer term, this should be a flag like above.
#
# Here's an example of how/what to run:
#   Rscript src/testR.R --json_path nextflow/lib.json
#     --samples_tsv data/reference/bobaseq_barseq_mgenitalium_samples.tsv
#     --codes_dir data/output/output-2025-03-21-230345/barseq
#     --r_image new.rimage
#     --output_tsv new.tsv
#
# The json/csv file should have the following columns:
#   - lib, the library name (e.g. Mgenitalium)
#   - path, the path to the mapping output directory (typically the path to get to '05-BC_and_genes_dfs')
#   - feature_table_path, path to the feature table (e.g. data/reference/GCF_000027325.1_ASM2732v1_feature_table.txt)

library(optparse)
# library(jsonlite)

skip_candidate_similarity = TRUE

# Define command-line arguments
option_list <- list(
  make_option("--bobaseq_path", type="character", default=file.path("/app/bobaseq.R"), help="Path to the script bobaseq.R (default: bobaseq.R)"),
  # make_option("--json_path", type="character", help="Path to library JSON file (alternatively, use --csv_path)"),
  # make_option("--csv_path", type="character", help="Path to library CSV file (alternatively, use --json_path)"),
  # make_option("--lib_name", type="character", help="Library name"),
  make_option("--mapping_dir", type="character", help="Path to mapping directory (generally has same name as library"),
  make_option("--feature_table_path", type="character", help="Path to feature table"),
  make_option("--samples_tsv", type="character", help="Path to samples TSV file"),
  make_option("--bc2best_pos", type="character", help="List of paths to BC2best_pos files"),
  make_option("--bc_gene_mapped", type="character", help="List of paths to BC_gene_mapped files"),
  make_option("--r_image", type="character", default=file.path("fitness.Rimage"), help="Path to R image output"),
  make_option("--output_tsv", type="character", default=file.path("fitness.tsv"), help="Path to output top proteins TSV file")
  # make_option("--output_pdf", type="character", help="Path to generate a PDF output; if none given, this step will be skipped")
)
# Parse arguments
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)

# if (is.null(args$csv_path) && !is.null(args$json_path)) {
#   libs <- fromJSON(args$json_path)
# } else if (!is.null(args$csv_path) && is.null(args$json_path)) {
#   libs <- read.csv(args$csv_path)
# } else {
#   stop("Error: You must provide either --json_path or --csv_path, but not both.")
# }

libs <- data.frame(
  lib = args$mapping_dir,
  dir = file.path(args$mapping_dir, '05-BC_and_genes_dfs'),
  feature_table_path = args$feature_table_path
)

# Load bobaseq.R
source(args$bobaseq_path);

#Proteins and similarity

## Set up libaries and barcode mappings
maps = lapply(libs$dir, function(dir) read.delim(paste0(dir, "/BC2best_pos.tsv"), as.is=T));

names(maps) = libs$lib;
for(l in names(maps)) maps[[l]]$lib = l;
mapsAll = do.call(rbind, maps);
names(mapsAll)[1] = "barcode";
names(mapsAll)[names(mapsAll)=="tstart"] = "start";
names(mapsAll)[names(mapsAll)=="tend"] = "end";
cat("Barcodes with at least 1 mapping: ", length(unique(mapsAll$barcode)), "\n");
maps = subset(mapsAll, is.unique(barcode));
cat("Uniquely-mapped barcodes: ", nrow(maps), "\n");

## My plans don't involve a conjugation donor
# DONOR_STRAIN_NAME is the strain name of the conjugation donor
# DONOR_BC_PATH is a tsv file, for mapping barcode count files to sources, it should have the following:
#     barseq.file - path to a barseq count file
#     strain - string that matches DONOR_STRAIN_NAME above for relevant rows
#     sample - sample name (can just be Rep1, Rep2, ...)
# # Subset the maps to barcodes seen in the conjugation donor
# divs = read.delim(DONOR_BC_PATH,as.is=T);
# divsWM = subset(divs, grepl(DONOR_STRAIN_NAME, strain));
# bcWM = lapply(divsWM$barseq.file, read.delim, as.is=T);
# names(bcWM) = divsWM$sample;
# for (n in names(bcWM)) bcWM[[n]]$sample = n;
# bcWM = do.call(rbind, bcWM);
# bcWM = merge(bcWM, divsWM[,c("lib","sample")]);
# cat("Barcodes seen in ", DONOR_STRAIN_NAME, ": ", nrow(bcWM), "\n");
# maps = subset(maps,  barcode %in% bcWM$barcode);
# cat("Mapped barcodes seen in ", DONOR_STRAIN_NAME, ": ", nrow(maps), "\n");

## Set up proteins and similarity
named_ft_paths <- setNames(as.list(libs$feature_table_path), libs$lib)
ft = lapply(named_ft_paths, function(nftp) read.delim(nftp, as.is=T));
for (n in names(ft)) ft[[n]]$lib = n;
ft = do.call(rbind, ft);
prot = subset(ft, X..feature=="CDS" & class=="with_protein")[,c("lib","genomic_accession","start","end","strand","locus_tag","product_accession","name")];
names(prot)[names(prot)=="genomic_accession"] = "contig";
names(prot)[names(prot)=="product_accession"] = "protId";
names(prot)[names(prot)=="name"] = "desc";
cat("Loaded " , nrow(prot), " proteins\n");

if (!skip_candidate_similarity){
  psim = proteinSimilarity(prot,
    paste0(details$paths$assembly_dir, "/faa.hits"),
    paste0(details$paths$assembly_dir, "/faa.len"));
  cat("Loaded " , nrow(psim), " pairs of similar proteins (counted both ways)\n");
}
## Gene mappings and barcodes to proteins
g = lapply(libs$dir, function(dir) read.delim(paste0(dir, "/BC_gene_mapped.tsv"), as.is=T));
names(g) = libs$lib;
for(l in names(g)) g[[l]]$lib = l;
g = do.call(rbind, g);
names(g)[1] = "barcode";
g = subset(g, barcode %in% maps$barcode);
glist = aggregate(g[,"locus_tag",drop=F],  g[,"barcode",drop=F], function(x) paste(sort(unique(x)), collapse=","));
# Add (comma-delimited) list of locus tags and #genes covered to maps
maps = merge(maps, glist);
names(maps)[names(maps)=="locus_tag"] = "locus_tags";
maps = merge(maps, unique(g[,c("barcode","gene_count")]));
p = merge(maps[,words("barcode start end contig lib gene_count")], unique(g[,c("barcode","locus_tag","directionality")]));
# note p$strand is for the locus
p = merge(p, prot, by=c("contig","lib","locus_tag"), suffixes=c(".insert",""));

## Set up samples
samples = read.delim(args$samples_tsv, as.is=T);
samples$name = paste0(samples$set, samples$index);
# samples$file = findCodesFiles(samples, baseDir=args$codes_dir);
samples$file <- paste0(samples$sampleId, "_", samples$index, ".codes");

## Set up counts per sample
counts = getCounts(samples, maps);
samples = cbind(samples, totalCounts(counts, samples, maps));
samples$used = samples$totLib >= 10000 & samples$fLib >= 0.9 &
  ifelse(samples$Group=="Time0", samples$minStrainFrac > 0.05, TRUE);

## Log ratios and high barcodes
countTime0 = time0Counts(samples, counts, maps);
lr = logRatios(samples, counts, countTime0, maps);
hi = highBarcodes(samples, counts, countTime0, lr, maps);#, minLogRatio=2, minZ=2);
hi$confirmedByOverlap = confirmedByOverlap(hi, hi$name);

## Replicates (within the same day)
repSamples = merge(
  samples[samples$used & samples$Group != "Time0", c("name","background","lib","t0set","Group","desc")],
  samples[samples$used & samples$Group != "Time0", c("name","background","lib","t0set","Group","desc")],
  by=c("background","lib","t0set","Group","desc"), suffixes=1:2);
repSamples = subset(repSamples, name1 < name2);
names(repSamples)[names(repSamples)=="lib"] = "libSpec";
repSamples = repSamples[order(repSamples$background, repSamples$Group, repSamples$desc, repSamples$libSpec, repSamples$name1, repSamples$name2),];
repSamples = subset(repSamples, name1 %in% names(lr) & name2 %in% names(lr));

safe_analysis <- function() {
  print("Trying the part that breaks now!!!")
  hiRep <- highBarcodesReplicates(samples, lr, hi)
  topRegion <- topHitPerRegion(hiRep, paste(hiRep$t0set, hiRep$background, hiRep$Group, hiRep$desc), hiRep$avg)
  topCand <- regionsToCandidates(topRegion, samples, lr, p)
  if (!skip_candidate_similarity) topCandSim <- candidateSimilarity(topCand, psim)
  topCand <- topCandAddInserts(topCand, hi)
  if (!skip_candidate_similarity) topCand <- topCandAddSim(topCand, topCandSim)
  if (skip_candidate_similarity) topCand$nSim <- 0
  topProt <- topProtPerRegion(topCand, prot)
  writeDelim(topProt[, words("background Group expDesc t0set locus_tag nHi avgHi confirmedByOverlap nRep avgProt protId protDesc")], args$output_tsv)
  cat("Wrote", args$output_tsv, "\n")
}

result <- tryCatch(safe_analysis(), error = function(e) {
  message("Warning: analysis failed — likely due to no high barcode replicates.");
  file.create(args$output_tsv);
})

save.image(args$r_image);
cat("Wrote", args$r_image, "\n");

# ## Plots
# if (!is.null(args$output_pdf)) {
#   pdf(args$output_pdf, width=5, height=5, pointsize=10);
#   par(mgp=2:0, mar=c(3,3,3,1));
#   mapply(function(n1,n2,label) cmpPlot(n1, n2, lr, hi, main=label),
#     repSamples$name1, repSamples$name2,
#     paste(repSamples$background, repSamples$Group, repSamples$libSpec, "\n", repSamples$desc));
#   dev.off();

#   cat("Wrote", args$output_pdf, "\n");
# }
