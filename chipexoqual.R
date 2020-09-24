# Read samples from input.
library('optparse')
option_list = list(
  make_option(c('-p', '--prefix'), type='character', default='',
              help='file name prefix', metavar='character'),
  make_option(c('-s', '--sampling'), type='integer', default=3000000,
              help='number of reads to keep in sampling', metavar='number'),
  make_option(c('-t', '--threads'), type='integer', default=1,
              help='number of threads', metavar='number')
);
opt_parser = OptionParser(option_list=option_list);
arguments = parse_args(opt_parser, positional_arguments = TRUE);
opt = arguments$options
bams = arguments$args
samples = lapply(bams, tools::file_path_sans_ext)

# Chip-exo quality
cat('Loading library ChIPexoQual\n')
suppressPackageStartupMessages(library('ChIPexoQual'))

cat('Reading bam files:', bams, '\n')
bam_reads = lapply(bams, readGAlignments)
cat('Creating ExoData objects\n')
computeExoData <- function() {
  for (i in 1:30) {
    cat('  computing exodata, iteration', i, '\n')
    rand_bam_reads = lapply(bam_reads, function(x) sample(x, min(length(x), opt$sampling)))
    try(return(lapply(rand_bam_reads, function(x) ExoData(reads=x, mc.cores = opt$threads, verbose = FALSE))))
  }
}
exoData <- computeExoData()
if (is.na(exoData)) stop("Could not compute exo data")

library_complexity_output = paste(opt$prefix, 'library_complexity.pdf', sep='')
cat('Creating library complexity plots in', library_complexity_output, '\n')
pdf(library_complexity_output)
ARCvURCplot(exoData, names.input = samples)

unique_library_complexity_output = paste(opt$prefix, 'unique_library_complexity.pdf', sep='')
cat('Creating unique library complexity plots in', unique_library_complexity_output, '\n')
pdf(unique_library_complexity_output)
ARCvURCplot(lapply(exoData, subset, uniquePos > 10), names.input = samples)

strand_balance_output = paste(opt$prefix, 'strand_balance.pdf', sep='')
cat('Creating strand balance plots in', strand_balance_output, '\n')
pdf(strand_balance_output)
p1 = regionCompplot(exoData, names.input = samples, depth.values = seq_len(50))
p2 = FSRDistplot(exoData, names.input = samples, quantiles = c(.25,.5,.75), depth.values = seq_len(100))
gridExtra::grid.arrange(p1, p2, nrow = 1)

quality_evaluation_output = paste(opt$prefix, 'quality_evaluation.pdf', sep='')
cat('Creating quality evaluation plots in', quality_evaluation_output, '\n')
pdf(quality_evaluation_output)
p1 = paramDistBoxplot(exoData, which.param = "beta1", names.input = samples)
p2 = paramDistBoxplot(exoData, which.param = "beta2", names.input = samples)
gridExtra::grid.arrange(p1, p2, nrow = 1)

while (!is.null(dev.list()))  dev.off()
