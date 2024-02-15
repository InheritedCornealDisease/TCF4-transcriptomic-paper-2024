library(maser)
library(rtracklayer)
library(httpgd)

# path to data
path <- "../output/fastq_start/rmats/run1_CvE"

# create maser object
events <- maser(path, c("Control", "Expansion+"), ftype = "JC")  # change according to the disease conditions you want to look at (Control, Exp+ or Exp-)

events

head(summary(events, type = "SE")[, 1:9]) # look at Skipped Exons, for example

# subset significant events by fdr and dpsi filtering
sig_events <- topEvents(events, fdr = 0.05, deltaPSI = 0.1)

sig_events # see how many total significant events for this particular comparison, and a break down by splice types.
