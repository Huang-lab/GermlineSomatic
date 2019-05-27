##### somatic_germline_overlap.R #####
# Kuan-lin Huang @ WashU 2018
# Find overlap of genes/variants for somatic/germline variants

### dependencies ###
# setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/germline_somatic_analysis/somatic_germline_overlap")
source("../global_aes_out.R")
source("../dependency_files.R")
source("../load_somatic.R")
library(eulerr)

# check if there is statistical enrichment of overlaps
# $ awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic_merged
# 49586385
exonSize = 49586385
nPath = nrow(pathVarP)

### somatic mutation
# $ gzcat /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz | cut -f1,4 | sort | uniq -c | awk '$1 >2 && $3 != "."' | wc -l
# 68537
numSomaticOverlap = nrow(pathVarP[pathVarP$SomaticMutation_peptide_location_count != 0,])
somaticMutRate = 68537/exonSize  
poisson.test(numSomaticOverlap, T = nPath, r = somaticMutRate, conf.level = 0.95, alternative = "greater")

numSomaticOverlapRec = nrow(pathVarP[pathVarP$SomaticMutation_peptide_location_count > 2,])

# plot
pathVarP_hot = pathVarP[pathVarP$SomaticMutation_peptide_location_count != 0,]
genes = names(table(pathVarP_hot$HUGO_Symbol)[table(pathVarP_hot$HUGO_Symbol) > 2])
pathVarP_hot$somatic_count_plot = pathVarP_hot$SomaticMutation_peptide_location_count
pathVarP_hot$HGVSp_short_plot = gsub("p.","",pathVarP_hot$HGVSp_short)
germline_recurrence = data.frame(table(pathVarP_hot$HGVSp_short))
colnames(germline_recurrence) = c("HGVSp_short","GermlineRecurrence")
pathVarP_hot_m = merge(pathVarP_hot,germline_recurrence, by = "HGVSp_short")
#pathVarP_hot$somatic_count_plot[pathVarP_hot$somatic_count_plot> 100 ]  = 100
pathVarP_hot_m = pathVarP_hot_m[!duplicated(pathVarP_hot_m$HGVSp_short),]

p = ggplot(pathVarP_hot_m[pathVarP_hot_m$HUGO_Symbol %in% genes,],aes(y=HUGO_Symbol, x =somatic_count_plot))
#p = p + facet_grid(PCGP~Gene_Classification,drop=T,scale="free",space="free")
p = p + facet_grid(Gene_Classification~ .,drop=T,scale="free_y",space="free_y")
p = p + geom_point(alpha = 0.5, stroke=0,aes(size = GermlineRecurrence)) + theme_bw()  #+ guides(color=FALSE)
#p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) #+ geom_density2d(alpha=0.5)
p = p + geom_text_repel(aes(label=ifelse(duplicated(HGVSp_short) | somatic_count_plot<3,NA,HGVSp_short_plot)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + scale_x_log10()
p = p + expand_limits(x = 0)
#p = p + coord_equal() + getLOHColorScale()
p = p + labs(x = "Recurrence of co−localizing somatic mutations", y = "Gene")
p
fn = "out/pathVarP_somatic_potlight.pdf"
ggsave(file=fn, width=7, h =7, useDingbats=FALSE)


# plotting intersects
fn = "out/germline_somatic_gene_overlap.pdf"
pdf(fn)
fit = euler(c("Germline" = length(unique(pathVarP$HUGO_Symbol)) - length(intersect(unique(pathVarP$HUGO_Symbol),somaticDriver299)), "SomaticDriver" =length(somaticDriver299) - length(intersect(unique(pathVarP$HUGO_Symbol),somaticDriver299)), "Germline&SomaticDriver" = length(intersect(unique(pathVarP$HUGO_Symbol),somaticDriver299))))
plot(fit, quantities = TRUE,fills = list(fill = c("red", "steelblue4"), alpha = 0.5))
dev.off()

fn = "out/germline_somatic_variant_overlap.pdf"
pdf(fn)
fit = euler(c("Germline" = nPath - numSomaticOverlap, "Somatic" =68537 - numSomaticOverlap, "Germline&Somatic" = numSomaticOverlap))
plot(fit, quantities = TRUE,fills = list(fill = c("red", "steelblue4"), alpha = 0.5))
dev.off()
