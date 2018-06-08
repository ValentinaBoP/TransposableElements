# libraries required
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# read and subset the main input file
DATA = fread(file = "tuatara_30Sep2015_rUdWx.fasta.align.k2p.noCpG.size", header = F, sep = "\t", stringsAsFactors = F)
TUATARA = TUATARA[,c(5,6,7,8,10,11,17,16)]
names(TUATARA) = c("Scaffold", "Begin", "End", "Length", "Strand", "Element", "Divergence", "ID")

# delete elements with divergence greater than 1 
TUATARA = TUATARA[TUATARA$Divergence < 1,]

# DNA TEs
DNA = TUATARA[grepl(pattern = 'DNA', x = TUATARA$Element),]
# get the names of the subfamilies - the string before the '#'
DNA$Subfamily = sapply(strsplit(as.character(DNA$Element), "#"), "[[", 1)
# get the names of the superfamily/family
DNA$Superfamily = sapply(strsplit(as.character(DNA$Element), "#"), "[[", 2)
# round the divergence values
DNA$Divergence = DNA$Divergence * 100
DNA$RoundDiv = floor(DNA$Divergence)
# create a factor with the name of the subfamily and the divergence associated to it
# so we can get the number of bps associated to that particular subfamily at that particular divergence
DNA$Factor = paste(DNA$Subfamily, DNA$RoundDiv, sep = "$")

# general landscape - bps occupied
DNA_bps = aggregate(Length ~ Factor, DNA, sum)
DNA_bps$Subfamily = sapply(strsplit(DNA_bps$Factor, "\\$"), "[[", 1)
DNA_bps$Divergence = sapply(strsplit(DNA_bps$Factor, "\\$"), "[[", 2)

# conversion in megabases
DNA_bps$Mb = DNA_bps$Length / 1000000
DNA_bps = DNA_bps[DNA_bps$Mb > 0.05,]

# here I aggregate some of the subfamilies to obtain a more tidy plot
PATTERN = c("Charlie1-L_tua", "Chompy-2a_tua", "DNA-1_Gav", "Harbinger-2-L_tua", "Harbinger-3-L_tua", "Harbinger-7_CPB", "hAT-11_AMi", "hAT-12_Crp", "hAT-12B_Crp", "hAT-1598-L_tua", "hAT-17_Croc", "hAT-17B_Croc", "hAT-2_Crp", "hAT-29-L_tua", "hAT-39_tua", "hAT-4_AMi", "hAT-4B_AMi", "hAT-6_AMi", "hAT-66-L-tua", "hAT-7-L-tua", "hAT-N11-L_tua", "hAT-N18-Lb_tua", "hAT-N18_tua", "hAT6-N1_Croc", "hATN-3-La_tua", "hATN-3-Lb_tua", "hATN-3-Lc_tua", "hATN-3-Ld_tua", "MER45B-L_tua", "MER45R-La_tua", "MER45R-Lb_tua", "MER45R-Lc_tua", "nhAT-3-La_tua", "piggyBac-1-Lb_tua", "Polinton-1_AMi", "Polinton-1_tua.inc", "Polinton-2_AMi", "REP131", "tuaDNA10", "tuaDNA11", "tuaDNA12", "tuaDNA13", "tuaDNA14", "tuaDNA17", "tuaDNA18", "tuaDNA19", "tuaDNA1a", "tuaDNA1b", "tuaDNA1c", "tuaDNA1d", "tuaDNA20", "tuaDNA21", "tuaDNA22", "tuaDNA23", "tuaDNA2a", "tuaDNA2b", "tuaDNA3a", "tuaDNA3b", "tuaDNA4", "tuaDNA5", "tuaDNA6", "tuaDNA7", "tuaDNA8", "tuaDNA9", "tuaDNAL-1")
REPLACEMENT = c("Charlie1_tua", "Chompy-2_tua", "DNA-1_Gav", "Harbinger_tua", "Harbinger_tua", "Harbinger_CPB", "hAT_AMi", "hAT_Crp", "hAT_Crp", "hAT_tua", "hAT_Croc", "hAT_Croc", "hAT_Crp", "hAT_tua", "hAT_tua", "hAT_AMi", "hAT_AMi", "hAT_AMi", "hAT_tua", "hAT_tua", "hAT_tua", "hAT_tua", "hAT_tua", "hAT_Croc", "hATN_tua", "hATN_tua", "hATN_tua", "hATN_tua", "MER45_tua", "MER45_tua", "MER45_tua", "MER45_tua", "nhAT_tua", "piggyBac-1_tua", "Polinton_AMi", "Polinton_tua", "Polinton_AMi", "REP131", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA", "tuaDNA")
for(i in 1:length(PATTERN)){
  
  DNA_bps$Superfamily[DNA_bps$Subfamily == PATTERN[i]] = REPLACEMENT[i]
}

# assign colors to the subfamilies
coldna = character()
colfunc <- colorRampPalette(c("#BCBDDC", "#3F007D"))
coldna = c(coldna, colfunc(6)[c(1,2,3,4,5)])
colfunc <- colorRampPalette(c("lightgreen", "darkgreen"))
coldna = c(coldna, colfunc(4))
colfunc <- colorRampPalette(c("#6BAED6", "#08306B"))
coldna = c(coldna, colfunc(4))
colfunc <- colorRampPalette(c("gold", "red3"))
coldna = c(coldna, colfunc(8)[c(1,3,5,7)])

# plot
ggplot(data = DNA_bps, aes(x = as.integer(Divergence), y = Mb, fill = factor(Superfamily))) + geom_bar(stat = "identity") + scale_fill_manual(name = "DNA Subfamilies", values = as.character(coldna)) + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + xlab("Divergence")

# save the plot
ggsave(filename = 'Tuatara_DNA_landscape_bps_white.png', plot = DNA_general_bps, device = "png", width = 50, height = 27, units = "cm", scale = .5)
