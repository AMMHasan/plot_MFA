# Written by Dr. A M Mahedi Hasan, Post-doctoral Research Associate, The Leach Lab, ICB, University of Edinburgh, UK.
# This script is under GPLv3 licensing criteria.

# This script is a visualisation tool for the plotting replication profile upon Marker Frequency Analysis on the depth of coverage, calculated from Illumina data. This R script is the next step to (MFA) project as an integrated workflow for comparing replication profiles of different Escherichia coli strains.
# Some important genomic loci (landmarks) are displayed with automated annotation.
# The potential is limitless, as you can plot as many graphs you want. The number, at present, is limited by the number of colours (seven) defined in the list. The number can theoritically be increased. However, one should keep in mind that too many overlaying plots will hamper the asthetic value of the plot.
# Another limitation of this script is - it has to be executed line by line.


# Load required packages ... ...
if (!require(seqinr)) install.packages('seqinr')
suppressPackageStartupMessages(library(seqinr))

if (!require(ggplot2)) install.packages('ggplot2')
suppressPackageStartupMessages(library(ggplot2))


# The reference genome file should given in .fasta format and the name of the file should start with "genome", like "genome_ref.fasta".
# mapping some landmarks in the genome ... ...

genome_file <- list.files(pattern = "^genome")
genome <- read.fasta(genome_file)[[1]]
genome_seq <- paste(genome[1:length(genome)], collapse = "")

oriC <- words.pos("GATCTATTTATTTAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGCACTGCCCTGTGGATAACAAGGATCCGGCTTTTAAGATCAACAACCTGGAAAGGATCATTAACTGTGAATGATCGGTGATCCTGGACCGTATAAGCTGGGATCAGAATGAGGGGTTATACACAACTCAAAAACTGAACAACAGTTGTTCTTTGGATAACTACCGGTTGATCCAAGCTTCCTGA", genome_seq, ignore.case=T)
rrsH <- words.pos("AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAA",genome_seq,ignore.case = T)[1]
pal246 <- words.pos("ATGGTCATAGCTGTTTCC",genome_seq,ignore.case = T)
terB <- words.pos("AATAAGTATGTTGTAACTAAAGT",genome_seq,ignore.case = T)
terC <- words.pos("ATATAGGATGTTGTAACTAATAT", genome_seq, ignore.case = T)
terA <- words.pos(c2s(rev(comp(s2c("AATTAGTATGTTGTAACTAAAGT")))), genome_seq, ignore.case = T)
terD <- words.pos(c2s(rev(comp(s2c("CATTAGTATGTTGTAACTAAATG")))), genome_seq, ignore.case = T)

# listing important loci ... ...
landmarks <- c("oriC","terA","terB","terC","terD", "pal246","rrsH")
landmarks_list = data.frame(loci=landmarks,position=c(oriC,terA,terB,terC,terD,pal246,rrsH),row.names = 1)

# Follow the instruction carefully ... ...
print(paste("Several important loci to choose from:", paste(landmarks,collapse = ", ")))
loci_chosen <- readline(prompt = "Please type in the loci name(s) to annotate the graph (comma separated, like - a,b,c):")

individual_locus <- trimws(strsplit(loci_chosen,",")[[1]])
individual_locus_pos <- landmarks_list[individual_locus,]


# colour_list for the plot... ... ...
colour_list <- c("red","blue","green","brown","black","orange", "bluish green")

# read the data files ... ... ...
FW_data <- read.table("normalisation_output/1Kb_fixed_window data.txt")
loess_data <- read.table("normalisation_output/loess_1Kb_fixed_winodw_data.txt")

# print out the options to choose from
print(colnames(loess_data)[2:length(colnames(loess_data))])

# take the input from the keyboard
col_names <- readline(prompt = "Choose from the options given below (with commas, like a,b,c):")
col_names <- trimws(strsplit(col_names,",")[[1]])

# put the plot command in a list as a serise of strings
z=list()
for(i in seq(length(col_names))){
  z[i] <- paste("geom_point(data=FW_data, aes_string(x='position',y='",
                toString(col_names[i]),
                "'),colour='",
                toString(colour_list[i]),
                "', alpha=0.02) + ",
                "geom_line(data = loess_data, aes_string(x='position',y='",
                toString(col_names[i]),
                "'),colour='",
                toString(colour_list[i]),
                "', size=2, alpha=0.5)",
                sep ="")

}

########### plot on the screen and save the combined MFA graph


# y values for the lines marking important loci
ymin <- seq(2,length.out = length(col_names), by=.1)
ymax <- seq(2.02,length.out = length(col_names), by=.1)

# coordinates for the labels for the loci lines
x_text <- 8.5e5
y_text <- ymax





# a vector for ggplot
plt <- eval(parse(text = paste("ggplot() + ",
                               "geom_vline(xintercept=individual_locus_pos,linetype=2) +",
                               paste(z, collapse = "+"),
                               " + xlab('Chromosomal position') + ylab('Normalised depth')+ ylim(0,3) + theme_bw()",
                               "+annotate('rect',xmin=0,xmax=5e5,ymin = ymin,ymax = ymax, fill=colour_list[seq(length(col_names))])",
                               "+annotate('text', x=x_text, y= y_text, label = col_names)",
                               "+annotate('text', x=(individual_locus_pos+1e5), y= seq(2.3,length.out = length(individual_locus), by=.1), label = individual_locus)"
)))

# suppressing the warnings for the missing values. Keeping them makes the users worried for no reason.
suppressWarnings(print(plt))

# same here
suppressWarnings(ggsave("combined_MFA_graph.tiff"))


