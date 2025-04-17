# convert the output of bbmap::readlength.sh to a histogram that is readable
library(argparse, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(scales, quietly = TRUE, warn.conflicts = FALSE)

# supress warnings - ok with doing bc the scripts works as expected, turn off when modifying
options(warn=-1)

# Define the argument parser
parser <- ArgumentParser()

# Add named arguments
parser$add_argument("--distrib_file", type = "character", help = "output file of bbmap::readlength.sh")
parser$add_argument("--output_file", type = "character", help = "name of output png")
parser$add_argument("--hmin", type = "character", help = "minimum read length to highlight")
parser$add_argument("--hmax", type = "character", help = "maximum read length to highlight")

# Parse the command line arguments
args <- parser$parse_args()


# Define your function
read_distrib_plot = function(distrib_file, hmin, hmax, output_file){
  
  
  df = read.table(distrib_file)
  colnames(df) = c("Length", "reads",	"pct_reads",	"cum_reads",	"cum_pct_reads",	"bases",
                   "pct_bases",	"cum_bases",	"cum_pct_bases")
  
  df %>% 
    ggplot(aes(x=Length, y = reads)) + 
    geom_line(color="gray40", linewidth=0.3) + 
    geom_point(alpha = 0.5, size = 0.3) + 
    scale_x_continuous(trans="log10", breaks = log_breaks(12), labels=comma) +
    scale_y_continuous(labels=comma,
                       sec.axis = sec_axis(~ (. / sum(df$reads))*100, 
                                           name = "Percentage of reads")
    ) + 
    annotate(geom = "rect",
             xmin = hmin, xmax = hmax,
             ymin = -Inf, ymax = Inf,
             color = "transparent",
             fill = "blue",
             alpha = .2) +
    labs(y = "Number of reads", x = "Read length") + 
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank())
  
  ggsave(plot=last_plot(), filename = output_file, width = 7, height = 5)
}



# Call the function with named parameters
read_distrib_plot(distrib_file = args$distrib_file, 
                  hmin = as.numeric(args$hmin), 
                  hmax = as.numeric(args$hmax), 
                  output_file = args$output_file)