
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)


# supress warnings - ok with doing bc the scripts works as expected, turn off when modifying
options(warn=-1)



# function to parse primer file for cutadapt, fasta id must have the work Forward or Reverse in it
parse_primer_for_cutadapt = function(fasta_file){
  
  # Function to reverse letters in a string
  reverse_string <- function(x) {
    paste(rev(strsplit(x, "")[[1]]), collapse = "")
  }
  
  # get df of primer seqs, and RC seqs
  lines_temp = readLines(fasta_file, warn=FALSE)
  lines_temp <- lines_temp[lines_temp != ""] # rm blank lines
  
  primer_df = 
    lines_temp %>% 
    matrix(nrow = length(.)/2, byrow = TRUE) %>% 
    as_tibble() %>% setNames(c("name", "sequence")) %>% 
    mutate(name = str_remove(name, "^>") %>% str_trim()) %>%
    mutate(sequence_rc = chartr("ATGCYM","TACGRK",sequence) %>% sapply(reverse_string) ) %>%
    mutate(for_rev = case_when(str_detect(name, regex("forward", ignore_case=TRUE)) ~ "Forward",
                            str_detect(name, regex("reverse", ignore_case=TRUE)) ~ "Reverse")
           )
  
  
  # get all combinations of cutadapt
  
  ## df to long
  primer_df_long = 
    primer_df %>%
    pivot_longer(cols = c("sequence", "sequence_rc"), names_to = "primer_type", values_to = "sequence")
  
  
  ## unlinked
  unlinked_df = 
    primer_df_long %>% 
    mutate(cutadapt_param = 
             case_when(primer_type == "sequence_rc" ~ "-a",
                       primer_type == "sequence" ~ "-g")) %>% 
    rename(cutadapt_seq = sequence) %>% 
    mutate(link_type = "unlinked")
  
  
  ## linked
  linked_combos1 = 
    expand.grid(primer_df %>% filter(for_rev=="Forward") %>% pull(name),
                primer_df %>% filter(for_rev=="Reverse") %>% pull(name))
  
  linked_combos2 = 
    expand.grid(primer_df %>% filter(for_rev=="Reverse") %>% pull(name),
                primer_df %>% filter(for_rev=="Forward") %>% pull(name))
  
  
  linked_df = 
    bind_rows(
      expand.grid(primer_df %>% filter(for_rev=="Forward") %>% pull(sequence), 
                  primer_df %>% filter(for_rev=="Reverse") %>% pull(sequence_rc)),
      
      expand.grid(primer_df %>% filter(for_rev=="Reverse") %>% pull(sequence), 
                  primer_df %>% filter(for_rev=="Forward") %>% pull(sequence_rc))
    ) %>% 
    setNames(c("sequence", "sequence_rc")) %>%
    mutate(name_first_primer =c(linked_combos1[,1] %>% as.character(),
                  linked_combos2[,1] %>% as.character())) %>% 
    mutate(name_second_primer =c(linked_combos1[,2] %>% as.character(),
                                linked_combos2[,2] %>% as.character())) %>%   
    mutate(cutadapt_param = "-g") %>% 
    mutate(cutadapt_seq = paste0(sequence, "...", sequence_rc)) %>% 
    mutate(link_type = "linked")
  
  
  # combine into final df
  linked_unlinked = bind_rows(unlinked_df %>% select(cutadapt_param, cutadapt_seq, link_type), 
                              linked_df %>% select(cutadapt_param, cutadapt_seq, link_type))
  
  return(linked_unlinked)
}




# Construct Cutadapt command function
construct_cutadapt <- function(input_file, output_file, 
                         adapter_type, adapter_sequence, 
                         n_cores, cutadapt_error_rate #, os
                         ) {
  cutadapt_command <- paste(
    # ifelse(os =="linux", "singularity exec images/cutadapt.sif", ifelse(os=="mac", "cutadapt", "")),
    "cutadapt",                  # Cutadapt command
    adapter_type, adapter_sequence,   # Adapter sequence to trim
    "--cores",  n_cores,            # Number of cores to use
    "-e ", cutadapt_error_rate,  # Error rate
    "--no-indels", "--discard-untrimmed",  # Discard untrimmed reads
    "--output", output_file,        # Output file
    input_file,               # Input file
    collapse = " "
  )
  
  return(cutadapt_command)
}




# construct all cutadapt commands & run them & cat the output into linked and unlinked
run_cutadapt = function(input_file, fasta_file,
                       output_path_linked, output_path_unlinked,
                       n_cores, cutadapt_error_rate #, os
                       ){
  
  # parse primer file
  primer_df = parse_primer_for_cutadapt(fasta_file)
  
  ## make temp dir for linked and unlinked
  tempdir_linked = tempfile()
  tempdir_unlinked = tempfile()
  dir.create(tempdir_linked)
  dir.create(tempdir_unlinked)
  
  # construct cutadapt commands
  cutadapt_commands = list()
  for(i in 1:nrow(primer_df)){
    cutadapt_commands[[i]] =
    construct_cutadapt(input_file = input_file, 
                      output_file = case_when(primer_df$link_type[i] =="linked" ~ paste0(tempdir_linked, "/temp", i, ".fastq"), 
                                              primer_df$link_type[i]=="unlinked" ~ paste0(tempdir_unlinked, "/temp", i, ".fastq")), 
                      adapter_type = primer_df$cutadapt_param[i], 
                      adapter_sequence = primer_df$cutadapt_seq[i], 
                      n_cores = n_cores, 
                      cutadapt_error_rate = cutadapt_error_rate #,
                      # os = os
                      )
  }  
  
  # echo number of cutadapt commands and what they are
  cat("Number of cutadapt commands: ", length(cutadapt_commands), "\n")
  cat("Cutadapt commands: ", "\n")
  cat(unlist(cutadapt_commands), sep = "\n")
  
  # echo running cutadapt and the name of the two tempdirs
  cat("Running cutadapt on ", input_file, "\n")
  cat("Tempdir linked: ", tempdir_linked, "\n")
  cat("Tempdir unlinked: ", tempdir_unlinked, "\n")

  # Run all all cutadapt commands
  for(j in 1:length(cutadapt_commands)){
    system(cutadapt_commands[[j]], intern = TRUE)
  }


  # cat linked to output path
  system(paste0("cat ", tempdir_linked, "/* > ", output_path_linked), intern = TRUE)
  
  # cat unlinked to output path
  system(paste0("cat ", tempdir_unlinked, "/* > ", output_path_unlinked), intern = TRUE)

  
}



# Inputs from command line
# input_file, primer_fasta, output_path_linked, output_path_unlinked, cutadapt_error_rate, n_cores


# Define the argument parser
parser <- ArgumentParser()

# Add named arguments
parser$add_argument("--input_file", type = "character", help = "input fastq file to be trimmed")
parser$add_argument("--primer_fasta", type = "character", help = "fasta file with primer sequences")
parser$add_argument("--output_path_linked", type = "character", help = "output path for linked primers")
parser$add_argument("--output_path_unlinked", type = "character", help = "output path for unlinked primers")
parser$add_argument("--cutadapt_error_rate", type = "numeric", help = "error rate for cutadapt")
parser$add_argument("--n_cores", type = "numeric", help = "number of cores to use")
# parser$add_argument("--os", type = "character", help = "operating system, linux or mac")


# Parse the command line arguments
args <- parser$parse_args()



# actual call to function, run_cutadapt
run_cutadapt(input_file = args$input_file, 
             fasta_file = args$primer_fasta,
             output_path_linked = args$output_path_linked, 
             output_path_unlinked = args$output_path_unlinked,
             n_cores = args$n_cores, 
             cutadapt_error_rate = args$cutadapt_error_rate #,
            #  os = args$os
             )