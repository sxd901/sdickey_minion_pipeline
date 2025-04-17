# Function to process the output of blast_compile_results_CL.R
# Will: remove domesticated species, process positive controls, process negative controls

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(openxlsx, quietly = TRUE, warn.conflicts = FALSE)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)

# supress warnings - ok with doing bc the scripts works as expected, turn off when modifying
options(warn=-1)



# Define the argument parser
parser <- ArgumentParser()

# Add named arguments
parser$add_argument("--local_blast_sum_all", type="character", help="output of cluster level summary from blast_compile_results_CL.R for local blast")
parser$add_argument("--global_blast_sum_all", type="character", help="output of cluster level summary from blast_compile_results_CL.R for global blast")
parser$add_argument("--barcode_key", type = "character", help = "seqkit tabulate before clustering")
parser$add_argument("--sample_key", type="character", help="excel sheet of sample key: barcode code to sample name and description")
parser$add_argument("--primer_name", type="character", help="name of primer that matches sample key")
parser$add_argument("--positive_control_list", type="character", help="Path to the positive control list txt file")
parser$add_argument("--contaminant_list", type="character", help="Path to the contaminant species list txt file")
parser$add_argument("--tag_jump_rate", type="numeric", help="user supplied tag jump rate to use (overrides calculated tag jump rate)")
parser$add_argument("--output_dir", type="character", help="dir of csv output files")


# Parse the command line arguments
args <- parser$parse_args()


post_process_blast_compiled_results = 
  function(local_blast_sum_all, # output of cluster level summary from blast_compile_results_CL.R for local blast
           global_blast_sum_all, # output of cluster level summary from blast_compile_results_CL.R for global blast
           barcode_key, # seqkit tabulate
           sample_key, # excel sheet of sample key: barcode code to sample name and description
           primer_name, # name of primer that matches sample key
           positive_control_list, # Path to the positive control list txt file
           contaminant_list, # Path to the contaminant species list txt file
           tag_jump_rate = NULL, # user supplied tag jump rate to use (overrides calculated tag jump rate)
           output_dir # dir of csv output files
          )
  { 


  # Cleaning sequences ---- !!!! DO THIS ON COMBINED LOCAL AND GLOBAL BLAST COMBINED 
  # 1. Remove domestic animal species
  # 2. account for positive controls
  # 3. account for negative controls
  # 4. re export summary by cluster and summary by taxa
  
    
  # what to do if global blast empty? maybe make sure it still writes a file even if empty
  
  # ---------------------------------------------------------------------------- #
    
  # make output dirs for pos/neg control information
  contam_dir = paste0(output_dir, "/contamination")
  positive_control_dir = paste0(output_dir, "/positive_control")
  negative_control_dir = paste0(output_dir, "/negative_control")
  read_distrib_dir = paste0(output_dir, "/sample_read_distribs")
  dir.create(contam_dir, showWarnings = F)
  dir.create(positive_control_dir, showWarnings = F)
  dir.create(negative_control_dir, showWarnings = F)
  dir.create(read_distrib_dir, showWarnings = F)
  
  # ---------------------------------------------------------------------------- #
    
  # combine 98 and 95 cluster summaries as well as local and global blast results ----
  df1 = read.csv(local_blast_sum_all)
  df2 = read.csv(global_blast_sum_all)
  blast_sum9895 = bind_rows(df1 %>% mutate(max_pident = as.character(max_pident),
                                           align_len = as.character(align_len),
                                           qstart = as.character(qstart),
                                           qend = as.character(qend),
                                           taxid = as.character(taxid)), 
                            df2 %>% mutate(max_pident = as.character(max_pident),
                                           align_len = as.character(align_len),
                                           qstart = as.character(qstart),
                                           qend = as.character(qend),
                                           taxid = as.character(taxid))
                            )  

  
  # load in sample key ----
  sample_key_df = read.xlsx(sample_key, na.strings = c("-", "NA"))
  barcode_code_key = get_barcode_code_to_number_key()
  sample_barcode_key =
    sample_key_df %>%
    pivot_longer(cols = starts_with("Primer", ignore.case = TRUE), names_to = "primer", values_to="code") %>%
    filter(!is.na(code)) %>%
    left_join(barcode_code_key, by = c("code"="key"))

  ## subset key to just have barcodes pertaining to the primer/run
  sample_barcode_key_primer =
    sample_barcode_key %>%
    filter(str_detect(primer, regex(paste0("Primer\\.*", primer_name), ignore_case = T))) %>%
    # add in unclassified as a 'sample'
    bind_rows(tibble(barcode = "unclassified", Sample.ID = "unclassified", Sample.Type = "unclassified"))
  
  
  # get read distrib by sample stats 
  read_distrib_by_sample1 = get_read_distrib_by_sample(blast_sum9895, Sample.Type,
                                                       "01_read_distrib_by_sample_pre_post_processing_with_pos_neg", 
                                                       read_distrib_dir, 
                                                       include_pos_neg = TRUE, 
                                                       na_sample_export = TRUE, barcode_code_key)  
  
  read_distrib_by_sample2 = get_read_distrib_by_sample(blast_sum9895, Sample.Type,
                                                       "02_read_distrib_by_sample_pre_post_processing", 
                                                       read_distrib_dir)
  
  
  # ---------------------------------------------------------------------------- #
  
  # remove domestic/contaminate animal species ----
  
  ## load in domestic animal species  (adjust as needed per location)
  add_word_boundary <- function(word) {
    paste0("\\b", word, "\\b")
  }
  
  # load in contaminant species list
  contam_sp = readLines(contaminant_list, warn = FALSE)

  # get genus names of contaminant species
  contam_genus_names = unique(str_extract(contam_sp, "^[A-z]+(?=\\s)"))

  # add word boundaries
  domesticated_contam_sp = paste0(paste0(add_word_boundary(contam_genus_names)), collapse = "|")
  
  
  ## subset clusters that contain these species
  blast_sum9895b_1 = 
    blast_sum9895 %>% 
    filter(!str_detect(blast_scinames, domesticated_contam_sp))
  
  
  ## export report of clusters/taxa removed
  
  ### which clusters removed
  contams_cluster = 
    blast_sum9895 %>% 
    filter(str_detect(blast_scinames, domesticated_contam_sp)) 
  
  ### summary by n reads/samples
  contams_taxa = 
    blast_sum9895 %>% 
    filter(str_detect(blast_scinames, domesticated_contam_sp)) %>% 
    group_by(sscinames_sp, Sample.ID) %>%
    summarise(n_reads = sum(n_reads_barcode), .groups = 'drop') %>% 
    pivot_wider(names_from = Sample.ID, values_from = n_reads) %>% 
    arrange(sscinames_sp) %>% 
    # add in total n reads
    mutate(sum_reads = rowSums(across(setdiff(colnames(.), "sscinames_sp")), na.rm = T), .after=sscinames_sp) %>% 
    # fill na with 0
    replace(is.na(.), 0)
  
  
  
  ### export
  write.csv(contams_cluster, paste0(contam_dir, "/", primer_name, "_contaminants_clusters.csv"), row.names = F)
  write.csv(contams_taxa, paste0(contam_dir, "/", primer_name, "_contaminants_taxa.csv"), row.names = F)
  
  
  # get read distrib by sample stats
  read_distrib_by_sample3 = get_read_distrib_by_sample(blast_sum9895b_1, Sample.Type,
                                                    "03_read_distrib_by_sample_post_contamination", 
                                                    read_distrib_dir)
  
  
  # ---------------------------------------------------------------------------- #
  
  # # assign any NA samples to unclassified ----
  # # BUT THEN need to re-summarize by cluster-sample bc there could now be repeat clusters inside unclassifieds

  # blast_sum9895b = blast_sum9895b_1
  blast_sum9895b = 
    blast_sum9895b_1 %>%
    mutate(Sample.ID = ifelse(is.na(Sample.ID), "unclassified", Sample.ID)) %>% 
    group_by(Sample.ID, id) %>% 
    summarize(n_reads_barcode = sum(n_reads_barcode), 
              barcode = paste0(unique(barcode), collapse = "; "),
              .groups = 'drop') %>% 
    ungroup() %>% 
    # join in original cols
    left_join(blast_sum9895b_1 %>% 
                mutate(Sample.ID = ifelse(is.na(Sample.ID), "unclassified", Sample.ID)) %>% 
                select(-n_reads_barcode, -barcode) %>% 
                distinct(id, Sample.ID, .keep_all = TRUE), 
              by = c("Sample.ID", "id"))
  ## original col order
  blast_sum9895b = blast_sum9895b[, colnames(blast_sum9895b_1)]
  

  
  # sum(blast_sum9895b_1$n_reads_barcode, na.rm = T) == sum(blast_sum9895b$n_reads_barcode, na.rm = T)
  
  ## read distrib
  read_distrib_by_sample4 = get_read_distrib_by_sample(blast_sum9895b, Sample.Type,
                                                       "04_read_distrib_by_sample_mv_na_barcode_unclassified", 
                                                       read_distrib_dir)
  
  
  
  # ---------------------------------------------------------------------------- #
  
  # process positive control ----
  
  # load in positive control names
  positive_control_species = readLines(positive_control_list, warn = FALSE)

  # get genus names from positive control species list
  positive_genus_names = unique(str_extract(positive_control_species, "^[A-z]+(?=\\s)"))
  
  # load in the genus name(s) of positive control (genus not found in amazonas) 
  # for now manually input, but will be on a sep sheet in excel sample barcode key (could be a parameter?)
  # note this is used to search through blast_sciname column
  positive_control_names =  paste0(paste0(add_word_boundary(positive_genus_names), collapse = "|"))
  
  ## load barcode key to get total n reads in a barcode/sample
  # load in barcode key (tabulated from seqkit) 
  # could pre process in command line if too slow here (by converting space into \t?)
  barcode_key1_raw = fread(barcode_key, sep="\t", header = F)
  barcode_key1 = barcode_key1_raw
  barcode_key1 = barcode_key1 %>% select(V1)
  barcode_key1$id = str_extract(barcode_key1$V1, "^[\\w-]+(?=\\s)")
  barcode_key1$barcode = str_extract(barcode_key1$V1, "(?<=barcode\\=)[A-z0-9=]+")
  # barcode_key1$start_time = str_extract(barcode_key1$V1, "(?<=start_time\\=).+(?=\\s)")
  barcode_key1 = barcode_key1 %>% 
    select(-V1)
  
  # # check dupes
  # barcode_key1 %>% group_by(id) %>% filter(n()>1)
  # barcode_key1 %>% count(barcode)
  
  rm(barcode_key1_raw)
  

  ## run positive control correction
  blast_sum9895c = positive_control_correction(blast_sum9895b, 
                                               barcode_key1,
                                               sample_barcode_key_primer,
                                               positive_control_names,
                                               10000,
                                               tag_jump_rate,
                                               primer_name,
                                               positive_control_dir)
  
  
  # # need to rm pos control spp for export (only?)
  # blast_sum9895c_no_pos_sp = 
  #   blast_sum9895c %>% 
  #   filter(!str_detect(blast_scinames, positive_control_names))
  
  
  ## read distrib
  read_distrib_by_sample5 = get_read_distrib_by_sample(blast_sum9895c, Sample.Type.Original,
                                                       "05_read_distrib_by_sample_post_positive_control", 
                                                       read_distrib_dir)
  
  
  # ---------------------------------------------------------------------------- #
  
  
  # process negative controls ----
  blast_sum9895e = negative_control_correction(blast_sum9895c, 
                                              sample_barcode_key_primer, 
                                              primer_name,
                                              positive_control_names,
                                              negative_control_dir)
  
  ## add in unique ID for MOTUS
  blast_sum9895f = blast_sum9895e %>% 
    # add 100000 zero padded id
    group_by(id) %>% 
    mutate(temp1 = cur_group_id()) %>% 
    ungroup() %>% 
    mutate(temp2 = sprintf("%06d", temp1)) %>%
    mutate(motu_id = paste0(primer_name, "_motu_", temp2), .before = Sample.ID) %>% 
    select(-temp1, -temp2)
  
  
  ## read distribs
  read_distrib_by_sample6 = get_read_distrib_by_sample(blast_sum9895f, Sample.Type.Original,
                                                       "06_read_distrib_by_sample_post_negative_control", 
                                                       read_distrib_dir)  
  
  
  
  # ---------------------------------------------------------------------------- #
  
  # summarize by blast_sciname_species ----
  
  ## split by threshold
  blast_sum_all = 
    blast_sum9895f %>% 
    group_by(blast_threshold) %>% 
    group_split()
  
  
  ## Harmonize MOL taxonomy here
  
  
  ## summarize by threshold
  sum_taxa = map_df(blast_sum_all, ~blast_summarize_by_taxa(.x)) %>% 
    arrange(class, order, family, blast_scinames_species) %>% 
    mutate(temp1 = row_number()) %>% 
    mutate(temp2 = sprintf("%05d", temp1)) %>%
    mutate(id = paste0(primer_name, "_taxon_", temp2), .before = primer) %>% 
    select(-temp1, -temp2, -motu_id) %>% 
    # fill in zeros where NA for sample columns (all columns after species)
    mutate(across(all_of(
      colnames(.)[(which(colnames(.) == "species")+1):length(colnames(.))]
      ), ~replace_na(., 0))
    ) # %>% 
    # # add in col for tree png
    # mutate(tree_file_name = paste0(id, ".png"), .after=species)
  
  ## now convert to proportion by sample instead of read count
  sample_totals = 
    sum_taxa %>% 
    pivot_longer(cols = colnames(.)[(which(colnames(.) == "species")+1):length(colnames(.))], 
                 names_to = "sample", values_to = "n") %>% 
    group_by(sample) %>%
    summarize(n_reads_in_sample = sum(n), .groups = 'drop')
  
  sum_taxa_prop = 
    sum_taxa %>% 
    pivot_longer(cols = colnames(.)[(which(colnames(.) == "species")+1):length(colnames(.))], 
                 names_to = "sample", values_to = "n") %>%
    left_join(sample_totals, by = "sample") %>% 
    mutate(prop = n/n_reads_in_sample) %>% 
    select(-n_reads_in_sample, -n) %>% 
    pivot_wider(names_from = sample, values_from = prop)

  
  # ---------------------------------------------------------------------------- #
  
  # summarize by lca name ----
  sum_lca = sum_taxa %>% 
    group_by(lca_sciname, lca_rank) %>%
    summarise(n_reads_corrected = sum(n_reads_total_corrected), 
              blast_refdb = paste0(unique(sort(unlist(str_split(blast_refdb, "; ")), decreasing = TRUE)), collapse = "; "),
              .groups = 'drop') %>% 
    left_join(sum_taxa %>% filter(!is.na(lca_sciname)) %>% 
                                    distinct(lca_sciname, class, order, family), 
              by = c("lca_sciname" = "lca_sciname")) %>% 
    select(class:family, lca_sciname, lca_rank, n_reads_corrected, everything()) %>%
    arrange(class, order, family, lca_sciname) %>% 
    select(blast_refdb, everything())
  
  
  
  # ---------------------------------------------------------------------------- #
  
  # export ---- 
  write.csv(blast_sum9895f, paste0(output_dir, "/", primer_name, "_blast_clusters_corrected.csv"), row.names = F)
  write.csv(sum_taxa, paste0(output_dir, "/", primer_name, "_sum_taxa.csv"), row.names = F)
  write.csv(sum_taxa_prop, paste0(output_dir, "/", primer_name, "_sum_taxa_prop.csv"), row.names = F)
  write.csv(sum_lca, paste0(output_dir, "/", primer_name, "_sum_lca.csv"), row.names = F)
  


  # create excel friendly versions of sum_taxa, sum_taxa_prop bc sometimes ids are too long for excel cells

  ## add / to end of output_dir if not there - OK to overwrite bc end of script
  if(!str_detect(output_dir, "/$")){
    output_dir = paste0(output_dir, "/")
  }

  ## sum_taxa
  fix_excel_overflow(csv_filepath = paste0(output_dir, "/", primer_name, "_sum_taxa.csv"),
                     output_dir = output_dir)

  ## sum_taxa_prop
  fix_excel_overflow(csv_filepath = paste0(output_dir, "/", primer_name, "_sum_taxa_prop.csv"),
                     output_dir = output_dir,
                     export_keys = FALSE)
  
  # ---------------------------------------------------------------------------- #
}



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# function to create excel friendly versions of sum_taxa, sum_taxa_prop bc sometimes ids are too long for excel cells
fix_excel_overflow = function(csv_filepath, output_dir, export_keys = TRUE){
  
  # Read in the data
  df = read.csv(csv_filepath)
  
  # create key for ids column 
  df_ids = df %>% select(id, ids) %>% 
    separate_longer_delim(ids, delim = "; ")
  
  # create key for motu_ids col
  df_motu_ids = df %>% select(id, motu_ids) %>% 
    separate_longer_delim(motu_ids, delim = "; ")
  
  # replace ids with key
  df2 = 
    df %>% 
    select(-ids, -motu_ids)
  
  
  # get filename, not including extension
  filename = str_remove(basename(csv_filepath), ".csv")

  # write to csv
  write.csv(df2, file.path(output_dir, paste0(filename, "_fixed.csv")), row.names = FALSE)
  if(export_keys){
    write.csv(df_ids, file.path(output_dir, paste0(filename, "_ids_key.csv")), row.names = FALSE)
    write.csv(df_motu_ids, file.path(output_dir, paste0(filename, "_motu_ids_key.csv")), row.names = FALSE)
  }
  
}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# function to convert barcode code to number
get_barcode_code_to_number_key = function(){
  bind_rows(
    tibble(key = paste0("A", 1:12), barcode = paste0("barcode", sprintf("%02d", 1:12))),
    tibble(key = paste0("B", 1:12), barcode = paste0("barcode", sprintf("%02d", 13:24))),
    tibble(key = paste0("C", 1:12), barcode = paste0("barcode", sprintf("%02d", 25:36))),
    tibble(key = paste0("D", 1:12), barcode = paste0("barcode", sprintf("%02d", 37:48))),
    tibble(key = paste0("E", 1:12), barcode = paste0("barcode", sprintf("%02d", 49:60))),
    tibble(key = paste0("F", 1:12), barcode = paste0("barcode", sprintf("%02d", 61:72))),
    tibble(key = paste0("G", 1:12), barcode = paste0("barcode", sprintf("%02d", 73:84))),
    tibble(key = paste0("H", 1:12), barcode = paste0("barcode", sprintf("%02d", 85:96)))
  )
}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# function to get read distrib by sample and export
get_read_distrib_by_sample = function(df, Sample.Type_col, name, output_dir, 
                                      include_pos_neg = FALSE, 
                                      na_sample_export = FALSE, barcode_code_key = NULL){
  
  
  ## read distrib by sample
  if(!include_pos_neg){
    df = 
      df %>% 
      filter(!str_detect({{Sample.Type_col}}, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T)) |
               is.na({{Sample.Type_col}})) 
  }
  
  sum_all_reads = 
    df %>%
    summarise(n_reads_total = sum(n_reads_barcode, na.rm = T)) %>% 
    pull(n_reads_total)
  
  read_distrib_by_sample =
    df %>% 
    group_by(Sample.ID) %>% 
    summarise(n_reads_total = sum(n_reads_barcode, na.rm = T), 
              prop_reads = round(n_reads_total/sum_all_reads, 4),
              .groups = 'drop') %>% 
    arrange(-n_reads_total) %>% 
    bind_rows(tibble(Sample.ID = "Total", n_reads_total = sum(.$n_reads_total)))
  
  
  ## export
  write.csv(read_distrib_by_sample, 
            paste0(output_dir, "/", name, ".csv"), row.names = F)
  
  
  if(na_sample_export){
    ## read distrib of barcodes where sample = NA
    read_distrib_na =
      df %>% 
      filter(is.na(Sample.ID)) %>% 
      group_by(barcode) %>% 
      summarise(n_reads_total = sum(n_reads_barcode, na.rm = T)) %>% 
      left_join(barcode_code_key, by = c("barcode")) %>% 
      select(key, barcode, everything()) %>% 
      arrange(-n_reads_total)
    
    ## export
    write.csv(read_distrib_na, 
              paste0(output_dir, "/", name, "_NA.csv"), row.names = F)
  }
  
  
  
  # return for debugging
  return(read_distrib_by_sample)
  
}



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #



# function to correct for positive controls
positive_control_correction = function(blast_sum9895b, 
                                       barcode_key1,
                                       sample_barcode_key_primer,
                                       positive_control_names, 
                                       n_reads_min_sample, # n reads min for a sample to be considered in pos. control calc
                                       tag_jump_rate = NULL, # user supplied tag jump rate to use (overrides calculated tag jump rate)
                                       primer_name,
                                       positive_control_dir){
  
  ## how positive control works is :
  ### Get the % of reads within each sample for each cluster
  ### See what the largest % is across all samples (except negative/positive control & unclassified) for all positive control species
  ### Move any cluster that is =< that % & move those n reads to unclassified.
  ### Remove pos controls from data
  ### Re-summarize by cluster-sample bc potentially multiple of the same clusters could be moved to unclassified
  
  
  # 1. make df of prop of reads per sample for each taxa
  n_reads_per_sample = 
    barcode_key1 %>% 
    group_by(barcode) %>% 
    summarize(n_reads_entire_barcode = n_distinct(id), .groups='drop') %>% 
    left_join(sample_barcode_key_primer %>% select(barcode, Sample.ID), by = "barcode") %>% 
    filter(!is.na(Sample.ID)) %>% 
    select(-barcode)
  
  blast_sum_prop = 
    blast_sum9895b %>% 
    left_join(n_reads_per_sample, 
              by = "Sample.ID") %>% 
    # get the percent of reads per cluster-sample
    group_by(Sample.ID, id) %>% 
    mutate(prop_reads = n_reads_barcode/n_reads_entire_barcode) %>% 
    ungroup() 
  
  
  # 2. get the max prop of pos control species across all Sample.ID (excluding pos/negs controls)
  blast_sum_pc_prop = 
    blast_sum_prop %>% 
    # add in if a positive control species
    mutate(positive_control = ifelse(str_detect(blast_scinames, positive_control_names), 
                                     "positive_control", "not_positive_control")) %>% 
    # filter out positive control samples & Negative controls
    filter(!str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T))) %>% 
    # filter out unclassified for positive control max prop consideration
    filter(Sample.ID != "unclassified") %>%
    # get the max prop_reads for any positive control species
    filter(positive_control == "positive_control") %>% 
    select(barcode, primer, blast_threshold, blast_refdb, id, n_reads_entire_barcode, 
           n_reads_barcode, prop_reads, everything())
  
  ##  get max prop for each Sample.ID
  blast_sum_pc_max_prop = 
    blast_sum_pc_prop %>% 
    group_by(Sample.ID) %>%
    summarize(prop_reads = max(prop_reads), .groups='drop') %>% 
    left_join(blast_sum_pc_prop %>% 
                select(Sample.ID, prop_reads, id, n_reads_entire_barcode, n_reads_barcode, everything()) %>% 
                distinct(), 
              by = c("Sample.ID", "prop_reads")) %>% 
    # remove any ties for max prop in a sample
    group_by(Sample.ID) %>% 
    slice(1) %>% 
    ungroup() %>% 
    arrange(prop_reads)
  
  
  ## get max overall prop AKA tag jump rate

  ### if user supplies value, use that
  if(!is.null(tag_jump_rate)){

    max_pos_prop = tag_jump_rate

  ### otherwise try to calc
  }else{

    if(nrow(blast_sum_pc_max_prop) > 0){
      max_pos_prop = 
        blast_sum_pc_max_prop %>% 
        filter(n_reads_entire_barcode >= n_reads_min_sample) %>% 
        summarise(max_prop_reads = max(prop_reads, na.rm = T)) %>%
        pull(max_prop_reads)
    }else{
      # !!! Hard code in a tag jumping rate when not enough data to get one
      max_pos_prop = 6.65e-5
      message("Not enough data to calculate tag jumping rate, using default value of 6.65e-5")
    }
    
    # !!! Hard code in a tag jumping rate when not enough data to get one
    if(is.infinite(max_pos_prop)){
      max_pos_prop = 6.65e-5
      message("Tag jumping rate is infinite - AKA not enough data to calculate it, using default value of 6.65e-5")
    }
    


  }


  
  # 3. Make positive control correction, move any clusters under this % to unclassified (by sample)
  # but not from pos/neg controls 
  # and rm pos control species clusters after done
  blast_sum_pc_corrected = 
    blast_sum_prop %>% 
    mutate(Sample.ID.Original = Sample.ID) %>% 
    mutate(Sample.ID = 
             case_when(
               !str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T)) & 
                 prop_reads <= max_pos_prop ~ "unclassified",
               TRUE ~ Sample.ID)) %>%
    mutate(under_pos_control_threshold = case_when(
      !str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T)) & 
        prop_reads <= max_pos_prop ~ "yes",
      TRUE ~ "no"))
  
  ## only select original cols (bc sample type needs to be re-categorized)
  blast_sum_pc_corrected2 = blast_sum_pc_corrected[, c(colnames(blast_sum_prop), "Sample.ID.Original", "under_pos_control_threshold")] %>% 
    select(-n_reads_entire_barcode, -prop_reads)
  
  
  # 4. redo unqiue cluster-sample combos now that more (potentially repeated clusters) have moved into unclassified
  blast_sum_pc_corrected3 = 
    blast_sum_pc_corrected2 %>% 
    group_by(Sample.ID, id) %>% 
    summarize(n_reads_barcode = sum(n_reads_barcode), 
              barcode = paste0(unique(sort(unlist(str_split(barcode, "; ")))), collapse = "; "), # need split here bc i combined barcodes above, code good i tested on test df
              Sample.ID.Original = paste0(unique(sort(Sample.ID.Original)), collapse = "; "),
              Sample.Type.Original = paste0(unique(sort(Sample.Type)), collapse = "; "),
              .groups = 'drop') %>% 
    # join in original cols
    left_join(blast_sum_pc_corrected2 %>% 
                select(-n_reads_barcode, -barcode, -Sample.Type, -Sample.ID.Original) %>% 
                distinct(id, Sample.ID, .keep_all = TRUE), 
              by = c("Sample.ID", "id")) %>% 
    select(
      Sample.ID, primer:blast_refdb, id, blast_scinames, blast_comnames,
      blast_n_species, superkingdom:max_pident, 
      n_reads_barcode, align_len:taxid, 
      Sample.ID.Original, Sample.Type.Original, under_pos_control_threshold, barcode
    )
  
  # ---------------------------------------------------------------------------- #
  
  # Write exports

  ## max proportion for each sample
  write.csv(blast_sum_pc_max_prop, 
            paste0(positive_control_dir, "/", primer_name, "_max_pos_prop_per_sample.csv"), row.names = F)
  ## all positive controls in samples
  write.csv(blast_sum_pc_prop, 
            paste0(positive_control_dir, "/", primer_name, "_all_pos_controls_in_samples.csv"), row.names = F)
  ## which sample clusters were moved to unclassifieds
  write.csv((blast_sum_pc_corrected2 %>% 
               filter(under_pos_control_threshold == "yes") %>% 
               filter(Sample.ID.Original != "unclassified") %>% 
               mutate(row = row_number(), .before = Sample.ID.Original)
             ), 
            paste0(positive_control_dir, "/", primer_name, "_clusters_moved_to_unclassified.csv"), row.names = F)
  ## summary by sample, n clusters, n reads
  write.csv((blast_sum_pc_corrected2 %>% 
               filter(under_pos_control_threshold == "yes") %>% 
               filter(Sample.ID.Original != "unclassified") %>% 
              group_by(Sample.ID.Original) %>% 
              summarise(n_clusters = n_distinct(id), 
                        n_reads = sum(n_reads_barcode), .groups = 'drop') %>% 
               arrange(-n_reads) %>% 
               mutate(row = row_number(), .before = Sample.ID.Original)
             ),
            paste0(positive_control_dir, "/", primer_name, "_clusters_moved_to_unclassified_by_sample.csv"), row.names = F)
  ## summary of which taxa moved to unclassified
  write.csv((blast_sum_pc_corrected2 %>% 
               filter(under_pos_control_threshold == "yes") %>% 
               filter(Sample.ID.Original != "unclassified") %>% 
               group_by(sscinames_sp) %>%
               summarise(n_clusters = n_distinct(id), 
                         n_reads = sum(n_reads_barcode), .groups = 'drop') %>% 
               arrange(-n_reads) %>% 
               mutate(row = row_number(), .before = sscinames_sp)
             ),
            paste0(positive_control_dir, "/", primer_name,  "_taxa_moved_to_unclassified.csv"), row.names = F)
  

  
  ## summary of species that moved to unclassifieds and are no longer found in samples
  species_moved_to_unclassified =
    blast_sum_pc_corrected2 %>% 
    filter(under_pos_control_threshold == "yes") %>% 
    filter(Sample.ID.Original != "unclassified") %>% 
    filter(lca_rank == "species") %>% 
    filter(blast_threshold == 98)
  
  species_moved_to_unclassified_not_in_samples_sum = 
    species_moved_to_unclassified  %>% 
       filter(!sscinames_sp %in% unique(blast_sum_pc_corrected2 %>% 
                                          filter(Sample.ID.Original != "unclassified") %>% 
                                          filter(under_pos_control_threshold == "no") %>% 
                                          filter(lca_rank == "species") %>% 
                                          filter(blast_threshold == 98) %>% 
                                          pull(sscinames_sp))) %>% 
       group_by(sscinames_sp) %>%
       summarise(n_clusters = n_distinct(id), 
                 n_reads = sum(n_reads_barcode), .groups = 'drop') %>% 
       arrange(-n_reads) %>% 
       mutate(row = row_number(), .before = sscinames_sp)
    
  
  write.csv(species_moved_to_unclassified_not_in_samples_sum,
  paste0(positive_control_dir, "/", primer_name,  "_species_moved_to_unclassified_not_in_samples.csv"), row.names = F)
  
  ## max proportion & basic stats
  write(paste0("The maximum proportion of a positive control species found in a sample. AKA Tag Jump Rate.\n  (If user supplied this value, it will be that value, not the calculated value)\n  (If == 6.65e-5 then not enough data was generated to calculate it):", max_pos_prop, 
               "\n", "Number of clusters from sample(s) moved to unclassified, as a result: ", 
               (blast_sum_pc_corrected %>% filter(under_pos_control_threshold == "yes") %>% nrow()),
               "\n", "Number of reads from sample(s) moved to unclassified, as a result: ", 
               (blast_sum_pc_corrected %>% filter(under_pos_control_threshold == "yes") %>% 
                  summarise(n_reads = sum(n_reads_barcode)) %>% pull()),
               "\n", "Number of species moved to unclassified and now not found in any sample, as a result: ",
               (species_moved_to_unclassified_not_in_samples_sum %>% nrow())
              ), 
  file = paste0(positive_control_dir, "/", primer_name, "_positive_control_info.txt"))
  
  
  # return df to end fun
  return(blast_sum_pc_corrected3)

}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# function to correct for negative controls
negative_control_correction = function(blast_sum9895c, 
                                       sample_barcode_key_primer, 
                                       primer_name,
                                       positive_control_names,
                                       negative_control_dir){
  
  ## how negative control works is :
  ### 1. Each sample has 3+ negative controls associated with it (field, extraction, pcr)
  ###    Each negative control type (field, extraction, pcr) will have their own lookup table with the SampleID that are associated with them
  ### 2. For each sample, we will check the 3+ negative controls and see if any of the clusters are present in the negative controls
  ### 3. If they are present, select the one negative control with the highest number of reads and substract this number from the sample-cluster
  ### 4. if new read count is < 0 , remove ... should i track somehow? -- perhaps make new sample name for these? for now yes
  
  # 1. Get which negatives are associated with each sample
  negatives_by_sample =
    sample_barcode_key_primer %>%
    filter(!str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T)))
  
  ## field
  negatives_field =
    negatives_by_sample %>%
    separate_longer_delim(Field.Negative, delim = regex(";\\s*")) %>%
    select(Sample.ID, Field.Negative) %>% 
    filter(!is.na(Field.Negative))
  
  ## extraction
  negatives_extraction =
    negatives_by_sample %>%
    separate_longer_delim(Extraction.Negative, delim = regex(";\\s*")) %>%
    select(Sample.ID, Extraction.Negative) %>% 
    filter(!is.na(Extraction.Negative))
  
  ## pcr
  ### get list of cols to filter by (if )
  pcr_cols = negatives_by_sample %>% select(matches(paste0("PCR\\.Negative\\.", primer_name))) %>% colnames()
  
  negatives_pcr =
    negatives_by_sample %>%
    separate_longer_delim(matches(paste0("PCR\\.Negative\\.", primer_name)), delim = regex(";\\s*")) %>%
    select(Sample.ID, matches(paste0("PCR\\.Negative\\.", primer_name))) %>% 
    filter(rowSums(is.na(select(., all_of(pcr_cols)))) == 0)
  
  
  ## combine all into one df
  negatives_all_long =
    bind_rows(
      negatives_field %>% pivot_longer(cols = -Sample.ID, names_to = "type", values_to = "negative_id"),
      negatives_extraction %>% pivot_longer(cols = -Sample.ID, names_to = "type", values_to = "negative_id"),
      negatives_pcr %>% pivot_longer(cols = -Sample.ID, names_to = "type", values_to = "negative_id")
    ) %>% 
    # add in unclassifieds as having all of the above negatives
    bind_rows(tibble(Sample.ID = "unclassified", 
                     type = "unclassified", 
                     negative_id = setdiff(unique(.$negative_id), "unclassified")
                     ))
  
  
  # 2. get sample-cluster combos w.read count in negative controls (should already be unique)
  reads_per_cluster_negative =
    blast_sum9895c %>%
    filter(Sample.ID %in% unique(negatives_all_long$negative_id)) %>%
    rename(negative_id = Sample.ID, n_reads_negative = n_reads_barcode) %>% 
    select(id, negative_id, n_reads_negative) 
  
  
  # 3. for each sample-cluster in samples, join all negative controls and select the one with the max reads
  sample_cluster_neg_read_counts = 
    blast_sum9895c %>% 
    # filter out all negatives and positives controls
    filter(!Sample.ID %in% 
             (sample_barcode_key_primer %>%
                filter(str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T))) %>% 
                pull(Sample.ID))
    ) %>% 
    # join to the sample-clusters in the negative controls
    # many to many bc blast_sum can have multiple clusters (if cluster is found across multiple barcodes)
    # ............ and because reads_per_cluster_negative can have multiple clusters found across multiple negative controls
    left_join(reads_per_cluster_negative, by = c("id"), relationship = "many-to-many") %>% 
    # if n_reads_negative is NA, it means the cluster was not found in any negative control set to 0
    mutate(n_reads_negative = ifelse(is.na(n_reads_negative), 0, n_reads_negative)) %>%
    # now get the max n_reads_negative within each sample cluster combo. IE when more than one negative control tied to the sample-cluster had the cluster present
    group_by(Sample.ID, id) %>%
    summarise(max_n_reads_negative = max(n_reads_negative, na.rm = T), .groups = 'drop') %>% 
    # note if sample is NA that means there were reads in the barcode that were not assigned to any sample or pos/neg control
    # meaning they can be removed here bc they are not tied to a sample or control
    filter(!is.na(Sample.ID))
  
  
  # 4. now subtract the max reads from the sample-cluster in the samples
  blast_sum9895d = 
    blast_sum9895c %>% 
    # filter out all negatives and positives controls
    filter(!Sample.ID %in% 
             (sample_barcode_key_primer %>%
                filter(str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T))) %>% 
                pull(Sample.ID))
    ) %>% 
    # filter out positive control species
    filter(!str_detect(blast_scinames, positive_control_names)) %>% 
    # join in the max neg number of reads per sample-cluster
    left_join(sample_cluster_neg_read_counts, by = c("Sample.ID", "id")) %>% 
    # filter(!is.na(Sample.ID)) %>% # no longer needed bc NA samples are now unclassified
    # subtract the max neg reads from the sample-cluster
    mutate(n_reads_corrected = n_reads_barcode - max_n_reads_negative, .after = n_reads_barcode) %>%
    # if new read count is < 0, set to 0
    mutate(n_reads_corrected = ifelse(n_reads_corrected < 0, 0, n_reads_corrected))
  
  
  # if new read count is 0, then remove them but track taxon names and check if exists elsewhere in other samples
  blast_sum9895e = 
    blast_sum9895d %>% 
    filter(n_reads_corrected > 0) %>% 
    select(-max_n_reads_negative)
  
  
  # ## double check read count is correct ... it is.
  # blast_sum9895c %>% 
  #   filter(!Sample.ID %in% 
  #            (sample_barcode_key_primer %>%
  #               filter(str_detect(Sample.Type, regex("\\bPositive\\b|\\bNegative\\b",ignore_case = T))) %>% 
  #               pull(Sample.ID))
  #   ) %>% 
  #   filter(!str_detect(blast_scinames, positive_control_names)) %>% 
  #   summarise(n_reads_total = sum(n_reads_barcode, na.rm = T)) %>% 
  #   pull() -
  #   blast_sum9895d %>% 
  #   mutate(n = ifelse(max_n_reads_negative < 0, 0, max_n_reads_negative)) %>% 
  #   summarise(n_reads_total = sum(n, na.rm = T)) %>% 
  #   pull()
  # 
  # sum(blast_sum9895e$n_reads_corrected, na.rm = T)
  
  
  # ---------------------------------------------------------------------------- #
  
  # Summary/Reporting Negative Controls
  
  ## negative controls
  
  ### which clusters were in negative controls
  clusters_in_negative_controls = 
    blast_sum9895c %>%
    filter(Sample.ID %in% unique(negatives_all_long$negative_id)) 
  
  ### negative summary by taxa-control
  negative_control_summary_by_taxa = 
    clusters_in_negative_controls %>% 
    group_by(sscinames_sp, Sample.ID) %>% 
    summarise(n_reads = sum(n_reads_barcode), .groups = 'drop')
  
  ### negative control summary by control
  negative_control_summary = 
    clusters_in_negative_controls %>% 
    group_by(Sample.ID) %>% 
    summarise(n_reads = sum(n_reads_barcode), .groups = 'drop')
  
  ## samples 
  
  ## which clusters got subtracted and by how much
  clusters_in_samples_that_got_negatived = 
    blast_sum9895d %>% 
    filter(max_n_reads_negative > 0)
  
  ## which clustered got zeroed out
  clusters_got_zeroed_out = 
    blast_sum9895d %>% 
    filter(n_reads_corrected == 0)
  
  
  ## which taxa got zeroed out that are no longer present in the data
  taxa_zeroed_out_not_in_samples =
    clusters_got_zeroed_out %>% 
    filter(!sscinames_sp %in% unique(blast_sum9895e$sscinames_sp)) %>%
    select(id, barcode, Sample.ID, sscinames_sp, blast_comnames, n_reads_barcode, max_n_reads_negative) %>% 
    distinct()

  
  ### overall summary of negatives
  negative_control_txt_summary = paste0(
    
    ## about negative controls themselves
    "Which negative controls had reads: ", paste0(negative_control_summary$Sample.ID, collapse = ", "),
    "\n", "How many total reads in the negatives: ", sum(negative_control_summary$n_reads),
    "\n", "Average number of reads in the negatives: ", round(mean(negative_control_summary$n_reads, na.rm = T), 1),
    "\n", "Number of clusters in the negatives: ", nrow(clusters_in_negative_controls),
    "\n", "Average number of reads per cluster in negatives: ", round(mean(clusters_in_negative_controls$n_reads_barcode, na.rm = T), 1),
    "\n", "How many total taxa in the negatives: ", length(unique(negative_control_summary_by_taxa$sscinames_sp)),
    # "\n", "How many total taxa in the negatives that are also in the Samples: ",
    # (negative_control_summary_by_taxa %>% filter(sscinames_sp %in% unique(blast_sum9895d$sscinames_sp)) %>% 
    #    pull(sscinames_sp) %>% unique() %>% length()),
    "\n", "Average number of reads per taxa in negatives: ", round(mean(negative_control_summary_by_taxa$n_reads, na.rm = T), 1),
    "\n", "How many total species were in the negatives: ", 
    (negative_control_summary_by_taxa %>% filter(!str_detect(sscinames_sp, ";")) %>% 
       pull(sscinames_sp) %>% unique() %>% length()),
    "\n", "Average number of reads per species in negatives: ", 
    (negative_control_summary_by_taxa %>% filter(!str_detect(sscinames_sp, ";")) %>% 
       pull(n_reads) %>% mean(na.rm = T) %>% round(1)),
    
    ### about samples (what happened to them)
    "\n\n", "How many cluster-SAMPLE combos needed to be corrected for bc a cluster was in negative: ", 
            nrow(clusters_in_samples_that_got_negatived),
    "\n", "Average number of reads corrected for in a cluster-SAMPLE combo (ie avg n reads subtracted): ", 
            round(mean(clusters_in_samples_that_got_negatived$max_n_reads_negative, na.rm = T), 1),
    "\n", "How many clusters got zeroed out: ", nrow(clusters_got_zeroed_out),
    "\n", "Average read count in clusters that got zeroed out: ", 
            round(mean(clusters_got_zeroed_out$n_reads_barcode, na.rm = T), 1),
    "\n", "How many taxa got zeroed out that are no longer present in the data: ",
            length(unique(taxa_zeroed_out_not_in_samples$sscinames_sp)),
    "\n", "Names of taxa that got zeroed out: ",
            paste0(unique(taxa_zeroed_out_not_in_samples$sscinames_sp), collapse = "; ")
  )
  
  
  # export report
  write.csv(clusters_in_negative_controls, paste0(negative_control_dir, "/", primer_name, "_clusters_in_negative_controls.csv"), row.names = F)
  write.csv(negative_control_summary_by_taxa, paste0(negative_control_dir, "/", primer_name, "_negative_control_summary_by_taxa.csv"), row.names = F)
  write.csv(negative_control_summary, paste0(negative_control_dir, "/", primer_name, "_negative_control_summary_by_control.csv"), row.names = F)
  write.csv(clusters_in_samples_that_got_negatived, paste0(negative_control_dir, "/", primer_name, "_clusters_in_samples_that_got_negatived.csv"), row.names = F)
  write.csv(clusters_got_zeroed_out, paste0(negative_control_dir, "/", primer_name, "_clusters_got_zeroed_out.csv"), row.names = F)
  write.csv(taxa_zeroed_out_not_in_samples, paste0(negative_control_dir, "/", primer_name, "_taxa_zeroed_out_not_in_samples.csv"), row.names = F)
  write(negative_control_txt_summary, file = paste0(negative_control_dir, "/", primer_name, "_negative_control_info.txt"))
  
  
  
   
  # ---------------------------------------------------------------------------- #
  
  # return
  return(blast_sum9895e)
  
  
}



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# summarize by blast sciname-sample combo
blast_summarize_by_taxa = function(sum_blast2){
  
  # Summarize by n reads of each taxa by Sample
  sum_taxa_sample = 
    sum_blast2 %>% 
    group_by(Sample.ID, sscinames_sp) %>% 
    summarize(n_reads_corrected = sum(n_reads_corrected), 
              .groups = 'drop') %>% 
    pivot_wider(names_from = Sample.ID, values_from = n_reads_corrected) %>%
    replace(is.na(.), 0)
  
  # get the 'best' (highest pident) consensus sequence for each taxa (sscinames_sp)
  top_con_seq_id = 
    sum_blast2 %>% 
    mutate(row = row_number()) %>% 
    separate_rows(max_pident, sep = "; ") %>%
    mutate(max_pident = as.numeric(max_pident)) %>%
    group_by(sscinames_sp) %>% 
    filter(max_pident == max(max_pident)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    pull(row)
  
  top_con_seq = 
    sum_blast2 %>% 
    mutate(row = row_number()) %>% 
    select(-barcode, -id, -n_reads_corrected, -blast_refdb, -subspecies, -n_reads_barcode) %>%
    filter(row %in% top_con_seq_id)
  
  
  # summarize by taxa
  
  # Notes: Across local/global at the same threshold is fine. They will be the same
  # taxa call (if in same threshold). I would just need to take care to summarize in 
  # a way that would re parse the max pident, sequence, etc. 
  
  # BUT this would not be OK across blast thresholds (98/95) because they could
  # be called to different taxa. E.g. Species A @98 call != Species A @95 call (one is species and one is genus)
  # Even if they are the same, like in the case of a genus call, the 98 could get revised
  # by an expert to change to species. So, we want to keep them separate to account for this.
  
  # Perhaps best way is to parse data by threshold (98 of local and global; 95 of local and global)
  # and run this script twice for each df and combine after
  # it makes sense because we are considering them as separate taxa calls
  
  sum_taxa = 
    sum_blast2 %>% 
    group_by(Sample.ID, sscinames_sp) %>%
    mutate(name_id = cur_group_id()) %>% 
    ungroup() %>% 
    group_by(sscinames_sp) %>%
    summarize(n_clusters = n_distinct(id),
              n_reads_total_corrected = sum(n_reads_corrected),
              ids = paste0(unique(id), collapse = "; "),
              motu_ids = paste0(unique(motu_id), collapse = "; "),
              blast_refdb = paste0(unique(blast_refdb), collapse = "; "),
              .groups = 'drop')  %>% 
    # join back in taxonomy & other info
    left_join(top_con_seq, by = "sscinames_sp") %>% 
    # join in reads per sample
    left_join(sum_taxa_sample, by = "sscinames_sp") %>% 
    rename(blast_scinames_species = sscinames_sp) %>% 
    select(-blast_scinames, -row) %>%
    arrange(class, order, family, blast_scinames_species) %>% 
    select(primer, blast_threshold, blast_refdb,
           class, order, family, 
           lca_sciname, lca_rank, blast_scinames_species, blast_comnames, blast_n_species, 
           max_pident, 
           n_clusters, n_reads_total_corrected,
           sequence, qstart, qend, align_len, sseqid, taxid, ids, motu_ids,
           superkingdom, phylum, genus, species,
           everything()) %>% 
    # mutate(id = row_number(), .before = primer) %>% 
    select(-Sample.ID, -Sample.ID.Original, -Sample.Type.Original, -under_pos_control_threshold) %>%
    # fill zeros for sample cols
    mutate(across(all_of(setdiff(colnames(sum_taxa_sample), "sscinames_sp")), ~ replace_na(.x, 0)))
    
    
    return(sum_taxa)
  
}


# # check to see if any local taxa are in global, there are ... what does this means
# # this means there are seqs in global that should be in local
# # also means that when we summarize by sscinames_sp we are losing this level of granularity
# MIKE IS INTERESTED IN THIS
# localsppinglobal = 
#   blast_sum9895e %>% 
#     filter(blast_refdb == "local") %>% 
#     distinct(sscinames_sp) %>% 
#     filter(sscinames_sp %in% 
#              (blast_sum9895e %>% 
#              filter(blast_refdb == "global") %>% 
#                distinct(sscinames_sp) %>% 
#                pull(sscinames_sp))
#            ) %>% 
#   pull()
# 
# # get their accession #s
# df_temp = 
#   blast_sum9895e %>% 
#   filter(sscinames_sp %in% localsppinglobal) %>% 
#   filter(blast_refdb == "global") %>% 
#   select(blast_refdb, max_pident, sscinames_sp, sseqid) %>% 
#   distinct() %>% 
#   arrange(sscinames_sp)
# 
# write.csv(df_temp, "~/Downloads/local_in_global.csv", row.names = F)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# Call the function with named parameters ----
post_process_blast_compiled_results(args$local_blast_sum_all, 
                                    args$global_blast_sum_all,
                                    args$barcode_key,
                                    args$sample_key, 
                                    args$primer_name, 
                                    args$positive_control_list,
                                    args$contaminant_list,
                                    args$tag_jump_rate,
                                    args$output_dir)