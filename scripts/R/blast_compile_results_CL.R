# process pipeline v2 w/ initial length filter, seqkit tab, vsearch clustering, blastn, and summary
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)
library(openxlsx, quietly = TRUE, warn.conflicts = FALSE)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)
# library(scales, quietly = TRUE, warn.conflicts = FALSE)

# supress warnings - ok with doing bc the scripts works as expected, turn off when modifying
options(warn=-1)



# Define the argument parser
parser <- ArgumentParser()

# Add named arguments
parser$add_argument("--cluster_key", type = "character", help = "vsearch output")
parser$add_argument("--barcode_key", type = "character", help = "seqkit tabulate before clustering")
parser$add_argument("--conseq_key", type = "character", help = "key for cluster id to consensus sequence")
parser$add_argument("--blast_results", type = "character", help = "local blastn output")
parser$add_argument("--sample_key", type = "character", help = "excel sheet of sample key: barcode code to sample name and description")
parser$add_argument("--primer_name", type = "character", help = "name of primer that matches sample key")
parser$add_argument("--min_alignment_length", type = "numeric", help = "minimum alignment length to be considered a match")
parser$add_argument("--local_global", type = "character", help = "string of either 'local' or 'global' depending on the refdb used")
parser$add_argument("--output_cluster_filepath", type = "character", help = "filepath of xlsx for cluster summary")
parser$add_argument("--output_taxa_filepath", type = "character", help = "filepath of xlsx for taxa summary")
parser$add_argument("--clusters_matched_filepath", type = "character", help = "filepath of txt file for clusters that had any local blast matches")
parser$add_argument("--taxdump_filepath", type = "character", help = "filepath of taxdump csv file")
parser$add_argument("--custom_taxa_filepath", type = "character", help = "filepath of custom taxa csv file")

# Parse the command line arguments
args <- parser$parse_args()


pipeline_v2 = function(cluster_key, # vsearch output
                       barcode_key, # seqkit tabulate
                       conseq_key, # key for cluster id to consensus sequence
                       blast_results, # local blastn output
                       sample_key, # excel sheet of sample key: barcode code to sample name and description
                       primer_name, # name of primer that matches sample key
                       min_alignment_length, # minimum alignment length to be considered a match
                       local_global, # string of either "local" or "global" depending on the refdb used
                       output_cluster_filepath, # filepath of csv for cluster summary
                       output_taxa_filepath, # filepath of csv for taxa summary
                       clusters_matched_filepath = FALSE, # filepath of txt file for clusters that had any local blast matches
                       taxdump_filepath, # filepath of taxdump csv file
                       custom_taxa_filepath = FALSE# filepath of custom taxa csv file
                       )
  { 

  
  # ---------------------------------------------------------------------------- #

  # check files exist
  check_file_exists <- function(file_path) {
    if (file.exists(file_path)) {
      print(paste("File exists:", file_path))
    } else {
      print(paste("File does not exist:", file_path))
    }
  }

  check_file_exists(cluster_key)
  check_file_exists(barcode_key)
  check_file_exists(conseq_key)
  check_file_exists(blast_results)
  check_file_exists(sample_key)
  check_file_exists(taxdump_filepath)
  if (is.character(custom_taxa_filepath)) {
    check_file_exists(custom_taxa_filepath)
  }


  # ---------------------------------------------------------------------------- #
  
  # load in cluster key info -----
  # cluster key contains the name of the cluster (named after the centroid sequence)
  # along with the ids of the seqs that belong in the cluster
  message("Reading in cluster key...")
  cluster_key_raw = fread(cluster_key, sep = "\t")
  # https://drive5.com/usearch/manual/opt_uc.html
  colnames(cluster_key_raw) = c("record_type", "cluster_no", "centroid_len", 
                                "p_simil_cen", "match_orientation", "not_used",
                                "not_used2", "compress_align", "query_lab", "con_lab")
  # cluster_key_raw %>% count(record_type)
  ## rm redundant rows
  cluster_key_raw = cluster_key_raw %>% filter(record_type != "S")
  
  
  ## create key
  cluster_key1 = 
    cluster_key_raw %>% 
    select(cluster_no, query_id = query_lab, con_lab) %>% 
    mutate(con_id = ifelse(con_lab == "*", query_id, 
                           con_lab)) %>% 
    select(con_id, everything(), -con_lab)
  
  rm(cluster_key_raw)
  
  # # check dupes
  # cluster_key1 %>% group_by(query_id) %>% filter(n()>1)
  
  # ---------------------------------------------------------------------------- #
  
  # load in consencus sequence key (one seq per cluster)
  message("Reading in conseq_key...")
  conseq_key = read.csv(conseq_key, sep = "\t", header = F)
  colnames(conseq_key) = c("id", "sequence")
  
  # ---------------------------------------------------------------------------- #
  
  # load in barcode key (tabulated from seqkit) ----
  # could pre process in command line if too slow here (by converting space into \t?)
  message("Reading in barcode key...")
  barcode_key1_raw = fread(barcode_key, sep="\t", header = F)
  barcode_key1 = barcode_key1_raw
  barcode_key1 = barcode_key1 %>% select(V1)
  barcode_key1$id = str_extract(barcode_key1$V1, "^[\\w-]+(?=\\s)")
  barcode_key1$barcode = str_extract(barcode_key1$V1, "(?<=barcode\\=)[A-z0-9=]+")
  barcode_key1$start_time = str_extract(barcode_key1$V1, "(?<=start_time\\=).+(?=\\s)")
  barcode_key1 = barcode_key1 %>% 
    select(-V1)
  
  # # check dupes
  # barcode_key1 %>% group_by(id) %>% filter(n()>1)
  # barcode_key1 %>% count(barcode)
  
  rm(barcode_key1_raw)
  
  
  # ---------------------------------------------------------------------------- #
  
  
  # get n reads per barcode w/ con seq ----
  barcodes_per_cluster = 
    barcode_key1 %>% 
    left_join(cluster_key1, by = c("id" = "query_id")) %>%
    group_by(con_id, barcode) %>%
    summarise(n_reads_barcode = n(), .groups = 'drop') %>% 
    ungroup() %>%
    ## join in con seq (only clusters that got blasted will have seqs)
    left_join(conseq_key, by = c("con_id"="id")) %>% 
    mutate(seq_len = nchar(sequence))
  
  
  
  # ---------------------------------------------------------------------------- #
  
  # load in taxononmy keys ----
  taxa_key = read.csv(taxdump_filepath)
  taxa_key$tax_id = as.character(taxa_key$tax_id)
  if(is.character(custom_taxa_filepath)){
    message("Reading in custom_taxa_filepath...")
    custom_taxa_key = read.csv(custom_taxa_filepath, na.strings = c("NA"))
    custom_taxa_key$taxid = as.character(custom_taxa_key$taxid)
  }else{
    custom_taxa_key = data.frame()
  }
  
  # load in sample key ----
  message("Reading in sample key...")
  sample_key_df = read.xlsx(sample_key)
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

  
  # ---------------------------------------------------------------------------- #
  
  # load in blast data & pre process ----
  message("Reading in blast results...")
  blastn_out = read.csv(blast_results, header = F, na.strings = c("N/A", "NA"))
  colnames(blastn_out) = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                           'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen',
                           'staxids', 'sscinames', 'scomnames', 'sskingdoms')
  blastn_out$staxids = as.character(blastn_out$staxids)
  
  ## get cluster size, per cluster that was blasted
  blastn_out2 = 
    blastn_out %>% 
    separate(qseqid, sep = ";", into = c("id", "size")) %>% 
    mutate(id = str_extract(id, "(?<==).+")) %>%
    mutate(size = as.integer(str_extract(size, "(?<==)[0-9]+"))) %>% 
    filter(!is.na(id))
  
  
  ## add in ncbi taxonomy and rank + custom taxonomy if not in ncbi or if blast didn't assign
  blastn_out2t = add_taxonomy(blastn_out2, taxa_key, custom_taxa_key)
  
  
  # ---------------------------------------------------------------------------- #
 
  # summarize blast output ----
  
  # summarize by cluster
  
  ## summarize blast output at >98% pident and 95-98% pident (per cluster)
  blast_sum98 = blast_summarize_by_cluster(blastn_out2t, barcodes_per_cluster, taxa_key, min_alignment_length, 98, 100)
  blast_sum95 = blast_summarize_by_cluster(
                    blastn_out2t %>% filter(!id %in% blast_sum98$id), 
                    barcodes_per_cluster, taxa_key, min_alignment_length, 95, 97.999999) 
  
  ##  95% pident df to change species matches to genus
  if(nrow(blast_sum95) > 0){
    blast_sum95 =
      blast_sum95 %>%
      mutate(lca_sciname = ifelse(lca_rank == "species", str_extract(lca_sciname, "^[A-z]+(?=\\s)"), lca_sciname),
             lca_rank = ifelse(lca_rank == "species", "genus", lca_rank)
      )
  }
  
  ## combine 95 and 98 files for export & join in sample names
  blast_sum_all = 
    bind_rows(blast_sum98 %>% mutate(blast_threshold = "98", .before = barcode) , 
              blast_sum95 %>% mutate(blast_threshold = "95", .before = barcode)
              ) %>% 
    mutate(blast_refdb = local_global, .after = blast_threshold) %>%
    mutate(primer = primer_name, .before = blast_threshold) %>% 
    left_join(sample_barcode_key_primer %>% 
                select(barcode, Sample.ID, Sample.Type), 
              by = c("barcode" = "barcode")) %>% 
    filter(!is.na(id))
  

 
   # ---------------------------------------------------------------------------- #
  
  # export summaries ----
  
  # by cluster
  write.csv(x = blast_sum_all, file = output_cluster_filepath, row.names = F)
  
  # # by taxa
  # write.csv(x = sum_taxa, file = output_taxa_filepath, row.names = F)
  
  
  # Export IDs of clusters that had any blast matches
  
  # if clusters_matched_filepath is not provided, do not export
  # if clusters_matched_filepath is a string, export
  if(is.character(clusters_matched_filepath)){
    
    ## get clusters that did have any blast matches
    clusters_matched = bind_rows(blast_sum98, blast_sum95) %>% 
      distinct(id) %>% 
      pull(id)
    
    write.table(clusters_matched, file = clusters_matched_filepath,
                row.names = F, col.names = F, quote = F)
  }
  
  
  # ---------------------------------------------------------------------------- #

  
} # end pipeline_v2


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# function to add in ncbi taxonomy and rank + custom taxonomy if not in ncbi or if blast didn't assign
add_taxonomy = function(blastn_out2, taxa_key, custom_taxa_key){
  
  # create new custom taxa key based on sseqid, but use taxid to get taxonomy from ncbi
  
  if(nrow(custom_taxa_key) > 0){
    ## for all taxa in blast_out
    custom_taxa_key2a = 
      custom_taxa_key %>% 
      select(-any_of(c("superkingdom", "phylum", "class", "order", "family"))) %>% # rm's superkingdom + family in ref db csv if exists, bc also exists in ncbi taxa dump
      filter(seqID %in% blastn_out2$sseqid) %>% 
      left_join(taxa_key %>% select(tax_id:species, subspecies), by = c("taxid" = "tax_id")) %>% 
      mutate(rank = case_when(str_detect(species.x, "\\[") ~ "genus",
                              TRUE ~ "species")) %>% 
      left_join(blastn_out2 %>% select(staxids, scomnames) %>% distinct(), by = c("taxid" = "staxids")) %>% 
      mutate(sscinames = species.x) %>% 
      mutate(species = case_when(str_detect(species.x, "\\[") ~ NA_character_,
                                 TRUE ~ species.y)) %>% 
      # no longer need taxid bc we have taxonomy and we use seqID for join
      select(seqID, taxid, rank, sscinames, scomnames, superkingdom:family, genus.x, species) %>% 
      rename(genus = genus.x)
    
    ## repeat except for species/taxa not in NCBI (taxid == NA)
    custom_taxa_key2b = 
      custom_taxa_key2a %>% 
      select(seqID, taxid, rank, sscinames, scomnames, genus) %>%
      filter(is.na(taxid)) %>% 
      left_join(taxa_key %>% filter(genus %in% .$genus) %>% 
                  select(superkingdom:genus) %>% distinct(), 
                by = c("genus" = "genus")) %>% 
      mutate(species = ifelse(rank == "species", sscinames, NA_character_))
    
    ## combine into one key
    if(nrow(custom_taxa_key2b) > 0 ){
      custom_taxa_key2 = 
        bind_rows(custom_taxa_key2a %>% filter(!is.na(taxid)), 
                  custom_taxa_key2b)
    }else{
      custom_taxa_key2 = custom_taxa_key2a
    }
    
  # if no custom taxa key provided then make empty df in its place
  }else{
    custom_taxa_key2 = data.frame(seqID = character())
  }
  
  
  # now make a key for all blast taxids that are not in the custom key (why they are not idk)
  taxa_key_subset =
    taxa_key %>%
    filter(tax_id %in% blastn_out2$staxids) %>%
    # filter(!tax_id %in% custom_taxa_key2$taxid) %>%
    select(tax_id:species, subspecies) %>%
    # make blanks na
    mutate(across(c(superkingdom:subspecies), ~na_if(., ""))) %>%
    # add rank
    mutate(rank = case_when(!is.na(subspecies) ~ "subspecies",
                            !is.na(species) ~ "species",
                            !is.na(genus) ~ "genus",
                            !is.na(family) ~ "family",
                            !is.na(order) ~ "order",
                            !is.na(class) ~ "class",
                            !is.na(phylum) ~ "phylum",
                            TRUE ~ NA_character_)) %>%
    left_join(blastn_out2 %>% select(staxids, scomnames) %>% distinct(), by = c("tax_id" = "staxids")) %>%
    mutate(sscinames = case_when(!is.na(subspecies) ~ subspecies,
                                 !is.na(species) ~ species,
                                 !is.na(genus) ~ genus,
                                 !is.na(family) ~ family,
                                 !is.na(order) ~ order,
                                 !is.na(class) ~ class,
                                 !is.na(phylum) ~ phylum))
  
  
  
  
  
  # join (have to do in 2 steps, one for each key)
  blastn_out2t1 = 
    blastn_out2 %>% 
    select(-scomnames, -sscinames, -sskingdoms) %>% 
    filter(sseqid %in% custom_taxa_key2$seqID) %>%
    left_join(custom_taxa_key2 %>% distinct(), by = c("sseqid" = "seqID"))
  
  blastn_out2t2 = 
    blastn_out2 %>% 
    select(-scomnames, -sscinames, -sskingdoms) %>% 
    filter(!sseqid %in% custom_taxa_key2$seqID) %>%
    left_join(taxa_key_subset, by = c("staxids" = "tax_id"))

  # set dtype sscinames of to character to avoid issues with dplyr::
  blastn_out2t1$sscinames = as.character(blastn_out2t1$sscinames)
  blastn_out2t2$sscinames = as.character(blastn_out2t2$sscinames)

  
  
  blastn_out2t = bind_rows(blastn_out2t1, blastn_out2t2)
  
  return(blastn_out2t)
  
}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# function to summarize blast output by pident threshold and min alignment length

# First fun is to summarize by cluster/barcode
blast_summarize_by_cluster = function(blastn_out2t, 
                                barcodes_per_cluster, taxa_key,
                                min_alignment_len, 
                                min_pident_threshold, max_pident_threshold){
  
  

  
  # filter by min. alignment length & by pident
  blastn_out2a =
    blastn_out2t %>%
    filter(as.numeric(`length`) >= min_alignment_len, # to re-visit
           pident >= min_pident_threshold & pident <= max_pident_threshold) %>% 
    mutate(row = row_number(), .before=id)
  
  
  # if none remain
  if(nrow(blastn_out2a) == 0){
    return(data.frame(barcode = NA_character_))
  }

  
  # summarize by cluster/barcode combo (at the species level - ie no subspecies)
  cluster_taxa_combos =
    blastn_out2a %>%
    # change global blast output names to match local - only needed to do when changed to Oliver's ref db format ... not sure why
    {if ( "staxids" %in% colnames(.) & !("taxid" %in% colnames(.)) ) mutate(., taxid = staxids) else .} %>% 
    group_by(id) %>%
    # filter(n_distinct(species, na.rm = T) > 1) %>%
    ungroup() %>% 
    # temporarily ovderride sciname for subspecies so i can summarize by species
    mutate(sscinames = ifelse(!is.na(subspecies), species, sscinames)) %>%
    group_by(id, sscinames) %>% 
    arrange(-pident) %>% 
    slice(1) %>% # should take the highest pident for each id-sciname combo
    select(row, id, pident, length, qstart, qend, sseqid, sscinames_sp = sscinames, scomnames, taxid, species) %>% 
    ungroup() %>% 
    # select(-sscinames) %>% 
    # join back in correct ssciname
    left_join(blastn_out2a %>% select(row, sscinames) %>% distinct(), by = "row") %>%
    group_by(id) %>% 
    # arrange(-pident) %>%
    arrange(sscinames) %>% 
    summarise(blast_n_species = n_distinct(sscinames, na.rm = T),
              sscinames_sp = paste0(sscinames_sp, collapse = "; "),
              sscinames = paste0(unique(sscinames), collapse = "; "),
              scomnames = paste0(unique(scomnames), collapse = "; "),
              max_pident = paste0(pident, collapse = "; "),
              align_len = paste0(length, collapse = "; "),
              qstart = paste0(qstart, collapse = "; "),
              qend = paste0(qend, collapse = "; "),
              sseqid = paste0(sseqid, collapse = "; "),
              taxid = paste0(taxid, collapse = "; "),
              .groups = 'drop') %>% 
    ## join in barcode info, incl conseq
    ## note will increase rows bc some clusters are across multiple barcodes
    left_join(barcodes_per_cluster, by = c("id"="con_id"))
  
  
  # # summarize by cluster/barcode combo
  # sum_blast = 
  #   ## summarize blast output by cluster id
  #   blastn_out2a %>% 
  #   group_by(id) %>%
  #   arrange(sscinames) %>% 
  #   summarize(scinames = paste(unique(sscinames), collapse = "; "),
  #             comnames = paste(unique(scomnames), collapse = "; "),
  #             n_taxa = n_distinct(sscinames),
  #             pident = paste(unique(as.numeric(pident)), collapse = "; "), 
  #             n_reads_cluster = unique(size),
  #             qstart = paste(unique(qstart), collapse = "; "), 
  #             qend = paste(unique(qend), collapse = "; "),
  #             avg_alignment_length = mean(length) %>% round(0),
  #             sseqids = paste(unique(sseqid), collapse = "; "),
  #             staxids = paste(unique(staxids), collapse = "; ")
  #             ) %>% 
  #   ungroup() %>% 
  #   ## join in barcode info, including consensus sequence + some clusters are across multiple barcodes
  #   left_join(barcodes_per_cluster, by = c("id"="con_id")) %>% 
  #   arrange(barcode, -n_reads_barcode, scinames, .by_group = TRUE) %>% 
  #   rename(consensus_sequence = sequence) %>% 
  #   select(barcode, id:pident, n_reads_barcode, n_reads_cluster, 
  #          avg_alignment_length,
  #          con_seq_len = seq_len, consensus_sequence, everything()) 
  
  # run lca
  sum_blast2 = add_lca_to_blast_output(blastn_out2a, cluster_taxa_combos) %>% 
    rename(lca_sciname = lca, 
           lca_rank = rank, 
           blast_scinames = sscinames, 
           blast_comnames = scomnames) %>%
    select(-seq_len) %>% 
    select(barcode, id, 
           blast_scinames, blast_comnames, blast_n_species, 
           superkingdom:subspecies, 
           lca_sciname, lca_rank,
           max_pident,
           n_reads_barcode, # n_reads_cluster, 
           align_len, # con_seq_len, 
           qstart, qend, 
           sequence, 
           everything()) %>% 
    arrange(-n_reads_barcode) %>% 
    # if species name has sp. then it should be genus and give lca_sciname genus name (and sometimes family)
    mutate(lca_rank = case_when(((lca_rank == "species") & str_detect(lca_sciname, "\\bsp\\.")) ~ "genus",
                                TRUE ~ lca_rank)) %>% 
    mutate(lca_sciname = case_when(str_detect(lca_sciname, "\\bsp\\.") & !is.na(genus) ~ genus,
                                   str_detect(lca_sciname, "\\bsp\\.") & !is.na(family) ~ family,
                                   TRUE ~ lca_sciname))
  
  return(sum_blast2)
  
}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# function to add on LCA to blast output
# input params: blast summary by barcode (df), csv of taxdump
add_lca_to_blast_output = function(blastn_out2a, cluster_taxa_combos){
  
  
  # join back in taxa info
  df2 = 
    cluster_taxa_combos %>% 
    mutate(row = row_number()) %>%
    separate_rows(sseqid, sep = "; ") %>%
    left_join(blastn_out2a %>% select(sseqid, rank, superkingdom:subspecies) %>% distinct(), 
              by = c("sseqid" = "sseqid"))
  
  # sum(cluster_taxa_combos$id %in% df2$id) == nrow(cluster_taxa_combos)
  
  
  # sum n taxa by ranks by row, to set up LCA
  df3 = 
    df2 %>% 
    group_by(row) %>% 
    summarise(n_na_species = sum(is.na(species)),
              n_species = n_distinct(species, na.rm = T),
              n_genus = n_distinct(genus, na.rm = T),
              n_family = n_distinct(family, na.rm = T),
              n_order = n_distinct(order, na.rm = T),
              n_class = n_distinct(class, na.rm = T),
              n_phylum = n_distinct(phylum, na.rm = T)) %>% 
    left_join(df2 %>% select(row, superkingdom:subspecies) %>% group_by(row) %>% slice(1), 
              by = "row") %>% 
    ungroup()
  
  
  
  # Run logic for LCA
  df4 = 
    df3 %>% 
    mutate(lca = case_when(
      # when one+ genus-level [] seq matches with one+ species of different order(s), same class
      n_order >= 2 & n_na_species >= 1 & n_class == 1 ~ class,
      # when one+ genus-level [] seq matches with one+ species of different familie(s), same order
      n_family >= 2 & n_na_species >= 1 & n_order == 1 ~ order,
      # when one+ genus-level [] seq matches with one+ species of different genus, same family
      n_genus >= 2 & n_na_species >= 1 & n_family == 1 ~ family,
      # when one genus level [] seq matches with one species of same genus
      n_genus == 1 & n_na_species >= 1 ~ genus,
      # normal matching
      n_species == 1 ~ species,
      n_genus == 1 ~ genus,
      n_family == 1 ~ family,
      n_order == 1 ~ order,
      n_class == 1 ~ class,
      n_phylum == 1 ~ phylum,
      TRUE ~ NA_character_
    )) %>% 
    mutate(rank = case_when(
      n_order >= 2 & n_na_species >= 1 & n_class == 1 ~ "class",
      n_family >= 2 & n_na_species >= 1 & n_order == 1 ~ "order",
      n_genus >= 2 & n_na_species >= 1 & n_family == 1 ~ "family",
      n_genus == 1 & n_na_species >= 1 ~ "genus",
      n_species == 1 ~ "species",
      n_genus == 1 ~ "genus",
      n_family == 1 ~ "family",
      n_order == 1 ~ "order",
      n_class == 1 ~ "class",
      n_phylum == 1 ~ "phylum",
      TRUE ~ NA_character_
    )) %>%
    # fix taxonomy to match lca taxonomy (ie add in NA to uppuer taxonomy, if rank is higher than species)
    mutate(class = ifelse(rank %in% c("phylum"), NA_character_, class),
           order = ifelse(rank %in% c("phylum", "class"), NA_character_, order),
           family = ifelse(rank %in% c("phylum", "class", "order"), NA_character_, family),
           genus = ifelse(rank %in% c("phylum", "class", "order", "family"), NA_character_, genus),
           species = ifelse(rank %in% c("phylum", "class", "order", "family", "genus"), NA_character_, species),
           subspecies = ifelse(rank %in% c("phylum", "class", "order", "family", "genus", "species"), NA_character_, subspecies)) %>%
    # rm no lca matches
    # filter(!is.na(lca)) %>%
    select(row, lca, rank, superkingdom:subspecies)  
  
  
  # add lca results back in the original data 
  df5 = 
    cluster_taxa_combos %>% 
    mutate(row = row_number()) %>% 
    left_join(df4, by = "row") %>% 
    select(-row) %>% 
    select(barcode, id, lca:subspecies, n_reads_barcode, everything()) 
  
  
  return(df5)
  
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





# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# Call the function with named parameters ----
pipeline_v2(cluster_key = args$cluster_key, 
            barcode_key = args$barcode_key, 
            blast_results = args$blast_results, 
            conseq_key = args$conseq_key,
            sample_key = args$sample_key,
            primer_name = args$primer_name,
            min_alignment_length = args$min_alignment_length,
            local_global = args$local_global,
            output_cluster_filepath = args$output_cluster_filepath,
            output_taxa_filepath = args$output_taxa_filepath,
            clusters_matched_filepath = args$clusters_matched_filepath,
            taxdump_filepath = args$taxdump_filepath,
            custom_taxa_filepath = args$custom_taxa_filepath)