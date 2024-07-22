parse_data_summary <- function(path){
  # This parses the data summary file which lists how many of each organisms reads went into each file, and it adds in the 0 count rows for those files without spikes
  # this makes merging much eaasier fown the road
  spikes <- c("Haloarcula_hispanica", "Salinibacter_ruber", "Trichoderma_reesei")
  spikes_pri <- paste0(spikes, "_primary") 
  
  data_summary_most <- read.csv(path, sep="\t") %>%
    #filter(!grepl("spike0.00000", file)) %>% 
    mutate(
      true_count=as.numeric(gsub(",", "", num_seqs)),
      sample=gsub(".*?_(\\d+_depth.*?_spike.*?)_.*", "\\1", file),
      sample=gsub(".*?_(\\d+_evencoverage.*?)_.*", "\\1", sample),
      Genome = gsub(".*?_(\\d+_depth.*?_spike.*?)_(.*).fasta.*", "\\2", file),
      Genome = gsub(".*?_(\\d+_evencoverage.*?)_(.*).fasta.*", "\\2", Genome),
      is_spike =Genome %in% c(spikes, spikes_pri)
    ) %>% 
    filter(!grepl("fastq.gz", sample)) %>% # these are the sample totals, we already calculate these
    group_by(sample,is_spike, Genome) %>% 
    summarize(true_count = sum(true_count)) %>%
    ungroup() %>%
    group_by(sample) %>% 
    mutate(total = sum(true_count)) %>%
    ungroup() %>% 
    #  filter( is_spike)%>% 
    select(sample, Genome, is_spike, true_count, total) %>% 
    distinct()
  # parse the ones lacking any spikes so we still have totals
  # the left join with the dummy data makes the merging below more straightforward
  data_summary_nospike <- read.csv(path, sep="\t") %>% 
    filter(grepl("spike0.00000", file)) %>% 
    mutate(
      true_count=as.numeric(gsub(",", "", num_seqs)),
      sample=gsub(".*?_(\\d+_depth.*?_spike.*?)_.*", "\\1", file),
      sample=gsub(".*?_(\\d+_evencoverage.*?)_.*", "\\1", sample),
    ) %>% 
    filter(!grepl("fastq.gz", sample)) %>% # these are the sample totals, we already calculate these
    group_by(sample) %>% 
    summarize(total = sum(true_count)) %>% ungroup() %>%
    mutate(true_count=0) %>% 
    left_join(data.frame(true_count=rep(0, length(spikes)), Genome=spikes), relationship= "many-to-many", by="true_count") %>%
    mutate(is_spike =Genome %in% c(spikes, spikes_pri))
  
  
  return(bind_rows(data_summary_most,data_summary_nospike) %>% 
           rename("genome_base" = "Genome"))
}
