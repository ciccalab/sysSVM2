# This script contains functions to annotate standard file formats for use with sysSVM2
# Small somatic mutations (SSMs) can be annotated with annotate_ssms()
# Copy number variants (CNVs) can be annotated with annotate_cnvs()
# Requirements:
#  - annotate_ssms() requires ANNOVAR to be installed
#  - annotate_cnvs() requires bedtoools to be installed
#  - Both functions require data files from the sysSVM2 GitHub repository (gene_coords_hgXX.tsv, hotspots_tcga.tsv, gene_aliases.tsv, XXXXX)!!!!



# Function to map hotspots to annotated SSMs
# Used by annotate_ssms
map_hotspots = function(
  ssms,
  hotspots
){
  
  # Only consider missense SNVs
  ssms_orig = ssms
  ssms = ssms %>% 
    subset(ExonicFunc.refGene == "nonsynonymous SNV") %>%
    select(Chr:Alt, symbol, AAChange.refGene) %>%
    mutate(AAChange.refGene = as.character(AAChange.refGene))
  
  
  # Only consider hotspots that are significant and in genes where there are SNVs
  hotspots = hotspots %>% select(GENE, CLUST_COORDS, QVALUE) %>% subset(QVALUE <= 0.1 & GENE %in% ssms$symbol)
  
  
  # Exit if there aren't any matches at this point
  if (nrow(hotspots) == 0){
    ssm_hotspots = ssms_orig %>% mutate(is_hotspot = F)
    return(ssm_hotspots)
    next
  }
  
  
  # Function to extract the amino acid position of mutations
  parse_AAChange = function(x){
    
    # Separate different transcripts
    x = strsplit(x, split = ",") %>% unlist
    
    # Get the amino acid position in each transcript
    x = sapply(x, function(y){
      strsplit(y, split = ":") %>%
        unlist %>%
        subset(grepl("p.", ., fixed = T)) %>%
        gsub("p.", "", ., fixed = T) %>%
        gsub("[^0-9.-]", "", .) %>%
        as.numeric
    })
    
    # Concatenate the unique positions
    x = unique(x) %>%
      as.numeric %>%
      paste(collapse = ";")
    return(x)
  }
  
  
  # Annotate the mutated AA positions for each missense SNV
  ssms$AA_position = sapply(ssms$AAChange.refGene, parse_AAChange, USE.NAMES = F)
  
  
  # Function to annotate if missense SNVs are in a given hotspot
  map_one_hotspot = function(oncodrive_row){
    
    # Only consider SNVs in this gene
    df = ssms %>% subset(symbol == oncodrive_row["GENE"])
    if (nrow(df) == 0) next
    
    # Initialise hotspot column
    df$is_hotspot = F
    
    # Extract hotspots in this gene
    hotspots_thisGene = strsplit(oncodrive_row["CLUST_COORDS"], split = ",[", fixed = T)[[1]] %>%
      gsub("[", "", ., fixed = T) %>%
      paste0("[", .)
    
    # Loop through hotspots and match to coordinates of SNVs
    for (hotspot in hotspots_thisGene){
      # Extract AA coordinates of this hotspot
      aa_coords = strsplit(hotspot, split = ":")[[1]][1] %>%
        strsplit(., split = ",") %>%
        unlist %>%
        gsub("[", "", ., fixed = T) %>%
        gsub("]", "", ., fixed = T) %>%
        as.numeric
      
      # Match to AA coordinates of SNVs in this gene
      is_hotspot = lapply(strsplit(df$AA_position, split = ";"), function(x){
        any(between(as.numeric(x), aa_coords[1], aa_coords[2]))
      }) %>%
        unlist
      
      # Update hotspot column
      df$is_hotspot = (df$is_hotspot == T) | is_hotspot
    }
    return(df)
  }
  
  
  # Match hotspots with SNVs
  ssm_hotspots = apply(hotspots, 1, map_one_hotspot)
  ssm_hotspots = do.call(rbind, ssm_hotspots)

  
  # Join back to full SSM table
  ssm_hotspots = ssm_hotspots %>% select(Chr:Alt, is_hotspot)
  ssm_hotspots = ssms_orig %>% left_join(ssm_hotspots, by = c("Chr", "Start", "End", "Ref", "Alt"))
  ssm_hotspots$is_hotspot[is.na(ssm_hotspots$is_hotspot)] = F
  return(ssm_hotspots)
}



# Function to annotate small somatic mutations from a VCF
# gene_aliases_entrez and hotspots can be passed as data frames or file names
annotate_ssms = function(
  vcf_fn,                    # File name of somatic VCF to annotate
  sample,                    # Name of this sample (NB can't process multi-sample VCFs)
  annovar_dir,               # Directory where ANNOVAR is installed (should contain table_annovar.pl)
  genome_version = "hg19",   # Version of the human genome to use for ANNOVAR - hg19 or hg38
  gene_aliases_entrez,       # gene_aliases_entrez.tsv from the sysSVM2 GitHub repository
  hotspots,                  # tcga_pancancer_hotspots_oncodriveclust.tsv from the sysSVM2 GitHub repository
  temp_dir = tempdir()       # Directory for temporary files to be created
){
  
  # Packages
  require(readr)
  require(tidyr)
  require(dplyr)
  
  
  # Run ANNOVAR on VCF
  annovar_output_fn = tempfile(tmpdir = temp_dir)
  annovar_cmd = paste0(
    "perl ", annovar_dir, "/table_annovar.pl -vcfinput ", vcf_fn, 
    " ", annovar_dir, "/humandb/ -buildver ", genome_version, " -out ", annovar_output_fn, " -remove ", 
    "-protocol refGene,dbnsfp35a,dbscsnv11 ", 
    "-operation g,f,f"
    )
  message("Running ANNOVAR...")
  system(annovar_cmd)
  message("Done")
  
  
  # Load annotated file and delete temporary files
  annovar_output = suppressWarnings(read_tsv(paste0(annovar_output_fn, ".", genome_version, "_multianno.txt"), col_types = cols()))
  cleanup_cmd = paste0("rm ", annovar_output_fn, "*")
  system(cleanup_cmd)
  
  
  # Clean up the annovar output
  if ("GERP++_RS" %in% colnames(annovar_output)){
    annovar_clean = annovar_output %>% rename(GERPpp_RS = `GERP++_RS`)
  } else if ("GERP.._RS" %in% colnames(annovar_output)){
    annovar_clean = annovar_output %>% rename(GERPpp_RS = `GERP.._RS`)
  } else {
    stop("GERP++ annotations missing")
  }
  annovar_clean = annovar_clean %>%
    select(
      Chr:Alt, 
      Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene,
      SIFT_pred, LRT_pred, FATHMM_pred, MutationTaster_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, MutationAssessor_pred,
      phyloP100way_vertebrate, GERPpp_RS, SiPhy_29way_logOdds,
      dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE
    ) %>%
    mutate_at(
      c("phyloP100way_vertebrate", "GERPpp_RS", "SiPhy_29way_logOdds", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE"),
      function(x) suppressWarnings(as.numeric(x))
    ) %>%
    mutate(Gene.refGene = strsplit(Gene.refGene, split = ";")) %>%
    tidyr::unnest(Gene.refGene) %>%
    unique
  
  
  # Get exonic/splicing variants within a consistently-named set of 19549 human genes
  annovar_exonic = annovar_clean %>% 
    subset(Func.refGene == "exonic" | Func.refGene %in% c("splicing", "exonic;splicing"))
  if (is.character(gene_aliases_entrez)) gene_aliases_entrez = read_tsv(gene_aliases_entrez, col_types = cols())
  gene_aliases_entrez = gene_aliases_entrez %>% select(symbol, alias, entrez)
  annovar_exonic = inner_join(
    annovar_exonic,
    gene_aliases_entrez,
    by = c("Gene.refGene" = "alias")
  )
  
  
  # Parse truncating and predicted damaging missense/splicing mutations
  damaging_trunc_ntdam = annovar_exonic %>%
    # Truncating
    mutate(dam_trunc = ExonicFunc.refGene %in% c("stopgain", "stoploss", "frameshift deletion", "frameshift insertion", "frameshift substitution")) %>%
    mutate(dam_trunc = replace_na(dam_trunc, F)) %>%
    # Functional 5/7
    mutate(
      dam_sift = SIFT_pred == "D",
      dam_lrt = LRT_pred == "D",
      dam_fathmm = FATHMM_pred == "D",
      dam_mutationtaster = MutationTaster_pred %in% c("D", "A"),
      dam_polyphen_hdiv = Polyphen2_HDIV_pred %in% c("D", "P"),
      dam_polyphen_hvar = Polyphen2_HVAR_pred %in% c("D", "P"),
      dam_mutationassessor = MutationAssessor_pred %in% c("H", "M")
    ) %>%
    mutate_at(c("dam_sift", "dam_lrt", "dam_fathmm", "dam_mutationtaster", "dam_polyphen_hdiv", "dam_polyphen_hvar", "dam_mutationassessor"),
              function(x) replace_na(x, F)) %>%
    mutate(dam_func = !is.na(ExonicFunc.refGene) & ExonicFunc.refGene == "nonsynonymous SNV" &
             dam_sift + dam_lrt + dam_fathmm + dam_mutationtaster + dam_polyphen_hdiv + dam_polyphen_hvar + dam_mutationassessor >= 5) %>%
    # Conservation 2/3
    mutate(
      dam_phylop = phyloP100way_vertebrate > 1.6,
      dam_gerp = GERPpp_RS > 4.4,
      dam_siphy = SiPhy_29way_logOdds > 12.17
    ) %>%
    mutate_at(c("dam_phylop", "dam_gerp", "dam_siphy"), function(x) replace_na(x, F)) %>%
    mutate(dam_cons = !is.na(ExonicFunc.refGene) & ExonicFunc.refGene == "nonsynonymous SNV" &
             dam_phylop + dam_gerp + dam_siphy >= 2) %>%
    # Splicing 1/2
    mutate(
      dam_ada = dbscSNV_ADA_SCORE > 0.6,
      dam_rf = dbscSNV_RF_SCORE > 0.6
    ) %>%
    mutate_at(c("dam_ada", "dam_rf"), function(x) replace_na(x, F)) %>%
    mutate(dam_splicing = dam_ada + dam_rf >= 1 & Func.refGene %in% c("exonic;splicing", "splicing"))
  
  
  # Map hotpots
  if (is.character(hotspots)) hotspots = read_tsv(hotspots, col_types = cols())
  damaging_trunc_ntdam_gof = map_hotspots(damaging_trunc_ntdam, hotspots)

  
  # Clean table
  damaging_ssms = damaging_trunc_ntdam_gof %>%
    unite(
      "variant_id", c(Chr, Start, End, Ref, Alt), sep = ";"
    ) %>%
    mutate(
      TRUNC_mut = dam_trunc,
      NTDam_mut = dam_func | dam_cons | dam_splicing,
      GOF_mut = is_hotspot,
      damaging = TRUNC_mut | NTDam_mut | GOF_mut,
      sample = sample
    ) %>%
    select(
      sample, variant_id, symbol, entrez, Func.refGene, ExonicFunc.refGene,
      TRUNC_mut, NTDam_mut, GOF_mut, damaging
    )
  
  
  # Output
  return(damaging_ssms)
}



# Function to annotate copy number variants (CNVs)
# Input files can be passed as data frames or file names
annotate_cnvs = function(
  # Arguments
  cnv_segments,            # Table with the following columns: sample; chromosome; start; end; and copy_number or segment_mean
  ploidy = NULL,           # Table with the following columns: sample; ploidy. Leave null if unavailable (assumes diploidy)
  ploidy_threshold = 2,    # Threshold for determining gene amplifications: CN >= ploidy_threshold*ploidy
  gene_coords,             # gene_coords_hg19.tsv or gene_coords_hg38.tsv from the sysSVM2 GitHub repository
  bedtools_bin_dir = NULL, # Directory where the bedtools binary executable is located, if not in $PATH
  temp_dir = tempdir()     # Directory for temporary files to be created
  ){
  
  
  # Packages
  require(readr)
  require(dplyr)
  
  
  # Load data if required
  if (is.character(cnv_segments)) cnv_segments = read_tsv(cnv_segments, col_types = cols())
  if (is.null(ploidy)){
    warning(paste0("No ploidy values provided, assuming diploidy (amplifications will have CN >= ", 2*ploidy_threshold, ")"))
    ploidy = data.frame(sample = unique(cnv_segments$sample), ploidy = 2, stringsAsFactors = F) # Assume diploidy in the absence of other data
  } else if (is.character(ploidy)){
    ploidy = read_tsv(ploidy, col_types = cols())
  }
  if (is.character(gene_coords)) gene_coords = read_tsv(gene_coords, col_types = cols())
  
  
  # Make sure we have the correct columns
  if (!("copy_number" %in% colnames(cnv_segments))){ # Can provide either segment mean or copy number (preferentially copy number)
    if ("segment_mean" %in% colnames(cnv_segments)){
      cnv_segments$copy_number = round(2^cnv_segments$segment_mean * 2)
    } else {
      stop("Must provide either a copy_number or segment_mean column in cnv_segments")
    }
  }
  cnv_segments = cnv_segments %>% select(sample, chromosome, start, end, copy_number)
  ploidy = ploidy %>% select(sample, ploidy)
  gene_coords = gene_coords %>% select(entrez, symbol, chromosome, tstart, tend)
  
  
  # Make sure cnv_segments chromosome names start with "chr" for consistency
  missing_chr_ind = !grepl("^chr", cnv_segments$chromosome)
  cnv_segments$chromosome[missing_chr_ind] = paste0("chr", cnv_segments$chromosome[missing_chr_ind])
  
  
  # Make temporary bed files
  # CNV segments
  cnv_segments_bed = cnv_segments %>% select(chromosome, start, end, sample)
  cnv_segments_bed_fn = tempfile(tmpdir = temp_dir, fileext = ".bed")
  write_tsv(cnv_segments_bed, path = cnv_segments_bed_fn, col_names = F)
  # Gene coordinates
  gene_coords_bed = gene_coords %>% select(chromosome, tstart, tend)
  gene_coords_bed_fn = tempfile(tmpdir = temp_dir, fileext = ".bed")
  write_tsv(gene_coords_bed, path = gene_coords_bed_fn, col_names = F)
  
  
  # Run bedtools intersect
  if (!is.null(bedtools_bin_dir)) bedtools_bin_dir = paste0(bedtools_bin_dir, "/")
  intersected_bed_fn = tempfile(tmpdir = temp_dir, fileext = ".bed")
  bedtools_command = paste0(bedtools_bin_dir, "bedtools intersect -wo -a ", cnv_segments_bed_fn, " -b ",  gene_coords_bed_fn, " > ", intersected_bed_fn)
  message("Running bedtools...")
  system(bedtools_command)
  message("Done")
  
  
  # Read results and remove temporary files
  bedtools_res = read_tsv(intersected_bed_fn, col_names = c("cnv_chromosome", "cnv_start", "cnv_end", "sample", "gene_chromosome", "gene_start", "gene_end", "overlap"), col_types = cols())
  cleanup_command = paste("rm", cnv_segments_bed_fn, gene_coords_bed_fn, intersected_bed_fn)
  system(cleanup_command)
  
  
  # Annotate genes and copy numbers
  cnv_genes = bedtools_res %>%
    left_join(cnv_segments, by = c("sample", "cnv_chromosome" = "chromosome", "cnv_start" = "start", "cnv_end" = "end")) %>%
    left_join(gene_coords, by = c("gene_chromosome" = "chromosome", "gene_start" = "tstart", "gene_end" = "tend")) %>%
    mutate(gene_length = gene_end - gene_start, overlap_percent = overlap / gene_length * 100)
  
  
  # Subset for damaging CNVs:
  #  - Overlap a gene with >=25% coverage
  #  - Either CN<=1 or CN>=ploidy_threshold*ploidy
  damaging_cnvs = cnv_genes %>%
    subset(overlap_percent >= 25) %>%
    left_join(ploidy, by = "sample") %>%
    mutate(
      CNVGain = case_when(copy_number >= ploidy_threshold*ploidy ~ 1, T ~ 0),
      CNVLoss = case_when(copy_number <= 1 ~ 1, T ~ 0)
      ) %>%
    subset(CNVGain == 1 | CNVLoss == 1)
  
  
  # If multiple segments overlap a gene, take the CN from the one with the largest overlap
  damaging_cnvs_noDup = damaging_cnvs %>%
    group_by(sample, entrez) %>%
    top_n(1, overlap_percent) %>%
    slice(1) %>% # If there is an exact tie, just take the first one
    ungroup
  
  
  # Prettify output
  # NB will often have an unreasonable amount of CN=1 variants, but these are only retained later if they have a damaging mutation
  damaging_cnvs_noDup = damaging_cnvs_noDup %>%
    select(
      sample, entrez, symbol,
      Copy_number = copy_number, ploidy, CNVGain, CNVLoss
    )
  message("CNVs annotated")
  return(damaging_cnvs_noDup)
}



# Function to combine annotated SSMs and CNVs into a totalTable
# This contains the molecular properties used by sysSVM2
make_totalTable = function(
  ssm_anno,         # Output of annotate_ssms()
  cnv_anno,         # Output of annotate_cnvs()
  canonical_drivers # canonical_drivers.rds from the sysSVM2 GitHub repository 
){
  
  # Make sure we have the correct columns in canonical_drivers
  if (is.character(canonical_drivers)) canonical_drivers = readRDS(canonical_drivers)
  canonical_drivers = canonical_drivers %>% select(entrez, geneType)
  
  
  # Summarise mutations at the level of sample-gene pairs
  sample_gene_muts = ssm_anno %>%
    group_by(sample, entrez) %>%
    summarise(
      no_ALL_muts = n(), 
      no_TRUNC_muts = sum(TRUNC_mut), 
      no_NTDam_muts = sum(NTDam_mut), 
      no_GOF_muts = sum(GOF_mut)
    ) %>%
    ungroup
  
  
  # Summarise CNVs at the level of sample-gene pairs
  sample_gene_cnvs = cnv_anno %>%
    select(sample, entrez, Copy_number, CNVGain, CNVLoss) %>%
    unique
  
  
  # Join
  totalTable = full_join(
    sample_gene_muts,
    sample_gene_cnvs,
    by = c("sample", "entrez")
  ) %>%
    replace_na(list(
      no_ALL_muts = 0, no_TRUNC_muts = 0, no_NTDam_muts = 0, no_GOF_muts = 0,
      Copy_number = 2, CNVGain = 0, CNVLoss = 0
    ))
  
  
  # Subset for damaged genes (including GoF in oncogenes and LoF in TSGs)
  totalTable = totalTable %>% 
    subset(no_TRUNC_muts + no_NTDam_muts + no_GOF_muts > 0 | Copy_number == 0 | CNVGain == 1) %>%
    left_join(canonical_drivers, by = "entrez") %>%
    subset(
      is.na(geneType) |
        (geneType == "oncogene" & no_NTDam_muts + no_GOF_muts >= 1 | CNVGain == 1) |
        (geneType == "tumourSuppressor" & no_TRUNC_muts + no_NTDam_muts >= 1 | Copy_number == 0) |
        (geneType == "TP53" & no_TRUNC_muts + no_NTDam_muts + no_GOF_muts >= 1 | Copy_number == 0)
    )
  
  
  # Output
  totalTable = totalTable %>% select(-geneType)
  return(totalTable)
}

