# created with the help of chat.ai using model Anthropic Claude Sonnet 5

# python code was run before

################################################################################
### step 2: Load and reshape merged sqlite ###
################################################################################
library(DBI)
library(RSQLite)
library(dplyr)
library(MsCoreUtils)
library(digest)

load_combined <- function(sqlite_path, source_label) {
  # load data
  con <- dbConnect(SQLite(), sqlite_path)
  df <- dbReadTable(con, "spectra")
  dbDisconnect(con)
  
  # reshape mz and intensity and make sure library source is provided
  df$mz <- lapply(strsplit(df$mz, ";"), as.numeric)
  df$intensity <- lapply(strsplit(df$intensity, ";"), as.numeric)
  df$library_source <- ifelse(is.na(df$library_source) | df$library_source == "",
                              source_label, df$library_source)
  
  df
}

pos_df <- load_combined("combined_pos.sqlite", "MoNA_GNPS")
neg_df <- load_combined("combined_neg.sqlite", "MoNA_GNPS")

################################################################################
### step 3: get MassBank from AnnotationHub and reshape ###
################################################################################
library(AnnotationHub)
library(CompoundDb) # version 1.14.2
library(Spectra)

ah <- AnnotationHub()
query(ah, "MassBank")
mb <- ah[["AH119519"]]

extract_massbank <- function (cdb, polarity_value, polarity) {
  # transform massbank into spectra and filter polarity mode
  sp <- Spectra(cdb)
  sp <- filterPolarity(sp, polarity = polarity_value)
  
  # transform to dataframe and create metadata
  #print(colnames(as.data.frame(compounds(cdb))))
  cmp <- as.data.frame(compounds(cdb, columns = c("compound_id", "name", "formula", 
                                                  "exactmass", "inchikey")))
  #print(colnames(as.data.frame(spectraData(sp))))
  meta <- as.data.frame(spectraData(sp, columns = c("compound_id", "precursorMz", 
                                                    "polarity", "collisionEnergy",
                                                    "adduct")))
  
  # rename precursor mz and collision energy
  colnames(meta) <- c("compound_id", "precursor_mz", "polarity", "collision_energy", "adduct")
  
  # extract peak data
  pd <- peaksData(sp)
  meta$mz <- lapply(pd, function(m) m[,"mz"])
  meta$intensity <- lapply(pd, function(m) m[,"intensity"])
  meta$library_source <- "MassBank"
  
  # modify some columns for merging in step 4
  meta$polarity <- polarity
  meta$collision_energy <- as.character(meta$collision_energy)
  
  # merge all information together into one data frame
  merge(meta, cmp, by = "compound_id", all.x = TRUE)
}

mb_pos <- extract_massbank(mb, 1L, "P")
mb_neg <- extract_massbank(mb, 0L, "N")

################################################################################
### step 4 - harmonize DB together ###
################################################################################
harmonize_dbs <- function (mona_gnps, massbank, sign) {
  dfs <- list(mona_gnps, massbank)
  
  # make sure the columns are identical before binding
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function (d) {
    missing_cols <- setdiff(all_cols, names(d))
    for (col in missing_cols) d[[col]] <- NA
    d[all_cols]
  })
  
  df <- bind_rows(dfs)
  
  # clean InChIKeys and unify missing values
  df$inchikey <- trimws(toupper(df$inchikey)) # removes whitespace at the beginning or end
  df$inchikey[df$inchikey %in% c("", "N/A", "NA", "None")] <- NA
  
  # standardize adduct annotation
  clean_adducts <- function (a, sign = "+"){
    a <- trimws(a) # remove whitespace
    ifelse(
      is.na(a) | a == "",
      NA_character_,
      {
        wrapped <- ifelse(startsWith(a, "["), a, paste0("[", a, "]"))
        has_sign <- grepl("[+-]$", wrapped)
        ifelse(has_sign, wrapped, paste0(wrapped, sign))
      }
    )
  }
  df$adduct <- clean_adducts(df$adduct, sign)
  
  # add source to compound_id to avoid collision
  df$compound_id <- paste(df$library_source, df$compound_id, sep="_")
  
  # use formula to add exactmass, where missing
  mass_missing <- is.na(df$exactmass) & !is.na(df$formula)
  df$exactmass[mass_missing] <- vapply(df$formula[mass_missing], function (f) {
    tryCatch(calculateMass(f), error = function(e) NA_real_)
  }, numeric(1))
  
  df
}

pos_all <- harmonize_dbs(pos_df, mb_pos, "+")
neg_all <- harmonize_dbs(neg_df, mb_neg, "-")

################################################################################
### step 5 - deduplication of the combined DB ###
# differentiate between true duplicates (same compound and spectrum found in multiple DB)
# and distinct spectra of the same compound
# deduplication key = inchikey + adduct + collision energy + spectrum content hash
################################################################################
# round mz and intensity to 4 digits after the comma
spectrum_hash <- function (mz, intensity, digits = 4){
  if (length(mz) == 0) return(NA_character_)
  ord <- order(mz)
  mz_r <- round(mz[ord], digits)
  int_r <- round(intensity[ord], digits)
  digest::digest(paste(mz_r, int_r, collapse = ";"))
}

# deduplicate spectra
deduplicate_spectra <- function (df) {
  # set spectrum hash for each spectrum
  df$spec_hash <- mapply(spectrum_hash, df$mz, df$intensity)
  
  # group spectra by inchikey, adduct, collision energy and spec_hash and select first
  df %>%
    group_by ( inchikey, adduct, collision_energy, spec_hash) %>%
    summarise (
      compound_id = compound_id[1],
      name = name[which.max(nchar(coalesce(name, "")))],
      formula = dplyr::first(na.omit(formula)),
      exactmass = dplyr::first(na.omit(exactmass)),
      precursor_mz = dplyr::first(na.omit(precursor_mz)),
      polarity = dplyr::first(na.omit(polarity)),
      mz = mz[1],
      intensity = intensity[1],
      library_source = paste(sort(unique(library_source)), collapse = ";"),
      n_source = n_distinct(library_source), # cross-library confidence signal
      .groups = "drop"
    ) %>%
    ungroup()
}

pos_dedup <- deduplicate_spectra(pos_all)
neg_dedup <- deduplicate_spectra(neg_all)

################################################################################
### step 6 - assign canonical compound IDs ###
# link multiple spectra correctly to one compound record
################################################################################
assign_canonical_id <- function (df) {
  df$compound_id <- ifelse(
    !is.na(df$inchikey),
    paste0("CMP_", substr(df$inchikey, 1, 14)), # only select first 14 letters of inchikey -> connectivity layer, no stereochemistry
    paste0("CMP_", df$compound_id)
  )
  df
}

pos_dedup <- assign_canonical_id(pos_dedup)
neg_dedup <- assign_canonical_id(neg_dedup)

################################################################################
### step 5 - build compound table and final CompDb for later usage ###
################################################################################
# build compound table only containing entries for unique compound_ids
build_compound_table <- function (df) {
  df %>%
    group_by (compound_id) %>%
    summarise (
      name = name[which.max(nchar(coalesce(name, "")))],
      formula = dplyr::first(na.omit(formula)),
      exactmass = dplyr::first(na.omit(exactmass)),
      inchikey = dplyr::first(na.omit(inchikey)),
      .groups = "drop"
    )
  # add additional columns required for createCompDb
  df$inchi <- NA_character_
  df$synonyms <- NA_character_
  select(df, -c("mz", "intensity"))
}
# create compound table with unique compounds
pos_compounds <- build_compound_table(pos_dedup)
neg_compounds <- build_compound_table(neg_dedup)

# add additional columns required for createCompID
# add spectrum id to all entries of deduplicated tables
# add predicted (whether spectrum was predicted or measured, we suppose measured for all)
# add columns named splash, instrument_type instrument containing NA
# positive mode
pos_dedup$spectrum_id <- seq_len(nrow(pos_dedup))
pos_dedup$predicted <- FALSE
pos_dedup$splash <- NA_character_
pos_dedup$instrument_type <- NA_character_
pos_dedup$instrument <- NA_character_
pos_dedup$polarity <- as.integer(1) # set polarity to 1 for CompoundDb
# negative mode
neg_dedup$spectrum_id <- seq_len(nrow(neg_dedup))
neg_dedup$predicted <- FALSE
neg_dedup$splash <- NA_character_
neg_dedup$instrument_type <- NA_character_
neg_dedup$instrument <- NA_character_
neg_dedup$polarity <- as.integer(0) # set polarity to 0 for CompoundDb

# create general metadata table as required for CompDb, this table should contain
# two columns (name and value) -> usage of make_metadata
meta_df <- make_metadata(
  source = "MoNA;GNPS;MassBank",
  url = NA_character_,
  source_version = 1.0,
  source_date = as.character(Sys.Date()),
  organism = NA_character_
)

# create CompDB object from the provided data
db_file_pos <- createCompDb (
  x = pos_compounds,
  metadata = meta_df,
  msms_spectra = pos_dedup,
  path = ".",
  dbFile = "combined_positive_libray"
)
db_file_neg <- createCompDb (
  x = neg_compounds,
  metadata = meta_df, 
  msms_spectra = neg_dedup,
  path = ".",
  dbFile = "combined_negative_library"
)

cdb_pos <- CompDb(db_file_pos)
cdb_neg <- CompDb(db_file_neg)
