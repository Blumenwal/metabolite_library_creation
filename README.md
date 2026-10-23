# Metabolite_library_creation
Create harmonized and deduplicated spectral library using different library files

This repository can be used to combine, harmonize and deduplicate the mass spectra included in multiple spectral libraries to create on big library file to perform MS1 and MS2 annotation. 
The spectra in the libraries are only separated by their polarity, to create one final library file only containing positive spectra and the second containing negative spectra.

# What data was used, where and when was it downloaded
**GNPS** - downloaded ALL_GNPS_NO_PROPAGATED as msp on 24.08.2026. I opted for the data with no preprocessing. Downloaded from https://external.gnps2.org/gnpslibrary

**MoNA** - downloaded All_LC-MS-MS Orbitrap data on 24.08.2026 as msp from https://massbank.us/downloads

**MassBank** - obtained on 24.08.2026 using the following R-code (that can be found in harmonize_libraries.R):<br>
```
library(AnnotationHub)
library(CompoundDb)

ah <- AnnotationHub()
mb <- ah["AH119519"]]
```

# How to create the combined library:
1. Use the **msp_to_sqlite.ipynb** script to transform the downloaded msp files into sqlite files. For each library, two sqlite files are created. One for each polarity. The sqlite files contain the table spectra with the needed parameters (compound_id, name, inchikey, formula, adduct, prcursor_mz, polarity, collision_energy, mz, intensity and library_source). 
2. Use the **merge_libraries.ipynb** script to merge the library files with identical polarity.
3. Use the **harmonize_libraries.R** to include the MassBank entries, harmonize and deduplicate the spectra, separated on their polarity. The sqlite-files are first read and merged from different libraries, then the MassBank entries are loaded and split by their polarity. Then the entries from all libraries are bound together, before making sure some values are cast correctly ([M+H]+ instead of M+H as adducts). Following this, duplicated entries are removed. These duplicated entries correspond to spectra included in multiple libraries. In the last step, the final CompDb is created, to use it later for the annotation of mass spectra.

# Used R libraries
- DBI - 1.3.0
- RSQLite - 3.53.3
- dplyr - 1.2.1
- MsCoreUtils
- digest - 0.6.39
- AnnotationHub - 4.0.0
- CompoundDb - 1.14.2
- Spectra - 1.20.1

# Used Python libraries
- sqlite3
- re

# Disclaimer regarding AI usage
https://chat-ai.academiccloud.de with the model Anthropic Claude Sonnet 5 was used during the creation of the code in this repository. The code was tested and modified to work for the used versions of packages using the respective manuals.
