# Loading -- Libraries ----------------------------------------------------
library(tidyverse)
library(broom)
library(data.table)
library(DescTools)
library(readxl)
library(vegan)
library(ape)
library(lme4)
library(MuMIn)

library(RColorBrewer)

select <- dplyr::select

# Loading -- datasets -----------------------------------------------------
#Raw metabolomics data
quantRaw <- read_csv('~/Documents/GitHub/greeneMaui/data/raw/quant.csv')


canopus <- read_tsv('~/Documents/GitHub/greeneMaui/data/raw/canopus_structure_summary.tsv')%>%
  rename(superclass = `ClassyFire#superclass`,
         class = `ClassyFire#class`,
         subclass =  `ClassyFire#subclass`,
         featureNumber = mappingFeatureId)%>%
  filter(`ClassyFire#superclass probability` >= 0.7) # this leaves 1350 /1488 feature annotations

## Raw SIRIUS exports
siriusFormula <- read_tsv('~/Documents/GitHub/greeneMaui/data/raw/formula_identifications.tsv')%>%
  rename(featureNumber = 'mappingFeatureId')%>%
  group_by(featureNumber)%>%
  nest()%>%
  mutate(elementalComposition = map(data, ~CHNOSZ::count.elements(.x$molecularFormula)%>%
                                      as.data.frame()%>%
                                      rename(`number` = 1)%>%
                                      rownames_to_column(var = 'element')%>%
                                      pivot_wider(names_from = 'element', values_from = 'number')))%>%
  unnest(c(data, elementalComposition))%>%
  mutate(C = replace_na(C,0),
         H = replace_na(H,0),
         O = replace_na(O,0),
         S = replace_na(S,0),
         N = replace_na(N,0),
         P = replace_na(P,0))

#GNPS outputs -- library Matches and networking information
libraryMatches <- read_tsv('~/Documents/GitHub/greeneMaui/data/raw/libraryMatches.tsv')%>%
  rename(featureNumber = `#Scan#`)

libraryMatchClassyFire <- read_csv('~/Documents/GitHub/greeneMaui/data/raw/greeneLibraryClassyFire.csv')%>%
  rename(featureNumber = `#Scan#`)%>%
  select(Compound_Name, featureNumber, superclass, class, subclass)%>% 
  unique()

nodeInfo <- read_tsv('~/Documents/GitHub/greeneMaui/data/raw/nodeInfo.tsv')%>%
  rename(featureNumber = `cluster index`,
         network = componentindex)


# Metadata 
rawMetadata <- read_csv('~/Documents/HIMB-PRF/projects/mauiDonahue/For Zach/MetaData_08312021.csv')

msMetadata <- read_csv('~/Documents/HIMB-PRF/projects/mauiDonahue/For Zach/MSMS_Metadata.csv')%>%
  mutate(sample = gsub('-', '_', `Sample name`))

landUseRaw <- read_csv('~/Documents/GitHub/greeneMaui/data/raw/CAH_3km_Intersection.csv')

#libray Match Sources
libraryMatchSources <- read_csv('~/Documents/GitHub/greeneMaui/data/analysis/libraryMatchTableCleaned.csv')

#coral cover data
# coralCoverRaw <- read_csv('~/Documents/HIMB-PRF/projects/mauiDonahue/siteMatching/matched_observations.csv')
coralCoverRaw <- read_csv('~/Documents/GitHub/greeneMaui/data/raw/maui_cramp_coralcover.csv')%>%
  mutate(siteGroup = case_when(site_id == 'MaHoS_3' ~ 'Site-01',
                               site_id == 'MaKah_3' ~ 'Site-03',
                               site_id == 'MaOlo_3' ~ 'Site-08',
                               site_id == 'MaMaa_3' ~ 'Site-11',
                               site_id == 'MaMol_8' ~ 'Site-15'))

## The below lines were done to extract InChIKeys from the libraryMatch file and use the ClassyFire Batch file from The Fiehn Lab
##   to then run concise using classyfire instead of NP Classifier

# libInchiKey <- libraryMatches%>%
#   filter(InChIKey != 'N/A')%>%
#   pull(InChIKey)
# 
# write_lines(libInchiKey, '~/Downloads/greeneInChiKeys.txt')
## Pulling classyfire IDs back in and making the library file

# classyFireIDs <- read_csv('~/Downloads/classyfire_20241031084550.csv')
# 
# libraryWithClassyFire <- libraryMatches%>%
#   select(-c(superclass:subclass))%>%
#   left_join(classyFireIDs, by = 'InChIKey')%>%
#   rename(`superclass` = 'Superclass',
#          `class` = 'Class',
#          `subclass` = 'Subclass',
#          `#Scan#` = 'featureNumber')
# #
# write_csv(libraryWithClassyFire, '~/Downloads/greeneLibraryClassyFire.csv')

# Writing functions ------------------------------------------------------
# This function is used while cleaning the data so that we can identify molecular subnetworks that are dominated by background features
# The inputs are:
# 1) dataframe
# 2) the column name of the subnetwork number, 
# 3) the background column that will say either 'background' or 'real' for each feature
# 4) Percent of the subnetwork that has to be real (e.g. minConsensus = 0.5 means greater than 50% of the subnetwork must be real)
flagBackgroundNetworks <- function(data, networkColumn, backgroundColumn, minConsensus = 0.5) {
  require("tidyverse")
  
  flagNetworks <- data%>%
    mutate(count = 1)%>%
    group_by({{networkColumn}}, {{backgroundColumn}})%>%
    summarize(count = sum(count))%>%
    ungroup()%>%
    group_by({{networkColumn}})%>%
    mutate(subnetworkPercentReal = count/sum(count),
           backgroundNetworks = case_when(subnetworkPercentReal > minConsensus ~ 'real',
                                          TRUE ~ 'background'))%>%
    filter({{backgroundColumn}} == 'real')%>%
    select({{networkColumn}}, subnetworkPercentReal, backgroundNetworks)
  
  join <- enquo(networkColumn)
  
  flagExport <- data%>%
    left_join(flagNetworks, by = quo_name(join))%>%
    mutate(backgroundNetworks = case_when(is.na(backgroundNetworks) ~ 'background',
                                          TRUE ~ backgroundNetworks))
  
}

zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

# This is the generic theme I use for plotting so that I dont have to copy and paste it continually throughout the script
genTheme <- function(x) {
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "transparent"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.border = element_rect(color = 'black', fill = "transparent"))
    # panel.border.y = element_rect(color = 'black', fill = "transparent"))
  # panel.grid.major = element_blank(), # get rid of major grid
  # panel.grid.minor = element_blank()) # get rid of minor grid
}

hcaColors <- c("#56B7E9", "#B485D8", "#E56C2F")

# Cleaning -- raw quant --------------------------------------------------
quantClean <- quantRaw%>%
  select(-c(2:13, X787))%>%
  pivot_longer(2:ncol(.), names_to = 'sample', values_to = 'xic')%>%
  mutate(sample = gsub('.mzML Peak area', '', sample))%>%
  rename(featureNumber = `row ID`)


# Sample 158, 186, 94 were all re-run at some point. The re-run samples were removed
# Montipora dataframe
montWm <- quantClean%>% 
  filter(sample %like% '%Batch1%',
         !sample %like% '%stopmethod%')%>%
  separate(sample, c('island', 'species', 'batch', 'sampleNum', 'extra', 'extra2'), sep = '_')%>%
  unite(sampleNum, c('sampleNum', 'extra', 'extra2'), sep = '')%>%
  mutate(sampleNum = gsub('NA', '', sampleNum))%>%
  filter(!sampleNum %like any% c('9420210425141536', '15820210425135754', '186re'))%>%
  left_join(msMetadata%>%
              mutate(sampleNum = as.character(`Extraction Vial No.`))%>%
              select(sampleNum, `Sample name`), 
            by = 'sampleNum')%>%
  rename(SampleID = 'Sample name')%>%
  left_join(rawMetadata,
            by = 'SampleID')%>%
  left_join(nodeInfo%>%
              select(featureNumber, network),
            by = 'featureNumber')
#Porites dataframe
porWm <- quantClean%>% 
  filter(sample %like% '%Batch3%',
         !sample %like% '%stopmethod%')%>%
  separate(sample, c('island', 'species', 'batch', 'sampleVal', 'extra', 'extra2', 'extra3', 'speciesCode'), sep = '_')%>%
  unite(sampleNum, c('sampleVal', 'extra', 'extra2'), sep = '', remove = FALSE)%>%
  unite(SampleID, c('extra', 'extra2', 'extra3', 'speciesCode'), sep = '-')%>%
  mutate(sampleNum = gsub('NA', '', sampleNum))%>%
  filter(!sampleNum %like any% c('9420210425141536', '15820210425135754', '186re'))%>%
  left_join(msMetadata%>%
              mutate(sampleNum = as.character(`Extraction Vial No.`))%>%
              select(sampleNum, `Sample name`), 
            by = 'sampleNum')%>%
  # rename(SampleID = 'Sample name')%>%
  left_join(rawMetadata,
            by = 'SampleID')%>%
  left_join(nodeInfo%>%
              select(featureNumber, network),
            by = 'featureNumber')
  

# Cleaning -- background remoaval -- DEFINING WHICH SPECIES TO LOOK AT-----------------------------------------
# Here we are flagging all features that are background 
# (e.g. half the average peak height in the sample is less than max peak height in the blanks)
# Then we look at each subnetwork, if more than 50% of the features in a subnetwork are background, we remove the whole subnetwork

## Change the original dataframe below depending on the species to look at
blanksFlagged <- montWm%>%
  mutate(netVals = case_when(network == -1 ~ as.numeric(featureNumber)*-1,
                             TRUE ~ network))%>%
  mutate(sampleType = case_when(sampleNum %like any% c('%Blank%', '%Metha%') ~ 'Blank',
                                TRUE ~ 'sample'))%>%
  group_by(featureNumber)%>%
  mutate(sampleXic = case_when(!sampleType %like% '%Blank' ~ xic,
                               TRUE ~ NA_real_),
         meanSample = mean(sampleXic, na.rm = TRUE),
         blankXic = case_when(sampleType %like% '%Blank%' ~ xic,
                              TRUE ~ NA_real_),
         maxBlank = max(blankXic, na.rm = TRUE),
         # transientNum = sum(sampleXic > 5E4, na.rm = TRUE), 
         # Real columns
         background = case_when(meanSample*0.5 > maxBlank ~ 'real',
                                TRUE ~ 'background'))%>%
  ungroup()%>%
  flagBackgroundNetworks(netVals, background, minConsensus = 0.3)
 
#Filtering everything so that we are only looking at how many features are real 
noBlanks <- blanksFlagged%>%
  filter(background == 'real',
         backgroundNetworks == 'real')


#Counting how many features we defined as background features
# 2142 features in 615 subnetworks/single-loop nodes before background removal
# Montipora 1426 features in 400 subnetworks/single-loop nodes after background removal
# Porites 1315 features in 386 subnetworks/single=loop nodes after background removal
blanksFlagged$featureNumber%>%
  unique()%>%
  length()

blanksFlagged%>%
  pull(netVals)%>%
  unique()%>%
  length()

noBlanks$featureNumber%>%
  unique()%>%
  length()

noBlanks%>%
  mutate(netVals = case_when(network == -1 ~ as.numeric(featureNumber)*-1,
                             TRUE ~ network))%>%
  pull(netVals)%>%
  unique()%>%
  length()

# Cleaning -- minimum peak height -----------------------------------------
# First we need to look at a histogram of the XICs in the data to identify what is a rare xic
# By looking at the max abundance of each ion feature we can identify what xic first quartile
maxIntensities <- noBlanks%>%
  filter(sampleType != 'Blank')%>%
  select(featureNumber, xic)%>%
  group_by(featureNumber)%>%
  summarize_if(is.numeric, max)
  
maxIntensities%>%
  summary()

# Anything below the red line we are defining as rare and will remove them
# We removed some specific sites do to sampling issues or were clear outliers:
# Sites were filtered out because they either had too few samples (WM-3 & WM-4), were deep samples (WM-04-Deep, WM-19) or were 
#    obvious outliers according to the HCA analysis (WM-S09-CO2-MC)
# Additionally sites that did not have a Site designation were removed
pdf('~/Documents/GitHub/greeneMaui/data/plots/rareRemoval.pdf', width = 10, height = 10)
noBlanks%>%
  filter(sampleType != 'Blank')%>%
  ggplot(aes(reorder(featureNumber, xic), xic)) +
  geom_point() +
  geom_hline(yintercept = 86813, color = 'red') +
  labs(x = 'Feature', y = 'Max Ion Intensity (XIC)')
dev.off()

noRare <- noBlanks%>%
  mutate(siteGroup = recode(Site, 'WM-01' = "Site-14", 'WM-02' = "Site-11", 'WM-05' = "Site-10", 'WM-06' = "Site-12", 'WM-08' = "Site-08", 
                            'WM-09' = "Site-04", 'WM-10' = "Site-03", 'WM-11' = "Site-01", 'WM-12' = "Site-06", 'WM-13' = "Site-09",
                            'WM-14' = "Site-07", 'WM-15' = "Site-15", 'WM-16' = "Site-02", 'WM-17' = "Site-16", 'WM-18' = "Site-05",
                            'WM-20' = "Site-13"))%>%
  filter(sampleType != 'Blank',
         !Site %in% c("WM-03", "WM-04", "WM-04-Deep", "WM-19", 'WM-07'),
         SampleID != "WM-S09-C02-MC",
         !is.na(siteGroup))%>%
  group_by(featureNumber)%>%
  filter(max(xic, na.rm = TRUE) > 86813)

# Counting how many features we have remaining after removing rare features
# Montipora 966 features belonging to 266 subnetworks/single loop nodes remain
# Porites 786 features bolong to 240 subnetworks/single loop nodes remain
noRare%>%
  pull(featureNumber)%>%
  unique()%>%
  length()

noRare%>%
  mutate(netVals = case_when(network == -1 ~ as.numeric(featureNumber)*-1,
                      TRUE ~ network))%>%
  pull(netVals)%>%
  unique()%>%
  length()


# Cleaning -- making the stats working data frame -------------------------
# This is the working data frame (WDF) which we will use for all stats and plotting
# It is all the raw feature xic and log10 transformed xic for each sample with metadata included
# The HCAGrouping is calculated from a cluster dendogram lower in this code 
# For the sake of simplicity, those groupings were added to the working data frame here
statsWdf <- noRare%>%
  ungroup()%>%
  select(-c(sampleType:backgroundNetworks))%>%
  group_by(SampleID)%>%
  mutate(log10 = log10(xic + 1),
         ra = xic/sum(xic),
         asin = asin(sqrt(ra)),
         HCAgrouping = recode(siteGroup, 'Site-15' = 1, 'Site-06' = 1, 'Site-07' = 1, 'Site-04' = 1, 'Site-09' = 1, 
                              'Site-12' = 0, #Montipora Groups
                              'Site-01' = 2, #Montipora Groups
                              # 'Site-12' = 2, #Porites Groups
                              # 'Site-01' = 1, #Porites Groups
                              'Site-10' = 2, 'Site-08' = 2, 'Site-11' = 2, 'Site-14' = 2,
                              'Site-16' = 3, 'Site-05' = 3, 'Site-13' = 3, 'Site-02' = 3, 'Site-03' = 3),
         HCAgrouping = as.factor(HCAgrouping),
         HCAgrouping = fct_relevel(HCAgrouping, '1','2','3'),
         netVals = case_when(network == -1 ~ as.numeric(featureNumber)*-1,
                                    TRUE ~ network))

statsWdf%>% 
  select(-c(network, netVals, log10, ra, asin))%>%
  pivot_wider(names_from = 'featureNumber', values_from = 'xic')%>% 
  select(SampleID, siteGroup)%>%
  mutate(n = 1)%>%
  group_by(siteGroup)%>%
  summarize_if(is.numeric, sum)%>%
  summary()

# SET SEED ----------------------------------------------------------------
set.seed(23068)



# STATS -- LMER -- Metabolites predicted by SITE ----------------------------------
lmerHCAGroups <- statsWdf%>%
  group_by(netVals)%>%
  nest()%>%
  # mutate(data = map(data, ~lmer(asin ~ HCAgrouping + (1|featureNumber) + (1|siteGroup), data = .x, 
  #                               control =  lmerControl(check.nlev.gtr.1 = "ignore",
  #                                                      check.conv.singular = 'ignore',
  #                                                      check.nobs.vs.nRE = 'ignore'))%>%
  mutate(data = map(data, ~lmer(asin ~ siteGroup + (1|featureNumber), data = .x,
                                control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                       check.conv.singular = 'ignore',
                                                       check.nobs.vs.nRE = 'ignore'))%>%
                      car::Anova()%>%
                      .[['Pr(>Chisq)']]))

significantValues <- lmerHCAGroups%>%
  unnest(data)%>%
  ungroup()%>%
  mutate(pAdjusted = p.adjust(data, method = 'BH'))%>%
  filter(pAdjusted < 0.05)


# 242 out of 260 subnetworks were signicant using asin ra for MONTIPORA FOR SITE 
# 189 out of 260 subnetworks were signicant using asin ra for PORITES
nrow(significantValues)

hcaSigSubnets <- significantValues$netVals



# Figure 1a --  HCA of metabolomes -----------------------------------------------
siteAverageMetabolomes <- statsWdf%>%
  left_join(canopus%>%
              select(superclass, class, subclass, featureNumber), by = 'featureNumber')%>%
  filter(netVals %in% hcaSigSubnets)%>%
  group_by(siteGroup, netVals, SampleID, superclass, class, subclass)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  mutate(ra = asin(sqrt(ra)))%>%
  unite(annotation, c(superclass, class, subclass), sep = ';')%>%
  group_by(siteGroup, annotation)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(siteGroup, annotation, ra)%>%
  group_by(annotation)%>%
  mutate(ra = zscore(ra))%>%
  ungroup()%>%
  separate(annotation, c('superclass', 'class', 'subclass'), sep = ';')%>%
  mutate(deepestAnnotation = case_when(subclass != 'NA' ~ subclass,
                                       class != 'NA' & subclass == 'NA' ~ class,
                                       TRUE ~ superclass))%>%
  select(-c(superclass, class, subclass))%>%
  pivot_wider(names_from = 'deepestAnnotation', values_from = 'ra')%>%
  column_to_rownames(var = 'siteGroup')


siteAverageMetabolomes%>%
  pheatmap::pheatmap()

#dendo grouping for montipora
dendo <- pheatmap::pheatmap(siteAverageMetabolomes)$tree_row
dendo <- dendextend::rotate(dendo, order = c('Site-12','Site-15', 'Site-06', 'Site-07', 'Site-04', 'Site-09', 
                                                 'Site-01', 'Site-10', 'Site-08', 'Site-11', 'Site-14',
                                                 'Site-16', 'Site-05', 'Site-13', 'Site-02', 'Site-03'))

annotationList <- list(Site = c('Site-15' = '#E56C2F', 'Site-06' = '#E56C2F', 'Site-07' = '#E56C2F', 'Site-04' = '#E56C2F', 'Site-09' = '#E56C2F', 
                       'Site-12' ='#B485D8', 'Site-01' = '#B485D8', 'Site-10' = '#B485D8','Site-08' = '#B485D8', 'Site-11' = '#B485D8', 'Site-14' = '#B485D8',
                       'Site-16' = '#56B7E9', 'Site-05' = '#56B7E9', 'Site-13' = '#56B7E9', 'Site-02' = '#56B7E9', 'Site-03' = '#56B7E9'))

annotationColorDf <- data.frame(Site = c('Site-15', 'Site-06', 'Site-07', 'Site-04', 'Site-09', 
                                         'Site-12', 'Site-01', 'Site-10', 'Site-08', 'Site-11', 'Site-14',
                                         'Site-16', 'Site-05', 'Site-13', 'Site-02', 'Site-03'),
                                row.names = c('Site-15', 'Site-06', 'Site-07', 'Site-04', 'Site-09', 
                                              'Site-12', 'Site-01', 'Site-10', 'Site-08', 'Site-11', 'Site-14',
                                              'Site-16', 'Site-05', 'Site-13', 'Site-02', 'Site-03'))

pdf('~/Documents/GitHub/greeneMaui/data/plots/biclusterSubnetworkAverages.pdf', width = 20, height = 12)  
pheatmap::pheatmap(siteAverageMetabolomes, cluster_rows = dendo, color = brewer.pal(n = 9, name = "Greys"),
                   angle_col = 45)
                   # annotation_color = annotationList, annotation_col = annotationColorDf)
dev.off()

labels <- c('Site-15', 'Site-06', 'Site-07', 'Site-04', 'Site-09', 
            'Site-12', 'Site-01', 'Site-10', 'Site-08', 'Site-11', 'Site-14',
            'Site-16', 'Site-05', 'Site-13', 'Site-02', 'Site-03')

#dendo grouping for porites
dendoPor <- pheatmap::pheatmap(siteAverageMetabolomes)$tree_row
dendoPor <- dendextend::rotate(dendo, order = c('Site-03', 'Site-05', 'Site-16','Site-02', 'Site-13',
                                                'Site-14', 'Site-08','Site-10', 'Site-11','Site-12',
                                                'Site-15', 'Site-09', 'Site-06', 'Site-07', 'Site-01', 'Site-04'))

pdf('~/Documents/GitHub/greeneMaui/data/plots/poritesBiclusterSubnetworkAverages.pdf', width = 20, height = 12)  
pheatmap::pheatmap(siteAverageMetabolomes, cluster_rows = dendoPor, color = brewer.pal(n = 9, name = "Greys"),
                   angle_col = 45)
dev.off()


# Visual -- Environmental variability across hca clusters -----------------
landUse <- landUseRaw%>%
  select(Site, Hab_Status, Shape_Area)%>%
  group_by(Site, Hab_Status)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  add_row(Site = 'WM-15', Hab_Status = 'Heavily Disturbed', Shape_Area = 0)%>%
  group_by(Site)%>%
  mutate(habPercent = Shape_Area/sum(Shape_Area, na.rm = TRUE))%>%
  ungroup()%>%
  filter(Hab_Status == 'Heavily Disturbed')

siteLatLong <- rawMetadata%>%
  left_join(statsWdf%>%
              ungroup()%>%
              select(Site, siteGroup, HCAgrouping), by = 'Site')%>%
  mutate(HCAgrouping = as.factor(HCAgrouping),
         HCAgrouping = fct_relevel(HCAgrouping, '3','2','1'))%>%
  select(siteGroup, HCAgrouping, TNC.Lat, TNC.Long)%>%
  filter(!is.na(siteGroup))%>%
  unique()

write_csv(siteLatLong, '~/Downloads/MauiContaminantSiteLatLong.csv')

enviroVariability <- rawMetadata%>%
  left_join(statsWdf%>%
              ungroup()%>%
              select(Site, siteGroup, HCAgrouping), by = 'Site')%>%
  left_join(landUse, by = 'Site')%>%
  select(TNC.Coral.Cover, TNC.Coral.Disease, TNC.Phosphate.Benthic, 
         TNC.Silicate.Benthic, TNC.NN.Benthic, TNC.Ammonia.Benthic,
         TNC.Coral.Diversity, TNC.Raw.Resilience.Score,
         TNC.Herbivore.Biomass, habPercent, HCAgrouping, Site, siteGroup)%>%
  unique()%>%
  filter(!is.na(TNC.Coral.Cover),
         siteGroup %in% labels,
         siteGroup != 'Site-12'
         )%>%
  mutate(TNC.Phosphate.Benthic.log = log10(TNC.Phosphate.Benthic), 
         TNC.Silicate.Benthic.log = log10(TNC.Silicate.Benthic), 
         TNC.NN.Benthic.log = log10(TNC.NN.Benthic),
         TNC.Ammonia.Benthic.log = log10(TNC.Ammonia.Benthic),
         N2P = TNC.NN.Benthic/TNC.Phosphate.Benthic,
         site = siteGroup,
         # HCAgrouping = recode(siteGroup, 'Site-15' = 1, 'Site-06' = 1, 'Site-07' = 1, 'Site-04' = 1, 'Site-09' = 1, 
         #                      'Site-12' = 2, 'Site-01' = 2, 'Site-10' = 2, 'Site-08' = 2, 'Site-11' = 2, 'Site-14' = 2,
         #                      'Site-16' = 3, 'Site-05' = 3, 'Site-13' = 3, 'Site-02' = 3, 'Site-03' = 3),
         HCAgrouping = as.factor(HCAgrouping),
         HCAgrouping = fct_relevel(HCAgrouping, '3','2','1'),
         site = as.factor(site),
         site = fct_relevel(site, 'Site-15', 'Site-06', 'Site-07', 'Site-04', 'Site-09',
                            'Site-01', 'Site-10', 'Site-08', 'Site-11', 'Site-14',
                            'Site-16', 'Site-05', 'Site-13', 'Site-02', 'Site-03')) # Montipora order


# Plotting the N and P values of benthic systems
pdf('~/Documents/GitHub/greeneMaui/data/plots/environmentalParameterViolins.pdf', width = 6, height = 10)
enviroVariability%>%
  ggplot(aes(TNC.NN.Benthic.log, HCAgrouping, color = HCAgrouping)) +
  geom_point(stat = 'identity', alpha= 0.5, size = 3) +
  # scale_y_discrete(limits=rev) +
  geom_boxplot() +
  scale_color_manual(values = rev(hcaColors)) +
  genTheme() +
  labs(y = 'Cluster', x = 'Benthic Nitrogen') +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        legend.position = 'None')

enviroVariability%>%
  ggplot(aes(TNC.Ammonia.Benthic.log, HCAgrouping, color = HCAgrouping)) +
  geom_boxplot() +
  geom_point(stat = 'identity', alpha= 0.5, size = 3) +
  # scale_y_discrete(limits=rev) +
  scale_color_manual(values = rev(hcaColors)) +
  genTheme() +
  labs(y = 'Cluster', x = 'Benthic Ammonium') +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        legend.position = 'None')

enviroVariability%>%
  ggplot(aes(TNC.Phosphate.Benthic.log, HCAgrouping, color = HCAgrouping)) +
  geom_boxplot() +
  geom_point(stat = 'identity', alpha= 0.5, size = 3) +
  # scale_y_discrete(limits=rev) +
  genTheme() +
  scale_color_manual(values = rev(hcaColors)) +
  labs(y = 'Cluster', x = 'Benthic Phosphorous') +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        legend.position = 'None')

enviroVariability%>%
  ggplot(aes(N2P, HCAgrouping, color = HCAgrouping)) +
  geom_boxplot() +
  geom_point(stat = 'identity', alpha= 0.5, size = 3) +
  # scale_y_discrete(limits=rev) +
  genTheme() +
  scale_color_manual(values = rev(hcaColors)) +
  labs(y = 'Cluster', x = 'Benthic Nitrogen/Benthic Phosphorous') +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        legend.position = 'None')

enviroVariability%>%
  ggplot(aes(TNC.Silicate.Benthic.log, HCAgrouping, color = HCAgrouping)) +
  geom_boxplot() +
  geom_point(stat = 'identity', alpha= 0.5, size = 3) +
  # scale_y_discrete(limits=rev) +
  genTheme() +
  scale_color_manual(values = rev(hcaColors)) +
  labs(y = 'Cluster', x = 'Benthic Silicate') +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        legend.position = 'None')

enviroVariability%>%
  ggplot(aes(habPercent*100, HCAgrouping, color = HCAgrouping)) +
  geom_boxplot() +
  geom_point(stat = 'identity', alpha= 0.5, size = 3) +
  # scale_y_discrete(limits=rev) +
  genTheme() +
  scale_color_manual(values = rev(hcaColors)) +
  labs(y = 'Cluster', x = 'Terrestrial percent disturbance') +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        legend.position = 'None')

enviroVariability%>%
  select(HCAgrouping, siteGroup, TNC.Ammonia.Benthic.log, TNC.Phosphate.Benthic.log, TNC.Silicate.Benthic.log, habPercent)%>%
  unite(site, c(HCAgrouping, siteGroup), sep = ' ')%>%
  pivot_longer(3:ncol(.), names_to = 'responseVar', values_to = 'values')%>%
  group_by(responseVar)%>%
  mutate(values = zscore(values),
         site = as.factor(site),
         # site = fct_relevel(site, '0 Site-12', '1 Site-15', '1 Site-06', '1 Site-07', '1 Site-04', '1 Site-09', 
         #                    '2 Site-01', '2 Site-10', '2 Site-08', '2 Site-11', '2 Site-14',
         #                    '3 Site-16', '3 Site-05', '3 Site-13', '3 Site-02', '3 Site-03'))%>% # Montipora order
         site = fct_relevel(site, '3 Site-03', '3 Site-05', '3 Site-16', '3 Site-02', '3 Site-13', 
                            '2 Site-14', '2 Site-08', '2 Site-10', '2 Site-11', '2 Site-12', 
                            '1 Site-15', '1 Site-09', '1 Site-06', '1 Site-07', '1 Site-01', '1 Site-04'
                            ))%>% # Porites order
  ungroup()%>%
  arrange(site)%>%
  pivot_wider(names_from = 'responseVar', values_from = 'values')%>%
  column_to_rownames(var = 'site')%>%
  pheatmap::pheatmap(cluster_rows = FALSE, color = brewer.pal(n = 9, name = "BuGn"))


dev.off()


pdf('~/Documents/GitHub/greeneMaui/data/plots/EnviroBiplotsExpanded.pdf', width = 15, height = 10)
enviroVariability%>%
  select(HCAgrouping, siteGroup, TNC.Ammonia.Benthic.log, TNC.Phosphate.Benthic.log, TNC.Silicate.Benthic.log, habPercent)%>%
  unite(site, c(HCAgrouping, siteGroup), sep = ' ')%>%
  pivot_longer(3:ncol(.), names_to = 'responseVar', values_to = 'values')%>%
  group_by(responseVar)%>%
  mutate(values = zscore(values),
         site = as.factor(site),
         # site = fct_relevel(site, '0 Site-12', '1 Site-15', '1 Site-06', '1 Site-07', '1 Site-04', '1 Site-09', 
         #                    '2 Site-01', '2 Site-10', '2 Site-08', '2 Site-11', '2 Site-14',
         #                    '3 Site-16', '3 Site-05', '3 Site-13', '3 Site-02', '3 Site-03'))%>% # Montipora order
         site = fct_relevel(site, '3 Site-03', '3 Site-05', '3 Site-16', '3 Site-02', '3 Site-13', 
                            '2 Site-14', '2 Site-08', '2 Site-10', '2 Site-11', '2 Site-12', 
                            '1 Site-15', '1 Site-09', '1 Site-06', '1 Site-07', '1 Site-01', '1 Site-04'
         ))%>% # Porites order
  ungroup()%>%
  arrange(site)%>%
  pivot_wider(names_from = 'responseVar', values_from = 'values')%>%
  column_to_rownames(var = 'site')%>%
  pheatmap::pheatmap(color = brewer.pal(n = 9, name = "BuGn"))

enviroVariability%>%
  select(HCAgrouping, siteGroup, TNC.Coral.Disease, TNC.Coral.Cover, TNC.Raw.Resilience.Score, TNC.Herbivore.Biomass, TNC.NN.Benthic.log, TNC.Ammonia.Benthic.log, TNC.Phosphate.Benthic.log, TNC.Silicate.Benthic.log, habPercent)%>%
  unite(site, c(HCAgrouping, siteGroup), sep = ' ')%>%
  pivot_longer(3:ncol(.), names_to = 'responseVar', values_to = 'values')%>%
  group_by(responseVar)%>%
  mutate(values = zscore(values),
         site = as.factor(site),
         # site = fct_relevel(site, '0 Site-12', '1 Site-15', '1 Site-06', '1 Site-07', '1 Site-04', '1 Site-09', 
         #                    '2 Site-01', '2 Site-10', '2 Site-08', '2 Site-11', '2 Site-14',
         #                    '3 Site-16', '3 Site-05', '3 Site-13', '3 Site-02', '3 Site-03'))%>% # Montipora order
         site = fct_relevel(site, '3 Site-03', '3 Site-05', '3 Site-16', '3 Site-02', '3 Site-13', 
                            '2 Site-14', '2 Site-08', '2 Site-10', '2 Site-11', '2 Site-12', 
                            '1 Site-15', '1 Site-09', '1 Site-06', '1 Site-07', '1 Site-01', '1 Site-04'
         ))%>% # Porites order
  ungroup()%>%
  arrange(site)%>%
  pivot_wider(names_from = 'responseVar', values_from = 'values')%>%
  column_to_rownames(var = 'site')%>%
  pheatmap::pheatmap(color = brewer.pal(n = 9, name = "BuGn"))

enviroVariability%>%
  select(HCAgrouping, siteGroup, TNC.Coral.Disease, TNC.Coral.Cover, TNC.Raw.Resilience.Score, TNC.Herbivore.Biomass)%>%
  unite(site, c(HCAgrouping, siteGroup), sep = ' ')%>%
  pivot_longer(3:ncol(.), names_to = 'responseVar', values_to = 'values')%>%
  group_by(responseVar)%>%
  mutate(values = zscore(values),
         site = as.factor(site),
         # site = fct_relevel(site, '0 Site-12', '1 Site-15', '1 Site-06', '1 Site-07', '1 Site-04', '1 Site-09', 
         #                    '2 Site-01', '2 Site-10', '2 Site-08', '2 Site-11', '2 Site-14',
         #                    '3 Site-16', '3 Site-05', '3 Site-13', '3 Site-02', '3 Site-03'))%>% # Montipora order
         site = fct_relevel(site, '3 Site-03', '3 Site-05', '3 Site-16', '3 Site-02', '3 Site-13', 
                            '2 Site-14', '2 Site-08', '2 Site-10', '2 Site-11', '2 Site-12', 
                            '1 Site-15', '1 Site-09', '1 Site-06', '1 Site-07', '1 Site-01', '1 Site-04'
         ))%>% # Porites order
  ungroup()%>%
  arrange(site)%>%
  pivot_wider(names_from = 'responseVar', values_from = 'values')%>%
  column_to_rownames(var = 'site')%>%
  pheatmap::pheatmap(color = brewer.pal(n = 9, name = "BuGn"))
dev.off()


# Stats -- enviroVariability -- ANOVA -------------------------------------
enviroVariability%>%
  select(HCAgrouping, siteGroup, TNC.Ammonia.Benthic.log, TNC.Phosphate.Benthic.log, TNC.Silicate.Benthic.log, habPercent)%>%
  pivot_longer(3:ncol(.), names_to = 'responseVar', values_to = 'values')%>%
  group_by(responseVar)%>%
  nest()%>%
  mutate(data = map(data, ~car::Anova(lm(values~HCAgrouping, data = .x), type = 'II')%>%
                      tidy()
                      ))%>%
  unnest(data)


# Visual -- TIC -----------------------------------------------------------
# Have to figure out how they did the sampling. Were all the samples chunks?
# Might have to calculate RA so that we are not dealing with concentration issues
statsWdf%>%
  unite(local, c(Site, HCAgrouping), sep = '_')%>%
  group_by(SampleID, local)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  ggplot(aes(local, xic)) +
  geom_boxplot()


# Supplemental Figure 1 -- PCoA -- grouped at the subnetwork consensus level four clusters ---------------------------------
pcoaSuperclassDf <- statsWdf%>%
  filter(netVals %in% hcaSigSubnets,
         siteGroup != 'Site-12'
         )%>%
  left_join(canopus%>%
              select(superclass, class, subclass, featureNumber), by = 'featureNumber')%>%
  # left_join(conciseRaw, by = 'featureNumber')%>%
  unite(annotation, c(superclass, class, subclass), sep = ';')%>%
  group_by(annotation, SampleID, Site, ColonyID, HCAgrouping)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  # mutate(log10 = log10(xic + 1))%>%
  mutate(asin = asin(sqrt(ra)))%>%
  select(annotation, SampleID, Site, ColonyID, HCAgrouping, asin)%>%
  unite(sample, c('SampleID', 'Site', 'ColonyID',  'HCAgrouping'), sep = ';')%>%
  pivot_wider(names_from = 'annotation', values_from = 'asin')%>%
  column_to_rownames(var = 'sample')


pcoaSuperclass <- pcoaSuperclassDf%>%
  vegdist()%>%
  pcoa()

#Plotting PCoAs
pdf('~/Documents/GitHub/greeneMaui/data/plots/pcoaSubclass3Clusters.pdf', width = 15, height = 10)
pcoaSuperclass$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('SampleID', 'Site', 'ColonyID',  'HCAgrouping'), sep = ';')%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = HCAgrouping)) +
  geom_point(size = 8) +
  # ggforce::geom_mark_ellipse(aes(fill = HCAgrouping)) +
  scale_color_manual(values = hcaColors) +
  # scale_fill_manual(values = hcaColors4) +
  # ggforce::geom_mark_ellipse(aes(fill = HCAgrouping)) +
  # geom_text(aes(label = Day), size = 12) +
  # scale_color_manual(values = healthColors) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25)) +
  xlab(str_c("Axis 1", " (", round(pcoaSuperclass$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(pcoaSuperclass$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))
dev.off()


# Stats Fig. S1 -- PERMANOVA for subclass PcOA ------------------------------------
pcoaSuperclassDf%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('SampleID', 'Site', 'ColonyID',  'HCAgrouping'), sep = ';')%>%
  adonis2(.[5:ncol(.)]~HCAgrouping, data = .)

pairwiseDf <- pcoaSuperclassDf%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('SampleID', 'Site', 'ColonyID',  'HCAgrouping'), sep = ';')

pairwiseAdonis::pairwise.adonis2(pairwiseDf[5:ncol(pairwiseDf)]~HCAgrouping, data = pairwiseDf)



# Stats -- Prechecks Multiple Regression -- pcoa Axis 1 vs environmental variables  -------------------------------------------
pcoaAxis <- pcoaSuperclass$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = 'sample')%>%
  select(sample, Axis.1)%>%
  separate(sample, c('SampleID', 'Site', 'ColonyID',  'HCAgrouping'), sep = ';')%>%
  group_by(Site, HCAgrouping)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

pcoaMulReg <- pcoaAxis%>%
  left_join(enviroVariability, by = c('Site', 'HCAgrouping'))

# Checking to see if there is too much correlation between any of the 
pcoaMulReg%>%
  ungroup()%>%
  select(Axis.1, TNC.Ammonia.Benthic.log, TNC.NN.Benthic.log, TNC.Phosphate.Benthic.log, TNC.Silicate.Benthic.log, habPercent)%>%
  cor()%>%
  corrplot::corrplot()
  

# Stats -- Multiple Regression AIC ----------------------------------------
lm_test <- lm(Axis.1 ~ TNC.Ammonia.Benthic.log + TNC.NN.Benthic.log + TNC.Phosphate.Benthic.log + TNC.Silicate.Benthic.log + habPercent, 
              data = pcoaMulReg)
# Other Variables:
#  TNC.Coral.Cover + TNC.Coral.Disease + TNC.Coral.Diversity + TNC.Raw.Resilience.Score + TNC.Herbivore.Biomass +

#Dredge
options(na.action = na.fail)
dredge_n_select <- MuMIn::dredge(lm_test)
sink('./analysis/dredge_model_selection.txt')
head(dredge_n_select)
sink()

# Stats -- MRegression ----------------------------------------------------
#Montipora
lmLimited <- lm(Axis.1 ~ TNC.Phosphate.Benthic.log + TNC.Silicate.Benthic.log + TNC.Ammonia.Benthic.log + habPercent, data = pcoaMulReg)


#Porites
# lmLimited <- lm(Axis.1 ~ TNC.Silicate.Benthic.log + TNC.Ammonia.Benthic.log + habPercent, data = pcoaMulReg)
#Summary
summary(lmLimited)

# Montipora 56.77%
# Porites 50.03%





# Stats -- lm a with coral cover -----------------------------------
coralCoverClean <- coralCoverRaw%>%
  mutate(yearGroup = ifelse(Year < 2015, 'pre-Bleaching', 'post-Bleaching'),
         siteGroup = as.factor(siteGroup),
         siteGroup = fct_relevel(siteGroup, 'Site-15', 'Site-01', 'Site-08', 'Site-11', 'Site-03'))

coralGroupRegression <- coralCoverClean%>%
  group_by(siteGroup, yearGroup)%>%
  nest()%>%
  mutate(data = map(data, ~ lm(Coral ~ Year, data = .x)%>%
                      summary()),
         p.value = map(data, ~tidy(.x)%>%
                         filter(term == 'Year')%>%
                         pull(p.value)),
         estimate = map(data, ~tidy(.x)%>%
                          filter(term == 'Year')%>%
                          pull(estimate)),
         r2 = map(data, ~.x$r.squared))%>%
  unnest(c(p.value, estimate, r2))%>%
  ungroup()%>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))


coralGroupRegression

# Visualization -- Regression with coral cover ------------------------------------
coralCoverPlotting <- coralCoverClean%>%
  filter(!siteGroup %in% c('Site-08', 'Site-15'))

#Have to separate out everything so that I can make each line appropriately dashed based on significance
plottingSite15 <- coralCoverClean%>%
  filter(siteGroup == 'Site-15')%>%
  mutate(line = 'dashed')

plottingSite08Pre <- coralCoverClean%>%
  filter(siteGroup == 'Site-08',
         yearGroup == 'pre-Bleaching')%>%
  mutate(line = ifelse(yearGroup == 'post-Bleaching', 'solid', 'dashed'))

plottingSite08Post <- coralCoverClean%>%
  filter(siteGroup == 'Site-08',
         yearGroup == 'post-Bleaching')%>%
  mutate(line = ifelse(yearGroup == 'post-Bleaching', 'solid', 'dashed'))


pdf('~/Documents/GitHub/greeneMaui/data/plots/coralCoverChangeSiteBreakdown.pdf', width = 10, height = 40)
coralCoverPlotting%>%
  ggplot(aes(Year, Coral, color = yearGroup)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_point(data = plottingSite15, size = 4) +
  geom_smooth(data = plottingSite15, method = 'lm', linetype = 'dashed', se = FALSE, size = 3) +
  geom_point(data = plottingSite08Pre, size = 4) +
  geom_smooth(data = plottingSite08Pre, method = 'lm', linetype = 'solid', size = 3) +
  geom_point(data = plottingSite08Post, size = 4) +
  geom_smooth(data = plottingSite08Post, method = 'lm', linetype = 'dashed', se = FALSE, size = 3) +
  facet_wrap(~siteGroup, ncol = 1) +
  genTheme() +
  scale_color_manual(values = c('#9E9B97', '#5D3A00')) + 
  theme(panel.spacing = unit(10, "lines"))
dev.off()


# Stats -- benthic nutrients ----------------------------------------------
# enviroStatsDf <- enviroVariability%>%
#   filter(siteGroup != 'Site-12')%>%
#   select(TNC.Phosphate.Benthic.log, TNC.Ammonia.Benthic.log,
#          TNC.Silicate.Benthic.log, TNC.Herbivore.Biomass, habPercent, site, HCAgrouping, siteGroup)%>%
#   pivot_longer(TNC.Phosphate.Benthic.log:habPercent, names_to = 'responseVar', values_to = 'values')%>%
#   group_by(responseVar)%>%
#   nest()%>%
#   mutate(data = map(data, ~aov(values ~ HCAgrouping , data = .x)%>%
#                       tidy()))%>%
#   unnest(data)
# 

# Figure 2a -- Beta diversity Shannon ----------------------------------------
betaDiv <- statsWdf%>%
  filter(siteGroup != 'Site-12')%>%
  filter(netVals %in% hcaSigSubnets)%>%
  select(siteGroup, HCAgrouping, SampleID, featureNumber, ra)%>%
  pivot_wider(names_from = 'featureNumber', values_from = 'ra')%>%
  group_by(siteGroup, SampleID, HCAgrouping)%>%
  nest()%>%
  mutate(shannon = map(data, ~ diversity(.x, index = 'shannon')),
         specN = map(data, ~ rowSums(.x > 0, na.rm = TRUE)))%>%
  unnest(c(shannon, specN))%>%
  mutate(eveness = shannon/log(specN))

pdf('~/Documents/GitHub/greeneMaui/data/plots/poritesshannonDiversity.pdf', width = 12, height = 10)
betaDiv%>%
  group_by(siteGroup)%>%
  mutate(HCAgrouping = as.factor(HCAgrouping),
         HCAgrouping = fct_relevel(HCAgrouping, '1', '2', '3', '4')) %>%
  ggplot(aes(HCAgrouping, shannon, fill = HCAgrouping)) +
  # geom_boxplot() +
  geom_violin() +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               colour = "white") +
  genTheme() +
  scale_fill_manual(values = hcaColors) +
  labs(x = 'HCA grouping', y = 'Shannon diversity')
dev.off()


# Stats -- Shannon Diversity ----------------------------------------------
# Testing if there is a significant difference between HCA groups in terms of diversity of the molecules
# 4.056368e-62 for mont
# 2.701403e-13 for porites
betaDiv%>%
  lmer(shannon~HCAgrouping + (1|siteGroup), data = .)%>%
  car::Anova()%>%
  .[['Pr(>Chisq)']]

# Tukey p-value < 0.001 for all Tukey differences
shannonPvalues <- betaDiv%>%
  lmer(shannon~HCAgrouping + (1|siteGroup), data = .)%>%
  multcomp::glht(., linfct = multcomp::mcp(HCAgrouping = "Tukey"))%>%
  summary(., test = multcomp::adjusted('BH'))%>%
  tidy()%>%
  mutate(test = 'shannon')

# Figure 2b-c -- elemental Ratios -- N:C vs P:C and NOSC plots -------------------------------------------------------
## N:C, P:C were weighted by relative abundance of the molecules
## NOSC and gCox were weighted by (xic*C) /sum(xic*C)
## NOSC is the nominal oxidation state of carbon for each molecule
## gCox is the potential energy of each molecule from catabolism
## N:C, P:C, NOSC and gCox are further summed within each sample to represent the total N:C
## Summed HCA group NOSC, gCox, N:C, P:C across all samples and metabolites 
elementalComposition <- statsWdf%>%
  filter(siteGroup != 'Site-12')%>%
  ungroup()%>%
  left_join(siriusFormula%>%
              select(featureNumber, C:P), by = 'featureNumber')%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  mutate(p2c = P/C,
         n2c = N/C,
         propXic = xic/sum(xic),
         NOSC = (-((4*C + H - 3*N - 2*O + 5*P - 2*S)/C)+4),
         NOSC = ifelse(NOSC <= -4 | NOSC >= 4 , NA_real_, NOSC),
         # x = xic*C,
         # weight = x/sum(x, na.rm = TRUE),
         gCOX = (60.3) - 28.5*NOSC,
         weightN2c = n2c*propXic,
         weightP2c = p2c*propXic,
         weightNOSC = NOSC*propXic,
         weightGcox = gCOX*propXic)%>%
  ungroup()%>%
  filter(netVals %in% hcaSigSubnets)%>%
  filter(n2c < 1)

# This dataframe allows us to calculate the standard errors to plot directly on top of the scatter plot
elementalAverages <- elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(HCAgrouping)%>%
  mutate(n = 1,
         n = sum(n),
         p2cErr = sd(weightP2c)/sqrt(n),
         n2cErr = sd(weightN2c)/sqrt(n))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

pdf('~/Documents/GitHub/greeneMaui/data/plots/n2cVp2cAverageSite.pdf', width = 15, height = 10)
elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(siteGroup, HCAgrouping)%>%
  mutate(n = 1,
         n = sum(n),
         yerr = sd(weightP2c)/sqrt(n),
         xerr = sd(weightN2c)/sqrt(n))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(weightN2c, weightP2c, color = HCAgrouping)) +
  # geom_boxplot(alpha = 1) +
  geom_errorbar(aes(ymin = weightP2c - yerr, ymax = weightP2c + yerr), linewidth = 1) +
  geom_errorbarh(aes(xmin = weightN2c - xerr, xmax = weightN2c + xerr), linewidth = 1) +
  geom_point(stat = 'identity', size = 8, alpha = 0.6) +
  ggforce::geom_mark_ellipse(aes(fill = HCAgrouping)) +
  scale_fill_manual(values = hcaColors) +
  scale_color_manual(values = hcaColors) +
  genTheme() +
  labs(x = 'Weighted N:C', y = 'Weighted P:C')

elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  ggplot(aes(HCAgrouping, weightP2c, color = HCAgrouping)) +
  geom_violin() +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  scale_fill_manual(values = hcaColors) +
  scale_color_manual(values = hcaColors) +
  genTheme() +
  labs(x = 'Site cluster', y = 'N:C per molecule')

elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  ggplot(aes(HCAgrouping, weightN2c, color = HCAgrouping)) +
  geom_violin() +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  scale_fill_manual(values = hcaColors) +
  scale_color_manual(values = hcaColors) +
  genTheme() +
  labs(x = 'Site cluster', y = 'N:C per molecule')
dev.off()

# "For instance, the energy potential gained during catabolism is greater in compounds with low C oxidation states" Wegley Kelly 2022 -- La rowe paper
# "We hypothesize that the production of more reduced, energy-rich compounds is one mechanism driving higher microbial growth rates"
# More reduced metabolomes could support the growth of faster growing microbial communities on coral tissue, leading to disbiosis of the coral metabolomes
# However these metabolites could also be intra-cellular meaning the corals are increasingly storing more reduced forms of carbon
# Surprsingly N:C does not vary similarly. We would expect N:C to increase with more reduced metabolites
# This could mean that there is more potential energy but it is more locked down in benzene rings and harder to degrade.
pdf('~/Documents/GitHub/greeneMaui/data/plots/poritesnosc.pdf', width = 15, height = 10)
elementalComposition%>%
  # group_by(SampleID, siteGroup, HCAgrouping)%>%
  # summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  # ungroup()%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  ggplot(aes(HCAgrouping, weightNOSC, fill = HCAgrouping)) +
  # geom_boxplot() +
  geom_violin() +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               colour = "white") +
  scale_fill_manual(values = hcaColors) +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  genTheme() +
  labs(x = 'HCA grouping', y = 'Weighted NOSC')

elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  ggplot(aes(HCAgrouping, weightGcox, fill = HCAgrouping)) +
  geom_violin() +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               colour = "white") +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  scale_fill_manual(values = hcaColors) +
  genTheme() +
  labs(x = 'HCA grouping', y = 'Weighted gCox')
dev.off()

# Stats -- NOSC LMER --------------------------------------------------------------
## There was a significant difference between HCA groups and average sample NOSC score
## p = 0.03002356 for Mont
## p = 2.429769e-07 for porites
elementalStats <- elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()

elementalStats%>%
  lmer(weightNOSC ~ HCAgrouping + (1|siteGroup), data = .)%>%
  car::Anova()%>%
  .[['Pr(>Chisq)']]

#Significant differences between each group
noscPvalues <- elementalStats%>%
  lmer(weightNOSC ~ HCAgrouping + (1|siteGroup), data = .)%>%
  multcomp::glht(., linfct = multcomp::mcp(HCAgrouping = "Tukey"))%>%
  summary(., test = multcomp::adjusted('BH'))%>%
  tidy()%>%
  mutate(test = 'NOSC')



# Stats -- Weighted N:C LMER ------------------------------------------------------
#  3.156069e-23 for mont
#  4.231315e-07 for por
elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  lmer(weightN2c ~ HCAgrouping + (1|siteGroup), data = .)%>%
  car::Anova()%>%
  .[['Pr(>Chisq)']]

ncPvalues <- elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  lmer(weightN2c ~ HCAgrouping + (1|siteGroup), data = .)%>%
  multcomp::glht(., linfct = multcomp::mcp(HCAgrouping = "Tukey"))%>%
  summary(., test = multcomp::adjusted('BH'))%>%
  tidy()%>%
  mutate(test = 'nc')

ncPvalues

# Stats -- Weighted P:C LMER ------------------------------------------------------
# No significant difference between the values for mont
#  5.332998e-13 for porites
elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  lmer(weightP2c ~ HCAgrouping + (1|siteGroup), data = .)%>%
  car::Anova()%>%
  .[['Pr(>Chisq)']]

pcPvalues <- elementalComposition%>%
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  lmer(weightP2c ~ HCAgrouping + (1|siteGroup), data = .)%>%
  multcomp::glht(., linfct = multcomp::mcp(HCAgrouping = "Tukey"))%>%
  summary(., test = multcomp::adjusted('BH'))%>%
  tidy()%>%
  mutate(test = 'pc')

pcPvalues


# Visualization -- elemental ratios vs.  environmental variables ----------
# Possibly delte this section  
# Maybe not helpful for interpretation of the data
elEnviro <- elementalComposition%>%
  filter(Site != 'Site-12')%>% # for Montipora
  group_by(SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  # group_by(siteGroup, HCAgrouping)%>%
  # mutate(n = 1,
  #        n = sum(n),
  #        serr = sd(weightN2c)/sqrt(n))%>%
  # summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  # ungroup()%>%
  select(siteGroup, HCAgrouping, weightN2c, weightP2c, weightNOSC)%>%
  pivot_longer(3:5, names_to = 'responseVar', values_to = 'values')%>%
  left_join(enviroVariability, by = c('siteGroup', 'HCAgrouping'))

# testLmData <- elEnviro%>%
#   filter(responseVar == 'weightN2c')
# 
# testLm <- lm(values ~ TNC.Phosphate.Benthic.log + TNC.Ammonia.Benthic.log + TNC.Silicate.Benthic.log + habPercent,
#    data = testLmData)
# 
# dredge(testLm)

n2cLm <- elEnviro%>%
  group_by(responseVar)%>%
  nest()%>%
  mutate(data = map(data, ~lmer(values ~ TNC.Phosphate.Benthic.log + TNC.Ammonia.Benthic.log + TNC.Silicate.Benthic.log + habPercent + (1|HCAgrouping),
                              data = .x)),
         r2 = map(data, ~ r.squaredGLMM(.x))
         )

# Porites conditional r2
# 0.5013139 nitrogen
# 0.4309536 NOSC
# 0.4890079 phosphorus 

# Montipora conditional r2
# 0.3724013 nitrogen

elEnviro%>%
  ggplot(aes(TNC.Phosphate.Benthic.log, values)) +
  geom_point() +
  facet_wrap(~responseVar, scales = 'free_y') +
  # geom_errorbar(aes(ymin = weightN2c - serr, ymax = weightN2c + serr)) +
  geom_smooth(method = 'lm') 

elEnviro%>%
  ggplot(aes(TNC.Ammonia.Benthic.log, values)) +
  geom_point() +
  facet_wrap(~responseVar, scales = 'free_y') +
  # geom_errorbar(aes(ymin = weightN2c - serr, ymax = weightN2c + serr)) +
  geom_smooth(method = 'lm') 

elEnviro%>%
  ggplot(aes(TNC.Silicate.Benthic.log, values)) +
  geom_point() +
  facet_wrap(~responseVar, scales = 'free_y') +
  # geom_errorbar(aes(ymin = weightN2c - serr, ymax = weightN2c + serr)) +
  geom_smooth(method = 'lm') 
  
elEnviro%>%
  ggplot(aes(habPercent, values)) +
  geom_point() +
  facet_wrap(~responseVar, scales = 'free_y') +
  # geom_errorbar(aes(ymin = weightN2c - serr, ymax = weightN2c + serr)) +
  geom_smooth(method = 'lm') 

# analysis -- p strucutres ------------------------------------------------
pContent <- elementalComposition%>%
  select(featureNumber, netVals, C:P)%>%
  unique()%>%
  filter(P > 0)%>%
  left_join(canopus%>%
              select(superclass, class, subclass, featureNumber), by = 'featureNumber')%>%
  left_join(libraryMatchClassyFire%>%
              select(-c(superclass, class, subclass)), by = 'featureNumber')%>%
  left_join(libraryMatchSources%>%
              select(`Compound_Name`, Source), by = 'Compound_Name')

length(pContent%>%
         filter(subclass %like% 'Glycerophospho%')%>%
         pull(P))/length(pContent$P)

# Visual -- Supplemental Figure S2/3  -- Significant compounds that contain library matches  --------
libraryBiclusterDf <- statsWdf%>%
  filter(netVals %in% hcaSigSubnets)%>%
  left_join(libraryMatchClassyFire, by = 'featureNumber')%>%
  left_join(libraryMatchSources, by = 'Compound_Name')%>%
  # select(-c(superclass, class, subclass))%>%
  # left_join(conciseRaw, by = c('network', 'featureNumber'))%>%
  filter(!is.na(Compound_Name),
         siteGroup != 'Site-12'
         )%>%
  select( 'Compound cleaned name', SampleID, Source, siteGroup, HCAgrouping, ra)%>%
  # select( 'Compound cleaned name', SampleID, Source, `Specific industry or broad chemical class`, siteGroup, HCAgrouping, ra)%>%
  # filter(Source == 'Cosmetic, agricultural, pharmaceutical or industrial product')%>%
  # mutate(Compound_Name = gsub('Spectral Match to ', '', Compound_Name),
  #        Compound_Name = gsub(' from NIST14', '', Compound_Name))%>%
  unite(sample, c('HCAgrouping', 'siteGroup'), sep = ';')%>%
  unite(annotation, c('Compound cleaned name', 'Source'), sep = '; ')%>%
  # unite(annotation, c('Compound cleaned name', 'Specific industry or broad chemical class'), sep = '; ')%>%
  group_by(sample, annotation, SampleID)%>% # Find the sum abundance of each subnetwork within each sample
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  # group_by(annotation)%>%
  # filter(max(ra)>0.01)%>%
  # ungroup()%>%
  group_by(sample, annotation)%>% # Find the average sum abundance of each subnetwork within each SITE
  summarize_if(is.numeric, mean)%>%
  # mutate(xic = log10(xic + 1))%>%
  # filter(max(ra) >= 0.01)%>%
  mutate(asin = asin(sqrt(ra)))%>%
  ungroup()%>%
  group_by(annotation)%>%
  filter(sum(ra) != 0)%>% # have to filter out any features that are not present in any sites. This is because we removed Site 12
  mutate(asin = zscore(asin))%>% # zscoring within the metabolite will allow for better visualization of metabolite variability across sites
  ungroup()%>%
  select(-ra)%>%
  pivot_wider(names_from = 'sample', values_from = 'asin')
  

rowColors <- libraryBiclusterDf%>%
  separate(annotation, c('Compound cleaned name', 'Source'), sep = '; ', remove = FALSE)%>%
  select(annotation, Source)%>%
  mutate(Source = case_when(Source %like any% c('Marine fungi%', 'MS mi%', 'Plant metabolite') ~ 'Contaminant, or other marine source', 
                                                TRUE ~ Source))%>%
  column_to_rownames(var = 'annotation')

colorList <- list(Source = c(
  `Contaminant, or other marine source` = '#EBEE2B',
  `Coral holobiont metabolite` = '#2B5CEE',
  `Coral holobiont stress response` = '#CF571F', 
  `Cosmetic, agricultural, pharmaceutical or industrial product` = '#8E51BD',
  # `NA` = 'black',
  `Unknown` = '#E3E3E3'
  ))
# colorList <- list(Source = colorList)

# colorList <- list(Source = c(
#   `Pharmaceutical` = '#EBEE2B',
#   `Industrial product` = '#2B5CEE',
#   `Agricultural/Pesticide` = '#CF571F', 
#   `Cosmetics` = '#8E51BD',
#   `Performance enhancer` = '#FFFFFF',
#   # `NA` = 'black',
#   `NA` = '#E3E3E3'
# ))


pdf('~/Documents/GitHub/greeneMaui/data/plots/biclusterLibraryMatches.pdf', width = 25, height = 35)
libraryBiclusterDf%>%
  separate(annotation, c('Compound cleaned name', 'Source'), sep = '; ')%>%
  select(-Source)%>%
  column_to_rownames(var = 'Compound cleaned name')%>%
  pheatmap::pheatmap(cluster_rows = TRUE,
                     fontsize_row = 10,
                     color = brewer.pal(n = 9, name = "Greys"),
                     annotation_row = rowColors,
                     annotation_colors = colorList)
dev.off()

# pdf('~/Documents/GitHub/greeneMaui/data/plots/anthropogenicContaminants.pdf', width = 25, height = 35)
# libraryBiclusterDf%>%
#   separate(annotation, c('Compound cleaned name', 'Source'), sep = '; ')%>%
#   select(-Source)%>%
#   column_to_rownames(var = 'Compound cleaned name')%>%
#   pheatmap::pheatmap(cluster_rows = TRUE,
#                      fontsize_row = 10,
#                      color = brewer.pal(n = 9, name = "Greys"),
#                      annotation_row = rowColors,
#                      annotation_colors = colorList)
# dev.off()

# Supplemental Table -- Unique library Matches ----------------------------
libMatchTable <- statsWdf%>%
  filter(netVals %in% hcaSigSubnets)%>%
  left_join(libraryMatchClassyFire, by = 'featureNumber')%>%
  # select(-c(superclass, class, subclass))%>%
  # left_join(conciseRaw, by = c('network', 'featureNumber'))%>%
  ungroup()%>%
  select(Compound_Name, superclass, class, subclass)%>%
  unique()%>%
  filter(!is.na(Compound_Name))

write_csv(libMatchTable, '~/Documents/GitHub/greeneMaui/data/analysis/libraryMatchTableRaw.csv')


# analysis Table -- Cytoscape export --------------------------------------
cytoTable <-  statsWdf%>%
  filter(netVals %in% hcaSigSubnets)%>%
  left_join(libraryMatchClassyFire, by = 'featureNumber')%>%
  # left_join(conciseRaw, by = c('network', 'featureNumber'))%>%
  ungroup()%>%
  select(Compound_Name, featureNumber, SampleID, siteGroup, HCAgrouping, ra)%>%
  group_by(featureNumber,Compound_Name, SampleID, siteGroup, HCAgrouping)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(featureNumber,Compound_Name, HCAgrouping)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  pivot_wider(names_from = 'HCAgrouping', values_from = 'ra')

write_csv(cytoTable, '~/Documents/GitHub/greeneMaui/data/analysis/cytoTable.csv')



# Visual -- library matches by source -- violin plots ---------------------
# Polyethylene glycols were all marked as coral metabolites because pentaethylene Glycol is produced by glycol
# Watson and Jones et al., 1977; Huang et al., 2005 are all 
librarySources <- statsWdf%>%
  filter(siteGroup != 'Site-12')%>%
  filter(netVals %in% hcaSigSubnets)%>%
  left_join(libraryMatchClassyFire, by = 'featureNumber')%>%
  left_join(libraryMatchSources, by = 'Compound_Name')%>%
  group_by(Source, HCAgrouping, siteGroup, SampleID)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  filter(Source %in% c('Coral holobiont metabolite', 'Cosmetic, agricultural, pharmaceutical or industrial product', 'Coral holobiont stress response'))%>%
  mutate(asin = asin(sqrt(ra)),
         log10 = log10(xic + 1))

pdf('~/Documents/GitHub/greeneMaui/data/plots/poritessourceViolins.pdf', width = 10, height = 30)
librarySources%>%
  ggplot(aes(HCAgrouping, ra, fill = HCAgrouping)) +
  geom_violin() +
  geom_point(stat = 'identity', alpha= 0.3, size = 3) +
  genTheme() +
  scale_fill_manual(values = hcaColors) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               colour = "white") +
  labs(x = 'HCA grouping', y = 'Sum feature intensity') +
  scale_y_log10() +
  facet_wrap(~Source, scales = 'free_y', nrow = 3)
dev.off()

# Stats -- library matches by source -- lmer ------------------------------
librarySourceLmer <- librarySources%>%
  group_by(Source)%>%
  nest()%>%
  mutate(lmer = map(data, ~lmer(asin ~ HCAgrouping + (1|siteGroup), data = .x)%>%
                      car::Anova()%>%
                      .[['Pr(>Chisq)']]),
         tukey = map(data, ~lmer(asin ~ HCAgrouping + (1|siteGroup), data = .x)%>%
                       multcomp::glht(., linfct = multcomp::mcp(HCAgrouping = "Tukey"))),
         tukeyPvalues = map(tukey, ~summary(., test = multcomp::adjusted('BH'))))%>%
  unnest(c(lmer))%>%
  ungroup()%>%
  mutate(lmerFDR = p.adjust(lmer, method = 'BH'))

sourcePvalues <- librarySourceLmer%>% 
  filter(lmerFDR < 0.05)%>%
  select(Source, tukeyPvalues)%>% 
  mutate(tukeyPvalues = map(tukeyPvalues, ~tidy(.x)))%>% 
  unnest(tukeyPvalues)


# Stats -- FDR correction on Tukey values ---------------------------------
statsCombined <- sourcePvalues%>%
  rename(test = Source)%>%
  bind_rows(shannonPvalues, 
            noscPvalues,
            pcPvalues, # only include for Porites
            ncPvalues)%>%
  mutate(fdr = p.adjust(p.value, method = 'BH'))
            

# Analysis -- source percent changes --------------------------------------
percentChangeSources <- librarySources%>%
  group_by(Source, HCAgrouping)%>%
  # mutate(sd = sd(ra))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(Source, HCAgrouping, ra)%>%
  filter(HCAgrouping != 2)%>%
  pivot_wider(names_from = 'HCAgrouping', values_from = 'ra')%>%
  mutate(difference = (`3` - `1`)/`1`)
  
# Porites percent changes:
# coral: -16.7% stress: 113.0% cosmetic: 32.8%

# Montipora percent changes:
# coral: -50.9% stress: 36.4% cosmetic: 26.0%


# Figure 3 and Stats --BoxPlot and LMER/ANOVA -- Relationship in HCA to Ra with superclasses      ----------------------
# For the ANOVA we need a dataframe that has summed the superclasses within each sample
anovaSuperclassesDf <- statsWdf%>%
  ungroup()%>%
  filter(netVals %in% hcaSigSubnets,
         siteGroup != 'Site-12'
         )%>%
  left_join(canopus%>%
              select(superclass, class, subclass, featureNumber), by = 'featureNumber')%>%
  select(network, superclass, SampleID, siteGroup, HCAgrouping, ra)%>%
  unite(sample, c('HCAgrouping', 'siteGroup'), sep = ';')%>%
  group_by(sample, superclass, SampleID)%>% # Find the sum abundance of each subnetwork within each sample
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  separate(sample, c('HCAgrouping', 'siteGroup'), sep = ';')
  # mutate(HCAgrouping = as.factor(HCAgrouping),
  #        HCAgrouping = fct_relevel(HCAgrouping, '1', '2', '4'))

## This linear model tests the sum abundance of each superclass between HCA groups
## Here we control for the site using it as a random effect
## from the mixed model Every superclass except lipids, nucleosides, and oxygen compounds had a significant relationship with HCA groups 
anovaSuperclasses <- anovaSuperclassesDf%>%
  mutate(asin = asin(sqrt(ra)))%>%
  group_by(superclass)%>%
  nest()%>%
  mutate(lmer = map(data, ~lmer(asin ~ HCAgrouping + (1|siteGroup), data = .x)%>%
                     car::Anova()%>%
                     .[['Pr(>Chisq)']]),
         # aov = map(data, ~ aov(asin~HCAgrouping, data = .x)%>%
         #             tidy()%>%
         #             filter(term != 'Residuals')%>%
         #             select(-c('term', 'df', 'sumsq', 'meansq', 'statistic'))))%>%
         tukey = map(data, ~lmer(asin ~ HCAgrouping + (1|siteGroup), data = .x)%>%
                       multcomp::glht(., linfct = multcomp::mcp(HCAgrouping = "Tukey"))),
         tukeyPvalues = map(tukey, ~summary(., test = multcomp::adjusted('BH'))))%>%
  unnest(c(lmer))%>%
  ungroup()%>%
  mutate(lmerFDR = p.adjust(lmer, method = 'BH'))%>%
  filter(lmerFDR < 0.05)

significantSuperclasses <- anovaSuperclasses%>%
  pull(superclass)
  
pdf('~/Documents/GitHub/greeneMaui/data/plots/MontiporasuperclassBoxplots.pdf', width = 15, height = 10)
anovaSuperclassesDf%>%
  filter(superclass %in% significantSuperclasses)%>%
  unite(sample, c('HCAgrouping', 'siteGroup'), sep = ' ')%>%
  group_by(sample, superclass)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  select(-c(network))%>%
  group_by(superclass)%>%
  mutate(ra = asin(sqrt(ra)),
         ra = zscore(ra))%>%
  pivot_wider(names_from = 'superclass', values_from = 'ra')%>%
  column_to_rownames(var = 'sample')%>%
  pheatmap::pheatmap(cluster_rows = FALSE, color = brewer.pal(n = 9, name = "BuGn"))
  # ggplot(aes(HCAgrouping, ra, fill = HCAgrouping)) +
  # geom_violin(color = 'white') +
  # geom_point(stat = 'identity', size = 2, alpha = 0.8) +
  # scale_fill_manual(values = hcaColors) +
  # facet_wrap(~superclass, scales = 'free_y') +
  # genTheme() +
  # labs(x = 'HCA group', y = 'Superclass sum relative abundance')
dev.off()  

# Stats -- linear model -- anthropogenic contamination DOES NOT explain coral cover, etc. ---------------
sourceVTNC <- librarySources%>%
  filter(Source == "Cosmetic, agricultural, pharmaceutical or industrial product")%>%
  select(siteGroup, TNC.Coral.Cover, TNC.Coral.Diversity, TNC.Coral.Disease, TNC.Herbivore.Biomass, log10)%>%
  group_by(siteGroup)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  left_join(coverChange%>%
              select(siteGroup, estimate),
            by = 'siteGroup')%>%
  filter(!is.na(estimate))

lm(estimate ~ log10, data = sourceVTNC)%>%
  summary()

lm(TNC.Coral.Disease ~ log10, data = sourceVTNC)%>%
  summary()

lm(TNC.Coral.Diversity ~ log10, data = sourceVTNC)%>%
  summary()

lm(TNC.Herbivore.Biomass ~ log10, data = sourceVTNC)%>%
  summary()
  

