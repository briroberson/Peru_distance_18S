# Load Packages ----
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(lme4)
library(car)
library(bbmle)
library(lmtest)
library(ape)
library(pairwiseAdonis)
library(LDM)
library(indicspecies)
library(MASS)
library(ecole)
library(plyr)
library(mirlyn)
library(ANCOMBC)
library(ggrepel)
library(patchwork)

# plotting theme

theme<-theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
             axis.title.y = element_text(face="bold", size = 18, margin = margin(r = 10)), 
             axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16),
             axis.title.x = element_text(size = 18, face = "bold",color = "black"),
             plot.margin = unit(c(0.1,0.1,0,0.1),"cm") ,
             strip.text = element_text(face = "bold", size = 16), 
             legend.position = "none")

dist_cols <- colorRampPalette(c("#9C6EB0", "cyan3"))(6)

# Load Metadata ----

###### 1. Load data. Alternative: once you run this once, you can then save it as an R file
#and load it directly. the code for this is in section 2d and 3e

### 1a. Metadata and elevation
#the metadata
metadata<-readr::read_tsv("peru-dec24-metadata18S.tsv")

# make distance a factor
metadata$distance<- as.factor(metadata$distance)

#this is the elevation file. this is already loaded into the 16S data but can also be loaded and joined separately
#waypoints<- read.csv("F:\\Research\\waypoints.csv")

### 1b. Other files, loaded into a phyloseq
# Load Other Data ----

#load it into a phyloseq object
phy <- qza_to_phyloseq("PeruDec24_18s_table.qza", 
                       "PeruDec24_18S_rooted-tree.qza", 
                       "PeruDec24_18S_taxonomy.qza")

#add the metadata to the phyloseq because for some reason I couldn't load it directly
#order the samples in the metadata to match the order of samples in phyloseq
metadata<-metadata[ order(match(metadata$sampleID, colnames(phy@otu_table))), ]
metadata$distance<- as.factor(metadata$distance)
#add the metadata to phyloseq
phy@sam_data<- sample_data(metadata)
#make the row names of the sample data match the actual sample ID
row.names(phy@sam_data)<- metadata$sampleID

#check that metadata didn't chop off names using sample_variables(phy). this should be the header names, not data
sample_variables(phy)


# Initial Processing ----
####### 2. Filtering

### 2a. format the heading names
## check if taxa need to be formatted
taxa_data<- as.data.frame(tax_table(phy)) #pull out taxa table
#looks good

### 2b. Filter out singletons 
pruned_phy<-prune_taxa(taxa_sums(phy)>1, phy) #this is similar to subset but it keeps only the taxa that had more than 1 occurence using the taxa sums function
pruned_phy
phy
#I printed them to compare 

### 2c. Filter out other things

#filter out positive and negative controls 
pruned_filtered_phy <-subset_samples(pruned_phy, !is.na(latrine)) 

#filter out NA at the Kingdom level
table(taxa_data$Kingdom)
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, Kingdom != "Unassigned")

#filter out archea and bacteria, plus mitochondria and plastids 
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !(Kingdom %in% c("Archaea", "Bacteria", "Eukaryota:mito", "Eukaryota:plas")))

#exclude vertebrates (family Craniata)
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !(Family %in% c("Craniata")))

#remove cephalapods
pruned_filtered_phy <- subset_taxa(pruned_filtered_phy, !(Genus %in% c('Cephalopoda')))

### 2d. Filter out the 3 samples that we chose to drop based on # of ASVs because 3 sites were sampled twice for quality
final_filtered_phy<-subset_samples(pruned_filtered_phy, 
                                   !sampleID %in% c("Dec-24-18S-30", 'Dec-24-18S-31', "Dec-24-18S-32") 
                                   &   !is.na(distance)) #this removes NAs (the pos and neg controls)
final_filtered_phy #print the two to compare
pruned_filtered_phy

#save it as R file so it can be easily loaded. at this point I recommend continuing
#through the rarefying step and save the rarefied file instead
saveRDS(final_filtered_phy, file="F:\\Research\\Dec24_18S\\final_filt_phy.rds") #use whatever file path for where you want to save it
final_filt_phy<-readRDS("F:\\Research\\Dec24_18S\\final_filt_phy.rds")


######## 3. Rarefying. if you want to skip to the rarefying and not do all of these steps (3a-3d) 
#that visualize the number of reads and determine the rarefy value, skip to section 3e

###3a. Plot histogram of number of reads per sample
reads_per_sample<- data.frame(sum=sample_sums(final_filtered_phy))
ggplot(reads_per_sample, aes(x=sum))+
  geom_histogram(binwidth=2500)


### 3b. Determine the minimum number of reads
smin<- min(sample_sums(final_filtered_phy))
smin #this was used in tutorials as the value to rarefy to but obviously 0 won't work


### 3c. Plot the rarefaction curves to determine sampling depth

#For all the data
otu.matrix = otu_table(final_filtered_phy) #make data into data frame
otu.matrix = as.data.frame(t(otu.matrix))

#plot
otu.rarecurve = rarecurve(otu.matrix, step = 50, label = F)
abline(v=14746)
abline(v=24869)
abline(v=12800)


### 3d. Do the rarefying using the chosen value. rngseed sets the seed for us within the function

#####rarefy using mirlyn
################do for 500 reps
#rarefy data
mirl_object_500<- mirl(final_filtered_phy, libsize=12800, set.seed=200, trimOTUs=T, replace=F, rep=500) #we chose this value with KBKs help. previously had tested lower and higher values and there was no difference so using higher

#make an empty object to put the ASV tables in
mirl_otu_500 <- vector("list", length(mirl_object_500))

#extract otu tables from each rarefied phyloseq and add to the empty object above
for (i in 1:length(mirl_object_500)){
  colnames(mirl_object_500[[i]]@otu_table) <- paste0(colnames(mirl_object_500[[i]]@otu_table))
  (mirl_otu_500[[i]] <- mirl_object_500[[i]]@otu_table)
}



#make metadata file with the correct samples (remove ones dropped during rarefying)
sample_id<- data.frame(final_filtered_phy@sam_data) 
sample_id$Samples<- row.names(sample_id)
sample_id<- sample_id %>% 
  filter(!Samples %in% c('Dec-24-18S-22','Dec-24-18S-23','Dec-24-18S-21','Dec-24-18S-6'))


sample_id <- sample_id$Samples

#make empty list for each sample
average_counts_500 <- vector("list", length(sample_id))

#give how many reps you will do
rep_500<-1:500
#make empty list to hold 5 dataframes
iter_list_500<- vector('list', length(rep_500))

#rewrite loop to select columns from each rep, then average them and put them in new otu table
for (i in 1:length(sample_id) ){
  for (j in rep_500){
    iter_list_500[[j]]<-dplyr::select(as.data.frame(mirl_otu_500[[j]]),i) #this selects each individual iteration's otu table and 
    iter_list_500[[j]]$ASVname<- row.names(iter_list_500[[j]])
  }
  
  sample_df_500<- reduce(iter_list_500[rep_500], full_join, by='ASVname')
  sample_df_500[is.na(sample_df_500)]<-0
  row.names(sample_df_500)<- sample_df_500$ASVname
  sample_df_500<- sample_df_500[,c(1, 3:(1+length(rep_500)))]
  sample_average_500 <- data.frame(rowMeans(sample_df_500))
  colnames(sample_average_500) <- sample_id[[i]]
  average_counts_500[[i]] <- sample_average_500
}
average_count_df_500 <- do.call(cbind, average_counts_500)

write.csv(x=average_count_df_500, file="D:\\Soil\\Dec24_18S\\500rep_averaged_OTUtable.csv")
write.csv(x=mirl_object_500, file="D:\\Soil\\Dec24_18S\\500rep_mirlobj.csv")

######## do some checks
#check that they all have the rarefied number of ASVs
# colSums(average_count_df_500)
# 
# #is this close to the number for the whole data frame?
# sum(iter_list_500[[1]]$`99`!=0)
# sum(average_count_df_500$`99` !=0)

#add to phyloseq
mirl_phyloseq <- final_filtered_phy
mirl_phyloseq@otu_table@.Data <- as.matrix(average_count_df_500)

rowSums(mirl_phyloseq@otu_table)==rowSums(average_count_df_500)  #should print a bunch of "TRUE"

#compare the two phyloseqs just to see and confirm that the expected number of samples are present
final_filtered_phy
mirl_phyloseq

#save to final phyloseq name used for analyses 
filt_rare_phy<- mirl_phyloseq

saveRDS(filt_rare_phy, file="D:\\Soil\\Dec24_18S\\filt_rare_phy_18s.rds") #use whatever file path for where you want to save it


##### Read in FINAL phyloseq----
filt_rare_phy<-readRDS("filt_rare_phy_18s_low.rds")

#order metadata to match phyloseq
metadata<-metadata[ order(match(metadata$sampleID, colnames(filt_rare_phy@otu_table))), ]




#### Alpha diversity 
#calculate shannon diversity
#all_shan_div<-estimate_richness(filt_rare_phy, measures='Shannon')
#add sample ID as column for the left join
#all_shan_div$sampleID<- row.names(all_shan_div)

#merge with the metadata so we can run a model and filter out the stuff already filtered out
#metadata_filt<- metadata %>% 
#  left_join(all_shan_div, by='sampleID') %>% 
#  filter(!is.na(Shannon))




# richness using specnumber() for non-integer counts from mirlyn rarefaction
otu_mat <- as(otu_table(filt_rare_phy), "matrix") 
if(taxa_are_rows(filt_rare_phy)) {otu_mat <- t(otu_mat)} #transpose

richness <- specnumber(otu_mat) #calculate richness
sample_ids <- sample_names(filt_rare_phy) #grab samp IDs
all_richness <- data.frame(`sampleID` = sample_ids, Observed = richness)

# merge with metadata
metadata_filt <- metadata %>% 
  left_join(all_richness, by='sampleID') %>%
  filter(!is.na(Observed))

#add elevation and format columns 
latrine_elevation <- read.csv("latrine_elevation.csv") %>%
  dplyr::select(latrine, elevation) %>%
  distinct(latrine, .keep_all = TRUE)

metadata_filt <- metadata_filt %>% 
  mutate(type=ifelse(grepl('L', latrine), 'latrine','veg_patch')) %>% 
  mutate(inside=ifelse(distance %in% c(1,2), 'inside', 'outside')) %>% 
  mutate(type_ins=paste(type, inside, sep='_')) %>%
  left_join(latrine_elevation, by = "latrine")
  

### Elevation models for manuscript

## Distance 
metadata_dist <- 
  metadata_filt %>%
  filter(
    type == 'latrine'
  ) %>%
  mutate(
    elevation_sc = 
      scale(elevation)
  ) %>% 
  mutate(
    dis_name = recode_factor(distance, '1'='Latrine center','2'='Latrine edge','3'='0.5m', '4'='1m', '5'='2m','6'='4m')) 

m_rich_dis <- lmer(Observed ~ distance + (1 | latrine), data = metadata_dist)
summary(m_rich_dis)
Anova(m_rich_dis)

m_rich_dis_elev_add <- lmer(Observed ~ distance + elevation_sc + (1 | latrine), data = metadata_dist, REML = FALSE)
summary(m_rich_dis_elev_add)
Anova(m_rich_dis_elev_add)

m_rich_dis_elev <- lmer(Observed ~ distance * elevation_sc + (1 | latrine), data = metadata_dist)
summary(m_rich_dis_elev)
Anova(m_rich_dis_elev)

ICtab(m_rich_dis, m_rich_dis_elev_add, m_rich_dis_elev)

# plot for manuscript 
ggplot(metadata_dist, aes(x = dis_name, y = Observed)) +
  geom_boxplot(
    aes(fill = dis_name),
    outliers = FALSE
  ) +
  geom_point(
    color = "black",
    alpha = 0.5,
    size = 3,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  scale_fill_manual(values = dist_cols, guide = "none") +
  labs(
    x = NULL,
    y = "Eukaryote ASV richness",
    tag = "(b)"
  ) +
  theme_bw() +
  theme +
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.01, 0.98),
    legend.position = "none")



### Elevation models for manuscript
metadata_filt <-
  metadata_filt %>%
  mutate(
    elevation_sc = 
      scale(elevation)
  ) 

m_rich_veg <- lmer(Observed ~ type_ins + (1 | latrine), data = metadata_filt)
summary(m_rich_veg)
Anova(m_rich_veg)

m_rich_veg_elev_add <- lmer(Observed ~ type_ins + elevation_sc + (1 | latrine), data = metadata_filt, REML = FALSE)
summary(m_rich_veg_elev_add)
Anova(m_rich_veg_elev_add)

m_rich_veg_elev <- lmer(Observed ~ type_ins * elevation_sc + (1 | latrine), data = metadata_filt)
summary(m_rich_veg_elev)
Anova(m_rich_veg_elev)

ICtab(m_rich_dis, m_rich_dis_elev_add, m_rich_dis_elev)

# plot for manuscript 
ggplot(metadata_filt, aes(type_ins, Observed)) +
  geom_boxplot(
    aes(fill = type_ins),
    outliers = FALSE
  ) +
  geom_point(
    color = "black",
    alpha = 0.5,
    size = 3,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  labs(
    x = NULL,
    y = "Eukaryote ASV richness",
    tag = "(b)"
  ) +
  scale_fill_manual(
    values = c(
      "latrine_inside" = "#9C6EB0",
      "veg_patch_inside" = "#4ACD66",
      "latrine_outside" = "#4DD7CE"
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c(
      "latrine_inside" = "Inside Latrine",
      "veg_patch_inside" = "Vegetation Patch",
      "latrine_outside" = "Outside Latrine"
    )
  ) +
  theme_bw() +
  theme +
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.01, 0.98),
    legend.position = "none"
  )

ggplot(metadata_filt, aes(x = elevation, y = Observed, color = type_ins)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    x = "Elevation (m)",
    y = "Eukaryote ASV richness",
    color = NULL,
    tag = "(c)"
  ) +
  scale_color_manual(
    values = c(
      "latrine_inside" = "#9C6EB0",
      "veg_patch_inside" = "#4ACD66",
      "latrine_outside" = "#4DD7CE"
    ),
    labels = c(
      "latrine_inside" = "Inside Latrine",
      "veg_patch_inside" = "Vegetation Patch",
      "latrine_outside" = "Outside Latrine"
    )
  ) +
  theme_bw() +
  theme +
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.01, 0.98),
    legend.position = "right"
  )


######### Beta diveristy
### NECESSARY Reroot the tree.----
# 5a. Reroot the tree
#It has to be binary but now it is not since we trimmed it
phy_tree<- phy_tree(filt_rare_phy) #put tree into an object
is.binary(phy_tree) #asking if it is binary. if false, go to next step

phy_tree(filt_rare_phy)<-multi2di(phy_tree) #fix the tree and put it back in the phyloseq
is.binary(phy_tree(filt_rare_phy)) #check if it's binary, should be true

## Permanova ----
### 5e. Permanova test----
#add these columns to the metadata in the phyloseq
sam_data<- filt_rare_phy@sam_data
sam_data<- sam_data %>% 
  mutate(type=factor(ifelse(grepl('L', latrine), 'latrine','veg_patch'))) %>% 
  mutate(inside=factor(ifelse(distance %in% c(1,2), 'inside', 'outside'))) %>% 
  mutate(type_ins=factor(paste(type, inside, sep='_')))

filt_rare_phy@sam_data<- sample_data(sam_data)

#subset for latrines
filt_phy_lat<-subset_samples(filt_rare_phy, type=='latrine')
metadata_lat<- metadata_filt %>% 
  filter(type=='latrine')


set.seed(200) ###VERY IMPORTANT, always keep the same

#######run permanova for distance
permanova_dis<- adonis2(distance(filt_phy_lat, method='wunifrac')~distance, data=metadata_lat, by='terms')
permanova_dis

#pairwise permanova 
permanova_pairwise(distance(filt_phy_lat, method='wunifrac'), grp=metadata_lat$distance, padj='holm')

#######permanova with distance and veg patch
metadata_filt$distance<- as.character((metadata_filt$distance))
metadata_filt$distance[c(25:27, 29:32)]<- 'veg'
metadata_filt$distance<- as.factor((metadata_filt$distance))
permanova_pairwise(distance(filt_rare_phy, method='wunifrac'), grp=metadata_filt$distance, padj='holm')


########## 2 category permanova
permanova_ins<- adonis2(distance(filt_phy_lat, method='wunifrac')~inside, data=metadata_lat, by='terms')
permanova_ins

#pairwise permanova 
permanova_pairwise(distance(filt_phy_lat, method='wunifrac'), grp=metadata_lat$inside, padj='holm')


######## 3 category with vegetation
permanova_veg<- adonis2(distance(filt_rare_phy, method='wunifrac')~type_ins, data=metadata_filt, by='terms')
permanova_veg

#pairwise permanova 
permanova_pairwise(distance(filt_rare_phy, method='wunifrac'), grp=metadata_filt$type_ins, padj='holm')


## Permanova with separate 1 & 2 groups 
metadata_filt_grouped <- metadata_filt %>%
  mutate(group = case_when(
    type == "veg_patch" ~ "veg_patch",
    type == "latrine" & distance == 1 ~ "1",
    type == "latrine" & distance == 2 ~ "2",
    type == "latrine" & distance %in% 3:6 ~ "outside",
    TRUE ~ NA_character_ ))
metadata_filt_grouped$group <- as.factor(metadata_filt_grouped$group)

permanova_grouped<- adonis2(distance(filt_rare_phy, method='wunifrac')~group, data=metadata_filt_grouped, by='terms')
permanova_grouped

#pairwise permanova 
permanova_pairwise(distance(filt_rare_phy, method='wunifrac'), grp=metadata_filt_grouped$group, padj='holm')


#### PCOA ----
## WUNIFRAC
#get asv table and transpose for just latrine samples
asvLat<- as.data.frame(otu_table(filt_phy_lat))
tasvLat <- data.frame(t(asvLat), check.names = F)

#calculate the pcoa
pcoaLat<-cmdscale(d=distance(filt_phy_lat, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorLat<-as.data.frame(wascores(x = pcoaLat$points, w = tasvLat))

#add the scores to the metadata
metadata_lat$axis01<- vegan::scores(pcoaLat)[,1]
metadata_lat$axis02<- vegan::scores(pcoaLat)[,2]

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hulls <- ddply(metadata_lat, "distance", find_hull)

#plot it for distance
ggplot(metadata_lat, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = distance, fill = distance), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = distance)) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

#plot it for inside outside
micro.hulls <- ddply(metadata_lat, "inside", find_hull)
ggplot(metadata_lat, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = inside, fill = inside), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = inside)) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

## for veg and latrine
#get asv table and transpose for all samples
asvAll<- as.data.frame(otu_table(filt_rare_phy))
tasvAll <- data.frame(t(asvAll), check.names = F)

#calculate the pcoa
pcoaAll<-cmdscale(d=distance(filt_rare_phy, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorAll<-as.data.frame(wascores(x = pcoaAll$points, w = tasvAll))

#add the scores to the metadata
metadata_filt$axis01<- vegan::scores(pcoaAll)[,1]
metadata_filt$axis02<- vegan::scores(pcoaAll)[,2]

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hulls <- ddply(metadata_filt, "type_ins", find_hull)

#plot it for distance
ggplot(metadata_filt, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = type_ins, fill = type_ins), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = type_ins)) +
  scale_color_manual(labels=c('Inside Latrine','Vegetation Patch','Outside Latrine'),
                     values=c('purple1','#74e374', 'cyan2'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))


# Simper ----

# run simper for Distance----
simper_dis<- simper(tasvLat, metadata_lat$distance)
simper_dis

#see the top 20
s_dis<- summary(simper_dis)
top10_1_2<-head(s_dis$`1_2`, n = 10)
top10_1_3<-head(s_dis$`1_3`, n = 10)
top10_1_4<-head(s_dis$`1_4`, n = 10)
top10_1_5<-head(s_dis$`1_5`, n = 10)

simpdis_asv1_2<- row.names(top10_1_2)
simpdis_asv1_3<- row.names(top10_1_3)
simpdis_asv1_4<- row.names(top10_1_4)
simpdis_asv1_5<- row.names(top10_1_5)

#get actual taxa info
taxa_dis <- as.data.frame(tax_table(filt_phy_lat)) #taxonomy
simper1_2_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_2,]
simper1_3_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_3,]
simper1_4_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_4,]
simper1_5_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_5,]


# run simper on Inside/Outside----
simper_ins<- simper(tasvLat, metadata_lat$inside)
simper_ins

#see the top 20
s_ins<- summary(simper_ins)
top10_ins<-head(s_ins$inside_outside, n = 10)

simpdis_asv_ins<- row.names(top10_ins)

#get actual taxa info
simper_ins_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv_ins,]


# run simper on 3 categories----
simper_veg<- simper(tasvAll, metadata_filt$type_ins)
simper_veg

#see the top 20
s_veg<- summary(simper_veg)
top10_veg_insL<-head(s_veg$latrine_inside_veg_patch_inside, n = 10)
top10_veg_ousL<-head(s_veg$latrine_outside_veg_patch_inside, n = 10)
top10_insL_ousL<-head(s_veg$latrine_inside_latrine_outside, n = 10)

simp_asv_veg_insL<- row.names(top10_veg_insL)
simp_asv_veg_ousL<- row.names(top10_veg_ousL)
simp_asv_insL_ousL<- row.names(top10_insL_ousL)

#get actual taxa info
simper_veg_insL_taxa<-taxa_dis[row.names(taxa_dis) %in% simp_asv_veg_insL,]
simper_veg_ousL_taxa<-taxa_dis[row.names(taxa_dis) %in% simp_asv_veg_ousL,]
simper_insL_ousL_taxa<-taxa_dis[row.names(taxa_dis) %in% simp_asv_insL_ousL,]

names(simper_veg)

#see only significant species
comparisons <- c("latrine_inside_latrine_outside", "latrine_inside_veg_patch_inside", "latrine_outside_veg_patch_inside") 

simper.results <- purrr::map_dfr(comparisons, function(comp) {
  df <- as.data.frame(simper_veg[[comp]]) %>%
    tibble::rownames_to_column("ASV")
  df %>%
    dplyr::mutate(
      Comparison = comp,
      Position = seq_len(nrow(df)))
})

#filter for significant 
sig_asvs_veg <- simper.results %>%
  filter(p <= 0.05) 

#create a df of significant ASVs with taxonomy 
taxaveg <- as.data.frame(tax_table(filt_rare_phy)) %>%
  tibble::rownames_to_column("ASV")
simper_taxaveg <- sig_asvs_veg %>%
  left_join(taxaveg, by = c("ASV" = "ASV"))

#grab top 10 overall only (choosing max avg within ASVs)
simper_taxaveg_top10all <- simper_taxaveg %>%
  group_by(ASV) %>%
  slice_max(
    order_by = average,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(desc(average)) %>%
  slice_head(n = 10) 

write.csv(simper_taxaveg_top10, "18S_simper_veg_top10.csv", row.names = FALSE)
saveRDS(simper_taxaveg_top5, file = "18s_simper_taxaveg_top5.rds")
simper_taxaveg_top5 <- readRDS("18s_simper_taxaveg_top5.rds")



###pcoa with simper labels 

## for veg and latrine
#get asv table and transpose for all samples
asvAll<- as.data.frame(otu_table(filt_rare_phy))
tasvAll <- data.frame(t(asvAll), check.names = F)

#calculate the pcoa
pcoaAll<-cmdscale(d=distance(filt_rare_phy, method='wunifrac'), eig=T)
#add the scores to the metadata
metadata_filt$axis01<- vegan::scores(pcoaAll)[,1]
metadata_filt$axis02<- vegan::scores(pcoaAll)[,2]

#retrieve species scores for it
spscorAll<-as.data.frame(wascores(x = pcoaAll$points, w = tasvAll))
#subset to top 10 significant asvs from simper to plot 
spscorAll$ASV <- rownames(spscorAll)
spscorAll_top10 <- spscorAll %>%
  dplyr::filter(ASV %in% simper_taxaveg_top10all$species)
#add taxonomy 
tax_df <- as.data.frame(phyloseq::tax_table(filt_rare_phy)) %>%
  tibble::rownames_to_column("ASV")
spscorAll_top10 <- spscorAll_top10 %>%
  dplyr::left_join(tax_df, by = "ASV")
#select lowest assigned taxonomy
spscorAll_top10 <- spscorAll_top10all %>%
  mutate(
    tax_label = case_when(
      !is.na(Species) & Species != "" &
        !Species %in% c("Embryophyceae_XX", "Microthamniales_X", "Rotifera_XX",
                        "Annelida_XX", "Chlamydomonadales_X", 	
                        "Chromadorea_X") ~ Species,
      !is.na(Genus) & Genus != "" &
        !Genus %in% c("Embryophyceae_X", "Heterotrichea_X",
                      "Annelida_X", "Rotifera_X") ~ Genus,
      !is.na(Family) & Family != "" ~ Family,
      !is.na(Order) & Order != "" ~ Order,
      !is.na(Class) & Class != "" ~ Class,
      !is.na(Phylum) & Phylum != "" ~ Phylum, TRUE ~ Kingdom,))

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hulls <- ddply(metadata_filt, "type_ins", find_hull)
#var explained labels 
var_exp <- round(100 * pcoaAll$eig / sum(pcoaAll$eig[pcoaAll$eig > 0]), 1)

#plot it 
ggplot(metadata_filt, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = type_ins, fill = type_ins), alpha = 0.4, linejoin = "round", show.legend = F) +
  geom_point(size = 3, aes(colour = type_ins)) +
  scale_color_manual(labels=c('Inside Latrine','Vegetation Patch','Outside Latrine'),
                     values=c("#9C6EB0","#4ACD66", "#4DD7CE")) + 
  scale_fill_manual(labels=c('Inside Latrine','Vegetation Patch','Outside Latrine'),
                    values=c("#9C6EB0","#4ACD66", "#4DD7CE")) + 
  geom_segment(aes(x=0, xend=V1, y=0, yend=V2), data=spscorAll_top10, arrow=arrow())+
  ggrepel::geom_text_repel(
    data = spscorAll_top10,
    aes(V1, V2, label = tax_label),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = NA)  + 
  xlab(paste0("PCoA 1 (", var_exp[1], "%)")) +
  ylab(paste0("PCoA 2 (", var_exp[2], "%)")) + 
  ggtitle("18S Beta Diversity") + 
  labs(color = NULL, fill = NULL) + 
  theme_bw()  +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black") , 
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm") ,
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5),
    legend.key = element_rect(
      colour = "black",
      fill = "white" )) +
  guides(
    colour = guide_legend(override.aes = list(size = 2.75)),
    fill = guide_legend(ncol = 3, override.aes = list(alpha = 1))) +
  labs(tag = "(b)") +
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98))


#export 1000 x 700


###export these taxa to BLAST

#read in sequences 
rep_seqs <- read_qza("PeruDec24_18s_rep-seqs.qza")$data
seq_df <- data.frame(
  ASV = names(rep_seqs),
  Sequence = as.character(rep_seqs),
  stringsAsFactors = FALSE)

# Extract PCoA axis scores for each ASV
axis_df <- spscorAll %>%
  dplyr::select(ASV, axis01 = V1, axis02 = V2)

# Join sequences and PCoA axis scores
taxaveg_top10_blast_18s <- simper_taxaveg_top10all %>%
  left_join(seq_df, by = "ASV") %>%
  left_join(axis_df, by = "ASV")

# Export
write.csv(
  taxaveg_top10_blast_18s,
  "taxaveg_blast_18s_TOP10_OVERALL.csv",
  row.names = FALSE
)






##### Differential Abundance ----
cat3_DA<-ancombc2(data = filt_rare_phy, tax_level = "Genus",
                 fix_formula = "type_ins", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.0, lib_cut = 0, s0_perc = 0.05,
                 group = "type_ins", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                     verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

cat3_DA_prim<- cat3_DA$res
saveRDS(cat3_DA_prim, file='F:\\Research\\Dec24_18S\\differentialAbundancePrim.rds')

#view pairwise
cat3_pair<- cat3_DA$res_pair %>% 
  mutate_if(is.numeric, function(x) round(x, 2))

saveRDS(cat3_pair, file='F:\\Research\\Dec24_18S\\differentialAbundancePair.rds')

# view Inside latrine and Vegetation Patch
res_insL_veg<- cat3_pair %>% 
  filter(q_type_insveg_patch_inside<.05 & passed_ss_type_insveg_patch_inside==T)

#inside and outside latrines
res_insL_ousL<- cat3_pair %>% 
  filter(q_type_inslatrine_outside<.05 & passed_ss_type_inslatrine_outside==T)

#outside latrines and veg patch
res_ousL_veg<- cat3_pair %>% 
  filter(q_type_insveg_patch_inside_type_inslatrine_outside<.05 & passed_ss_type_insveg_patch_inside_type_inslatrine_outside==T)
