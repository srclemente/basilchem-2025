# Multivariate analyses for basil VOC data 
  # created by/on: S Clemente / 7 October 2023
  # last updated by/on: S Clemente / 7 August 2026

# 1: Setup #### 
if(TRUE){
# opening packages
library(tidyverse) # for data organization, manipulation, & visualization
library(corrr) # for correlation plots
library(ggcorrplot) # also for correlation plots
library(vegan) # for multivariate analyses and visualizations
library(ggrepel) # for neat labels 
library(flextable) # for exporting into .docx tables
library(officer) # for declaring .docx format
library(rcompanion) # for getting compact letter display in SIMPER table
library(ggpubr) # for concatenating NMDS plots
library(readr) # for exporting into .csv files
library(factoextra) # for pca loadings


# set wd, upload data
chems <- read.csv('basil_chemicals_all.csv',
                  check.names = FALSE)
#re-sorting by group
chems <- chems[order(chems$type),]

#subsetting numerical data only
chems.n <- chems[,8:84]

#removing RIs from main chems dataset
chems <- chems[-1,]

#then, standardize everything so that the 

#removing unresolved compounds / contaminants
colnames(chems.n) # check!
chems.n <- chems.n[,-75] # removing 2,6,10,14-tetramethyl...pentadecane (unresolved)
chems.n <- chems.n[,-61] # removing benzenecycloheptane (contaminant)
chems.n <- chems.n[,-20] # removing cuminal (contaminant)
chems.n <- chems.n[,-1] # removing acetic acid butyl ester (contaminant)
colnames(chems.n) # check again~

# create vector for cultivar names 
vars <- c("Aroma 2", "Cinnamon", "Genovese", "Greek Dwarf", "Lemon",
          "Lettuce", "Lime", "Mammolo", "Dark Purple Opal", "Sweet Nufar")

# create vector for aroma names
aros <- c("Citrus", "Clove", "Eucalyptol","Pungent","Sweet Basil")

#creating a matrix to check for correlations
#log transforming data
chems.p <- chems.n[-1,]
min.val.check <- min(abs(chems.p[chems.p!=0])) # finds minimum value in dataset, divides by 10
chems.log.check <- log10((chems.p + sqrt(chems.p^2 + min.val.check^2))/10) # log transforming

corr.matrix.check <- cor(chems.log.check)
ggcorrplot(corr.matrix.check,
           tl.cex = 6)

# creating a smaller matrix to see later relationships
chems.z <- chems.log.check[,50:73]
corr.matrix.z <- cor(chems.z)
ggcorrplot(corr.matrix.z,
           tl.cex = 6)

if(TRUE){
#summing up very correlated peaks
#(use colnames(chems.n) to get column location)
#correlated peaks #1: epicubenol + tau_cadinol
chems.n$epic <- chems.n$Epicubenol + chems.n$`tau-Cadinol`
chems.n <- chems.n %>% relocate(epic, .before = Epicubenol) #reorder
chems.n <- chems.n[,-c(65,67)] # remove
colnames(chems.n)[64] <- "Epicubenol"
colnames(chems.n) ##check

#correlated peaks #2: Unknown 71/43 and Unknown 108/81
chems.n$Unknown_109_81_43 <- chems.n$`Unknown 109, 81, 43` + chems.n$`Unknown 71, 43, 79`
chems.n <- chems.n %>% relocate(`Unknown_109_81_43`, .before = `Unknown 71, 43, 79`) #reorder
chems.n <- chems.n[,-c(61,63)] # remove
colnames(chems.n)[60] <- "Unknown 109, 81, 43"
colnames(chems.n) ##check

#correlated peaks #3: b_caryophyllene
chems.n$`beta-Caryophyllene` <- chems.n$`β-Caryophyllene` + chems.n$`β-Ylangene`
chems.n <- chems.n %>% relocate(`beta-Caryophyllene`, .before = `β-Ylangene`) #reorder
chems.n <- chems.n[,-c(36:37)] # remove
colnames(chems.n)[35] <- "β-Caryophyllene"
colnames(chems.n) ##check

#correlated peaks #4: methyl eugenol
chems.n$`Methyl eugenol` <- chems.n$Methyleugenol + chems.n$`β-Elemene`
chems.n <- chems.n %>% relocate(`Methyl eugenol`, .before = Methyleugenol) #reorder new column
chems.n <- chems.n[,-c(29:30)] #remove old columns
colnames(chems.n) ##check

#correlated peaks #5: beta ocimene
chems.n$`B-Ocimene` <- chems.n$`β-Ocimene` + chems.n$Benzeneacetaldehyde #sum
chems.n <- chems.n %>% relocate(`B-Ocimene`, .before = `β-Ocimene`) #reorder
chems.n <- chems.n[,-c(7:8)] # remove
colnames(chems.n)[6] <- "β-Ocimene"
colnames(chems.n) ##check

#changing troublesome names
names(chems.n)[21] <- "Unknown 97, 111, 70, 71"
names(chems.n)[30] <- "Unknown 95, 161, 107, 93"
names(chems.n)[46] <- "Unknown 81, 105, 119, 161"
names(chems.n)[57] <- "Unknown 109, 81, 93, 95"
names(chems.n)[62] <- "Unknown 161, 204, 105, 162"
names(chems.n)[63] <- "Unknown 105, 204, 71, 81"
names(chems.n)[65] <- "Unknown 93, 107, 81, 91"
names(chems.n)[66] <- "Unknown 95, 97, 81, 110"
names(chems.n)[67] <- "Unknown 161, 204, 91, 105"
names(chems.n)[68] <- "Unknown 161, 204, 95, 81"
names(chems.n)

# standardizing names
names(chems.n)[6] <- "(Z)-β-Ocimene"
names(chems.n)[8] <- "(Z)-linalool oxide"
names(chems.n)[16] <- "(E)-geraniol"
names(chems.n)[17] <- "(Z)-geraniol"
names(chems.n)[29] <- "(Z)-ɑ-Bergamotene"
names(chems.n)[34] <- "(E)-ɑ-Bergamotene"
names(chems.n)[41] <- "epi-Bicyclosesquiphellandrene"
names(chems.n)[54] <-"[S-(Z)]- 3,7,11-trimethyl-1,6,10-Dodecatrien-3-ol"


names(chems.n)

#reassign compound RIs
chems.n[1,60] <- 1640
chems.n[1,57] <- 1618
chems.n[1,33] <- 1444
chems.n[1,27] <- 1406
chems.n[1,6]  <- 1048
}

#extract RIs
RIs <- (chems.n[1,])
chems.n <- chems.n[-1,]
}

# 2: Flowers v. Leaves ####
# 2.1 difference tests 
# create distance matrix
dist.lvf <- vegdist(chems.n, method = "bray")

#dispersion test
(disptest.lvf <- anova(betadisper(dist.lvf, chems$type)))
#statistically significant. PERMANOVA preferred over ANOSIM, indicate in paper
  # see: http://doi.org/10.1890/12-2010.1

# assessing statistical differences between groups
(basil.perm <- adonis2(dist.lvf ~ chems$type,
                      method = "bray"))
 # significantly different!

# 2.2 SIMPER 
simper.lvf <- simper(comm = chems.n, 
                   group = chems$type,
                   permutations = 10000)
simper.output <- summary(simper.lvf)
simper.output <- simper.output$flower_leaf

# Row names as a column
simper.output$`Compound Name` <- rownames(simper.output)
rownames(simper.output) <- NULL

# Make trait/compound names first column
simper.output <- simper.output %>% 
  relocate(`Compound Name`, 
           .before = average)

# renaming
colnames(simper.output) <- c("Compound Name", "Average", "SD", "Ratio",
                             "Avg. Flower", "Avg. Leaf", "Cum. Sum", "p")


# 2.3 NMDS
# first standardize data 
dissim <- chems.n %>%
  decostand(method = "hellinger")

# set seed for reproduciblity. 
  # different seeds will simply change the orientation of polygons.
  # but polygon shape and distances will remain the same.
set.seed(2025)
nmds.lvf <- metaMDS(dissim, 
                    autotransform = FALSE,
                    distance = "bray",
                    k = 7,
                    model = "global",
                    try = 30,
                    trymax = 200)

# k = 7 allows for n > 2k+1 for both leaf & flower data
  # while minimizing stress. Although stress is >0.1 from k=3 onwards

# Get scores
nmds.lvf.scores <- as_tibble(scores(nmds.lvf)$sites)

# Add sample type
nmds.lvf.scores <- nmds.lvf.scores %>% 
  mutate(type = chems$type)

# Identify extreme points of convex hull
# convex hull values - flowers
hull.flowers<- 
  nmds.lvf.scores[nmds.lvf.scores$type == "flower", ][chull(nmds.lvf.scores[
    nmds.lvf.scores$type == "flower", c("NMDS1", "NMDS2")]), ] 

# convex hull values - leaves
hull.leaves <- 
  nmds.lvf.scores[nmds.lvf.scores$type == "leaf", ][chull(nmds.lvf.scores[
    nmds.lvf.scores$type == "leaf", c("NMDS1", "NMDS2")]), ]  

# Combine
hull.scores <- rbind(hull.flowers, 
                       hull.leaves)  

# get centroids 
centroids1 <- nmds.lvf.scores %>% 
                group_by(type) %>%
                  summarize(NMDS1a = mean(NMDS1),
                            NMDS2a = mean(NMDS2)) 
centroids1$NMDS1a <- centroids1$NMDS1a + c(0.5, 0)
centroids1$NMDS2a <- centroids1$NMDS2a + c(-0.41, 0)

names(centroids1) <- c("type","NMDS1", "NMDS2")

# Plot - IGNORE ggrepel warnings
plot.nmds.lvf <- ggplot(nmds.lvf.scores) + 
  geom_point(aes(x = NMDS1, # points
                 y = NMDS2, 
                 color = type, 
                 shape = type), 
             size = 3) + 
  geom_polygon(data = hull.scores, # shaded area
               aes(x = NMDS1,
                   y = NMDS2, 
                   fill = type, 
                   col = type),
               alpha = 0.3) + 
  scale_color_manual(name = "Type",
                     values = c("#93B7BE", "#554348"), 
                     labels = c("flower", "leaf")) +
  scale_fill_manual(name = "Type",
                    values = c("#93B7BE", "#554348"),
                    labels = c("flower", "leaf")) +
  scale_shape_manual(name = "Type",
                     values = c(18,19),
                     labels = c("flower", "leaf")) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") 

# adding labels
plot.nmds.lvf <-
  plot.nmds.lvf + geom_text_repel(data = centroids1, 
                                mapping = aes(x = NMDS1,
                                              y = NMDS2,
                                              label = c("flower", 
                                                        "leaf")),
                                fontface = c("bold"),
                                color = c("#93B7BE", "#554348"),
                                min.segment.length = 10,
                                direction = "both",
                                size = 8,
                                box.padding = 0,
                                point.padding = 0,
                                max.overlaps = 10) +
  annotate(
    "text", label = paste("Stress = ", 
                          round(nmds.lvf$stress, 
                                digits = 3)),
    x = 0.45, y = -1.121, size = 8, fontface = "italic",
    color = "#ffffff")

plot.nmds.lvf
  
# 3 Varieties ####
# 3.1: Setting Up
# creating a flower-only dataset
flowers <- chems.n[c(1:37),]

#re-adding important information
chems.info <- chems[c(1:37),c(1:5)]
flowers <- cbind(chems.info, flowers)

# rearranging by variety
flowers <- flowers[order(flowers$variety),]

# create distance matrix
dist.f <- flowers %>% 
  select(where(~is.numeric(.x))) %>% 
    vegdist(method = "bray")

# 3.2 Difference tests
# dispersion test
(disptest.f <- anova(betadisper(dist.f, flowers$variety)))
#ns - ANOSIM and PERMANOVA OK

# anosim
(anosim19 <-  anosim(dist.f, flowers$variety, permutations = 10000))
# sig dif

# permanova
(basil.perm.f <- adonis2(dist.f ~ flowers$variety, method = "bray"))
#sig dif

# 3.3 SIMPER 
# WARNING - computationally taxing. leave as TRUE if not in a rush.
if(TRUE){
simper.f <- simper(comm = flowers %>% 
                     select(where(~is.numeric(.x))), 
                     group = flowers$variety,
                     permutations = 10000)
simper.f.output <- summary(simper.f)

# extracting just p-values from large simper object
#taken from Konrad Rudolph at
  #https://stackoverflow.com/questions/23758858/how-can-i-extract-elements-from-lists-of-lists-in-r
p.vals <- lapply(simper.f.output, `[[`, "p") #added an extra "["
p.table <- as.data.frame(do.call(rbind, p.vals))
p.table <- p.table %>% 
            mutate_all(~replace(., is.na(.), 0.9999)) #replacing NA with 1
            names(p.table) <- names(chems.n)

# Create "ID" indicating pairwise comparisons
p.table$Comparison <- gsub("_", "-", names(p.vals))

# Make trait/compound names first column
p.table <- p.table %>% 
  relocate(`Comparison`, 
           .before = `Sabinene`)

# now get tukey letters for every chemical 
p.list <- as.list(p.table[2:69])
t.lets <- lapply(p.list, function(x)
                cldList(p.value = x, 
                        comparison = p.table$Comparison))
t.list <- lapply(t.lets, `[[`, "Letter")
f.simper.cld <- as.data.frame(do.call(rbind, t.list))
names(f.simper.cld ) <- vars
f.simper.cld $Compound <- rownames(f.simper.cld )

rownames(f.simper.cld) <- NULL

f.simper.cld <- f.simper.cld %>% 
  relocate(`Compound`, 
           .before = `Aroma 2`)
}

# 3.4 NMDS: VARIETIES
# first standardize data 
dissim.f <- flowers %>%
  select(where(~is.numeric(.))) %>%
    decostand(method = "hellinger")

# run nmds
set.seed(2025) # so results are reproducible
nmds.f <- metaMDS(dissim.f, 
                    autotransform = FALSE,
                    distance = "bray",
                    k = 3,
                    model = "global",
                    try = 30,
                    trymax = 200)

# Get scores
nmds.f.scores <- as_tibble(scores(nmds.f)$sites)

# Add variety
nmds.f.scores <- nmds.f.scores %>% 
  mutate(variety = flowers$variety)

# Identify extreme points of convex hull for each variety
hull.a2 <- 
  nmds.f.scores[nmds.f.scores$variety == "aroma 2", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "aroma 2", c("NMDS1", "NMDS2")]), ] 

hull.cn <- 
  nmds.f.scores[nmds.f.scores$variety == "cinnamon", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "cinnamon", c("NMDS1", "NMDS2")]), ] 

hull.gv <- 
  nmds.f.scores[nmds.f.scores$variety == "genovese", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "genovese", c("NMDS1", "NMDS2")]), ] 

hull.gd <- 
  nmds.f.scores[nmds.f.scores$variety == "grk dwarf", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "grk dwarf", c("NMDS1", "NMDS2")]), ] 

hull.lm <- 
  nmds.f.scores[nmds.f.scores$variety == "lemon", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "lemon", c("NMDS1", "NMDS2")]), ] 

hull.lt <- 
  nmds.f.scores[nmds.f.scores$variety == "lettuce", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "lettuce", c("NMDS1", "NMDS2")]), ] 

hull.li <- 
  nmds.f.scores[nmds.f.scores$variety == "lime", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "lime", c("NMDS1", "NMDS2")]), ] 

hull.mm <- 
  nmds.f.scores[nmds.f.scores$variety == "mammolo", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "mammolo", c("NMDS1", "NMDS2")]), ] 

hull.po <- 
  nmds.f.scores[nmds.f.scores$variety == "purp opal", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "purp opal", c("NMDS1", "NMDS2")]), ] 

hull.sn <- 
  nmds.f.scores[nmds.f.scores$variety == "swt nufar", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "swt nufar", c("NMDS1", "NMDS2")]), ] 


# Combine
hullf.scores <- rbind(hull.a2, hull.cn, hull.gv, hull.gd, hull.lm, 
                     hull.lt, hull.li, hull.mm, hull.po, hull.sn)

# Plot
#necessary values for plot
cols <- c("#781c81", "#43328d", "#99bd5c", "#d92120","#70b484", 
          "#416fb8", "#c3ba45", "#e0a239", "#e66b2d","#519cb8")
null.cols <- c("#ffffff", "#43328d", "#ffffff", "#ffffff","#70b484", 
          "#ffffff", "#c3ba45", "#ffffff", "#ffffff","#000000")
null.vars <- c("All others", "Cinnamon", "All others", "All others", "Lemon",
           "All others", "Lime", "All others", "All others", "All others")
centroids <- nmds.f.scores %>%
  group_by(variety) %>%
  summarize(NMDS1 = mean(NMDS1),
            NMDS2 = mean(NMDS2))

# plotting - IGNORE ggrepel warnings
if(TRUE){
plot.nmds.f <- ggplot(nmds.f.scores) + 
  geom_point(aes(x = NMDS1, # points
                 y = NMDS2, 
                 color = variety), 
             size = 3) + 
  geom_polygon(data = hullf.scores, # shaded area
               aes(x = NMDS1,
                   y = NMDS2, 
                   fill = variety, 
                   col = variety),
               alpha = 0.3) + 
  scale_color_manual(name = "Variety",
                     values = cols, 
                     labels = vars) +
  scale_fill_manual(name = "Variety",
                    values = cols,
                    labels = vars) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") 

# adding labels
plot.nmds.f <-
  plot.nmds.f + geom_text_repel(data = centroids, 
                                mapping = aes(x = NMDS1,
                                              y = NMDS2,
                                             label = null.vars),
                                xlim = c(-Inf, NA),
                                fontface = c("italic", "bold", "italic","italic","bold",
                                             "italic","bold","italic","italic","italic"),
                                nudge_x = c(0,0.6,0,0,0.1,0,0,0,0,0.6),
                                nudge_y = c(0,0.5,0,0,0.1,0,0.5,0,0,0.1),
                                color = null.cols,
                                min.segment.length = 3,
                                direction = "both",
                                size = c(1,8,1,1,8,1,8,1,1,8),
                                box.padding = 2.5,
                                point.padding = 0.5,
                                max.overlaps = 10) +
              annotate("text", 
                         label = paste("S = ",round(nmds.f$stress, 
                              digits = 2)),
                         x = 1, 
                         y = -1, 
                         size = 8, 
                         fontface = "italic",
                         color = "#000000") 

plot.nmds.f
}

# 4: Additional Plots & Analyses for Reviewer Comments #### 
# 4.1 NMDS C: AROMA PROFILES
# dispersion test
(disptest.aro<- anova(betadisper(dist.f, flowers$flavor)))
#ns p = 0.318- ANOSIM and PERMANOVA OK

# anosim
(anosim19a <-  anosim(dist.f, flowers$flavor, permutations = 10000))
#sig dif p = 0.002

# permanova
(basil.perm.aro <- adonis2(dist.f ~ flowers$flavor, method = "bray"))
#sig dif p = 0.002

# Recall general NMDS scores
nmds.ap.scores <- as_tibble(scores(nmds.f)$sites)

# Add aroma profile
nmds.ap.scores <- nmds.ap.scores %>% 
  mutate(aroma = flowers$flavor)
# Identify extreme points of convex hull for each variety
hull.ci <- 
  nmds.ap.scores[nmds.ap.scores$aroma == "citrus", ][chull(nmds.ap.scores[
    nmds.ap.scores$aroma == "citrus", c("NMDS1", "NMDS2")]), ] 

hull.cl <- 
  nmds.ap.scores[nmds.ap.scores$aroma == "clove", ][chull(nmds.ap.scores[
    nmds.ap.scores$aroma == "clove", c("NMDS1", "NMDS2")]), ] 

hull.eu <- 
  nmds.ap.scores[nmds.ap.scores$aroma == "eucalyptol", ][chull(nmds.ap.scores[
    nmds.ap.scores$aroma == "eucalyptol", c("NMDS1", "NMDS2")]), ] 

hull.pu <- 
  nmds.ap.scores[nmds.ap.scores$aroma == "pungent", ][chull(nmds.ap.scores[
    nmds.ap.scores$aroma == "pungent", c("NMDS1", "NMDS2")]), ] 

hull.sb <- 
  nmds.ap.scores[nmds.ap.scores$aroma == "swt bas", ][chull(nmds.ap.scores[
    nmds.ap.scores$aroma == "swt bas", c("NMDS1", "NMDS2")]), ] 

# Combine
hulla.scores <- rbind(hull.ci, hull.cl, hull.eu, hull.pu, hull.sb)

# Plot
#necessary values for plot
col.a <- c("#ff0029", "#377eb8", "#66a61e", "#984ea3","#00d2d5")

centroids.a <- nmds.ap.scores %>%
  group_by(aroma) %>%
  summarize(NMDS1 = mean(NMDS1),
            NMDS2 = mean(NMDS2))



# plotting - IGNORE ggrepel warnings
if(TRUE){
  plot.nmds.a <- ggplot(nmds.ap.scores) + 
    geom_point(aes(x = NMDS1, # points
                   y = NMDS2, 
                   color = aroma), 
               size = 3) + 
    geom_polygon(data = hulla.scores, # shaded area
                 aes(x = NMDS1,
                     y = NMDS2, 
                     fill = aroma, 
                     col = aroma),
                 alpha = 0.3) + 
    scale_color_manual(name = "Variety",
                       values = col.a, 
                       labels = aros) +
    scale_fill_manual(name = "Variety",
                      values = col.a,
                      labels = aros) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none") 
  
  # adding labels
  plot.nmds.a <-
    plot.nmds.a + geom_text_repel(data = centroids.a, 
                                  mapping = aes(x = NMDS1,
                                                y = NMDS2,
                                                label = aros),
                                  xlim = c(0, NA),
                                  fontface = c("bold", "bold", "bold","bold","bold"),
                                  nudge_x = c(0,0,0,0,0),
                                  nudge_y = c(0.4,0,0,0.7,0),
                                  color = col.a,
                                  min.segment.length = 2,
                                  direction = "both",
                                  size = 8,
                                  box.padding = 3,
                                  point.padding = 2.5,
                                  max.overlaps = 10) +
    annotate("text", 
             label = paste("S = ",round(nmds.f$stress, 
                                             digits = 2)),
             x = 1, 
             y = -1, 
             size = 8, 
             fontface = "italic",
             color = "#000000") 
  
  plot.nmds.a
}

# 4.2 NMDS C: REMOVE LEMON, LIME, & CINNAMON
# taking flower dataset and removing varieties of interest
flowers.red <- flowers[flowers$variety != "lemon" ,]
flowers.red <- flowers.red[flowers.red$variety != "lime" ,]
flowers.red <- flowers.red[flowers.red$variety != "cinnamon",]

# create distance matrix
dist.fred <- flowers.red %>% 
  select(where(~is.numeric(.x))) %>% 
  vegdist(method = "bray")

# dispersion test
(disptest.fred <- anova(betadisper(dist.fred, flowers.red$variety)))
#ns p = 0.867- ANOSIM and PERMANOVA OK

# anosim
(anosim19 <-  anosim(dist.fred, flowers.red$variety, permutations = 10000))
# sig dif p = 0.017

# permanova
(basil.perm.fred <- adonis2(dist.fred ~ flowers.red$variety, method = "bray"))
#sig dif p = 0.003


# first standardize data 
dissim.fred <- flowers.red %>%
  select(where(~is.numeric(.))) %>%
  decostand(method = "hellinger")

# run nmds
set.seed(2025) # so results are reproducible
nmds.fred <- metaMDS(dissim.fred, 
                  autotransform = FALSE,
                  distance = "bray",
                  k = 3,
                  model = "global",
                  try = 30,
                  trymax = 200)

# Get scores
nmds.fred.scores <- as_tibble(scores(nmds.fred)$sites)

# Add variety
nmds.fred.scores <- nmds.fred.scores %>% 
  mutate(variety = flowers.red$variety)

# Identify extreme points of convex hull for each variety
hullf.a2 <- 
  nmds.fred.scores[nmds.fred.scores$variety == "aroma 2", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "aroma 2", c("NMDS1", "NMDS2")]), ] 

hullf.gv <- 
  nmds.fred.scores[nmds.fred.scores$variety == "genovese", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "genovese", c("NMDS1", "NMDS2")]), ] 

hullf.gd <- 
  nmds.fred.scores[nmds.fred.scores$variety == "grk dwarf", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "grk dwarf", c("NMDS1", "NMDS2")]), ] 

hullf.lt <- 
  nmds.fred.scores[nmds.fred.scores$variety == "lettuce", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "lettuce", c("NMDS1", "NMDS2")]), ] 

hullf.mm <- 
  nmds.fred.scores[nmds.fred.scores$variety == "mammolo", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "mammolo", c("NMDS1", "NMDS2")]), ] 

hullf.po <- 
  nmds.fred.scores[nmds.fred.scores$variety == "purp opal", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "purp opal", c("NMDS1", "NMDS2")]), ] 

hullf.sn <- 
  nmds.fred.scores[nmds.fred.scores$variety == "swt nufar", ][chull(nmds.fred.scores[
    nmds.fred.scores$variety == "swt nufar", c("NMDS1", "NMDS2")]), ] 


# Combine
hullfred.scores <- rbind(hullf.a2, hullf.gv, hullf.gd, 
                      hullf.lt, hullf.mm, hullf.po, hullf.sn)

# Plot
#necessary values for plot
cols.fr <- c("#781c81", "#99bd5c", "#d92120",
          "#416fb8", "#e0a239", "#e66b2d","#519cb8")
vars.fr <- c("Aroma 2", "Genovese", "G. Dwarf", "Lettuce Lf.", 
             "Mammolo", "P. Opal", "S. Nufar")
centroids.fr <- nmds.fred.scores %>%
  group_by(variety) %>%
  summarize(NMDS1 = mean(NMDS1),
            NMDS2 = mean(NMDS2))

# plotting - IGNORE ggrepel warnings
if(TRUE){
  plot.nmds.fr <- ggplot(nmds.fred.scores) + 
    geom_point(aes(x = NMDS1, # points
                   y = NMDS2, 
                   color = variety), 
               size = 3) + 
    geom_polygon(data = hullfred.scores, # shaded area
                 aes(x = NMDS1,
                     y = NMDS2, 
                     fill = variety, 
                     col = variety),
                 alpha = 0.3) + 
    scale_color_manual(name = "Variety",
                       values = cols.fr, 
                       labels = vars.fr) +
    scale_fill_manual(name = "Variety",
                      values = cols.fr,
                      labels = vars.fr) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none") 
  
  #adjusting axes
  plot.nmds.fr <- 
      plot.nmds.fr +
      xlim(-1, 1) +
      ylim(-0.6, 0.6) 
  
  # adding labels
  plot.nmds.fr <-
    plot.nmds.fr + geom_text_repel(data = centroids.fr, 
                                  mapping = aes(x = NMDS1,
                                                y = NMDS2,
                                                label = vars.fr),
                                  #xlim = c(-1, NA),
                                  fontface = "bold",
                                  nudge_x = c(-0.15,0,0.2,0,0.6,-0.5,0),
                                  nudge_y = c(0.3,0.2,-0.1,-0.35,0.4,0.0,-0.1), 
                                  color = cols.fr,
                                  min.segment.length = 2,
                                  direction = "both",
                                  size = 8,
                                  box.padding = 1,
                                  point.padding = 1,
                                  max.overlaps = 10) +
    annotate("text", 
             label = paste("S = ",round(nmds.fred$stress, 
                                             digits = 2)),
             x = 0.55, 
             y = -0.6, 
             size = 8, 
             fontface = "italic",
             color = "#000000") 
  
  plot.nmds.fr
}


# 5: Creating Tables, Exporting Tables & Figs ####
# 5.1 TABLES 2 / A1: ALL FLORAL CHEMISTRY BY CULTIVAR
  # FOR TABLE 2, RUN AS IS.
  # TO MAKE TABLE A1, DON'T RUN LINES THAT ARE INDICATED
    # and be sure to rename .docx filenames at all export prompts 

# get means for each variety
fleurs <- flowers[,-c(1:2,4:5)]
flowers.bvo <- fleurs %>% 
  group_by(variety) %>%
  summarize_all(list(average = mean))
flowers.bvar <- flowers.bvo[,-1]
names(flowers.bvar) <- names(chems.n)
flowers.bv <- rbind(RIs, flowers.bvar )

# keep only values that are at least 1% of one variety
flower.sums <- rowSums(flowers.bv) # calculate total chemicals per variety
# create a logical df where TRUE means a chemical is at least 1%
flower.logic <- as.data.frame((apply(flowers.bv, 
                     2, function(x) x/flower.sums)) > 0.01)
# replace logic with numerical values
flower.logic[flower.logic == FALSE] <- 0
flower.logic[flower.logic == TRUE] <- 1
flower.logic[1,] <- rep(1, times = ncol(flower.logic)) # so RIs are always 1

# index original dataset so that only columns with at least 
  # one "1" in flower.logic (except the RIs, which are always 1) are kept
  # TO CREATE TABLE A2, DON'T RUN NEXT LINE, proceed with rest. 
flowers.bv<- flowers.bv[colSums(flower.logic) > 1] 

#extract RIs of retained compounds, remove from dataset for now.
retainedRIs <- flowers.bv[1,]
flowers.bv <- flowers.bv[-1,]

# replace zeroes with NA for next steps
flowers.bv <- flowers.bv %>%
  select(where(~is.numeric(.x))) %>%
  round(digits = 0)
flowers.bv[flowers.bv < 20] <- NA

# get standard errors for each variety
flowers.sd <- fleurs %>%
  group_by(variety) %>%
   summarize_all(list(se = sd))
    
samples <- fleurs %>%
  group_by(variety) %>%
  count()

flowers.se <- flowers.sd %>%
  select(where(~is.numeric(.x))) %>%
    mutate(./samples$n) %>%
     round(digits = 0)
names(flowers.se) <- names(chems.n)
flowers.se.full <- flowers.se

# keep only se values for remaining chemicals in flowers.bv
# DON'T RUN NEXT TWO LINES IF YOU WANT TABLE A2
flowers.se<- flowers.se[colSums(flower.logic[-1,]) > 0] 
names(flowers.bv) == names(flowers.se) #check

#convert numbers to characters
flowers.bv <- flowers.bv %>%
  select(where(~is.numeric(.x))) %>%
  mutate(across(everything(), as.character))

flowers.se <- flowers.se %>%
  select(where(~is.numeric(.x))) %>%
  mutate(across(everything(), as.character))

# merge dataframes together
# first convert mean and sd to matrices
flowers.bv <- as.matrix(flowers.bv)
flowers.se <- as.matrix(flowers.se)

# then merge 
  # (thanks to Aaron at:https://stackoverflow.com/questions/6408016/paste-together-two-dataframes-element-by-element)

flowers.mse <- as.data.frame(matrix(paste(flowers.bv, flowers.se, sep = " ± "),
                                    nrow = nrow(flowers.bv), 
                                    dimnames = dimnames(flowers.bv)))
# replace anything containing "NA±" with N.D.
flowers.mse <- as.data.frame(apply(flowers.mse, 
                                   2, # apply to columns not rows
                                   function(x) gsub(".*NA.*", "N.D.", x))) 

# how above works (for future reference)
  # apply() applies function(x) (in this case, gsub) to the whole data frame flowers.mse
    # we specify "2" to indicate columns (rows would be 1)
  # in gsub(), specifying ".*NA.*" would replace anything containing "NA" to "N.D."  
    # just using "NA" would replace the string "NA" with "N.D." and not the whole value.

flowers.mse <- rbind(retainedRIs, flowers.mse)
flowers.mse <- as.data.frame(t(flowers.mse)) #transposing

abbrvars <- c("A2", "CN", "GV", "GD", "LE",
             "LL", "LM", "MM", "DP", "NF")

names(flowers.mse) <- c("RI", abbrvars)
flowers.mse$Compound <- rownames(flowers.mse)
flowers.mse <- flowers.mse%>% 
  relocate(`Compound`, 
           .before = `RI`)

# exporting to .docx
# custom code to make autofit work
fit_flex_to_page <- function(ft, pgwidth = x){
  
  ft_out <- ft %>% autofit()
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}
# specify desired page properties
horiz_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 8.5, height = 11
  ),
  page_margins = page_mar()
)

# export
table2 <- as_flextable(flowers.mse, 
                             max_row = nrow(flowers.mse), 
                             show_coltype = F) %>% 
  fontsize(size =10 ) %>%
  fontsize(size = 10, part = "header") %>%
  colformat_double(digits = 0, big.mark = "") %>% 
  theme_vanilla() %>% 
  fit_flex_to_page(pgwidth = 9) %>%
  save_as_docx(path = "Table_2.docx", pr_section = horiz_properties) 

# 5.2 APPENDIX TABLE A3: LEAF CHEMISTRY SUMMARY
# get means for each variety
feuilles <- chems.n[c(38:47),]
feuilles <- cbind(chems$variety[c(38:47)], feuilles)
names(feuilles)[1] <- "variety"

leaves.bvo <- feuilles %>% 
  group_by(variety) %>%
  summarize_all(list(average = mean))
leaves.bvar <- leaves.bvo[,-1]
names(leaves.bvar) <- names(chems.n)
leaves.bv <- rbind(RIs, leaves.bvar )


# keep only values that are at least 1% of one variety
leaf.sums <- rowSums(leaves.bv) # calculate total chemicals per variety
leaf.logic <- as.data.frame((apply(leaves.bv, 
                                     2, function(x) x/leaf.sums)) > 0.01)

# replace logic with numerical values
leaf.logic[leaf.logic == FALSE] <- 0
leaf.logic[leaf.logic == TRUE] <- 1
leaf.logic[1,] <- rep(1, times = ncol(leaf.logic)) # so RIs are kept

# drop chems where TRUE < 1
leaves.bv <- leaves.bv[colSums(leaf.logic) > 1] 

# replace zeroes with NA for next steps
leaves.bv <- leaves.bv %>%
  select(where(~is.numeric(.x))) %>%
  round(digits = 0)
leaves.bv[leaves.bv < 20] <- NA

# replace anything containing "NA" with N.D.
leaves.bv <- as.data.frame(leaves.bv)
leaves.bv[is.na(leaves.bv) == TRUE] <- "N.D."

#final touches before export
leaf.means <- as.data.frame(t(leaves.bv)) #transposing
names(leaf.means) <- c("RI", abbrvars[-10])
leaf.means$Compound <- rownames(leaf.means)
leaf.means <- leaf.means%>% 
  relocate(`Compound`, 
           .before = `RI`)
rownames(leaf.means) <- NULL

#export: leaf means
tableA1 <- as_flextable(leaf.means, 
                       max_row = nrow(leaf.means), 
                       show_coltype = F) %>% 
  fontsize(size =9 ) %>%
  fontsize(size = 9, part = "header") %>%
  colformat_double(digits = 0, big.mark = "") %>% 
  theme_vanilla() %>% 
  fit_flex_to_page(pgwidth = 6.5) %>%
  save_as_docx(path = "Table_A3.docx", pr_section = horiz_properties) 

# 5.3 Creating Datafile for Regression Analyses
# merge flower and leaf means together
names(flowers.bvo) <- c("variety", paste(names(chems.n), 
                            rep("floral",
                                times = length(names(chems.n)))))
names(leaves.bvo) <- c("variety", names(chems.n))
regressions<- cbind(flowers.bvo[-10,-1], leaves.bvo[,-1])
  # note that there are no leaf samples for swt nufar
regressions$variety <- vars[-10]
regressions <- regressions %>% 
  relocate(`variety`, 
           .before = `Sabinene floral`)

# export
write_csv(regressions, "regressions.csv")

# 5.4 Figure 1 & A1: NMDS Ordinations
# concatenate plots for fig 1
plot.nmds <- ggarrange(plot.nmds.f, plot.nmds.fr, plot.nmds.a, 
                    labels = c("A", "B", "C"),
                    font.label = list(size = 20, color = "black", 
                                      face = "bold", 
                                      family = NULL),
                    ncol = 3, nrow = 1)

# export fig 1 - IGNORE ggrepel warnings
ggsave(plot = plot.nmds,
       filename = "Fig_1.pdf",
       width = 12.5*3,
       height = 12, 
       units = "cm")

# exporting fig a1 
ggsave(plot = plot.nmds.lvf,
       filename = "Fig_A1.pdf",
       width = 12.5*1,
       height = 12,
       units = "cm")

# 5.5 APPENDIX TABLE A4: LEAF VS. FLOWERS SIMPER TABLE
colnames(simper.output)
simper.lvf.ft <- as_flextable(simper.output, 
                            max_row = nrow(simper.output), 
                            show_coltype = F) %>% 
    colformat_double(digits = 3) %>% 
    theme_vanilla() %>% 
    bold(part = "body", # bold only significant 
         ~p < 0.05) %>%
    fit_flex_to_page(pgwidth = 9.25) %>%
    fontsize(size = 12) %>%
    fontsize(size = 12, part = "header") %>%
    save_as_docx(path = "Table_A4.docx", pr_section = horiz_properties) #save

# 5.6 APPENDIX TABLE A2: FLOWER VARIETY SIMPER ANALYSIS LETTERS
# Convert to flextable
names(f.simper.cld) <- c("Compound", "A2", "CN", "GV", "GD", "LE", "LL", "LM", "MM",
                         "DP", "NF")
simper.f.ft <- as_flextable(f.simper.cld, 
                           max_row = nrow(f.simper.cld), 
                           show_coltype = F) %>% 
  theme_vanilla() %>% 
  fit_flex_to_page(pgwidth = 9) %>%
  fontsize(size = 12) %>%
  fontsize(size = 12, part = "header") %>%
  save_as_docx(path = "Table_A2.docx", pr_section = horiz_properties) #save

# 5.7 EXTRA STATS: ANOVA & PERMANOVA WHILE REMOVING <1% 
nova.f <- flowers[,-c(1:5)]
nova.f <- nova.f[colSums(flower.logic) > 1]
dist.nova <- vegdist(nova.f, method = "bray")

# dispersion test
(disptest.nova <- anova(betadisper(dist.nova, flowers$variety)))
#ns - ANOSIM and PERMANOVA OK

# anosim
(anosimnova <-  anosim(dist.nova, flowers$variety, permutations = 10000))
# sig dif

# permanova
(basil.perm.n <- adonis2(dist.nova ~ flowers$variety, method = "bray"))
#sig dif

# 6: Unused - PCA Code and Simple Stats. WILL NOT RUN :)####
if(FALSE){
library(FactoMineR)
library(factoextra)   
  
# 6.1 Leaf v. Flower PCA
#redo log transformation
min.val.lvf <- min(abs(chems.n[chems.n!=0])) # finds minimum value in dataset, divides by 10
chems.log.lvf <- log10((chems.n + sqrt(chems.n^2 + min.val.lvf^2))/10) # arcsinh transforming


#creating the distance matrix
chems.pca.lvf <- prcomp(chems.log.lvf, center = TRUE, scale. = FALSE)
cpca.lvf <- summary(chems.pca.lvf)

#scree plot
fviz_eig(chems.pca.lvf, addlabels=TRUE)

#biplot - shows important variables
fviz_pca_var(chems.pca.lvf, repel = TRUE)

# importance of variables
cos21 <- fviz_cos2(chems.pca.lvf, choice = "var", axes =1:2)
cos21 <- cos21 + theme(text = element_text(size = 20)) 
cos21


# 6.2: Plotting LvF
#creating groups
chems %>% count(type) # use these values below
gp <- c(rep("Flower", times = 37), rep("Leaf", times = 10))

#plotting: method 1
plot1 <- fviz_pca_ind(chems.pca.lvf, repel = TRUE, pointsize = 7, habillage = chems$type, 
                      axes.linetype = NA, geom = "point", col.ind = gp,
                      palette = c("#93B7BE", "#554348"), addEllipses = TRUE,
                      ellipse.level = 0.95, invisible = "quali") +
  theme_minimal() +
  xlab("PC 1 (32.9%)") +
  ylab("PC 2 (11.5%)")

plot1 + theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              text = element_text(size = 20))

#extracting PCs for later analyses
PCs <- as.data.frame(chems.pca.lvf$x)
PCs$sample <- paste(chems$type, chems$variety, sep = "_")
PCs$sample[43] <- "leaf_grk dwarf"

#alternative plot
pc.plot <- ggplot(PCs, aes(x = PC1, y = PC2, color = Sample_Type))
pc.plot <- pc.plot + geom_point()
pc.plot <- pc.plot + stat_ellipse()

# getting mean values for each sample type
pc.bt <- PCs %>% 
  group_by(sample)

pc.means <- pc.bt %>%
  summarize_all(list(average = mean))

# 6.3 Plotting Vars 1
#log transforming data
min.val.f <- min(abs(flowers.bvar[flowers.bvar!=0])) # finds minimum value in dataset, divides by 10
chems.log.f <- log10((flowers.bvar + sqrt(flowers.bvar^2 + min.val.f^2))/10) # log transforming

#creating the PCA
chems.pca.f <- prcomp(chems.log.f, center = TRUE, scale. = FALSE)
cpca.f <- summary(chems.pca.f)

#scree plot
fviz_eig(chems.pca.f, addlabels=TRUE)

#biplot - shows important variables
fviz_pca_var(chems.pca.f, repel = TRUE)

# importance of variables
cos22<- fviz_cos2(chems.pca.f, choice = "var", axes =1:2)
cos22 <- cos22 + theme(text = element_text(size = 20))
pca.info <- cos22$data %>% 
  arrange(-cos2)
names(pca.info) <- c("Chemical", "importance")

# 6.4 Plotting Vars 2
#creating groups
chems.fl <- chems[c(1:37),]
chems.fl %>% count(variety) # use these values below
gp.f <- c(rep("Aroma 2", times = 3), rep("Cinnamon", times = 3),
          rep("Genovese", times = 4), rep("Greek Dwarf", times = 4),
          rep("Lemon", times = 3), rep("Lettuce", times = 5),
          rep("Lime", times = 3), rep("Mammolo", times = 5),
          rep("Dark Purple Opal", times = 4), rep("Nufar", times = 3))

#plotting ellipses: ugly, but will be useful for future plots
bb <-fviz_pca_ind(chems.pca.f, repel = TRUE, pointsize = 3, 
                  labelsize = 5, col.ind = gp.f,
                  habillage = chems.fl$variety,
                  addEllipses = FALSE, ellipse.level = 0.95,
                  palette = "Spectral")

# extracting x-y coordinates for the centroids
pc.x<- bb[["data"]][["x"]]
pc.y<- bb[["data"]][["y"]]
pc.df <- data.frame(pc.x, pc.y)
pc.df$variety <- chems.fl$variety

# renaming shortcutted varieties
pc.df <- pc.df %>%
  mutate(variety =
           replace(variety, variety == "grk dwarf", "greek dwarf")) 

pc.df <- pc.df %>%
  mutate(variety =
           replace(variety, variety == "purp opal", "purple opal")) 

pc.df <- pc.df %>%
  mutate(variety =
           replace(variety, variety == "swt nufar", "sweet nufar")) 

# grouping x-y coords by variety
pc.df_by <- pc.df %>%
  group_by(variety) 

# rearranges x-y coords so polygons are formed around points
pc.df_pol <- pc.df_by%>%
  do(.[chull(.[1:2]), ]) 

# calculating x-y means for centroids
x.y.list <- pc.df_by %>%
  summarize_all(list(average = mean))

# setting up ggplot, plotting samples & polygons
plot2 <- ggplot(pc.df_pol, aes(x = pc.x, y = pc.y,
                               col = factor(variety),
                               fill = factor(variety)))

plot2 <- plot2 + geom_polygon(alpha = 0.3
                              #, fill = NA
) +
  geom_point(data = pc.df_pol, aes(x = pc.x, y = pc.y, col = factor(variety))) +
  xlim(-10,16) +
  ylim(-7.5,15)

# plotting the centroids
plot2 <- plot2 + geom_point(data = x.y.list, aes(x = pc.x_average, pc.y_average,
                                                 color = factor(variety)), size = 5, shape = 18)

#recoloring 
coolors <- c("#781c81", "#43328d", "#99bd5c",  "#d92120",
             "#70b484", "#416fb8", "#c3ba45", "#e0a239",
             "#e66b2d", "#519cb8")
#guide:
#"lime" = "#781c81", "cinnamon" = "#43328d", "genovese" = "#99bd5c",
#"lettuce" = "#519cb8", "lemon" = "#70b484", "aroma 2" = "#416fb8",
#"mammolo" = "#c3ba45", "grk dwarf" = "#e0a239", "purp opal" = "#e66b2d",
#"swt nufar" = "#d92120"

plot2 <- plot2 + scale_color_manual(values = coolors) +
  scale_fill_manual(values = coolors)

plot2 <- plot2 + geom_text_repel(data = x.y.list, 
                                 mapping = aes(x = pc.x_average,
                                               pc.y_average,
                                               label = variety),
                                 fontface = "bold",
                                 nudge_y = c(2,1,-2,-5,1,-4,-2,-3,2,4),
                                 nudge_x = c(5,4,-2,-1,4,10,-2,3,5,1),
                                 color = coolors,
                                 min.segment.length = 0.1,
                                 size = 8)
#box.padding = unit(2, "lines"))

plot2 <- plot2 + labs(x = "PC 1 (24.8%)", 
                      y = "PC 2 (15.2%)")
plot2 <- plot2 + theme(axis.title = element_text(size = 20))
plot2 <- plot2 + theme(axis.text = element_text(size = 20))


plot2 <- plot2 + theme(legend.position = "none")
plot2 <- plot2 + theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank())
plot2

# some simple stats
chem.means$total <- rowSums(chem.means[,-1])
feuilles$total <- rowSums(feuilles[,-1])
AOVdata <- data.frame(chem.means[,1], chem.means[,73])
AOVdata[c(11:20),] <- data.frame(feuilles[,1], feuilles[,73])
AOVdata$source <- c(rep("flower", times = 10), rep("leaf", times = 10))

summary(aov(total ~ variety*source, data = AOVdata))

#old data tables
#mins and maxs
chem.min.max <- flowers.bv %>%
  summarize_all(list(min, max))

#mean
chem.means <- flowers.bv %>%
  summarize_all(list(average = mean))

# replace zeroes with NA for next steps
flowers.all <- flowers.full %>%
  select(where(~is.numeric(.x))) %>%
  round(digits = 0)
flowers.all[flowers.all < 20] <- NA

# create a logical df where TRUE means a chemical ≠ NA
flower.logic1 <- as.data.frame(is.na(flowers.all))

#replacing logic with numerical values
flower.logic1[flower.logic1 == TRUE] <- 0
flower.logic1[flower.logic1 == FALSE] <- 1


}

# Multivariate analyses for basil VOC data 
  # created by/on: S Clemente / 7 October 2023
  # last updated by/on: S Clemente / 7 January 2025

# 1: Setup #### 
if(TRUE){
# opening packages
library(tidyverse) # for data organization, manipulation, & visualization
library(corrr) # for correlation plots
library(ggcorrplot) # also for correlation plots
library(vegan) # for multivariate analyses and visualizations
library(ggrepel) # for neat labels 
library(flextable) # for exporting into .docx tables
library(officer) # for declaring .docx format
library(rcompanion) # for getting compact letter display in SIMPER table
library(ggpubr) # for concatenating NMDS plots
library(readr) # for exporting into .csv files


# set wd, upload data
chems <- read.csv('basil_chemicals_all.csv',
                  check.names = FALSE)
#re-sorting by group
chems <- chems[order(chems$type),]

#subsetting numerical data only
chems.n <- chems[,6:90]

#removing RIs from main chems dataset
chems <- chems[-1,]

#removing unresolved compounds / contaminants
colnames(chems.n) # check!
chems.n <- chems.n[,-85] # removing 2,6,10,14-tetramethyl...octadecane (unresolved)
chems.n <- chems.n[,-82] # removing 2,6,10,14-tetramethyl...pentadecane (unresolved)
chems.n <- chems.n[,-68] # removing benzenecycloheptane (contaminant)
chems.n <- chems.n[,-29] # removing eicosane (contaminant)
chems.n <- chems.n[,-26] # removing cuminal (contaminant)
chems.n <- chems.n[,-20] # removing methylindene (unresolved)
chems.n <- chems.n[,-19] # removing borneol (unresolved)
chems.n <- chems.n[,-12] # removing trans-oxiran-3-undecylmethyl ester (contaminant)
chems.n <- chems.n[,-10] # removing 2-hexen-1-yl dodecanoic acid ester (all 0)
chems.n <- chems.n[,-1] # removing acetic acid butyl ester (contaminant)
colnames(chems.n) # check again~

# create vector for cultivar names 
vars <- c("Aroma 2", "Cinnamon", "Genovese", "Greek Dwarf", "Lemon",
          "Lettuce", "Lime", "Mammolo", "Dark Purple Opal", "Sweet Nufar")

#creating a matrix to check for correlations
#log transforming data
chems.p <- chems.n[-1,]
min.val.check <- min(abs(chems.p[chems.p!=0])) # finds minimum value in dataset, divides by 10
chems.log.check <- log10((chems.p + sqrt(chems.p^2 + min.val.check^2))/10) # log transforming

corr.matrix.check <- cor(chems.log.check)
ggcorrplot(corr.matrix.check,
           tl.cex = 6)

# creating a smaller matrix to see later relationships
chems.z <- chems.log.check[,48:75]
corr.matrix.z <- cor(chems.z)
ggcorrplot(corr.matrix.z,
           tl.cex = 6)

if(TRUE){
#summing up very correlated peaks
#(use colnames(chems.n) to get column location)
#correlated peaks #1: epicubenol + tau_cadinol
chems.n$epic <- chems.n$Epicubenol + chems.n$`tau-Cadinol`
chems.n <- chems.n %>% relocate(epic, .before = Epicubenol) #reorder
chems.n <- chems.n[,-c(67,69)] # remove
colnames(chems.n)[66] <- "Epicubenol"
colnames(chems.n) ##check

#correlated peaks #2: Unknown 71/43 and Unknown 108/81
chems.n$Unknown_109_81_43 <- chems.n$`Unknown 109, 81, 43` + chems.n$`Unknown 71, 43, 79`
chems.n <- chems.n %>% relocate(`Unknown_109_81_43`, .before = `Unknown 71, 43, 79`) #reorder
chems.n <- chems.n[,-c(63,65)] # remove
colnames(chems.n)[62] <- "Unknown 109, 81, 43"

#correlated peaks #3: b_caryophyllene
chems.n$`beta-Caryophyllene` <- chems.n$`β-Caryophyllene` + chems.n$`β-Ylangene`
chems.n <- chems.n %>% relocate(`beta-Caryophyllene`, .before = `β-Ylangene`) #reorder
chems.n <- chems.n[,-c(38:39)] # remove
colnames(chems.n)[37] <- "β-Caryophyllene"

#correlated peaks #4: methyl eugenol
chems.n$`Methyl eugenol` <- chems.n$Methyleugenol + chems.n$`β-Elemene`
chems.n <- chems.n %>% relocate(`Methyl eugenol`, .before = Methyleugenol) #reorder new column
chems.n <- chems.n[,-c(31:32)] #remove old columns

#correlated peaks #5: beta ocimene
chems.n$`B-Ocimene` <- chems.n$`β-Ocimene` + chems.n$Benzeneacetaldehyde #sum
chems.n <- chems.n %>% relocate(`B-Ocimene`, .before = `β-Ocimene`) #reorder
chems.n <- chems.n[,-c(8:9)] # remove
colnames(chems.n)[7] <- "β-Ocimene"

#check that all is fine
colnames(chems.n) 

#changing troublesome names
names(chems.n)
names(chems.n)[23] <- "Unknown 97, 111, 70, 71"
names(chems.n)[32] <- "Unknown 95, 161, 107, 93"
names(chems.n)[48] <- "Unknown 81, 105, 119, 161"
names(chems.n)[59] <- "Unknown 109, 81, 93, 95"
names(chems.n)[64] <- "Unknown 161, 204, 105, 162"
names(chems.n)[65] <- "Unknown 105, 204, 71, 81"
names(chems.n)[67] <- "Unknown 93, 107, 81, 91"
names(chems.n)[68] <- "Unknown 95, 97, 81, 110"
names(chems.n)[69] <- "Unknown 161, 204, 91, 105"
names(chems.n)[70] <- "Unknown 161, 204, 95, 81"

#reassign compound RIs
chems.n[1,62] <- 1640
chems.n[1,59] <- 1618
chems.n[1,35] <- 1444
chems.n[1,29] <- 1406
chems.n[1,7]  <- 1048
}

#extract RIs
RIs <- (chems.n[1,])
chems.n <- chems.n[-1,]
}

# 2: Flowers v. Leaves ####
# 2.1 difference tests 
# create distance matrix
dist.lvf <- vegdist(chems.n, method = "bray")

#dispersion test
(disptest.lvf <- anova(betadisper(dist.lvf, chems$type)))
#statistically significant. PERMANOVA preferred over ANOSIM, indicate in paper
  # see: http://doi.org/10.1890/12-2010.1

# assessing statistical differences between groups
(basil.perm <- adonis2(dist.lvf ~ chems$type,
                      method = "bray"))
 # significantly different!

# 2.2 SIMPER 
simper.lvf <- simper(comm = chems.n, 
                   group = chems$type,
                   permutations = 10000)
simper.output <- summary(simper.lvf)
simper.output <- simper.output$flower_leaf

# Row names as a column
simper.output$`Compound Name` <- rownames(simper.output)
rownames(simper.output) <- NULL

# Make trait/compound names first column
simper.output <- simper.output %>% 
  relocate(`Compound Name`, 
           .before = average)

# renaming
colnames(simper.output) <- c("Compound Name", "Average", "SD", "Ratio",
                             "Avg. Flower", "Avg. Leaf", "Cum. Sum", "p")


# 2.3 NMDS
# first standardize data 
dissim <- chems.n %>%
  decostand(method = "hellinger")

# set seed for reproduciblity. 
  # different seeds will simply change the orientation of polygons.
  # but polygon shape and distances will remain the same.
set.seed(2025)
nmds.lvf <- metaMDS(dissim, 
                    autotransform = FALSE,
                    distance = "bray",
                    k = 7,
                    model = "global",
                    try = 30,
                    trymax = 200)

# k = 7 allows for n > 2k+1 for both leaf & flower data
  # while minimizing stress. Although stress is >0.1 from k=3 onwards

# Get scores
nmds.lvf.scores <- as_tibble(scores(nmds.lvf)$sites)

# Add sample type
nmds.lvf.scores <- nmds.lvf.scores %>% 
  mutate(type = chems$type)

# Identify extreme points of convex hull
# convex hull values - flowers
hull.flowers<- 
  nmds.lvf.scores[nmds.lvf.scores$type == "flower", ][chull(nmds.lvf.scores[
    nmds.lvf.scores$type == "flower", c("NMDS1", "NMDS2")]), ] 

# convex hull values - leaves
hull.leaves <- 
  nmds.lvf.scores[nmds.lvf.scores$type == "leaf", ][chull(nmds.lvf.scores[
    nmds.lvf.scores$type == "leaf", c("NMDS1", "NMDS2")]), ]  

# Combine
hull.scores <- rbind(hull.flowers, 
                       hull.leaves)  

# get centroids 
centroids1 <- nmds.lvf.scores %>% 
                group_by(type) %>%
                  summarize(NMDS1a = mean(NMDS1),
                            NMDS2a = mean(NMDS2)) 
centroids1$NMDS1a <- centroids1$NMDS1a + c(0.5, 0)
centroids1$NMDS2a <- centroids1$NMDS2a + c(-0.41, 0)

names(centroids1) <- c("type","NMDS1", "NMDS2")

# Plot - IGNORE ggrepel warnings
plot.nmds.lvf <- ggplot(nmds.lvf.scores) + 
  geom_point(aes(x = NMDS1, # points
                 y = NMDS2, 
                 color = type, 
                 shape = type), 
             size = 3) + 
  geom_polygon(data = hull.scores, # shaded area
               aes(x = NMDS1,
                   y = NMDS2, 
                   fill = type, 
                   col = type),
               alpha = 0.3) + 
  scale_color_manual(name = "Type",
                     values = c("#93B7BE", "#554348"), 
                     labels = c("flower", "leaf")) +
  scale_fill_manual(name = "Type",
                    values = c("#93B7BE", "#554348"),
                    labels = c("flower", "leaf")) +
  scale_shape_manual(name = "Type",
                     values = c(18,19),
                     labels = c("flower", "leaf")) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") 

# adding labels
plot.nmds.lvf <-
  plot.nmds.lvf + geom_text_repel(data = centroids1, 
                                mapping = aes(x = NMDS1,
                                              y = NMDS2,
                                              label = c("flower", 
                                                        "leaf")),
                                fontface = c("bold"),
                                color = c("#93B7BE", "#554348"),
                                min.segment.length = 10,
                                direction = "both",
                                size = 8,
                                box.padding = 0,
                                point.padding = 0,
                                max.overlaps = 10) +
  annotate("rect",
           xmin = -0.1, xmax = 1,
           ymin = -1, ymax = -1.23,
           alpha = 1, fill = "#000000") +
  annotate(
    "text", label = paste("Stress = ", 
                          round(nmds.lvf$stress, 
                                digits = 3)),
    x = 0.45, y = -1.121, size = 8, fontface = "bold",
    color = "#ffffff")

plot.nmds.lvf
  
# 3 Varieties ####
# 3.1: Setting Up
# creating a flower-only dataset
flowers <- chems.n[c(1:37),]

#re-adding important information
chems.info <- chems[c(1:37),c(1:5)]
flowers <- cbind(chems.info, flowers)

# rearranging by variety
flowers <- flowers[order(flowers$variety),]

# create distance matrix
dist.f <- flowers[,-5] %>% 
  select(where(~is.numeric(.x))) %>% 
    vegdist(method = "bray")

# 3.2 Difference tests
# dispersion test
(disptest.f <- anova(betadisper(dist.f, flowers$variety)))
#ns - ANOSIM and PERMANOVA OK

# anosim
(anosim19 <-  anosim(dist.f, flowers$variety, permutations = 10000))
# sig dif

# permanova
(basil.perm.f <- adonis2(dist.f ~ flowers$variety, method = "bray"))
#sig dif

# 3.3 SIMPER 
# WARNING - computationally taxing. leave as TRUE if not in a rush.
if(TRUE){
simper.f <- simper(comm = flowers[-5] %>% 
                     select(where(~is.numeric(.x))), 
                     group = flowers$variety,
                     permutations = 10000)
simper.f.output <- summary(simper.f)

# extracting just p-values from large simper object
#taken from Konrad Rudolph at
  #https://stackoverflow.com/questions/23758858/how-can-i-extract-elements-from-lists-of-lists-in-r
p.vals <- lapply(simper.f.output, `[[`, "p") #added an extra "["
p.table <- as.data.frame(do.call(rbind, p.vals))
p.table <- p.table %>% 
            mutate_all(~replace(., is.na(.), 0.9999)) #replacing NA with 1
            names(p.table) <- names(chems.n)

# Create "ID" indicating pairwise comparisons
p.table$Comparison <- gsub("_", "-", names(p.vals))

# Make trait/compound names first column
p.table <- p.table %>% 
  relocate(`Comparison`, 
           .before = `Sabinene`)

# now get tukey letters for every chemical 
p.list <- as.list(p.table[2:71])
t.lets <- lapply(p.list, function(x)
                cldList(p.value = x, 
                        comparison = p.table$Comparison))
t.list <- lapply(t.lets, `[[`, "Letter")
f.simper.cld <- as.data.frame(do.call(rbind, t.list))
names(f.simper.cld ) <- vars
f.simper.cld $Compound <- rownames(f.simper.cld )

rownames(f.simper.cld) <- NULL

f.simper.cld <- f.simper.cld %>% 
  relocate(`Compound`, 
           .before = `Aroma 2`)
}

# 3.4 NMDS 
# first standardize data 
dissim.f <- flowers[-5] %>%
  select(where(~is.numeric(.))) %>%
    decostand(method = "hellinger")

# run nmds
set.seed(2025) # so results are reproducible
nmds.f <- metaMDS(dissim.f, 
                    autotransform = FALSE,
                    distance = "bray",
                    k = 3,
                    model = "global",
                    try = 30,
                    trymax = 200)

# Get scores
nmds.f.scores <- as_tibble(scores(nmds.f)$sites)

# Add variety
nmds.f.scores <- nmds.f.scores %>% 
  mutate(variety = flowers$variety)

# Identify extreme points of convex hull for each variety
hull.a2 <- 
  nmds.f.scores[nmds.f.scores$variety == "aroma 2", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "aroma 2", c("NMDS1", "NMDS2")]), ] 

hull.cn <- 
  nmds.f.scores[nmds.f.scores$variety == "cinnamon", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "cinnamon", c("NMDS1", "NMDS2")]), ] 

hull.gv <- 
  nmds.f.scores[nmds.f.scores$variety == "genovese", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "genovese", c("NMDS1", "NMDS2")]), ] 

hull.gd <- 
  nmds.f.scores[nmds.f.scores$variety == "grk dwarf", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "grk dwarf", c("NMDS1", "NMDS2")]), ] 

hull.lm <- 
  nmds.f.scores[nmds.f.scores$variety == "lemon", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "lemon", c("NMDS1", "NMDS2")]), ] 

hull.lt <- 
  nmds.f.scores[nmds.f.scores$variety == "lettuce", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "lettuce", c("NMDS1", "NMDS2")]), ] 

hull.li <- 
  nmds.f.scores[nmds.f.scores$variety == "lime", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "lime", c("NMDS1", "NMDS2")]), ] 

hull.mm <- 
  nmds.f.scores[nmds.f.scores$variety == "mammolo", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "mammolo", c("NMDS1", "NMDS2")]), ] 

hull.po <- 
  nmds.f.scores[nmds.f.scores$variety == "purp opal", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "purp opal", c("NMDS1", "NMDS2")]), ] 

hull.sn <- 
  nmds.f.scores[nmds.f.scores$variety == "swt nufar", ][chull(nmds.f.scores[
    nmds.f.scores$variety == "swt nufar", c("NMDS1", "NMDS2")]), ] 


# Combine
hullf.scores <- rbind(hull.a2, hull.cn, hull.gv, hull.gd, hull.lm, 
                     hull.lt, hull.li, hull.mm, hull.po, hull.sn)

# Plot
#necessary values for plot
cols <- c("#781c81", "#43328d", "#99bd5c", "#d92120","#70b484", 
          "#416fb8", "#c3ba45", "#e0a239", "#e66b2d","#519cb8")
null.cols <- c("#ffffff", "#43328d", "#ffffff", "#ffffff","#70b484", 
          "#ffffff", "#c3ba45", "#ffffff", "#ffffff","#000000")
null.vars <- c("All others", "Cinnamon", "All others", "All others", "Lemon",
           "All others", "Lime", "All others", "All others", "All others")
centroids <- nmds.f.scores %>%
  group_by(variety) %>%
  summarize(NMDS1 = mean(NMDS1),
            NMDS2 = mean(NMDS2))

# plotting - IGNORE ggrepel warnings
if(TRUE){
plot.nmds.f <- ggplot(nmds.f.scores) + 
  geom_point(aes(x = NMDS1, # points
                 y = NMDS2, 
                 color = variety), 
             size = 3) + 
  geom_polygon(data = hullf.scores, # shaded area
               aes(x = NMDS1,
                   y = NMDS2, 
                   fill = variety, 
                   col = variety),
               alpha = 0.3) + 
  scale_color_manual(name = "Variety",
                     values = cols, 
                     labels = vars) +
  scale_fill_manual(name = "Variety",
                    values = cols,
                    labels = vars) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") 

# adding labels
plot.nmds.f <-
  plot.nmds.f + geom_text_repel(data = centroids, 
                                mapping = aes(x = NMDS1,
                                              y = NMDS2,
                                             label = null.vars),
                                xlim = c(-Inf, NA),
                                fontface = c("italic", "bold", "italic","italic","bold",
                                             "italic","bold","italic","italic","italic"),
                                nudge_x = c(0,0.6,0,0,0.1,0,0,0,0,0.6),
                                nudge_y = c(0,0.4,0,0,0,0,0.5,0,0,0.1),
                                color = null.cols,
                                min.segment.length = 0,
                                direction = "both",
                                size = 8,
                                box.padding = 2.5,
                                point.padding = 0.5,
                                max.overlaps = 10) +
              annotate("rect",
                       xmin = 0.02, xmax = 1.07,
                       ymin = -0.6, ymax = -0.8,
                       alpha = 1, fill = "#000000") +
              annotate("text", 
                         label = paste("Stress = ",round(nmds.f$stress, 
                              digits = 3)),
                         x = 0.55, 
                         y = -0.7, 
                         size = 8, 
                         fontface = "bold",
                         color = "#ffffff") 

plot.nmds.f
}

# 4: Creating Tables, Exporting Tables & Figs ####
# 4.1 TABLES 2 & S2: ALL FLORAL CHEMISTRY BY CULTIVAR
  # FOR TABLE 2, RUN AS IS.
  # TO MAKE TABLE S2, DON'T RUN LINES 508, 539-540
    # and be sure to rename .docx filename at L612

# get means for each variety
fleurs <- flowers[,-c(1:2,4:5)]
flowers.bvo <- fleurs %>% 
  group_by(variety) %>%
  summarize_all(list(average = mean))
flowers.bvar <- flowers.bvo[,-1]
names(flowers.bvar) <- names(chems.n)
flowers.bv <- rbind(RIs, flowers.bvar )

# keep only values that are at least 1% of one variety
flower.sums <- rowSums(flowers.bv) # calculate total chemicals per variety
# create a logical df where TRUE means a chemical is at least 1%
flower.logic <- as.data.frame((apply(flowers.bv, 
                     2, function(x) x/flower.sums)) > 0.01)
# replace logic with numerical values
flower.logic[flower.logic == FALSE] <- 0
flower.logic[flower.logic == TRUE] <- 1
flower.logic[1,] <- rep(1, times = ncol(flower.logic)) # so RIs are always 1

# index original dataset so that only columns with at least 
  # one "1" in flower.logic (except the RIs, which are always 1) are kept
  # TO CREATE TABLE S2, DON'T RUN L508, proceed with rest. 
flowers.bv<- flowers.bv[colSums(flower.logic) > 1] 

#extract RIs of retained compounds, remove from dataset for now.
retainedRIs <- flowers.bv[1,]
flowers.bv <- flowers.bv[-1,]

# replace zeroes with NA for next steps
flowers.bv <- flowers.bv %>%
  select(where(~is.numeric(.x))) %>%
  round(digits = 0)
flowers.bv[flowers.bv < 20] <- NA

# get standard errors for each variety
flowers.sd <- fleurs %>%
  group_by(variety) %>%
   summarize_all(list(se = sd))
    
samples <- fleurs %>%
  group_by(variety) %>%
  count()

flowers.se <- flowers.sd %>%
  select(where(~is.numeric(.x))) %>%
    mutate(./samples$n) %>%
     round(digits = 0)
names(flowers.se) <- names(chems.n)
flowers.se.full <- flowers.se

# keep only se values for remaining chemicals in flowers.bv
# DON'T RUN 538-539 IF YOU WANT TABLE S2
flowers.se<- flowers.se[colSums(flower.logic[-1,]) > 0] 
names(flowers.bv) == names(flowers.se) #check

#convert numbers to characters
flowers.bv <- flowers.bv %>%
  select(where(~is.numeric(.x))) %>%
  mutate(across(everything(), as.character))

flowers.se <- flowers.se %>%
  select(where(~is.numeric(.x))) %>%
  mutate(across(everything(), as.character))

# merge dataframes together
# first convert mean and sd to matrices
flowers.bv <- as.matrix(flowers.bv)
flowers.se <- as.matrix(flowers.se)

# then merge 
  # (thanks to Aaron at:https://stackoverflow.com/questions/6408016/paste-together-two-dataframes-element-by-element)

flowers.mse <- as.data.frame(matrix(paste(flowers.bv, flowers.se, sep = " ± "),
                                    nrow = nrow(flowers.bv), 
                                    dimnames = dimnames(flowers.bv)))
# replace anything containing "NA±" with N.D.
flowers.mse <- as.data.frame(apply(flowers.mse, 
                                   2, # apply to columns not rows
                                   function(x) gsub(".*NA.*", "N.D.", x))) 

# how above works (for future reference)
  # apply() applies function(x) (in this case, gsub) to the whole data frame flowers.mse
    # we specify "2" to indicate columns (rows would be 1)
  # in gsub(), specifying ".*NA.*" would replace anything containing "NA" to "N.D."  
    # just using "NA" would replace the string "NA" with "N.D." and not the whole value.

flowers.mse <- rbind(retainedRIs, flowers.mse)
flowers.mse <- as.data.frame(t(flowers.mse)) #transposing

abbrvars <- c("A2", "CN", "GV", "GD", "LE",
             "LL", "LM", "MM", "DP", "NF")

names(flowers.mse) <- c("RI", abbrvars)
flowers.mse$Compound <- rownames(flowers.mse)
flowers.mse <- flowers.mse%>% 
  relocate(`Compound`, 
           .before = `RI`)

# exporting to .docx
# custom code to make autofit work
fit_flex_to_page <- function(ft, pgwidth = x){
  
  ft_out <- ft %>% autofit()
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}
# specify desired page properties
horiz_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 8.5, height = 11
  ),
  page_margins = page_mar()
)
# export
table2 <- as_flextable(flowers.mse, 
                             max_row = nrow(flowers.mse), 
                             show_coltype = F) %>% 
  fontsize(size =10 ) %>%
  fontsize(size = 10, part = "header") %>%
  colformat_double(digits = 0, big.mark = "") %>% 
  theme_vanilla() %>% 
  fit_flex_to_page(pgwidth = 9) %>%
  save_as_docx(path = "Table_2.docx", pr_section = horiz_properties) 

# 4.2 TABLE 3: LEAF CHEMISTRY SUMMARY
# get means for each variety
feuilles <- chems.n[c(38:47),]
feuilles <- cbind(chems$variety[c(38:47)], feuilles)
names(feuilles)[1] <- "variety"

leaves.bvo <- feuilles %>% 
  group_by(variety) %>%
  summarize_all(list(average = mean))
leaves.bvar <- leaves.bvo[,-1]
names(leaves.bvar) <- names(chems.n)
leaves.bv <- rbind(RIs, leaves.bvar )


# keep only values that are at least 1% of one variety
leaf.sums <- rowSums(leaves.bv) # calculate total chemicals per variety
leaf.logic <- as.data.frame((apply(leaves.bv, 
                                     2, function(x) x/leaf.sums)) > 0.01)

# replace logic with numerical values
leaf.logic[leaf.logic == FALSE] <- 0
leaf.logic[leaf.logic == TRUE] <- 1
leaf.logic[1,] <- rep(1, times = ncol(leaf.logic)) # so RIs are kept

# drop chems where TRUE < 1
leaves.bv <- leaves.bv[colSums(leaf.logic) > 1] 

# replace zeroes with NA for next steps
leaves.bv <- leaves.bv %>%
  select(where(~is.numeric(.x))) %>%
  round(digits = 0)
leaves.bv[leaves.bv < 20] <- NA

# replace anything containing "NA" with N.D.
leaves.bv <- as.data.frame(leaves.bv)
leaves.bv[is.na(leaves.bv) == TRUE] <- "N.D."

#final touches before export
leaf.means <- as.data.frame(t(leaves.bv)) #transposing
names(leaf.means) <- c("RI", abbrvars[-10])
leaf.means$Compound <- rownames(leaf.means)
leaf.means <- leaf.means%>% 
  relocate(`Compound`, 
           .before = `RI`)
rownames(leaf.means) <- NULL

#export
table3 <- as_flextable(leaf.means, 
                       max_row = nrow(leaf.means), 
                       show_coltype = F) %>% 
  fontsize(size =9 ) %>%
  fontsize(size = 9, part = "header") %>%
  colformat_double(digits = 0, big.mark = "") %>% 
  theme_vanilla() %>% 
  fit_flex_to_page(pgwidth = 6.5) %>%
  save_as_docx(path = "Table_3.docx", pr_section = horiz_properties) 

# 4.3 Creating Datafile for Regression Analyses
# merge flower and leaf means together
names(flowers.bvo) <- c("variety", paste(names(chems.n), 
                            rep("floral",
                                times = length(names(chems.n)))))
names(leaves.bvo) <- c("variety", names(chems.n))
regressions<- cbind(flowers.bvo[-10,-1], leaves.bvo[,-1])
  # note that there are no leaf samples for swt nufar
regressions$variety <- vars[-10]
regressions <- regressions %>% 
  relocate(`variety`, 
           .before = `Sabinene floral`)

# export
write_csv(regressions, "regressions.csv")

# 4.4 Figure 1: NMDS Ordination
# concatenate plots 
plot.nmds <- ggarrange(plot.nmds.f, plot.nmds.lvf,
                    labels = c("A", "B"),
                    font.label = list(size = 20, color = "black", 
                                      face = "bold", 
                                      family = NULL),
                    ncol = 2, nrow = 1)

# export - IGNORE ggrepel warnings
ggsave(plot = plot.nmds,
       filename = "Fig_1.pdf",
       width = 25,
       height = 12, 
       units = "cm")

# 4.5 Supplemental S1: LEAF VS. FLOWERS SIMPER TABLE
colnames(simper.output)
simper.lvf.ft <- as_flextable(simper.output, 
                            max_row = nrow(simper.output), 
                            show_coltype = F) %>% 
    colformat_double(digits = 3) %>% 
    theme_vanilla() %>% 
    bold(part = "body", # bold only significant 
         ~p < 0.05) %>%
    fit_flex_to_page(pgwidth = 9.25) %>%
    fontsize(size = 12) %>%
    fontsize(size = 12, part = "header") %>%
    save_as_docx(path = "Supplemental_S1.docx", pr_section = horiz_properties) #save

# 4.6 Supplemental S3: FLOWER VARIETY SIMPER ANALYSIS LETTERS
# Convert to flextable
simper.f.ft <- as_flextable(f.simper.cld, 
                           max_row = nrow(f.simper.cld), 
                           show_coltype = F) %>% 
  theme_vanilla() %>% 
  fit_flex_to_page(pgwidth = 9) %>%
  fontsize(size = 12) %>%
  fontsize(size = 12, part = "header") %>%
  save_as_docx(path = "Supplemental_S3.docx", pr_section = horiz_properties) #save

# fin.
