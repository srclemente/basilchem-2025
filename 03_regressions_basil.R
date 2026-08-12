# Regressions for basil VOC data - all varieties
# originally created by SC / 4 April 2024
# last updated by SC / 12 August 2026

# 1 - SETUP ####
# opening packages, uploading data
library(ggplot2)
library(ggrepel)
library(ggcorrplot)
library(ggpubr)

# set wd, upload data
chems <- read.csv('regressions.csv', check.names = FALSE)

# changing names
chems[,1] <- c("A2", "CN", "GD", "GV", "LE", "LL", "LM", "MM", "DP")

#custom function for displaying R^2 values of regression lines
eq <- function(x,y) {
  m <- lm(y ~ x)
  as.numeric(format(summary(m)$r.squared, digits = 3))
}

# extract floral values 
flowers <- chems[,2:69]

# extract leaf values (even numbered columns)
leaves <- chems[,70:137]

#renaming
names(flowers) <- names(leaves)

# create correlation plot
corr.data <- ggcorrplot(cor(leaves,flowers))$data

# extract only values where same chemical is compared
  # between leaves and flowers
lvf <-corr.data[which(corr.data$Var1 == corr.data$Var2),]
  # note that not all 72 compounds are given
  # those where the compound is absent from all leaf samples
  # cannot be correlated (n = 11 compounds)

#order descending by absolute value of correlation
lvf<- lvf[order(-abs(lvf$value)), ]

# 2 - Getting the top 30 correlated compounds ####
lvf$Var1[1:50] 
  # note that many compounds are perfectly correlated
  # or close to so
  # because only one variety contains that compound

# for example, estragole:
  # which is only abundant in cv. Lime: 
if(TRUE){
summary(lm(`Estragole floral` ~ `Estragole`, data = chems)) #significant

estragole <- ggplot(chems, aes(x=`Estragole`, y = `Estragole floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE) +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$Estragole, 
                             chems$`Estragole floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                  mapping = aes(x = Estragole,
                                `Estragole floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 10,
                  size = 4,
                  max.overlaps = 20)
estragole <- estragole + labs(x = "Estragole", 
                              y = "")
estragole <- estragole  + theme(axis.title = element_text(size = 12, face = "bold"))
estragole <- estragole  + theme(axis.text = element_text(size = 10))


estragole <- estragole  + theme(legend.position = "none")
estragole <- estragole   + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank())

estragole }
 
# 3 - Building graphs #
# 3a ROW 1  Prominent compounds in basil####
# Linalool, eugenol, Methyl eugenol, # eucalyptol 
# note: a-bergamotene & estragole 2nd + 3rd, but only found in 1 var.

# linalool NS
if(TRUE){
  summary(lm(`Linalool floral` ~ Linalool, data = chems)) #ns
  
  linalool <- ggplot(chems, aes(x=Linalool, y = `Linalool floral`))+
    geom_point() +
    geom_smooth(method = lm, color = "blue", se = TRUE,
                linetype = "dashed") +  
    annotate(geom = "text", 
             label = paste0("Rs=", 
                            eq(chems$Linalool, 
                               chems$`Linalool floral`)),
             x = Inf,
             y = Inf,
             hjust = 1,
             vjust = 1,
             fontface = "bold", 
             size = 4,
             color = "red") + 
    geom_text_repel(data = chems, 
                    mapping = aes(x = Linalool,
                                  `Linalool floral`,
                                  label = variety),
                    fontface = "bold",
                    min.segment.length = 1,
                    size = 4)
  linalool <- linalool + labs(x = "Linalool", 
                              y = "")
  linalool <- linalool + theme(axis.title = element_text(size = 12, face = "bold"))
  linalool <- linalool + theme(axis.text = element_text(size = 10))
  
  
  linalool <- linalool + theme(legend.position = "none")
  linalool <- linalool + theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank())
  linalool
}
# eugenol NS
if(TRUE){
summary(lm(`Eugenol floral` ~ Eugenol, data = chems)) #ns

eugenol <- ggplot(chems, aes(x=Eugenol, y = `Eugenol floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$Eugenol, 
                             chems$`Eugenol floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                  mapping = aes(x = Eugenol,
                                y= `Eugenol floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
eugenol <- eugenol + labs(x = "Eugenol", 
                          y = "")
eugenol <- eugenol + theme(axis.title = element_text(size = 12, face = "bold"))
eugenol <- eugenol + theme(axis.text = element_text(size = 10))


eugenol <- eugenol + theme(legend.position = "none")
eugenol <- eugenol + theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())

eugenol}
# methyl eugenol NS
if(TRUE){
summary(lm(`Methyl eugenol floral` ~ `Methyl eugenol`, data = chems)) #ns

M_eugenol <- ggplot(chems, aes(x=`Methyl eugenol`, y = `Methyl eugenol floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`Methyl eugenol`, 
                             chems$`Methyl eugenol floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                  mapping = aes(x = `Methyl eugenol`,
                                `Methyl eugenol floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 10,
                  size = 4)
M_eugenol <- M_eugenol + labs(x = "Methyl eugenol", 
                              y = "")
M_eugenol <- M_eugenol + theme(axis.title = element_text(size = 12, face = "bold"))
M_eugenol <- M_eugenol + theme(axis.text = element_text(size = 10))


M_eugenol <- M_eugenol + theme(legend.position = "none")
M_eugenol <- M_eugenol + theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank())

M_eugenol}
# eucalyptol
if(TRUE){
  summary(lm(`Eucalyptol floral` ~ Eucalyptol, data = chems)) #ns
  
  eucalyptol <- ggplot(chems, aes(x=Eucalyptol, y = `Eucalyptol floral`))+
    geom_point() +
    geom_smooth(method = lm, color = "blue", se = TRUE,
                linetype = "dashed") +  
    annotate(geom = "text", 
             label = paste0("Rs=", 
                            eq(chems$Eucalyptol, 
                               chems$`Eucalyptol floral`)),
             x = Inf,
             y = Inf,
             hjust = 1,
             vjust = 1,
             fontface = "bold", 
             size = 4,
             color = "red") + 
    geom_text_repel(data = chems, 
                    mapping = aes(x = Eucalyptol,
                                  `Eucalyptol floral`,
                                  label = variety),
                    fontface = "bold",
                    min.segment.length = 1,
                    size = 4)
  eucalyptol <- eucalyptol + labs(x = "Eucalyptol", 
                                  y = "")
  eucalyptol <- eucalyptol + theme(axis.title = element_text(size = 12, face = "bold"))
  eucalyptol <- eucalyptol + theme(axis.text = element_text(size = 10))
  
  
  eucalyptol <- eucalyptol + theme(legend.position = "none")
  eucalyptol <- eucalyptol + theme(axis.line = element_line(colour = "black"),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.border = element_blank(),
                                   panel.background = element_blank())
  eucalyptol}






# 3b ROW 2: unk 97/111/70/71, caryophyllene, B ocimene, unk 95/161/107/93 ####
#unknown 97, 111, 70, 71
if(TRUE){
summary(lm(`Unknown 97, 111, 70, 71 floral` ~ `Unknown 97, 111, 70, 71`, data = chems)) #ns
  
unknown43a <- ggplot(chems, aes(x=`Unknown 97, 111, 70, 71`, y = `Unknown 97, 111, 70, 71 floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`Unknown 97, 111, 70, 71`, 
                             chems$`Unknown 97, 111, 70, 71 floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
          color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `Unknown 97, 111, 70, 71`,
                                `Unknown 97, 111, 70, 71 floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
unknown43a <- unknown43a + labs(x = "m/z 97, 111, 70, 71", 
                                y = "")
unknown43a <- unknown43a + theme(axis.title = element_text(size = 12, face = "bold"))
unknown43a <- unknown43a + theme(axis.text = element_text(size = 10))
  
unknown43a <- unknown43a + theme(legend.position = "none")
unknown43a <- unknown43a + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank())
unknown43a
}

#caryophyllene NS
if(TRUE){
summary(lm(`β-Caryophyllene floral` ~ `β-Caryophyllene`, data = chems)) #sig

caryophyllene <- ggplot(chems, aes(x=`β-Caryophyllene`, 
                                   y = `β-Caryophyllene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`β-Caryophyllene`, 
                             chems$`β-Caryophyllene floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `β-Caryophyllene`,
                                `β-Caryophyllene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 0,
                  size = 4)
caryophyllene <- caryophyllene + labs(x = "beta-Caryophyllene", 
                                      y = "")
caryophyllene <- caryophyllene + theme(axis.title = element_text(size = 12, face = "bold"))
caryophyllene <- caryophyllene + theme(axis.text = element_text(size = 10))


caryophyllene <- caryophyllene + theme(legend.position = "none")
caryophyllene <- caryophyllene + theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_blank(),
                                       panel.background = element_blank())
caryophyllene
}

# b-Ocimene NS
if(TRUE){
summary(lm(`(Z)-β-Ocimene floral` ~ `(Z)-β-Ocimene`, data = chems))# ns

b_ocimene <- ggplot(chems, aes(x=`(Z)-β-Ocimene`, y = `(Z)-β-Ocimene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`(Z)-β-Ocimene`, 
                             chems$`(Z)-β-Ocimene floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `(Z)-β-Ocimene`,
                                `(Z)-β-Ocimene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 10,
                  size = 4)
b_ocimene <- b_ocimene + labs(x = "(Z)-beta-Ocimene",
                              y = "")
b_ocimene <- b_ocimene + theme(axis.title = element_text(size = 12, face = "bold"))
b_ocimene <- b_ocimene + theme(axis.text = element_text(size = 10))


b_ocimene <- b_ocimene + theme(legend.position = "none")
b_ocimene <- b_ocimene + theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank())
b_ocimene
}

#unknown 95, 161, 107, 93 NS
if(TRUE){
summary(lm(`Unknown 95, 161, 107, 93 floral` ~ `Unknown 95, 161, 107, 93`, data = chems))# ns
  
unknown95 <- ggplot(chems, aes(x=`Unknown 95, 161, 107, 93`, y = `Unknown 95, 161, 107, 93 floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
         label = paste0("Rs=", 
                          eq(chems$`Unknown 95, 161, 107, 93`, 
                               chems$`Unknown 95, 161, 107, 93 floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
          color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `Unknown 95, 161, 107, 93`,
                                `Unknown 95, 161, 107, 93 floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
unknown95 <- unknown95 + labs(x = "m/z 95, 161, 107, 93", 
                              y = "")
unknown95 <- unknown95 + theme(axis.title = element_text(size = 12, face = "bold"))
unknown95 <- unknown95 + theme(axis.text = element_text(size = 10))
  
unknown95 <- unknown95 + theme(legend.position = "none")
unknown95 <- unknown95 + theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank())
unknown95}

# 3c ROW 3: epi-bicyclo, unk 105/205/71/81, a cubebene, TPEP ####
# epi-bicyclo NS
if(TRUE){
summary(lm(`epi-Bicyclosesquiphellandrene floral` ~ 
             `epi-Bicyclosesquiphellandrene`, data = chems))  
  
epib <- ggplot(chems, aes(x=`epi-Bicyclosesquiphellandrene`,
                          y = `epi-Bicyclosesquiphellandrene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
   annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`epi-Bicyclosesquiphellandrene`, 
                             chems$`epi-Bicyclosesquiphellandrene floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                   mapping = aes(x = `epi-Bicyclosesquiphellandrene`,
                                `epi-Bicyclosesquiphellandrene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
epib <- epib + labs(x = "epi-Bicyclosesquiphellandrene", 
                    y = "")
epib <- epib   + theme(axis.title = element_text(size = 9, face = "bold"))
epib <- epib  + theme(axis.text = element_text(size = 10))
  
  
epib <- epib + theme(legend.position = "none")
epib <- epib  + theme(axis.line = element_line(colour = "black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank())
epib}
# unk 105/205/71/81 NS
if(TRUE){
summary(lm(`Unknown 105, 204, 71, 81 floral` ~ 
             `Unknown 105, 204, 71, 81`, data = chems)) 
  
unknown105 <- ggplot(chems, aes(x=`Unknown 105, 204, 71, 81`, 
                                y = `Unknown 105, 204, 71, 81 floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`Unknown 105, 204, 71, 81`, 
                             chems$`Unknown 105, 204, 71, 81 floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `Unknown 105, 204, 71, 81`,
                                `Unknown 105, 204, 71, 81 floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
unknown105 <- unknown105 + labs(x = "m/z 105, 204, 71, 81", 
                                  y = "")
unknown105 <- unknown105 + theme(axis.title = element_text(size = 12, face = "bold"))
unknown105 <- unknown105 + theme(axis.text = element_text(size = 10))
  
unknown105 <- unknown105 + theme(legend.position = "none")
unknown105 <- unknown105 + theme(axis.line = element_line(colour = "black"),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.border = element_blank(),
                                   panel.background = element_blank())
unknown105
}
# cubebene NS
if(TRUE){
summary(lm(`ɑ-Cubebene floral` ~ `ɑ-Cubebene`, data = chems)) #ns 

cubebene <- ggplot(chems, aes(x=`ɑ-Cubebene`, y = `ɑ-Cubebene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`ɑ-Cubebene`, 
                             chems$`ɑ-Cubebene floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                  mapping = aes(x = `ɑ-Cubebene`,
                                `ɑ-Cubebene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
cubebene <- cubebene + labs(x = "alpha-Cubebene", 
                                y = "")
cubebene <- cubebene   + theme(axis.title = element_text(size = 12, face = "bold"))
cubebene <- cubebene   + theme(axis.text = element_text(size = 10))


cubebene <- cubebene + theme(legend.position = "none")
cubebene <- cubebene  + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank())
cubebene} 
# 2-(1-phenylethyl)-phenol NS
if(TRUE){
summary(lm(`2-(1-phenylethyl)-Phenol floral` ~ `2-(1-phenylethyl)-Phenol`, data = chems))# ns

TPEP <- ggplot(chems, aes(x=`2-(1-phenylethyl)-Phenol`, 
                          y = `2-(1-phenylethyl)-Phenol floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`2-(1-phenylethyl)-Phenol`, 
                             chems$`2-(1-phenylethyl)-Phenol floral`)),
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                  mapping = aes(x = `2-(1-phenylethyl)-Phenol`,
                                `2-(1-phenylethyl)-Phenol floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
TPEP <- TPEP + labs(x = "2-(1-phenylethyl)-phenol", 
                    y = "")
TPEP <- TPEP    + theme(axis.title = element_text(size = 10, face = "bold"))
TPEP <- TPEP   + theme(axis.text = element_text(size = 10))


TPEP <- TPEP   + theme(legend.position = "none")
TPEP <- TPEP   + theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank())
TPEP
}

# 3d ROW 4: d muurolene, a cadinene, bornyl acetate,  a guaiene####
# d-muurolene NS
if(TRUE){
summary(lm(`δ-Muurolene floral` ~ `δ-Muurolene`, data = chems)) 

d_muurolene <- ggplot(chems, aes(x=`δ-Muurolene`, y = `δ-Muurolene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`δ-Muurolene`, 
                             chems$`δ-Muurolene floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") + 
  geom_text_repel(data = chems, 
                  mapping = aes(x = `δ-Muurolene`,
                                `δ-Muurolene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 10,
                  size = 4)
d_muurolene <- d_muurolene + labs(x = "delta-Muurolene", 
                                  y = "")
d_muurolene <- d_muurolene + theme(axis.title = element_text(size = 12, face = "bold"))
d_muurolene <- d_muurolene + theme(axis.text = element_text(size = 10))

d_muurolene <- d_muurolene + theme(legend.position = "none")
d_muurolene <- d_muurolene + theme(axis.line = element_line(colour = "black"),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.border = element_blank(),
                                   panel.background = element_blank())
d_muurolene
}
#a-cadinene NS
if(TRUE){
summary(lm(`ɑ-Cadinene floral` ~ `ɑ-Cadinene`, data = chems))
  
a_cadinene <- ggplot(chems, aes(x=`ɑ-Cadinene`, y = `ɑ-Cadinene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`ɑ-Cadinene`, 
                             chems$`ɑ-Cadinene floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `ɑ-Cadinene`,
                                `ɑ-Cadinene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
a_cadinene <- a_cadinene + labs(x = "alpha-Cadinene", 
                                 y = "")
a_cadinene <- a_cadinene + theme(axis.title = element_text(size = 12, face = "bold"))
a_cadinene <- a_cadinene + theme(axis.text = element_text(size = 10))
  
a_cadinene <- a_cadinene + theme(legend.position = "none")
a_cadinene <- a_cadinene + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank())
a_cadinene}
# bornyl acetate NS
if(TRUE){
summary(lm(`Bornyl acetate floral` ~ `Bornyl acetate`, data = chems)) 
  
bornyl_acetate <- ggplot(chems, aes(x=`Bornyl acetate`, y = `Bornyl acetate floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`Bornyl acetate`, 
                             chems$`Bornyl acetate floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `Bornyl acetate`,
                                `Bornyl acetate floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 10,
                  size = 4)
bornyl_acetate <- bornyl_acetate + labs(x = "Bornyl acetate", 
                                        y = "")
bornyl_acetate <- bornyl_acetate   + theme(axis.title = element_text(size = 12, face = "bold"))
bornyl_acetate <- bornyl_acetate  + theme(axis.text = element_text(size = 12))
  
bornyl_acetate <- bornyl_acetate + theme(legend.position = "none")
bornyl_acetate <- bornyl_acetate + theme(axis.line = element_line(colour = "black"),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         panel.border = element_blank(),
                                         panel.background = element_blank())
bornyl_acetate}
#a-guaiene NS
if(TRUE){
summary(lm(`ɑ-Guaiene floral` ~ `ɑ-Guaiene`, data = chems)) 
  
a_guaiene <- ggplot(chems, aes(x=`ɑ-Guaiene`, y = `ɑ-Guaiene floral`))+
  geom_point() +
  geom_smooth(method = lm, color = "blue", se = TRUE,
              linetype = "dashed") +  
  annotate(geom = "text", 
           label = paste0("Rs=", 
                          eq(chems$`ɑ-Guaiene`, 
                             chems$`ɑ-Guaiene floral`)),
           x = Inf,
           y = Inf,
           hjust = 1,
           vjust = 1,
           fontface = "bold", 
           size = 4,
           color = "red") +
  geom_text_repel(data = chems, 
                  mapping = aes(x = `ɑ-Guaiene`,
                                `ɑ-Guaiene floral`,
                                label = variety),
                  fontface = "bold",
                  min.segment.length = 1,
                  size = 4)
a_guaiene <- a_guaiene + labs(x = "alpha-Guaiene", 
                                 y = "")
a_guaiene <- a_guaiene + theme(axis.title = element_text(size = 12, face = "bold"))
a_guaiene <- a_guaiene + theme(axis.text = element_text(size = 10))
  
a_guaiene <- a_guaiene + theme(legend.position = "none")
a_guaiene <- a_guaiene + theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank())
a_guaiene}

# 4- Exporting ####
# stitching it all together
regressions <- ggarrange(linalool, eugenol, M_eugenol, eucalyptol,
                         unknown43a, caryophyllene, b_ocimene, unknown95,
                         epib, unknown105, cubebene, TPEP, 
                         d_muurolene, a_cadinene, bornyl_acetate, a_guaiene,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                                    "M", "N", "O", "P"),
                         font.label = list(size = 16, color = "black", face = "bold", family = NULL),
                         ncol = 4, nrow = 4)

ggsave(plot = regressions,
       filename = "Fig_A2.pdf",
       width = 3000,
       height = 2700, 
       units = "px")
