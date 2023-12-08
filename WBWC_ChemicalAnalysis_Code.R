# AUTHOR: Yeraldi Loera, yeraldiloera@gmail.com
# CREATION DATE 20 Jun 2022
# PURPOSE: Wedgebilled woodcreeper feather metal analysis 
# DETAILS: Metal analysis by Dartmouth trace metal analysis lab
# QUESTION: Are the feathers from Tiputini and French Guiana contaminated despite their protection?

####SET UP####

#PACKAGES#
library(readxl)
library(dplyr)
library(tidyverse)
library(corrr)
library(corrplot)
library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(AICcmodavg)
library(broom)
library(FactoMineR)
library(factoextra)
library(rstatix)
library(sjstats)
library(survey)
library(sp)
library(amt)
library(rgdal)
library(raster)
library(rasterVis)
library(rgeos)

####DATA OVERVIEW####
#upload full database (both batches, meta data and results combined)
full_data<-as.data.frame(read.csv("Full_Metal_results_RAnalysis.csv"))

#take the mean of all heavy metal concentrations per individual
full_data$Full_load<-rowSums(full_data[,6:28])

#look over database
dim(full_data) 

#23 metals
metals<-colnames(full_data)[6:28]
metals
length(metals)

#metadata
metadata<-colnames(full_data)[29:42]
metadata
#note: I recovered the weights by looking across all dataframes and merging by Field ID. Also collected lat long or looked them up online (Nouragues and San Rafael)

#rename sites with long names
full_data[full_data=="Sangay National Park Macas"] <- "Sangay"
full_data[full_data=="San Rafael Falls"] <- "San Rafael"

full_data



####NORMALITY TEST FOR FULL LOAD####

#density plot not normal tail
ggdensity(full_data$Full_load, 
          main = "Density plot of Full HM Load",
          xlab = "Full HM Load")

#qq plot not notmal
ggqqplot(full_data$Full_load,
         main = "qqplot of Full HM Load",
         xlab = "Full HM Load")

#shapiro-wilkes tests sig different so not normal...
shapiro.test(full_data$Full_load)

#f test for unequal variances show EQUAL variances surprisingly...
var.test(Full_load ~ Experimental.Notes, full_data, 
         alternative = "two.sided")

####AIC TEST FOR AOV BY GROUP (test vs control) = NEEDS NON-PARAMETRIC KRUSKAL####
#by full metal load... 

one.way <- aov(Full_load ~ Experimental.Notes, data = full_data)
summary(one.way)

two.way <- aov(Full_load ~ Experimental.Notes + Feather.Weight, data = full_data)
summary(two.way)

interaction <- aov(Full_load ~ Experimental.Notes*Feather.Weight, data = full_data)
summary(interaction)

blocking <- aov(Full_load ~ Experimental.Notes + Feather.Weight + Country.code, data = full_data)
summary(blocking)

#AIC model test (lowest is best)
model.set <- list(one.way, two.way, interaction, blocking)
model.names <- c("one.way", "two.way feather", "interaction", "blocking")
aictab(model.set, modnames = model.names)

#blocking, which includes feather weights and country code best explains the variation... so separate EC and FG and use feather weights as a weighing measure
par(mfrow=c(2,2))
plot(interaction)
par(mfrow=c(1,1))

#Make a group variable to separate FG and EC controls vs tests...
#rename this column for future analyses
full_data$Group<-paste(full_data$Experimental.Notes, full_data$Country.code)



####PCA####

#normalization is done by subtracting its mean and dividing by its standard deviation
numerical_data <- full_data[,6:28]
data_normalized <- scale(numerical_data)
head(data_normalized)

#correlation matrix
corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)

#3. eigenvalues and eigenvectors
# eigenvalue, on the other hand, is a number representing the amount of variance present in the data for a given direction.
data.pca <- princomp(corr_matrix)
summary(data.pca)

#looking at the first two components (loadings)
data.pca$loadings[, 1:2]

#scree plot of % of variance explained by each dimension
fviz_eig(data.pca, addlabels = TRUE)

# Graph of the variables (metals)
fviz_pca_var(data.pca, col.var = "black")

#scree plot of quality of representation in dimensions 1 + 2 
fviz_cos2(data.pca, choice = "var", axes = 1:2)

#plot the quality of representation for each metal on dimensions 1 + 2 
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)#4. selection of PCS


##PCA other way (normalized)
res.pca <- prcomp(full_data[,6:28],  scale = TRUE)

#plot the first two PC's and group by sampling groups
fviz_pca_ind(res.pca, label="var", habillage=full_data$Group,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2")

#plot the eigenvector of each metal on the first two PCS and color by group
fviz_pca_biplot(res.pca, label = "var", habillage=full_data$Group,
                addEllipses=TRUE, ellipse.level=0.95,
                ggtheme = theme_minimal())

#test how each metal is increasing along the first two PCs...




##EC ONLY ##

EC_data<-subset(full_data, full_data$Group != "Contaminated French Guiana")

# #subset out metals not involved in oil (Sn, Sb, Mo, Tl, Cs, Ag)
# EC_data<-EC_data[ ,-c(19,21,22,23,26)]


numerical_data <- EC_data[,6:23]
data_normalized <- scale(numerical_data)
head(data_normalized)
data.pca <- princomp(corr_matrix)
summary(data.pca)
data.pca$loadings[, 1:2]

fviz_eig(data.pca, addlabels = TRUE)

fviz_pca_var(data.pca, col.var = "black")

fviz_cos2(data.pca, choice = "var", axes = 1:2)

fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)#4. selection of PCS

##PCA other way
res.pca <- prcomp(EC_data[,6:23],  scale = TRUE)

fviz_pca_ind(res.pca, label="var", habillage=EC_data$Group,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2")

fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("white", "blue", "red"),repel = TRUE,
             ggtheme = theme_minimal())


fviz_pca_biplot(res.pca, label = "var", alpha.var="contrib", habillage=EC_data$Group,
                addEllipses=TRUE, ellipse.level=0.95,palette=c('royalblue4','gold2'),
                ggtheme = theme_tufte())



##FG##
FG_data<-subset(full_data, full_data$Group != "Contaminated Ecuador")

colnames(FG_data)
# #subset out metals not involved in gold (Ba, Sr, Sn, Mo, Tl, Cs, Ag)
# FG_data<-FG_data[ ,-c(24, 17, 21, 18, 26, 23, 19)]


numerical_data <- FG_data[,6:21]
data_normalized <- scale(numerical_data)
head(data_normalized)
data.pca <- princomp(corr_matrix)
summary(data.pca)
data.pca$loadings[, 1:2]

fviz_eig(data.pca, addlabels = TRUE)

fviz_pca_var(data.pca, col.var = "black")

fviz_cos2(data.pca, choice = "var", axes = 1:2)

fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)#4. selection of PCS

##PCA other way
res.pca <- prcomp(FG_data[,6:21],  scale = TRUE)

fviz_pca_ind(res.pca, label="var", habillage=FG_data$Group,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2")

fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("white", "blue", "red"),repel = TRUE,
             ggtheme = theme_minimal())


pcaplot<-fviz_pca_biplot(res.pca, label = "var", alpha.var="contrib", habillage=FG_data$Group,
                addEllipses=TRUE, ellipse.level=0.95,palette=c('gold2','green4'),
                ggtheme = theme_minimal())



####TEST AND PLOTTING EC and FG SEPARATELY####

EC_data<-subset(full_data, full_data$Group != "Contaminated French Guiana")
FG_data<-subset(full_data, full_data$Group != "Contaminated Ecuador")


#####ECUADOR ONLY#####

#F test for unequal variances show EQUAL variances surprisingly...
var.test(Full_load ~ Group, EC_data, 
         alternative = "two.sided")


##ALL by GROUP#
#All weighted by featherweight 
EC_All_mann_weighted <- weighted_mannwhitney(EC_data, Full_load,Group, Feather.Weight)
EC_All_mann_weighted

pdf("Figures/All_Mangroup.pdf")
ggboxplot(EC_data, x = "Group", y = "Full_load",
          ylab = "Mean Heavy Metal Load", xlab = "Group", title = "Mean Heavy Metal Load Across Ecuador" , subtitle = EC_All_mann_weighted$p.value) + geom_jitter(color="black", size=0.4, alpha=0.9)
dev.off()


##TEST SIG > CONTROL##

##Co##
EC_Co_mann_weighted <- weighted_mannwhitney(EC_data, Co,Group, Feather.Weight)
EC_Co_mann_weighted

EC_Co_Mann<-ggboxplot(EC_data, x = "Group", y = "Co",
          ylab = "Co Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Co_mann_weighted$p.value)


##Fe##
EC_Fe_mann_weighted <- weighted_mannwhitney(EC_data, Fe,Group, Feather.Weight)
EC_Fe_mann_weighted

EC_Fe_Mann<-ggboxplot(EC_data, x = "Group", y = "Fe",
                      ylab = "Fe Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Fe_mann_weighted$p.value)


##Ni##
EC_Ni_mann_weighted <- weighted_mannwhitney(EC_data, Ni,Group, Feather.Weight)
EC_Ni_mann_weighted

EC_Ni_Mann<-ggboxplot(EC_data, x = "Group", y = "Ni",
                      ylab = "Ni Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Ni_mann_weighted$p.value)


##V##
EC_V_mann_weighted <- weighted_mannwhitney(EC_data, V,Group, Feather.Weight)
EC_V_mann_weighted

EC_V_Mann<-ggboxplot(EC_data, x = "Group", y = "V",
                      ylab = "V Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_V_mann_weighted$p.value)


##Cu##
EC_Cu_mann_weighted <- weighted_mannwhitney(EC_data, Cu,Group, Feather.Weight)
EC_Cu_mann_weighted

EC_Cu_Mann<-ggboxplot(EC_data, x = "Group", y = "Cu",
                      ylab = "Cu Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Cu_mann_weighted$p.value)


##Ba##
EC_Ba_mann_weighted <- weighted_mannwhitney(EC_data, Ba,Group, Feather.Weight)
EC_Ba_mann_weighted

EC_Ba_Mann<-ggboxplot(EC_data, x = "Group", y = "Ba",
                      ylab = "Ba Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Ba_mann_weighted$p.value)


##Cd##
EC_Cd_mann_weighted <- weighted_mannwhitney(EC_data, Cd,Group, Feather.Weight)
EC_Cd_mann_weighted

EC_Cd_Mann<-ggboxplot(EC_data, x = "Group", y = "Cd",
                      ylab = "Cd Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Cd_mann_weighted$p.value)


##Mn##
EC_Mn_mann_weighted <- weighted_mannwhitney(EC_data, Mn,Group, Feather.Weight)
EC_Mn_mann_weighted

EC_Mn_Mann<-ggboxplot(EC_data, x = "Group", y = "Mn",
                      ylab = "Mn Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Mn_mann_weighted$p.value)


##U##
EC_U_mann_weighted <- weighted_mannwhitney(EC_data, U,Group, Feather.Weight)
EC_U_mann_weighted

EC_U_Mann<-ggboxplot(EC_data, x = "Group", y = "U",
                      ylab = "U Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_U_mann_weighted$p.value)


##CONTROL SIG > TEST##

##Zn##
EC_Zn_mann_weighted <- weighted_mannwhitney(EC_data, Zn,Group, Feather.Weight)
EC_Zn_mann_weighted

EC_Zn_Mann<-ggboxplot(EC_data, x = "Group", y = "Zn",
                      ylab = "Zn Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Zn_mann_weighted$p.value)


##As##
EC_As_mann_weighted <- weighted_mannwhitney(EC_data, As,Group, Feather.Weight)
EC_As_mann_weighted

EC_As_Mann<-ggboxplot(EC_data, x = "Group", y = "As",
                      ylab = "As Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_As_mann_weighted$p.value)


##Se##
EC_Se_mann_weighted <- weighted_mannwhitney(EC_data, Se,Group, Feather.Weight)
EC_Se_mann_weighted

EC_Se_Mann<-ggboxplot(EC_data, x = "Group", y = "Se",
                      ylab = "Se Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Se_mann_weighted$p.value)


##Sb##
EC_Sb_mann_weighted <- weighted_mannwhitney(EC_data, Sb,Group, Feather.Weight)
EC_Sb_mann_weighted

EC_Sb_Mann<-ggboxplot(EC_data, x = "Group", y = "Sb",
                      ylab = "Sb Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Sb_mann_weighted$p.value)


##Tl##
EC_Tl_mann_weighted <- weighted_mannwhitney(EC_data, Tl,Group, Feather.Weight)
EC_Tl_mann_weighted

EC_Tl_Mann<-ggboxplot(EC_data, x = "Group", y = "Tl",
                      ylab = "Tl Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Tl_mann_weighted$p.value)


##Cs##
EC_Cs_mann_weighted <- weighted_mannwhitney(EC_data, Cs,Group, Feather.Weight)
EC_Cs_mann_weighted

EC_Cs_Mann<-ggboxplot(EC_data, x = "Group", y = "Cs",
                      ylab = "Cs Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Cs_mann_weighted$p.value)


##NOT SIG##

##Pb##
EC_Pb_mann_weighted <- weighted_mannwhitney(EC_data, Pb,Group, Feather.Weight)
EC_Pb_mann_weighted

EC_Pb_Mann<-ggboxplot(EC_data, x = "Group", y = "Pb",
                      ylab = "Pb Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Pb_mann_weighted$p.value)


##Al##
EC_Al_mann_weighted <- weighted_mannwhitney(EC_data, Al,Group, Feather.Weight)
EC_Al_mann_weighted

EC_Al_Mann<-ggboxplot(EC_data, x = "Group", y = "Al",
                      ylab = "Al Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Al_mann_weighted$p.value)


##Hg##
EC_Hg_mann_weighted <- weighted_mannwhitney(EC_data, Hg,Group, Feather.Weight)
EC_Hg_mann_weighted

EC_Mn_Mann<-ggboxplot(EC_data, x = "Group", y = "Mn",
                      ylab = "Mn Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Mn_mann_weighted$p.value)


##Cr##
EC_Cr_mann_weighted <- weighted_mannwhitney(EC_data, Cr,Group, Feather.Weight)
EC_Cr_mann_weighted

EC_Cr_Mann<-ggboxplot(EC_data, x = "Group", y = "Cr",
                      ylab = "Cr Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Cr_mann_weighted$p.value)


##Sr##
EC_Sr_mann_weighted <- weighted_mannwhitney(EC_data, Sr,Group, Feather.Weight)
EC_Sr_mann_weighted

EC_Sr_Mann<-ggboxplot(EC_data, x = "Group", y = "Sr",
                      ylab = "Sr Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Sr_mann_weighted$p.value)


##Mo##
EC_Mo_mann_weighted <- weighted_mannwhitney(EC_data, Mo,Group, Feather.Weight)
EC_Mo_mann_weighted

EC_Mo_Mann<-ggboxplot(EC_data, x = "Group", y = "Mo",
                      ylab = "Mo Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Mo_mann_weighted$p.value)


##Ag##
EC_Ag_mann_weighted <- weighted_mannwhitney(EC_data, Ag,Group, Feather.Weight)
EC_Ag_mann_weighted

EC_Ag_Mann<-ggboxplot(EC_data, x = "Group", y = "Ag",
                      ylab = "Ag Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Ag_mann_weighted$p.value)


##Sn##
EC_Sn_mann_weighted <- weighted_mannwhitney(EC_data, Sn,Group, Feather.Weight)
EC_Sn_mann_weighted

EC_Sn_Mann<-ggboxplot(EC_data, x = "Group", y = "Sn",
                      ylab = "Sn Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = EC_Sn_mann_weighted$p.value)



#####FG AND FG CONTROL ONLY#####

#f test for unequal variances show EQUAL variances surprisingly...
var.test(Full_load ~ Group, FG_data, 
         alternative = "two.sided")


##ALL by GROUP#
#All weighted by featherweight 
FG_All_mann_weighted <- weighted_mannwhitney(FG_data, Full_load,Group, Feather.Weight)
FG_All_mann_weighted

pdf("Figures/All_MangroupFG.pdf")
ggboxplot(FG_data, x = "Group", y = "Full_load",
          ylab = "Mean Heavy Metal Load", xlab = "Group", title = "Mean Heavy Metal Load" , subtitle = FG_All_mann_weighted$p.value) + geom_jitter(color="black", size=0.4, alpha=0.9)
dev.off()


##TEST SIG > CONTROL##

##Ni##
FG_Ni_mann_weighted <- weighted_mannwhitney(FG_data, Ni,Group, Feather.Weight)
FG_Ni_mann_weighted

FG_Ni_Mann<-ggboxplot(FG_data, x = "Group", y = "Ni",
                      ylab = "Ni Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Ni_mann_weighted$p.value)


##Pb##
FG_Pb_mann_weighted <- weighted_mannwhitney(FG_data, Pb,Group, Feather.Weight)
FG_Pb_mann_weighted

FG_Pb_Mann<-ggboxplot(FG_data, x = "Group", y = "Pb",
                      ylab = "Pb Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Pb_mann_weighted$p.value)


##Cu##
FG_Cu_mann_weighted <- weighted_mannwhitney(FG_data, Cu,Group, Feather.Weight)
FG_Cu_mann_weighted

FG_Cu_Mann<-ggboxplot(FG_data, x = "Group", y = "Cu",
                      ylab = "Cu Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Cu_mann_weighted$p.value)


##Zn##
FG_Zn_mann_weighted <- weighted_mannwhitney(FG_data, Zn,Group, Feather.Weight)
FG_Zn_mann_weighted

FG_Zn_Mann<-ggboxplot(FG_data, x = "Group", y = "Zn",
                      ylab = "Zn Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Zn_mann_weighted$p.value)


##Mn##
FG_Mn_mann_weighted <- weighted_mannwhitney(FG_data, Mn,Group, Feather.Weight)
FG_Mn_mann_weighted

FG_Mn_Mann<-ggboxplot(FG_data, x = "Group", y = "Mn",
                      ylab = "Mn Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Mn_mann_weighted$p.value)


##Fe##
FG_Fe_mann_weighted <- weighted_mannwhitney(FG_data, Fe,Group, Feather.Weight)
FG_Fe_mann_weighted

FG_Fe_Mann<-ggboxplot(FG_data, x = "Group", y = "Fe",
                      ylab = "Fe Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Fe_mann_weighted$p.value)


##Cr##
FG_Cr_mann_weighted <- weighted_mannwhitney(FG_data, Cr,Group, Feather.Weight)
FG_Cr_mann_weighted

FG_Cr_Mann<-ggboxplot(FG_data, x = "Group", y = "Cr",
                      ylab = "Cr Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Cr_mann_weighted$p.value)


##Co##
FG_Co_mann_weighted <- weighted_mannwhitney(FG_data, Co,Group, Feather.Weight)
FG_Co_mann_weighted

FG_Co_Mann<-ggboxplot(FG_data, x = "Group", y = "Co",
                      ylab = "Co Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Co_mann_weighted$p.value)


##U##
FG_U_mann_weighted <- weighted_mannwhitney(FG_data, U,Group, Feather.Weight)
FG_U_mann_weighted

FG_U_Mann<-ggboxplot(FG_data, x = "Group", y = "U",
                      ylab = "U Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_U_mann_weighted$p.value)


##Se##
FG_Se_mann_weighted <- weighted_mannwhitney(FG_data, Se,Group, Feather.Weight)
FG_Se_mann_weighted

FG_Se_Mann<-ggboxplot(FG_data, x = "Group", y = "Se",
                      ylab = "Se Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Se_mann_weighted$p.value)


##V##
FG_V_mann_weighted <- weighted_mannwhitney(FG_data, V,Group, Feather.Weight)
FG_V_mann_weighted

FG_V_Mann<-ggboxplot(FG_data, x = "Group", y = "V",
                      ylab = "V Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_V_mann_weighted$p.value)


##Hg##
FG_Hg_mann_weighted <- weighted_mannwhitney(FG_data, Hg,Group, Feather.Weight)
FG_Hg_mann_weighted

FG_Hg_Mann<-ggboxplot(FG_data, x = "Group", y = "Hg",
                      ylab = "Hg Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Hg_mann_weighted$p.value)


##Cs##
FG_Cs_mann_weighted <- weighted_mannwhitney(FG_data, Cs,Group, Feather.Weight)
FG_Cs_mann_weighted

FG_Cs_Mann<-ggboxplot(FG_data, x = "Group", y = "Cs",
                      ylab = "Cs Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Cs_mann_weighted$p.value)



##CONTROL SIG > TEST##

##Sr##
FG_Sr_mann_weighted <- weighted_mannwhitney(FG_data, Sr,Group, Feather.Weight)
FG_Sr_mann_weighted

FG_Sr_Mann<-ggboxplot(FG_data, x = "Group", y = "Sr",
                      ylab = "Sr Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Sr_mann_weighted$p.value)

##Ba##
FG_Ba_mann_weighted <- weighted_mannwhitney(FG_data, Ba,Group, Feather.Weight)
FG_Ba_mann_weighted

FG_Ba_Mann<-ggboxplot(FG_data, x = "Group", y = "Ba",
                      ylab = "Ba Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Ba_mann_weighted$p.value)


##Tl##
FG_Tl_mann_weighted <- weighted_mannwhitney(FG_data, Tl,Group, Feather.Weight)
FG_Tl_mann_weighted

FG_Tl_Mann<-ggboxplot(FG_data, x = "Group", y = "Tl",
                      ylab = "Tl Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Tl_mann_weighted$p.value)


##Sb##
FG_Sb_mann_weighted <- weighted_mannwhitney(FG_data, Sb,Group, Feather.Weight)
FG_Sb_mann_weighted

FG_Sb_Mann<-ggboxplot(FG_data, x = "Group", y = "Sb",
                      ylab = "Sb Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Sb_mann_weighted$p.value)


##As##
FG_As_mann_weighted <- weighted_mannwhitney(FG_data, As,Group, Feather.Weight)
FG_As_mann_weighted

FG_As_Mann<-ggboxplot(FG_data, x = "Group", y = "As",
                      ylab = "As Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_As_mann_weighted$p.value)


##Cd##
FG_Cd_mann_weighted <- weighted_mannwhitney(FG_data, Cd,Group, Feather.Weight)
FG_Cd_mann_weighted

FG_Cd_Mann<-ggboxplot(FG_data, x = "Group", y = "Cd",
                      ylab = "Cd Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Cd_mann_weighted$p.value)


##NOT SIG ##

##Al##
FG_Al_mann_weighted <- weighted_mannwhitney(FG_data, Al,Group, Feather.Weight)
FG_Al_mann_weighted

FG_Al_Mann<-ggboxplot(FG_data, x = "Group", y = "Al",
                      ylab = "Al Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Al_mann_weighted$p.value)


##Mo##
FG_Mo_mann_weighted <- weighted_mannwhitney(FG_data, Mo,Group, Feather.Weight)
FG_Mo_mann_weighted

FG_Mo_Mann<-ggboxplot(FG_data, x = "Group", y = "Mo",
                      ylab = "Mo Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Mo_mann_weighted$p.value)


##Ag##
FG_Ag_mann_weighted <- weighted_mannwhitney(FG_data, Ag,Group, Feather.Weight)
FG_Ag_mann_weighted
FG_Ag_Mann<-ggboxplot(FG_data, x = "Group", y = "Ag",
                      ylab = "Ag Concentration", xlab = "Group") +
  labs(title = "Weighted Mann-Whitney" ,subtitle = FG_Ag_mann_weighted$p.value)



####SUMMARY PLOTTING####

#####PLOTTING 3 GROUPS BY METAL #####
# Visualization: box plots with Mann wallis p-values and wilcox test
# Visualization: box plots with Mann wallis p-values 

#get full load info per indiv
full_data$Group<-paste(full_data$Experimental.Notes, full_data$Country.code)

#full_data$Full_load<-rowSums(full_data[,6:28])

#make a new column of "Sampling Group" with control ecuador renamed as reference...

full_data$Sampling_Group<-full_data$Group

full_data$Sampling_Group[full_data$Sampling_Group == "Control Ecuador"] <- "Reference Ecuador"

full_data$Sampling_Group

###KRUSKAL TEST PLOT###

#Mean HM load (note that "weighted mann whitney" funciton takes into account 3 groups = kruskal wallis!)
Full_mann_weighted <- weighted_mannwhitney(full_data, Full_load,Sampling_Group, Feather.Weight)
Full_mann_weighted

# Pairwise comparisons
Full_wilcox <- full_data %>% 
  wilcox_test(Full_load ~ Sampling_Group) 
Full_wilcox

Full_wilcox <- Full_wilcox %>% add_xy_position(x = "Sampling_Group")

full_data$Group = factor(full_data$Sampling_Group, levels=c("Contaminated Ecuador", "Reference Ecuador", "Contaminated French Guiana"))

FULL_PLOT<-ggboxplot(full_data, x = "Group", y = "Full_load", color="Group",palette=c('royalblue4','gold2','green4'),
          ylab = "Mean Heavy Metal Load (ppm)", xlab = "Group")+ geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Full_wilcox, hide.ns = TRUE) +
  labs( title = "Mean Heavy Metal Load by Group" ,
        subtitle = "Weighted Kruskal-Wallis, p,5.82e-10",
        caption = get_pwc_label(Full_wilcox)
  )+theme(legend.position='upper right')
FULL_PLOT




#Al#
#Al weighted by featherweight = NOTSIG
Al_mann_weighted <- weighted_mannwhitney(full_data, Al,Group, Feather.Weight)
Al_mann_weighted

# Pairwise comparisons
Al_wilcox <- full_data %>% 
  wilcox_test(Al ~ Group) 
Al_wilcox

pdf("Figures/Al_Manngroup.pdf")
Al_wilcox <- Al_wilcox %>% add_xy_position(x = "Group")
AL_PLOT<-ggboxplot(full_data, x = "Group", y = "Al",
          ylab = "Al Concentration (ppm)", xlab = "Group") + geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Al_wilcox, hide.ns = TRUE) +
  labs( title =Al_mann_weighted$method ,
        subtitle = Al_mann_weighted$p.value,
        caption = get_pwc_label(Al_wilcox)
  )
dev.off()


#Mo weighted by featherweight =  NOT SIG
Mo_mann_weighted <- weighted_mannwhitney(full_data, Mo,Group, Feather.Weight)
Mo_mann_weighted

# Pairwise comparisons NOT SIG
Mo_wilcox <- full_data %>% 
  wilcox_test(Mo ~ Group) 
Mo_wilcox

pdf("Figures/Mo_Manngroup.pdf")
Mo_wilcox <- Mo_wilcox %>% add_xy_position(x = "Group")
MO_PLOT<-ggboxplot(full_data, x = "Group", y = "Mo",
          ylab = "Mo Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Mo_wilcox, hide.ns = TRUE) +
  labs( title =Mo_mann_weighted$method ,
        subtitle = Mo_mann_weighted$p.value,
        caption = get_pwc_label(Mo_wilcox)
  )
dev.off()


#Ag weighted by featherweight = NOT SIG
Ag_mann_weighted <- weighted_mannwhitney(full_data, Ag,Group, Feather.Weight)
Ag_mann_weighted


# Pairwise comparisons
Ag_wilcox <- full_data %>% 
  wilcox_test(Ag ~ Group) 
Ag_wilcox

pdf("Figures/Ag_Manngroup.pdf")
Ag_wilcox <- Ag_wilcox %>% add_xy_position(x = "Group")
AG_PLOT<-ggboxplot(full_data, x = "Group", y = "Ag",
          ylab = "Ag Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Ag_wilcox, hide.ns = TRUE) +
  labs( title =Ag_mann_weighted$method ,
        subtitle = Ag_mann_weighted$p.value,
        caption = get_pwc_label(Ag_wilcox)
  )
dev.off()


##rest are sig##


#V weighted by featherweight = SIG
V_mann_weighted <- weighted_mannwhitney(full_data, V,Group, Feather.Weight)
V_mann_weighted

# Pairwise comparisons
V_wilcox <- full_data %>% 
  wilcox_test(V ~ Group) 
V_wilcox

pdf("Figures/V_Manngroup.pdf")
V_wilcox <- V_wilcox %>% add_xy_position(x = "Group")
V_PLOT<-ggboxplot(full_data, x = "Group", y = "V",
          ylab = "V Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(V_wilcox, hide.ns = TRUE) +
  labs( title =V_mann_weighted$method ,
        subtitle = V_mann_weighted$p.value,
        caption = get_pwc_label(V_wilcox)
  )
dev.off()


#Cr weighted by featherweight = SIG
Cr_mann_weighted <- weighted_mannwhitney(full_data, Cr,Group, Feather.Weight)
Cr_mann_weighted

# Pairwise comparisons
Cr_wilcox <- full_data %>% 
  wilcox_test(Cr ~ Group) 
Cr_wilcox

pdf("Figures/Cr_Manngroup.pdf")
Cr_wilcox <- Cr_wilcox %>% add_xy_position(x = "Group")
CR_PLOT<-ggboxplot(full_data, x = "Group", y = "Cr",
          ylab = "Cr Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Cr_wilcox, hide.ns = TRUE) +
  labs( title =Cr_mann_weighted$method ,
        subtitle = Cr_mann_weighted$p.value,
        caption = get_pwc_label(Cr_wilcox)
  )
dev.off()


#Mn weighted by featherweight = SIG
Mn_mann_weighted <- weighted_mannwhitney(full_data, Mn,Group, Feather.Weight)
Mn_mann_weighted


# Pairwise comparisons
Mn_wilcox <- full_data %>% 
  wilcox_test(Mn ~ Group) 
Mn_wilcox

pdf("Figures/Mn_Manngroup.pdf")
Mn_wilcox <- Mn_wilcox %>% add_xy_position(x = "Group")
MN_PLOT<-ggboxplot(full_data, x = "Group", y = "Mn",
          ylab = "Mn Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Mn_wilcox, hide.ns = TRUE) +
  labs( title =Mn_mann_weighted$method ,
        subtitle = Mn_mann_weighted$p.value,
        caption = get_pwc_label(Mn_wilcox)
  )
dev.off()

#Fe weighted by featherweight = SIG
Fe_mann_weighted <- weighted_mannwhitney(full_data, Fe,Group, Feather.Weight)
Fe_mann_weighted

# Pairwise comparisons
Fe_wilcox <- full_data %>% 
  wilcox_test(Fe ~ Group) 
Fe_wilcox

pdf("Figures/Fe_Manngroup.pdf")
Fe_wilcox <- Fe_wilcox %>% add_xy_position(x = "Group")
FE_PLOT<-ggboxplot(full_data, x = "Group", y = "Fe",
          ylab = "Fe Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Fe_wilcox, hide.ns = TRUE) +
  labs( title =Fe_mann_weighted$method ,
        subtitle = Fe_mann_weighted$p.value,
        caption = get_pwc_label(Fe_wilcox)
  )
dev.off()


#Co weighted by featherweight =  SIG
Co_mann_weighted <- weighted_mannwhitney(full_data, Co,Group, Feather.Weight)
Co_mann_weighted

# Pairwise comparisons
Co_wilcox <- full_data %>% 
  wilcox_test(Co ~ Group) 
Co_wilcox

pdf("Figures/Co_Manngroup.pdf")
Co_wilcox <- Co_wilcox %>% add_xy_position(x = "Group")
CO_PLOT<-ggboxplot(full_data, x = "Group", y = "Co",
          ylab = "Co Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Co_wilcox, hide.ns = TRUE) +
  labs( title =Co_mann_weighted$method ,
        subtitle = Co_mann_weighted$p.value,
        caption = get_pwc_label(Co_wilcox)
  )
dev.off()

#Ni weighted by featherweight =  SIG
Ni_mann_weighted <- weighted_mannwhitney(full_data, Ni,Group, Feather.Weight)
Ni_mann_weighted

# Pairwise comparisons
Ni_wilcox <- full_data %>% 
  wilcox_test(Ni ~ Group) 
Ni_wilcox

pdf("Figures/Ni_Manngroup.pdf")
Ni_wilcox <- Ni_wilcox %>% add_xy_position(x = "Group")
NI_PLOT<-ggboxplot(full_data, x = "Group", y = "Ni",
          ylab = "Ni Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Ni_wilcox, hide.ns = TRUE) +
  labs( title =Ni_mann_weighted$method ,
        subtitle = Ni_mann_weighted$p.value,
        caption = get_pwc_label(Ni_wilcox)
  )
dev.off()


#Cu weighted by featherweight =  SIG
Cu_mann_weighted <- weighted_mannwhitney(full_data, Cu,Group, Feather.Weight)
Cu_mann_weighted

# Pairwise comparisons
Cu_wilcox <- full_data %>% 
  wilcox_test(Cu ~ Group) 
Cu_wilcox

pdf("Figures/Cu_Manngroup.pdf")
Cu_wilcox <- Cu_wilcox %>% add_xy_position(x = "Group")
CU_PLOT<-ggboxplot(full_data, x = "Group", y = "Cu",
          ylab = "Cu Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Cu_wilcox, hide.ns = TRUE) +
  labs( title =Cu_mann_weighted$method ,
        subtitle = Cu_mann_weighted$p.value,
        caption = get_pwc_label(Cu_wilcox)
  )
dev.off()


#Zn weighted by featherweight =  SIG
Zn_mann_weighted <- weighted_mannwhitney(full_data, Zn,Group, Feather.Weight)
Zn_mann_weighted

# Pairwise comparisons
Zn_wilcox <- full_data %>% 
  wilcox_test(Zn ~ Group) 
Zn_wilcox

pdf("Figures/Zn_Manngroup.pdf")
Zn_wilcox <- Zn_wilcox %>% add_xy_position(x = "Group")
ZN_PLOT<-ggboxplot(full_data, x = "Group", y = "Zn",
          ylab = "Zn Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Zn_wilcox, hide.ns = TRUE) +
  labs( title =Zn_mann_weighted$method ,
        subtitle = Zn_mann_weighted$p.value,
        caption = get_pwc_label(Zn_wilcox)
  )
dev.off()


#Se weighted by featherweight =  SIG
Se_mann_weighted <- weighted_mannwhitney(full_data, Se,Group, Feather.Weight)
Se_mann_weighted

# Pairwise comparisons
Se_wilcox <- full_data %>% 
  wilcox_test(Se ~ Group) 
Se_wilcox

pdf("Figures/Se_Manngroup.pdf")
Se_wilcox <- Se_wilcox %>% add_xy_position(x = "Group")
SE_PLOT<-ggboxplot(full_data, x = "Group", y = "Se",
          ylab = "Se Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Se_wilcox, hide.ns = TRUE) +
  labs( title =Se_mann_weighted$method ,
        subtitle = Se_mann_weighted$p.value,
        caption = get_pwc_label(Se_wilcox)
  )
dev.off()


#Sr weighted by featherweight =  SIG
Sr_mann_weighted <- weighted_mannwhitney(full_data, Sr,Group, Feather.Weight)
Sr_mann_weighted

# Pairwise comparisons
Sr_wilcox <- full_data %>% 
  wilcox_test(Sr ~ Group) 
Sr_wilcox

pdf("Figures/Sr_Manngroup.pdf")
Sr_wilcox <- Sr_wilcox %>% add_xy_position(x = "Group")
SR_PLOT<-ggboxplot(full_data, x = "Group", y = "Sr",
          ylab = "Sr Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Sr_wilcox, hide.ns = TRUE) +
  labs( title =Sr_mann_weighted$method ,
        subtitle = Sr_mann_weighted$p.value,
        caption = get_pwc_label(Sr_wilcox)
  )
dev.off()


#Cd weighted by featherweight =  SIG
Cd_mann_weighted <- weighted_mannwhitney(full_data, Cd,Group, Feather.Weight)
Cd_mann_weighted

# Pairwise comparisons
Cd_wilcox <- full_data %>% 
  wilcox_test(Cd ~ Group) 
Cd_wilcox

pdf("Figures/Cd_Manngroup.pdf")
Cd_wilcox <- Cd_wilcox %>% add_xy_position(x = "Group")
CD_PLOT<-ggboxplot(full_data, x = "Group", y = "Cd",
          ylab = "Cd Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Cd_wilcox, hide.ns = TRUE) +
  labs( title =Cd_mann_weighted$method ,
        subtitle = Cd_mann_weighted$p.value,
        caption = get_pwc_label(Cd_wilcox)
  )
dev.off()


#Sn weighted by featherweight =  SIG
Sn_mann_weighted <- weighted_mannwhitney(full_data, Sn,Group, Feather.Weight)
Sn_mann_weighted

# Pairwise comparisons
Sn_wilcox <- full_data %>% 
  wilcox_test(Sn ~ Group) 
Sn_wilcox

pdf("Figures/Sn_Manngroup.pdf")
Sn_wilcox <- Sn_wilcox %>% add_xy_position(x = "Group")
SN_PLOT<-ggboxplot(full_data, x = "Group", y = "Sn",
          ylab = "Sn Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Sn_wilcox, hide.ns = TRUE) +
  labs( title =Sn_mann_weighted$method ,
        subtitle = Sn_mann_weighted$p.value,
        caption = get_pwc_label(Sn_wilcox)
  )
dev.off()


#Sb weighted by featherweight =  SIG
Sb_mann_weighted <- weighted_mannwhitney(full_data, Sb,Group, Feather.Weight)
Sb_mann_weighted

# Pairwise comparisons
Sb_wilcox <- full_data %>% 
  wilcox_test(Sb ~ Group) 
Sb_wilcox

pdf("Figures/Sb_Manngroup.pdf")
Sb_wilcox <- Sb_wilcox %>% add_xy_position(x = "Group")
SB_PLOT<-ggboxplot(full_data, x = "Group", y = "Sb",
          ylab = "Sb Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Sb_wilcox, hide.ns = TRUE) +
  labs( title =Sb_mann_weighted$method ,
        subtitle = Sb_mann_weighted$p.value,
        caption = get_pwc_label(Sb_wilcox)
  )
dev.off()


#Cs weighted by featherweight =  SIG
Cs_mann_weighted <- weighted_mannwhitney(full_data, Cs,Group, Feather.Weight)
Cs_mann_weighted

# Pairwise comparisons
Cs_wilcox <- full_data %>% 
  wilcox_test(Cs ~ Group) 
Cs_wilcox

pdf("Figures/Cs_Manngroup.pdf")
Cs_wilcox <- Cs_wilcox %>% add_xy_position(x = "Group")
CS_PLOT<-ggboxplot(full_data, x = "Group", y = "Cs",
          ylab = "Cs Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Cs_wilcox, hide.ns = TRUE) +
  labs( title =Cs_mann_weighted$method ,
        subtitle = Cs_mann_weighted$p.value,
        caption = get_pwc_label(Cs_wilcox)
  )
dev.off()


#Ba weighted by featherweight =  SIG
Ba_mann_weighted <- weighted_mannwhitney(full_data, Ba,Group, Feather.Weight)
Ba_mann_weighted

# Pairwise comparisons
Ba_wilcox <- full_data %>% 
  wilcox_test(Ba ~ Group) 
Ba_wilcox

pdf("Figures/Ba_Manngroup.pdf")
Ba_wilcox <- Ba_wilcox %>% add_xy_position(x = "Group")
BA_PLOT<-ggboxplot(full_data, x = "Group", y = "Ba",
          ylab = "Ba Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Ba_wilcox, hide.ns = TRUE) +
  labs( title =Ba_mann_weighted$method ,
        subtitle = Ba_mann_weighted$p.value,
        caption = get_pwc_label(Ba_wilcox)
  )
dev.off()

#Hg weighted by featherweight =  SIG
Hg_mann_weighted <- weighted_mannwhitney(full_data, Hg,Group, Feather.Weight)
Hg_mann_weighted

# Pairwise comparisons
Hg_wilcox <- full_data %>% 
  wilcox_test(Hg ~ Group) 
Hg_wilcox

pdf("Figures/Hg_Manngroup.pdf")
Hg_wilcox <- Hg_wilcox %>% add_xy_position(x = "Group")
HG_PLOT<-ggboxplot(full_data, x = "Group", y = "Hg",
          ylab = "Hg Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Hg_wilcox, hide.ns = TRUE) +
  labs( title =Hg_mann_weighted$method ,
        subtitle = Hg_mann_weighted$p.value,
        caption = get_pwc_label(Hg_wilcox)
  )
dev.off()


#Tl weighted by featherweight =  SIG
Tl_mann_weighted <- weighted_mannwhitney(full_data, Tl,Group, Feather.Weight)
Tl_mann_weighted

# Pairwise comparisons
Tl_wilcox <- full_data %>% 
  wilcox_test(Tl ~ Group) 
Tl_wilcox

pdf("Figures/Tl_Manngroup.pdf")
Tl_wilcox <- Tl_wilcox %>% add_xy_position(x = "Group")
TL_PLOT<-ggboxplot(full_data, x = "Group", y = "Tl",
          ylab = "Tl Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Tl_wilcox, hide.ns = TRUE) +
  labs( title =Tl_mann_weighted$method ,
        subtitle = Tl_mann_weighted$p.value,
        caption = get_pwc_label(Tl_wilcox)
  )
dev.off()


#Pb weighted by featherweight =  SIG
Pb_mann_weighted <- weighted_mannwhitney(full_data, Pb,Group, Feather.Weight)
Pb_mann_weighted


# Pairwise comparisons
Pb_wilcox <- full_data %>% 
  wilcox_test(Pb ~ Group) 
Pb_wilcox

pdf("Figures/Pb_Manngroup.pdf")
Pb_wilcox <- Pb_wilcox %>% add_xy_position(x = "Group")
PB_PLOT<-ggboxplot(full_data, x = "Group", y = "Pb",
          ylab = "Pb Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(Pb_wilcox, hide.ns = TRUE) +
  labs( title =Pb_mann_weighted$method ,
        subtitle = Pb_mann_weighted$p.value,
        caption = get_pwc_label(Pb_wilcox)
  )
dev.off()

#U weighted by featherweight =  SIG
U_mann_weighted <- weighted_mannwhitney(full_data, U,Group, Feather.Weight)
U_mann_weighted

# Pairwise comparisons
U_wilcox <- full_data %>% 
  wilcox_test(U ~ Group) 
U_wilcox

pdf("Figures/U_Manngroup.pdf")
U_wilcox <- U_wilcox %>% add_xy_position(x = "Group")
U_PLOT<-ggboxplot(full_data, x = "Group", y = "U",
          ylab = "U Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_pvalue_manual(U_wilcox, hide.ns = TRUE) +
  labs( title =U_mann_weighted$method ,
        subtitle =U_mann_weighted$p.value,
        caption = get_pwc_label(U_wilcox)
  )
dev.off()

# #As weighted by featherweight = MISSING DATA
# As_mann_weighted <- weighted_mannwhitney(full_data, As,Group, Feather.Weight)
# As_mann_weighted
# 
# # Pairwise comparisons
# As_wilcox <- full_data %>% 
#   wilcox_test(As ~ Group) 
# As_wilcox
# 
# pdf("Figures/As_Manngroup.pdf")
# As_wilcox <- As_wilcox %>% add_xy_position(x = "Group")
# AS_PLOT<-ggboxplot(full_data, x = "Group", y = "As",
#           ylab = "As Concentration (ppm)", xlab = "Group") +geom_jitter(color="black", size=0.4, alpha=0.9) +
#   stat_pvalue_manual(As_wilcox, hide.ns = TRUE) +
#   labs( title =As_mann_weighted$method ,
#         subtitle = As_mann_weighted$p.value,
#         caption = get_pwc_label(As_wilcox)
#   )
# dev.off()



####LOOK OVER RESULTS####
metals

#weighted kruskal tests with wilcoxon pairwise analyses with holm-bonferroni corrections
#ignoring relatinship between tiputini and fg

AL_PLOT #notsig
V_PLOT  #sig kruskal, FG>Control
CR_PLOT #sig kruskal, Control>Tip
MN_PLOT #sig kruskal, FG>Control
FE_PLOT #sig kruskal, Tip>Control, FG>control
CO_PLOT #sig kruskal, FG>control 
NI_PLOT #sig kruskal, FG>control Tip> Control
CU_PLOT #sig kruskal, FG>control Tip> Control
ZN_PLOT #sig kruskal, FG>control Control>tip
AS_PLOT #NOT ENOUGH DATA
SE_PLOT #sig kruskal, FG>control Control>tip
SR_PLOT #sig kruskal, Tip> Control
MO_PLOT #not sig
AG_PLOT #not sig
CD_PLOT #sig kruskal, Tip> Control, controll>FG
SN_PLOT #sig kruskal, Tip> Control,
SB_PLOT #sig kruskal, control > tip and Fg
CS_PLOT #sig kruskal, FG> control, control > tip
BA_PLOT #sig kruskal, Tip>control, control>FG
HG_PLOT #sig kruskal, control> tip
TL_PLOT #sig kruskal, contro>Fg
PB_PLOT #sig kruskal, control>tip
U_PLOT #sig kruskal, fg>control

##in summary, 
#tiputini is greater than control for (Fe, Ni, Cu, Sr, Cd, Sn, Ba) 
#control greater than tiputini for (Cr, zn, se, sb, cs, hg, pb)

#Fg greater than control for (v,mn, fe, co, ni, cu, zn, se, cs, u)
#control greater than FG for (cd, sb, ba, tl)


#####PLOTTING ALL EC and FG SEPARATELY#####


##ECUADOR dataset where the means and sd was calculated for each metal per group
EC_plottingdf<-read.csv("ECUADORplottingdataset.csv")
       
##split up into two graphs for better axis...###

sig1<-c("Zn","Mn","Fe","Ba", "Cu")

sig2<-c("Ni","Cd","Se","Co","V","As","U")

sig3<-c("Sb","Tl","Cs")

EC_sigplottingdf1<-filter(EC_plottingdf, metals %in% sig1)
EC_sigplottingdf2<-filter(EC_plottingdf, metals %in% sig2)
EC_sigplottingdf3<-filter(EC_plottingdf, metals %in% sig3)

#level the bars
EC_sigplottingdf1$metals <- factor(EC_sigplottingdf1$metals, levels = sig1)
EC_sigplottingdf2$metals <- factor(EC_sigplottingdf2$metals, levels = sig2)
EC_sigplottingdf3$metals <- factor(EC_sigplottingdf3$metals, levels = sig3)


##first graph
ggboxplot(FG_data, x = "Group", y = "Zn",
          ylab = "Zn Concentration", xlab = "Group") +
  stat_pvalue_manual(FG_Zn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = FG_Zn_kruskal_weighted$p.value,
       caption = get_pwc_label(FG_Zn_wilcox)
  )





pbar1<- ggplot(EC_sigplottingdf1, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(title="Ecuador Wedgebilled Woodcreeper Feather Metal Concentrations (Kruskal-Wallis p<0.05)", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2'))+
  coord_flip()

pbar1


#####Second graph


pbar2<- ggplot(EC_sigplottingdf2, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2'))+
  coord_flip()

pbar2


#####Third graph


pbar3<- ggplot(EC_sigplottingdf3, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2'))+
  coord_flip()

pbar3


##PLOT TOGETHER##

pdf("Figures/EC metalplotseparate.pdf")
ggarrange(pbar1, pbar2, pbar3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()

ggarrange(pbar1, pbar2, pbar3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

###COORD FLIP###

##first graph

pbar1<- ggplot(EC_sigplottingdf1, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(title="Ecuador Wedgebilled Woodcreeper Feather Metal Concentrations (Kruskal-Wallis p<0.05)", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2'))

pbar1


#####Second graph


pbar2<- ggplot(EC_sigplottingdf2, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2'))

pbar2


#####Third graph


pbar3<- ggplot(EC_sigplottingdf3, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2'))

pbar3



ggarrange(pbar1, pbar2, pbar3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()






##FRENCH GUIANA##

FG_plottingdf<-read.csv("FRENCHGUIANAplottingdataset.csv")

##split up into two graphs for better axis...###
sig1<-c("Zn","Mn","Fe","Cu","Ni")

sig2<-c("Cr","Pb","Se","Co","Cd")

sig3<-c("Hg","V","As","U")

sig4<-c("Ba","Sr","Sn","Sb","Tl","Cs")

FG_sigplottingdf1<-filter(FG_plottingdf, metals %in% sig1)
FG_sigplottingdf2<-filter(FG_plottingdf, metals %in% sig2)
FG_sigplottingdf3<-filter(FG_plottingdf, metals %in% sig3)
FG_sigplottingdf4<-filter(FG_plottingdf, metals %in% sig4)


#level the bars
FG_sigplottingdf1$metals <- factor(FG_sigplottingdf1$metals, levels = sig1)
FG_sigplottingdf2$metals <- factor(FG_sigplottingdf2$metals, levels = sig2)
FG_sigplottingdf3$metals <- factor(FG_sigplottingdf3$metals, levels = sig3)
FG_sigplottingdf4$metals <- factor(FG_sigplottingdf4$metals, levels = sig4)


##first graph

pbar1<- ggplot(FG_sigplottingdf1, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(title="French Guiana Wedgebilled Woodcreeper Feather Metal Concentrations (Kruskal-Wallis p<0.05)", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))+
  coord_flip()

pbar1



#####Second graph


pbar2<- ggplot(FG_sigplottingdf2, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))+
  coord_flip()

pbar2

#####Third graph


pbar3<- ggplot(FG_sigplottingdf3, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))+
  coord_flip()

pbar3



pbar4<- ggplot(FG_sigplottingdf4, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))+
  coord_flip()

pbar4

##PLOT TOGETHER##

pdf("Figures/FG metalplotseparate.pdf")
ggarrange(pbar1, pbar2, pbar3, pbar4,
          labels = c("A", "B", "C","D"),
          ncol = 1, nrow = 4)
dev.off()

ggarrange(pbar1, pbar2, pbar3, pbar4,
          labels = c("A", "B", "C","D"),
          ncol = 1, nrow = 4)



##FLIP COORDS##

##first graph

pbar1<- ggplot(FG_sigplottingdf1, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(title="French Guiana Wedgebilled Woodcreeper Feather Metal Concentrations (Kruskal-Wallis p<0.05)", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))

pbar1



#####Second graph


pbar2<- ggplot(FG_sigplottingdf2, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))

pbar2

#####Third graph


pbar3<- ggplot(FG_sigplottingdf3, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))

pbar3



pbar4<- ggplot(FG_sigplottingdf4, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="metals", y="Mean Concentration (ug/g)")+
  theme_classic() +
  scale_fill_manual(values=c('green4','gold2'))

pbar4

##PLOT TOGETHER##

ggarrange(pbar1, pbar2, pbar3, pbar4,
          labels = c("A", "B", "C","D"),
          ncol = 1, nrow = 4)


dev.off()






####COMBINED EC AND FG SUMMARY FOR PAPER####

EC_plottingdf<-read.csv("ECUADORplottingdataset.csv")
FG_plottingdf<-read.csv("FRENCHGUIANAplottingdataset.csv")


EC_plottingdf$Country<-"Ecuador"
EC_plottingdf$Group<-paste(EC_plottingdf$Country,EC_plottingdf$Group)


FG_plottingdf$Country<-"French Guiana"
FG_plottingdf$Group<-paste(FG_plottingdf$Country,FG_plottingdf$Group)
FG_plottingdf_new<-FG_plottingdf[24:46,]
FG_plottingdf_new

combineddf<-rbind(FG_plottingdf_new, EC_plottingdf)


##GRAPH##
#make a new column of "Sampling Group" with control ecuador renamed as reference...
combineddf$Sampling_Group<-combineddf$Group
combineddf$Sampling_Group[combineddf$Sampling_Group == "Ecuador control"] <- "Reference Ecuador"
combineddf$Sampling_Group[combineddf$Sampling_Group == "Ecuador contaminated"] <- "Contaminated Ecuador"
combineddf$Sampling_Group[combineddf$Sampling_Group == "French Guiana contaminated"] <- "Contaminated French Guiana"
combineddf$Sampling_Group

combineddf$Sampling_Group = factor(combined$Sampling_Group, levels=c("Contaminated Ecuador", "Reference Ecuador", "Contaminated French Guiana"))


pbar1<- ggplot(combineddf, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(title="Ecuador Wedgebilled Woodcreeper Feather Metal Concentrations (Kruskal-Wallis p<0.05)", y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))

pbar1



##split up into two graphs for better axis...###
sig1<-c("Zn","Mn","Fe","Al","Cu")

sig2<-c("Ba","Sr","Ni","Cr","Pb")

sig3<-c("Se","Sn","Cd","Co","Hg","Sb")

sig4<-c("V","Mo","Tl","Cs","Ag","As","U")

combineddf1<-filter(combineddf, metals %in% sig1)
combineddf2<-filter(combineddf, metals %in% sig2)
combineddf3<-filter(combineddf, metals %in% sig3)
combineddf4<-filter(combineddf, metals %in% sig4)


# #level the bars
combineddf1$metals <- factor(combineddf1$metals, levels = sig1)
combineddf2$metals <- factor(combineddf2$metals, levels = sig2)
combineddf3$metals <- factor(combineddf3$metals, levels = sig3)
combineddf4$metals <- factor(combineddf4$metals, levels = sig4)


##first graph

pbar1<- ggplot(combineddf1, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))
pbar1



#####Second graph


pbar2<- ggplot(combineddf2, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))
pbar2


#####Third graph


pbar3<- ggplot(combineddf3, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))
pbar3


#####Fourth graph


pbar4<- ggplot(combineddf4, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))
pbar4


ggarrange(pbar1, pbar2, pbar3, pbar4,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

dev.off()




###POINT PLOT

ppoint1<- ggplot(combineddf1, aes(x=metals, y=means, fill=Group)) + 
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) + geom_point(shape=21, size=2, position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))

ppoint1

ppoint2<- ggplot(combineddf2, aes(x=metals, y=means, fill=Group)) + 
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) + geom_point(shape=21, size=2, position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))

ppoint2

ppoint3<- ggplot(combineddf3, aes(x=metals, y=means, fill=Group)) + 
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) + geom_point(shape=21, size=2, position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))

ppoint3

ppoint4<- ggplot(combineddf4, aes(x=metals, y=means, fill=Group)) + 
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) + geom_point(shape=21, size=2, position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('royalblue4','gold2','green4'))

ppoint4


ggarrange(ppoint1, ppoint2, ppoint3, ppoint4,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

dev.off()


#####BW COLORBLIND#####

##first graph

pbar1<- ggplot(combineddf1, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('black','white','darkgrey'))

pbar1


#####Second graph


pbar2<- ggplot(combineddf2, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs( y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('black','white','darkgrey'))
pbar2


#####Third graph


pbar3<- ggplot(combineddf3, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('black','white','darkgrey'))
pbar3


#####Fourth graph


pbar4<- ggplot(combineddf4, aes(x=metals, y=means, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=means-se, ymax=means+se), width=.2,
                position=position_dodge(.9)) +
  labs(y="Mean Concentration (ug/g)", x="Metals")+
  theme_classic() +
  scale_fill_manual(values=c('black','white','darkgrey'))
pbar4


ggarrange(pbar1, pbar2, pbar3, pbar4,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

dev.off()




####SPATIAL DATA FOR GROUPING...TIPUTINI IS HIGH RISK COMPARED TO CONTROLS####


#READ IN OIL DATA FROM OUR WORLD IN DATA (SPECIAL PERMISSION FOR THIS DATA- ALREADY PRUNED FOR ECUADOR ONLY)
oil <- read_excel('~/Library/CloudStorage/GoogleDrive-yl0518@princeton.edu/.shortcut-targets-by-id/1LJjchwY-EphSn7Gv624wVwwe8gEtAHn_/Ecuador Oil Pollution Project/Raw Data/Oil Data/Oil_data_Ecuador_GEM.xlsx')

plot(Latitude ~ Longitude, oil)


#Subset Ecuador
oil_aoi <- oil[oil$Country == 'Ecuador',]

summary(oil_aoi)

oil_aoi<-oil_aoi[!(is.na(oil_aoi$Latitude)),]

summary(oil_aoi)

#rename
oil<-oil_aoi

# transform the dataframe to a "track" which is preferred by amt package 
#EPSG:4326, also known as the WGS84 projection is a coordinate system used in Google Earth and GSP systems
oil_amt <- make_track(oil, .x = Longitude, .y = Latitude, crs = 'EPSG:4326' )

summary(oil_amt)

# project data from WGS1984 to UTM Zone 17S (Ecuador-based mapping)
oil_amt_utm <- transform_coords(oil_amt, crs_from = 4326, crs_to = 32717)

# make a kerenel density map from the oil lat and longs
oil_kde <- hr_kde(oil_amt_utm)

#plot the "utilization distribution"
plot(oil_kde$ud)


# read in samples csv
samples <- read.csv("Full_Metal_results_RAnalysis.csv")
head(samples)

#subset by Ecuador samples
samples<-samples[samples$`Country.code`=="Ecuador",]
  
# clean
# samples <- samples[c(1:71),]
# summary(samples)

#indicate lat and long
samples$Lat<-as.numeric(samples$Lat)
samples$Lon<-as.numeric(samples$Lon)

colnames(samples)
# make spatial data in coordiinate system used before
samples_sp <- SpatialPointsDataFrame(samples[,c(39,38)],
                                     samples,
                                     proj4string = CRS('EPSG:4326'))
plot(samples_sp)

# transform to UTM Zone 17S Ecador
samples_utm <- spTransform(samples_sp, CRS('EPSG:32717'))
plot(samples_utm)


# transform kde to raster
oil_rast <- raster(oil_kde$ud)

# extract value of the kernel denisty estimation for the sample points
samples_utm$kde_values <- raster::extract(oil_rast, samples_utm[,c(17,16)])
summary(samples_utm)
samples_utm@data

listofkde<-as.data.frame(samples_utm$kde_values)
listofkde

plot(oil_rast)
points(samples_utm)

#add Ecuador outline


#read in ecuador map
ecmap<-readOGR("~/Library/CloudStorage/GoogleDrive-yl0518@princeton.edu/.shortcut-targets-by-id/1LJjchwY-EphSn7Gv624wVwwe8gEtAHn_/Ecuador Oil Pollution Project/Code & Results/EC_Shapefile/sb978rn7613.shp")

#what system is it in
proj4string(ecmap)

#transform to CRS zone 17 used before
ec_projected<- spTransform(ecmap, CRS('EPSG:32717'))

#now in zone 17 UTM ecuador projected
proj4string(ec_projected)

#plot heatmap then ecuador outline then points
plot(oil_rast)
plot(ec_projected, add=TRUE)
plot(samples_utm, add=TRUE)

#save
pdf("~/Library/CloudStorage/GoogleDrive-yl0518@princeton.edu/.shortcut-targets-by-id/1LJjchwY-EphSn7Gv624wVwwe8gEtAHn_/Ecuador Oil Pollution Project/Code & Results/Figures/Rasterizedplot.pdf")
plot(oil_rast )
plot(ec_projected, add=TRUE)
plot(samples_utm, add=TRUE)
dev.off()

#instead of plotting in r, write out the files and load into acrmap
#writeRaster(oil_rast, "oil_raster.tif")
#shapefile(samples_utm, "sampledata.shp")
#shapefile(ec_projected, "ecmapdata.shp")



#####PLOTTING KDE MAP#####
#plotting in r using rastervis

pdf("~/Library/CloudStorage/GoogleDrive-yl0518@princeton.edu/.shortcut-targets-by-id/1LJjchwY-EphSn7Gv624wVwwe8gEtAHn_/Ecuador Oil Pollution Project/Code & Results/Figures/KDPLOT.pdf")

masked_oil<- mask(oil_rast, ec_projected)
p <- levelplot(masked_oil, layers = 1, margin=T, contour=TRUE, colorkey = list(title = "kde") , main = "Kernel Density Plot of Oil Fields in Ecuador",  xlab="Longitude (UTM 17N)", ylab="Latitude (UTM 17N)")
p + latticeExtra::layer(sp.lines(ec_projected, col="white", lwd=0.5)) + latticeExtra::layer(sp.lines(samples_utm, col="white", lwd=0.5)) 


dev.off()

##calculating distance from kernel center
test<-as.data.frame(oil_kde$data)
test2<-as.data.frame(oil_kde$ud)
test3<-cbind(test,test2)
test4<-test3[which.max(test3$lyr.1),]

#largest kernel ud value
test4

#x and y location for that max ud value
x<-as.numeric(960578.9)
y<-as.numeric(9964177)

#make spatial point for that coordinate 
largestkde<-SpatialPoints(coords=test4[,1:2],CRS('EPSG:32717'))
largestkde

#gdistance from each site to kernel peak
distances<-gDistance(largestkde,samples_utm, byid=TRUE)

#get the distances for each location, then do an anova across all sites by full heavy metal load (boxplot) thena. linear model of averaged heavy metal load within a pop and the KDE distance values as a scatterplot. 
samples_utm$distances<-distances

samples_df<-samples_utm@data

colnames(samples_df)

#look at city distances

citydistances<-samples_df %>% distinct(`City.Town`, .keep_all=TRUE)
citydistances<-cbind(citydistances$`City.Town`, citydistances$distances)
citydistances

# "Sangay National Park Macas" "259783.704607078" meters
# "Miazal"                     "276488.84296328"  meters
# "Tiputini"                   "72177.4837829124" meters
# "Hollin River"               "104363.594230283" meters
# "San Rafael Falls"           "83560.0101596611" meters



####KDE AND METAL ANALYSIS####

#rename since we will be adding more values for each individual
all_data<-samples_df

#calculate full load per inividual
all_data$Full_load<-rowSums(all_data[,6:28])

#rename long names
all_data[all_data=="Sangay National Park Macas"] <- "Sangay"
all_data[all_data=="San Rafael Falls"] <- "San Rafael"
all_data[all_data=="Hollin River"] <- "Hollin"


###Analysis###

# mean full conamination load by site...
# #summarize mean"full contamination load" by each site"
# all_data$Full_load
# 
# tiputini_df<- all_data[ which(all_data$City.Town=='Tiputini'), ]
# tiputini_df
# 
# hollin_df<-all_data[ which(all_data$City.Town=='Hollin'), ]
# sanrafael_df<-all_data[ which(all_data$City.Town=='San Rafael'), ]
# sangay_df<-all_data[ which(all_data$City.Town=='Sangay'), ]
# miazal_df<-all_data[ which(all_data$City.Town=='Miazal'), ]
# 
# #calculate mean full loads...
# tiputini_mean<-mean(tiputini_df$Full_load)
# hollin_mean<-mean(hollin_df$Full_load)
# sanrafael_mean<-mean(sanrafael_df$Full_load)
# sangay_mean<-mean(sangay_df$Full_load)
# miazal_mean<-mean(miazal_df$Full_load)
# 
# #new dataframe of just the cities, kde values, and full load
# Sites<-c("Tiputini","Hollin","San Rafael","Sangay","Miazal")
# KDE<-c(72177.48, 104363.59,83560.01,259783.70,276488.84)
# MeanHMLoad<-c(tiputini_mean, hollin_mean, sanrafael_mean, sangay_mean, miazal_mean)
# 
# mean_df<-cbind(Sites, KDE, MeanHMLoad)
# mean_df<-as.data.frame(mean_df)
# 
# mean_df$MeanHMLoad<-as.numeric(mean_df$MeanHMLoad)
# mean_df$KDE<-as.numeric(mean_df$KDE)
# 
# ##Pearsons correlation coefficient of how KDE and full load are related- p=0.5606, cor=-0.352
# cortest <- cor.test(mean_df$MeanHMLoad, mean_df$KDE,
#                     method = "pearson")
# cortest
# 
# #linear model... p=0.5606
# linearmodel<- lm(MeanHMLoad ~ KDE, data = mean_df)
# summary(linearmodel)


##Linear mixed effects model of how KDE affects full load by site
all_data
#lm of full indiv load sig to KDE (p=0.0254)
linearmodel<- lm(Full_load ~ kde_values, data = all_data)
summary(linearmodel)
linearmodel
#Residual standard error: 179.3 on 71 degrees of freedom
#Multiple R-squared:  0.06842,	Adjusted R-squared:  0.0553 
#F-statistic: 5.214 on 1 and 71 DF,  p-value: 0.02539


##Pearsons correlation coefficient of how KDE and full load are related
cortest <- cor.test(all_data$Full_load, all_data$kde_values, 
                    method = "pearson")
cortest

#t = 2.296, df = 71, p-value = 0.02463
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:  0.0349494 0.4648532
#sample estimates: cor= 0.2629028 


#plot#
all_data$City.Town <- factor(all_data$City.Town , levels=c("Tiputini","Hollin","San Rafael","Sangay","Miazal"))
unique(all_data$City.Town)


#scatterplot
p <- ggplot(all_data, aes(x=kde_values, y=Full_load, color = City.Town)) +
  geom_point() + 
 # scale_color_manual(values=c("red3","violetred3","orchid4","darkorchid4", "black"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+ 
  labs(title = "LMM Model of Full Contamination Load by KDE Values per Site" ,  subtitle = "r2=0.0553, p-value=0.025", x="KDE Value", y="Full HM Load")

#manually add the lmm regression line
p + xlim(0,2.352535e-11) +geom_abline(intercept = 3.969e+02, slope = 0.0553,, color="grey", linetype="dashed") 


#####ALL by POP / KDE####
all_data$Full_load
all_data$City.Town
all_data$kde_values
all_data$distances

#perform Bartlett's test for unequal variances...shows sig different variances!
bartlett.test(all_data$Full_load ~ all_data$City.Town)

#welch's anova for unequal variances...
All_aov <- oneway.test(Full_load ~ City.Town, data=all_data, var.equal = FALSE)
All_aov

#posthoc test!
#Games-Howell multiple comparisons method.The Games-Howell post hoc test, like Welch’s analysis of variance, does not require the groups to have equal standard deviations.
all_data %>% games_howell_test(Full_load ~ City.Town)

#regular anova...not meaninful here...
All_aov <- aov(Full_load ~ City.Town, data=all_data)
All_aov

# Pairwise comparisons
All_post<-all_data %>% games_howell_test(Full_load ~ City.Town)
All_post

All_post <-  All_post %>% mutate(y.position = "y.position")

All_plot<-ggboxplot(all_data, x = "City.Town", y = "Full_load", color="kde_values",
          xlab = "Site",ylab = "Mean Heavy Metal Load (ppm)") + stat_pvalue_manual(All_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+ labs(title = "Mean Heavy Metal Load by Site (KDE values)" ,subtitle = "AOV, p=0.0254",
       caption = get_pwc_label(All_post))+theme(legend.position="right") + gradient_color(c("purple4", "red3"))

All_plot
#no sig differences in full load across pops



#####DO THIS FOR ALL METALS##
##do this for each metals to see if they are related to kde values
metals<-colnames(all_data[6:28])
metals


##V sig (tip-hol, tip-sang, tip-mia)
#welch's anova for unequal variances...
V_aov <- oneway.test(V ~ City.Town, data=all_data, var.equal = FALSE)
V_aov
V_post<-all_data %>% games_howell_test(V ~ City.Town)
V_post

max(all_data$V)

V_post <- V_post  %>% mutate(y.position = c(0.6,0.62,0.64,0.68,0.70,0.72,0.74,0.76, 0.78, 0.80))
V_plot<-ggboxplot(all_data, x = "City.Town", y = "V", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "V (ppm)") + stat_pvalue_manual(V_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(V_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Co sig (tip-sang, tip-mia)
#welch's anova for unequal variances...
Co_aov <- oneway.test(Co ~ City.Town, data=all_data, var.equal = FALSE)
Co_aov
Co_post<-all_data %>% games_howell_test(Co ~ City.Town)
Co_post

max(all_data$Co)

Co_post <- Co_post  %>% mutate(y.position = c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9))
Co_plot<-ggboxplot(all_data, x = "City.Town", y = "Co", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Co (ppm)") + stat_pvalue_manual(Co_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Co_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Al not sig
#welch's anoAla for unequal Alariances...
Al_aov <- oneway.test(Al ~ City.Town, data=all_data, var.equal = FALSE)
Al_post<-all_data %>% games_howell_test(Al ~ City.Town)
Al_post

max(all_data$Al)

Al_post <- Al_post  %>% mutate(y.position = c(140))
Al_plot<-ggboxplot(all_data, x = "City.Town", y = "Al", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Al (ppm)") + stat_pvalue_manual(Al_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Al_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Cr not sig
#welch's anova for unequal variances...
Cr_aov <- oneway.test(Cr ~ City.Town, data=all_data, var.equal = FALSE)
Cr_post<-all_data %>% games_howell_test(Cr ~ City.Town)
Cr_post

max(all_data$Cr)

Cr_post <- Cr_post  %>% mutate(y.position = c(3))
Cr_plot<-ggboxplot(all_data, x = "City.Town", y = "Cr", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Cr (ppm)") + stat_pvalue_manual(Cr_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Cr_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Mn tip-sanr, tip-sang
#welch's anova for unequal variances...
Mn_aov <- oneway.test(Mn ~ City.Town, data=all_data, var.equal = FALSE)
Mn_post<-all_data %>% games_howell_test(Mn ~ City.Town)
Mn_post

max(all_data$Mn)

Mn_post <- Mn_post  %>% mutate(y.position = c(610,620,630,640,650,660,670,680,690,700))
Mn_plot<-ggboxplot(all_data, x = "City.Town", y = "Mn", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Mn (ppm)") + stat_pvalue_manual(Mn_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Mn_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Fe not sig
#welch's anova for unequal variances...
Fe_aov <- oneway.test(Fe ~ City.Town, data=all_data, var.equal = FALSE)
Fe_post<-all_data %>% games_howell_test(Fe ~ City.Town)
Fe_post

max(all_data$Fe)

Fe_post <- Fe_post  %>% mutate(y.position = c(200))
Fe_plot<-ggboxplot(all_data, x = "City.Town", y = "Fe", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Fe (ppm)") + stat_pvalue_manual(Fe_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Fe_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Ni sig (tip-hol, tip-sanr, tip-sang)
#welch's anova for unequal variances...
Ni_aov <- oneway.test(Ni ~ City.Town, data=all_data, var.equal = FALSE)
Ni_post<-all_data %>% games_howell_test(Ni ~ City.Town)
Ni_post

max(all_data$Ni)

Ni_post <- Ni_post  %>% mutate(y.position = c(6.5,7,7.5,8,8.5,9,9.5,10,10.5,11))
Ni_plot<-ggboxplot(all_data, x = "City.Town", y = "Ni", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Ni (ppm)") + stat_pvalue_manual(Ni_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Ni_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Cu sig (tip-hol)
#welch's anova for unequal variances...
Cu_aov <- oneway.test(Cu ~ City.Town, data=all_data, var.equal = FALSE)
Cu_post<-all_data %>% games_howell_test(Cu ~ City.Town)
Cu_post

max(all_data$Cu)

Cu_post <- Cu_post  %>% mutate(y.position = c(40))
Cu_plot<-ggboxplot(all_data, x = "City.Town", y = "Cu", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Cu (ppm)") + stat_pvalue_manual(Cu_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Cu_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Zn sig (tip-hol, tip-sanr)
#welch's anova for unequal variances...
Zn_aov <- oneway.test(Zn ~ City.Town, data=all_data, var.equal = FALSE)
Zn_post<-all_data %>% games_howell_test(Zn ~ City.Town)
Zn_post

max(all_data$Zn)

Zn_post <- Zn_post  %>% mutate(y.position = c(310,320,330,340,350,360,370,380,390,400))
Zn_plot<-ggboxplot(all_data, x = "City.Town", y = "Zn", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Zn (ppm)") + stat_pvalue_manual(Zn_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Zn_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


# ##As= NA
# #welch's anova for unequal variances...
# As_aov <- oneway.test(As ~ City.Town, data=all_data, var.equal = FALSE)
# As_post<-all_data %>% games_howell_test(As ~ City.Town)
# As_post
# 
# As_post <- As_post  %>% mutate(y.position = c(0.5,0.52,0.54,0.58,0.60,0.62,0.64,0.66, 0.68, 0.70))
# As_plot<-ggboxplot(all_data, x = "City.Town", y = "As", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
#           xlab = "KDE Values",ylab = "As (ppm)") + stat_pvalue_manual(As_post, label = "p.adj.signif", hide.ns = TRUE) +
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   labs(title = "Full HM load by site (KDE values)" ,
#        caption = get_pwc_label(As_post)
#   )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Se sig (tip-hol)
#welch's anova for unequal variances...
Se_aov <- oneway.test(Se ~ City.Town, data=all_data, var.equal = FALSE)
Se_post<-all_data %>% games_howell_test(Se ~ City.Town)
Se_post

max(all_data$Se)

Se_post <- Se_post  %>% mutate(y.position = c(1.5))
Se_plot<-ggboxplot(all_data, x = "City.Town", y = "Se", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Se (ppm)") + stat_pvalue_manual(Se_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Se_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Sr sig (tip-hol,hol-mia)
#welch's anova for unequal variances...
Sr_aov <- oneway.test(Sr ~ City.Town, data=all_data, var.equal = FALSE)
Sr_post<-all_data %>% games_howell_test(Sr ~ City.Town)
Sr_post

max(all_data$Sr)

Sr_post <- Sr_post  %>% mutate(y.position = c(75,77,79,81,83,85,87,89,91,93))
Sr_plot<-ggboxplot(all_data, x = "City.Town", y = "Sr", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Sr (ppm)") + stat_pvalue_manual(Sr_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Sr_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Mo not sig
#welch's anova for unequal variances...
Mo_aov <- oneway.test(Mo ~ City.Town, data=all_data, var.equal = FALSE)
Mo_post<-all_data %>% games_howell_test(Mo ~ City.Town)
Mo_post

max(all_data$Mo)

Mo_post <- Mo_post  %>% mutate(y.position = c(0.5,0.52,0.54,0.58,0.60,0.62,0.64,0.66, 0.68, 0.70))
Mo_plot<-ggboxplot(all_data, x = "City.Town", y = "Mo", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Mo (ppm)") + stat_pvalue_manual(Mo_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Mo_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Ag not sig
#welch's anova for unequal variances...
Ag_aov <- oneway.test(Ag ~ City.Town, data=all_data, var.equal = FALSE)
Ag_post<-all_data %>% games_howell_test(Ag ~ City.Town)
Ag_post

max(all_data$Ag)

Ag_post <- Ag_post  %>% mutate(y.position = c(0.5,0.52,0.54,0.58,0.60,0.62,0.64,0.66, 0.68, 0.70))
Ag_plot<-ggboxplot(all_data, x = "City.Town", y = "Ag", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Ag (ppm)") + stat_pvalue_manual(Ag_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Ag_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Cd sig (tip-hol, tip-sanr, tip-sang)
#welch's anova for unequal variances...
Cd_aov <- oneway.test(Cd ~ City.Town, data=all_data, var.equal = FALSE)
Cd_post<-all_data %>% games_howell_test(Cd ~ City.Town)
Cd_post

max(all_data$Cd)

Cd_post <- Cd_post  %>% mutate(y.position = c(5,5.5,6,6.5,7,7.5,8,8.5,9,9.5))
Cd_plot<-ggboxplot(all_data, x = "City.Town", y = "Cd", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Cd (ppm)") + stat_pvalue_manual(Cd_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Cd_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Sn sig (hol-mia)
#welch's anova for unequal variances...
Sn_aov <- oneway.test(Sn ~ City.Town, data=all_data, var.equal = FALSE)
Sn_post<-all_data %>% games_howell_test(Sn ~ City.Town)
Sn_post

max(all_data$Sn)

Sn_post <- Sn_post  %>% mutate(y.position = c(5,5.5,6,6.5,7,7.5,8,8.5,9,9.5))
Sn_plot<-ggboxplot(all_data, x = "City.Town", y = "Sn", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Sn (ppm)") + stat_pvalue_manual(Sn_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Sn_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Sb sig (tip-hol, tip-sanr, tip-sang, hol-mia, sanr-mia, sang-mia)
#welch's anova for unequal variances...
Sb_aov <- oneway.test(Sb ~ City.Town, data=all_data, var.equal = FALSE)
Sb_post<-all_data %>% games_howell_test(Sb ~ City.Town)
Sb_post

max(all_data$Sb)

Sb_post <- Sb_post  %>% mutate(y.position = c(0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7))
Sb_plot<-ggboxplot(all_data, x = "City.Town", y = "Sb", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Sb (ppm)") + stat_pvalue_manual(Sb_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Sb_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Cs sig (tip-hol, tip-sanr, tip-sang, hol-mia,sanr-mia, sang-mia)
#welch's anova for unequal variances...
Cs_aov <- oneway.test(Cs ~ City.Town, data=all_data, var.equal = FALSE)
Cs_post<-all_data %>% games_howell_test(Cs ~ City.Town)
Cs_post

max(all_data$Cs)

Cs_post <- Cs_post  %>% mutate(y.position = c(0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19))
Cs_plot<-ggboxplot(all_data, x = "City.Town", y = "Cs", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Cs (ppm)") + stat_pvalue_manual(Cs_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Cs_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Ba sig (tip-sanrl, tip-sang, sang-hol)
#welch's anova for unequal variances...
Ba_aov <- oneway.test(Ba ~ City.Town, data=all_data, var.equal = FALSE)
Ba_post<-all_data %>% games_howell_test(Ba ~ City.Town)
Ba_post

max(all_data$Ba)

Ba_post <- Ba_post  %>% mutate(y.position = c(200,205,210,215,220,225,230,235,240,245))
Ba_plot<-ggboxplot(all_data, x = "City.Town", y = "Ba", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Ba (ppm)") + stat_pvalue_manual(Ba_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Ba_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Hg sig (hol-mia, sang-mia)
#welch's anova for unequal variances...
Hg_aov <- oneway.test(Hg ~ City.Town, data=all_data, var.equal = FALSE)
Hg_post<-all_data %>% games_howell_test(Hg ~ City.Town)
Hg_post

max(all_data$Hg)

Hg_post <- Hg_post  %>% mutate(y.position = c(0.5,0.52,0.54,0.58,0.60,0.62,0.64,0.66, 0.68, 0.70))
Hg_plot<-ggboxplot(all_data, x = "City.Town", y = "Hg", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Hg (ppm)") + stat_pvalue_manual(Hg_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Hg_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Tl sig ((tip-hol), tip-sang, sang-mia)
#welch's anova for unequal variances...
Tl_aov <- oneway.test(Tl ~ City.Town, data=all_data, var.equal = FALSE)
Tl_post<-all_data %>% games_howell_test(Tl ~ City.Town)
Tl_post

max(all_data$Tl)

Tl_post <- Tl_post  %>% mutate(y.position = c(0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81))
Tl_plot<-ggboxplot(all_data, x = "City.Town", y = "Tl", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Tl (ppm)") + stat_pvalue_manual(Tl_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Tl_post)
  )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


##Pb not sig
#welch's anova for unequal variances...
Pb_aov <- oneway.test(Pb ~ City.Town, data=all_data, var.equal = FALSE)
Pb_post<-all_data %>% games_howell_test(Pb ~ City.Town)
Pb_post

max(all_data$Pb)

Pb_post <- Pb_post  %>% mutate(y.position = c(10))
Pb_plot<-ggboxplot(all_data, x = "City.Town", y = "Pb", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
          xlab = "KDE Values",ylab = "Pb (ppm)") + stat_pvalue_manual(Pb_post, label = "p.adj.signif", hide.ns = TRUE) +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  labs(title = "Full HM load by site (KDE values)" ,
       caption = get_pwc_label(Pb_post)
  ) +theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


# ##U sig NA
# #welch's anova for unequal variances...
# U_aov <- oneway.test(U ~ City.Town, data=all_data, var.equal = FALSE)
# U_post<-all_data %>% games_howell_test(U ~ City.Town)
# U_post
# 
# 
# U_post <- U_post  %>% mutate(y.position = c(0.5,0.52,0.54,0.58,0.60,0.62,0.64,0.66, 0.68, 0.70))
# U_plot<-ggboxplot(all_data, x = "City.Town", y = "U", color="kde_values",title = "T test of full contamination load by KDE values per site" ,
#           xlab = "KDE Values",ylab = "U (ppm)") + stat_pvalue_manual(U_post, label = "p.adj.signif", hide.ns = TRUE) +
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   labs(title = "Full HM load by site (KDE values)" ,
#        caption = get_pwc_label(U_post)
#   )+theme(legend.position="right")+ gradient_color(c("purple4", "red3"))


###ALL PLOTS

All_plot

Al_plot
V_plot
Cr_plot
Mn_plot
Fe_plot
Co_plot
Ni_plot
Cu_plot
Zn_plot
As_plot # too many NA's
Se_plot
Sr_plot
Mo_plot
Ag_plot
Cd_plot
Sn_plot
Sb_plot
Cs_plot
Ba_plot
Hg_plot
Tl_plot
Pb_plot
U_plot # too many NA's


# #quickplots
# metals
# 
# ggboxplot(all_data, x = "City.Town", y = "Al", color="City.Town") #kde trend
# ggboxplot(all_data, x = "City.Town", y = "V", color="City.Town") #san rafael is high then according to kde
# ggboxplot(all_data, x = "City.Town", y = "Cr", color="City.Town") #opposite of kde?
# ggboxplot(all_data, x = "City.Town", y = "Mn", color="City.Town") #tip>miazal>hollin>sangay>sanrafael
# ggboxplot(all_data, x = "City.Town", y = "Fe", color="City.Town") #same as mn
# ggboxplot(all_data, x = "City.Town", y = "Co", color="City.Town") #tiputini is super high then the rest lowww
# ggboxplot(all_data, x = "City.Town", y = "Ni", color="City.Town") #tip, miazzal, sangay, sanrafael, hollin
# ggboxplot(all_data, x = "City.Town", y = "Cu", color="City.Town") #tiputini and miazal, then the rest
# ggboxplot(all_data, x = "City.Town", y = "Zn", color="City.Town") #sanrafael, then hollin and sngay then tip and miz
# ggboxplot(all_data, x = "City.Town", y = "As", color="City.Town") #high san rafael then hollin and sangay and the rest
# ggboxplot(all_data, x = "City.Town", y = "Se", color="City.Town") #rafael then hollin and miazal sangau then tiputini
# ggboxplot(all_data, x = "City.Town", y = "Sr", color="City.Town") #tip and hollin then the rest down kde
# ggboxplot(all_data, x = "City.Town", y = "Mo", color="City.Town") #san rafael then hollin, tip, sanga, miaz
# ggboxplot(all_data, x = "City.Town", y = "Ag", color="City.Town") #sang, rafael and miazal, then tip and hollin
# ggboxplot(all_data, x = "City.Town", y = "Cd", color="City.Town") #tip and miazal, then the rest
# ggboxplot(all_data, x = "City.Town", y = "Sn", color="City.Town") #all about the same
# ggboxplot(all_data, x = "City.Town", y = "Sb", color="City.Town") #hollin san rafael an sangau high and tip and miaz low
# ggboxplot(all_data, x = "City.Town", y = "Cs", color="City.Town") #hollin sanaf and sang high then tip an dmiaz
# ggboxplot(all_data, x = "City.Town", y = "Ba", color="City.Town") #hollin, tip, then miazal and rest
# ggboxplot(all_data, x = "City.Town", y = "Hg", color="City.Town") #hollin,sanraf, sangay, miaz
# ggboxplot(all_data, x = "City.Town", y = "Tl", color="City.Town") #all low
# ggboxplot(all_data, x = "City.Town", y = "Pb", color="City.Town") #all low
# ggboxplot(all_data, x = "City.Town", y = "U", color="City.Town") #tiputini high then rest




#####EACH METAL BY RISK GROUPS (high medium low)#####

##Add Risk category by KDE
#add a new group by Kernel density distances (hollin and sanrafael as "medium risk" and miazal and sangay as "low risk")
City.Town<-c("Tiputini", "Hollin", "San Rafael","Sangay", "Miazal")
Risks<-c("High","Medium","Medium","Low","Low")
df<-data.frame(City.Town,Risks)
full_data<-merge(all_data, df, by="City.Town")
##THIS IS ONLY FOR ECUADOR

#All weighted by featherweight 
All_kruskal_weighted <- weighted_mannwhitney(all_data, Full_load,Risks, Feather.Weight)
All_kruskal_weighted

All_kruskal <- all_data %>% kruskal_test(Full_load ~ Risks)
All_kruskal

# Pairwise comparisons
All_wilcox <- all_data %>% 
  wilcox_test(Full_load ~ Risks, p.adjust.method = "bonferroni") 
All_wilcox

ggboxplot(all_data, x = "Risks", y = "Full_load")

pdf("Figures/All_Kruskalgroup.pdf")

All_wilcox <-  All_wilcox %>% mutate(y.position = c(4000, 4000, 4000))

ggboxplot(all_data, x = "Risks", y = "Full_load",
          ylab = "Full Heavy Metal Load", xlab = "Risks") + stat_pvalue_manual(All_wilcox, label = "p.adj", hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = All_kruskal_weighted$p.value,
       caption = get_pwc_label(All_wilcox)
       
  )
dev.off()



###ELEMENTS BY RISK GROUPS###

##TEST SIG > CONTROL##

##Co##
EC_Co_kruskal_weighted <- weighted_mannwhitney(EC_data, Co,Risks, Feather.Weight)
EC_Co_kruskal_weighted

EC_Co_kruskal <- EC_data %>% kruskal_test(Co ~ Risks)
EC_Co_kruskal

# Pairwise comparisons
EC_Co_wilcox <- EC_data %>% 
  wilcox_test(Co ~ Risks, p.adjust.method = "bonferroni") 
EC_Co_wilcox

pdf("Figures/EC_Co_KruskalRisks.pdf")
EC_Co_wilcox <- EC_Co_wilcox %>% mutate(y.position = 1,2,3)
ggboxplot(EC_data, x = "Risks", y = "Co",
          ylab = "Co Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Co_wilcox, hide.ns = TRUE,  y.position = "y.position") +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Co_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Co_wilcox)
  )
dev.off()

EC_Co_wilcox <- EC_Co_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Co",
          ylab = "Co Concentration", xlab = "Risks") +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Co_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Co_wilcox)
  )


##Fe##
EC_Fe_kruskal_weighted <- weighted_mannwhitney(EC_data, Fe,Risks, Feather.Weight)
EC_Fe_kruskal_weighted

EC_Fe_kruskal <- EC_data %>% kruskal_test(Fe ~ Risks)
EC_Fe_kruskal

# Pairwise comparisons
EC_Fe_wilcox <- EC_data %>% 
  wilcox_test(Fe ~ Risks, p.adjust.method = "bonferroni") 
EC_Fe_wilcox

pdf("Figures/EC_Fe_KruskalRisks.pdf")
EC_Fe_wilcox <- EC_Fe_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Fe",
          ylab = "Fe Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Fe_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Fe_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Fe_wilcox)
  )
dev.off()

EC_Fe_wilcox <- EC_Fe_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Fe",
          ylab = "Fe Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Fe_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Fe_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Fe_wilcox)
  )


##Ni##
EC_Ni_kruskal_weighted <- weighted_mannwhitney(EC_data, Ni,Risks, Feather.Weight)
EC_Ni_kruskal_weighted

EC_Ni_kruskal <- EC_data %>% kruskal_test(Ni ~ Risks)
EC_Ni_kruskal

# Pairwise comparisons
EC_Ni_wilcox <- EC_data %>% 
  wilcox_test(Ni ~ Risks, p.adjust.method = "bonferroni") 
EC_Ni_wilcox

pdf("Figures/EC_Ni_KruskalRisks.pdf")
EC_Ni_wilcox <- EC_Ni_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Ni",
          ylab = "Ni Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Ni_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Ni_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Ni_wilcox)
  )
dev.off()

EC_Ni_wilcox <- EC_Ni_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Ni",
          ylab = "Ni Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Ni_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Ni_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Ni_wilcox)
  )

##V##
EC_V_kruskal_weighted <- weighted_mannwhitney(EC_data, V,Risks, Feather.Weight)
EC_V_kruskal_weighted

EC_V_kruskal <- EC_data %>% kruskal_test(V ~ Risks)
EC_V_kruskal

# Pairwise comparisons
EC_V_wilcox <- EC_data %>% 
  wilcox_test(V ~ Risks, p.adjust.method = "bonferroni") 
EC_V_wilcox

pdf("Figures/EC_V_KruskalRisks.pdf")
EC_V_wilcox <- EC_V_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "V",
          ylab = "V Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_V_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_V_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_V_wilcox)
  )
dev.off()

EC_V_wilcox <- EC_V_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "V",
          ylab = "V Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_V_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_V_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_V_wilcox)
  )


##Cu##
EC_Cu_kruskal_weighted <- weighted_mannwhitney(EC_data, Cu,Risks, Feather.Weight)
EC_Cu_kruskal_weighted

EC_Cu_kruskal <- EC_data %>% kruskal_test(Cu ~ Risks)
EC_Cu_kruskal

# Pairwise comparisons
EC_Cu_wilcox <- EC_data %>% 
  wilcox_test(Cu ~ Risks, p.adjust.method = "bonferroni") 
EC_Cu_wilcox

pdf("Figures/EC_Cu_KruskalRisks.pdf")
EC_Cu_wilcox <- EC_Cu_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cu",
          ylab = "Cu Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cu_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cu_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cu_wilcox)
  )
dev.off()

EC_Cu_wilcox <- EC_Cu_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cu",
          ylab = "Cu Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cu_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cu_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cu_wilcox)
  )



##Ba##
EC_Ba_kruskal_weighted <- weighted_mannwhitney(EC_data, Ba,Risks, Feather.Weight)
EC_Ba_kruskal_weighted

EC_Ba_kruskal <- EC_data %>% kruskal_test(Ba ~ Risks)
EC_Ba_kruskal

# Pairwise comparisons
EC_Ba_wilcox <- EC_data %>% 
  wilcox_test(Ba ~ Risks, p.adjust.method = "bonferroni") 
EC_Ba_wilcox

pdf("Figures/EC_Ba_KruskalRisks.pdf")
EC_Ba_wilcox <- EC_Ba_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Ba",
          ylab = "Ba Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Ba_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Ba_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Ba_wilcox)
  )
dev.off()

EC_Ba_wilcox <- EC_Ba_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Ba",
          ylab = "Ba Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Ba_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Ba_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Ba_wilcox)
  )


##Cd##
EC_Cd_kruskal_weighted <- weighted_mannwhitney(EC_data, Cd,Risks, Feather.Weight)
EC_Cd_kruskal_weighted

EC_Cd_kruskal <- EC_data %>% kruskal_test(Cd ~ Risks)
EC_Cd_kruskal

# Pairwise comparisons
EC_Cd_wilcox <- EC_data %>% 
  wilcox_test(Cd ~ Risks, p.adjust.method = "bonferroni") 
EC_Cd_wilcox

pdf("Figures/EC_Cd_KruskalRisks.pdf")
EC_Cd_wilcox <- EC_Cd_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cd",
          ylab = "Cd Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cd_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cd_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cd_wilcox)
  )
dev.off()

EC_Cd_wilcox <- EC_Cd_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cd",
          ylab = "Cd Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cd_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cd_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cd_wilcox)
  )


##Mn##
EC_Mn_kruskal_weighted <- weighted_mannwhitney(EC_data, Mn,Risks, Feather.Weight)
EC_Mn_kruskal_weighted

EC_Mn_kruskal <- EC_data %>% kruskal_test(Mn ~ Risks)
EC_Mn_kruskal

# Pairwise comparisons
EC_Mn_wilcox <- EC_data %>% 
  wilcox_test(Mn ~ Risks, p.adjust.method = "bonferroni") 
EC_Mn_wilcox

pdf("Figures/EC_Mn_KruskalRisks.pdf")
EC_Mn_wilcox <- EC_Mn_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Mn",
          ylab = "Mn Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Mn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Mn_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Mn_wilcox)
  )
dev.off()

EC_Mn_wilcox <- EC_Mn_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Mn",
          ylab = "Mn Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Mn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Mn_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Mn_wilcox)
  )



##U##
EC_U_kruskal_weighted <- weighted_mannwhitney(EC_data, U, Risks, Feather.Weight)
EC_U_kruskal_weighted

EC_U_kruskal <- EC_data %>% kruskal_test(U ~ Risks)
EC_U_kruskal

# Pairwise comparisons
EC_U_wilcox <- EC_data %>% 
  wilcox_test(U ~ Risks, p.adjust.method = "bonferroni") 
EC_U_wilcox

pdf("Figures/EC_U_KruskalRisks.pdf")
EC_U_wilcox <- EC_U_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "U",
          ylab = "U Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_U_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_U_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_U_wilcox)
  )
dev.off()

EC_U_wilcox <- EC_U_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "U",
          ylab = "U Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_U_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_U_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_U_wilcox)
  )



##CONTROL SIG > TEST##

##Zn##
EC_Zn_kruskal_weighted <- weighted_mannwhitney(EC_data, Zn,Risks, Feather.Weight)
EC_Zn_kruskal_weighted

EC_Zn_kruskal <- EC_data %>% kruskal_test(Zn ~ Risks)
EC_Zn_kruskal

# Pairwise comparisons
EC_Zn_wilcox <- EC_data %>% 
  wilcox_test(Zn ~ Risks, p.adjust.method = "bonferroni") 
EC_Zn_wilcox

pdf("Figures/EC_Zn_KruskalRisks.pdf")
EC_Zn_wilcox <- EC_Zn_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Zn",
          ylab = "Zn Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Zn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Zn_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Zn_wilcox)
  )
dev.off()

EC_Zn_wilcox <- EC_Zn_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Zn",
          ylab = "Zn Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Zn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Zn_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Zn_wilcox)
  )


##As##
EC_As_kruskal_weighted <- weighted_mannwhitney(EC_data, As,Risks, Feather.Weight)
EC_As_kruskal_weighted

EC_As_kruskal <- EC_data %>% kruskal_test(As ~ Risks)
EC_As_kruskal

# Pairwise comparisons
EC_As_wilcox <- EC_data %>% 
  wilcox_test(As ~ Risks, p.adjust.method = "bonferroni") 
EC_As_wilcox

pdf("Figures/EC_As_KruskalRisks.pdf")
EC_As_wilcox <- EC_As_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "As",
          ylab = "As Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_As_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_As_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_As_wilcox)
  )
dev.off()

EC_As_wilcox <- EC_As_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "As",
          ylab = "As Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_As_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_As_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_As_wilcox)
  )



##Se##
EC_Se_kruskal_weighted <- weighted_mannwhitney(EC_data, Se,Risks, Feather.Weight)
EC_Se_kruskal_weighted

EC_Se_kruskal <- EC_data %>% kruskal_test(Se ~ Risks)
EC_Se_kruskal

# Pairwise comparisons
EC_Se_wilcox <- EC_data %>% 
  wilcox_test(Se ~ Risks, p.adjust.method = "bonferroni") 
EC_Se_wilcox

pdf("Figures/EC_Se_KruskalRisks.pdf")
EC_Se_wilcox <- EC_Se_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Se",
          ylab = "Se Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Se_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Se_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Se_wilcox)
  )
dev.off()

EC_Se_wilcox <- EC_Se_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Se",
          ylab = "Se Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Se_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Se_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Se_wilcox)
  )


##Sb##
EC_Sb_kruskal_weighted <- weighted_mannwhitney(EC_data, Sb,Risks, Feather.Weight)
EC_Sb_kruskal_weighted

EC_Sb_kruskal <- EC_data %>% kruskal_test(Sb ~ Risks)
EC_Sb_kruskal

# Pairwise comparisons
EC_Sb_wilcox <- EC_data %>% 
  wilcox_test(Sb ~ Risks, p.adjust.method = "bonferroni") 
EC_Sb_wilcox

pdf("Figures/EC_Sb_KruskalRisks.pdf")
EC_Sb_wilcox <- EC_Sb_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Sb",
          ylab = "Sb Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Sb_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Sb_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Sb_wilcox)
  )
dev.off()

EC_Sb_wilcox <- EC_Sb_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Sb",
          ylab = "Sb Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Sb_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Sb_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Sb_wilcox)
  )



##Tl##
EC_Tl_kruskal_weighted <- weighted_mannwhitney(EC_data, Tl,Risks, Feather.Weight)
EC_Tl_kruskal_weighted

EC_Tl_kruskal <- EC_data %>% kruskal_test(Tl ~ Risks)
EC_Tl_kruskal

# Pairwise comparisons
EC_Tl_wilcox <- EC_data %>% 
  wilcox_test(Tl ~ Risks, p.adjust.method = "bonferroni") 
EC_Tl_wilcox

pdf("Figures/EC_Tl_KruskalRisks.pdf")
EC_Tl_wilcox <- EC_Tl_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Tl",
          ylab = "Tl Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Tl_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Tl_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Tl_wilcox)
  )
dev.off()

EC_Tl_wilcox <- EC_Tl_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Tl",
          ylab = "Tl Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Tl_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Tl_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Tl_wilcox)
  )


##Cs##
EC_Cs_kruskal_weighted <- weighted_mannwhitney(EC_data, Cs,Risks, Feather.Weight)
EC_Cs_kruskal_weighted

EC_Cs_kruskal <- EC_data %>% kruskal_test(Cs ~ Risks)
EC_Cs_kruskal

# Pairwise comparisons
EC_Cs_wilcox <- EC_data %>% 
  wilcox_test(Cs ~ Risks, p.adjust.method = "bonferroni") 
EC_Cs_wilcox

pdf("Figures/EC_Cs_KruskalRisks.pdf")
EC_Cs_wilcox <- EC_Cs_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cs",
          ylab = "Cs Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cs_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cs_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cs_wilcox)
  )
dev.off()

EC_Cs_wilcox <- EC_Cs_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cs",
          ylab = "Cs Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cs_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cs_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cs_wilcox)
  )


##NOT SIG##


##Pb##
EC_Pb_kruskal_weighted <- weighted_mannwhitney(EC_data, Pb,Risks, Feather.Weight)
EC_Pb_kruskal_weighted

EC_Pb_kruskal <- EC_data %>% kruskal_test(Pb ~ Risks)
EC_Pb_kruskal

# Pairwise comparisons
EC_Pb_wilcox <- EC_data %>% 
  wilcox_test(Pb ~ Risks, p.adjust.method = "bonferroni") 
EC_Pb_wilcox

pdf("Figures/EC_Pb_KruskalRisks.pdf")
EC_Pb_wilcox <- EC_Pb_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Pb",
          ylab = "Pb Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Pb_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Pb_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Pb_wilcox)
  )
dev.off()

EC_Pb_wilcox <- EC_Pb_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Pb",
          ylab = "Pb Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Pb_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Pb_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Pb_wilcox)
  )


##Al##
EC_Al_kruskal_weighted <- weighted_mannwhitney(EC_data, Al,Risks, Feather.Weight)
EC_Al_kruskal_weighted

EC_Al_kruskal <- EC_data %>% kruskal_test(Al ~ Risks)
EC_Al_kruskal

# Pairwise comparisons
EC_Al_wilcox <- EC_data %>% 
  wilcox_test(Al ~ Risks, p.adjust.method = "bonferroni") 
EC_Al_wilcox

pdf("Figures/EC_Al_KruskalRisks.pdf")
EC_Al_wilcox <- EC_Al_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Al",
          ylab = "Al Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Al_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Al_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Al_wilcox)
  )
dev.off()

EC_Al_wilcox <- EC_Al_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Al",
          ylab = "Al Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Al_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Al_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Al_wilcox)
  )


##Hg##
EC_Hg_kruskal_weighted <- weighted_mannwhitney(EC_data, Hg,Risks, Feather.Weight)
EC_Hg_kruskal_weighted

EC_Hg_kruskal <- EC_data %>% kruskal_test(Hg ~ Risks)
EC_Hg_kruskal

# Pairwise comparisons
EC_Hg_wilcox <- EC_data %>% 
  wilcox_test(Hg ~ Risks, p.adjust.method = "bonferroni") 
EC_Hg_wilcox

pdf("Figures/EC_Hg_KruskalRisks.pdf")
EC_Hg_wilcox <- EC_Hg_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Hg",
          ylab = "Hg Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Hg_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Hg_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Hg_wilcox)
  )
dev.off()

EC_Hg_wilcox <- EC_Hg_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Hg",
          ylab = "Hg Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Hg_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Hg_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Hg_wilcox)
  )



##Cr##
EC_Cr_kruskal_weighted <- weighted_mannwhitney(EC_data, Cr,Risks, Feather.Weight)
EC_Cr_kruskal_weighted

EC_Cr_kruskal <- EC_data %>% kruskal_test(Cr ~ Risks)
EC_Cr_kruskal

# Pairwise comparisons
EC_Cr_wilcox <- EC_data %>% 
  wilcox_test(Cr ~ Risks, p.adjust.method = "bonferroni") 
EC_Cr_wilcox

pdf("Figures/EC_Cr_KruskalRisks.pdf")
EC_Cr_wilcox <- EC_Cr_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cr",
          ylab = "Cr Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cr_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cr_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cr_wilcox)
  )
dev.off()

EC_Cr_wilcox <- EC_Cr_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Cr",
          ylab = "Cr Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Cr_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Cr_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Cr_wilcox)
  )



##Sr##
EC_Sr_kruskal_weighted <- weighted_mannwhitney(EC_data, Sr,Risks, Feather.Weight)
EC_Sr_kruskal_weighted

EC_Sr_kruskal <- EC_data %>% kruskal_test(Sr ~ Risks)
EC_Sr_kruskal

# Pairwise comparisons
EC_Sr_wilcox <- EC_data %>% 
  wilcox_test(Sr ~ Risks, p.adjust.method = "bonferroni") 
EC_Sr_wilcox

pdf("Figures/EC_Sr_KruskalRisks.pdf")
EC_Sr_wilcox <- EC_Sr_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Sr",
          ylab = "Sr Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Sr_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Sr_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Sr_wilcox)
  )
dev.off()

EC_Sr_wilcox <- EC_Sr_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Sr",
          ylab = "Sr Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Sr_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Sr_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Sr_wilcox)
  )


##Mo##
EC_Mo_kruskal_weighted <- weighted_mannwhitney(EC_data, Mo,Risks, Feather.Weight)
EC_Mo_kruskal_weighted

EC_Mo_kruskal <- EC_data %>% kruskal_test(Mo ~ Risks)
EC_Mo_kruskal

# Pairwise comparisons
EC_Mo_wilcox <- EC_data %>% 
  wilcox_test(Mo ~ Risks, p.adjust.method = "bonferroni") 
EC_Mo_wilcox

pdf("Figures/EC_Mo_KruskalRisks.pdf")
EC_Mo_wilcox <- EC_Mo_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Mo",
          ylab = "Mo Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Mo_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Mo_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Mo_wilcox)
  )
dev.off()

EC_Mo_wilcox <- EC_Mo_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Mo",
          ylab = "Mo Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Mo_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Mo_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Mo_wilcox)
  )



##Ag##
EC_Ag_kruskal_weighted <- weighted_mannwhitney(EC_data, Ag,Risks, Feather.Weight)
EC_Ag_kruskal_weighted

EC_Ag_kruskal <- EC_data %>% kruskal_test(Ag ~ Risks)
EC_Ag_kruskal

# Pairwise comparisons
EC_Ag_wilcox <- EC_data %>% 
  wilcox_test(Ag ~ Risks, p.adjust.method = "bonferroni") 
EC_Ag_wilcox

pdf("Figures/EC_Ag_KruskalRisks.pdf")
EC_Ag_wilcox <- EC_Ag_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Ag",
          ylab = "Ag Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Ag_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Ag_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Ag_wilcox)
  )
dev.off()

EC_Ag_wilcox <- EC_Ag_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Ag",
          ylab = "Ag Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Ag_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Ag_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Ag_wilcox)
  )


##Sn##
EC_Sn_kruskal_weighted <- weighted_mannwhitney(EC_data, Sn,Risks, Feather.Weight)
EC_Sn_kruskal_weighted

EC_Sn_kruskal <- EC_data %>% kruskal_test(Sn ~ Risks)
EC_Sn_kruskal

# Pairwise comparisons
EC_Sn_wilcox <- EC_data %>% 
  wilcox_test(Sn ~ Risks, p.adjust.method = "bonferroni") 
EC_Sn_wilcox

pdf("Figures/EC_Sn_KruskalRisks.pdf")
EC_Sn_wilcox <- EC_Sn_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Sn",
          ylab = "Sn Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Sn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Sn_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Sn_wilcox)
  )
dev.off()

EC_Sn_wilcox <- EC_Sn_wilcox %>% add_xy_position(x = "Risks")
ggboxplot(EC_data, x = "Risks", y = "Sn",
          ylab = "Sn Concentration", xlab = "Risks") +
  stat_pvalue_manual(EC_Sn_wilcox, hide.ns = TRUE) +
  labs(title = "Weighted Kruskal-Wallis" ,        subtitle = EC_Sn_kruskal_weighted$p.value,
       caption = get_pwc_label(EC_Sn_wilcox)
  )



#all wilcox tests...
EC_Zn_wilcox
EC_Mn_wilcox
EC_Fe_wilcox
EC_Al_wilcox
EC_Cu_wilcox
EC_Ba_wilcox
EC_Sr_wilcox
EC_Ni_wilcox
EC_Cr_wilcox
EC_Pb_wilcox
EC_Se_wilcox
EC_Sn_wilcox
EC_Cd_wilcox
EC_Co_wilcox
EC_Hg_wilcox
EC_Sb_wilcox
EC_V_wilcox
EC_Mo_wilcox
EC_Tl_wilcox
EC_Cs_wilcox
EC_Ag_wilcox
EC_As_wilcox
EC_U_wilcox


