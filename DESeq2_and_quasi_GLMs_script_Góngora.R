setwd("<path_to_data_folder>")
data=read.csv("Table_S1.csv", header=T)
#Alpha is the proportion of littoral carbon on the diet of individual murres
#d15N and d13C values for the littoral and pelagic prey items were obtained from Gongora et al. 2017 (Marine Pollution Bulletin)
data$alpha <- ((data$d13C-(-23.09))/((-17.00)-(-23.09)))
#Baseline is the expectation of d15N at the base of the ecosystem
data$baseline_d15N <- ((data$alpha*9.54)+((1-data$alpha)*9.18))
#tpos is the trophic position of individual murres
data$tpos <- 2+((1/3.4)*(data$d15N-data$baseline_d15N))
#See Bolnick et al., 2014 (Nature Communications) for more information on how to calculate alpha and tpos
#Rearrange order of columns to make the rest of the script easier to execute
data <- data[c(1,2,3,4,5,7,9,11,12,13,6,8,15,10,14,16)]

#Create a separate data frame for the filtered relative frequency table in order to merge it with the isotope data
rel_freq=read.csv("filtered-ASVs.csv", header=T)
#Remove columns that contain ASVs with only zeros
rel_freq <- rel_freq[, colSums(rel_freq !=0) > 0]
#Merge the two data frames
data <- merge(data, rel_freq, by="SampleID")
#Saved a copy of the data file
write.csv(data, "data_GLMs.csv", row.names = FALSE)


#Statistical tests for differences in isotope ratios between sexes and seasons
females = subset(data, Sex == "F")
males = subset(data, Sex == "M")
sex = data[!(data$Sex==""),]

shapiro.test(females$alpha)
shapiro.test(females$tpos)
shapiro.test(females$d34S)
shapiro.test(males$alpha)
shapiro.test(males$tpos)
shapiro.test(males$d34S)

wilcox.test(alpha ~ Sex, data = sex)
wilcox.test(d34S ~ Sex, data = sex)
bartlett.test(tpos ~ Sex, data = sex)
t.test(tpos ~ Sex, data = sex)

chick = subset(data, Status == "Chick")
egg = subset(data, Status == "Incubating")

shapiro.test(chick$alpha)
shapiro.test(chick$tpos)
shapiro.test(chick$d34S)
shapiro.test(egg$alpha)
shapiro.test(egg$tpos)
shapiro.test(egg$d34S)

wilcox.test(alpha ~ Status, data = data)
wilcox.test(d34S ~ Status, data = data)
bartlett.test(tpos ~ Status, data = data)
wilcox.test(tpos ~ Status, data = data)


#Quasibinomial GLMs

#Define ASVs to be analyzed (outcome)
out_start=17
out_end= 440
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se=rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)
out_AIC=rep(NA, out_nvar)
out_AICc=rep(NA, out_nvar)

#Define variables to be used (exposure)
exp_start=14
exp_end=16
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se=rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)
exp_AIC=rep(NA, exp_nvar)
exp_AICc=rep(NA, exp_nvar)

number=1

#Create "hacked" quasibinomial function so that AIC values from the binomial model are used
#Quasibinomial models cannot produce AIC values so this is the best alternative
hacked.quasibinomial <- function(...){
  res <- quasibinomial(...)
  res$aic <- binomial(...)$aic
  res
}
library(MuMIn)

#Loop to run a model for each variable and ASV combination
for (i in out_start:out_end){
  outcome = colnames(data)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(data)[j]
    model <- glm(formula = get(outcome) ~ get(exposure), family = hacked.quasibinomial,
                 maxit = 100,
                 na.action = na.exclude,
                 data=data)
    
    out_beta[number] = as.numeric(coef(model)[2])
    out_se[number] = as.numeric(coef(summary(model))[2,2])
    out_pvalue[number] = as.numeric(coef(summary(model))[2,4])
    out_AIC[number] = AIC(model)
    out_AICc[number] = AICc(model)
    out_variable[number] = outcome
    
    exp_beta[number] = as.numeric(coef(model)[2])
    exp_se[number] = as.numeric(coef(summary(model))[2,2])
    exp_pvalue[number] = as.numeric(coef(summary(model))[2,4])
    exp_AIC[number] = AIC(model)
    exp_AICc[number] = AICc(model)
    exp_variable[number] = exposure
    number = number + 1
  }
}

outcome = data.frame(out_variable, out_beta, out_se, out_pvalue, out_AIC, out_AICc)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue, exp_AIC, exp_AICc)

#Rename the columns of the results and merge them into a single dataset
library(tidyverse)
outcome = outcome %>% 
  rename(
    ASV = out_variable,
    estimate = out_beta,
    se = out_se,
    pvalue = out_pvalue,
    AIC = out_AIC,
    AICc = out_AICc
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    estimate = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue,
    AIC = exp_AIC,
    AICc = exp_AICc
  )
library(plyr); library(dplyr)
glms = join(outcome,exposure)
detach("package:plyr", unload = TRUE)
library(tidyverse)
glms = glms %>%
  distinct()
glms$qvalue <- p.adjust(glms$pvalue, method = "fdr", n = length(glms$pvalue))
glms = select(glms, ASV, variable, estimate, se, pvalue, qvalue, AIC, AICc)

#Load taxonomy file
ASVs <- read.delim("taxonomy.tsv", na.strings="")
#Remove ASVs for which no GLMs were performed on
ASVs = ASVs %>%
  filter(ASV %in% glms$ASV)
#Save modified taxonomy file for backup
write.csv(ASVs, "ASVs.csv", row.names = FALSE)
#Merge GLM results with the taxonomic information of each ASV
glms = merge(glms, ASVs, by="ASV")
#Save unified results dataset for backup
write.csv(glms, "glms.csv", row.names = FALSE)

#Quasibinomial GLMs with quadratic effect of d34S + linear effect of alpha

#Define ASVs to be analyzed (quadratic model outcome)
q_out_variable=rep(NA, out_nvar)
out_beta_alpha=rep(NA, out_nvar)
out_se_alpha=rep(NA, out_nvar)
out_pvalue_alpha=rep(NA, out_nvar)
out_beta_d34S1=rep(NA, out_nvar)
out_se_d34S1=rep(NA, out_nvar)
out_pvalue_d34S1=rep(NA, out_nvar)
out_beta_d34S2=rep(NA, out_nvar)
out_se_d34S2=rep(NA, out_nvar)
out_pvalue_d34S2=rep(NA, out_nvar)
q_out_AIC=rep(NA, out_nvar)
q_out_AICc=rep(NA, out_nvar)

number=1

#Loop to run a quadratic model for each variable and ASV combination
for (i in out_start:out_end){
  outcome = colnames(data)[i]
  model2 <- glm(formula = get(outcome) ~ alpha + poly(d34S, degree = 2), family = hacked.quasibinomial,
                maxit = 700,
                na.action = na.exclude,
                data=data)
  
  out_beta_alpha[number] = as.numeric(coef(model2)[2])
  out_se_alpha[number] = as.numeric(coef(summary(model2))[2,2])
  out_pvalue_alpha[number] = as.numeric(coef(summary(model2))[2,4])
  out_beta_d34S1[number] = as.numeric(coef(model2)[3])
  out_se_d34S1[number] = as.numeric(coef(summary(model2))[3,2])
  out_pvalue_d34S1[number] = as.numeric(coef(summary(model2))[3,4])
  out_beta_d34S2[number] = as.numeric(coef(model2)[4])
  out_se_d34S2[number] = as.numeric(coef(summary(model2))[4,2])
  out_pvalue_d34S2[number] = as.numeric(coef(summary(model2))[4,4])
  q_out_AIC[number] = AIC(model2)
  q_out_AICc[number] = AICc(model2)
  q_out_variable[number] = outcome
  number = number + 1
}

quadratic = data.frame(q_out_variable, out_beta_alpha, out_se_alpha, out_pvalue_alpha, out_beta_d34S1, out_se_d34S1, out_pvalue_d34S1, out_beta_d34S2, out_se_d34S2, out_pvalue_d34S2, q_out_AIC, q_out_AICc)
quadratic = quadratic %>% 
  rename(
    ASV = q_out_variable,
    alpha_estimate = out_beta_alpha,
    alpha_se = out_se_alpha,
    alpha_pvalue = out_pvalue_alpha,
    d34S1_estimate = out_beta_d34S1,
    d34S1_se = out_se_d34S1,
    d34S1_pvalue = out_pvalue_d34S1,
    d34S2_estimate = out_beta_d34S2,
    d34S2_se = out_se_d34S2,
    d34S2_pvalue = out_pvalue_d34S2,
    AIC = q_out_AIC,
    AICc = q_out_AICc
  )

quadratic$alpha_qvalue <- p.adjust(quadratic$alpha_pvalue, method = "fdr", n = length(quadratic$alpha_pvalue))
quadratic$d34S1_qvalue <- p.adjust(quadratic$d34S1_pvalue, method = "fdr", n = length(quadratic$d34S1_pvalue))
quadratic$d34S2_qvalue <- p.adjust(quadratic$d34S2_pvalue, method = "fdr", n = length(quadratic$d34S2_pvalue))
quadratic = select(quadratic, ASV, alpha_estimate, alpha_se, alpha_pvalue, alpha_qvalue, d34S1_estimate, d34S1_se, d34S1_pvalue, d34S1_qvalue, d34S2_estimate, d34S2_se, d34S2_pvalue, d34S2_qvalue, AIC, AICc)
quadratic = merge(quadratic, ASVs, by="ASV")
#Save quadratic GLM results for backup
write.csv(quadratic, "quadratic.csv", row.names = FALSE)

#Create contingency tables for chi-squared tests
variables = c("d34S", "alpha", "tpos")
a = 0.05
c = 0.95

name_var = rep(NA, 4)
n_sig = rep(NA, 4)
n_nosig = rep(NA, 4)
prop_sig = rep(NA, 4)
p_value = rep(NA, 4)
n_sig_q = rep(NA, 4)
n_nosig_q = rep(NA, 4)
prop_sig_q = rep(NA, 4)
q_value = rep(NA, 4)

number = 1

for(num in 1:3){
  var = variables[num]
  b = nrow(glms[glms$variable == var & glms$pvalue <= 0.05,])
  d = nrow(glms[glms$variable == var & glms$pvalue > 0.05,])
  e = nrow(glms[glms$variable == var & glms$qvalue <= 0.05,])
  f = nrow(glms[glms$variable == var & glms$qvalue > 0.05,])
  contingency = c(b,d)
  expected = c(a,c)
  prop = b/(b+d)
  chtest = chisq.test(contingency, p = expected)
  p_val = chtest$p.value
  contingency_q = c(e,f)
  prop_q = e/(e+f)
  chtest_q = chisq.test(contingency_q, p = expected)
  q_val = chtest_q$p.value
  
  name_var[number] = var
  n_sig[number] = b
  n_nosig[number] = d
  prop_sig[number] = prop
  p_value[number] = p_val
  n_sig_q[number] = e
  n_nosig_q[number] = f
  prop_sig_q[number] = prop_q
  q_value[number] = q_val
  
  number = number + 1
}

var = "aS2"
b = nrow(quadratic[quadratic$alpha_pvalue <= 0.05 & quadratic$d34S1_pvalue <= 0.05 & quadratic$d34S2_pvalue <= 0.05,])
d = nrow(quadratic[quadratic$alpha_pvalue > 0.05 | quadratic$d34S1_pvalue > 0.05 | quadratic$d34S2_pvalue > 0.05,])
e = nrow(quadratic[quadratic$alpha_qvalue <= 0.05 & quadratic$d34S1_qvalue <= 0.05 & quadratic$d34S2_qvalue <= 0.05,])
f = nrow(quadratic[quadratic$alpha_qvalue > 0.05 | quadratic$d34S1_qvalue > 0.05 | quadratic$d34S2_qvalue > 0.05,])
contingency = c(b,d)
expected = c(a,c)
prop = b/(b+d)
chtest = chisq.test(contingency, p = expected)
p_val = chtest$p.value
contingency_q = c(e,f)
prop_q = e/(e+f)
chtest_q = chisq.test(contingency_q, p = expected)
q_val = chtest_q$p.value

name_var[number] = var
n_sig[number] = b
n_nosig[number] = d
prop_sig[number] = prop
p_value[number] = p_val
n_sig_q[number] = e
n_nosig_q[number] = f
prop_sig_q[number] = prop_q
q_value[number] = q_val

#Combine results
contingency_tests = data.frame(name_var, n_sig, n_nosig, prop_sig, p_value, n_sig_q, n_nosig_q, prop_sig_q, q_value)
contingency_tests = contingency_tests %>%
  rename(
    variable = name_var,
    significant = n_sig,
    non_significant = n_nosig,
    significant_proportion = prop_sig,
    pvalue = p_value,
    significant_FDR = n_sig_q,
    non_significant_FDR = n_nosig_q,
    significant_proportion_FDR = prop_sig_q,
    pvalue_FDR = q_value
  )
#Save file for backup
write.csv(contingency_tests, "contingency.csv", row.names = FALSE)

#Create table and figure to observe which isotope models have positive or negative effects on ASV abundance

#For single variable GLMs
Effect = rep(NA, nrow(glms))

for (i in 1:nrow(glms)){
  if (glms$pvalue[i] > 0.05){
    Effect[i] = 0
  } else {
    if (glms$estimate[i] > 0){
      Effect[i] = 1
    } else {
      Effect[i] = -1
    }
  }
} 

glm_effects = data.frame(glms$ASV, glms$Class, glms$variable, Effect)
glm_effects = glm_effects %>%
  rename(ASV = glms.ASV,
         Class = glms.Class,
         Variable = glms.variable)

#for quadratic GLMs
q_effect_alpha = rep(NA, nrow(quadratic))
alpha_variable = rep("alpha_q", nrow(quadratic))
q_effect_d34S1 = rep(NA, nrow(quadratic))
d34S1_variable = rep("d34S1_q", nrow(quadratic))
q_effect_d34S2 = rep(NA, nrow(quadratic))
d34S2_variable = rep("d34S2_q", nrow(quadratic))

for (i in 1:nrow(quadratic)){
  if (quadratic$alpha_pvalue[i] > 0.05 || quadratic$d34S1_pvalue[i] > 0.05 || quadratic$d34S2_pvalue[i] > 0.05){
    q_effect_alpha[i] = 0
    q_effect_d34S1[i] = 0
    q_effect_d34S2[i] = 0
  } else {
    if (quadratic$alpha_estimate[i] > 0){
      q_effect_alpha[i] = 1
    } else {
      q_effect_alpha[i] = -1
    }
  }
}

for (i in 1:nrow(quadratic)){
  if (quadratic$alpha_pvalue[i] > 0.05 || quadratic$d34S1_pvalue[i] > 0.05 || quadratic$d34S2_pvalue[i] > 0.05){
    q_effect_alpha[i] = 0
    q_effect_d34S1[i] = 0
    q_effect_d34S2[i] = 0
  } else {
    if (quadratic$d34S1_estimate[i] > 0){
      q_effect_d34S1[i] = 1
    } else {
      q_effect_d34S1[i] = -1
    }
  }
}

for (i in 1:nrow(quadratic)){
  if (quadratic$alpha_pvalue[i] > 0.05 || quadratic$d34S1_pvalue[i] > 0.05 || quadratic$d34S2_pvalue[i] > 0.05){
    q_effect_alpha[i] = 0
    q_effect_d34S1[i] = 0
    q_effect_d34S2[i] = 0
  } else {
    if (quadratic$d34S2_estimate[i] > 0){
      q_effect_d34S2[i] = 1
    } else {
      q_effect_d34S2[i] = -1
    }
  }
}


glm_effects_alpha = data.frame(quadratic$ASV, quadratic$Class, alpha_variable, q_effect_alpha)
glm_effects_alpha = glm_effects_alpha %>%
  rename(ASV = quadratic.ASV,
         Class = quadratic.Class,
         Variable = alpha_variable,
         Effect = q_effect_alpha)

glm_effects_d34S1 = data.frame(quadratic$ASV, quadratic$Class, d34S1_variable, q_effect_d34S1)
glm_effects_d34S1 = glm_effects_d34S1 %>%
  rename(ASV = quadratic.ASV,
         Class = quadratic.Class,
         Variable = d34S1_variable,
         Effect = q_effect_d34S1)

glm_effects_d34S2 = data.frame(quadratic$ASV, quadratic$Class, d34S2_variable, q_effect_d34S2)
glm_effects_d34S2 = glm_effects_d34S2 %>%
  rename(ASV = quadratic.ASV,
         Class = quadratic.Class,
         Variable = d34S2_variable,
         Effect = q_effect_d34S2)

#Combine the linear and quadratic GLMs
glm_effects = rbind(glm_effects, glm_effects_alpha, glm_effects_d34S1, glm_effects_d34S2)
glm_effects = glm_effects[order(glm_effects$ASV),]
#Rename variables so that the three components of the quadratic model are grouped next to each other when plotting the results
glm_effects$Variable <- as.character(glm_effects$Variable)
glm_effects$Variable[glm_effects$Variable == "tpos"] <- "btpos"
glm_effects$Variable[glm_effects$Variable == "d34S"] <- "cd34S"
glm_effects$Variable[glm_effects$Variable == "alpha_q"] <- "d1alpha_q"
#Save file for backup
write.csv(glm_effects, "glm_effects.csv", row.names = FALSE)

#Rename variables for the figure and add number of models with an effect
facet_labels = c("alpha (23*)", "tpos (40*)", "d34S (49*)", "alpha (88**)", "d34S (88**)", "d34S^2 (88**)")
names(facet_labels) = c("alpha", "btpos", "cd34S", "d1alpha_q", "d34S1_q", "d34S2_q")
addline_format = function(x,...){
  gsub('\\s', '\n', x)
}
Effects <- ggplot(data = glm_effects, mapping = aes(x = ASV,
                                                    y = reorder(Class, desc(Class)),
                                                    fill = Effect)) +
  geom_tile() +
  ylab(label = "Class") +
  facet_grid(~ Variable, switch = "x", labeller = labeller(Variable = addline_format(facet_labels))) +
  scale_fill_gradient2(Effect, low = "blue", mid = "white", high = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
Effects


#DESeq analysis

#Copy the sample-metadata-full.txt, filtered-table.qza, taxonomy-16S-silva-132.qza, rooted-tree-fasttree.qza into a data folder
#Qiime2 artifacts were converted into phyloseq format using the qiime2R package
library(qiime2R)
library(readr)
metadata = read_tsv("sample-metadata-full.txt")
metadata
SV = read_qza("filtered-table.qza")
names(SV)
taxonomy = read_qza("taxonomy-16S-silva-132.qza")
taxonomy$uuid
taxtable = taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxtable
tree = read_qza("rooted-tree-fasttree.qza")
tree$uuid
tree$data

library(phyloseq)
library(DESeq2)

#Create a phyloseq object with our dataset
pso_full = phyloseq(
  otu_table(SV$data, taxa_are_rows = T),
  phy_tree(tree$data),
  tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()),
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("#SampleID"))
)

#Subset samples to only include those from the first replicate of each bird
pso = subset_samples(pso_full, Replicate == 1)

#DESeq for differences between sexes
#Subset the phyloseq object to remove the individuals that were not sexed
pso_sex = subset_samples(pso, Sex == "F" | Sex == "M")
pso_sex

#Convert phyloseq object to DESeq format and run DESeq
dso_s = phyloseq_to_deseq2(pso_sex, ~ Sex)
dso_s = DESeq(dso_s, test = "Wald", fitType = "parametric")

#Filter results to include only those with p value of less than 0.05
res = results(dso_s, cooksCutoff = FALSE)
alph = 0.05
sigtab_sex = res[which(res$padj < alph), ]
sigtab_sex = cbind(as(sigtab_sex, "data.frame"), as(tax_table(pso_sex)[rownames(sigtab_sex), ], "matrix"))
head(sigtab_sex)
dim(sigtab_sex)
#Fix Silva Taxonomy
sigtab_sex$Genus[sigtab_sex$Genus == "D_5__uncultured bacterium"] <- as.character(NA)
sigtab_sex$Genus[sigtab_sex$Genus == "D_5__uncultured"] <- as.character(NA)
levels(sigtab_sex$Kingdom) = str_replace_all(levels(sigtab_sex$Kingdom), "D_[0-9]__", "")
levels(sigtab_sex$Phylum) = str_replace_all(levels(sigtab_sex$Phylum), "D_[0-9]__", "")
levels(sigtab_sex$Class) = str_replace_all(levels(sigtab_sex$Class), "D_[0-9]__", "")
levels(sigtab_sex$Order) = str_replace_all(levels(sigtab_sex$Order), "D_[0-9]__", "")
levels(sigtab_sex$Family) = str_replace_all(levels(sigtab_sex$Family), "D_[0-9]__", "")
levels(sigtab_sex$Genus) = str_replace_all(levels(sigtab_sex$Genus), "D_[0-9]__", "")
levels(sigtab_sex$Species) = str_replace_all(levels(sigtab_sex$Species), "D_[0-9]__", "")
head(sigtab_sex)
dim(sigtab_sex)

#Create figure
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_sex$log2FoldChange, sigtab_sex$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_sex$Phylum = factor(as.character(sigtab_sex$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_sex$log2FoldChange, sigtab_sex$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_sex$Genus = factor(as.character(sigtab_sex$Genus), levels=names(x))
DE_sex = ggplot(sigtab_sex, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
  ) +
  labs (y = "log2 Fold Change (Male vs Female)")
DE_sex

#DESeq for differences between reproductive stages
#Convert phyloseq object to DESeq format and run DESeq
dso_stat = phyloseq_to_deseq2(pso, ~ Status)
dso_stat = DESeq(dso_stat, test = "Wald", fitType = "parametric")

#Filter results to include only those with p value of less than 0.05
res_stat = results(dso_stat, cooksCutoff = FALSE)
sigtab_status = res_stat[which(res_stat$padj < alph), ]
sigtab_status = cbind(as(sigtab_status, "data.frame"), as(tax_table(pso)[rownames(sigtab_status), ], "matrix"))
head(sigtab_status)
dim(sigtab_status)
#Fix Silva Taxonomy
sigtab_status$Genus[sigtab_status$Genus == "D_5__uncultured bacterium"] <- as.character(NA)
sigtab_status$Genus[sigtab_status$Genus == "D_5__uncultured"] <- as.character(NA)
levels(sigtab_status$Kingdom) = str_replace_all(levels(sigtab_status$Kingdom), "D_[0-9]__", "")
levels(sigtab_status$Phylum) = str_replace_all(levels(sigtab_status$Phylum), "D_[0-9]__", "")
levels(sigtab_status$Class) = str_replace_all(levels(sigtab_status$Class), "D_[0-9]__", "")
levels(sigtab_status$Order) = str_replace_all(levels(sigtab_status$Order), "D_[0-9]__", "")
levels(sigtab_status$Family) = str_replace_all(levels(sigtab_status$Family), "D_[0-9]__", "")
levels(sigtab_status$Genus) = str_replace_all(levels(sigtab_status$Genus), "D_[0-9]__", "")
levels(sigtab_status$Species) = str_replace_all(levels(sigtab_status$Species), "D_[0-9]__", "")
head(sigtab_status)
dim(sigtab_status)

#Create figure
# Phylum order
x = tapply(sigtab_status$log2FoldChange, sigtab_status$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_status$Phylum = factor(as.character(sigtab_status$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_status$log2FoldChange, sigtab_status$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_status$Genus = factor(as.character(sigtab_status$Genus), levels=names(x))
DE_status = ggplot(sigtab_status, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(size=3) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
  ) +
  labs (y = "log2 Fold Change (Chick vs Incubating)")
DE_status

#Combine the two figures into a single object
library(gridExtra)
grid.arrange(DE_sex, DE_status, ncol=1)


#NMDS Ordination plots

#NMDS ordination with unweighted UniFrac metric
murre_unw_ord = ordinate(pso, "NMDS", "unifrac")
stress_unw = murre_unw_ord$stress
stress_unw

#Plot the unweighted UniFrac NMDS ordination
group_colors = c("#08f002", "#f0027f")
unifrac_nmds = plot_ordination(pso, murre_unw_ord, type = "samples", color = "Sex", shape = "Status") +
  geom_point(size=4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()
  ) +
  scale_color_manual(values = group_colors, na.value = "#02f0ea", labels = c("Female", "Male", "Unknown")) +
  scale_shape_discrete(name = "Stage", labels = c("Chick-rearing", "Incubating")) +
  annotate("text", x=0.32, y=-0.32, label="Stress = 0.21", size = 4)
unifrac_nmds

#NMDS ordination with weighted UniFrac metric
murre_w_ord = ordinate(pso, "NMDS", "wunifrac")
stress_w = murre_w_ord$stress
stress_w

#Plot the weighted UniFrac NMDS ordination
wunifrac_nmds = plot_ordination(pso, murre_w_ord, type = "samples", color = "Sex", shape = "Status") +
  geom_point(size=4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()
  ) +
  scale_color_manual(values = group_colors, na.value = "#02f0ea", labels = c("Female", "Male", "Unknown")) +
  scale_shape_discrete(name = "Stage", labels = c("Chick-rearing", "Incubating")) +
  annotate("text", x=0.32, y=-0.11, label="Stress = 0.06", size = 4)
wunifrac_nmds

#Combine the two figures into a single object
grid.arrange(unifrac_nmds, wunifrac_nmds, ncol=1)