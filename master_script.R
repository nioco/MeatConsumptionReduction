install.packages("here")
# Declare the location of the current script
here::i_am("master_script.R")

library(here)
# Check if location was declared correctly
here()

install.packages("haven")
install.packages("plyr")
install.packages("dplyr")
install.packages("psych")
install.packages("corrplot")
install.packages("lavaan")
install.packages("labelled")
install.packages("readr")
install.packages("QuantPsyc")
install.packages("RColorBrewer")
install.packages("semTools")
install.packages("devtools")
install.packages("htmlTable")
library(haven)
library(plyr)
library(dplyr)
library(psych) 
library(corrplot)
library(ggplot2)
library(lavaan) 
library(labelled)
library(readr)
library(QuantPsyc)
library(RColorBrewer)
library(semTools)
library(devtools)
library(htmlTable)

install("lavaan.mi")
library(lavaan.mi)

install("dynamic/dynamic-master")
library(dynamic)

# Read data, convert to csv and read again to deal with Stata format issues
data <- read_dta("data/appended_clean.dta")
class(data)
write.csv(data, "data/data.csv")
data_csv <- read_csv("data/data.csv")

################################################################################
### Data Cleaning ###
################################################################################
# Only keep baseline and wave 1
str(data_csv$wave)
data_csv <- subset(data_csv, wave <= 1)
str(data_csv$wave)

# Only keep control group subjects
str(data_csv$treat)
summary(data_csv$treat)
data_csv <- subset(data_csv, treat == 4)
str(data_csv$treat)

# Subdatensatz mit Indikatoren_0 erstellen
indicators <- data_csv %>% 
                dplyr::select(wave,
                              meat_freq, 
                              meat, 
                              attitude1, 
                              attitude2, 
                              attitude3, 
                              intention1, 
                              intention2, 
                              intention3, 
                              injnorm1, 
                              injnorm2, 
                              injnorm3,
                              pbc1, 
                              pbc2, 
                              pbc3, 
                              identity1, 
                              identity2,
                              loc_1_rev,
                              loc_2_rev,
                              loc_3_rev,
                              loc_4_rev,
                              loc_5_rev,
                              loc_6,
                              loc_7,
                              key)

# Nur Cases mit wave = 0 behalten
str(indicators$wave)
indicators <- subset(indicators, wave == 0)
str(indicators$wave)
indicators <- indicators %>% 
  dplyr::select(-wave)

# Add Meat-Consumption Behavior from Wave 1
indicators_1 <- subset(data_csv, wave == 1)
indicators_1 <- indicators_1 %>% 
  dplyr::select(meat_freq, 
                meat,
                key,
                wave) 
str(indicators_1$wave)

indicators_1 <- indicators_1 %>% 
  dplyr::select(-wave)
names(indicators_1) <- c("meat_freq_1", "meat_1", "key")

# Merge datasets
indicators <- merge(indicators, indicators_1, by = "key")

# Descriptive Statistics
summary(indicators)
sapply(indicators[,2:26], sd, na.rm = TRUE)

# Grafische Veranschaulichung der Verteilung
multi.hist(indicators, 
           global = FALSE)

# Korrelationen zwischen Indikatoren_0
corr1 <- cor(indicators, 
             use = "complete.obs")
corrplot(corr1, 
         method="number",
         type = "lower",
         col=brewer.pal(n=8, name="RdYlBu"),
         tl.col="black",
         tl.srt=45,
         diag=FALSE)

# Percentage of missing cases
missing_cells <- sum(is.na(indicators))
total_cells <- prod(dim(indicators)) 
percent_missing <- (missing_cells * 100 )/(total_cells) 
percent_missing
################################################################################
################################################################################
#### Confirmatory Factor Analyses ####
################################################################################
################################################################################

# Check for Estimator Bias for IJN

inj_cfa <- 'IJN =~ injnorm1 + injnorm2 + injnorm3
'

fit_inj_cfa <- cfa(model = inj_cfa,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
parameterEstimates(fit_inj_cfa)
standardizedSolution(fit_inj_cfa)

indicators_cat <- indicators %>% 
  dplyr::select(-key)
fit_inj_cfa_cat <- semTools::cfa.mi(model = inj_cfa,
                   data = indicators_cat,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   ordered = TRUE,
                   estimator = "ULSMVS",
                   m = 50,
                   miPackage = "mice",
                   seed = 1234)
parameterEstimates.mi(fit_inj_cfa_cat)
standardizedSolution.mi(fit_inj_cfa_cat)

### TPB Model

tpb_cfa <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc1 + pbc2 + pbc3
          INT =~ intention1 + intention2 + intention3
'
fit_tpb_cfa <- cfa(model = tpb_cfa,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_tpb_cfa, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_tpb_cfa,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_tpb_cfa)
# Convert to HTML
res_tpb_cfa <- format(round(lavResiduals(fit_tpb_cfa)$cov,
                            3), 
                      nsmall = 3)  
res_tpb_cfa[upper.tri(res_tpb_cfa)] <- ""
res_tpb_cfa %>% htmlTable() # correlation residuals > .10
lavInspect(fit_tpb_cfa, "sampstat")
res_tpb_cfa_sam <- format(round(lavInspect(fit_tpb_cfa,
                                           "sampstat")$cov, 
                                3), 
                          nsmall = 3)  
res_tpb_cfa_sam[upper.tri(res_tpb_cfa)] <- ""
res_tpb_cfa_sam %>% htmlTable()
modindices(fit_tpb_cfa, sort. = TRUE, op = "~~")

# CFA ohne pbc1
tpb_cfa_re <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention1 + intention2 + intention3
'
fit_tpb_cfa_re <- cfa(model = tpb_cfa_re,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_tpb_cfa_re, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_tpb_cfa_re,
      data = indicators,
      estimator = "MLR")

# Residuals
lavInspect(fit_tpb_cfa_re, "sampstat")
lavResiduals(fit_tpb_cfa_re) # no correlation residuals > .10 remain
# Convert to HTML
res_tpb_cfa_re <- format(round(lavResiduals(fit_tpb_cfa_re)$cov, 3), nsmall = 3)  
res_tpb_cfa_re[upper.tri(res_tpb_cfa_re)] <- ""
res_tpb_cfa_re %>% htmlTable()
modindices(fit_tpb_cfa_re, sort. = TRUE)

################################################################################

### TPB Model with Meat-Eater Identity

mei_cfa <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention1 + intention2 + intention3
          MEI =~ identity1 + identity2
'
fit_mei_cfa <- cfa(model = mei_cfa,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_mei_cfa, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mei_cfa,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mei_cfa)
res_mei_cfa <- lavResiduals(fit_mei_cfa)$cov 
View(res_mei_cfa)
# correlation residuals > .10 between intention1 and identity2
# Convert to HTML
res_mei_cfa <- format(round(lavResiduals(fit_mei_cfa)$cov, 3), nsmall = 3)  
res_mei_cfa[upper.tri(res_mei_cfa)] <- ""
res_mei_cfa %>% htmlTable()

lavInspect(fit_mei_cfa, "sampstat") 
# intention2 and 3 have higher and more similar covariances with MEI items than 
# intention1
res_mei_cfa_sam <- format(round(lavInspect(fit_mei_cfa, 
                                           "sampstat")$cov,
                                3), 
                          nsmall = 3)  
res_mei_cfa_sam[upper.tri(res_mei_cfa_sam)] <- ""
res_mei_cfa_sam %>% htmlTable()
modindices(fit_mei_cfa, sort. = TRUE)

# Respecification without intention1

mei_cfa_re <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~  intention2 + intention3
          MEI =~ identity1 + identity2
'
fit_mei_cfa_re <- cfa(model = mei_cfa_re,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_mei_cfa_re, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mei_cfa_re,
      data = indicators,
      estimator = "MLR")

# Residuals
lavInspect(fit_mei_cfa_re, 
           "sampstat")
lavResiduals(fit_mei_cfa_re) # no correlation residuals > .10 remain
# Convert to HTML
res_mei_cfa_re <- format(round(lavResiduals(fit_mei_cfa_re)$cov, 
                               3), 
                         nsmall = 3)  
res_mei_cfa_re[upper.tri(res_mei_cfa_re)] <- ""
res_mei_cfa_re %>% htmlTable()

################################################################################

### TPB Model with Meat-Eater Identity and Locus of Control

loc_cfa <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_3_rev + loc_4_rev + loc_5_rev + 
                 loc_6 + loc_7
'
fit_loc_cfa <- cfa(model = loc_cfa,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_loc_cfa, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_loc_cfa,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_loc_cfa) # 0.246 correlation residual between loc_6 and loc_7
# Convert to HTML
res_loc_cfa <- format(round(lavResiduals(fit_loc_cfa)$cov,
                            3), 
                      nsmall = 3)  
res_loc_cfa[upper.tri(res_loc_cfa)] <- ""
res_loc_cfa %>% htmlTable()
lavInspect(fit_loc_cfa, "sampstat") 
modindices(fit_loc_cfa, sort. = TRUE)

# Exploratory factor analysis
# Check LOC for dimensionality

loc_efa <- 'efa("efa")*f1 + 
            efa("efa")*f2 =~ loc_1_rev + loc_2_rev + loc_3_rev + loc_4_rev +
                             loc_5_rev + loc_6 + loc_7
'
fit_loc_efa <- cfa(loc_efa,
                   indicators, 
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_loc_efa,
        fit.measures=TRUE, 
        standardized= TRUE)
# loc_6 loads strongly on the second factor with some cross-loadings of loc_7

# Standardized estimates with standard errors and p-values
options(max.print = 2000)
standardizedSolution(fit_loc_efa)


# CFA without loc_6

loc_cfa_re <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_3_rev + loc_4_rev + loc_5_rev 
                + loc_7
'
fit_loc_cfa_re <- cfa(model = loc_cfa_re,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_loc_cfa_re, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_loc_cfa_re,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_loc_cfa_re)
# loc_3_rev and loc_4_rev have correlation residuals > .10
# Convert to HTML
res_loc_cfa_re <- format(round(lavResiduals(fit_loc_cfa_re)$cov,
                            3), 
                      nsmall = 3)  
res_loc_cfa_re[upper.tri(res_loc_cfa_re)] <- ""
res_loc_cfa_re %>% htmlTable()

lavInspect(fit_loc_cfa_re, "sampstat") 
# The covariances of loc_3_rev and loc_4_rev with items of other constructs
# differ from the rest of the LOC items
# Convert to HTML
res_loc_cfa_re_sam <- format(round(lavInspect(fit_loc_cfa_re, 
                                          "sampstat")$cov,
                               3), 
                         nsmall = 3)  
res_loc_cfa_re_sam[upper.tri(res_loc_cfa_re_sam)] <- ""
res_loc_cfa_re_sam %>% htmlTable()
modindices(fit_loc_cfa_re, sort. = TRUE)

# CFA without loc_6 and loc_4_rev

loc_cfa_re2 <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_3_rev + loc_5_rev + loc_7
'
fit_loc_cfa_re2 <- cfa(model = loc_cfa_re2,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_loc_cfa_re2, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)
# Barely non-significant X^2

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_loc_cfa_re2,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_loc_cfa_re2)
# loc_3_rev has a correlation residual > .10
# Convert to HTML
res_loc_cfa_re2 <- format(round(lavResiduals(fit_loc_cfa_re2)$cov,
                               3), 
                         nsmall = 3)  
res_loc_cfa_re2[upper.tri(res_loc_cfa_re2)] <- ""
res_loc_cfa_re2 %>% htmlTable()
modindices(fit_loc_cfa_re2, sort. = TRUE)

# CFA without loc_6, loc_4_rev and loc_3_rev

loc_cfa_re3 <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
'
fit_loc_cfa_re3 <- cfa(model = loc_cfa_re3,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_loc_cfa_re3, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)
# Non-significant X^2

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_loc_cfa_re3,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_loc_cfa_re3) # no correlation residuals > .10 remain
# Convert to HTML
res_loc_cfa_re3 <- format(round(lavResiduals(fit_loc_cfa_re3)$cov,
                                3), 
                          nsmall = 3)  
res_loc_cfa_re3[upper.tri(res_loc_cfa_re3)] <- ""
res_loc_cfa_re3 %>% htmlTable()

################################################################################
################################################################################
### Confirmatory Factor Analyses for Repeated MCB Measurements ###
################################################################################
################################################################################

### Configural MI

mcb_cfa_con <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + meat
          MCB_1 =~ 1*meat_freq_1 + meat_1
          
          # Error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
'
fit_mcb_cfa_con <- cfa(model = mcb_cfa_con,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_mcb_cfa_con, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mcb_cfa_con,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mcb_cfa_con) # no correlation residuals > .10
# Convert to HTML
res_mcb_con_cfa <- format(round(lavResiduals(fit_mcb_cfa_con)$cov, 
                                3), 
                          nsmall = 3)  
res_mcb_con_cfa[upper.tri(res_mcb_con_cfa)] <- ""
res_mcb_con_cfa %>% htmlTable()

################################################################################

### Weak MI

mcb_cfa_met <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
'
fit_mcb_cfa_met <- cfa(model = mcb_cfa_met,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_mcb_cfa_met, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mcb_cfa_met,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mcb_cfa_met) # no correlation residuals > .10
# Convert to HTML
res_mcb_met_cfa <- format(round(lavResiduals(fit_mcb_cfa_met)$cov, 
                                3), 
                          nsmall = 3)  
res_mcb_met_cfa[upper.tri(res_mcb_met_cfa)] <- ""
res_mcb_met_cfa %>% htmlTable()

################################################################################

### Strong MI

mcb_cfa_sca <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
'
fit_mcb_cfa_sca <- cfa(model = mcb_cfa_sca,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_mcb_cfa_sca, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mcb_cfa_sca,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mcb_cfa_sca) # no correlation residuals > .10
# Convert to HTML
res_mcb_sca_cfa <- format(round(lavResiduals(fit_mcb_cfa_sca)$cov, 
                                3), 
                          nsmall = 3)  
res_mcb_sca_cfa[upper.tri(res_mcb_sca_cfa)] <- ""
res_mcb_sca_cfa %>% htmlTable()

################################################################################

### Strict MI

mcb_cfa_res <-'ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
          
          # Equal error variances across waves
          meat_freq ~~ errVar1*meat_freq
          meat_freq_1 ~~ errVar1*meat_freq_1
          meat ~~ errVar2*meat
          meat_1 ~~ errVar2*meat_1
'
fit_mcb_cfa_res <- cfa(model = mcb_cfa_res,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_mcb_cfa_res, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mcb_cfa_res,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mcb_cfa_res) # no correlation residuals > .10
# Convert to HTML
res_mcb_res_cfa <- format(round(lavResiduals(fit_mcb_cfa_res)$cov, 
                                3), 
                          nsmall = 3)  
res_mcb_res_cfa[upper.tri(res_mcb_res_cfa)] <- ""
res_mcb_res_cfa %>% htmlTable()

################################################################################
################################################################################
### Structural Equation Models ###
################################################################################
################################################################################

### Model with effects of TPB predictors

tpb_sem <-'# Measurment Model
          ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
          
          # Equal error variances across waves
          meat_freq ~~ errVar1*meat_freq
          meat_freq_1 ~~ errVar1*meat_freq_1
          meat ~~ errVar2*meat
          meat_1 ~~ errVar2*meat_1
          
          # Structural Model
          INT ~ attInt*ATT + ijnInt*IJN + pbcInt*PBC
          MCB_1 ~ intMcb*INT + MCB
          
          # Indirect Effects
          ind_att_int_mcb := attInt * intMcb
          ind_ijn_int_mcb := ijnInt * intMcb
          ind_pbc_int_mcb := pbcInt * intMcb
'
fit_tpb_sem <- sem(model = tpb_sem,
                      data = indicators,
                      effect.coding = "loadings", # Coding Effects Method Constraints
                      estimator = "MLR",
                      missing = "FIML")
summary(fit_tpb_sem, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_tpb_sem,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_tpb_sem) # no correlation residuals > .10
# Convert to HTML
res_tpb_sem <- format(round(lavResiduals(fit_tpb_sem)$cov,
                            3), 
                      nsmall = 3)  
res_tpb_sem[upper.tri(res_tpb_sem)] <- ""
res_tpb_sem %>% htmlTable()

################################################################################

### Model with effects for TPB predictors and Meat-Eater Identity

mei_sem <-'# Measurment Model
          ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
          
          # Equal error variances across waves
          meat_freq ~~ errVar1*meat_freq
          meat_freq_1 ~~ errVar1*meat_freq_1
          meat ~~ errVar2*meat
          meat_1 ~~ errVar2*meat_1
          
          # Structural Model
          ATT ~ meiATT*MEI
          IJN ~ meiIjn*MEI
          PBC ~ meiPbc*MEI
          INT ~ attInt*ATT + ijnInt*IJN + pbcInt*PBC + meiInt*MEI
          MCB_1 ~ intMcb*INT + MCB
          
          # Indirect Effects
          ind_mei_att_int := meiATT * attInt
          ind_mei_ijn_int := meiIjn * ijnInt
          ind_mei_pbc_int := meiPbc * pbcInt
          ind_att_int_mcb := attInt * intMcb
          ind_ijn_int_mcb := ijnInt * intMcb
          ind_pbc_int_mcb := pbcInt * intMcb
          ind_mei_att_int_mcb := meiATT * attInt * intMcb
          ind_mei_ijn_int_mcb := meiIjn * ijnInt * intMcb
          ind_mei_pbc_int_mcb := meiPbc * pbcInt * intMcb
'
fit_mei_sem <- sem(model = mei_sem,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_mei_sem, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mei_sem,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mei_sem) # several absolute correlation residuals > .10
# Convert to HTML
res_mei_sem <- format(round(lavResiduals(fit_mei_sem)$cov,
                            3), 
                      nsmall = 3)  
res_mei_sem[upper.tri(res_mei_sem)] <- ""
res_mei_sem %>% htmlTable()

lavInspect(fit_mei_sem, "sampstat") 
# Model underestimates covariances between items --> MEI doesn't explain the 
# covariances between ATT, IJN and PBC well
# Convert to HTML
res_mei_sem_sam <- format(round(lavInspect(fit_mei_sem,
                                           "sampstat")$cov,
                            3), 
                      nsmall = 3)  
res_mei_sem_sam[upper.tri(res_mei_sem_sam)] <- ""
res_mei_sem_sam %>% htmlTable()
modindices(fit_mei_sem_wol, 
           sort. = TRUE)

# respecification with disturbance covariances
mei_sem_re <-'# Measurment Model
          ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
          
          # Equal error variances across waves
          meat_freq ~~ errVar1*meat_freq
          meat_freq_1 ~~ errVar1*meat_freq_1
          meat ~~ errVar2*meat
          meat_1 ~~ errVar2*meat_1
          
          # Structural Model
          ATT ~ meiATT*MEI
          IJN ~ meiIjn*MEI
          PBC ~ meiPbc*MEI
          INT ~ attInt*ATT + ijnInt*IJN + pbcInt*PBC + meiInt*MEI
          MCB_1 ~ intMcb*INT + MCB
          
          # Disturbance Covariances
          ATT ~~ IJN + PBC
          IJN ~~ PBC
          
          # Indirect Effects
          ind_mei_att_int := meiATT * attInt
          ind_mei_ijn_int := meiIjn * ijnInt
          ind_mei_pbc_int := meiPbc * pbcInt
          ind_att_int_mcb := attInt * intMcb
          ind_ijn_int_mcb := ijnInt * intMcb
          ind_pbc_int_mcb := pbcInt * intMcb
          ind_mei_att_int_mcb := meiATT * attInt * intMcb
          ind_mei_ijn_int_mcb := meiIjn * ijnInt * intMcb
          ind_mei_pbc_int_mcb := meiPbc * pbcInt * intMcb
'
fit_mei_sem_re <- sem(model = mei_sem_re,
                       data = indicators,
                       effect.coding = "loadings", # Coding Effects Method Constraints
                       estimator = "MLR",
                       missing = "FIML")
summary(fit_mei_sem_re, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_mei_sem_re,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_mei_sem_re) # no absolute correlation residuals > .10 remain
# Convert to HTML
res_mei_sem_re <- format(round(lavResiduals(fit_mei_sem_re)$cov,
                            3), 
                      nsmall = 3)  
res_mei_sem_re[upper.tri(res_mei_sem_re)] <- ""
res_mei_sem_re %>% htmlTable()
lavInspect(fit_mei_sem_re, "sampstat") 

################################################################################

### Full model with effects for TPB predictors and Meat-Eater Identity, 

full_sem <-'# Measurment Model
          ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          LOC =~ loc_1_rev + loc_2_rev + loc_5_rev + loc_7
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
          
          # Equal error variances across waves
          meat_freq ~~ errVar1*meat_freq
          meat_freq_1 ~~ errVar1*meat_freq_1
          meat ~~ errVar2*meat
          meat_1 ~~ errVar2*meat_1
          
          # Structural Model
          ATT ~ meiATT*MEI + locATT*LOC
          IJN ~ meiIjn*MEI
          PBC ~ meiPbc*MEI + locPbc*LOC
          INT ~ attInt*ATT + ijnInt*IJN + pbcInt*PBC + meiInt*MEI + locInt*LOC
          MCB_1 ~ intMcb*INT + MCB

          # Disturbance Covariances
          ATT ~~ IJN + PBC
          IJN ~~ PBC
          
          # Indirect Effects
          ind_mei_att_int := meiATT * attInt
          ind_mei_ijn_int := meiIjn * ijnInt
          ind_mei_pbc_int := meiPbc * pbcInt
          ind_loc_att_int := locATT * attInt
          ind_loc_pbc_int := locPbc * pbcInt
          ind_att_int_mcb := attInt * intMcb
          ind_ijn_int_mcb := ijnInt * intMcb
          ind_pbc_int_mcb := pbcInt * intMcb
          ind_mei_att_int_mcb := meiATT * attInt * intMcb
          ind_mei_ijn_int_mcb := meiIjn * ijnInt * intMcb
          ind_mei_pbc_int_mcb := meiPbc * pbcInt * intMcb
          ind_loc_att_int_mcb := locATT * attInt * intMcb
          ind_loc_pbc_int_mcb := locPbc * pbcInt * intMcb
'
fit_full_sem <- sem(model = full_sem,
                   data = indicators,
                   effect.coding = "loadings", # Coding Effects Method Constraints
                   estimator = "MLR",
                   missing = "FIML")
summary(fit_full_sem, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Adjusted R-squared (Using code from the rsquareCalc function by Hayes (2021))
  # Adjusted R-squared for INT
    # Parameter estimates
    pe_full <- parameterEstimates(fit_full_sem, 
                                standardized = TRUE, 
                                rsquare = TRUE)
    # Grab unadjusted R-Squared
    Rsq_full_INT <- pe_full[pe_full$lhs == "INT" & pe_full$op == "r2", 
                    "est"]
    #Retrieve number of observations used in the analysis.
    n_full <- lavInspect(fit_full_sem, 
                  what = "nobs")
    #Number of predictors in the full model.
    p_full_INT <- nrow(pe_full[pe_full$lhs == "INT" & pe_full$op == "~",])
    #Adjusted R-square calculations.
    multiplier_full_INT <- (n_full-1)/(n_full - p_full_INT - 1)
    adj_Rsq_full_INT <- 1 - multiplier_full_INT*(1 - Rsq_full_INT)
    adj_Rsq_full_INT
    
  # Adjusted R-squared for MCB_1
    # Parameter estimates
    pe_full <- parameterEstimates(fit_full_sem, 
                                  standardized = TRUE, 
                                  rsquare = TRUE)
    # Grab unadjusted R-Squared
    Rsq_full_MCB_1 <- pe_full[pe_full$lhs == "MCB_1" & pe_full$op == "r2", 
                            "est"]
    #Retrieve number of observations used in the analysis.
    n_full <- lavInspect(fit_full_sem, 
                         what = "nobs")
    #Number of predictors in the full model.
    p_full_MCB_1 <- nrow(pe_full[pe_full$lhs == "MCB_1" & pe_full$op == "~",])
    #Adjusted R-square calculations.
    multiplier_full_MCB_1 <- (n_full-1)/(n_full - p_full_MCB_1 - 1)
    adj_Rsq_full_MCB_1 <- 1 - multiplier_full_MCB_1*(1 - Rsq_full_MCB_1)
    adj_Rsq_full_MCB_1

# Standardized estimates with standard errors and p-values
options(max.print = 2000)
standardizedSolution(fit_full_sem)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_full_sem,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_full_sem)  # no absolute correlation residuals > .10
# Convert to HTML
res_full_sem <- format(round(lavResiduals(fit_full_sem)$cov,
                               3), 
                         nsmall = 3)  
res_full_sem[upper.tri(res_full_sem)] <- ""
res_full_sem %>% htmlTable()
lavInspect(fit_full_sem, "sampstat") 

################################################################################
################################################################################
### Test Model ###
################################################################################
################################################################################

### Test Model without LOC (MEI only)

test_sem <-'# Measurment Model
          ATT =~ attitude1 + attitude2 + attitude3 
          IJN =~ injnorm1 + injnorm2 + injnorm3
          PBC =~ pbc2 + pbc3
          INT =~ intention2 + intention3
          MEI =~ identity1 + identity2
          MCB =~ 1*meat_freq + equf2*meat    # equal fator loadings for items
          MCB_1 =~ 1*meat_freq_1 + equf2*meat_1
          
          # Equal error covariances across waves
          meat_freq ~~ errCov1*meat_freq_1
          meat ~~ errCov2*meat_1
          
          # Equal intercepts across waves
          meat_freq ~ int1*1
          meat_freq_1 ~ int1*1
          meat ~ int2*1
          meat_1 ~ int2*1
          
          # Equal error variances across waves
          meat_freq ~~ errVar1*meat_freq
          meat_freq_1 ~~ errVar1*meat_freq_1
          meat ~~ errVar2*meat
          meat_1 ~~ errVar2*meat_1
          
          # Structural Model
          ATT ~ meiATT*MEI
          IJN ~ meiIjn*MEI
          PBC ~ meiPbc*MEI
          INT ~ attInt*ATT + ijnInt*IJN + pbcInt*PBC + meiInt*MEI
          MCB_1 ~ intMcb*INT + MCB
          
          # Disturbance Covariances
          ATT ~~ IJN + PBC
          IJN ~~ PBC
          
          # Indirect Effects
          ind_mei_att_int := meiATT * attInt
          ind_mei_ijn_int := meiIjn * ijnInt
          ind_mei_pbc_int := meiPbc * pbcInt
          ind_att_int_mcb := attInt * intMcb
          ind_ijn_int_mcb := ijnInt * intMcb
          ind_pbc_int_mcb := pbcInt * intMcb
          ind_mei_int_mcb := meiInt * intMcb
          ind_mei_att_int_mcb := meiATT * attInt * intMcb
          ind_mei_ijn_int_mcb := meiIjn * ijnInt * intMcb
          ind_mei_pbc_int_mcb := meiPbc * pbcInt * intMcb
'
fit_test_sem <- sem(model = test_sem,
                    data = indicators,
                    effect.coding = "loadings", # Coding Effects Method Constraints
                    estimator = "MLR",
                    missing = "FIML")
summary(fit_test_sem, 
        fit.measures=TRUE, 
        standardized= TRUE, 
        rsquare=TRUE)

# Dynamic Fit Index Cutoffs
DDDFI(model = fit_test_sem,
      data = indicators,
      estimator = "MLR")

# Residuals
lavResiduals(fit_test_sem)  # no absolute correlation residuals > .10
# Convert to HTML
res_test_sem <- format(round(lavResiduals(fit_test_sem)$cov,
                             3), 
                       nsmall = 3)  
res_test_sem[upper.tri(res_test_sem)] <- ""
res_test_sem %>% htmlTable()
lavInspect(fit_test_sem, "sampstat")

################################################################################
################################################################################
### R-Squared Difference Calculations ###
################################################################################
################################################################################

source('R2functions.R', chdir = TRUE)

# R-Squared difference between TPB and TPB + MEI with INT as outcome
rsquareCalc(fit_test_sem,
            "INT",
            "MEI",
            adj = TRUE)
rsquareCalc.Boot(test_sem,
                 indicators,
                 "INT",
                 "MEI",
                 adj = TRUE,
                 seed = 1234) # set seed for reproducibility

# R-Squared difference between TPB and TPB + MEI with MCB_1 as outcome
rsquareCalc(fit_test_sem,
            "MCB_1",
            "MEI",
            adj = TRUE)
rsquareCalc.Boot(test_sem,
                 indicators,
                 "MCB_1",
                 "MEI",
                 adj = TRUE,
                 seed = 1234)

# R-Squared difference between TPB + MEI and TPB + MEI + LOC with INT as outcome
rsquareCalc(fit_full_sem,
            "INT",
            "LOC",
            adj = TRUE)
rsquareCalc.Boot(full_sem,
                 indicators,
                 "INT",
                 "LOC",
                 adj = TRUE,
                 seed = 1234)

# R-Squared difference between TPB + MEI and TPB + MEI + LOC with MCB_1 as outcome
rsquareCalc(fit_full_sem,
            "MCB_1",
            "LOC",
            adj = TRUE)
rsquareCalc.Boot(full_sem,
                 indicators,
                 "MCB_1",
                 "LOC",
                 adj = TRUE,
                 seed = 1234)
