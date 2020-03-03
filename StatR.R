cat("\014") 
rm(list = ls())
library('lme4')
library('nlme')

FileName = '/Users/elhamb/Documents/Codes/Git/LetterVernier-SSVEP/ResultData/Letter_forRANOVA.txt'
T1 = read.table(FileName,header = TRUE)

#df[col_names] <- lapply(df[col_names] , factor)

lmeModel = lmer(Score ~ RC*logMAR + (1|ids), data=T1)
lmeModel2 = lme(Score ~ RC*logMAR , random=~1 | ids, data=T1)
A = anova(lmeModel2)


AOV = aov(Score ~ RC*logMAR + (1|ids), data=T1)
summary(AOV)


# the MANOVA for SNR comparison
#factor1 = S/N
#dependent variables = 5 different logMAR
# We do separately for letter and vernier, what about RC1 and RC2?

