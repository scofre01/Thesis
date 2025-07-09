#Script for Meta-analysis

##1 Artículo análisis de confiabilidad de métodos basados en inteligencia artificial
#Instalar paquete
install.packages("meta")
library(meta)
library(metafor)
calories_cor$r <-as.numeric(calories_cor$r) #verificar que el R este como estadistico numerico

### análisis de correlaciones para calorias
cal <-metacor(data=calories_cor, r_energy, total, sm="COR" , random=T, fixed=F, studlab = ID) #Permite calcular el meta-analisis
cal
forest(cal, layout = "RevMan5") ## Para hacer el forest plot

eg_test_cal <- regtest(cal, model = "rma", predictor = "sei")

### análisis de correlaciones para cho
cho <-metacor(data=calories_cor, r_cho, total, sm="COR" , random=T, fixed=F, studlab = ID) #Permite calcular el meta-analisis
cho
forest(cho, layout = "RevMan5") ## Para hacer el forest plot

### análisis de correlaciones para prot
prot <-metacor(data=calories_cor, r_prot, total, sm="COR" , random=T, fixed=F, studlab = ID) #Permite calcular el meta-analisis
prot
forest(prot, layout = "RevMan5") ## 

### análisis de correlaciones para fat
fat <-metacor(data=calories_cor, r_fat, total, sm="COR" , random=T, fixed=F, studlab = ID) #Permite calcular el meta-analisis
fat
forest(fat, layout = "RevMan5") ## 


## 2. Artículo 2 meta-analisis UPF vs CH.
library(meta)
library(metadat)
library(metafor)
attach(metaanalisis)
settings.meta('RevMan5')
m<-metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = author,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic",
  title= "UPFs and Liver Cancer"## method to calculate estimator's CI
)
summary(m)
## Forest plot of the estimates
forest(m,xlim = c(0.3,4.5))
m1<-forest(m,main = "Riesgo de Alto Consumo de UP vs Bajo Consumo")


## Funnel plot of the meta-analysis
funnel(m, studlab = TRUE)

## test de egger
sesgo_publicacion <- metabias(m,  method.bias = "rank")

### Meta analysis of cumulative incidence
m3<- metaprop(cases, total, studlab= author, sm="PLOGIT", data=metaanalisis, method="GLMM",  
              fixed = FALSE,
              random = TRUE,
              method.random.ci = "HK",
              title = "Cumulative incidence of HCC")
summary(m3)

forest(m3, layout="RevMan5", xlab="Proportion", comb.r=T, comb.f=F, xlim = c(0,1), fontsize=10, digits=3)

meta::forest(m3, layout = "JAMA")

m2 <- metaprop(event = cases,
                       n = total,
                       studylab = author, # Etiquetas para los estudios (opcional)
                       sm = "PFT", # Especifica que la medida de efecto es la proporción transformada (útil para estabilizar la varianza)
                       method = "Inverse", # Especifica el método para estimar la varianza entre estudios (REML es recomendado)
                       random = TRUE, # Indica que se utilizará un modelo de efectos aleatorios
                       hakn = TRUE) # Aplica la corrección de Hartung-Knapp para los intervalos de confianza (recomendado para pocos estudios)

# 4. Resumen de los resultados del Meta-Análisis
summary(m2)

# También puedes ver una representación gráfica del meta-análisis (forest plot)
forest(m2)


#########################################
#Meta analisis dose - response 
library(dplyr)
library(dosresmeta)
dsma <- dsma %>%
  mutate(
    log_rr = log(hr),
    se_log_rr = (log(ci_upper) - log(ci_lower)) / (2 *qnorm(0.975)))

###Graficar
library(ggplot2)
dsma$inver_se <- 1/dsma$se_log_rr

ggplot(dsma, aes(dose, log_rr, size=inver_se)) + geom_point(shape=1, colour= "black") + scale_size_area (max_size=20)

### Modelo lineal
lin_bin <- dosresmeta(formula=log_rr ~ dose, 
                      id = id, 
                      type= type,
                      se = se_log_rr, 
                      cases= cases,
                      n= n,
                      data = dsma)

summary(lin_bin)
predict(lin_bin, delta=1, exp =TRUE)

dosex_bin <- data.frame(dose=seq(0, 100, 1))
with(predict(lin_bin, dosex_bin, order=TRUE, exp=TRUE), {plot(dose, pred, type="l", col="blue", ylim=c(0, 2), ylab= "liver cancer relative risk", xlab= "UPF consumption, %/day")
  lines(dose, ci.lb, lty=2)
  lines(dose, ci.ub, lty=2)})


### cubic splines
library("rms")
k<-with(dsma, quantile(dose, probs=c(0.1, 0.5, 0.9), na.rm=TRUE)) #knots
spl <- dosresmeta(log_rr ~ rcs(dose, k), id = id, se = se_log_rr, type = type, cases = cases, n = n, data = dsma, method="mm")
summary(spl)

xref_bin <- 0
with(predict(spl, dosex_bin, xref_bin, exp=TRUE),{plot(get("rcs(dose, k)dose"), pred, type="l", ylim=c(0.4, 10), ylab="liver cancer relative risk", xlab="UPF consumption, %/day", log="y", bty="l", las=1) matlines (get("rcs(dose, k)dose"), cbind(ci.ub, ci.lb), col=1, lty="dashed")})
points(dosex_bin$dose, predict(lin_bin, dosex_bin, xref_bin, exp=TRUE)$pred, type="l", lty=3, col="blue")

