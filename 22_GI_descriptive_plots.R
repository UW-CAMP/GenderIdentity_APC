## Code for descriptive plots, Figures 1, S1, and S2
## of the following manuscript:
##
##    Barry MP, Godwin J, Lucas R, Tordoff DM, Koenig LJ, Goodreau SM. 2025.
##      Demographic trends in gender identity among adults in the United States,
##      2014-2021. International Journal of Transgender Health. 
##      Online first: https://www.tandfonline.com/doi/full/10.1080/26895269.2025.2537874.
## 
##  Script authors: Godwin J
##
##  This script creates 6-panel period-cohort (Figure 1) and
##  age-period (Figure S1) plots for the non cisgender responses to 
##  the BRFSS GI module. An additional 2-panel plot (Figure S2) is
##  created for period-cohort and age-period trends for cisgender respondents
##
##  Inputs:
##     -  data - clean/brfss_final.rds
##     -  data - clean/gi_props_A.rds
##
##  Outputs:
##     -  plots/descriptive/Fig1_GI_pc_all.png 
##     -  plots/descriptive/FigS1_GI_ap_all.png
##     -  plots/descriptive/FigS2_GI_cis.png

# Setup ####
rm(list = ls())

## Folders ####
if(!dir.exists("plots/descriptive/")){
  if(!dir.exists("plots/")){
    dir.create("plots/")
  }
  dir.create("plots/descriptive/")
}

## Functions ####
tableNA <- function(x, ...){
  table(x, useNA = "ifany", ...)  
}

# Load data ####
brfss_dat <- readRDS("data - clean/brfss_final.rds")
gi_props_A <- readRDS("data - clean/gi_props_A.rds")

## Get important objects ####
gi_values <- levels(as.factor(brfss_dat$gender))
sex_values <- levels(as.factor(brfss_dat$sex))
cohort_n <- length(min(brfss_dat$cohort):max(brfss_dat$cohort))
period_n <- length(min(brfss_dat$year):max(brfss_dat$year))
cohort_v <- min(brfss_dat$cohort):max(brfss_dat$cohort) 
period_v <- min(brfss_dat$year):max(brfss_dat$year)

# Plots & Visual Analyses ####

## Plot pars ####
ids <- c("Transgender woman", "Transgender man", "Non-binary/Gender non-conforming", 
         "Cisgender", "Don't know/Not sure",
         "Declined to Answer")

cols <- rainbow_hcl(cohort_n) 
cols_pds <- rainbow_hcl(period_n+1)[1:period_n] 
n_ids <- length(ids)

lwd <- 0.3
cex <- 0.3
loess.lwd <- 1.5
cex.axis <- 0.8
lwd.ticks <- 0.5
cex.axis.labs <- 0.6
cex.main <- 0.8
yline <- 2
xline <- 2
loess.span <- 0.5

## Fig 1: GI period-cohort ####
## Sexual minority identity prevalences by birth cohort & survey year (period)

png(filename = paste0("plots/descriptive/Fig1_GI_pc_all.png"),
    height = 9, width = 6.5, units = "in", res=300)
par(mfrow = c(3,2), lend = 1)
for (id in c(1:3, 5:6)) {
  id_no <- match(id, c(1:3, 5:6))
  plot(x = cohort_v,
       y = gi_props_A[,1,id],
       type = "o",
       col = cols_pds[1],
       xlim = c(cohort_v[1], cohort_v[cohort_n]),
       ylim = (id==4)*c(95, 100) + (id!=4)*c(0,3),
       cex = cex,
       xlab = "", ylab = "",
       xaxs = "i", yaxs = "i",
       axes = FALSE,
       lwd = 0.3)
  
  temp <- loess(gi_props_A[,1,id] ~ cohort_v, span = loess.span)
  lines(temp$x, temp$fitted, col = cols_pds[1], lwd = loess.lwd)
  
  for(period in 2:length(period_v)) {
    lines(x = cohort_v,
          y = gi_props_A[,period,id],
          type = "o", cex = cex,
          col = cols_pds[period], lwd = lwd)
    temp <- loess(gi_props_A[,period,id] ~ cohort_v, span = loess.span)
    lines(temp$x, temp$fitted, col = cols_pds[period], lwd = loess.lwd)
  }
  
  axis(1, col = "grey40", col.axis = "grey20", 
       at = seq(1940,2000,10), cex.axis=cex.axis, lwd=0, lwd.ticks=lwd.ticks,
       mgp=c(3,0.5,0))
  axis(2, col = "grey40", col.axis = "grey20", 
       at = seq(0, 3, 0.25), cex.axis = cex.axis, lwd=0, lwd.ticks=lwd.ticks,
       mgp=c(3,0.5,0))
  box(col = "grey60")
  mtext("Birth cohort", side = 1, outer = FALSE, cex = cex.axis.labs, 
        line = xline, col = "grey20")
  mtext("% of respondents", side = 2, outer = FALSE, cex = cex.axis.labs, 
        line = yline, col = "grey20")
  mtext(paste0(LETTERS[id_no], ") ", ids[id]), side = 3, outer = FALSE, 
        cex = cex.main, col = "grey20", line = 0.5, adj = 0)
  
}

plot(NA, xlim = c(0,1),
     ylim = c(0,1), axes = F, xlab = "", ylab = "")
legend(x = "center",inset = 0,
       border = "black",
       title = "Period (survey year)",
       legend = c(seq(min(period_v), max(period_v))),
       col = cols_pds[seq(1, length(period_v))],
       lwd = 2,
       lty = 1.5,
       cex = 1,
       ncol = 2,
       bty = "n",
       bg = "gray99")
dev.off()

## Fig S1: GI age-period ####
##  Sexual minority identities by age and survey year (period)

png(filename = paste0("plots/descriptive/FigS1_GI_ap_all.png"),
    height = 9, width = 6.5, units = "in", res=300)
par(mfrow = c(3,2), lend = 1)
for (id in c(1:3, 5:6)) {
  id_no <- match(id, c(1:3, 5:6))
  plot(x = rev(period_v[1] - cohort_v),
       y = rev(gi_props_A[,1,id]),
       type = "o",
       col = cols_pds[1],
       xlim= c(18, 80),
       ylim = (id==4)*c(95, 100) + (id!=4)*c(0,3),
       cex = cex,
       xlab = "", ylab = "",
       xaxs = "i", yaxs = "i",
       axes = FALSE,
       lwd = 0.3)
  
  temp <- loess(rev(gi_props_A[,1,id]) ~ rev(period_v[1] - cohort_v), 
                span = loess.span)
  lines(temp$x, temp$fitted, col = cols_pds[1], lwd = loess.lwd)
  
  for(period in 2:length(period_v)) {
    lines(x = rev(period_v[period] - cohort_v),
          y = rev(gi_props_A[,period,id]),
          type = "o", cex = cex,
          col = cols_pds[period], lwd=lwd)
    temp <- loess(rev(gi_props_A[,period,id]) ~ 
                    rev(period_v[period] - cohort_v), span = loess.span)
    lines(temp$x, temp$fitted, col = cols_pds[period], lwd = loess.lwd)
    
  }
  axis(1, col = "grey40", col.axis = "grey20", 
       at = seq(20, 80, 5), cex.axis = cex.axis, lwd = 0, 
       lwd.ticks = lwd.ticks,
       mgp=c(3,0.5,0))
  axis(2, col = "grey40", col.axis = "grey20", 
       at = seq(0, 3, 0.25), cex.axis = cex.axis, lwd = 0, 
       lwd.ticks = lwd.ticks,
       mgp=c(3,0.5,0))
  box(col = "grey60")
  mtext("Age", side = 1, outer = FALSE, cex = cex.axis.labs,
        line = xline, col = "grey20")
  mtext("% of respondents", side = 2, outer = FALSE, cex = cex.axis.labs, 
        line = yline, col = "grey20")
  mtext(paste0(LETTERS[id_no], ") ", ids[id]), side = 3, outer = FALSE, 
        cex = cex.main, col = "grey20", line = 0.5, adj = 0)
}

plot(NA, xlim = c(0,1),
     ylim = c(0,1), axes = F, xlab = "", ylab = "")
legend(x = "center",inset = 0,
       border = "black",
       title = "Period (survey year)",
       legend = c(seq(min(period_v), max(period_v))),
       col = cols_pds[seq(1, length(period_v))],
       lwd = 2,
       lty = 1.5,
       cex = 1,
       ncol = 2,
       bty = "n",
       bg="gray99")
dev.off()


## Fig S2: Cis #### 
## cisgender by cohort-period and age-period
### Pars ####
lwd <- 0.3
cex <- 0.3
loess.lwd <- 1.5
cex.axis <- 0.5
lwd.ticks <- 0.5
cex.axis.labs <- 0.6
cex.main <- 0.8
yline <- 2
xline <- 2
loess.span <- 0.5

png(filename = paste0("plots/descriptive/FigS2_GI_cis.png"),
    height = 4, width = 6.5, units = "in", res=300)
par(mfrow = c(1,2), lend = 1)

id <- 4
plot(x = cohort_v,
     y = gi_props_A[,1,id],
     type = "o",
     col = cols_pds[1],
     xlim = c(cohort_v[1], cohort_v[cohort_n]),
     ylim = (id==4)*c(94, 100) + (id!=4)*c(0,3),
     cex = cex,
     xlab = "", ylab = "",
     xaxs = "i", yaxs = "i",
     axes = FALSE,
     lwd = 0.3
)

temp <- loess(gi_props_A[,1,id] ~ cohort_v, span = loess.span)
lines(temp$x, temp$fitted, col = cols_pds[1], lwd = loess.lwd)

for(period in 2:length(period_v)) {
  lines(x = cohort_v,
        y = gi_props_A[,period,id],
        type = "o", cex=cex,
        col = cols_pds[period], lwd=lwd)
  temp <- loess(gi_props_A[,period,id] ~ cohort_v, span = loess.span)
  lines(temp$x, temp$fitted, col = cols_pds[period], lwd = loess.lwd)
}

axis(1, col = "grey40", col.axis = "grey20", 
     at = seq(1940,2000,10), cex.axis = cex.axis, lwd = 0, 
     lwd.ticks = lwd.ticks, mgp = c(3,0.5,0))
axis(2, col = "grey40", col.axis = "grey20", 
     at = seq(94, 100, 0.5), cex.axis = cex.axis, lwd = 0,
     lwd.ticks = lwd.ticks, mgp = c(3,0.5,0))
box(col = "grey60")
mtext("Birth cohort", side = 1, outer = FALSE, cex = cex.axis.labs, 
      line = xline, col = "grey20")
mtext("% of respondents", side = 2, outer = FALSE, cex = cex.axis.labs, 
      line = yline, col = "grey20")
mtext(paste0(ids[id], ", by cohort"), side = 3, outer = FALSE, 
      cex = cex.main, col = "grey20", line = 0.5, adj = 0)

legend(x = 1955, #"center",
       y = 95.75,
       border = "black",
       title = "Period (survey year)",
       legend = c(seq(min(period_v), max(period_v))),
       col = cols_pds[seq(1, length(period_v))],
       lwd = 2,
       lty = 1.5,
       cex = 0.5,
       ncol = 2,
       bty = "n",
       bg = "gray99")


plot(x= rev(period_v[1] - cohort_v),
     y= rev(gi_props_A[,1,id]),
     type = "o",
     col = cols_pds[1],
     xlim= c(18, 80),
     ylim = (id==4)*c(94, 100) + (id!=4)*c(0,3),
     cex = cex,
     xlab = "", ylab = "",
     xaxs = "i", yaxs = "i",
     axes = FALSE,
     lwd = 0.3
)

temp <- loess(rev(gi_props_A[,1,id]) ~ rev(period_v[1] - cohort_v),
              span= loess.span)
lines(temp$x, temp$fitted, col = cols_pds[1], lwd = loess.lwd)

for(period in 2:length(period_v)) {
  lines(x= rev(period_v[period] - cohort_v),
        y= rev(gi_props_A[,period,id]),
        type = "o", cex=cex,
        col = cols_pds[period], lwd=lwd)
  temp <- loess(rev(gi_props_A[,period,id]) ~ 
                  rev(period_v[period] - cohort_v), span = loess.span)
  lines(temp$x, temp$fitted, col = cols_pds[period], lwd = loess.lwd)
  
}

axis(1, col = "grey40", col.axis = "grey20", 
     at = seq(20, 80, 5), cex.axis = cex.axis, lwd = 0,
     lwd.ticks = lwd.ticks, mgp = c(3,0.5,0))
axis(2, col = "grey40", col.axis = "grey20", 
     at = seq(94, 100, 0.5), cex.axis = cex.axis, lwd = 0,
     lwd.ticks = lwd.ticks, mgp = c(3,0.5,0))
box(col = "grey60")
mtext("Age", side = 1, outer = FALSE, cex = cex.axis.labs, 
      line = xline, col = "grey20")
mtext("% of respondents", side = 2, outer = FALSE, cex = cex.axis.labs, 
      line = yline, col = "grey20")
mtext(paste0(ids[id], ", by age"), side = 3, outer = FALSE, 
      cex = cex.main, col = "grey20", line = 0.5, adj = 0)

dev.off()