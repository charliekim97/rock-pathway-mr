## Task 9 — Power / minimum detectable effect (null-paper shield, tone DOWN).
## Reports, under stated assumptions, the smallest effect the ROCK2 instrument could detect at 80% power.
## NOT a claim that larger effects are excluded with certainty — assumptions are listed in the legend.
source("R/config.R"); source("R/helpers.R"); suppressMessages(library(ggplot2))

## instrument variance explained (R2) for ROCK2, summed over the r2<0.01 set (approx, from frozen betas/eaf)
inst <- load_instruments()[gene=="ROCK2" & tier %in% c("ROCK2_blood_lead","ROCK2_blood_r2.01")]
R2_x <- sum(2*inst$eaf*(1-inst$eaf)*inst$beta^2)          # SD-units expression; approximate
F_x  <- 409                                               # lead F (reported)

## continuous outcome MDE (per SD): se(beta_MR) ~ 1/sqrt(N * R2_x); MDE = (z_a/2 + z_b)/sqrt(N*R2_x)
mde_cont <- function(N, R2=R2_x, alpha=0.05, power=0.8)
  (qnorm(1-alpha/2)+qnorm(power))/sqrt(N*R2)

## binary outcome: effective N = 4 / (1/ncase + 1/nctrl); MDE on log-OR scale
mde_bin <- function(ncase,nctrl,R2=R2_x,alpha=0.05,power=0.8){
  Neff <- 4/(1/ncase + 1/nctrl); (qnorm(1-alpha/2)+qnorm(power))/sqrt(Neff*R2) }

tab <- rbind(
  data.table(outcome="BMI (cont, N=806834)",        MDE=mde_cont(806834),        scale="SD"),
  data.table(outcome="Hepatic fat (cont, N=32974)", MDE=mde_cont(32974),         scale="SD"),
  data.table(outcome="NAFLD (8434/770180)",         MDE=mde_bin(8434,770180),    scale="logOR"),
  data.table(outcome="Cirrhosis (5545/494803)",     MDE=mde_bin(5545,494803),    scale="logOR"),
  data.table(outcome="T2D DIAMANTE (80154/853816)", MDE=mde_bin(80154,853816),   scale="logOR")
)
tab[scale=="logOR", OR_MDE := exp(MDE)]
fwrite(tab, file.path(CONFIG$out,"task9_MDE_table.csv")); print(tab)

## power curves
grid <- seq(0,0.6,0.01)
pw <- function(eff,N,R2=R2_x,alpha=0.05) pnorm(eff*sqrt(N*R2)-qnorm(1-alpha/2))
curves <- rbindlist(lapply(list(c("BMI",806834),c("HepaticFat",32974)), function(x)
  data.table(outcome=x[1], effect=grid, power=pw(grid, as.numeric(x[2])))))
g <- ggplot(curves, aes(effect,power,color=outcome))+geom_line(linewidth=1)+
  geom_hline(yintercept=0.8,linetype=2)+theme_bw(base_size=11)+
  labs(x="True effect (per-SD expression)", y="Power",
       title="Minimum detectable effect (continuous outcomes)",
       caption=paste0("Assumptions: instrument R2=",signif(R2_x,3),
       " (ROCK2 r2<0.01 set), alpha=0.05, two-sided; binary outcomes use effective N=4/(1/ncase+1/nctrl);",
       " log-OR converted to OR. Power is a shield, not a sharp exclusion."))
ggsave(file.path(CONFIG$out,"task9_power_curves.pdf"), g, width=7, height=4.2)
ggsave(file.path(CONFIG$out,"task9_power_curves.png"), g, width=7, height=4.2, dpi=300)
cat("\nLegend MUST state: R2 derivation, effective N, alpha, effect scale, OR conversion.\n")
