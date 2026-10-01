# STGB33 - Labb 3 - Kapitel 8
# Regression models for quantitative and qualitative predictors
# Kutner et al., Applied Linear Statistical Models, §§8.1-8.5
# Graphs in this script are saved as PNG files in the plots/ folder.
dir.create("plots", showWarnings = FALSE)

# ============================================================
# DEL A: POWER CELLS - polynomial regression (§8.1)
# ============================================================
power <- read.csv("power_cells.csv")
head(power)
str(power)
summary(power)

# Centered and scaled variables according to Kutner (8.15)
power$x1 <- (power$CHARGE - 1.0) / 0.4
power$x2 <- (power$TEMP - 20) / 10
power$x1sq <- power$x1^2
power$x2sq <- power$x2^2
power$x1x2 <- power$x1 * power$x2

power

# Why coding helps: compare correlation X with X^2
cor(power$CHARGE, power$CHARGE^2)
cor(power$x1, power$x1^2)
cor(power$TEMP, power$TEMP^2)
cor(power$x2, power$x2^2)

# Second-order model, Kutner (8.13)
power2 <- lm(CYCLES ~ x1 + x2 + x1sq + x2sq + x1x2, data = power)
summary(power2)
anova(power2)

# Figure 8.5: four diagnostic plots for the second-order Power Cells model
png("plots/power-diagnostics.png", width=1500, height=1050, res=160)
par(mfrow = c(2,2), mar=c(4,4,2,1))
plot(fitted(power2), resid(power2), pch=19,
     xlab="Fitted values", ylab="Residuals",
     main="Residuals versus fitted values")
abline(h=0, lty=2, col="grey40")
plot(power$x1, resid(power2), pch=19,
     xlab="x1 (coded charge)", ylab="Residuals",
     main="Residuals versus x1")
abline(h=0, lty=2, col="grey40")
plot(power$x2, resid(power2), pch=19,
     xlab="x2 (coded temperature)", ylab="Residuals",
     main="Residuals versus x2")
abline(h=0, lty=2, col="grey40")
qqnorm(resid(power2), pch=19, main="Normal probability plot")
qqline(resid(power2), col="grey40", lty=2)
dev.off()
par(mfrow = c(1,1))

# Power Cells scatterplot: observed cycles by charge rate and temperature
png("plots/power-observations.png", width=1500, height=850, res=160)
plot(CYCLES ~ CHARGE, data=power,
     pch=c(1,16,17)[match(TEMP, c(10,20,30))],
     col=c("#092235", "#0a8e82", "#587383")[match(TEMP, c(10,20,30))],
     xlab="Charge rate (A)", ylab="Cycles before failure",
     main="Power Cells: observed cycles")
legend("topright", legend=c("10°C", "20°C", "30°C"),
       pch=c(1,16,17), col=c("#092235", "#0a8e82", "#587383"),
       bty="n")
dev.off()

# Step 4A: lack-of-fit test for the second-order Power Cells model
# H0: The model has no lack-of-fit; differences are due to random variation.
# Ha: The model has lack-of-fit at one or more design points.
# 1. Select the three repeated measurements at the center of the design.
center_cycles <- power$CYCLES[power$CHARGE == 1.0 & power$TEMP == 20]

# Pure-error sum of squares: variation among those three measurements.
SSPE <- sum((center_cycles - mean(center_cycles))^2)

# Error sum of squares from the fitted second-order model.
SSE <- sum(resid(power2)^2)

# There are 11 observations, 9 distinct CHARGE/TEMP combinations,
# and 6 estimated model coefficients (including the intercept).
n <- nrow(power)
m <- 9
p <- length(coef(power2))

df_PE <- n - m
df_LOF <- m - p
SS_LOF <- SSE - SSPE

# Lack-of-fit F test and its upper-tail p-value.
F_LOF <- (SS_LOF / df_LOF) / (SSPE / df_PE)
p_value <- pf(F_LOF, df_LOF, df_PE, lower.tail = FALSE)

c(SSPE=SSPE, SS_LOF=SS_LOF, df_PE=df_PE, df_LOF=df_LOF,
  F=F_LOF, p_value=p_value)

# Is a first-order model sufficient?
power1 <- lm(CYCLES ~ x1 + x2, data = power)
anova(power1, power2)
summary(power1)

# Fitted first-order model in original variables, Kutner (8.19)
# x1=(CHARGE-1)/0.4 and x2=(TEMP-20)/10
b <- coef(power1)
b0_original <- b[1] - b[2]*(1/0.4) - b[3]*(20/10)
b_charge <- b[2]/0.4
b_temp <- b[3]/10
c(Intercept=b0_original, CHARGE=b_charge, TEMP=b_temp)

# 90% Bonferroni family confidence intervals for the two slopes
alpha <- 0.10
g <- 2
B <- qt(1 - alpha/(2*g), df=df.residual(power1))
se <- coef(summary(power1))[,"Std. Error"]
ci_charge <- b_charge + c(-1,1)*B*(se["x1"]/0.4)
ci_temp <- b_temp + c(-1,1)*B*(se["x2"]/10)
B
ci_charge
ci_temp

# ============================================================
# DEL B: BODY FAT - interactions among quantitative predictors (§8.2)
# ============================================================
bodyfat <- read.csv("body_fat.csv")

# Center the three predictors exactly as in Kutner's example
bodyfat$x1 <- bodyfat$TRICEPS - mean(bodyfat$TRICEPS)  # mean 25.305
bodyfat$x2 <- bodyfat$THIGH   - mean(bodyfat$THIGH)    # mean 51.170
bodyfat$x3 <- bodyfat$MIDARM  - mean(bodyfat$MIDARM)   # mean 27.620

additive <- lm(BODYFAT ~ x1 + x2 + x3, data=bodyfat)
interact <- lm(BODYFAT ~ x1 + x2 + x3 + x1:x2 + x1:x3 + x2:x3,
               data=bodyfat)
summary(interact)
anova(additive, interact)

# Body Fat: textbook recommends screening residuals from the additive model
# against candidate interaction terms when deciding which interactions to examine.
png("plots/body-fat-interaction-screen.png", width=1650, height=700, res=160)
par(mfrow=c(1,3), mar=c(4,4,3,1))
body.resid <- resid(additive)
plot(bodyfat$x1 * bodyfat$x2, body.resid, pch=19,
     xlab="x1 * x2", ylab="Additive-model residuals",
     main="Residuals versus x1*x2")
abline(h=0, lty=2, col="grey40")
plot(bodyfat$x1 * bodyfat$x3, body.resid, pch=19,
     xlab="x1 * x3", ylab="Additive-model residuals",
     main="Residuals versus x1*x3")
abline(h=0, lty=2, col="grey40")
plot(bodyfat$x2 * bodyfat$x3, body.resid, pch=19,
     xlab="x2 * x3", ylab="Additive-model residuals",
     main="Residuals versus x2*x3")
abline(h=0, lty=2, col="grey40")
dev.off()
par(mfrow=c(1,1))

# Compare correlations before/after centering for selected interaction terms
raw_x1x2 <- bodyfat$TRICEPS * bodyfat$THIGH
ctr_x1x2 <- bodyfat$x1 * bodyfat$x2
cor(bodyfat$TRICEPS, raw_x1x2)
cor(bodyfat$x1, ctr_x1x2)

# ============================================================
# DEL C: INSURANCE INNOVATION - indicator variables (§§8.3-8.5)
# ============================================================
insurance <- read.csv("insurance_innovation.csv")
insurance$TYPE <- factor(insurance$TYPE, levels=c("Mutual","Stock"))
head(insurance)
str(insurance)

# Model without interaction: parallel lines, Kutner (8.33)
parallel <- lm(MONTHS ~ SIZE + TYPE, data=insurance)
summary(parallel)
confint(parallel, "TYPEStock", level=0.95)

# Figure 8.12: data and fitted parallel lines from model (8.33)
png("plots/insurance-parallel-lines.png", width=1500, height=900, res=160)
plot(MONTHS ~ SIZE, data=insurance,
     pch=ifelse(TYPE=="Stock", 19, 1),
     col=ifelse(TYPE=="Stock", "#0a8e82", "#092235"),
     xlab="Size of firm (million dollars)",
     ylab="Months elapsed",
     main="Insurance Innovation: fitted parallel lines")
xx <- seq(min(insurance$SIZE), max(insurance$SIZE), length.out=100)
lines(xx, predict(parallel, newdata=data.frame(SIZE=xx,
      TYPE=factor("Mutual", levels=levels(insurance$TYPE)))), col="#092235")
lines(xx, predict(parallel, newdata=data.frame(SIZE=xx,
      TYPE=factor("Stock", levels=levels(insurance$TYPE)))),
      lty=2, col="#0a8e82")
legend("topright", legend=c("Mutual","Stock"), pch=c(1,19),
       lty=c(1,2), col=c("#092235","#0a8e82"), bty="n")
dev.off()

# Add quantitative x qualitative interaction, Kutner (8.49)
interaction.model <- lm(MONTHS ~ SIZE * TYPE, data=insurance)
summary(interaction.model)
anova(parallel, interaction.model)

# The interaction coefficient tests whether slopes differ.
# If it is unnecessary, retain the simpler parallel-lines model.
