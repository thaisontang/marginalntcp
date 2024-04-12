res.dir <- here("simulated_illustration", "figures", "Figure 1")

alph <- 0.12
color_labs <- c("blue", "orange")

# Matrices of DVHs
blad.dvh <- long.data %>% dplyr::select(ID, Dose, Volume_Bladder) %>%
  pivot_wider(names_from = Dose, values_from = Volume_Bladder) %>%
  dplyr::select(-ID) %>% data.matrix
skin.dvh <- long.data %>% dplyr::select(ID, Dose, Volume_Skin) %>%
  pivot_wider(names_from = Dose, values_from = Volume_Skin) %>%
  dplyr::select(-ID) %>% data.matrix

# DVHs Stratified by Age at Diagnosis
strat.var <- "Age65_numeric"
strat.vals <- levels(wide.data$Age65)
name_labs <- paste0(strat.vals, " (N = ", table(wide.data[,strat.var]), "; ", round(prop.table(table(wide.data[,strat.var]))*100, 0), "%)")


# FIGURE 1A: BLADDER DVHs (by Age65+)
pdf(here(res.dir, "Figure1a.pdf"), width = 8, height =5)
op <- par(las = 1, mar = c(4.5, 4.5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
sub.data <- long.data[which(long.data$ID == 1),]
plot(sub.data$Dose, sub.data$Volume_Bladder, type = "l", col = "white",
     ylab = "Volume (%)", xlab = "Dose (Gy)")
for(id in wide.data$ID) {
  sub.data <- long.data[which(long.data$ID == id),]
  lines(sub.data$Dose, sub.data$Volume_Bladder, col = ggplot2::alpha(color_labs, alph)[pull(wide.data[id,strat.var]) + 1])
}
lines(dvals, colMeans(blad.dvh[which(wide.data$Age65_numeric == 1),]), lwd = 3, col = "orange")
lines(dvals, colMeans(blad.dvh[which(wide.data$Age65_numeric == 0),]), lwd = 3, col = "blue")
legend("topright", col = color_labs, name_labs, lty = rep(1, 2), lwd = rep(2, 2))
par(op)
dev.off()

# FIGURE 1B: SKIN DVHs (by Age65+)
pdf(here(res.dir, "Figure1b.pdf"), width = 8, height =5)
op <- par(las = 1, mar = c(4.5, 4.5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
sub.data <- long.data[which(long.data$ID == 1),]
plot(sub.data$Dose, sub.data$Volume_Skin, type = "l", col = "white",
     ylab = "Volume (%)", xlab = "Dose (Gy)")
for(id in wide.data$ID) {
  sub.data <- long.data[which(long.data$ID == id),]
  lines(sub.data$Dose, sub.data$Volume_Skin, col = ggplot2::alpha(color_labs, alph)[pull(wide.data[id,strat.var]) + 1])
}
lines(dvals, colMeans(skin.dvh[which(wide.data$Age65_numeric == 1),]), lwd = 3, col = "orange")
lines(dvals, colMeans(skin.dvh[which(wide.data$Age65_numeric == 0),]), lwd = 3, col = "blue")
legend("topright", col = color_labs, name_labs, lty = rep(1, 2), lwd = rep(2, 2))
par(op)
dev.off()


# DVHs Stratified by Toxicities

# FIGURE 1C: BLADDER DVHs (by Genitourinary (GU) Toxicity)
strat.var <- "GU Toxicity"
strat.vals <- c("No", "Yes")
name_labs <- paste0(strat.vals, " (N = ", table(wide.data[,strat.var]), "; ", round(prop.table(table(wide.data[,strat.var]))*100, 0), "%)")

pdf(here(res.dir, "Figure1c.pdf"), width = 8, height =5)
op <- par(las = 1, mar = c(4.5, 4.5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
sub.data <- long.data[which(long.data$ID == 1),]
plot(sub.data$Dose, sub.data$Volume_Bladder, type = "l", col = "white",
     ylab = "Volume (%)", xlab = "Dose (Gy)")
for(id in wide.data$ID) {
  sub.data <- long.data[which(long.data$ID == id),]
  lines(sub.data$Dose, sub.data$Volume_Bladder, col = ggplot2::alpha(color_labs, alph)[pull(wide.data[id,strat.var]) + 1])
}
lines(dvals, colMeans(blad.dvh[which(wide.data$Age65_numeric == 1),]), lwd = 3, col = "orange")
lines(dvals, colMeans(blad.dvh[which(wide.data$Age65_numeric == 0),]), lwd = 3, col = "blue")
legend("topright", col = color_labs, name_labs, lty = rep(1, 2), lwd = rep(2, 2))
par(op)
dev.off()

# FIGURE 1D: SKIN DVHs (by Gastrointestinal (GI) Toxicity)
strat.var <- "GI Toxicity"
strat.vals <- c("No", "Yes")
name_labs <- paste0(strat.vals, " (N = ", table(wide.data[,strat.var]), "; ", round(prop.table(table(wide.data[,strat.var]))*100, 0), "%)")

pdf(here(res.dir, "Figure1d.pdf"), width = 8, height =5)
op <- par(las = 1, mar = c(4.5, 4.5, 1.0, 1.0), oma = c(0, 0, 0, 0), mgp = c(3, 1, 0))
sub.data <- long.data[which(long.data$ID == 1),]
plot(sub.data$Dose, sub.data$Volume_Skin, type = "l", col = "white",
     ylab = "Volume (%)", xlab = "Dose (Gy)")
for(id in wide.data$ID) {
  sub.data <- long.data[which(long.data$ID == id),]
  lines(sub.data$Dose, sub.data$Volume_Skin, col = ggplot2::alpha(color_labs, alph)[pull(wide.data[id,strat.var]) + 1])
}
lines(dvals, colMeans(skin.dvh[which(wide.data$Age65_numeric == 1),]), lwd = 3, col = "orange")
lines(dvals, colMeans(skin.dvh[which(wide.data$Age65_numeric == 0),]), lwd = 3, col = "blue")
legend("topright", col = color_labs, name_labs, lty = rep(1, 2), lwd = rep(2, 2))
par(op)
dev.off()

