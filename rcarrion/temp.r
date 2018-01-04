library(ggplot2)
size <- c(240, 960, 2160, 4000, 14400)
cpu1 <- cbind(size, "cpu1", c(0.1582022315938957,1.9797110311648187,9.435577033172983,30.4298446176302,406.666469568567))
cpu24 <- cbind(size, "cpu24", c(0.028694471441364538,0.19484387850388885,0.8395908067235723,2.368085369638478,26.41598431600723))
cpu48 <- cbind(size, "cpu48", c(0.032188582773475596,0.15417953062181672,0.6157196523000796,1.7285729438493338,20.5429633437190))
gpu8 <- cbind(size, "gpu8", c(2.5673562765121463,2.7172582308451334,3.192317875226339,4.595608019828797,32.2447277069091))
DP <- c(0.011091450184640557,0.07538209549066997,0.24890434449030704,1.1237477278976828,21.07282838271291,0.0058534082505159275,0.007686879825102305,0.026293040203117207,0.09115285716522209,1.2297147445218948,0.009643587515246639,0.00873484003978532,0.03255008131287857,0.03862431941229486,0.2595471596150877,0.0004160932314931688,0.0002170553367192258,0.007127573675714469,0.007996297398413028,0.12217453281608205)
DF <- data.frame(rbind(cpu1,cpu24,cpu48,gpu8))
names(DF) <- c("Size", "Version", "Mean")
DF$Size <- factor(DF$Size, levels = c("240", "960", "2160", "4000", "14400"))
DF$Mean <- as.numeric(as.character(DF$Mean))
DP <- as.numeric(as.character(DP))

Graph <- ggplot(DF, aes(x = Size, y = Mean, color = Version, group = Version)) + 
  geom_point(size = 2.5,mapping = aes(shape=Version)) + 
  geom_line(aes(linetype = Version), size=1.25) +
  theme_bw() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(sides = "l")  +
  geom_errorbar(aes(ymax = Mean + 1.96*DP/sqrt(30), ymin = Mean - 1.96*DP/sqrt(30)), width=0.4) +
  ylab("Mean (seconds)") + 
  xlab(expression(paste("Mesh size"))) +
  theme(plot.title = element_text(family = "Times", face="bold", size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(family = "Times", face="bold", size=12)) +
  theme(axis.text  = element_text(family = "Times", face="bold", size=10, colour = "Black")) +
  theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
  theme(legend.text  = element_text(family = "Times", face="bold", size=12)) +
  theme(legend.key.size = unit(5, "cm")) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key=element_rect(size=0),
        legend.key.size = unit(1, "lines")) +
  guides(col = guide_legend(nrow = 1))
ggsave(paste("~/Giuliano-plots.pdf",sep=""), Graph,  height=10, width=15, units="cm")
