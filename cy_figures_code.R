##### FIGURE 1: DIET AND HEALTH RESPONSES ####
diet.figure <- cy.metadata[,c("Location","d13C","d15N","Vol.Anthro","VOL.PREY")]

# (a) Stable isotope isoscape ####
diet.isoscape <- diet.figure %>%
  dplyr::group_by(Location) %>%
  dplyr::summarise(meanC=mean(d13C),stdC=sd(d13C),meanN=mean(d15N),stdN=sd(d15N))

diet.figure.a <- ggplot() +
  geom_point(data=diet.figure, aes(x=d13C, y=d15N, color=Location, shape=Location), size=1.5) + 
  ggforce::geom_ellipse(data=diet.isoscape, aes(x0=meanC, y0=meanN, a=stdC, b=stdN, angle=0, color=Location, linetype=Location), size=0.75) +
  geom_path() +
  theme_bw() +
  labs(x=expression(paste(delta^{13}, "C (\u2030)")), y=expression(paste(delta^{15}, "N (\u2030)")), tag="a") +
  guides(color=guide_legend(title="Location"), linetype=guide_legend(title="Location"), shape=guide_legend(title="Location")) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_y_continuous(limits=c(7,11)) +
  scale_x_continuous(limits=c(-24.5,-18.9), breaks=c(-24, -23, -22, -21, -20)) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=11),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        legend.position="none")

# (b) Stable isotope mixing model results ####
simmr_plot <- summary(simmr_groups_indv_out, type=c('statistics'), group=c(1:2))$statistics
simmr_plot <- cbind(simmr_plot[[1]],
                    simmr_plot[[2]])
colnames(simmr_plot) <- c("rural_mean", "rural_sd", "urban_mean", "urban_sd")
simmr_plot <- subset(simmr_plot, rownames(simmr_plot) %in% c("anthro","prey","fruit"))
simmr_plot <- 100*simmr_plot

simmr_plot <- as.data.frame(simmr_plot)
simmr_plot$Source <- rownames(simmr_plot)
simmr_plot <- reshape2::melt(simmr_plot)
simmr_plot <- tidyr::separate(simmr_plot, variable, into=c("Location", "Measure"),
                                   sep = "_", remove = TRUE)
temp1 <- subset(simmr_plot, Measure=="mean")
temp2 <- subset(simmr_plot[c("Measure","value")], Measure=="sd")
colnames(temp2) <- c("Measure2", "sd")
simmr_plot <- cbind(temp1, temp2)
rm(temp1, temp2)
simmr_plot$Source <- as.factor(simmr_plot$Source)
levels(simmr_plot$Source)[levels(simmr_plot$Source)=="anthro"] <- "Anthro\nfood"
levels(simmr_plot$Source)[levels(simmr_plot$Source)=="fruit"] <- "Fruit"
levels(simmr_plot$Source)[levels(simmr_plot$Source)=="prey"] <- "Prey"

plot.simmr <- ggplot(simmr_plot, aes(x=Source, y=value, fill=Location)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=value-sd, ymax=value+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() +
  scale_fill_manual(values=c("forestgreen", "purple")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,80)) +
  theme(axis.text.y=element_text(size=10, color="black"),
        axis.text.x=element_text(size=10, color="black"),
        axis.title=element_text(size=10, color="black"),
        legend.text=element_text(size=11, color="black"),
        legend.position="none",
        plot.tag=element_text(face="bold"),
        axis.line = element_line(colour = "black", size=0.5),
        panel.border=element_blank(),
        panel.grid=element_blank()) +
  labs(x="Diet item", y="Proportion (%)", tag="c")

# (c) Stomach contents ####
diet.figure.data <- melt(diet.figure)
levels(diet.figure.data$variable)[levels(diet.figure.data$variable)=="Vol.Anthro"] <- "Anthropogenic\nFood"
levels(diet.figure.data$variable)[levels(diet.figure.data$variable)=="VOL.PREY"] <- "Prey"
levels(diet.figure.data$Location)[levels(diet.figure.data$Location)=="Rural"] <- "rural"
diet.figure.data <- subset(diet.figure.data, variable !="d13C" & variable !="d15N")
diet.figure.data$variable <- factor(diet.figure.data$variable)

diet.figure.b <- ggplot(diet.figure.data, aes(x=Location, y=log(value), fill=Location)) + geom_boxplot() +
  facet_wrap(~variable) +
  scale_fill_manual(values=c("forestgreen", "purple")) +
  labs(x="\nLocation", y="log Volume (mL)", tag="b") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        axis.text.x=element_blank(), 
        strip.background=element_blank(),
        strip.text=element_text(color="black", size=10),
        legend.position="none",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

# (d) Health measures ####
health.figure.data <- cy.metadata[,c("Location", "KFI", "Spleen.to.Body.Ratio","index3.PC1", "index3.PC1.z", "Em.PCR.Overall","mass.size.residual",
                                     "mass.size.residual.z")]
health.figure.data <- melt(health.figure.data)
levels(health.figure.data$variable)[levels(health.figure.data$variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass (g/kg)"
levels(health.figure.data$variable)[levels(health.figure.data$variable)=="index3.PC1"] <- "Health score"
levels(health.figure.data$Location)[levels(health.figure.data$Location)=="peri-urban"] <- "rural"

health.plot <- ggplot(subset(health.figure.data, variable %in% c("KFI", "Spleen mass\n(g/kg)", "Health score")),
       aes(x=Location, y=value, fill=Location)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~variable, scales="free") +
  scale_fill_manual(values=c("forestgreen", "purple")) +
  labs(x="\nLocation", y="Value", tag="d") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=10),
        axis.text.x=element_blank(), 
        strip.background=element_blank(),
        strip.text=element_text(color="black", size=10),
        legend.position="none",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

# (e) E. multilocularis prevalence ####
cy.metadata %>%
  dplyr::group_by(Location) %>%
  dplyr::count(Em.PCR.Overall)

em.table <- data.frame(Location=c("rural", "urban"),
                       Em.Perc=c(35.38462,57.14286))

em.plot <- ggplot(em.table, aes(x=Location, y=Em.Perc, fill=Location)) + geom_bar(stat="identity",position="dodge") +
  scale_fill_manual(values=c("forestgreen", "purple")) +
  labs(x="\nLocation", y="E. multi prevalence (%)", tag="e") +
  theme_bw() +
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=10),
        axis.text.x=element_blank(), 
        strip.background=element_blank(),
        strip.text=element_text(color="black", size=10),
        legend.position="none",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

# - Figure design ####
grid.arrange(grobs=list(diet.figure.a, diet.figure.b, plot.simmr, health.plot, em.plot),
             layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),c(NA,NA,4,4,4,4,4,4,4,5,5,5,NA,NA)))


Figure_1_legend <- cowplot::get_legend(ggplot(em.table, aes(x=Location, y=Em.Perc, fill=Location)) + geom_bar(stat="identity",position="dodge") +
                                         scale_fill_manual(values=c("forestgreen", "purple")) +
                                         labs(x="\nLocation", y="E. multi prevalence (%)", tag="e") +
                                         theme_bw() +
                                         scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
                                         theme(panel.grid=element_blank(),
                                               axis.text.y=element_text(color="black", size=10),
                                               axis.title=element_text(color="black", size=11),
                                               axis.text.x=element_blank(), 
                                               strip.background=element_blank(),
                                               strip.text=element_text(color="black", size=10),
                                               legend.position="right",
                                               legend.title=element_text(color="black", size=10),
                                               legend.text=element_text(color="black", size=9)))
             

grid.arrange(grobs=list(diet.figure.a, plot.simmr+labs(tag="b"), diet.figure.b+labs(tag="c"), health.plot, em.plot, Figure_1_legend),
             layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),c(NA,4,4,4,4,4,4,4,5,5,5,NA,6,NA)))

Figure_1_diet_health <- arrangeGrob(grobs=list(diet.figure.a, 
                                               plot.simmr+labs(tag="b"), 
                                               diet.figure.b+labs(tag="c"),
                                               health.plot, em.plot, Figure_1_legend),
                                    layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),c(NA,4,4,4,4,4,4,4,5,5,5,NA,6,NA)))
ggsave("~/health/Tables/Figure_1_diet_health.jpg", Figure_1_diet_health, dpi=300, width=6.5, height=5, units="in")


#### FIGURE 2: MICROBIOME STRUCTURE & COMPOSITION ####
# (a) Urban/rural alpha diversity boxplots ####
adiv.urb.rur <- picante.data[,c("Location","Observed","Shannon","PD","NTI")]
adiv.urb.rur <- reshape2::melt(adiv.urb.rur)

adiv.urb.rur.plot <- ggplot(adiv.urb.rur, aes(x=Location, y=value, fill=Location)) + 
  geom_boxplot(outlier.shape=NA) + 
  facet_wrap(~variable, scales="free", nrow=1) +
  scale_fill_manual(values=c("forestgreen", "purple")) +
  labs(x=NULL, y="Value", tag="a") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=9),
        axis.text.x=element_blank(), 
        strip.background=element_blank(),
        panel.spacing=unit(0.5, "lines"),
        strip.text=element_text(color="black", size=9),
        legend.position="none",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

# (b) Alpha diversity correlations with other measures ####
cor.feces <- spearman.correlations(picante.data)

for(i in 1:length(cor.feces)){
  cor.feces[[i]] <- subset(cor.feces[[i]], colnames(cor.feces[[i]]) %in% c("Mass", "Snout.to.Base", "Girth", "BMI", "Cementum.Age", "Spleen.to.Body.Ratio", "KFI", "d13C", "d15N"))
  cor.feces[[i]] <- as.data.frame(t(cor.feces[[i]]))
  cor.feces[[i]] <- subset(cor.feces[[i]], rownames(cor.feces[[i]]) %in% c("Observed", "Shannon", "PD", "NTI"))
  cor.feces[[i]]$Measure <- rownames(cor.feces[[i]])
}

for(i in 1:length(cor.feces)){
  cor.feces[[i]] <- subset(cor.feces[[i]], colnames(cor.feces[[i]]) %in% c("index3.PC1", "Cementum.Age", "Spleen.to.Body.Ratio", "d13C", "d15N", "Vol.Anthro","VOL.PREY"))
  cor.feces[[i]] <- as.data.frame(t(cor.feces[[i]]))
  cor.feces[[i]] <- subset(cor.feces[[i]], rownames(cor.feces[[i]]) %in% c("Observed", "Shannon", "PD"))
  cor.feces[[i]]$Measure <- rownames(cor.feces[[i]])
}


cor.feces.R.melt <- reshape2::melt(cor.feces[["r"]])
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="index3.PC1"] <- "Health index"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="Cementum.Age"] <- "Age"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="Vol.Anthro"] <- "Vol. anthro"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="VOL.PREY"] <- "Vol. prey"
cor.feces.R.melt$Measure[cor.feces.R.melt$Measure=="Observed"] <- "Richness"
cor.feces.R.melt$Measure <- factor(cor.feces.R.melt$Measure, levels=c("PD", "Shannon", "Richness"))

cor.feces.R.melt$variable <- factor(cor.feces.R.melt$variable, levels = c(levels(cor.feces.R.melt$variable), " "))
desired_Order <- c("Age", "Health index", "Spleen mass", "d13C", "Vol. anthro", "d15N", "Vol. prey")
cor.feces.R.melt$variable <- factor(as.character(cor.feces.R.melt$variable), levels=desired_Order)

Fig2_b_legend <- cowplot::get_legend(ggplot(subset(cor.feces.R.melt, variable !="BMI"), aes(x=variable, y=value)) + 
  geom_col(aes(fill=Measure), position="dodge", width=0.75) +
  #scale_x_discrete(labels=function(variable) str_wrap(variable, width=5)) +
  #scale_fill_manual(values=c("goldenrod", "deepskyblue3")) +
  scale_x_discrete(labels=c("Age", "Health score", "Spleen mass", expression(paste(delta^{13}, "C")),
                            "Vol. anthro", expression(paste(delta^{15}, "N")), "Vol. prey")) +
  guides(fill=guide_legend(reverse=TRUE, ncol=1)) +
  labs(x="\nMeasure", y="Spearman's R\n") +
  coord_flip() +
  labs(tag="a") +
  theme_bw() +
  theme(axis.title=element_text(color="black", size=10),
        axis.text.y=element_text(color="black", size=9, hjust=1),
        axis.text.x=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(color="black", size=9),
        panel.border=element_rect(color="black"),
        panel.grid=element_blank(),
        legend.key.size=unit(1,"line"),
        plot.title=element_text(size=10),
        panel.grid.minor=element_blank()) +
  ylim(c(-0.31,0.31)) +
  geom_hline(yintercept=0, color="black"))

Fig.Adiv.legend <- g_legend(Fig.Adiv)
rm(Fig.Adiv)

Fig.Adiv <- ggplot(subset(cor.feces.R.melt, variable !="BMI"), aes(x=variable, y=value)) + 
  geom_col(aes(fill=Measure), position="dodge", width=0.75) +
  #scale_x_discrete(labels=function(variable) str_wrap(variable, width=5)) +
  #scale_fill_manual(values=c("goldenrod", "deepskyblue3")) +
  scale_x_discrete(labels=c("Age", "Health score", "Spleen mass", expression(paste(delta^{13}, "C")),
                            "Vol. anthro", expression(paste(delta^{15}, "N")), "Vol. prey")) +
  guides(fill=guide_legend(reverse=TRUE, ncol=1)) +
  labs(x="\nMeasure", y="Spearman's R\n") +
  coord_flip() +
  labs(tag="a") +
  theme_bw() +
  theme(axis.title=element_text(color="black", size=9),
        axis.text.y=element_text(color="black", size=9, hjust=1),
        axis.text.x=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        legend.position="none",
        legend.title=element_blank(),
        legend.text=element_text(color="black", size=9),
        panel.border=element_rect(color="black"),
        panel.grid=element_blank(),
        legend.key.size=unit(1,"line"),
        plot.title=element_text(size=10),
        panel.grid.minor=element_blank()) +
  ylim(c(-0.31,0.31)) +
  geom_hline(yintercept=0, color="black")

# (c) Relative abundance bar charts ####
y1 <- tax_glom(pseq.rarefied, taxrank="Class")
y2 <- merge_samples(y1, "Location")

y3 <- as.data.frame(sample_data(y2))
y3[,10:94] <- NULL
y3$Location <- as.character(y3$Location)
y3$Location[y3$Location=="1"] <- "Rural"
y3$Location[y3$Location=="2"] <- "Urban"
sample_data(y2) <- y3

y2 <- transform_sample_counts(y2, function(x) 100*x/sum(x))

y4 <- psmelt(y2)
y4$Phylum <- as.character(y4$Phylum)
y4$Phylum[y4$Phylum %in% c("Acidobacteria","Candidatus_Saccharibacteria","Cyanobacteria/Chloroplast",
                           "Deferribacteres","Deinococcus-Thermus","Euryarchaeota","Fibrobacteres","Nitrospirae",
                           "Planctomycetes","Spirochaetes","Tenericutes","Thaumarchaeota","Verrucomicrobia")] <- "Other"
y4$Class <- as.character(y4$Class)
y4$Class[y4$Phylum=="Other"] <- "Taxa <1% Abundance"
y4$Class[y4$Phylum=="Proteobacteria"] <- "Proteobacteria"
y4$Class[y4$Class %in% c("Bacteroidia", "Flavobacteriia", "Sphingobacteriia", "Other Bacteroidetes")] <- "Bacteroidetes"

y4$Phylum <- as.factor(y4$Phylum)
y4$Class <- as.factor(y4$Class)

desired_order <- c("Other","Actinobacteria","Bacteroidetes","Proteobacteria","Fusobacteria","Firmicutes")
y4$Phylum <- factor(as.character(y4$Phylum), levels=desired_order)
y4 <- y4[order(y4$Phylum, y4$Class),]
colorCount <- length(unique(y4$Class))

desired_order <- c("Taxa <1% Abundance", "Actinobacteria",
                   "Bacteroidetes", "Proteobacteria", "Fusobacteriia",
                   "Bacilli", "Clostridia", "Erysipelotrichia", "Negativicutes")
y4 <- subset(y4, is.na(Class)=="FALSE")
y4$Class <- factor(as.character(y4$Class), levels=desired_order)


myColors <- c("#999999", #Other
              "#CC0000", #Actinobacteria
              "#0000FF", # Bacteroidetes - Bacteroidia, Other
              "#99CC00", # Proteobacteria
              "#FF6699", #
              "#FFCC00","#CC9900","#FFCC66","#996600") # Firmicutes
names(myColors) <- levels(y4$Class)

Fig_2A <- ggplot(y4, aes(x=Location, y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  #facet_rep_grid(~Sample.Type, repeat.tick.labels='left') +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  scale_x_discrete(labels=function(Location) str_wrap(Location, width=3)) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(panel.spacing=unit(1, "lines"), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y=element_text(color="black", size = 8, vjust=1),
        axis.text.x=element_text(color="black", size = 8, vjust=0.5, hjust=0.5, angle=90),
        axis.title=element_text(size=9, color="black"),
        strip.background=element_blank(), 
        strip.text=element_text(size=9, color="black"),
        plot.tag=element_text(face="bold"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=8),
        legend.key.size=unit(0.5, "line"),
        legend.justification=0.5) +
  ylab("Mean relative abundance (%)") +
  xlab(NULL) +
  labs(tag="a")

Fig_2_c_legend <- cowplot::get_legend(ggplot(y4, aes(x=Location, y=Abundance, fill=Class)) +
                                        geom_bar(stat="identity") +
                                        #facet_rep_grid(~Sample.Type, repeat.tick.labels='left') +
                                        scale_fill_manual(name="Class", values=myColors) +
                                        scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
                                        scale_x_discrete(labels=function(Location) str_wrap(Location, width=3)) +
                                        guides(fill=guide_legend(ncol=2)) +
                                        theme_bw() +
                                        theme(panel.spacing=unit(1, "lines"), 
                                              panel.border = element_blank(),
                                              panel.grid.major = element_blank(), 
                                              panel.grid.minor = element_blank(),
                                              axis.line = element_line(colour = "black"),
                                              axis.text.y=element_text(color="black", size = 8, vjust=1),
                                              axis.text.x=element_text(color="black", size = 8, vjust=0.5, hjust=0.5, angle=90),
                                              axis.title=element_text(size=9, color="black"),
                                              strip.background=element_blank(), 
                                              strip.text=element_text(size=9, color="black"),
                                              plot.tag=element_text(face="bold"),
                                              legend.title=element_blank(),
                                              legend.position="bottom",
                                              legend.text=element_text(size=10),
                                              legend.key.size=unit(0.75, "line"),
                                              legend.justification=0.5) +
                                        ylab("Mean relative abundance (%)") +
                                        xlab(NULL) +
                                        labs(tag="a"))


grid.arrange(Fig.Adiv, Fig4B, ncol=2)

gs <- list(Fig.Adiv, Fig4B)

FIGURE_1 <- arrangeGrob(grobs=gs, ncol=2)

ggsave("Fig.Adiv.jpg", plot=Fig.Adiv, width=3.25, height=3.5, units="in", dpi=300)

ggsave("Fig1.pdf", plot = FIGURE_1, width=6.5, height=3.5, units="in", dpi=300)
ggsave("Fig1.jpg", plot = FIGURE_1, width=6.5, height=3.5, units="in", dpi=300)
















# (d) Differential abundance ####
diff.abund.plot.data <- diff.abund.location[[5]]
diff.abund.plot.data <- subset(diff.abund.plot.data, Genus %in% rownames(rgcca.scores))
diff.abund.plot.data[,c(1:14,19,21)] <- NULL
diff.abund.plot.data <- diff.abund.plot.data[order(diff.abund.plot.data$X7),]

diff.abund.plot.data$p.adj <- p.adjust(diff.abund.plot.data$model.Locationurban.Pr...t.., method="BH")

diff.abund.plot.data <- subset(diff.abund.plot.data, abs(X7) > 0.3)

diff.abund.plot.data$Genus <- factor(diff.abund.plot.data$Genus)
diff.abund.plot.data$Genus <- forcats::fct_inorder(diff.abund.plot.data$Genus)

abund.means <- as.data.frame(colMeans(feces.strict.sum.rel$Genus, na.rm=TRUE))
abund.means$Genus <- rownames(abund.means)
colnames(abund.means) <- c("Abundance", "Genus")
diff.abund.plot.data <- merge(diff.abund.plot.data, abund.means, by="Genus", all=FALSE)

diff.abund.plot.data$Significance <- cut(diff.abund.plot.data$p.adj,
                                         breaks=c(-Inf,0.05,Inf), labels=c("p<0.05","p>0.05"))

diff.abund.plot.data$Phylum <- as.character(diff.abund.plot.data$Phylum)
diff.abund.plot.data$Class <- as.character(diff.abund.plot.data$Class)

diff.abund.plot.data$Phylum[diff.abund.plot.data$Class=="Bacilli"] <- "Bacilli"
diff.abund.plot.data$Phylum[diff.abund.plot.data$Class=="Clostridia"] <- "Clostridia"
diff.abund.plot.data$Phylum[diff.abund.plot.data$Class=="Erysipelotrichia"] <- "Erysipelotrichia"
diff.abund.plot.data$Phylum[diff.abund.plot.data$Class=="Negativicutes"] <- "Negativicutes"

diff.abund.plot.data$Phylum <- as.factor(diff.abund.plot.data$Phylum)
diff.abund.plot.data$Phylum <- factor(diff.abund.plot.data$Phylum, levels=c("Actinobacteria", "Bacilli", "Bacteroidetes", "Clostridia",
                                                                            "Erysipelotrichia", "Fusobacteria", "Negativicutes", "Proteobacteria"))

myColors2 <- c("#CC0000", # Actinobacteria
               "#FFCC00", # Bacilli 
               "#0000FF", # Bacteroidetes - Bacteroidia, Other
               "#CC9900", # Clostridia
               "#FFCC66", # Erysipelotrichia
               "#FF6699", # Fusobacteria
               "#996600", # Negativicutes
               "#99CC00") # Proteobacteria
names(myColors2) <- levels(diff.abund.plot.data$Phylum)

levels(diff.abund.plot.data$Genus)[levels(diff.abund.plot.data$Genus)=="Clostridium_sensu_stricto"] <- "Clostridium"

plot.diff.abund.genus <- ggplot(diff.abund.plot.data, aes(x=X7, y=Genus)) + 
  geom_point(aes(size=Abundance, color=Phylum, shape=Significance)) +
  scale_color_manual(values=myColors2) +
  scale_shape_manual(values=c(19,1)) +
  geom_vline(xintercept=0, color="grey") +
  geom_vline(xintercept=0.2, color="grey", linetype=2) +
  #geom_vline(xintercept=0.5, color="grey", linetype=2) +
  geom_vline(xintercept=-0.2, color="grey", linetype=2) +
  #geom_vline(xintercept=-0.5, color="grey", linetype=2) +
  theme_bw() +
  labs(x="Hedge's g", y="Genus") +
  scale_y_discrete(position="right") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        legend.position="none")

diff.abund.legend <- cowplot::get_legend( ggplot(diff.abund.plot.data, aes(x=X7, y=Family)) + 
                                   geom_point(aes(size=Abundance, color=Phylum, shape=Significance)) +
                                   scale_color_manual(values=myColors2, guide=FALSE) +
                                   scale_shape_manual(values=c(19,1)) +
                                   scale_size_continuous(name="Abundance (%)") +
                                   geom_vline(xintercept=0, color="grey") +
                                   geom_vline(xintercept=0.2, color="grey", linetype=2) +
                                   #geom_vline(xintercept=0.5, color="grey", linetype=2) +
                                   geom_vline(xintercept=-0.2, color="grey", linetype=2) +
                                   #geom_vline(xintercept=-0.5, color="grey", linetype=2) +
                                     #guides(color=NULL) +
                                   theme_bw() +
                                   labs(x="Hedge's g", y="Genus") +
                                   scale_y_discrete(position="right") +
                                   theme(panel.grid=element_blank(),
                                         axis.text=element_text(color="black", size=8),
                                         axis.title=element_text(color="black", size=10),
                                         legend.title=element_text(color="black", size=10),
                                         legend.text=element_text(color="black", size=10),
                                         plot.tag=element_text(face="bold")))

# (e) Aitchison distance-based PCA, with vectors ####

ord.PCA.AIT <- rda(as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr"))))
ord.PCA.AIT.envfit <- envfit(ord.PCA.AIT, 
                             temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                                "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                                "Em.PCR.Overall", "index3.PC1","Location")],
                             permutations=1000, na.rm=TRUE)

scores <- scores(ord.PCA.AIT, choices=c(1:3), display="sites")
scores <- cbind(scores, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs <- as.data.frame(scores(ord.PCA.AIT.envfit, display="vectors"))
spp.scrs <- cbind(spp.scrs, Rsq = ord.PCA.AIT.envfit[["vectors"]]$r,
                  pvals = ord.PCA.AIT.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs))
spp.scrs <- subset(spp.scrs, pvals < 0.05)
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="Cementum.Age"] <- "Age"
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="Vol.Anthro"] <- "Vol. anthro"
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="index3.PC1"] <- "Health"
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"

# Confidence ellipse.
plot.new()
ellipse <- ordiellipse(ord.PCA.AIT, 
                       cy.metadata.pseq$Location, 
                       display="sites", 
                       kind="sd", 
                       conf=0.95, 
                       label=T)
ellipse.data <- data.frame()
for(g in levels(scores$Location)){
  ellipse.data <- rbind(ellipse.data, cbind(as.data.frame(with(scores[scores$Location==g,],
                                                               veganCovEllipse(ellipse[[g]]$cov,
                                                                               ellipse[[g]]$center,
                                                                               ellipse[[g]]$scale)))
                                            ,Type=g))
}
rm(ellipse)

ord.plot <- ggplot(scores) + 
  geom_point(mapping = aes(x=PC1, y=PC2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs,
               aes(x=0, y=0, xend=11*PC1, yend=11*PC2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ellipse.data, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  #geom_text(data=subset(spp.scrs, Variable=="Health"), aes(x=11*PC1-0.2, y=11*PC2, label=Variable), size=3.175, hjust=1, vjust=0.5) +
  #geom_text(data=subset(spp.scrs, Variable=="Vol. anthro"), aes(x=11*PC1-0.1, y=11*PC2, label=Variable), size=3.175, hjust=1, vjust=0.5) +
  #geom_text(data=subset(spp.scrs, Variable=="Age"), aes(x=11*PC1+0.1, y=11*PC2, label=Variable), size=3.175, hjust=0, vjust=0.5) +
  #geom_text(data=subset(spp.scrs, Variable=="d13C"), aes(x=11*PC1+0.2, y=11*PC2, label=Variable), size=3.175, hjust=0, vjust=0.5) +
  #geom_text(data=subset(spp.scrs, Variable=="Spleen mass"), aes(x=11*PC1+0.2, y=11*PC2, label=Variable), size=3.175, hjust=0, vjust=0.5) +
  theme_bw() +
  labs(x="PC1 (11.0%)", y="PC2 (10.0%)", tag="b") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

# - Figure design ####
grid.arrange(adiv.urb.rur.plot, 
             Fig.Adiv+theme(legend.position="none")+labs(tag="b"),
             Fig_2A+theme(legend.position="none")+labs(tag="c"), 
             plot.diff.abund+theme(legend.position="none")+labs(tag="d"),
             ord.plot+theme(legend.position="none")+labs(tag="e"),
             layout_matrix=rbind(c(1,3,4),c(2,3,5)), widths=c(0.4,0.2,0.4))



FIGURE_2_REAL <- arrangeGrob(grobs=list(adiv.urb.rur.plot, 
                             Fig.Adiv+theme(legend.position="none")+labs(tag="b"),
                             Fig_2A+theme(legend.position="none")+labs(tag="c"), 
                             plot.diff.abund.genus+theme(legend.position="none")+labs(tag="d"),
                             ord.plot+theme(legend.position="none")+labs(tag="e")),
                             layout_matrix=rbind(c(1,3,4),c(2,3,5)), widths=c(0.4,0.2,0.4))

grid.arrange(plot.legend, Fig2_b_legend, Fig_2_c_legend, diff.abund.legend)
FIGURE_2_LEGENDS <- arrangeGrob(grobs=list(plot.legend, Fig2_b_legend, Fig_2_c_legend, diff.abund.legend))

ggsave("FIGURE_2_REAL.jpg", plot=FIGURE_2_REAL, width=6.5, height=4.5, units="in", dpi=300)
ggsave("FIGURE_2_LEGENDS.jpg", plot=FIGURE_2_LEGENDS, dpi=300)



#### FIGURE 3: TAXON SWOOPS WITH COEFFICIENTS ####
# (a) Spearman correlation heat map ####
temp <- comp.strict.clr[["Genus"]]

# Separate into taxon and diet/health categories.
temp <- temp[,sapply(temp, is.numeric)]
temp <- as.data.frame(t(temp))
temp1 <- subset(temp, rownames(temp) %in% colnames(feces.strict.sum.clr))
temp2 <- subset(temp, rownames(temp) %in% c("Mass","Snout.to.Base","Girth",
                                            "Cementum.Age","Spleen.to.Body.Ratio","KFI",
                                            "d13C", "d15N", "VOL.PREY", "Vol.Anthro", "VOL.VEG"))

temp <- rbind(temp1, temp2)
rm(temp1, temp2, prev)
temp <- as.data.frame(t(temp))

# Calculate correlations between each taxa and metadata variable.
temp.spearman <- spearman.correlations(temp)
temp.spearman.R <- temp.spearman[["r"]]

# Melt correlation coefficients and clean data frame for plotting.
temp.spearman.melt <- reshape2::melt(temp.spearman.R)
temp.spearman.melt <- subset(temp.spearman.melt, Var1 %in% c("Mass",
                                                             "Snout.to.Base",
                                                             "Girth",
                                                             "Cementum.Age",
                                                             "Spleen.to.Body.Ratio",
                                                             "KFI",
                                                             "d13C",
                                                             "d15N",
                                                             "VOL.PREY",
                                                             "Vol.Anthro",
                                                             "VOL.VEG"))
temp.spearman.melt <- subset(temp.spearman.melt, !(Var2 %in% c("Mass",
                                                               "Snout.to.Base",
                                                               "Girth",
                                                               "Cementum.Age",
                                                               "Spleen.to.Body.Ratio",
                                                               "KFI",
                                                               "d13C",
                                                               "d15N",
                                                               "VOL.PREY",
                                                               "Vol.Anthro",
                                                               "VOL.VEG")))
temp.spearman.melt$Var1 <- factor(temp.spearman.melt$Var1)
temp.spearman.melt$Var2 <- factor(temp.spearman.melt$Var2)

# Arrange families alphabetically by phylum.
tax.table <- as.data.frame(cbind(tax_table(temp.physeq)))
tax.table <- subset(tax.table, Genus %in% levels(temp.spearman.melt$Var2))
tax.table <- dplyr::arrange(tax.table, Phylum, Class, Family, Genus)
tax.table <- dplyr::mutate(tax.table, Genus = factor(Genus, Genus))
temp.spearman.melt$Var2 <- factor(as.character(temp.spearman.melt$Var2), levels=levels(tax.table$Genus))

# Assign colornames based on phylum.
colornames <- c("#CC0000", "#CC0000", # Actinobacteria x2
                "magenta3", "magenta3", # Bacteroidetes x2
                "#996600", "#996600", "#996600", "#996600", "#996600", "#996600", 
                "#996600", "#996600", "#996600", "#996600", "#996600", "#996600",
                "#996600", "#996600", "#996600", "#996600", # Firmicutes x13
                "darkgreen", # Fusobacteria x1
                "#0000FF", "#0000FF", "#0000FF", "#0000FF") # Proteobacteria x4
names(colornames) <- levels(temp.spearman.melt$Var2)

# Rename variables for plotting.
levels(temp.spearman.melt$Var1)[levels(temp.spearman.melt$Var1)=="Snout.to.Base"] <- "Body size"
levels(temp.spearman.melt$Var1)[levels(temp.spearman.melt$Var1)=="Cementum.Age"] <- "Age"
levels(temp.spearman.melt$Var1)[levels(temp.spearman.melt$Var1)=="Spleen.to.Body.Ratio"] <- "Spleen mass"
levels(temp.spearman.melt$Var1)[levels(temp.spearman.melt$Var1)=="VOL.PREY"] <- "Vol prey"
levels(temp.spearman.melt$Var1)[levels(temp.spearman.melt$Var1)=="Vol.Anthro"] <- "Vol anthro"
levels(temp.spearman.melt$Var1)[levels(temp.spearman.melt$Var1)=="VOL.VEG"] <- "Vol veg"

# Add an empty factor level to create space between diet and condition variables, then factor in order for plotting.
temp.spearman.melt$Var1 <- factor(temp.spearman.melt$Var1, levels = c(levels(temp.spearman.melt$Var1), " "))
desired_order <- c("KFI", "Spleen mass", "Age", "Girth", "Body size", "Mass", " ", "d15N", "d13C", "Vol veg", "Vol anthro", "Vol prey")
temp.spearman.melt$Var1 <- factor(as.character(temp.spearman.melt$Var1), levels=desired_order)

# Clean taxa names
levels(temp.spearman.melt$Var2)[levels(temp.spearman.melt$Var2)=="Clostridium_sensu_stricto"] <- "Clostridium"
levels(temp.spearman.melt$Var2)[levels(temp.spearman.melt$Var2)=="Ruminococcus2"] <- "Ruminococcus"

# Clean workspace.
rm(temp, tax.table, temp.spearman, temp.spearman.R, temp.physeq, desired_order)

# Generate plot
genus.heatmap <- ggplot(temp.spearman.melt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue",
                       high = "red",
                       mid = "white", 
                       midpoint = 0,
                       limit = c(-0.5,0.5),
                       space = "Lab", 
                       name="Spearman's R") +
  theme_bw() + 
  scale_y_discrete(limits = levels(temp.spearman.melt$Var1), 
                   labels=c('KFI', 'Spleen mass', 'Age', 'Girth', 'Body size',
                            'Mass', ' ', expression(paste(delta^{15}, "N")), expression(paste(delta^{13}, "C")), 
                            'Vol veg', 'Vol anthro', 'Vol prey', 'Diet diversity')) +
  scale_x_discrete(limits = levels(temp.spearman.melt$Var2), position="top") +
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 8, hjust = 0, color=colornames),
        axis.text.y = element_text(color="black", size=8),
        axis.title = element_text(color="black", size=9),
        axis.line = element_blank(),
        panel.spacing=unit(2, "lines"),
        panel.border = element_blank(),
        strip.text = element_text(color="black", size=9),
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(color="black", size=9),
        legend.position="right",
        legend.text=element_text(color="black", size=8),
        legend.key.size = unit(1,"lines"),
        plot.tag=element_text(color="black", face="bold")) +
  xlab("Genus") +
  ylab("Measure") + labs(tag="a") +
  coord_fixed()

# (b) Swoop with health score ####
swoop2.data <- coef.health.regression.ranks[["Genus"]]

swoop2.data <- swoop2.data[order(swoop2.data$Coefficient),]
swoop2.data$Rank <- rep(1:nrow(swoop2.data))
swoop2.data$Taxon <- forcats::fct_inorder(swoop2.data$Taxon)

levels(swoop2.data$Taxon)[levels(swoop2.data$Taxon)=="Clostridium_sensu_stricto"] <- "Clostridium"
levels(swoop2.data$Taxon)[levels(swoop2.data$Taxon)=="Ruminococcus2"] <- "Ruminococcus"

# Expand the axis to make the plot look nicer.
swoop2.data$Coefficient <- 10*swoop2.data$Coefficient
swoop2.data$Lower <- 10*swoop2.data$Lower
swoop2.data$Upper <- 10*swoop2.data$Upper

swoop2 <- ggplot(swoop2.data, aes(x=Coefficient, y=Taxon)) +
  theme_bw() +
  geom_errorbarh(aes(xmin=swoop2.data$Lower, xmax=swoop2.data$Upper, y=Taxon), position=position_dodge(0.9), height=0, color="darkgrey") +
  geom_point(aes(size=Abundance), color="black") +
  geom_vline(xintercept=0) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=9),
        legend.position="none",
        plot.tag=element_text(face="bold")) +
  labs(x="Coef: Health", y="Genus", tag="b")

# (c) Swoop with spleen size ####
swoop3.data <- coef.spleen.regression.ranks[["Genus"]]
swoop3.data <- swoop3.data[order(swoop3.data$Coefficient),]
swoop3.data$Rank <- rep(1:nrow(swoop2.data))
swoop3.data$Taxon <- forcats::fct_inorder(swoop3.data$Taxon)

levels(swoop3.data$Taxon)[levels(swoop3.data$Taxon)=="Clostridium_sensu_stricto"] <- "Clostridium"
levels(swoop3.data$Taxon)[levels(swoop3.data$Taxon)=="Ruminococcus2"] <- "Ruminococcus"

swoop3 <- ggplot(swoop3.data, aes(x=Coefficient, y=Taxon)) +
  theme_bw() +
  geom_errorbarh(aes(xmin=swoop3.data$Lower, xmax=swoop3.data$Upper, y=Taxon), position=position_dodge(0.9), height=0, color="darkgrey") +
  geom_point(aes(size=Abundance), color="black") +
  geom_vline(xintercept=0) +
  scale_y_discrete(position="right") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=9),
        legend.position="none",
        plot.tag=element_text(face="bold")) +
  labs(x="Coef: Spleen mass", y="Genus", tag="c")

# - Figure design ####
grid.arrange(grobs=list(genus.heatmap, swoop2, swoop3), layout_matrix=rbind(c(1,1),c(2,3)), heights=c(0.5,0.5))

Figure_3 <- arrangeGrob(grobs=list(genus.heatmap, swoop2, swoop3), layout_matrix=rbind(c(1,1),c(2,3)), heights=c(0.5,0.5))
ggsave("~/health/Tables/Figure_3_taxon_correlations.jpg", Figure_3, dpi=300, width=6.5, height=8, units="in")


#### FIGURE 4: RCCA CLUSTERED IMAGE MAP ####
# Extract dendrograms.
rgcca.taxa.dendro <- dendro_data(rgcca.cim$ddr, type="rectangle")
rgcca.data.dendro <- dendro_data(rgcca.cim$ddc, type="rectangle")

# Extract heatmap data.
rgcca.heatmap <- as.data.frame(rgcca.cim$mat)
rgcca.heatmap$Taxon <- rownames(rgcca.heatmap)
levels(rgcca.heatmap$Taxon) <- rgcca.cim$row.names
rgcca.heatmap$Taxon <- forcats::fct_inorder(rgcca.heatmap$Taxon)
rgcca.heatmap <- reshape2::melt(rgcca.heatmap)

levels(rgcca.heatmap$Taxon)[levels(rgcca.heatmap$Taxon)=="Clostridium_sensu_stricto"] <- "Clostridium"
levels(rgcca.heatmap$Taxon)[levels(rgcca.heatmap$Taxon)=="Ruminococcus2"] <- "Ruminococcus"

# Plot the heat map with dendrograms.
A <- ggplot(rgcca.heatmap, aes(x=variable, y=Taxon)) + geom_tile(aes(fill=value)) +
  scale_y_discrete(position="right") +
  scale_x_discrete(labels=c("Age", "Health", "Vol. prey", expression(paste(delta^{15}, "N")),
                            "Spleen mass", expression(paste(delta^{13}, "C")), "Vol. anthro")) +
  geom_segment(data=segment(rgcca.taxa.dendro), aes(x=-2*y+0.5, y=x, xend=-2*yend+0.5, yend=xend)) +
  geom_segment(data=segment(rgcca.data.dendro), aes(x=x, y=2*y+0.5+ncol(temp.otu.newgen), xend=xend, yend=2*yend+0.5+ncol(temp.otu.newgen))) +
  theme_bw() +
  scale_fill_gradientn(colors=color.spectral(25), limits=c(-0.5,0.52), name="Correlation\nScore") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle=90, color="black", size=9, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="black", size=9),
        axis.title = element_text(color="black", size=10),
        legend.title= element_text(color="black", size=10, hjust=0.5),
        legend.position="none",
        plot.tag=element_text(face="bold")) +
  labs(x="Measure", y="Genus")

# Create a bar graph showing the sum of absolute values of correlation scores, as an indicator of taxon importance.
rgcca.scores$Taxon <- as.factor(rownames(rgcca.scores))
levels(rgcca.scores$Taxon)[levels(rgcca.scores$Taxon)=="Clostridium_sensu_stricto"] <- "Clostridium"
levels(rgcca.scores$Taxon)[levels(rgcca.scores$Taxon)=="Ruminococcus2"] <- "Ruminococcus"
rgcca.scores$Taxon <- factor(rgcca.scores$Taxon, levels=levels(rgcca.heatmap$Taxon))

B <- ggplot(rgcca.scores, aes(x=Taxon, y=Important)) + geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  labs(y=NULL, x=NULL) +
  theme(axis.text=element_blank(),
        panel.grid=element_blank(),
        plot.margin=unit(c(50,5.5,73,5.5), unit="pt"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank())

grid.arrange(A, B, ncol=2, widths=c(0.8,0.2))

# Extract the legend for the heat map.
C <- cowplot::get_legend(ggplot(rgcca.heatmap, aes(x=variable, y=Taxon)) + geom_tile(aes(fill=value)) +
                           scale_y_discrete(position="right") +
                           scale_x_discrete(labels=c(expression(paste(delta^{15}, "N")), "Vol. prey", "Age", "Health", "Vol. anthro",
                                                     expression(paste(delta^{13}, "C")), "Spleen mass")) +
                           geom_segment(data=segment(rgcca.taxa.dendro), aes(x=-2*y+0.5, y=x, xend=-2*yend+0.5, yend=xend)) +
                           geom_segment(data=segment(rgcca.data.dendro), aes(x=x, y=2*y+0.5+ncol(temp.otu), xend=xend, yend=2*yend+0.5+ncol(temp.otu))) +
                           theme_bw() +
                           scale_fill_gradientn(colors=color.spectral(25), limits=c(-0.5,0.52), name="Correlation\nScore") +
                           theme(panel.grid = element_blank(),
                                 panel.border = element_blank(),
                                 axis.text.x = element_text(angle=90, color="black", size=9, hjust=1, vjust=0.5),
                                 axis.text.y = element_text(color="black", size=9),
                                 axis.title = element_text(color="black", size=10),
                                 legend.title= element_text(color="black", size=10, hjust=0.5),
                                 legend.position="right",
                                 plot.tag=element_text(face="bold")) +
                           labs(x="Measure", y="Family"))

# Figure assembly. #
FIGURE_4 <- arrangeGrob(grobs=list(A, B), ncol=2, widths=c(0.8, 0.2))
ggsave("~/health/Tables/FIGURE_4.jpg", FIGURE_4, width=6.5, height=6, units="in")
ggsave("~/health/Tables/Figure_4_legend.jpg", C, width=2, height=2, units="in")

#### FIGURE 5: STRUCTURAL EQUATION MODELS ####
# These models were exported to CytoScape for figure design.
#### FIG S1: EMPTY STOMACHS AND SPLEEN MASS, AND EM INFECTION VS. AGE ####
supp_plot_data <- cy.metadata

levels(supp_plot_data$PA.Empty)[levels(supp_plot_data$PA.Empty)==1] <- "Y"
levels(supp_plot_data$PA.Empty)[levels(supp_plot_data$PA.Empty)==0] <- "N"
levels(supp_plot_data$Location)[levels(supp_plot_data$Location)=="peri-urban"] <- "rural"

Fig_S2a <- ggplot(supp_plot_data, aes(x=PA.Empty, y=Spleen.to.Body.Ratio)) +
  geom_boxplot(aes(color=Location, fill=PA.Empty)) +
  facet_grid(~Location) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_fill_manual(values=c("antiquewhite", "darkgrey")) +
  guides(color=guide_legend(title.position="top"), fill=guide_legend(title.position="top", title="Empty stomach")) +
  labs(x="Empty stomach (Y/N)", y="Spleen mass (g/kg)", tag="a") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        strip.background=element_blank(),
        panel.spacing=unit(0.5, "lines"),
        strip.text=element_text(color="black", size=10),
        legend.position="bottom",
        legend.title=element_text(color="black", size=10, hjust=0.5),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

Fig_S2b <- ggplot(subset(supp_plot_data, is.na(Em.PCR.Overall)=="FALSE"), aes(x=Em.PCR.Overall, y=Cementum.Age)) +
  geom_boxplot(aes(color=Location, fill=Em.PCR.Overall)) +
  facet_grid(~Location) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_fill_manual(values=c("antiquewhite", "darkgrey")) +
  guides(color=guide_legend(title.position="top"), fill=guide_legend(title.position="top", title="E. multilocularis")) +
  labs(x="E. multilocularis infection (Y/N)", y="Age (yr)", tag="b") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        strip.background=element_blank(),
        panel.spacing=unit(0.5, "lines"),
        strip.text=element_text(color="black", size=10),
        legend.position="bottom",
        legend.title=element_text(color="black", size=10, hjust=0.5),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

grid.arrange(Fig_S2a, Fig_S2b, nrow=1)
Fig_S2 <- arrangeGrob(grobs=list(Fig_S2a, Fig_S2b), nrow=1)
ggsave("~/health/Tables/Figure_S2_spleen_stomach_and_age_Em.jpg", Fig_S2, dpi=300, width=6, height=3.5, units="in")

#### FIG S2: ALPHA DIVERSITY MEASURES VS. EM, CONTROLLING FOR LOCATION ####
supp_plot_data <- reshape2::melt(supp_plot_data)
supp_plot_data <- subset(supp_plot_data, variable %in% c("Observed", "Shannon", "PD", "NTI"))
supp_plot_data$variable <- factor(supp_plot_data$variable, levels=c("Observed", "Shannon", "PD", "NTI"))

levels(supp_plot_data$variable)[levels(supp_plot_data$variable)=="Observed"] <- "ASV Richness"
colnames(supp_plot_data)[colnames(supp_plot_data)=="Em.PCR.Overall"] <- "E. multi"

Fig_S3 <- ggplot(supp_plot_data, aes(x=Location, y=value, color=Location, fill=`E. multi`)) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free", nrow=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  #guides(color=FALSE) +
  scale_fill_manual(values=c("antiquewhite", "darkgrey")) +
  labs(x="\nLocation / E. multilocularis infection", y="Value") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        strip.background=element_blank(),
        panel.spacing=unit(0.5, "lines"),
        strip.text=element_text(color="black", size=10),
        legend.position="right",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

ggsave("~/health/Tables/Figure_S3_Em_alpha_diversity.jpg", Fig_S3, dpi=300, width=6, height=3, units="in")


#### FIG S3: ALPHA DIVERSITY DREDGED MODEL PLOTS ####
adiv.model.plots <- list()
adiv.model.plots[["Observed"]] <- ggplot(adiv.model.results[["Observed"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Location", "E. multi", "Age", "Health", "Sex",
                            expression(paste(delta^{15}, "N")), "Vol. anthro", expression(paste(delta^{13}, "C")),
                            "Spleen mass", "Vol. prey")) +
  coord_flip() +
  theme_bw() +
  labs(x="Coefficient", y="Predictor", tag="a", title="ASV richness") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="none",
        axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

adiv.model.plots[["Shannon"]] <- ggplot(adiv.model.results[["Shannon"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Location", "E. multi", expression(paste(delta^{15}, "N")), "Sex",
                            expression(paste(delta^{13}, "C")), "Spleen mass", "Health", "Age", "Vol. prey", "Vol. anthro")) +
  coord_flip() +
  theme_bw() +
  labs(x="Coefficient", y="Predictor", tag="b", title="Shannon diversity") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="none",
        axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

adiv.model.plots[["PD"]] <- ggplot(adiv.model.results[["PD"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Location", "Health", "E. multi", expression(paste(delta^{15}, "N")), "Sex",
                            "Age", expression(paste(delta^{13}, "C")), "Vol. anthro", "Vol. prey", "Spleen mass")) +
  coord_flip() +
  theme_bw() +
  labs(x="Coefficient", y="Predictor", tag="c", title="Faith's PD") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="none",
        axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

adiv.model.plots[["NTI"]] <- ggplot(adiv.model.results[["NTI"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Location", "Health", "Vol. prey", "Spleen mass", "E. multi", "Vol. anthro", "Age",
                            expression(paste(delta^{15}, "N")), expression(paste(delta^{13}, "C")), "Sex")) +
  coord_flip() +
  theme_bw() +
  labs(x="Coefficient", y="Predictor", tag="d", title="NTI") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="none",
        axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

Fig_S4 <- arrangeGrob(grobs=list(adiv.model.plots[["Observed"]],
                                 adiv.model.plots[["Shannon"]],
                                 adiv.model.plots[["PD"]],
                                 adiv.model.plots[["NTI"]]), nrow=2, ncol=2)

ggsave("~/health/Tables/Figure_S4_adiv_dredged_models.jpg", Fig_S4, dpi=300, width=6.5, height=6.5, units="in")

#### FIG S4: ALPHA DIVERSITY AND EMPTY STOMACHS ####
Fig_S5 <- ggplot(supp_plot_data, aes(x=Location, y=value, color=Location, fill=PA.Empty)) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free", nrow=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  #guides(color=FALSE) +
  scale_fill_manual(values=c("antiquewhite", "darkgrey"), name="Empty\nstomach?") +
  labs(x="\nLocation / Empty stomach", y="Value") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        strip.background=element_blank(),
        panel.spacing=unit(0.5, "lines"),
        strip.text=element_text(color="black", size=10),
        legend.position="right",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

ggsave("~/health/Tables/Figure_S5_empty_stomach_alpha_diversity.jpg", Fig_S5, dpi=300, width=6, height=3, units="in")

#### FIG S5: DIFFERENTIAL ABUNDANCE BASED ON AGE, SEX, ETC. ####
# (a) Sex only (controlled for age) ####
diff.abund.supp.data <- diff.abund.just.sex.and.age[[5]]
diff.abund.supp.data <- subset(diff.abund.supp.data, Genus %in% colnames(temp.otu.newgen))
diff.abund.supp.data <- diff.abund.supp.data[order(diff.abund.supp.data$X7),]

diff.abund.supp.data <- subset(diff.abund.supp.data, abs(X7) > 0.25)

diff.abund.supp.data$Genus <- factor(diff.abund.supp.data$Genus)
diff.abund.supp.data$Genus <- forcats::fct_inorder(diff.abund.supp.data$Genus)

abund.means <- as.data.frame(colMeans(feces.strict.sum.rel$Genus, na.rm=TRUE))
abund.means$Genus <- rownames(abund.means)
colnames(abund.means) <- c("Abundance", "Genus")
diff.abund.supp.data <- merge(diff.abund.supp.data, abund.means, by="Genus", all=FALSE)

diff.abund.supp.data$p.adj <- p.adjust(diff.abund.supp.data$model.SexM.Pr...t.., method="BH")
diff.abund.supp.data$shape <- cut(diff.abund.supp.data$p.adj,
                                  breaks=c(-Inf,0.1,Inf), labels=c("significant","not"))

diff.abund.supp.data$Phylum <- as.character(diff.abund.supp.data$Phylum)
diff.abund.supp.data$Class <- as.character(diff.abund.supp.data$Class)

diff.abund.supp.data$Phylum <- as.factor(diff.abund.supp.data$Phylum)

myColors2 <- c("#CC0000", # Actinobacteria
               "#0000FF", # Bacteroidetes - Bacteroidia, Other
               "#996600", # Firmicutes
               "#FF6699", # Fusobacteria
               "#99CC00") # Proteobacteria
names(myColors2) <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria")

plot.diff.abund.sex <- ggplot(diff.abund.supp.data, aes(x=X7, y=Genus)) + 
  geom_point(aes(size=Abundance, color=Phylum, shape=shape)) +
  scale_color_manual(values=myColors2) +
  scale_shape_manual(values=c(19)) +
  scale_size_continuous(range=c(2,6)) +
  geom_vline(xintercept=0, color="grey") +
  geom_vline(xintercept=0.2, color="grey", linetype=2) +
  geom_vline(xintercept=-0.2, color="grey", linetype=2) +
  theme_bw() +
  labs(x="Hedge's g", y="Genus", tag="a", title="Sex") +
  scale_y_discrete(position="left") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=10),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=10, color="black"),
        legend.position="none")

# (b) Em only (controlled for location) ####
diff.abund.supp.data.em <- diff.abund.loc.Em[[5]]
diff.abund.supp.data.em <- subset(diff.abund.supp.data.em, Genus %in% colnames(temp.otu.newgen))
diff.abund.supp.data.em <- diff.abund.supp.data.em[order(diff.abund.supp.data.em$X7),]

diff.abund.supp.data.em <- subset(diff.abund.supp.data.em, abs(X7) > 0.25)

diff.abund.supp.data.em$Genus <- factor(diff.abund.supp.data.em$Genus)
diff.abund.supp.data.em$Genus <- forcats::fct_inorder(diff.abund.supp.data.em$Genus)

abund.means <- as.data.frame(colMeans(feces.strict.sum.rel$Genus, na.rm=TRUE))
abund.means$Genus <- rownames(abund.means)
colnames(abund.means) <- c("Abundance", "Genus")
diff.abund.supp.data.em <- merge(diff.abund.supp.data.em, abund.means, by="Genus", all=FALSE)

diff.abund.supp.data.em$p.adj <- p.adjust(diff.abund.supp.data.em$model.Em.PCR.OverallY.Pr...t.., method="BH")
diff.abund.supp.data.em$shape <- cut(diff.abund.supp.data.em$p.adj,
                                     breaks=c(-Inf,0.1,Inf), labels=c("significant","not"))

diff.abund.supp.data.em$Phylum <- as.character(diff.abund.supp.data.em$Phylum)
diff.abund.supp.data.em$Class <- as.character(diff.abund.supp.data.em$Class)

diff.abund.supp.data.em$Phylum <- as.factor(diff.abund.supp.data.em$Phylum)

levels(diff.abund.supp.data.em$Family)[levels(diff.abund.supp.data.em$Family)=="Clostridiaceae_1"] <- "Clostridiaceae"
levels(diff.abund.supp.data.em$Family)[levels(diff.abund.supp.data.em$Family)=="Peptococcaceae_1"] <- "Peptococcaceae"

plot.diff.abund.em <- ggplot(diff.abund.supp.data.em, aes(x=X7, y=Genus)) + 
  geom_point(aes(size=Abundance, color=Phylum, shape=shape)) +
  scale_color_manual(values=myColors2) +
  scale_shape_manual(values=c(19)) +
  scale_size_continuous(range=c(2,6)) +
  geom_vline(xintercept=0, color="grey") +
  geom_vline(xintercept=0.2, color="grey", linetype=2) +
  geom_vline(xintercept=-0.2, color="grey", linetype=2) +
  theme_bw() +
  labs(x="Hedge's g", y="Genus", tag="b", title="E. multilocularis") +
  scale_y_discrete(position="left") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=10),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=10, color="black"),
        legend.position="none")

# (c) Age (controlled for sex) ####
diff.abund.supp.data.age <- diff.abund.just.sex.and.age[[5]]
diff.abund.supp.data.age <- subset(diff.abund.supp.data.age, Genus %in% colnames(temp.otu.newgen))
abund.means <- as.data.frame(colMeans(feces.strict.sum.rel$Genus, na.rm=TRUE))
abund.means$Genus <- rownames(abund.means)
colnames(abund.means) <- c("Abundance", "Genus")
diff.abund.supp.data.age <- merge(diff.abund.supp.data.age, abund.means, by="Genus", all=FALSE)

temp.corr <- spearman.correlations(temp.otu.newgen.comp)
temp.corr <- as.data.frame(temp.corr[["r"]])
temp.corr <- subset(temp.corr, rownames(temp.corr) %in% diff.abund.supp.data.age$Genus)
temp.corr <- data.frame(Genus = rownames(temp.corr),
                        Age.R = temp.corr[,"Cementum.Age"])

diff.abund.supp.data.age <- merge(diff.abund.supp.data.age, temp.corr, by="Genus", all=FALSE)
diff.abund.supp.data.age <- diff.abund.supp.data.age[order(diff.abund.supp.data.age$model.Age.Estimate),]

diff.abund.supp.data.age <- subset(diff.abund.supp.data.age, abs(Age.R) > 0.2)

diff.abund.supp.data.age$Genus <- factor(diff.abund.supp.data.age$Genus)
diff.abund.supp.data.age$Genus <- forcats::fct_inorder(diff.abund.supp.data.age$Genus)

diff.abund.supp.data.age$p.adj <- p.adjust(diff.abund.supp.data.age$model.Age.Pr...t.., method="BH")
diff.abund.supp.data.age$shape <- cut(diff.abund.supp.data.age$p.adj,
                                      breaks=c(-Inf,0.1,Inf), labels=c("significant","not"))

diff.abund.supp.data.age$Phylum <- as.character(diff.abund.supp.data.age$Phylum)
diff.abund.supp.data.age$Class <- as.character(diff.abund.supp.data.age$Class)

diff.abund.supp.data.age$Phylum <- as.factor(diff.abund.supp.data.age$Phylum)

plot.diff.abund.age <- ggplot(diff.abund.supp.data.age, aes(x=model.Age.Estimate, y=Genus)) + 
  geom_point(aes(size=Abundance, color=Phylum, shape=shape)) +
  scale_color_manual(values=myColors2) +
  scale_shape_manual(values=c(19, 19)) +
  scale_size_continuous(range=c(2,6)) +
  geom_vline(xintercept=0, color="grey") +
  geom_vline(xintercept=0.2, color="grey", linetype=2) +
  geom_vline(xintercept=-0.2, color="grey", linetype=2) +
  theme_bw() +
  labs(x="Coefficient", y="Genus", tag="c", title="Age") +
  scale_y_discrete(position="left") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=10),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=10),
        legend.position="none")

# - Figure design ####
diff.abund.legend.data <- subset(diff.abund.location[[5]], Phylum %in% c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria"))
diff.abund.legend.data$Phylum <- factor(diff.abund.legend.data$Phylum)
abund.means <- as.data.frame(colMeans(feces.strict.sum.rel$Genus, na.rm=TRUE))
abund.means$Genus <- rownames(abund.means)
colnames(abund.means) <- c("Abundance", "Genus")
diff.abund.legend.data <- merge(diff.abund.legend.data, abund.means, by="Genus", all=FALSE)
rm(abund.means)

Fig_S6_legend <- cowplot::get_legend(ggplot(diff.abund.legend.data, aes(x=X7, y=Genus)) + 
                                       geom_point(aes(size=Abundance, color=Phylum)) +
                                       scale_color_manual(values=myColors2) +
                                       scale_shape_manual(values=c(19)) +
                                       scale_size_continuous(range=c(2,6)) +
                                       guides(shape=FALSE) +
                                       geom_vline(xintercept=0, color="grey") +
                                       geom_vline(xintercept=0.2, color="grey", linetype=2) +
                                       geom_vline(xintercept=-0.2, color="grey", linetype=2) +
                                       theme_bw() +
                                       labs(x="Coefficient", y="Family", tag="c", title="Age", size="Relative\nAbundance (%)", color="Phylum") +
                                       scale_y_discrete(position="left") +
                                       theme(panel.grid=element_blank(),
                                             axis.text=element_text(color="black", size=8),
                                             axis.title=element_text(color="black", size=10),
                                             plot.tag=element_text(face="bold"),
                                             plot.title=element_text(size=10),
                                             legend.title=element_text(size=10, hjust=0.5),
                                             legend.position="right"))


grid.arrange(plot.diff.abund.sex, plot.diff.abund.em, plot.diff.abund.age, layout_matrix=rbind(c(1,1,2,2),c(NA,3,3,NA)))

Fig_S6 <- arrangeGrob(grobs=list(plot.diff.abund.sex, plot.diff.abund.em, plot.diff.abund.age, Fig_S6_legend),
                      layout_matrix=rbind(c(1,1,2,2),c(3,3,4,4)))

ggsave("~/health/Tables/Figure_S6_differential_abundance_tests.jpg", Fig_S6, dpi=300, width=6.5, height=6, units="in")

rm(diff.abund.supp.data, diff.abund.supp.data.age, diff.abund.supp.data.em, diff.abund.supp.data.health, temp.corr)

#### FIG S6: RANDOM FOREST RESULTS ####
# Sort graph results by decreasing Gini scores.
rf.plot.data <- rf.importance[order(-rf.importance$MeanDecreaseGini),]

# Insert ASVID as a row name.
rf.plot.data$ASVID <- rownames(rf.plot.data)

# Order plots by ASVID.
rf.plot.data$ASVID <- forcats::fct_inorder(rf.plot.data$ASVID)

# Determine which grouping (segment: small/large, individual: identity, site: intestinal site) each taxon is most abundant in.
rf.plot.data$Max <- as.factor(colnames(rf.plot.data[,8:9])[max.col(rf.plot.data[,8:9], ties.method="first")])
rf.plot.data$Max <- factor(rf.plot.data$Max, levels=c("peri-urban","urban"))
levels(rf.plot.data$Max)[levels(rf.plot.data$Max)=="peri-urban"] <- "rural"

rf.plot.data <- rf.plot.data[c(1:20),]
rf.plot.data$ASVID <- factor(rf.plot.data$ASVID)

levels(rf.plot.data$Genus)[levels(rf.plot.data$Genus)=="Clostridium_sensu_stricto"] <- "Clostridium"

Fig_S7 <- ggplot(rf.plot.data, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Overall)) + 
  scale_y_discrete(limits=rev(levels(rf.plot.data$ASVID)), labels=rev(as.character(rf.plot.data$Genus))) +
  scale_color_manual(name="Location", values=c("forestgreen", "purple")) +
  scale_size_continuous(name="Abundance (%)", range=c(3,6)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n85.2%", size=3.75) +
  theme_bw() +
  theme(legend.position=c("bottom"), 
        legend.justification=c("center"),
        legend.title=element_text(size=10, color="black"),
        legend.text=element_text(size=9, color="black"),
        plot.tag=element_text(face="bold"),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(size=9, color="black"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=10, color="black")) +
  guides(colour=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=1)) +
  guides(size=guide_legend(ncol=1, title.position="top", title.hjust=0.5, order=0)) +
  labs(y=NULL, x="Mean decrease Gini coefficient")


ggsave("~/health/Tables/Figure_S7_random_forest_Gini.jpg", Fig_S7, dpi=300, width=5, height=6, units="in")

#### FIG S7: RANDOM FOREST MISCLASSIFIED MODELS ####
rf.plot.b <- ggplot(rf.model.results) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                           size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("NTI", "PD", "Richness", "Vol. prey", "Spleen mass", "Shannon", "Health", "Age",
                            expression(paste(delta^{15}, "N")), expression(paste(delta^{13}, "C")), "Sex", "Vol. anthro")) +
  coord_flip() +
  theme_bw() +
  labs(x="Coefficient", y="Predictor", tag="b") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="none",
        axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

urban.data$VOL.PREY <- log(urban.data$VOL.PREY)
urban.data$Vol.Anthro <- log(urban.data$Vol.Anthro)
rf.differences.plot.data <- reshape2::melt(urban.data)
rf.differences.plot.data <- subset(rf.differences.plot.data, variable %in% c("NTI", "PD", "Observed", "VOL.PREY", "Spleen.to.Body.Ratio", "Shannon", "index1.PC1",
                                                 "Cementum.Age", "d15N", "d13C", "Sex", "Vol.Anthro"))

levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="index1.PC1"] <- "Health"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Cementum.Age"] <- "Age"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="VOL.PREY"] <- "Vol. prey"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Vol.Anthro"] <- "Vol. anthro"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Observed"] <- "ASV Richness"

rf.differences.plot.data$variable <- factor(rf.differences.plot.data$variable, levels=c("ASV Richness", "Shannon", "PD", "NTI",
                                                                                        "d13C", "d15N", "Vol. anthro", "Vol. prey",
                                                                                        "Age", "Health", "Spleen mass"))

rf.differences.plot.data$Equal <- as.factor(rf.differences.plot.data$Equal)
levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="FALSE"] <- "NO"
levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="TRUE"] <- "YES"

rf.plot.a <- ggplot(rf.differences.plot.data, aes(x=Equal, y=value)) + 
  geom_boxplot(aes(fill=Equal)) +
  scale_fill_manual(values=c("darkgoldenrod3", "purple")) +
  facet_wrap(~variable, scales="free") +
  theme_bw() +
  labs(x="\nCorrectly classified?", y="Value", tag="a") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        strip.background=element_blank(),
        panel.spacing=unit(0.5, "lines"),
        strip.text=element_text(color="black", size=9),
        legend.position="none",
        legend.title=element_text(color="black", size=10),
        legend.text=element_text(color="black", size=9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        plot.tag=element_text(face="bold"))

levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="NO"] <- "Incorrect"
levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="YES"] <- "Correct"

rf.plot.a.legend <- cowplot::get_legend(ggplot(rf.differences.plot.data, aes(x=Equal, y=value)) + 
                                          geom_boxplot(aes(fill=Equal)) +
                                          scale_fill_manual(values=c("darkgoldenrod3", "purple")) +
                                          guides(fill=guide_legend(title="Classification\naccuracy")) +
                                          facet_wrap(~variable, scales="free") +
                                          theme_bw() +
                                          labs(x="\nCorrectly classified?", y="Value", tag="a") +
                                          theme(panel.grid=element_blank(),
                                                axis.text=element_text(color="black", size=9),
                                                axis.title=element_text(color="black", size=10),
                                                strip.background=element_blank(),
                                                panel.spacing=unit(0.5, "lines"),
                                                strip.text=element_text(color="black", size=9),
                                                legend.position="right",
                                                legend.title=element_text(color="black", size=10),
                                                legend.text=element_text(color="black", size=9),
                                                panel.border = element_blank(),
                                                panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(),
                                                axis.line = element_line(colour = "black", size=0.5),
                                                plot.tag=element_text(face="bold")))

grid.arrange(rf.plot.a, rf.plot.b, ncol=1, heights=c(0.6,0.4))
  
Fig_S8 <- arrangeGrob(grobs=list(rf.plot.a, rf.plot.b), heights=c(0.6, 0.4), ncol=1)

ggsave("~/health/Tables/Figure_S8_rf_misclassified.jpg", Fig_S8, dpi=300, width=5, height=6.5, units="in")
ggsave("~/health/Tables/Figure_S8_legend.jpg", rf.plot.a.legend, dpi=300, width=1, height=1, units="in")


#### FIG S8: ORDINATION PANEL ####
ord.BC.NMDS <- metaMDS(dist.cy.list[["Bray-Curtis"]], trymax=1000, k=2)

temp.JAC <- as.data.frame(otu_table(pseq.rarefied))
temp.JAC <- temp.JAC %>% dplyr::mutate_if(is.numeric, ~1 * (. != 0))
temp.JAC <- vegdist(temp.JAC, method="jaccard")

ord.JAC.NMDS <- metaMDS(temp.JAC, trymax=1000, k=2)
ord.wUF.NMDS <- metaMDS(dist.cy.list[["wUF"]], trymax=1000, k=2)
ord.uwUF.NMDS <- metaMDS(dist.cy.list[["uwUF"]], trymax=1000, k=2)

ord.BC.envfit <- envfit(ord.BC.NMDS, 
                             temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                                "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                                "Em.PCR.Overall", "index3.PC1","Location")],
                             permutations=1000, na.rm=TRUE)
ord.JAC.envfit <- envfit(ord.JAC.NMDS, 
                        temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                           "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                           "Em.PCR.Overall", "index3.PC1","Location")],
                        permutations=1000, na.rm=TRUE)
ord.wUF.envfit <- envfit(ord.wUF.NMDS, 
                        temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                           "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                           "Em.PCR.Overall", "index3.PC1","Location")],
                        permutations=1000, na.rm=TRUE)
ord.uwUF.envfit <- envfit(ord.uwUF.NMDS, 
                        temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                           "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                           "Em.PCR.Overall", "index3.PC1","Location")],
                        permutations=1000, na.rm=TRUE)

scores.BC <- scores(ord.BC.NMDS, choices=c(1:3), display="sites")
scores.BC <- cbind(scores.BC, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.BC <- as.data.frame(scores(ord.BC.envfit, display="vectors"))
spp.scrs.BC <- cbind(spp.scrs.BC, Rsq = ord.BC.envfit[["vectors"]]$r,
                  pvals = ord.BC.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.BC))
spp.scrs.BC <- subset(spp.scrs.BC, pvals < 0.05)

scores.JAC <- scores(ord.JAC.NMDS, choices=c(1:3), display="sites")
scores.JAC <- cbind(scores.JAC, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.JAC <- as.data.frame(scores(ord.JAC.envfit, display="vectors"))
spp.scrs.JAC <- cbind(spp.scrs.JAC, Rsq = ord.JAC.envfit[["vectors"]]$r,
                     pvals = ord.JAC.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.JAC))
spp.scrs.JAC <- subset(spp.scrs.JAC, pvals < 0.05)

scores.wUF <- scores(ord.wUF.NMDS, choices=c(1:3), display="sites")
scores.wUF <- cbind(scores.wUF, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.wUF <- as.data.frame(scores(ord.wUF.envfit, display="vectors"))
spp.scrs.wUF <- cbind(spp.scrs.wUF, Rsq = ord.wUF.envfit[["vectors"]]$r,
                     pvals = ord.wUF.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.wUF))
spp.scrs.wUF <- subset(spp.scrs.wUF, pvals < 0.05)

scores.uwUF <- scores(ord.uwUF.NMDS, choices=c(1:3), display="sites")
scores.uwUF <- cbind(scores.uwUF, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.uwUF <- as.data.frame(scores(ord.uwUF.envfit, display="vectors"))
spp.scrs.uwUF <- cbind(spp.scrs.uwUF, Rsq = ord.uwUF.envfit[["vectors"]]$r,
                     pvals = ord.uwUF.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.uwUF))
spp.scrs.uwUF <- subset(spp.scrs.uwUF, pvals < 0.05)

plot.new()
temp <- ordiellipse(ord.BC.NMDS, temp.feces.data$Location,
                    display="sites", kind="sd", label=T, conf=0.95)
ord.BC.ellipse <- data.frame()
for(g in levels(temp.feces.data$Location)){
  ord.BC.ellipse <- rbind(ord.BC.ellipse, 
                                   cbind(as.data.frame(with(temp.data[temp.feces.data$Location==g,],
                                                            veganCovEllipse(temp[[g]]$cov,
                                                                            temp[[g]]$center,
                                                                            temp[[g]]$scale)))
                                         ,Type=g))
}
rm(temp)

plot.new()
temp <- ordiellipse(ord.JAC.NMDS, temp.feces.data$Location,
                    display="sites", kind="sd", label=T, conf=0.95)
ord.JAC.ellipse <- data.frame()
for(g in levels(temp.feces.data$Location)){
  ord.JAC.ellipse <- rbind(ord.JAC.ellipse, 
                          cbind(as.data.frame(with(temp.data[temp.feces.data$Location==g,],
                                                   veganCovEllipse(temp[[g]]$cov,
                                                                   temp[[g]]$center,
                                                                   temp[[g]]$scale)))
                                ,Type=g))
}
rm(temp)

plot.new()
temp <- ordiellipse(ord.wUF.NMDS, temp.feces.data$Location,
                    display="sites", kind="sd", label=T, conf=0.95)
ord.wUF.ellipse <- data.frame()
for(g in levels(temp.feces.data$Location)){
  ord.wUF.ellipse <- rbind(ord.wUF.ellipse, 
                          cbind(as.data.frame(with(temp.data[temp.feces.data$Location==g,],
                                                   veganCovEllipse(temp[[g]]$cov,
                                                                   temp[[g]]$center,
                                                                   temp[[g]]$scale)))
                                ,Type=g))
}
rm(temp)

plot.new()
temp <- ordiellipse(ord.uwUF.NMDS, temp.feces.data$Location,
                    display="sites", kind="sd", label=T, conf=0.95)
ord.uwUF.ellipse <- data.frame()
for(g in levels(temp.feces.data$Location)){
  ord.uwUF.ellipse <- rbind(ord.uwUF.ellipse, 
                          cbind(as.data.frame(with(temp.data[temp.feces.data$Location==g,],
                                                   veganCovEllipse(temp[[g]]$cov,
                                                                   temp[[g]]$center,
                                                                   temp[[g]]$scale)))
                                ,Type=g))
}
rm(temp)

levels(spp.scrs.BC$Variable)[levels(spp.scrs.BC$Variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"

plot.BC <- ggplot(scores.BC) + 
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.BC,
               aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.BC.ellipse, aes(x=NMDS1, y=NMDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", tag="a") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

plot.JAC <- ggplot(scores.JAC) + 
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.JAC,
               aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.JAC.ellipse, aes(x=NMDS1, y=NMDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", tag="b") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

plot.wUF <- ggplot(scores.wUF) + 
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.wUF,
               aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.wUF.ellipse, aes(x=NMDS1, y=NMDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", tag="c") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

plot.uwUF <- ggplot(scores.uwUF) + 
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.uwUF,
               aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.uwUF.ellipse, aes(x=NMDS1, y=NMDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", tag="d") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

grid.arrange(plot.BC, plot.JAC, plot.wUF, plot.uwUF, nrow=2, ncol=2)

levels(scores.BC$Location)[levels(scores.BC$Location)=="peri-urban"] <- "rural"
levels(ord.BC.ellipse$Type)[levels(ord.BC.ellipse$Type)=="peri-urban"] <- "rural"

plot.legend <- cowplot::get_legend(ggplot(scores.BC) + 
                                     geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Location, shape=Location)) + 
                                     geom_segment(data = spp.scrs.BC,
                                                  aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
                                                  arrow = arrow(length=unit(0.25, "cm")), color="black") +
                                     geom_path(data=ord.BC.ellipse, aes(x=NMDS1, y=NMDS2, colour=Type), size=0.5, linetype=1) +
                                     scale_color_manual(values=c("forestgreen", "purple")) +
                                     guides(color=guide_legend(title.position="top", title.hjust=0.5, ncol=2)) +
                                     scale_size_continuous(range=c(1,3)) +
                                     theme_bw() +
                                     labs(x="NMDS1", y="NMDS2", tag="a") +
                                     theme(axis.title=element_text(color="black", size=9),
                                           axis.text=element_text(color="black", size=8),
                                           plot.tag=element_text(face="bold"),
                                           panel.border=element_rect(color="black"),
                                           plot.title=element_text(size=10),
                                           legend.position="bottom",
                                           legend.title=element_text(size=10),
                                           legend.text=element_text(size=10),
                                           panel.grid.major=element_blank(),
                                           panel.grid.minor=element_blank(),
                                           plot.margin = unit(c(5.5,5.5,16,5.5), "pt")))


Fig_S9 <- arrangeGrob(grobs=list(plot.BC, plot.JAC, plot.wUF, plot.uwUF, plot.legend), layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
                      heights=c(0.4,0.4,0.1))
ggsave("~/health/Tables/Figure_S9_ordination_pane.jpg", Fig_S9, dpi=300, width=5, height=5, units="in")

rm(ord.BC.NMDS, ord.BC.envfit, ord.BC.ellipse, scores.BC, spp.scrs.BC,
   ord.JAC.NMDS, ord.JAC.envfit, ord.JAC.ellipse, scores.JAC, spp.scrs.JAC,
   ord.wUF.NMDS, ord.wUF.envfit, ord.wUF.ellipse, scores.wUF, spp.scrs.wUF,
   ord.uwUF.NMDS, ord.uwUF.envfit, ord.uwUF.ellipse, scores.uwUF, spp.scrs.uwUF)

#### FIG S9: SEX AND E. MULTI ORDINATIONS ####
ord.PCA.AIT <- rda(as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr"))))
scores <- scores(ord.PCA.AIT, choices=c(1:3), display="sites")
scores <- cbind(scores, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

# Confidence ellipse.
plot.new()
ellipse <- ordiellipse(ord.PCA.AIT, scores$Sex, display="sites", 
                       kind="sd", conf=0.95, label=T)
ellipse.data <- data.frame()
for(g in levels(scores$Sex)){
  ellipse.data <- rbind(ellipse.data, cbind(as.data.frame(with(scores[scores$Sex==g,],
                                                               veganCovEllipse(ellipse[[g]]$cov,
                                                                               ellipse[[g]]$center,
                                                                               ellipse[[g]]$scale)))
                                            ,Type=g))
}
rm(ellipse)

scores$Sex <- factor(scores$Sex, levels=c("M", "F"))

ord.plot.sex <- ggplot(scores) + 
  geom_point(mapping = aes(x=PC1, y=PC2, color=Sex, shape=Sex)) + 
  geom_path(data=ellipse.data, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("orange2", "maroon2")) +
  guides(color=guide_legend(ncol=2, title.position="top")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (11.0%)", y="PC2 (10.0%)", tag="a", shape="\nSex", color="\nSex") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        legend.position="bottom",
        legend.title=element_text(size=10, hjust=0.5),
        plot.title=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

plot.new()
ellipse <- ordiellipse(ord.PCA.AIT, scores$Em.PCR.Overall, display="sites", 
                       kind="sd", conf=0.95, label=T)
ellipse.data <- data.frame()
for(g in levels(scores$Em.PCR.Overall)){
  ellipse.data <- rbind(ellipse.data, cbind(as.data.frame(with(scores[scores$Em.PCR.Overall==g,],
                                                               veganCovEllipse(ellipse[[g]]$cov,
                                                                               ellipse[[g]]$center,
                                                                               ellipse[[g]]$scale)))
                                            ,Type=g))
}
rm(ellipse)

scores$Em.PCR.Overall <- factor(scores$Em.PCR.Overall, levels=c("Y", "N"))

ord.plot.em <- ggplot(scores) + 
  geom_point(mapping = aes(x=PC1, y=PC2, color=Em.PCR.Overall, shape=Em.PCR.Overall)) + 
  geom_path(data=ellipse.data, aes(x=PC1, y=PC2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("orangered1", "deepskyblue3")) +
  scale_size_continuous(range=c(1,3)) +
  guides(color=guide_legend(ncol=2, title.position="top")) +
  theme_bw() +
  labs(x="PC1 (11.0%)", y="PC2 (10.0%)", tag="b", shape="E.multilocularis\nInfection", color="E.multilocularis\nInfection") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=10, hjust=0.5),
        legend.position="bottom",
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

grid.arrange(ord.plot.sex, ord.plot.em, ncol=2)

Fig_S10 <- arrangeGrob(grobs=list(ord.plot.sex, ord.plot.em), ncol=2)
ggsave("~/health/Tables/Fig_S10_sex_em_ords.jpg", Fig_S10, dpi=300, width=6, height=4, units="in")

#### FIG S10, S11, S12: MODEL CONFIRMATIONS ####
rgcca.taxon.model.plots <- list()
for(i in 1:length(rgcca.taxon.models)){
  plot <- ggplot(rgcca.taxon.models[[i]]) +
    ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                              size=2, color="black") + 
    ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                              color="black", size=0.5) +
    geom_point(aes(y=merge, x=psd_Coef),
               fill="white", shape=21, stroke=1, show.legend=TRUE) +
    coord_flip() +
    theme_bw() +
    labs(x="Coefficient", y="Predictor", title=names(rgcca.taxon.models)[[i]]) +
    geom_vline(xintercept=0, linetype="dashed", color="black") +
    theme(plot.tag=element_text(face="bold"),
          panel.grid.major.y = element_line(linetype = "solid"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position="none",
          axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
          axis.text.y=element_text(size=9, color="black"),
          axis.title=element_text(size=11, color="black"),
          plot.title=element_text(size=9, color="black"))
  rgcca.taxon.model.plots[[i]] <- plot
  rm(plot)
}
names(rgcca.taxon.model.plots) <- names(rgcca.taxon.models)

Fig_S11_anthro.taxa <- arrangeGrob(grobs=list(rgcca.taxon.model.plots[["Streptococcus"]] +
                                                scale_y_discrete(labels=c('Vol. anthro', 'Spleen mass', 'Age', 'Health', 'E. multi', 
                                                                          'Sex', expression(paste(delta^{13}, "C")),
                                                                          expression(paste(delta^{15}, "N")), 'Vol. prey')) +
                                                labs(x=NULL, y=NULL),
                                              rgcca.taxon.model.plots[["Enterococcus"]] +
                                                scale_y_discrete(labels=c('Age', 'Spleen mass', 'Vol. anthro', 'E. multi',
                                                                          'Health', expression(paste(delta^{15}, "N")), "Vol. prey",
                                                                          expression(paste(delta^{13}, "C")), 'Sex')) +
                                                labs(x=NULL, y=NULL),
                                              rgcca.taxon.model.plots[["Lactobacillus"]] +
                                                scale_y_discrete(labels=c('Vol. anthro', 'Age', 'Sex', 'Vol. prey',
                                                                          'Spleen mass', expression(paste(delta^{13}, "C")), 'Health',
                                                                          expression(paste(delta^{15}, "N")), 'E. multi')) +
                                                labs(x=NULL, y=NULL)),
                                   ncol=3)


Fig_S12_protein.taxa <- arrangeGrob(grobs=list(rgcca.taxon.model.plots[["Fusobacterium"]] +
                                                 scale_y_discrete(labels=c(expression(paste(delta^{15}, "N")), 'Vol. anthro', 'Vol. prey', 
                                                                           expression(paste(delta^{13}, "C")), 'Age', 'Spleen mass',
                                                                           'Health', 'E. multi', 'Sex'))+
                                                 labs(x=NULL, y=NULL),
                                               rgcca.taxon.model.plots[["Anaerobiospirillum"]] +
                                                 scale_y_discrete(labels=c('Vol. anthro', expression(paste(delta^{15}, "N")), 'Vol. prey',
                                                                           expression(paste(delta^{13}, "C")), "Age", 
                                                                           "Sex", "Spleen mass", "E. multi", "Health"))+
                                                 labs(x=NULL, y=NULL),
                                               rgcca.taxon.model.plots[["Sutterella"]] +
                                                 scale_y_discrete(labels=c('Vol. prey',  expression(paste(delta^{15}, "N")), "Spleen mass",
                                                                           "Sex", expression(paste(delta^{13}, "C")), "Vol. anthro",
                                                                           "Age", "Health", "E. multi"))+
                                                 labs(x=NULL, y=NULL),
                                               rgcca.taxon.model.plots[["Bacteroides"]] +
                                                 scale_y_discrete(labels=c(expression(paste(delta^{15}, "N")), "Sex", "Spleen mass", "Age", "Vol. prey",
                                                                           "Health", "Vol. anthro", "E. multi", expression(paste(delta^{13}, "C"))))+
                                                 labs(x=NULL, y=NULL),
                                               rgcca.taxon.model.plots[["Alloprevotella"]] +
                                                 scale_y_discrete(labels=c('Age',  'E. multi', 'Spleen mass', "Vol. anthro", 
                                                                           expression(paste(delta^{13}, "C")), "Health",
                                                                           expression(paste(delta^{15}, "N")), "Sex", "Vol. prey"))+
                                                 labs(x=NULL, y=NULL)),
                                    layout_matrix=rbind(c(1,1,2,2,3,3),c(NA,4,4,5,5,NA)))

Fig_S13_spleen.taxa <- arrangeGrob(grobs=list(rgcca.taxon.model.plots[["Erysipelotrichaceae"]] + labs(x=NULL, y=NULL) +
                                                scale_y_discrete(labels=c(expression(paste(delta^{15}, "N")), 'Spleen mass',
                                                                          expression(paste(delta^{13}, "C")), 'Health', 'Vol. prey',
                                                                          'Vol. anthro', 'Age', 'E. multi', 'Sex')),
                                              rgcca.taxon.model.plots[["Coriobacteriaceae"]] + labs(x=NULL, y=NULL) +
                                                scale_y_discrete(labels=c('Spleen mass', expression(paste(delta^{15}, "N")), 
                                                                          expression(paste(delta^{13}, "C")), "Vol. prey", "Age", 
                                                                          "Health", "Vol. anthro", "Sex", "E. multi")),
                                              rgcca.taxon.model.plots[["Lachnospiraceae"]] + labs(x=NULL, y=NULL) +
                                                scale_y_discrete(labels=c('Spleen mass',  expression(paste(delta^{13}, "C")), 
                                                                          expression(paste(delta^{15}, "N")), "Health", "Vol. prey",
                                                                          "E. multi", "Age", "Sex", "Vol. anthro"))),
                                   ncol=3)


ggsave("~/health/Tables/Figure_S11_anthro_prediction_models.jpg", Fig_S11_anthro.taxa, dpi=300, height=3.25, width=9, units="in")
ggsave("~/health/Tables/Figure_S12_protein_prediction_models.jpg", Fig_S12_protein.taxa, dpi=300, height=6.5, width=9, units="in")
ggsave("~/health/Tables/Figure_S13_spleen_prediction_models.jpg", Fig_S13_spleen.taxa, dpi=300, height=3.25, width=9, units="in")

#### FIG S14: COLLINEARITY IN AGE AND HEALTH ####
plot.data.temp <- cy.metadata
plot.data.temp$Sex <- factor(plot.data.temp$Sex, levels=c("M","F"))

Fig_S14 <- ggplot(plot.data.temp, aes(x=Cementum.Age, y=index1.PC1, color=Sex)) + geom_point() + geom_smooth(method="lm", se=FALSE) +
  theme_bw() +
  scale_color_manual(values=c("orange2", "maroon2")) +
  scale_y_continuous(expand=c(0.05,0.05)) +
  scale_x_continuous(expand=c(0.05,0.05)) +
  labs(x="Age", y="Health score") +
  theme(axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        legend.title=element_text(size=10, color="black"),
        panel.grid=element_blank())

ggsave("~/health/Tables/Figure_S14_sex_and_age.jpg", Fig_S14, dpi=300, width=3.5, height=3, units="in")