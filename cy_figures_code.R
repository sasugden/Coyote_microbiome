##### FIGURE 1: DIET AND HEALTH RESPONSES ####
diet.figure <- cy.metadata[,c("Location","d13C","d15N","Vol.Anthro","VOL.PREY")]

# (a) Stable isotope isoscape ####
diet.figure <- cy.metadata[,c("FecesID", "Location", "d13C", "d15N", "Vol.Anthro", "VOL.PREY")]

diet.isoscape <- cy.metadata %>%
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
simmr_plot <- reshape2::melt(mixing.model.results)
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
health.figure.data <- cy.metadata[,c("Location", "KFI", "Spleen.to.Body.Ratio","index1.PC1", "index1.PC1.z", "Em.PCR.Overall")]
health.figure.data <- melt(health.figure.data)
levels(health.figure.data$variable)[levels(health.figure.data$variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass\n(g/kg)"
levels(health.figure.data$variable)[levels(health.figure.data$variable)=="index1.PC1"] <- "Health score"
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
  labs(x="\nLocation", y=expression(paste(italic("E. multi"), " prevalence (%)")), tag="e") +
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
ggsave("~/health/Figure_1_diet_health.jpg", Figure_1_diet_health, dpi=300, width=6.5, height=5, units="in")


#### FIGURE 2: MICROBIOME STRUCTURE & COMPOSITION ####
# (a) Urban/rural alpha diversity boxplots ####
adiv.urb.rur <- cy.metadata[,c("Location","Observed_Extrap","Shannon_Extrap","PD_Extrap","NTI_Raw")]
adiv.urb.rur <- reshape2::melt(adiv.urb.rur)

levels(adiv.urb.rur$variable)[levels(adiv.urb.rur$variable)=="Observed_Extrap"] <- "Richness"
levels(adiv.urb.rur$variable)[levels(adiv.urb.rur$variable)=="Shannon_Extrap"] <- "Shannon"
levels(adiv.urb.rur$variable)[levels(adiv.urb.rur$variable)=="PD_Extrap"] <- "PD"
levels(adiv.urb.rur$variable)[levels(adiv.urb.rur$variable)=="NTI_Raw"] <- "NTI"

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
cor.feces <- Hmisc::rcorr(as.matrix(cy.metadata[,sapply(cy.metadata, is.numeric)]), type="spearman")

for(i in 1:length(cor.feces)){
  cor.feces[[i]] <- subset(cor.feces[[i]], colnames(cor.feces[[i]]) %in% c("index1.PC1", "Cementum.Age", "Spleen.to.Body.Ratio", "d13C", "d15N", "Vol.Anthro","VOL.PREY"))
  cor.feces[[i]] <- as.data.frame(t(cor.feces[[i]]))
  cor.feces[[i]] <- subset(cor.feces[[i]], rownames(cor.feces[[i]]) %in% c("Observed_Extrap", "Shannon_Extrap", "PD_Extrap"))
  cor.feces[[i]]$Measure <- rownames(cor.feces[[i]])
}


cor.feces.R.melt <- reshape2::melt(cor.feces[["r"]])
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="index1.PC1"] <- "Health index"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="Cementum.Age"] <- "Age"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="Vol.Anthro"] <- "Vol. anthro"
levels(cor.feces.R.melt$variable)[levels(cor.feces.R.melt$variable)=="VOL.PREY"] <- "Vol. prey"
cor.feces.R.melt$Measure[cor.feces.R.melt$Measure=="Observed_Extrap"] <- "Richness"
cor.feces.R.melt$Measure[cor.feces.R.melt$Measure=="Shannon_Extrap"] <- "Shannon"
cor.feces.R.melt$Measure[cor.feces.R.melt$Measure=="PD_Extrap"] <- "PD"
cor.feces.R.melt$Measure <- factor(cor.feces.R.melt$Measure, levels=c("PD", "Shannon", "Richness"))

cor.feces.R.melt$value[cor.feces.R.melt$value > 0.2] <- 0.2

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

# (d) Differential abundance ####
diff.abund.plot.data <- diff.abund.location[[5]]
diff.abund.plot.data <- subset(diff.abund.plot.data, Genus %in% rownames(rgcca.scores))
diff.abund.plot.data[,c(1:14,19,21)] <- NULL
diff.abund.plot.data <- diff.abund.plot.data[order(diff.abund.plot.data$X7),]

diff.abund.plot.data$p.adj <- p.adjust(diff.abund.plot.data$model.Locationurban.Pr...t.., method="BH")

diff.abund.plot.data <- subset(diff.abund.plot.data, abs(X7) > 0.3)

diff.abund.plot.data$Genus <- factor(diff.abund.plot.data$Genus)
diff.abund.plot.data$Genus <- forcats::fct_inorder(diff.abund.plot.data$Genus)

abund.means <- prev.data[["Genus"]][,c("Taxon", "Abundance")]
colnames(abund.means) <- c("Genus", "Abundance")
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
temp.feces.data <- subset(cy.metadata, FecesID %in% sample_names(pseq.))

ord.PCA.AIT <- rda(as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr"))))
ord.PCA.AIT.envfit <- envfit(ord.PCA.AIT, 
                             temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                                "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                                "Em.PCR.Overall", "index1.PC1","Location")],
                             permutations=1000, na.rm=TRUE)

scores <- as.data.frame(scores(ord.PCA.AIT, choices=c(1:3), display="sites"))
scores$FecesID <- rownames(scores)
scores <- merge(scores, cy.metadata[,c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")], by="FecesID", all=FALSE)

spp.scrs <- as.data.frame(scores(ord.PCA.AIT.envfit, display="vectors"))
spp.scrs <- cbind(spp.scrs, Rsq = ord.PCA.AIT.envfit[["vectors"]]$r,
                  pvals = ord.PCA.AIT.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs))
spp.scrs <- subset(spp.scrs, pvals < 0.05)
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="Cementum.Age"] <- "Age"
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="Vol.Anthro"] <- "Vol. anthro"
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="index1.PC1"] <- "Health"
levels(spp.scrs$Variable)[levels(spp.scrs$Variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"

# Confidence ellipse.
plot.new()
ellipse <- ordiellipse(ord.PCA.AIT, 
                       temp.feces.data$Location, 
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
             plot.diff.abund.genus+theme(legend.position="none")+labs(tag="d"),
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
temp1 <- subset(temp, rownames(temp) %in% colnames(feces.strict.sum.clr[[5]]))
temp2 <- subset(temp, rownames(temp) %in% c("Mass","Snout.to.Base","Girth",
                                            "Cementum.Age","Spleen.to.Body.Ratio","KFI",
                                            "d13C", "d15N", "VOL.PREY", "Vol.Anthro", "VOL.VEG"))

temp <- rbind(temp1, temp2)
rm(temp1, temp2)
temp <- as.data.frame(t(temp))

# Calculate correlations between each taxa and metadata variable.
temp.spearman <- Hmisc::rcorr(as.matrix(temp[,sapply(temp, is.numeric)]), type="spearman")
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
tax.table <- as.data.frame(cbind(tax_table(tax_glom(subset_taxa(pseq.clr.sub, Genus %in% temp.spearman.melt$Var2), "Genus"))))
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
rm(temp, tax.table, temp.spearman, temp.spearman.R, desired_order)

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
  geom_errorbarh(aes(xmin=Lower, xmax=Upper, y=Taxon), position=position_dodge(0.9), height=0, color="darkgrey") +
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
#### FIGURE S1: AGE VS. E MULTI INFECTION, FACET BY LOCATION ####
supp_plot_data <- cy.metadata

levels(supp_plot_data$PA.Empty)[levels(supp_plot_data$PA.Empty)==1] <- "Y"
levels(supp_plot_data$PA.Empty)[levels(supp_plot_data$PA.Empty)==0] <- "N"
levels(supp_plot_data$Location)[levels(supp_plot_data$Location)=="peri-urban"] <- "rural"

Fig_S1 <- ggplot(subset(supp_plot_data, is.na(Em.PCR.Overall)=="FALSE"), aes(x=Em.PCR.Overall, y=Cementum.Age)) +
  geom_boxplot(aes(color=Location, fill=Em.PCR.Overall)) +
  facet_grid(~Location) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_fill_manual(values=c("antiquewhite", "darkgrey")) +
  guides(color=guide_legend(title.position="top", order=2), 
         fill=guide_legend(title.position="top", order=1,
                           title=expression(paste(italic("E. multilocularis")))),
         order=1) +
  labs(x=expression(paste(italic("E. multilocularis")," infection (Y/N)")), y="Age (yr)") +
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


supp_plot_data$Test <- paste0(supp_plot_data$Em.PCR.Overall, supp_plot_data$Location)

summary(aov(Observed_Extrap ~ Test, supp_plot_data))
TukeyHSD(aov(Observed_Extrap ~ Test, supp_plot_data))

summary(aov(Shannon_Extrap ~ Test, supp_plot_data))
TukeyHSD(aov(Shannon_Extrap ~ Test, supp_plot_data))

summary(aov(PD_Extrap ~ Test, supp_plot_data))
TukeyHSD(aov(PD_Extrap ~ Test, supp_plot_data))

summary(aov(NTI_Raw ~ Test, supp_plot_data))
TukeyHSD(aov(NTI_Raw ~ Test, supp_plot_data))

supp_plot_data <- subset(supp_plot_data, is.na(Em.PCR.Overall)=="FALSE")
summary(aov(Cementum.Age ~ Test, supp_plot_data))
TukeyHSD(aov(Cementum.Age ~ Test, supp_plot_data))

rm(supp_plot_data)

ggsave("~/health/Figure_S1_age_Emulti.jpg", Fig_S1, dpi=300, width=3, height=3.5, units="in")




#### FIG S2: ALPHA DIVERSITY MEASURES VS. EM, CONTROLLING FOR LOCATION ####
supp_plot_data <- reshape2::melt(supp_plot_data)
supp_plot_data <- subset(supp_plot_data, variable %in% c("Observed_Extrap", "Shannon_Extrap", "PD_Extrap", "NTI_Raw"))
supp_plot_data$variable <- factor(supp_plot_data$variable, levels=c("Observed_Extrap", "Shannon_Extrap", "PD_Extrap", "NTI_Raw"))

levels(supp_plot_data$variable)[levels(supp_plot_data$variable)=="Observed_Extrap"] <- "ASV Richness"
levels(supp_plot_data$variable)[levels(supp_plot_data$variable)=="Shannon_Extrap"] <- "Shannon"
levels(supp_plot_data$variable)[levels(supp_plot_data$variable)=="PD_Extrap"] <- "PD"
levels(supp_plot_data$variable)[levels(supp_plot_data$variable)=="NTI_Raw"] <- "NTI"

colnames(supp_plot_data)[colnames(supp_plot_data)=="Em.PCR.Overall"] <- "E. multi"

Fig_S2 <- ggplot(supp_plot_data, aes(x=Location, y=value, color=Location, fill=`E. multi`)) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free", nrow=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  guides(fill=guide_legend(title=expression(paste(italic("E. multi"))))) +
  scale_fill_manual(values=c("antiquewhite", "darkgrey")) +
  labs(x=expression(paste("\nLocation / ", italic("E. multilocularis"), " infection")), y="Value") +
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

ggsave("~/health/Figure_S2_Em_alpha_diversity.jpg", Fig_S2, dpi=300, width=6, height=3, units="in")

#### FIG S3: ALPHA DIVERSITY DREDGED MODEL PLOTS ####
adiv.model.plots <- list()
adiv.model.plots[["Observed"]] <- ggplot(adiv.model.results[["Observed_Extrap"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("E. multi", "Location", "Age", expression(paste(delta^{15}, "N")), "Health", 
                            "Spleen mass", "Vol. prey", "Vol. anthro", "Sex", expression(paste(delta^{13}, "C")))) +
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

adiv.model.plots[["Shannon"]] <- ggplot(adiv.model.results[["Shannon_Extrap"]]) +
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

adiv.model.plots[["PD"]] <- ggplot(adiv.model.results[["PD_Extrap"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Location", "E. multi", "Spleen mass", "Health", "Age", "Vol. anthro", expression(paste(delta^{13}, "C")),
                            expression(paste(delta^{15}, "N")), "Sex", "Vol. prey")) +
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

adiv.model.plots[["NTI"]] <- ggplot(adiv.model.results[["NTI_Raw"]]) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Location", "Sex", "Health", "Spleen mass", "E. multi", "Age", "Vol. anthro", "Vol. prey",
                            expression(paste(delta^{13}, "C")), expression(paste(delta^{15}, "N")))) +
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

ggsave("~/health/Figure_S4_adiv_dredged_models.jpg", Fig_S4, dpi=300, width=6.5, height=6.5, units="in")

#### FIG S4: ALPHA DIVERSITY AND EMPTY STOMACHS ####
Fig_S4 <- ggplot(supp_plot_data, aes(x=Location, y=value, color=Location, fill=PA.Empty)) +
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

ggsave("~/health/Figure_S4_empty_stomach_alpha_diversity.jpg", Fig_S4, dpi=300, width=6, height=3, units="in")

supp_plot_data <- cy.metadata
levels(supp_plot_data$PA.Empty)[levels(supp_plot_data$PA.Empty)==1] <- "Y"
levels(supp_plot_data$PA.Empty)[levels(supp_plot_data$PA.Empty)==0] <- "N"
levels(supp_plot_data$Location)[levels(supp_plot_data$Location)=="peri-urban"] <- "rural"

supp_plot_data$Test <- paste0(supp_plot_data$Location, supp_plot_data$PA.Empty)

summary(aov(Observed_Extrap ~ Test, supp_plot_data))
TukeyHSD(aov(Observed_Extrap ~ Test, supp_plot_data))

summary(aov(Shannon_Extrap ~ Test, supp_plot_data))
TukeyHSD(aov(Shannon_Extrap ~ Test, supp_plot_data))

summary(aov(PD_Extrap ~ Test, supp_plot_data))
TukeyHSD(aov(PD_Extrap ~ Test, supp_plot_data))

summary(aov(NTI_Raw ~ Test, supp_plot_data))
TukeyHSD(aov(NTI_Raw ~ Test, supp_plot_data))

#### FIG S5: DIFFERENTIAL ABUNDANCE BASED ON AGE, SEX, ETC. ####
# (a) Sex only (controlled for age) ####
diff.abund.supp.data <- diff.abund.just.sex.and.age[[5]]
diff.abund.supp.data <- subset(diff.abund.supp.data, Genus %in% rownames(rgcca.scores))
diff.abund.supp.data <- diff.abund.supp.data[order(diff.abund.supp.data$X7),]

diff.abund.supp.data <- subset(diff.abund.supp.data, abs(X7) > 0.25)

diff.abund.supp.data$Genus <- factor(diff.abund.supp.data$Genus)
diff.abund.supp.data$Genus <- forcats::fct_inorder(diff.abund.supp.data$Genus)

abund.means <- prev.data[["Genus"]][,c("Taxon", "Abundance")]
colnames(abund.means) <- c("Genus", "Abundance")
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
diff.abund.supp.data.em <- subset(diff.abund.supp.data.em, Genus %in% rownames(rgcca.scores))
diff.abund.supp.data.em <- diff.abund.supp.data.em[order(diff.abund.supp.data.em$X7),]

diff.abund.supp.data.em <- subset(diff.abund.supp.data.em, abs(X7) > 0.25)

diff.abund.supp.data.em$Genus <- factor(diff.abund.supp.data.em$Genus)
diff.abund.supp.data.em$Genus <- forcats::fct_inorder(diff.abund.supp.data.em$Genus)

abund.means <- prev.data[["Genus"]][,c("Taxon", "Abundance")]
colnames(abund.means) <- c("Genus", "Abundance")
diff.abund.supp.data.em <- merge(diff.abund.supp.data.em, abund.means, by="Genus", all=FALSE)

diff.abund.supp.data.em$p.adj <- p.adjust(diff.abund.supp.data.em$model.Em.PCR.OverallY.Pr...t.., method="BH")
diff.abund.supp.data.em$shape <- cut(diff.abund.supp.data.em$p.adj,
                                     breaks=c(-Inf,0.1,Inf), labels=c("significant","not"))

diff.abund.supp.data.em$Phylum <- as.character(diff.abund.supp.data.em$Phylum)
diff.abund.supp.data.em$Class <- as.character(diff.abund.supp.data.em$Class)

diff.abund.supp.data.em$Phylum <- as.factor(diff.abund.supp.data.em$Phylum)

plot.diff.abund.em <- ggplot(diff.abund.supp.data.em, aes(x=X7, y=Genus)) + 
  geom_point(aes(size=Abundance, color=Phylum, shape=shape)) +
  scale_color_manual(values=myColors2) +
  scale_shape_manual(values=c(19)) +
  scale_size_continuous(range=c(2,6)) +
  geom_vline(xintercept=0, color="grey") +
  geom_vline(xintercept=0.2, color="grey", linetype=2) +
  geom_vline(xintercept=-0.2, color="grey", linetype=2) +
  theme_bw() +
  labs(x="Hedge's g", y="Genus", tag="b", title=expression(paste(italic("E. multilocularis")))) +
  scale_y_discrete(position="left") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=8),
        axis.title=element_text(color="black", size=10),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=10, color="black"),
        legend.position="none")

# (c) Age (controlled for sex) ####
diff.abund.supp.data.age <- diff.abund.just.sex.and.age[[5]]
diff.abund.supp.data.age <- subset(diff.abund.supp.data.age, Genus %in% rownames(rgcca.scores))

abund.means <- prev.data[["Genus"]][,c("Taxon", "Abundance")]
colnames(abund.means) <- c("Genus", "Abundance")

diff.abund.supp.data.age <- merge(diff.abund.supp.data.age, abund.means, by="Genus", all=FALSE)

temp.corr <- Hmisc::rcorr(as.matrix(comp.strict.clr$Genus[,sapply(comp.strict.clr$Genus, is.numeric)]), type="spearman")
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
diff.abund.supp.data.age$Phylum <- factor(diff.abund.supp.data.age$Phylum)

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

Fig_S5_legend <- cowplot::get_legend(ggplot(diff.abund.legend.data, aes(x=X7, y=Genus)) + 
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
                                             legend.position="right",
                                             legend.box="horizontal"))


grid.arrange(plot.diff.abund.sex, plot.diff.abund.em, plot.diff.abund.age, layout_matrix=rbind(c(1,1,2,2),c(NA,3,3,NA)))

Fig_S5 <- arrangeGrob(grobs=list(plot.diff.abund.sex, plot.diff.abund.em, plot.diff.abund.age, Fig_S5_legend),
                      layout_matrix=rbind(c(1,1,2,2),c(3,3,4,4)))

ggsave("~/health/Figure_S5_differential_abundance_tests.jpg", Fig_S5, dpi=300, width=6.5, height=6, units="in")

rm(diff.abund.supp.data, diff.abund.supp.data.age, diff.abund.supp.data.em, diff.abund.supp.data.health,
   temp.corr)

#### FIG S6: RANDOM FOREST RESULTS ####
# Sort graph results by decreasing Gini scores.
rf.plot.data <- rf.importance[order(-rf.importance$MeanDecreaseGini),]

# Insert ASVID as a row name.
rf.plot.data$ASVID <- rownames(rf.plot.data)

# Order plots by ASVID.
rf.plot.data$ASVID <- forcats::fct_inorder(rf.plot.data$ASVID)

# Determine which grouping (segment: small/large, individual: identity, site: intestinal site) each taxon is most abundant in.
rf.plot.data$Max <- as.factor(colnames(rf.plot.data[,8:9])[max.col(rf.plot.data[,8:9], ties.method="first")])
rf.plot.data$Max <- factor(rf.plot.data$Max, levels=c("rural","urban"))

rf.plot.data <- rf.plot.data[c(1:20),]
rf.plot.data$ASVID <- factor(rf.plot.data$ASVID)

levels(rf.plot.data$Genus)[levels(rf.plot.data$Genus)=="Clostridium_sensu_stricto"] <- "Clostridium"

Fig_S6 <- ggplot(rf.plot.data, aes(x=MeanDecreaseGini, y=ASVID)) + 
  geom_point(aes(color=Max, size=Overall)) + 
  scale_y_discrete(limits=rev(levels(rf.plot.data$ASVID)), labels=rev(as.character(rf.plot.data$Genus))) +
  scale_color_manual(name="Location", values=c("forestgreen", "purple")) +
  scale_size_continuous(name="Abundance (%)", range=c(3,6)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.5, label="Model accuracy:\n86.4%", size=3.75) +
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

ggsave("~/health/FigS6_random_forest_Gini.jpg", Fig_S6, dpi=300, width=5, height=6, units="in")

#### FIG S7: RANDOM FOREST MISCLASSIFIED MODELS ####
# (a) Comparison of correctly and incorrectly classified
# Reproduce data frame of urban data.
rf.test <- cbind(subset(cy.metadata, FecesID %in% sample_names(pseq.clr)), misclassified$Equal)
colnames(rf.test) <- c(colnames(rf.test)[1:(ncol(rf.test)-1)], "Equal")
rownames(rf.test) <- rf.test$FecesID
urban.data <- subset(rf.test, Microbiome=="Y" & Location=="urban")
urban.data <- urban.data[,c("Equal","Sex","Em.PCR.Overall","Mass","Snout.to.Base", "Girth", "KFI","index1.PC1",
                            "Spleen.to.Body.Ratio", "Cementum.Age","d13C","d15N","VOL.PREY","Vol.Anthro",
                            "Observed", "Observed_Extrap", "Shannon", "Shannon_Extrap", "PD", "PD_Extrap", "NTI", "NTI_Raw")]

urban.data$VOL.PREY <- log(urban.data$VOL.PREY+0.01)
urban.data$Vol.Anthro <- log(urban.data$Vol.Anthro+0.01)

rf.differences.plot.data <- reshape2::melt(urban.data)
rf.differences.plot.data <- subset(rf.differences.plot.data, variable %in% c("NTI_Raw", "PD_Extrap", "Observed_Extrap",
                                                                             "VOL.PREY", "Spleen.to.Body.Ratio",
                                                                             "Shannon_Extrap", "index1.PC1",
                                                                             "Cementum.Age", "d15N", "d13C", 
                                                                             "Sex", "Vol.Anthro"))

levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="index1.PC1"] <- "Health"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Spleen.to.Body.Ratio"] <- "Spleen mass"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Cementum.Age"] <- "Age"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="VOL.PREY"] <- "Vol. prey"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Vol.Anthro"] <- "Vol. anthro"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Observed_Extrap"] <- "ASV Richness"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="Shannon_Extrap"] <- "Shannon"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="PD_Extrap"] <- "PD"
levels(rf.differences.plot.data$variable)[levels(rf.differences.plot.data$variable)=="NTI_Raw"] <- "NTI"

rf.differences.plot.data$variable <- factor(rf.differences.plot.data$variable, levels=c("ASV Richness", "Shannon", "PD", "NTI",
                                                                                        "d13C", "d15N", "Vol. anthro", "Vol. prey",
                                                                                        "Age", "Health", "Spleen mass"))

rf.differences.plot.data$Equal <- as.factor(rf.differences.plot.data$Equal)
levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="FALSE"] <- "N"
levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="TRUE"] <- "Y"

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

levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="NO"] <- "Incorrect (rural)"
levels(rf.differences.plot.data$Equal)[levels(rf.differences.plot.data$Equal)=="YES"] <- "Correct (urban)"

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

rf.plot.b <- ggplot(rf.model.results) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  scale_y_discrete(labels=c("Vol. prey", "Sex", "ASV Richness", expression(paste(delta^{15}, "N")),
                            "Age", "Health", "Shannon", "Vol. anthro", "Spleen mass", expression(paste(delta^{13}, "C")))) +
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

grid.arrange(rf.plot.a, rf.plot.b, ncol=1, heights=c(0.6,0.4))

Fig_S7 <- arrangeGrob(grobs=list(rf.plot.a, rf.plot.b), heights=c(0.6, 0.4), ncol=1)

ggsave("~/health/Figure_S7_rf_misclassified.jpg", Fig_S7, dpi=300, width=5, height=6.5, units="in")
ggsave("~/health/Figure_S7_legend.jpg", rf.plot.a.legend, dpi=300, width=1, height=1, units="in")


#### FIG S8: ORDINATION PANEL ####
cy.metadata.pseq <- subset(cy.metadata, FecesID %in% sample_names(pseq.clr))

ord.BC.PCoA <- capscale(dist.cy.list[["Bray-Curtis"]] ~ 1)

temp.JAC <- as.data.frame(otu_table(pseq.rarefied))
temp.JAC <- temp.JAC %>% dplyr::mutate_if(is.numeric, ~1 * (. != 0))
temp.JAC <- vegdist(temp.JAC, method="jaccard")

ord.JAC.PCoA <- capscale(temp.JAC ~ 1)
ord.wUF.PCoA <- capscale(dist.cy.list[["wUF"]] ~ 1)
ord.uwUF.PCoA <- capscale(dist.cy.list[["uwUF"]] ~ 1)

ord.BC.envfit <- envfit(ord.BC.PCoA, 
                        temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                           "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                           "Em.PCR.Overall", "index1.PC1","Location")],
                        permutations=1000, na.rm=TRUE)
ord.JAC.envfit <- envfit(ord.JAC.PCoA, 
                         temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                            "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                            "Em.PCR.Overall", "index1.PC1","Location")],
                         permutations=1000, na.rm=TRUE)
ord.wUF.envfit <- envfit(ord.wUF.PCoA, 
                         temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                            "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                            "Em.PCR.Overall", "index1.PC1","Location")],
                         permutations=1000, na.rm=TRUE)
ord.uwUF.envfit <- envfit(ord.uwUF.PCoA, 
                          temp.feces.data[,c("Sex","Cementum.Age","Spleen.to.Body.Ratio",
                                             "Vol.Anthro", "VOL.PREY", "d15N", "d13C",
                                             "Em.PCR.Overall", "index1.PC1","Location")],
                          permutations=1000, na.rm=TRUE)

scores.BC <- scores(ord.BC.PCoA, choices=c(1:3), display="sites")
scores.BC <- cbind(scores.BC, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.BC <- as.data.frame(scores(ord.BC.envfit, display="vectors"))
spp.scrs.BC <- cbind(spp.scrs.BC, Rsq = ord.BC.envfit[["vectors"]]$r,
                     pvals = ord.BC.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.BC))
spp.scrs.BC <- subset(spp.scrs.BC, pvals < 0.05)

scores.JAC <- scores(ord.JAC.PCoA, choices=c(1:3), display="sites")
scores.JAC <- cbind(scores.JAC, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.JAC <- as.data.frame(scores(ord.JAC.envfit, display="vectors"))
spp.scrs.JAC <- cbind(spp.scrs.JAC, Rsq = ord.JAC.envfit[["vectors"]]$r,
                      pvals = ord.JAC.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.JAC))
spp.scrs.JAC <- subset(spp.scrs.JAC, pvals < 0.05)

scores.wUF <- scores(ord.wUF.PCoA, choices=c(1:3), display="sites")
scores.wUF <- cbind(scores.wUF, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.wUF <- as.data.frame(scores(ord.wUF.envfit, display="vectors"))
spp.scrs.wUF <- cbind(spp.scrs.wUF, Rsq = ord.wUF.envfit[["vectors"]]$r,
                      pvals = ord.wUF.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.wUF))
spp.scrs.wUF <- subset(spp.scrs.wUF, pvals < 0.05)

scores.uwUF <- scores(ord.uwUF.PCoA, choices=c(1:3), display="sites")
scores.uwUF <- cbind(scores.uwUF, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID")])

spp.scrs.uwUF <- as.data.frame(scores(ord.uwUF.envfit, display="vectors"))
spp.scrs.uwUF <- cbind(spp.scrs.uwUF, Rsq = ord.uwUF.envfit[["vectors"]]$r,
                       pvals = ord.uwUF.envfit[["vectors"]]$pvals, Variable = rownames(spp.scrs.uwUF))
spp.scrs.uwUF <- subset(spp.scrs.uwUF, pvals < 0.05)

plot.new()
temp <- ordiellipse(ord.BC.PCoA, temp.feces.data$Location,
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
temp <- ordiellipse(ord.JAC.PCoA, temp.feces.data$Location,
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
temp <- ordiellipse(ord.wUF.PCoA, temp.feces.data$Location,
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
temp <- ordiellipse(ord.uwUF.PCoA, temp.feces.data$Location,
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

summary(ord.BC.PCoA)

plot.BC <- ggplot(scores.BC) + 
  geom_point(mapping = aes(x=MDS1, y=MDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.BC,
               aes(x=0, y=0, xend=3*MDS1, yend=3*MDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.BC.ellipse, aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (17.3%)", y="PC2 (12.7%)", tag="a") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

summary(ord.JAC.PCoA)

plot.JAC <- ggplot(scores.JAC) + 
  geom_point(mapping = aes(x=MDS1, y=MDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.JAC,
               aes(x=0, y=0, xend=3*MDS1, yend=3*MDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.JAC.ellipse, aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (8.6%)", y="PC2 (6.7%)", tag="b") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

summary(ord.wUF.PCoA)

plot.wUF <- ggplot(scores.wUF) + 
  geom_point(mapping = aes(x=MDS1, y=MDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.wUF,
               aes(x=0, y=0, xend=3*MDS1, yend=3*MDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.wUF.ellipse, aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (32.2%)", y="PC2 (14.1%)", tag="c") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        plot.title=element_text(size=10),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(5.5,5.5,16,5.5), "pt"))

summary(ord.uwUF.PCoA)

plot.uwUF <- ggplot(scores.uwUF) + 
  geom_point(mapping = aes(x=MDS1, y=MDS2, color=Location, shape=Location)) + 
  geom_segment(data = spp.scrs.uwUF,
               aes(x=0, y=0, xend=3*MDS1, yend=3*MDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="black") +
  geom_path(data=ord.uwUF.ellipse, aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
  scale_color_manual(values=c("forestgreen", "purple")) +
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (16.0%)", y="PC2 (8.8%)", tag="d") +
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
                                     geom_point(mapping = aes(x=MDS1, y=MDS2, color=Location, shape=Location)) + 
                                     geom_segment(data = spp.scrs.BC,
                                                  aes(x=0, y=0, xend=MDS1, yend=MDS2),
                                                  arrow = arrow(length=unit(0.25, "cm")), color="black") +
                                     geom_path(data=ord.BC.ellipse, aes(x=MDS1, y=MDS2, colour=Type), size=0.5, linetype=1) +
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


Fig_S8 <- arrangeGrob(grobs=list(plot.BC, plot.JAC, plot.wUF, plot.uwUF, plot.legend), layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
                      heights=c(0.45,0.45,0.1))
ggsave("~/health/Figure_S8_ordination_pane.jpg", Fig_S8, dpi=300, width=5, height=5, units="in")

rm(ord.BC.PCoA, ord.BC.envfit, ord.BC.ellipse, scores.BC, spp.scrs.BC, plot.BC,
   ord.JAC.PCoA, ord.JAC.envfit, ord.JAC.ellipse, scores.JAC, spp.scrs.JAC, plot.JAC,
   ord.wUF.PCoA, ord.wUF.envfit, ord.wUF.ellipse, scores.wUF, spp.scrs.wUF, plot.wUF,
   ord.uwUF.PCoA, ord.uwUF.envfit, ord.uwUF.ellipse, scores.uwUF, spp.scrs.uwUF, plot.uwUF, plot.legend)

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

Fig_S9 <- arrangeGrob(grobs=list(ord.plot.sex, ord.plot.em), ncol=2)
ggsave("~/health/Tables/Fig_S9_sex_em_ords.jpg", Fig_S9, dpi=300, width=6, height=4, units="in")

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

Fig_S10_anthro.taxa <- arrangeGrob(grobs=list(rgcca.taxon.model.plots[["Streptococcus"]] +
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


Fig_S11_protein.taxa <- arrangeGrob(grobs=list(rgcca.taxon.model.plots[["Fusobacterium"]] +
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

Fig_S12_spleen.taxa <- arrangeGrob(grobs=list(rgcca.taxon.model.plots[["Erysipelotrichaceae"]] + labs(x=NULL, y=NULL) +
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


ggsave("~/health/Tables/Figure_S10_anthro_prediction_models.jpg", Fig_S10_anthro.taxa, dpi=300, height=3.25, width=9, units="in")
ggsave("~/health/Tables/Figure_S11_protein_prediction_models.jpg", Fig_S11_protein.taxa, dpi=300, height=6.5, width=9, units="in")
ggsave("~/health/Tables/Figure_S12_spleen_prediction_models.jpg", Fig_S12_spleen.taxa, dpi=300, height=3.25, width=9, units="in")

#### FIG S13: READ COUNT DISTRIBUTIONS ####
temp <- as.data.frame(sample_sums(subset_samples(pseq.raw, Study=="study")))

# Read count histograms
Fig_S13_reads <- ggplot(temp, aes(x=temp[,1])) + geom_histogram(bins=35) + theme_bw() +
  labs(x="Number of Reads", y="Count", title="Read count histogram", tag="a") +
  scale_x_continuous(breaks=c(0,5000,10000,15000,20000, 25000, 30000, 35000, 40000, 45000, 50000)) +
  scale_y_continuous(limits=c(0,9), expand=c(0,0), breaks=c(0,2,4,6,8)) +
  theme(panel.grid=element_blank(),
        plot.tag=element_text(face="bold"),
        axis.text.y=element_text(color="black", size=9),
        plot.title=element_text(color="black", size=10),
        axis.text.x=element_text(color="black", size=9, angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(color="black", size=10))

# Sample completeness summary
temp <- as.data.frame(t(otu_table(subset_samples(pseq.raw, Study=="study"))))
iNext.result <- iNEXT(temp, q=0, datatype="abundance") # All samples, except the five low-read samples, have completion >98%
rm(temp)

Fig_S13_completeness <- ggplot(iNext.result$DataInfo, aes(x=100*SC)) + geom_histogram() + theme_bw() +
  labs(x="Sample completeness (%)", y="Count", title = "Sample completeness", tag="b") +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  theme(panel.grid=element_blank(),
        plot.tag=element_text(face="bold"),
        axis.text=element_text(color="black", size=9),
        plot.title=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=10))

# Contaminant prevalence distribution
Fig_S13_contaminants <- ggplot(data=subset(contam.data$prevalence, pa.pos > 0), aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  labs(x = "Prevalence (negative controls)", y = "Prevalence (true samples)", title = "Contaminants", tag="c") +
  theme_bw() +
  guides(color=guide_legend(title="Contaminant")) +
  scale_color_manual(values=c("grey", "red")) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(color="black", size=9),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=10),
        legend.position=c(0.95,0.05),
        legend.justification=c("right", "bottom"),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8))

# Contaminant abundances
Fig_S13_abundances <- ggplot(contam.data$abundance, aes(x=ASV, y=means)) + geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  labs(x="ASV", y="Mean relative abundance (%)", title = "Contaminant abundances", tag="d") +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
  theme(panel.grid=element_blank(),
        plot.tag=element_text(face="bold"),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(color="black", size=10),
        plot.title=element_text(color="black", size=10),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.position="bottom")

Fig_S13_reads <- ggplotGrob(Fig_S13_reads)
Fig_S13_completeness <- ggplotGrob(Fig_S13_completeness)
Fig_S13_contaminants <- ggplotGrob(Fig_S13_contaminants)
Fig_S13_abundances <- ggplotGrob(Fig_S13_abundances)

maxHeight = grid::unit.pmax(Fig_S13_reads$heights[1:12], Fig_S13_completeness$heights[1:12])
Fig_S13_reads$heights[1:12] <- as.list(maxHeight)
Fig_S13_completeness$heights[1:12] <- as.list(maxHeight)

maxHeight = grid::unit.pmax(Fig_S13_contaminants$heights[1:12], Fig_S13_abundances$heights[1:12])
Fig_S13_contaminants$heights[1:12] <- as.list(maxHeight)
Fig_S13_abundances$heights[1:12] <- as.list(maxHeight)

rm(maxHeight)

grid.arrange(Fig_S13_reads, Fig_S13_completeness,
             Fig_S13_contaminants, Fig_S13_abundances, ncol=2, nrow=2)

Fig_S13_reads_contaminants <- arrangeGrob(Fig_S13_reads, Fig_S13_completeness,
                                          Fig_S13_contaminants, Fig_S13_abundances, ncol=2, nrow=2)

ggsave("~/health/FigS13_reads_contaminants.jpg", Fig_S13_reads_contaminants, dpi=300, width=6.5, height=6.5, units="in")

#### FIG S14: MOCK COMMUNITY ANALYSIS ####
df <- as.data.frame(
  sample_sums(subset_samples(pseq.raw, sample_names(pseq.raw) %in% c("K1", "K2", "K3", "NCK456"))))
colnames(df) <- "Reads"
df$SampleID <- c("NC1", "NC2", "NC3", "NC4")

plot1 <- ggplot(df, aes(x=SampleID, y=Reads)) + geom_bar(stat="identity", fill="darkgrey") +
  labs(x="Sample ID", y="Number of reads") +
  scale_y_continuous(limits=c(0,260), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text=element_text(size=10, color="black"),
        panel.grid=element_blank())

temp <- tax_glom(pseq.mock.com, taxrank="Genus")
df2 <- as.data.frame(t(otu_table(transform_sample_counts(temp, function(x) 100*x/sum(x)))))
rownames(df2) <- tax_table(temp)[,6]
df2 <- subset(df2, rownames(df2) %in% c("Methylomicrobium", "Sphingomonas", "Escherichia/Shigella",
                                        "Proteus", "Vibrio", "Staphylococcus"))
df2 <- rbind(df2,
             100-sum(df2$NCMC1))
rownames(df2) <- c(rownames(df2)[c(1:6)], "Extra")
df2$Genus <- rownames(df2)
df2$Genus <- factor(df2$Genus, levels=c("Extra", "Methylomicrobium", "Sphingomonas", "Escherichia/Shigella",
                                        "Proteus", "Vibrio", "Staphylococcus"))
df2$Plot <- rep("Plot")

plot2 <- ggplot(df2, aes(x=Plot, y=NCMC1, fill=Genus)) + geom_bar(stat="identity") +
  labs(x=NULL, y="Relative abundance (%)") +
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(color="black"))

plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)
plot2$heights <- plot1$heights

grid.arrange(plot1, plot2, ncol=2, widths=c(1.5,1))
FIG_S14 <- arrangeGrob(plot1, plot2, ncol=2, widths=c(2,1))
ggsave("~/FIG_S14_mock_com.jpg", FIG_S14, width=6.5, height=5, units="in", dpi=300)

rm(temp, df1, df2)

#### FIG S15: COLLINEARITY IN AGE AND HEALTH ####
plot.data.temp <- cy.metadata
plot.data.temp$Sex <- factor(plot.data.temp$Sex, levels=c("M","F"))

Fig_S15 <- ggplot(plot.data.temp, aes(x=Cementum.Age, y=index1.PC1, color=Sex)) + geom_point() + geom_smooth(method="lm", se=FALSE) +
  theme_bw() +
  scale_color_manual(values=c("orange2", "maroon2")) +
  scale_y_continuous(expand=c(0.05,0.05)) +
  scale_x_continuous(expand=c(0.05,0.05)) +
  labs(x="Age", y="Health score") +
  theme(axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        legend.title=element_text(size=10, color="black"),
        panel.grid=element_blank())

ggsave("~/health/Tables/Figure_S15_sex_and_age.jpg", Fig_S14, dpi=300, width=3.5, height=3, units="in")
#### FIG S16: STORAGE EFFECTS ####
# Group samples by their capture month
capture_data <- cy.metadata %>% dplyr::group_by(Capture.Month) %>% dplyr::count()
capture_data$Percent <- 100*capture_data$n/sum(capture_data$n)

# Import average monthly temperatures
temp_data <- data.frame(
  Capture.Month=c("2017/08", "2017/09", "2017/10", "2017/11", "2017/12", "2018/01", "2018/02",
                  "2018/03", "2018/04", "2018/05", "2018/10", "2018/11", "2018/12", "2019/01",
                  "2019/02", "2019/03"),
  Temp=c(17.8, 12.7, 5.6, -7.0, -7.4, -9.3, -12.3, -4.7, 2.0, 15.7,
         4.6, -2.2, -6.4, -6.8, -19.4, -2.7))

capture_data <- merge(capture_data, temp_data, by="Capture.Month", all=TRUE)
capture_data <- melt(capture_data)
capture_data <- subset(capture_data, variable !="n")

scales_y <- list(
  `Percent` = scale_y_continuous(limits=c(0,40), expand=c(0,0)), 
  `Temp` = scale_y_continuous(limits=c(-20,20), expand=c(0,0)) 
)

library(facetscales)

capture_data$variable <- factor(capture_data$variable, levels=c("Temp", "Percent"))

# Seasonal temperatures and sample collection
storage_plot0a <- ggplot(capture_data) + 
  geom_bar(data=subset(capture_data, variable=="Percent"), aes(x=Capture.Month, y=value), stat="identity") +
  geom_line(data=subset(capture_data, variable=="Temp"), aes(x=Capture.Month, y=value), group=1) +
  theme_bw() +
  facet_wrap(~variable, scales="free_y") +
  geom_hline(yintercept=0, color="grey") +
  labs (x="Capture month", y="Percent of sample (bottom); Temperature (*C) (top)",
        title="Sample collection by month", tag="a") +
  facet_grid_sc(rows=vars(variable), scales=list(y = scales_y)) +
  theme(axis.title=element_text(color="black", size=9),
        axis.text.y=element_text(color="black", size=8),
        axis.text.x=element_text(color="black", size=8, angle=90, hjust=1, vjust=0.5),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=10),
        panel.border=element_rect(color="black"))

# b) Compare alpha diversity measurements
capture_data <- cy.metadata

behav_data <- subset(cy.metadata, Location=="urban")
behav_data$Behaviour <- factor(behav_data$Behaviour)
behav_data <- subset(behav_data, is.na(Observed)=="FALSE")

behav_diversity <- subset(cy.metadata, Location=="urban")[,c("FecesID","Behaviour", "Observed_Extrap", "Shannon_Extrap", "PD_Extrap", "NTI_Raw")]
behav_diversity <- reshape2::melt(behav_diversity)

levels(behav_diversity$variable)[levels(behav_diversity$variable)=="Observed_Extrap"] <- "ASV\nRichness"
levels(behav_diversity$variable)[levels(behav_diversity$variable)=="Shannon_Extrap"] <- "Shannon"
levels(behav_diversity$variable)[levels(behav_diversity$variable)=="PD_Extrap"] <- "PD"
levels(behav_diversity$variable)[levels(behav_diversity$variable)=="NTI_Raw"] <- "NTI"

levels(behav_diversity$Behaviour)[levels(behav_diversity$Behaviour)=="conflict"] <- "managed"

storage_plot1 <- ggplot(behav_diversity, aes(x=Behaviour, y=value)) + geom_boxplot() + facet_wrap(~variable, scales="free_y", ncol=4) + theme_bw() +
  labs(x="Behaviour", y="Value", tag="b", title="Alpha-diversity: managed v. roadkill") +
  theme(axis.title=element_text(color="black", size=9),
        axis.text.y=element_text(color="black", size=8),
        axis.text.x=element_text(color="black", size=8, angle=90, hjust=1, vjust=0.5),
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=10),
        panel.border=element_rect(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# c) Compare beta diversity: roadkill vs. managed 
scores <- scores(ord.PCA.AIT, choices=c(1:3), display="sites")
scores <- cbind(scores, cy.metadata.pseq[c("Sex", "Cementum.Age","Em.PCR.Overall", "Location", "FecesID", "Analysis.Season")])
scores <- merge(scores, cy.metadata[,c("FecesID", "Behaviour","Capture.Month")], by="FecesID", all=FALSE)

scores$Season <- rep("winter")
scores$Season[scores$Capture.Month %in% c("2017/08", "2017/09", "2017/10", "2018/05", "2018/10")] <- "fall/spring"

levels(scores$Behaviour)[levels(scores$Behaviour)=="conflict"] <- "managed"

storage_plot2 <- ggplot(subset(scores, Behaviour !="trapped")) + 
  geom_point(mapping = aes(x=PC1, y=PC2, color=Behaviour, shape=Behaviour), size=3) + 
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (11.0%)", y="PC2 (10.0%)", tag="c", title="Aitchison distance-based PCA") +
  scale_x_continuous(limits=c(-5.1,4.5), breaks=c(-5,-2.5,0,2.5)) +
  scale_y_continuous(limits=c(-6.5,2.75), breaks=c(-5,-2.5,0,2.5)) +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        legend.position="bottom",
        legend.title=element_blank(),
        plot.title=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# Compare shoulder season data to winter data for roadkill samples
levels(scores$Season)[levels(scores$Season)=="shoulder"] <- "fall/spring"

storage_plot2_alt <- ggplot(subset(scores, Behaviour=="roadkill")) + 
  geom_point(mapping = aes(x=PC1, y=PC2, color=Season, shape=Season), size=3) + 
  scale_size_continuous(range=c(1,3)) +
  theme_bw() +
  labs(x="PC1 (11.0%)", y="PC2 (10.0%)", tag="d", title="Aitchison distance-based PCA") +
  scale_x_continuous(limits=c(-5.1,4.5), breaks=c(-5,-2.5,0,2.5)) +
  scale_y_continuous(limits=c(-6.5,2.75), breaks=c(-5,-2.5,0,2.5)) +
  theme(axis.title=element_text(color="black", size=9),
        axis.text=element_text(color="black", size=8),
        plot.tag=element_text(face="bold"),
        panel.border=element_rect(color="black"),
        legend.position="bottom",
        legend.title=element_blank(),
        plot.title=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

grid.arrange(storage_plot0a, storage_plot1, storage_plot2, storage_plot2_alt,
             layout_matrix=rbind(c(1,1,1,2,2,2,2,2),c(3,3,3,3,4,4,4,4)))

FigS16_STORAGE <- arrangeGrob(storage_plot0a, storage_plot1, storage_plot2, storage_plot2_alt,
                           layout_matrix=rbind(c(1,1,1,2,2,2,2,2),c(3,3,3,3,4,4,4,4)), heights=c(0.6,0.4))

ggsave("~/FigS16_STORAGE.jpg", Fig_STORAGE, height=6.5, width=6.5, dpi=300, units="in")