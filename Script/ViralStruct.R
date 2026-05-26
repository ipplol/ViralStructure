library(tidyverse)
library(reshape2)
library(ggsignif)
library(tibble)
library(ggpubr)
library(rgl)
library(venneuler)
library(eulerr)
library(ggalluvial)
library(ggrepel)
library(ggfortify)
library(ggmap)
library(ggtern)
library(factoextra)
library(FactoMineR)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(segmented)
library(Hmisc)
library(pROC)
library(RColorBrewer)
library(MASS)
library(Seurat)
library(pls) 
library(SeuratData)
library(stats)
library(glmnet)
library(randomForest)


##【基于结构与病毒流行历史的直接预测】##
MergeFreqIndices <- read.delim("G:/ViralStruct/SARS2/MergeFreqIndices.txt")
RBD<-MergeFreqIndices[MergeFreqIndices$Region=="RBD",]
NTD<-MergeFreqIndices[MergeFreqIndices$Region=="NTD",]
Other<-MergeFreqIndices[MergeFreqIndices$Region=="Other",]

#单个突变
SARS2_Single_RBD <- read.delim("G:/ViralStruct/TradeOff/SARS2_Single_RBD_Cao.txt")
SARS2_Single_Spike <- read.delim("G:/ViralStruct/TradeOff/SARS2_Single_Spike_Bloom.txt")
SARS2_Single_RBD_TylerBloom <- read.delim("G:/ViralStruct/TradeOff/SARS2_Single_RBD_TylerBloom.txt")

#ESC增大ACE2降低
EscAbove0 <- SARS2_Single_Spike[SARS2_Single_Spike$SeraEscape>0,]
ggplot(EscAbove0,aes(SeraEscape,ACE2_binding))+geom_point(aes(color=Clade,shape=region),alpha=0.8)+geom_smooth(method="lm")+stat_cor()+facet_grid(,vars(Clade))
ggplot(EscAbove0,aes(SeraEscape,CellEntry))+geom_point(aes(color=Clade,shape=region),alpha=0.8)+geom_smooth(method="lm")+stat_cor()+facet_grid(,vars(Clade))
tradeofffenbu <- read.delim("G:/ViralStruct/TradeOff/tradeofffenbu.txt")
ggplot(tradeofffenbu,aes(Clade,Proportion,fill=Group))+geom_bar(stat = "identity")+facet_grid(,vars(ESC))

#整个Spike
EscAbove0 <- SARS2_Single_Spike[SARS2_Single_Spike$SeraEscape>0,]
ggplot(EscAbove0,aes(SeraEscape,ACE2_binding))+geom_point(aes(color=Clade))+geom_smooth(method="lm")+stat_cor()+facet_grid(vars(Clade),vars(region))
ggplot(EscAbove0,aes(ACE2_binding,CellEntry))+geom_point(aes(color=Clade))+geom_smooth(method="lm")+stat_cor()+facet_grid(vars(Clade),vars(region))

SpikeSeraESC_ACE2Correlation <- read.delim("G:/ViralStruct/TradeOff/SpikeSeraESC_ACE2Correlation.txt")
ggplot(SpikeSeraESC_ACE2Correlation,aes(Region,Correlation))+geom_bar(stat = "identity")+facet_grid(,vars(Clade))

#解析trade off 阈值
SARS2_RBD_Bloom_ESCMut <- SARS2_Single_Spike[SARS2_Single_Spike$SeraEscape>0,]
SARS2_RBD_Bloom_ESCMut <- SARS2_RBD_Bloom_ESCMut[SARS2_RBD_Bloom_ESCMut$region=="RBD",]
SARS2_RBD_tylerBloom_ESCMut <- SARS2_Single_RBD_TylerBloom[SARS2_Single_RBD_TylerBloom$SeraEscape>0,]
SARS2_RBD_Cao_ESCMut <- SARS2_Single_RBD[SARS2_Single_RBD$ESC>0,]
P1<-ggplot(SARS2_RBD_Cao_ESCMut,aes(delta_bind,ESC,group=BindingGroup,color=BindingGroup))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_x_continuous(limits = c(-4.5,2))
#P2<-ggplot(SARS2_RBD_tylerBloom_ESCMut,aes(deltaYeastBind,SeraEscape))+geom_point()+geom_smooth()+stat_cor()+scale_x_continuous(limits = c(-5,2))
P3<-ggplot(SARS2_RBD_Bloom_ESCMut,aes(ACE2_binding,SeraEscape,group=BindingGroup,color=BindingGroup))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_x_continuous(limits = c(-4.5,2))
cowplot::plot_grid(P1,P3,ncol = 1)
P1<-ggplot(SARS2_RBD_Cao_ESCMut,aes(delta_bind,ESC))+geom_point()+geom_smooth(method = "gam")+stat_cor()+scale_x_continuous(limits = c(-4.5,2))
P3<-ggplot(SARS2_RBD_Bloom_ESCMut,aes(ACE2_binding,SeraEscape))+geom_point()+geom_smooth(method = "gam")+stat_cor()+scale_x_continuous(limits = c(-4.5,2))
cowplot::plot_grid(P1,P3,ncol = 1)

tmp1<-MergeFreqIndices
tmp1$ACE2Group<-MergeFreqIndices$Spike_ACE2binding_Max<0.167
P1<-ggplot(tmp1,aes(Spike_ACE2binding_Max,EscapeMax,shape=Region))+geom_point()+geom_smooth(method = "gam")+stat_cor()+scale_x_continuous(limits = c(0,1))
P3<-ggplot(tmp1,aes(Spike_ACE2binding_Max,SeraEscapeMax,shape=Region))+geom_point()+geom_smooth(method = "gam")+stat_cor()+scale_x_continuous(limits = c(0,1))+geom_vline(xintercept = c(0.167))
cowplot::plot_grid(P1,P3,ncol = 1)

ggplot(SARS2_RBD_Bloom_ESCMut,aes(ACE2_binding,SeraEscape))+geom_point()+geom_smooth(method = "gam")+
  stat_cor()+scale_x_continuous(limits = c(-4.5,2))+facet_grid(vars(Clade),)


ggplot(SARS2_RBD_tylerBloom_ESCMut,aes(deltaPseudoBind,deltaYeastBind))+geom_point()+stat_cor()+geom_smooth()+geom_hline(aes(yintercept = c(-0.625)))+geom_vline(aes(xintercept = c(-0.625)))

#解析适应性阈值断点
x<-MergeFreqIndices$Spike_ACE2binding_Max
y<-MergeFreqIndices$SeraEscapeMax
lin_mod <- lm(y ~ x)
seg_mod_auto <- segmented(lin_mod, seg.Z = ~x, psi = list(x = c(0.1,0.25)))
summary(seg_mod_auto)

#结构与trade off的关联
df_sorted <- MergeFreqIndices[order(MergeFreqIndices$SeraEscapeMax, decreasing = FALSE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=SeraEscapeMax,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#E50049"))
df_sorted <- MergeFreqIndices[order(MergeFreqIndices$Spike_ACE2binding_Max, decreasing = FALSE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=Spike_ACE2binding_Max,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#EBEBEB", "#007CFF", "#F74E69"))
df_sorted <- MergeFreqIndices[order(MergeFreqIndices$WCN_8hri, decreasing = TRUE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=WCN_8hri,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#00FF8A", "#00AFBB", "#EBEBEB"))
df_sorted <- MergeFreqIndices[order(MergeFreqIndices$MutRate, decreasing = FALSE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=MutRate,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))

SARS2 <- read.delim("G:/ViralStruct/StructurePredictESC/SARS2.txt")
df_sorted <- SARS2[order(SARS2$ClashMerge, decreasing = FALSE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=log10(SASA_Sm*log10(ClashA+2*ClashC+1)*Dist_to_Receptor),shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))

RBD<-SARS2[SARS2$Region=="RBD",]
ggplot(RBD,aes(log10(SASA_Sm*log10(ClashA+2*ClashC+1)*Dist_to_Receptor),SeraEscapeMax))+geom_point()+geom_smooth(method="lm")+stat_cor()
#解析随着酵母展示ACE2的下降，有多少突变在假病毒实验中消失
#WithYeastNotPseudo <- read.delim("G:/ViralStruct/TradeOff/data/WithYeastNotPseudo.txt")
#ggplot(WithYeastNotPseudo,aes(deltaACE2Yeast,Proportion,color=Type))+geom_point()+geom_line()

#适应性平衡对病毒进化的影响
#筛选非ACE2驱动的突变
NoACE2<-SARS2_Single_RBD[SARS2_Single_RBD$delta_bind<0,]
NoACE2_ESC<-NoACE2[NoACE2$ESC>0.1,]
NoACE2_ESC_Happened<-NoACE2_ESC[NoACE2_ESC$MutRate>0,]
ggplot(NoACE2_ESC_Happened,aes(Dist_to_Receptor,MutRate))+geom_point()+scale_x_continuous(limits = c(0,35))+geom_smooth(method = "gam")
P1<-ggplot(NoACE2_ESC_Happened,aes(Dist_to_Receptor,weight=MutRate))+geom_density()+scale_x_continuous(limits = c(0,45))
P2<-ggplot(RBD,aes(Dist_to_Receptor,EscapeMax))+geom_point(aes(size=Spike_ACE2binding_Max,color=Spike_ACE2binding_Max))+scale_x_continuous(limits = c(0,45))+scale_color_gradientn(colors = c("#EBEBEB", "#007CFF", "#F74E69"))+geom_smooth(method="gam")
P3<-ggplot(NoACE2_ESC_Happened,aes(Dist_to_Receptor,weight=1))+geom_density()+scale_x_continuous(limits = c(0,45))
cowplot::plot_grid(P1,P2,P3,ncol = 1)
#clade
ggplot(NoACE2_ESC_Happened,aes(Dist_to_Receptor,weight=MutRate, color=Clade))+geom_density()+scale_x_continuous(limits = c(0,35))

#15-20唉小突起是什么
ggplot(RBD,aes(DistGroup2,WCN_8hri))+geom_boxplot()
#不同距离组对适应性的影响大小
ggplot(RBD,aes(DistGroup2,CellEntryMax))+geom_boxplot()
compaired <- list(c("D0-5", "D5-10"),c("D5-10","D10-15"),c("D10-15","D15-20"),c("D15-20","D20plus"))
ggplot(RBD,aes(DistGroup2,Spike_ACE2binding_Max))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
#免疫背景
compaired <- list(c("D0-5", "D5-10"),c("D5-10","D10-15"))
P1<-ggplot(RBD,aes(DistGroup2,ESC_WT))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P2<-ggplot(RBD,aes(DistGroup2,ESC_BA1))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P3<-ggplot(RBD,aes(DistGroup2,ESC_BA2))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P4<-ggplot(RBD,aes(DistGroup2,ESC_BA5))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P5<-ggplot(RBD,aes(DistGroup2,ESC_XBB))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
cowplot::plot_grid(P1,P2,P3,P4,P5,ncol = 5)

#不同dist区域的进化速度差异
LineagePicked <- read.delim("G:/ViralStruct/SARS2/RCoV19/LineagePicked.txt")
P1<-ggplot(LineagePicked,aes(Date,D0.5))+geom_point(aes(color=Clade))+geom_smooth(method = "lm")+scale_y_continuous(limits = c(0,12))
P2<-ggplot(LineagePicked,aes(Date,D5.10))+geom_point(aes(color=Clade))+geom_smooth(method = "lm")+scale_y_continuous(limits = c(0,12))
P3<-ggplot(LineagePicked,aes(Date,D10.15))+geom_point(aes(color=Clade))+geom_smooth(method = "lm")+scale_y_continuous(limits = c(0,12))
P4<-ggplot(LineagePicked,aes(Date,D15.20))+geom_point(aes(color=Clade))+geom_smooth(method = "lm")+scale_y_continuous(limits = c(0,12))
P5<-ggplot(LineagePicked,aes(Date,D20Plus))+geom_point(aes(color=Clade))+geom_smooth(method = "lm")
cowplot::plot_grid(P1,P2,P3,P4,P5,ncol = 3)
lm(D0.5~Date,data=LineagePicked)
lm(D5.10~Date,data=LineagePicked)
lm(D10.15~Date,data=LineagePicked)
lm(D15.20~Date,data=LineagePicked)
lm(D20Plus~Date,data=LineagePicked)

#比较连续进化lineage与跃变，在不同距离区域的突变集中度
LineageMutevent <- read.delim("G:/ViralStruct/SARS2/UShERCoV/LineageMutevent.txt")
ggplot(LineageMutevent,aes(Evolution,Proportion,fill=DistGroup))+geom_bar(stat = "identity")

NoNA<-SARS2_Single_RBD_Expand[SARS2_Single_RBD_Expand$YeastACE2Bind!="NA",]
NoNA<-rbind(NoNA,SARS2_Single_RBD[SARS2_Single_RBD$Happened=="N",])
NoNA<-NoNA[NoNA$ESC>0,]
ggplot(NoNA, aes(x = Clade,y=YeastACE2Bind,color=Happened)) +geom_boxplot()+stat_compare_means(method = "t.test")+geom_hline(yintercept = c(8.188,8.53,8.7,10.09,9.89))

SARS2_Single_RBD_Expand <- read.delim("G:/ViralStruct/TradeOff/data/SARS2_Single_RBD_Cao.txt.Expanded.txt")
tmp<-SARS2_Single_RBD_Expand[SARS2_Single_RBD_Expand$delta_bind<0,]
ggplot(SARS2_Single_RBD, aes(x = YeastACE2Bind,y=ESC)) +stat_cor()+
  geom_point(aes(size = MutRate))+geom_smooth(method = "gam")+facet_grid(vars(Clade),)
ggplot(SARS2_Single_RBD_Expand, aes(x = delta_bind,weight=ESC))+geom_density()+facet_grid(vars(Clade),)
ggplot(SARS2_Single_RBD, aes(x = delta_bind,weight=ESC))+geom_density()+facet_grid(vars(Clade),)
ggplot(SARS2_Single_RBD_Expand, aes(x = delta_bind,y=ESC)) +stat_cor()+
  geom_point(aes(color = Happened))+geom_smooth(method = "gam")
#解析适应性阈值断点
tmp<-SARS2_Single_RBD_Expand[SARS2_Single_RBD_Expand$Clade == "BA.2",]
tmp<-tmp[tmp$MutRate>0,]
x<-tmp$YeastACE2Bind
y<-tmp$ESC
lin_mod <- lm(y ~ x)
seg_mod_auto <- segmented(lin_mod, seg.Z = ~x, psi = list(x = c(8,10)))
summary(seg_mod_auto)

tmp$ACE2Group<-tmp$YeastACE2Bind<8.7
PBA2<-ggplot(tmp, aes(x = YeastACE2Bind,y=ESC,group=ACE2Group)) +
  geom_point()+geom_smooth(method = "lm")+stat_cor()

cowplot::plot_grid(PWT,PBA1,PBA2,PBA5,XBB15,ncol = 1)
#分段绘图
tmp<-SARS2_Single_RBD_Expand[SARS2_Single_RBD_Expand$delta_bind < -0.16,]
ggplot(tmp, aes(x = delta_bind,y=ESC)) +stat_cor()+
  geom_point(aes(color = Happened))+geom_smooth(method = "lm")

#解析DMS数据差异
tmp<-MergeFreqIndices[MergeFreqIndices$ACE2_bind_WT>0.4,]
P1<-ggplot(tmp,aes(SeraEscape_WT_yeastdisplay,WTEscape_WithNeu,color=ACE2_bind_WT))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P2<-ggplot(tmp,aes(SeraEscape_WT_yeastdisplay,WTEscape_NoNeu,color=ACE2_bind_WT))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P3<-ggplot(tmp,aes(SeraEscape_XBB,ESC_XBB,color=ACE2_bind_WT))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P4<-ggplot(tmp,aes(SeraEscapeMax,EscapeMax,color=ACE2_bind_WT))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P5<-ggplot(tmp,aes(SeraEscape_WT_yeastdisplay,SeraEscapeMax,color=ACE2_bind_WT))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
cowplot::plot_grid(P1,P2,P3,P4,ncol = 2)

P1<-ggplot(MergeFreqIndices,aes(SeraEscape_WT_yeastdisplay,WTEscape_WithNeu,color=CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P2<-ggplot(MergeFreqIndices,aes(SeraEscape_WT_yeastdisplay,WTEscape_NoNeu,color=CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P3<-ggplot(MergeFreqIndices,aes(SeraEscape_XBB,ESC_XBB,color=CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P4<-ggplot(MergeFreqIndices,aes(SeraEscapeMax,EscapeMax,color=CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P5<-ggplot(MergeFreqIndices,aes(SeraEscape_WT_yeastdisplay,SeraEscapeMax,color=CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
cowplot::plot_grid(P1,P2,P3,P4,P5,ncol = 3)

#最可以反映抗原性的指标
P1<-ggplot(RBD,aes(MutRate,EscapeMax))+geom_point()+geom_smooth(method = "lm")+ stat_cor(method = "spearman")
P2<-ggplot(RBD,aes(SiteOscillating,EscapeMax))+geom_point()+geom_smooth(method = "lm")+ stat_cor(method = "spearman")
P3<-ggplot(RBD,aes(FitnessEffect,EscapeMax))+geom_point()+geom_smooth(method = "lm")+ stat_cor(method = "spearman")
cowplot::plot_grid(P1,P2,P3,ncol = 3)
fm1<-aov(MutRate~ACE2_bind_WT+EscapeMax+CellEntryMax,data=RBD)
summary(fm1)
fm2<-aov(SiteOscillating~ACE2_bind_WT+EscapeMax+CellEntryMax,data=RBD)
summary(fm2)
fm3<-aov(FitnessEffect~ACE2_bind_WT+EscapeMax+CellEntryMax,data=RBD)
summary(fm3)#手动汇总计算方差解释度
MutIndexAOV <- read.delim("G:/ViralStruct/SARS2/Plot/MutIndexAOV.txt")
ggplot(MutIndexAOV,aes(MutationIndex,VarianceProp,fill=Function))+geom_bar(stat="identity")

#反复突变的结构分布
ggplot(MergeFreqIndices,aes(-1*MDS2,-1*MDS1,color=WCN_8hri,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
ggplot(MergeFreqIndices,aes(-1*MDS2,-1*MDS1,color=Spike_ACE2binding_Max,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#EBEBEB", "#AF8DC8", "#F74E69"))
ggplot(MergeFreqIndices,aes(-1*MDS2,-1*MDS1,color=SeraEscapeMax,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#90D42D"))

#相同WCN情况下，抗体偏向针对RBD？
ggplot(MergeFreqIndices,aes(WCN_8hri,Spike_ACE2binding_Max,color=SeraEscapeMax,shape=Region,size=SeraEscapeMax))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
cor.test(MergeFreqIndices$WCN_8hri,MergeFreqIndices$SeraEscapeMax)
cor.test(MergeFreqIndices$Spike_ACE2binding_Max,MergeFreqIndices$SeraEscapeMax)
compaired <- list(c("NTD", "RBD"))
ggplot(MergeFreqIndices,aes(Region,WCN_8hri,color=Region))+geom_boxplot() + geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = t.test)
ggplot(MergeFreqIndices,aes(Region,SeraEscapeMax,color=Region))+geom_boxplot() + geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = t.test)

#RBD区域的抗体结合与抗体中和
df_sorted <- RBD[order(RBD$WTEscape_NoNeu, decreasing = FALSE), ]
p2<-ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=WTEscape_WithNeu,size=Spike_ACE2binding_Max))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#E50049"))
p3<-ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=WTEscape_WithNeu,size=Spike_ACE2binding_Max))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#E50049"))
cowplot::plot_grid(p2,p3,ncol = 2)
df_sorted <- RBD[order(RBD$WTEscape_WithNeu, decreasing = FALSE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=WTEscape_WithNeu/WTEscape_NoNeu,size=Spike_ACE2binding_Max))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#00AFBB", "#EBEBEB", "#FC4E07"))

#突变率的空间分布
df_sorted <- RBD[order(RBD$SiteOscillating, decreasing = FALSE), ]
ggplot(df_sorted,aes(-1*MDS2,-1*MDS1,color=SiteOscillating,size=Spike_ACE2binding_Max))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))


#对Sera逃逸的解释度
fm1<-aov(SeraEscapeMax~WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=MergeFreqIndices)
summary(fm1)
fm2<-aov(SeraEscapeMax~WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=RBD)
summary(fm2)
fm3<-aov(SeraEscapeMax~WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=NTD)
summary(fm3)
fm4<-aov(SeraEscapeMax~WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=Other)
summary(fm4)
ANOVA <- read.delim("G:/ViralStruct/SARS2/DMS/ANOVA.txt")
p1<-ggplot(ANOVA,aes(1,VariPropSpike,fill=Index))+geom_bar(stat = "identity")
p2<-ggplot(ANOVA,aes(1,VariPropRBD,fill=Index))+geom_bar(stat = "identity")
p3<-ggplot(ANOVA,aes(1,VariPropNTD,fill=Index))+geom_bar(stat = "identity")
p4<-ggplot(ANOVA,aes(1,VariPropOther,fill=Index))+geom_bar(stat = "identity")
cowplot::plot_grid(p1,p2,p3,p4,ncol = 4)

fm2<-aov(EscapeMax~WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=RBD)
summary(fm2)

#与免疫逃逸的相关性
tmpSpike<-MergeFreqIndices[,c("SeraEscapeMax","WCN_8hri","Dist_to_Receptor","Dist_to_Glycan","B_factor_Z","Spike_ACE2binding_Max","CellEntryMax")]
tmpRBD<-RBD[,c("SeraEscapeMax","WCN_8hri","Dist_to_Receptor","Dist_to_Glycan","B_factor_Z","Spike_ACE2binding_Max","CellEntryMax")]
tmpNTD<-NTD[,c("SeraEscapeMax","WCN_8hri","Dist_to_Receptor","Dist_to_Glycan","B_factor_Z","Spike_ACE2binding_Max","CellEntryMax")]
tmpOther<-Other[,c("SeraEscapeMax","WCN_8hri","Dist_to_Receptor","Dist_to_Glycan","B_factor_Z","Spike_ACE2binding_Max","CellEntryMax")]
CorMatrix_Spike<-data.frame(cor(tmpSpike, use = "pairwise.complete.obs"))
CorMatrix_Spike$Region<-"Spike"
CorMatrix_RBD<-data.frame(cor(tmpRBD, use = "pairwise.complete.obs"))
CorMatrix_RBD$Region<-"RBD"
CorMatrix_NTD<-data.frame(cor(tmpNTD, use = "pairwise.complete.obs"))
CorMatrix_NTD$Region<-"NTD"
CorMatrix_Other<-data.frame(cor(tmpOther, use = "pairwise.complete.obs"))
CorMatrix_Other$Region<-"Other"
SeraEscapeCor<-rbind(CorMatrix_Spike,CorMatrix_RBD,CorMatrix_NTD,CorMatrix_Other)
SeraEscapeCor$indics<-(row.names(SeraEscapeCor))
ggplot(SeraEscapeCor,aes(indics,SeraEscapeMax,fill=Region))+geom_bar(stat="identity")
ggplot(MergeFreqIndices,aes(Dist_to_Receptor,SeraEscapeMax))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri,shape=Region))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
ggplot(RBD,aes(Dist_to_Receptor,EscapeMax))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri),shape=15)+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))


P1<-ggplot(MergeFreqIndices,aes(Dist_to_Receptor,SeraEscapeMax))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri,shape=Region))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P2<-ggplot(MergeFreqIndices,aes(Dist_to_Receptor,EscapeMax))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri),shape=15)+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P3<-ggplot(MergeFreqIndices,aes(Dist_to_Receptor,SiteOscillating))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri,shape=Region))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P4<-ggplot(MergeFreqIndices,aes(Dist_to_Receptor,Count))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri,shape=Region))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P3,P4,ncol = 1)


#不同突变指标的功能解释度
fm1<-aov(MutRate~EscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=MergeFreqIndices)
summary(fm1)
fm2<-aov(FitnessEffect~EscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=MergeFreqIndices)
summary(fm2)
fm3<-aov(SiteOscillating~EscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=MergeFreqIndices)
summary(fm3)
fm1<-aov(MutRate~EscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+ACE2_bind_WT+CellEntryMax,data=RBD)
summary(fm1)
fm2<-aov(FitnessEffect~EscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+ACE2_bind_WT+CellEntryMax,data=RBD)
summary(fm2)
fm3<-aov(SiteOscillating~EscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+ACE2_bind_WT+CellEntryMax,data=RBD)
summary(fm3)
ANOVA <- read.delim("G:/ViralStruct/SARS2/DMS/ANOVA.txt")
ggplot(ANOVA,aes(1,VariProp,fill=Index))+geom_bar(stat = "identity")+facet_grid(. ~ Type) 

fm1<-aov(MutRate~SeraEscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=NTD)
summary(fm1)
fm3<-aov(SiteOscillating~SeraEscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=NTD)
summary(fm3)

fm1<-aov(MutRate~SeraEscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=Other)
summary(fm1)
fm3<-aov(SiteOscillating~SeraEscapeMax+WCN_8hri+Dist_to_Receptor+Dist_to_Glycan+B_factor_Z+Spike_ACE2binding_Max+CellEntryMax,data=Other)
summary(fm3)

ANOVA_Other <- read.delim("G:/ViralStruct/SARS2/DMS/ANOVA_Other.txt")
ggplot(ANOVA_Other,aes(1,VariProp,fill=Index))+geom_bar(stat = "identity")+facet_grid(. ~ Type) 

#对Incidence的分析
#Sera最优逃逸对比理论最优逃逸
P1<-ggplot(RBD,aes(SeraEscapeMax,Count))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))+stat_cor()
P2<-ggplot(RBD,aes(EscapeMax,Count))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))+stat_cor()
P3<-ggplot(RBD,aes(SeraEscapeMax,SiteOscillating))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))+stat_cor()
P4<-ggplot(RBD,aes(EscapeMax,SiteOscillating))+geom_point(aes(size=Spike_ACE2binding_Max,color=WCN_8hri))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))+stat_cor()
cowplot::plot_grid(P1,P2,P3,P4,ncol = 2)

EscapeSite<-RBD[RBD$EscapeMut=="Yes",]
#EscapeSite<-RBD[RBD$EscapeMax>0.5,]
#免疫逃逸突变的基因组与结构分布
ggplot(RBD,aes(Pos,EscapeMax,color=EscapeMut,label=Pos))+geom_point(size=3)+geom_text()
ggplot(RBD,aes(EscapeMut,Dist_to_Receptor))+geom_boxplot()+stat_compare_means(aes(group = EscapeMut), method = "wilcox.test")
#按incidence分组逃逸突变
RBD.txt.Expanded <- read.delim("G:/ViralStruct/SARS2/RBD.txt.Expanded.txt")
EscapeSite_exp<-RBD.txt.Expanded[RBD.txt.Expanded$EscapeMut=="Yes",]
EscapeSite_exp$IfHighIncidence<-EscapeSite_exp$Incidence=="high"
EscapeSite_exp_high<-EscapeSite_exp[EscapeSite_exp$Incidence=="High",]
EscapeSite_exp_low<-EscapeSite_exp[EscapeSite_exp$Incidence=="Low",]

t.test(EscapeSite_exp_high$Dist_to_Receptor,EscapeSite_exp_low$Dist_to_Receptor)
P1<-ggplot(EscapeSite, aes(x = Dist_to_Receptor, weight = Count,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(7.18,6.11,6.645),col=c("red","blue","black"))

t.test(EscapeSite_exp_high$CellEntryMax,EscapeSite_exp_low$CellEntryMax)
P2<-ggplot(EscapeSite, aes(x = CellEntryMax, weight = Count,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(0.39,0.55,0.47),col=c("red","blue","black"))

t.test(EscapeSite_exp_high$ACE2_bind_WT,EscapeSite_exp_low$ACE2_bind_WT)
P3<-ggplot(EscapeSite, aes(x = ACE2_bind_WT, weight = Count,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(0.44,0.63,0.535),col=c("red","blue","black"))

cowplot::plot_grid(P1,P2,P3,ncol = 3)

#对反复突变的分析
EscapeSite_exp_OSC<-EscapeSite_exp[EscapeSite_exp$Oscillating=="Yes",]
EscapeSite_exp_Other<-EscapeSite_exp[EscapeSite_exp$Oscillating=="No",]
t.test(EscapeSite_exp_OSC$Dist_to_Receptor,EscapeSite_exp_Other$Dist_to_Receptor)
t.test(EscapeSite_exp_OSC$CellEntryMax,EscapeSite_exp_Other$CellEntryMax)
t.test(EscapeSite_exp_OSC$ACE2_bind_WT,EscapeSite_exp_Other$ACE2_bind_WT)

P1<-ggplot(EscapeSite, aes(x = Dist_to_Receptor, weight = SiteOscillating,fill = Oscillating)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(8.39,5.02,6.7),col=c("red","blue","black"))

P2<-ggplot(EscapeSite, aes(x = CellEntryMax, weight = SiteOscillating,fill = Oscillating)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(0.398,0.428,0.413),col=c("red","blue","black"))

P3<-ggplot(EscapeSite, aes(x = ACE2_bind_WT, weight = SiteOscillating,fill = Oscillating)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(0.353,0.645,0.499),col=c("red","blue","black"))

cowplot::plot_grid(P1,P2,P3,ncol = 3)


#RBD区域反复突变位点的功能特征
ggplot(RBD,aes(ACE2_bind_WT,EscapeMax,size=SiteOscillating,color=CellEntryMax))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P1<-ggplot(RBD,aes(Oscillating,EscapeMax,group=Oscillating,color=Oscillating))+geom_violin(draw_quantiles = 0.5)+stat_compare_means(label = "p.signif", method = "wilcox.test")
P2<-ggplot(RBD,aes(Oscillating,ACE2_bind_WT,group=Oscillating,color=Oscillating))+geom_violin(draw_quantiles = 0.5)+stat_compare_means(label = "p.signif", method = "wilcox.test")
P3<-ggplot(RBD,aes(Oscillating,CellEntryMax,group=Oscillating,color=Oscillating))+geom_violin(draw_quantiles = 0.5)+stat_compare_means(label = "p.signif", method = "wilcox.test")
cowplot::plot_grid(P1,P2,P3,ncol = 3)

#逃逸位点的结构分布
RBD$ESCColorOSC <- paste(RBD$EscapeMut, RBD$Oscillating, sep = " ")
RBD$ESCColorIncidence <- paste(RBD$EscapeMut, RBD$Incidence, sep = " ")
P1<-ggplot(RBD,aes(-1*MDS2,-1*MDS1,color=ESCColorOSC,shape=Region,size=ACE2_bind_WT))+geom_point(alpha=0.8)
P2<-ggplot(RBD,aes(-1*MDS2,-1*MDS1,color=ESCColorIncidence,shape=Region,size=ACE2_bind_WT))+geom_point(alpha=0.8)
cowplot::plot_grid(P1,P2,ncol = 2)
ggplot(MergeFreqIndices,aes(-1*MDS2,-1*MDS1,color=EscapeMax,shape=Region))+geom_point(alpha=0.8,size=3)+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))

#反复突变热力图
OscillatingSite.UK.detail <- read.delim("G:/ViralStruct/SARS2/UShERCoV/OscillatingSite.UK.tsv.detail", row.names=1)
OscillatingSite.UK.freq <- read.delim("G:/ViralStruct/SARS2/UShERCoV/OscillatingSite.UK.tsv.freq", row.names=1)
pheatmap(OscillatingSite.UK.freq, cluster_cols = FALSE,cluster_rows = FALSE, display_numbers = OscillatingSite.UK.detail)
OscillatingSite.USA.detail <- read.delim("G:/ViralStruct/SARS2/UShERCoV/OscillatingSite.USA.tsv.detail", row.names=1)
OscillatingSite.USA.freq <- read.delim("G:/ViralStruct/SARS2/UShERCoV/OscillatingSite.USA.tsv.freq", row.names=1)
pheatmap(OscillatingSite.USA.freq, cluster_cols = FALSE,cluster_rows = FALSE, display_numbers = OscillatingSite.USA.detail)
OSCESC_CellE <- read.delim("G:/ViralStruct/SARS2/OSCMut/OSCESC_CellE.txt", row.names=1)
OSCESC_ACE2 <- read.delim("G:/ViralStruct/SARS2/OSCMut/OSCESC_ACE2.txt", row.names=1)
OSCESC_ESC <- read.delim("G:/ViralStruct/SARS2/OSCMut/OSCESC_ESC.txt", row.names=1)
pheatmap(OSCESC_CellE, cluster_cols = FALSE,cluster_rows = FALSE,color = colorRampPalette(colors = c("white","#996DC2"))(100))
pheatmap(OSCESC_ACE2, cluster_cols = FALSE,cluster_rows = FALSE,color = colorRampPalette(colors = c("white","#F05F7C"))(100))
pheatmap(OSCESC_ESC, cluster_cols = FALSE,cluster_rows = FALSE,color = colorRampPalette(colors = c("white","#170F5C"))(100))

#RBD反复突变集中
a <- matrix(c(17, 18, 184, 1055), ncol = 2)
fisher.test(a)

#反复突变位点的结构分布
ggplot(MergeFreqIndices,aes(-1*MDS2,-1*MDS1,color=WCN_VirusOnly,size=SiteOscillating,shape=Region))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))

#反复突变位点的功能特征
P1<-ggplot(MergeFreqIndices,aes(Region,WCN_VirusOnly,color=Oscillating))+geom_boxplot()+stat_compare_means(method = "wilcox.test")
P2<-ggplot(MergeFreqIndices,aes(Region,SeraEscapeMax,color=Oscillating))+geom_boxplot()+stat_compare_means(method = "wilcox.test")
P3<-ggplot(MergeFreqIndices,aes(Region,EscapeMax,color=Oscillating))+geom_boxplot()+stat_compare_means(method = "wilcox.test")
P4<-ggplot(MergeFreqIndices,aes(Region,Spike_ACE2binding_Max,color=Oscillating))+geom_boxplot()+stat_compare_means(method = "wilcox.test")
P5<-ggplot(MergeFreqIndices,aes(Region,CellEntryMax,color=Oscillating))+geom_boxplot()+stat_compare_means(method = "wilcox.test")
cowplot::plot_grid(P1,P2,P3,P4,P5,ncol = 5)

#MDS1与受体距离的相关性
cor.test(MergeFreqIndices$MDS1,MergeFreqIndices$Dist_to_Receptor)
ggplot(MergeFreqIndices,aes(MDS1,Dist_to_Receptor))+geom_point(aes(shape=Region,color=WCN_VirusOnly))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_cor(method = "pearson")+ stat_smooth(method = lm, level = 0.99)

#抗原位点的结构分布
ESCMut<-RBD[RBD$EscapeMut=="Yes",]
NonEscMut<-RBD[RBD$EscapeMut=="No",]
ggplot(RBD,aes(CellEntryMax,Dist_to_Receptor,size=SiteOscillating,color=EscapeMax,shape=Oscillating))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(RBD,aes(ACE2_bind_WT,Dist_to_Receptor,size=Count,color=EscapeMax,shape=Incidence))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(ESCMut,aes(Incidence,Dist_to_Receptor))+geom_violin(draw_quantiles = 0.5)+stat_compare_means(method = "wilcox.test")
ggplot(ESCMut,aes(Incidence,ACE2_bind_WT))+geom_violin(draw_quantiles = 0.5)+stat_compare_means(method = "wilcox.test")
var.test(Dist_to_Receptor ~ KeyAntigenic, data = ESCMut)

#反复突变更有可能是逃逸突变而非感染力突变
RBDIndRank <- read.delim("G:/ViralStruct/SARS2/OSCMut/RBDIndRank.txt")
P1<-ggplot(RBDIndRank,aes(ESCRank,EscapeMax,color=Oscillating))+geom_point()
P2<-ggplot(RBDIndRank,aes(ACE2Rank,ACE2_bind_WT,color=Oscillating))+geom_point()
cowplot::plot_grid(P1,P2,ncol = 1)
tmp<-RBDIndRank[RBDIndRank$Oscillating=="Yes",]
wilcox.test(tmp$ESCRank,tmp$ACE2Rank)

RBDIncidenceRank <- read.delim("G:/ViralStruct/SARS2/OSCMut/RBDIncidenceRank.txt")
P1<-ggplot(RBDIncidenceRank,aes(ESCRank,EscapeMax,color=Incidence))+geom_point()
P2<-ggplot(RBDIncidenceRank,aes(ACE2Rank,ACE2_bind_WT,color=Incidence))+geom_point()
cowplot::plot_grid(P1,P2,ncol = 1)
tmp<-RBDIndRank[RBDIncidenceRank$Incidence=="B_High",]
wilcox.test(tmp$ESCRank,tmp$ACE2Rank)

#到受体的距离与进化速度
ggplot(RBD,aes(Dist_to_Receptor,Count))+geom_point()+geom_smooth()
compaired <- list(c("Contact", "Close"),c("Far", "Close"))
P1<-ggplot(RBD,aes(DistGroup,Count,color=DistGroup))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P2<-ggplot(RBD,aes(DistGroup,SiteOscillating,color=DistGroup))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
cowplot::plot_grid(P1,P2,ncol = 2)
ggplot(RBD,aes(Dist_to_Receptor,Count))+geom_point(aes(color=DistGroup))+ geom_smooth(method="gam")


#486解析
ggplot(RBD,aes(ACE2_bind_WT,Dist_to_Receptor,size=SiteOscillating,color=ESCColorOSC,shape=EscapeMut,label=Pos))+geom_point()
ggplot(RBD,aes(Spike_ACE2binding_Max,Dist_to_Receptor,size=SiteOscillating,color=ESCColorOSC,shape=EscapeMut,label=Pos))+geom_point()+geom_text()

#甜区区域拟合
KeyAntigenic<-RBD[RBD$KeyAntigenic=="Yes",]
hist(KeyAntigenic$Dist_to_Receptor)
shapiro.test(KeyAntigenic$Dist_to_Receptor)
shapiro.test(ESCMut$Dist_to_Receptor)
# 估计正态分布参数
mu <- mean(KeyAntigenic$Dist_to_Receptor)
sigma <- sd(KeyAntigenic$Dist_to_Receptor)
hist(KeyAntigenic$Dist_to_Receptor, prob=TRUE, col="lightblue", 
     main="正态性检验与拟合", xlab="Value", breaks=20)
curve(dnorm(x, mean=mu, sd=sigma), col="red", lwd=2, add=TRUE)
legend("topright", legend=c("理论正态曲线"), col=c("red"), lwd=2)
mu;sigma


#RBD区域分布情况
ggplot(RBD,aes(-1*MDS2,-1*MDS1,size=SiteOscillating,color=CellEntryMax,shape=KeyAntigenic))+geom_point()+scale_color_gradientn(colors = c("#F0F0F0", "#E7B800", "#FF0000"))
ggplot(RBD,aes(-1*MDS2,-1*MDS1,size=SiteOscillating,color=ACE2_bind_WT,shape=KeyAntigenic))+geom_point()+scale_color_gradientn(colors = c("#F0F0F0", "#E7B800", "#FF0000"))


#对比EVE
VSEVEscape <- read.delim("G:/ViralStruct/SARS2/VSEVEscape.txt")
P1<-ggplot(VSEVEscape,aes(KeyAntigenic,-1*fitness_eve,label=Pos))+geom_text()+geom_violin()+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_hline(yintercept = c(0.7,4.66))
P2<-ggplot(VSEVEscape,aes(KeyAntigenic,-1*evescape,label=Pos))+geom_text()+geom_violin()+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_hline(yintercept = c(0.91,1.61))
P3<-ggplot(VSEVEscape,aes(KeyAntigenic,WCN_VirusOnly,label=Pos))+geom_text()+geom_violin()+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_hline(yintercept = c(0.8,1.29))
P4<-ggplot(VSEVEscape,aes(KeyAntigenic,Dist_to_Receptor,label=Pos))+geom_text()+geom_violin()+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_hline(yintercept = c(3.64,18.76))
cowplot::plot_grid(P1,P2,P3,P4,ncol = 4)

#H1N1pdm
HA1.XYZ <- read.delim("G:/ViralStruct/OtherVirus/H1N1pdm/Data/PDB/HA1.XYZ.txt", header=FALSE, row.names=1)
MDS<-cmdscale(dist(HA1.XYZ))
write.csv(MDS,"G:/ViralStruct/OtherVirus/H1N1pdm/Data/PDB/HA1.XYZ.MDS.csv")
HA2.XYZ <- read.delim("G:/ViralStruct/OtherVirus/H1N1pdm/Data/PDB/HA2.XYZ.txt", header=FALSE, row.names=1)
MDS<-cmdscale(dist(HA2.XYZ))
write.csv(MDS,"G:/ViralStruct/OtherVirus/H1N1pdm/Data/PDB/HA2.XYZ.MDS.csv")

MergeH1N1pdm <- read.delim("G:/ViralStruct/OtherVirus/H1N1pdm/MergeH1N1pdm.txt")
MergeH1N1pdm$ESCColorOSC <- paste(MergeH1N1pdm$EscapeMut, MergeH1N1pdm$Oscillating, sep = " ")
MergeH1N1pdm$ESCColorIncidence <- paste(MergeH1N1pdm$EscapeMut, MergeH1N1pdm$Incidence, sep = " ")
HA1<-MergeH1N1pdm[MergeH1N1pdm$Region=="HA1",]
ggplot(MergeH1N1pdm,aes(-1*DMS_ENTROPY,Dist_to_Receptor,size=SiteOscillating,color=DMSEscape,shape=KeyAntigenic))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_hline(yintercept = c(8,44.9))
ggplot(MergeH1N1pdm,aes(-1*DMS_ENTROPY,EscapeMax,size=SiteOscillating,weight=SiteOscillating))+
  geom_point(aes(color=ESCColorOSC,shape=EscapeMut))
#结构降维
ggplot(HA1,aes(MDS2,MDS1,size=-1*DMS_ENTROPY,color=ESCColorOSC,shape=Region))+geom_point()
ggplot(HA1,aes(MDS2,MDS1,size=-1*DMS_ENTROPY,color=ESCColorIncidence,shape=Region))+geom_point()
ggplot(HA1,aes(MDS2,MDS1,size=-1*DMS_ENTROPY,color=Dist_to_Receptor,shape=Region))+geom_point()+scale_color_gradientn(colors = c("#EBEBEB", "#000000"))

compaired <- list(c("Other", "ESCNo"),c("ESCOSC", "ESCNo"))
P1<-ggplot(HA1, aes(TypeOSC, Dist_to_Receptor,color = TypeOSC)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P2<-ggplot(HA1, aes(TypeOSC, -1*DMS_ENTROPY,color = TypeOSC)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
cowplot::plot_grid(P1,P2,ncol = 2)

compaired <- list(c("Other", "ESCLow"),c("ESCHigh", "ESCLow"))
P1<-ggplot(HA1, aes(TypeEvent, Dist_to_Receptor,color = TypeEvent)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P2<-ggplot(HA1, aes(TypeEvent, -1*DMS_ENTROPY,color = TypeEvent)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
cowplot::plot_grid(P1,P2,ncol = 2)

#预测作用
OSC<-MergeH1N1pdm[MergeH1N1pdm$SiteOscillating>1,]
ggplot(OSC,aes(EscapeMut,WCN_VirusOnly))+geom_violin(draw_quantiles = 0.5)+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(OSC,aes(EscapeMut,EVEFitness))+geom_violin(draw_quantiles = 0.5)+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(OSC,aes(EscapeMut,-1*DMS_ENTROPY))+geom_violin(draw_quantiles = 0.5)+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))

#H3N2
#单突变位点分析
SingleMutEffect <- read.delim("G:/ViralStruct/OtherVirus/H3N2/H3N2_singleMut.txt")
ggplot(SingleMutEffect, aes(Occurred, Dist_to_Receptor)) + geom_boxplot()
ggplot(SingleMutEffect, aes(Occurred, MDCKSIAT1.cell.entry)) + geom_boxplot()

H3N2ESCAbove0<-SingleMutEffect[SingleMutEffect$sera.escape>0,]
df_sorted <- H3N2ESCAbove0[order(H3N2ESCAbove0$Dist_to_Receptor, decreasing = TRUE), ]
ggplot(df_sorted,aes(sera.escape,MDCKSIAT1.cell.entry))+geom_point(aes(color=Dist_Group),size=3,alpha=0.7)+geom_smooth(method="lm")+stat_cor()
SARS2_Single_Spike <- read.delim("G:/ViralStruct/TradeOff/SARS2_Single_Spike_Bloom.txt")
EscAbove0 <- SARS2_Single_Spike[SARS2_Single_Spike$SeraEscape>0,]
df_sorted <- EscAbove0[order(EscAbove0$Dist_to_Receptor, decreasing = TRUE), ]
ggplot(df_sorted,aes(SeraEscape,CellEntry))+geom_point(aes(color=Dist_Group),alpha=0.8)+geom_smooth(method="lm")+stat_cor()+facet_grid(,vars(Clade))


HA1.XYZ <- read.delim("G:/ViralStruct/OtherVirus/H3N2/Data/PDB/HA1.XYZ.txt", header=FALSE, row.names=1)
MDS<-cmdscale(dist(HA1.XYZ))
write.csv(MDS,"G:/ViralStruct/OtherVirus/H3N2/Data/PDB/HA1.XYZ.MDS.csv")
MergeH3N2 <- read.delim("G:/ViralStruct/OtherVirus/H3N2/MergeH3N2.txt")
MergeH3N2$TypeOSC <- paste(MergeH3N2$EscapeMut, MergeH3N2$Oscillating, sep = " ")
MergeH3N2$TypeEvent <- paste(MergeH3N2$EscapeMut, MergeH3N2$Incidence, sep = " ")
HA1<-MergeH3N2[MergeH3N2$Region=="HA1",]
HA1top<-HA1[HA1$Dist_to_Receptor<5,]
HA1HighRate<-HA1[HA1$Oscillating=="Yes",]
#DMS对比
HA1_sorted <- HA1[order(HA1$WCN_VirusOnly, decreasing = TRUE), ]
ggplot(HA1_sorted,aes(Dist_to_Receptor,SeraEscape))+geom_point(aes(size=-1*DMS_ENTROPY,color=WCN_VirusOnly))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
ggplot(HA1_sorted,aes(Dist_to_Receptor,SeraEscape))+geom_point(aes(size=CellEntry,color=WCN_VirusOnly))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
ggplot(HA1_sorted,aes(-1*DMS_ENTROPY,CellEntry))+geom_point()

ggplot(HA1,aes(Dist_Group,SeraEscape))+geom_boxplot()
ggplot(HA1,aes(Dist_Group,CellEntry))+geom_boxplot()
ggplot(HA1,aes(Dist_to_Receptor,SeraEscape))+geom_point(aes(size=-1*DMS_ENTROPY,color=WCN_VirusOnly))+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))+scale_x_continuous(limits = c(0,10))
ggplot(HA1HighRate,aes(Dist_to_Receptor,SeraEscape))+geom_point(aes(size=-1*DMS_ENTROPY,color=WCN_VirusOnly))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))


#突变发生率
compaired <- list(c("D0", "D1"),c("D1", "D2"),c("D2", "D3"),c("D3", "D4"))
HA1SeraESCMut<-HA1[HA1$SeraEscape>0.1,]
ggplot(HA1SeraESCMut, aes(Dist_Group, MutEvent)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
ggplot(HA1SeraESCMut,aes(Dist_Group,MutEvent))+geom_boxplot()
ggplot(HA1,aes(Dist_to_Receptor,weight=MutEvent))+geom_density()
ggplot(HA1SeraESCMut,aes(Dist_to_Receptor,weight=MutEvent))+geom_density()
ggplot(HA1,aes(Dist_Group,MutEvent))+geom_boxplot()
#+scale_y_continuous(limits = c(0,45))
cowplot::plot_grid(P1,P2,P3,ncol = 1)
#单个突变
compaired <- list(c("D0", "D1"),c("D1", "D2"),c("D2", "D3"),c("D3", "D4"))
H3N2_singleMut <- read.delim("G:/ViralStruct/OtherVirus/H3N2/H3N2_singleMut.txt")
ggplot(H3N2_singleMut,aes(Dist_Group,sera.escape))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
H3N2ESCSingle <- H3N2_singleMut[H3N2_singleMut$sera.escape>0.1,]
H3N2ESCSingle_NoCellentry <- H3N2ESCSingle[H3N2ESCSingle$MDCKSIAT1.cell.entry<0,]
ggplot(H3N2ESCSingle_NoCellentry,aes(Dist_Group,Incidence))+geom_boxplot(outliers = FALSE)

#进化速度
WindowMut <- read.delim("G:/ViralStruct/OtherVirus/H3N2/EvoSpeed/WindowMut.txt")
P1<-ggplot(WindowMut,aes(Window,Contact))+geom_point()+geom_smooth(method = "lm")
P2<-ggplot(WindowMut,aes(Window,Close10))+geom_point()+geom_smooth(method = "lm")
P3<-ggplot(WindowMut,aes(Window,Close15))+geom_point()+geom_smooth(method = "lm")
P4<-ggplot(WindowMut,aes(Window,Close20))+geom_point()+geom_smooth(method = "lm")
P5<-ggplot(WindowMut,aes(Window,Far))+geom_point()+geom_smooth(method = "lm")
cowplot::plot_grid(P1,P2,P3,P4,P5,ncol = 3)
lm(Contact~Window,data=WindowMut)
lm(Close10~Window,data=WindowMut)
lm(Close15~Window,data=WindowMut)
lm(Close20~Window,data=WindowMut)
lm(Far~Window,data=WindowMut)

ggplot(HA1,aes(Dist_to_Receptor,CellEntry))+geom_point(aes(size=-1*DMS_ENTROPY,color=WCN_VirusOnly))+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))

#对比新冠的距离分组逃逸
ggplot(HA1,aes(Dist_Group,SeraEscape))+geom_boxplot()
Dist <- read.delim("G:/ViralStruct/OtherVirus/H3N2/Dist.txt")
ggplot(Dist,aes(Region,Dist_to_Receptor,color=Incidence))+geom_boxplot()+stat_compare_means(aes(group = Region), method = "wilcox.test")
ggplot(DistExpand,aes(Region,Dist_to_Receptor))+geom_boxplot(outliers = FALSE)

DistExpand <- read.delim("G:/ViralStruct/OtherVirus/H3N2/Dist.txt.Expanded.txt")
ggplot(DistExpand,aes(Dist_to_Receptor,weight=1))+geom_boxplot()+stat_compare_means(aes(group = Region), method = "wilcox.test")

DistGroup <- read.delim("G:/ViralStruct/OtherVirus/H3N2/DistGroup.txt")
Escape<-DistGroup[DistGroup$EscapeMut=="Yes",]
Surface<-DistGroup[DistGroup$RSAType=="Outside",]
GContact<-Surface[Surface$Dist_to_Receptor<15,]
ggplot(GContact, aes(RankP, Dist_to_Receptor,group = Virus))+geom_line() + geom_point()#aes(size=Incidence,color=EscapeMut))
ggplot(Escape, aes(Virus, Dist_to_Receptor,group = Virus))+geom_boxplot()

SurfaceRegionProportion <- read.delim("G:/ViralStruct/OtherVirus/H3N2/SurfaceRegionProportion.txt")
ggplot(SurfaceRegionProportion,aes(Region,Proportion,fill=Type))+geom_bar(stat="identity")+facet_grid(vars(Virus),)
ggplot(SurfaceRegionProportion,aes(Region,log10(ComparetoBase),fill=Type))+geom_bar(stat="identity")+facet_grid(vars(Virus),)

#分布密度图的比较
ggplot(HA1, aes(x = Dist_to_Receptor, weight = MutEvent,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(15.17,15.43),col=c("red","blue"))
ggplot(Spi, aes(x = Dist_to_Receptor, weight = SiteOscillating,fill = Oscillating)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(15.17,15.43),col=c("red","blue"))
#分布密度图
ESCMut<-HA1[HA1$EscapeMut=="Yes",]
ESCMut_high<-ESCMut[ESCMut$Incidence=="Yes",]
ESCMut_low<-ESCMut[ESCMut$Incidence=="No",]
t.test(ESCMut_high$Dist_to_Receptor,ESCMut_low$Dist_to_Receptor)
P1<-ggplot(ESCMut, aes(x = Dist_to_Receptor, weight = MutEvent,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(15.17,15.43),col=c("red","blue"))

t.test(ESCMut_high$CellEntry,ESCMut_low$CellEntry)
P2<-ggplot(ESCMut, aes(x = CellEntry, weight = MutEvent,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(1.02,1.99),col=c("red","blue"))

t.test(ESCMut_high$pHstability,ESCMut_low$pHstability)
P3<-ggplot(ESCMut, aes(x = pHstability, weight = MutEvent,fill = Incidence)) + 
  geom_density(alpha = 0.6)+geom_vline(xintercept = c(0.0434,0.0904),col=c("red","blue"))
cowplot::plot_grid(P1,P2,P3,ncol = 3)
#功能
compaired <- list(c("No No", "No Yes"),c("No No", "Yes No"),c("No Yes", "Yes Yes"),c("Yes Yes", "Yes No"))
P1<-ggplot(HA1, aes(TypeOSC, Dist_to_Receptor,color = TypeOSC)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P2<-ggplot(HA1, aes(TypeOSC, pHstability,color = TypeOSC)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P3<-ggplot(HA1, aes(TypeOSC, CellEntry,color = TypeOSC)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
cowplot::plot_grid(P1,P2,P3,ncol = 3)

compaired <- list(c("No No", "No Yes"),c("No No", "Yes No"),c("No Yes", "Yes Yes"),c("Yes Yes", "Yes No"))
P1<-ggplot(HA1, aes(TypeEvent, Dist_to_Receptor,color = TypeEvent)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P2<-ggplot(HA1, aes(TypeEvent, pHstability,color = TypeEvent)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
P3<-ggplot(HA1, aes(TypeEvent, CellEntry,color = TypeEvent)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
cowplot::plot_grid(P1,P2,P3,ncol = 3)
#结构
ggplot(HA1,aes(-1*MDS2,MDS1,size=CellEntry,color=TypeEvent,shape=Region))+geom_point()
ggplot(HA1,aes(-1*MDS2,MDS1,size=CellEntry,color=TypeOSC,shape=Region))+geom_point()
ggplot(HA1,aes(-1*MDS2,MDS1,size=CellEntry,color=Dist_to_Receptor,shape=Region))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#EBEBEB", "#000000"))
ggplot(HA1,aes(MDS2,MDS1,size=-1*DMS_ENTROPY,color=Dist_to_Receptor,shape=Region))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(HA1,aes(-1*DMS_ENTROPY,Dist_to_Receptor,size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_hline(yintercept = c(8,44.9))
OSC<-HA1[HA1$SiteOscillating>3,]
P1<-ggplot(OSC,aes(EscapeMut,Dist_to_Receptor))+geom_violin(draw_quantiles = c(0.5))+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P2<-ggplot(OSC,aes(EscapeMut,WCN_VirusOnly))+geom_violin(draw_quantiles = c(0.5))+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
P3<-ggplot(OSC,aes(EscapeMut,-1*DMS_ENTROPY))+geom_violin(draw_quantiles = 0.5)+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=KeyAntigenic))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
cowplot::plot_grid(P1,P2,P3,ncol = 3)

#基于结构的预测
SARS2 <- read.delim("G:/ViralStruct/StructurePredictESC/SARS2.txt")
Surface<-SARS2[SARS2$RSAType=="Outside",]
RBD<-SARS2[SARS2$Region=="RBD",]
#Evescape与DMS
P1<-ggplot(SARS2,aes(evescape,EscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P2<-ggplot(SARS2,aes(evescape,SeraEscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P3<-ggplot(RBD,aes(evescape,SeraEscapeMax))+geom_point(size=3,shape=15,aes(color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P4<-ggplot(SARS2,aes(SASA_Sm*log10(ClashA+2*ClashC+1),EscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P5<-ggplot(SARS2,aes(SASA_Sm*log10(ClashA+2*ClashC+1),SeraEscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P6<-ggplot(RBD,aes(SASA_Sm*log10(ClashA+2*ClashC+1),SeraEscapeMax))+geom_point(size=3,shape=15,aes(color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P1,P2,P3,P4,P5,P6,ncol = 3)
#二者合并改进预测效果
P2<-ggplot(SARS2,aes((1/(1+exp(-1*evescape)))*(1/(1+exp(-1*(ClashA+2*ClashC+1)))),SeraEscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P1<-ggplot(SARS2,aes((1/(1+exp(-1*evescape)))*(1/(1+exp(-1*(ClashA+2*ClashC+1)))),EscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P1,P2,ncol = 2)

#预测结果的结构分布
ggplot(RBD,aes(-1*MDS2,-1*MDS1,color=SASA_Sm*log10(ClashA+2*ClashC+1),shape=Region,size=SeraEscapeMax))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#000000", "#860000", "#FF0000"))
ggplot(HA1,aes(-1*MDS2,MDS1,color=SASA_Sm*log10(ClashA+1),shape=Region,size=SeraEscape))+geom_point(alpha=0.8)+scale_color_gradientn(colors = c("#000000", "#860000", "#FF0000"))


#H3N2预测
H3N2 <- read.delim("G:/ViralStruct/StructurePredictESC/H3N2.txt")
HA1<-H3N2[H3N2$Region=="HA1",]
P1<-ggplot(H3N2,aes(Evescape,DMSEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P2<-ggplot(H3N2,aes(Evescape,SeraEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P3<-ggplot(HA1,aes(Evescape,SeraEscape))+geom_point(size=3,aes(color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P4<-ggplot(H3N2,aes(SASA_Sm*log10(ClashA+1),DMSEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P5<-ggplot(H3N2,aes(SASA_Sm*log10(ClashA+1),SeraEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P6<-ggplot(HA1,aes(SASA_Sm*log10(ClashA+1),SeraEscape))+geom_point(size=3,aes(color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P1,P2,P3,P4,P5,P6,ncol = 3)
#二者合并改进预测效果
P2<-ggplot(H3N2,aes((1/(1+exp(-1*Evescape)))*(1/(1+exp(-1*(ClashA+1)))),SeraEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P1<-ggplot(H3N2,aes((1/(1+exp(-1*Evescape)))*(1/(1+exp(-1*(ClashA+1)))),DMSEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P1,P2,ncol = 2)

#预测竞争性与进化
MergeFreqIndices <- read.delim("G:/ViralStruct/SARS2/MergeFreqIndices.txt")
RBD<-MergeFreqIndices[MergeFreqIndices$Region=="RBD",]
RBDsort<-RBD[order(RBD$Diversity, decreasing = FALSE), ]
MergeH3N2 <- read.delim("G:/ViralStruct/OtherVirus/H3N2/MergeH3N2.txt")
HA1<-MergeH3N2[MergeH3N2$Region=="HA1",]
HA1sort<-HA1[order(HA1$Diveristy, decreasing = FALSE), ]
P1<-ggplot(RBDsort,aes(SiteCompete,Count))+geom_point(size=3,aes(color=Diversity))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#FF0000"))
P2<-ggplot(HA1sort,aes(SiteCompete,MutEvent))+geom_point(size=3,aes(color=Diveristy))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#FF0000"))
cowplot::plot_grid(P1,P2,ncol = 2)

P3<-ggplot(RBD,aes(SiteCompete,EscapeMax))+geom_point(size=3,aes(color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P4<-ggplot(HA1,aes(SiteCompete,SeraEscape))+geom_point(size=3,aes(color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P3,P4,ncol = 2)

#SARS非典预测新冠
SARS2 <- read.delim("G:/ViralStruct/StructurePredictESC/SARS2.txt")
Surface<-SARS2[SARS2$RSAType=="Outside",]
RBD<-SARS2[SARS2$Region=="RBD",]
P4<-ggplot(SARS2,aes(SARS_SASA_Sm*log10(SARS_ClashSum+1),EscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P5<-ggplot(SARS2,aes(SARS_SASA_Sm*log10(SARS_ClashSum+1),SeraEscapeMax))+geom_point(size=3,aes(shape=Region,color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P6<-ggplot(RBD,aes(SARS_SASA_Sm*log10(SARS_ClashSum+1),SeraEscapeMax))+geom_point(size=3,shape=15,aes(color=WCN_8hri))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P4,P5,P6,ncol = 3)

P7<-ggplot(SARS2,aes(SARS_SASA_Sm*log10(SARS_ClashSum+1),SASA_Sm*log10(ClashA+2*ClashC+1)))+geom_point(size=3,aes(shape=Region,color=SeraEscapeMax))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#90D42D"))
P8<-ggplot(HA1,aes(H1_SASA_Lg*log10(H1_Clash+1),SASA_Sm*log10(ClashA+1)))+geom_point(size=3,aes(color=SeraEscape))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#90D42D"))
cowplot::plot_grid(P7,P8,ncol = 2)

PA<-ggplot(RBD,aes(Pos,SARS_SASA_Sm*log10(SARS_ClashSum+1),label=Pos))+geom_line()+geom_text()
PB<-ggplot(RBD,aes(Pos,SASA_Sm*log10(ClashA+2*ClashC+1),label=Pos))+geom_line()+geom_text()
cowplot::plot_grid(PB,PA,ncol = 1)
ggplot(RBD,aes(SASA_Sm*log10(ClashA+2*ClashC+1),SARS_SASA_Sm*log10(SARS_ClashSum+1)))+geom_point()+stat_cor()

#H1N1预测H3N2
H3N2 <- read.delim("G:/ViralStruct/StructurePredictESC/H3N2.txt")
HA1<-H3N2[H3N2$Region=="HA1",]
P4<-ggplot(H3N2,aes(H1_SASA_Lg*log10(H1_Clash+1),DMSEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P5<-ggplot(H3N2,aes(H1_SASA_Lg*log10(H1_Clash+1),SeraEscape))+geom_point(size=3,aes(shape=Region,color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
P6<-ggplot(HA1,aes(H1_SASA_Lg*log10(H1_Clash+1),SeraEscape))+geom_point(size=3,aes(color=WCN_VirusOnly))+geom_smooth(method = "lm")+stat_cor()+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
cowplot::plot_grid(P4,P5,P6,ncol = 3)

tmp1<-HA1[HA1$Pos>100,]
tmp2<-tmp1[tmp1$Pos<250,]
PA<-ggplot(tmp2,aes(Pos,H1_SASA_Lg*log10(H1_Clash+1),label=Pos))+geom_line()
PB<-ggplot(tmp2,aes(Pos,SASA_Sm*log10(ClashA+1),label=Pos))+geom_line()
cowplot::plot_grid(PB,PA,ncol = 1)
ggplot(tmp2,aes(H1_SASA_Lg*log10(H1_Clash+1),SASA_Sm*log10(ClashA+1)))+geom_point()+stat_cor()

#Poliovirus1
Poliovirus <- read.delim("G:/ViralStruct/StructurePredictESC/Poliovirus.txt")
Comp<-Poliovirus[Poliovirus$Competitiveness>0,]
ggplot(Comp,aes(Competitiveness,FoldChange,label=Strain))+geom_point(aes(color=Virus,size=3))+stat_cor()+geom_smooth(method = "lm")+geom_text()


P1.XYZ <- read.delim("G:/ViralStruct/OtherVirus/Poliovirus1/Data/PDB/P1.XYZ.txt", header=FALSE)
MDS<-cmdscale(dist(P1.XYZ))
write.csv(MDS,"G:/ViralStruct/OtherVirus/Poliovirus1/Data/PDB/P1.XYZ.MDS.CSV")
MergePoliovirus1 <- read.delim("G:/ViralStruct/OtherVirus/Poliovirus1/MergePoliovirus1.txt")
ggplot(MergePoliovirus1,aes(MDS1,MDS2,size=-1*Dist_to_Receptor,color=Type,shape=Region))+geom_point()
ggplot(MergePoliovirus1,aes(MDS1,MDS2,size=-1*Dist_to_Receptor,color=Dist_to_Receptor,shape=Region))+geom_point()+scale_color_gradientn(colors = c("#EBEBEB", "#000000"))
#结构展示
MergePoliovirus1 <- read.delim("G:/ViralStruct/OtherVirus/Poliovirus1/MergePoliovirus1.txt")
ggplot(MergePoliovirus1,aes(MDS1,MDS2,color=WCN_VirusOnly))+geom_point(size=3,alpha=0.8)+scale_color_gradientn(colors = c("#000000", "#E7B800", "#00919B"))
ggplot(MergePoliovirus1,aes(MDS1,MDS2,color=Dist_to_Receptor))+geom_point(size=3,alpha=0.8)+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#90D42D"))
ggplot(MergePoliovirus1,aes(MDS1,MDS2,color=Region))+geom_point(size=3,alpha=0.8)
ggplot(MergePoliovirus1,aes(MDS1,MDS2,color=EscapeMax,shape=EscapeMut,size=SiteOscillating))+geom_point()+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
#功能
compaired <- list(c("Other", "ESCOSC"),c("ESCOSC", "ESCNo"))
ggplot(MergePoliovirus1, aes(Type, Dist_to_Receptor,color = Type)) + geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
ggplot(MergePoliovirus1,aes(Dist_to_Receptor,EscapeMax))+geom_point()+geom_smooth(method = "loess")
#预测
OSC<-MergePoliovirus1[MergePoliovirus1$SiteOscillating>1,]
ggplot(OSC,aes(EscapeMut,WCN_VirusOnly))+geom_violin(draw_quantiles = c(0.5))+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=EscapeMut))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(OSC,aes(EscapeMut,Dist_to_Receptor))+geom_violin(draw_quantiles = c(0.5))+geom_point(aes(size=SiteOscillating,color=EscapeMax,shape=EscapeMut))+scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))

#不同病毒的合并研究
#逃逸热点到受体的距离
AllVirusESCMut <- read.delim("G:/ViralStruct/OtherVirus/AllVirusESCMut.txt")
compaired <- list(c("SARS2", "H3N2"),c("H1N1pdm", "H3N2"),c("H1N1pdm", "PV"))
ggplot(AllVirusESCMut,aes(Virus,Dist_to_Receptor))+geom_boxplot()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)

#计算突变率
WindowMut_SARS2 <- read.delim("G:/ViralStruct/SARS2/EvoSpeed/WindowMut.txt")
lm(Contact~Window,data=WindowMut_SARS2)
lm(Close~Window,data=WindowMut_SARS2)
lm(Far~Window,data=WindowMut_SARS2)
lm((Contact+Close+Far)~Window,data=WindowMut_SARS2)

WindowMut_H1N1pdm <- read.delim("G:/ViralStruct/OtherVirus/H1N1pdm/EvoSpeed/WindowMut.txt")
lm(Contact~Window,data=WindowMut_H1N1pdm)
lm(Close~Window,data=WindowMut_H1N1pdm)
lm(Far~Window,data=WindowMut_H1N1pdm)
lm((Contact+Close+Far)~Window,data=WindowMut_H1N1pdm)
ggplot(WindowMut_H1N1pdm,aes(Window,Contact))+geom_point()+geom_smooth(method = "lm")

WindowMut_H3N2 <- read.delim("G:/ViralStruct/OtherVirus/H3N2/EvoSpeed/WindowMut.txt")
lm(Contact~Window,data=WindowMut_H3N2)
lm(Close~Window,data=WindowMut_H3N2)
lm(Far~Window,data=WindowMut_H3N2)
lm((Contact+Close+Far)~Window,data=WindowMut_H3N2)
ggplot(WindowMut_H3N2,aes(Window,Contact))+geom_point()+geom_smooth(method = "lm")

WindowMut_PV1 <- read.delim("G:/ViralStruct/OtherVirus/Poliovirus1/EvoSpeed/WindowMut.txt")
lm(Contact~Window,data=WindowMut_PV1)
lm(Close~Window,data=WindowMut_PV1)
lm(Far~Window,data=WindowMut_PV1)
lm((Contact+Close+Far)~Window,data=WindowMut_PV1)
ggplot(WindowMut_PV1,aes(Window,Contact))+geom_point()+geom_smooth(method = "lm")

#进化速度
EvoSpeed <- read.delim("G:/ViralStruct/OtherVirus/EvoSpeed/EvoSpeed.txt")
cor.test(EvoSpeed$MedianDist,EvoSpeed$EvoSpeed)
ggplot(EvoSpeed,aes(MedianDist,EvoSpeed))+geom_point()
ggplot(EvoSpeed,aes(Virus,log10(EvoSpeed),fill=Region,group=Region))+geom_bar(stat = "identity",position = "dodge")
ggplot(EvoSpeed,aes(Region,log10(EvoSpeed),color=Virus,group=Virus))+geom_point(size=3)+geom_line()

#表面site
RSASurface <- read.delim("G:/ViralStruct/OtherVirus/EvoSpeed/RSASurface.txt")
Top100<-RSASurface[RSASurface$DistRank<=50,]
ggplot(tmp,aes(DistRank,Dist_to_Receptor,group=Virus))+geom_point(aes(color=Type))+geom_line(aes(color=Virus))

D5<-RSASurface[RSASurface$Dist_to_Receptor<=5,]
tmp<-RSASurface[RSASurface$Dist_to_Receptor<=15,]
D10<-tmp[tmp$Dist_to_Receptor>5,]
D15<-RSASurface[RSASurface$Dist_to_Receptor>15,]

compaired <- list(c("SARS", "PV"),c("SARS", "H1N1"),c("SARS", "H3N"))
ggplot(D5,aes(Virus,WCN_8hri,color=Virus))+geom_boxplot()+geom_point()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
ggplot(D10,aes(Virus,WCN_8hri,color=Virus))+geom_boxplot()+geom_point()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)

compaired <- list(c("SARS", "PV"),c("PV", "H1N1"),c("PV", "H3N"))
ggplot(D5,aes(Virus,B_Factor_Z,color=Virus))+geom_boxplot()+geom_point()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)
ggplot(D10,aes(Virus,B_Factor_Z,color=Virus))+geom_boxplot()+geom_point()+ geom_signif(comparisons = compaired, step_increase = 0.1, map_signif_level = F, test = wilcox.test)

ggplot(D5,aes(DistRank,Dist_to_Receptor,group=Virus))+geom_point(aes(size=3,alpha=0.5,color=Type))+geom_line()+facet_grid(vars(Virus),)
ggplot(D10,aes(DistRank,Dist_to_Receptor,group=Virus))+geom_point(aes(size=3,alpha=0.5,color=Type))+geom_line()+facet_grid(vars(Virus),)
ggplot(D15,aes(DistRank,Dist_to_Receptor,group=Virus))+geom_point(aes(size=3,alpha=0.5,color=Type))+geom_line()+facet_grid(vars(Virus),)

#新冠ESC热点位于内测
SARS2 <- read.delim("G:/ViralStruct/OtherVirus/EvoSpeed/SARS2.txt")
ggplot(SARS2,aes(Region,Proportion,fill=Type))+geom_bar(stat = "identity")

#非竞争性抗体对病毒进化的影响
RBDSurface<-RBD[RBD$RSAType=="Outside",]
ggplot(RBDSurface,aes(Dist_to_Receptor,SiteOscillating,color=WTEscape_NoNeu,size=WTEscape_NoNeu))+geom_point()+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#90D42D"))
ggplot(RBDSurface,aes(Dist_to_Receptor,MutRate,color=WTEscape_NoNeu,size=WTEscape_NoNeu))+geom_point()+geom_smooth(method = "gam")+scale_color_gradientn(colors = c("#EBEBEB", "#E7B800", "#90D42D"))

#比对两个驱动力的相互干扰
P3<-ggplot(HA1,aes(SeraEscape,CellEntry))+geom_point()+geom_smooth(method = "lm")+stat_cor()
P4<-ggplot(RBD,aes(SeraEscapeMax,CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()
P1<-ggplot(HA1,aes(DMSEscape,CellEntry))+geom_point()+geom_smooth(method = "lm")+stat_cor()
P2<-ggplot(RBD,aes(EscapeMax,CellEntryMax))+geom_point()+geom_smooth(method = "lm")+stat_cor()
cowplot::plot_grid(P1,P2,P3,P4,ncol = 2)