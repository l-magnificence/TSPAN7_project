setwd("~/bioinfo_mill/tspan7_proj")
library(dplyr)
library(survminer)
library(survival)
library(ggplot2)
library(ggstatsplot)
library(plyr)

###-----------------PRJNA482620 datasets-----------------
clinical<-readRDS(file = "~/bioinfo_mill/dataset/glioma/ICB_corhort/PRJNA482620_treat_before_meta.rds")
exprSet<-readRDS(file = "~/bioinfo_mill/dataset/glioma/ICB_corhort/PRJNA482620_treat_before_expr.rds")

identical(colnames(exprSet),rownames(clinical))

clinical$gene<-t(exprSet["TSPAN7",])[,1]
clinical$group<-ifelse(clinical[,"gene"]>median(clinical[,"gene"]),"High","Low")
clinical$response<-clinical$`Response to PD1?`
clinical[clinical$`Response to PD1?`=="Yes","response"]<-"R"
clinical[clinical$`Response to PD1?`=="No","response"]<-"NR"
clinical$response<-factor(clinical$response,levels = c("R","NR"))

##--------tspan7 response --------
a <- data.frame(table(clinical$group,clinical$response))
a<- ddply(a, .(Var1), transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")

fisher.test(table(clinical$group,clinical$response))$p.value

pdf(file="./result/immuotherapy/PRJNA482620_TSPAN7_immunotherapy.pdf",bg="white",width=3,height=4.5,pointsize=12)
ggplot(a,aes(x = Var1,y=percent,fill = Var2))+
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 16)+
  geom_text(mapping = aes(label = label),
            size = 5, colour = 'black', vjust = 3, hjust = .5, position = position_stack())+
  ggsci::scale_fill_nejm()+ 
  theme_classic()+
  labs(title="Anti-PD-1 immunotherapy", x="TSPAN7 expression", y="Percentage")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position="right",
        legend.title=element_blank())
dev.off()

## survival in risk group
dat_select<-clinical
mySurv<-Surv(dat_select$`OS (days) from PD1 Inhibitor`, dat_select$`Dead (1) or alive (0)`)
group <-dat_select$group
group <- factor(group, levels = c("Low", "High"))
survival_dat <- data.frame(group = group)
fit <- survfit(mySurv ~ group)
#计算HR以及95%CI
data.survdiff <- survdiff(mySurv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- dat_select[order(dat_select[,"gene"]),]

pdf(file="./result/immuotherapy/PRJNA482620_TSPAN7_survival.pdf",bg="white",width=5.5,height=6.5,pointsize=12)
ggsurvplot(fit, data = survival_dat ,
           #ggtheme = theme_bw(), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           censor = F, #不显示观察值所在的位置
           palette = "aaas", #线的颜色对应高、低
           pval = T,
           legend.title = paste0("Group"),#基因名写在图例题目的位置
           font.legend = 11,#图例的字体大小
           #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
           xlab = "Survival after anti-PD-1 therapy (Months)",
           #在图例上标出高低分界点的表达量，和组内sample数量
           legend.labs=c(paste0("TSPAN7 Low expression"),
                         paste0("TSPAN7 High expression")),
           risk.table=T,risk.table.y.text=F)
dev.off()

##--------PDL1 response --------
clinical$PDL1<-t(exprSet["CD274",])[,1]
clinical$group2<-ifelse(clinical[,"PDL1"]>median(clinical[,"PDL1"]),"High","Low")

a <- data.frame(table(clinical$group2,clinical$response))
a<- ddply(a, .(Var1), transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")

fisher.test(table(clinical$group2,clinical$response))$p.value

pdf(file="./result/immuotherapy/PRJNA482620_PDL1_immunotherapy.pdf",bg="white",width=3,height=4.5,pointsize=12)
ggplot(a,aes(x = Var1,y=percent,fill = Var2))+
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 16)+
  geom_text(mapping = aes(label = label),
            size = 5, colour = 'black', vjust = 3, hjust = .5, position = position_stack())+
  ggsci::scale_fill_nejm()+ 
  theme_classic()+
  labs(title="Anti-PD-1 immunotherapy", x="PDL1 expression", y="Percentage")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position="right",
        legend.title=element_blank())
dev.off()

## survival in risk group
dat_select<-clinical
mySurv<-Surv(dat_select$`OS (days) from PD1 Inhibitor`, dat_select$`Dead (1) or alive (0)`)
group <-dat_select$group2
group <- factor(group, levels = c("Low", "High"))
survival_dat <- data.frame(group = group)
fit <- survfit(mySurv ~ group)
#计算HR以及95%CI
data.survdiff <- survdiff(mySurv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- dat_select[order(dat_select[,"gene"]),]

pdf(file="./result/immuotherapy/PRJNA482620_PDL1_survival.pdf",bg="white",width=5.5,height=6.5,pointsize=12)
ggsurvplot(fit, data = survival_dat ,
           #ggtheme = theme_bw(), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           censor = F, #不显示观察值所在的位置
           palette = "aaas", #线的颜色对应高、低
           pval = T,
           legend.title = paste0("Group"),#基因名写在图例题目的位置
           font.legend = 11,#图例的字体大小
           #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
           xlab = "Survival after anti-PD-1 therapy (Months)",
           #在图例上标出高低分界点的表达量，和组内sample数量
           legend.labs=c(paste0("PDL1 Low expression"),
                         paste0("PDL1 High expression")),
           risk.table=T,risk.table.y.text=F)
dev.off()

##--------TSPAN7 PDL1 response --------
df<-clinical
df$TSPAN7_con <- cut(df$gene,breaks=c(-Inf,median(df$gene), Inf),
                      labels=c("0","1"))
df$PDL1_con <- cut(df$PDL1,breaks=c(-Inf,median(df$PDL1), Inf),
                   labels=c("0","1"))
df$medGroup <- paste0(df[,"TSPAN7_con"],df[,"PDL1_con"])
df$medGroup <- ifelse(df$medGroup=='10','1',ifelse(df$medGroup=='00','2',ifelse(df$medGroup=="01","3","4")))
df[df$medGroup==1,"medGroup"]<-"TSPAN7 high & PD-L1 low"
df[df$medGroup==2,"medGroup"]<-"TSPAN7 low & PD-L1 low"
df[df$medGroup==3,"medGroup"]<-"TSPAN7 low & PD-L1 high"
df[df$medGroup==4,"medGroup"]<-"TSPAN7 high & PD-L1 high"

#TSPAN7 PDL1 survival
fit <- survfit(mySurv ~ df$medGroup)

pdf(file="./result/immuotherapy/PRJNA482620_PDL1_TSPAN7_survival.pdf",bg="white",width=5.5,height=6.5,pointsize=12)
ggsurvplot(fit, data = df,pval = T,palette="lancet",
           legend =  c(0.24,0.3),legend.title="Group",
           legend.labs=c("TSPAN7 high & PD-L1 high","TSPAN7 high & PD-L1 low","TSPAN7 low & PD-L1 high","TSPAN7 low & PD-L1 low"),
           risk.table=T,risk.table.y.text=F,
           xlab = "Survival after anti-PD-1 therapy (Months)",
           pval.coord=c(10,0.05))
dev.off()

#TSPAN7 PDL1 Group
a <- data.frame(table(df$medGroup,df$response))
a<- ddply(a, .(Var1), transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")

fisher.test(table(df$medGroup,df$response))$p.value

pdf(file="./result/immuotherapy/PRJNA482620_PDL1_TSPAN7_immunotherapy.pdf",bg="white",width=4.5,height=4.5,pointsize=12)
ggplot(a,aes(x = Var1,y=percent,fill = Var2))+
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 16)+
  geom_text(mapping = aes(label = label),
            size = 5, colour = 'black', vjust = 3, hjust = .5, position = position_stack())+
  ggsci::scale_fill_nejm()+ 
  theme_classic()+
  labs(title="Anti-PD-1 immunotherapy", x="", y="Percentage")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =8),
        plot.title = element_text(hjust = 0.5),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position="right",
        legend.title=element_blank())
dev.off()

