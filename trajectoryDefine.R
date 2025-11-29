rm(list=ls())
library(haven)
library(openxlsx)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(lcmm)        
library(nlme)    
library(ggplot2)
options(scipen = 200)
#数据整理
z_score <- read_sav("E:z_score250528.sav")
data1 <- z_score %>%
  dplyr::mutate(time=case_when(agemos==0~1,
                               agemos>=5 & agemos<=7~2,
                               agemos>=11 & agemos<=13~3,
                               agemos>=17 & agemos<=19~4,
                               agemos>=23 & agemos<=25~5,
                               agemos>=29 & agemos<=31~6,
                               agemos>=35 & agemos<=37~7)) %>%
  drop_na(time) 
# 处理函数_保留最靠近年龄的一次
filter_closest_age <- function(data) {
  data %>%
    group_by(IDchild, time) %>%
    mutate(
      age_diff = case_when(
        time == 2 ~ abs(agemos - 6),
        time == 3 ~ abs(agemos - 12),
        time == 4 ~ abs(agemos - 18),
        time == 5 ~ abs(agemos - 24),
        time == 6 ~ abs(agemos - 30),
        time == 7 ~ abs(agemos - 36),
        TRUE ~ Inf  
      )
    ) %>%
    # 按age_diff排序，选择最小的(最接近目标年龄)
    arrange(IDchild, time, age_diff) %>%
    group_by(IDchild, time) %>%
    slice(1) %>%
    ungroup() %>%
    select(-age_diff) 
}  
data2 <- filter_closest_age(data1) 
# 计算每个ID的测量次数
id_counts <- data2 %>%
  group_by(IDchild) %>%
  summarise(n_visits = n()) %>%
  filter(n_visits > 2)  %>%# 保留至少3次测量的ID
  dplyr::distinct(IDchild) %>%
  unlist()
data4 <- data2 %>%
  dplyr::filter(IDchild %in% id_counts) %>%
  dplyr::mutate(IDchild=as.numeric(IDchild))
#组轨迹模型(Group-based trajectory modelling, GBTM)
#线性(未纳入表格)
#m1_1<-hlme(zwhz~1+time, random=~1+time, ng=1, data=data4, subject="IDchild")
#m2_1<-hlme(zwhz~1+time, mixture=~1+time, random=~1+time, ng=2, data=data4, subject="IDchild",B=m1_1)
#m3_1<-hlme(zwhz~1+time, mixture=~1+time, random=~1+time, ng=3, data=data4, subject="IDchild",B=m1_1)
#m4_1<-hlme(zwhz~1+time, mixture=~1+time, random=~1+time, ng=4, data=data4, subject="IDchild",B=m1_1)
#m5_1<-hlme(zwhz~1+time, mixture=~1+time, random=~1+time, ng=5, data=data4, subject="IDchild",B=m1_1)
#Quadratic
m1_2<-hlme(zwhz~1+time+I(time^2), random=~1+time, ng=1, data=data4, subject="IDchild")
m2_2<-hlme(zwhz~1+time+I(time^2), mixture=~1+time+I(time^2), random=~1+time, ng=2, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_2)
m3_2<-hlme(zwhz~1+time+I(time^2), mixture=~1+time+I(time^2), random=~1+time, ng=3, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_2)
m4_2<-hlme(zwhz~1+time+I(time^2), mixture=~1+time+I(time^2), random=~1+time, ng=4, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_2)
m5_2<-hlme(zwhz~1+time+I(time^2), mixture=~1+time+I(time^2), random=~1+time, ng=5, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_2)
#Cubic
m1_3 <- hlme(zwhz~1+time+I(time^2)+I(time^3), random=~1+time, ng=1, data=data4, subject="IDchild")
m2_3 <- hlme(zwhz~1+time+I(time^2)+I(time^3), mixture=~1+time+I(time^2)+I(time^3), random=~1+time, ng=2, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_3)
m3_3 <- hlme(zwhz~1+time+I(time^2)+I(time^3), mixture=~1+time+I(time^2)+I(time^3), random=~1+time, ng=3, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_3)
m4_3 <- hlme(zwhz~1+time+I(time^2)+I(time^3), mixture=~1+time+I(time^2)+I(time^3), random=~1+time, ng=4, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_3)
m5_3 <- hlme(zwhz~1+time+I(time^2)+I(time^3), mixture=~1+time+I(time^2)+I(time^3), random=~1+time, ng=5, nwg=TRUE, idiag=FALSE, data = data4, subject = "IDchild", B=m1_3)
#proportion per class%
model_summary1 <- as.data.frame(summarytable(m1_2,m2_2,m3_2,m4_2,m5_2,m1_3,m2_3,m3_3,m4_3,m5_3)) %>%
  rownames_to_column(.,var="Model") %>%
  dplyr::rename(class1="%class1",
                class2="%class2",
                class3="%class3",
                class4="%class4",
                class5="%class5",
                Group=G) %>%
    dplyr::mutate(across(6:10, ~ round(., 1)),
                  Class = apply(cbind(class1, class2, class3, class4, class5), 1, 
                                function(x) paste(na.omit(x), collapse = "-")),
                  Parameter=case_when(grepl("_2",Model)~"Quadratic",
                                      grepl("_3",Model)~"Cubic")
                  ) 
  
# APP
get_avg_post_prob <- function(model, model_name) {
  post_probs <- model$pprob
  n_classes <- max(post_probs$class)
  
  post_probs %>%
    group_by(class) %>%
    summarise(across(paste0("prob", 1:n_classes), mean)) %>%
    mutate(Model = model_name)
}

model_names <- c("m1_2","m2_2","m3_2","m4_2","m5_2","m1_3","m2_3","m3_3","m4_3","m5_3")
model_list <- list(m1_2,m2_2,m3_2,m4_2,m5_2,m1_3,m2_3,m3_3,m4_3,m5_3)

model_summary2 <- map2_dfr(model_list, model_names, get_avg_post_prob) %>%
  dplyr::select(3,1,everything()) %>%
  dplyr::mutate(across(3:7, ~ round(., 2)),
                prob=case_when(class==1~prob1,
                               class==2~prob2,
                               class==3~prob3,
                               class==4~prob4,
                               class==5~prob5),
                prob_class=paste0("prob",class)) %>%
  dplyr::select(1,9,8) %>%
  pivot_wider(names_from = prob_class,values_from = prob) %>%
  dplyr::mutate(APP=apply(cbind(prob1, prob2, prob3, prob4, prob5), 1, 
                            function(x) paste(na.omit(x), collapse = "-"))) 

# Ek
entropy <- function(m) {
  if (!"pprob" %in% names(m)) stop("Model does not contain posterior probabilities!")
  p <- as.matrix(m$pprob[, grep("prob", names(m$pprob))])  
  if (ncol(p) < 1) stop("No probability columns found!")
  entropy_val <- -sum(p * log(p + 1e-12), na.rm = TRUE)
  n <- nrow(p)
  k <- ncol(p)
  rel_entropy <- 1 - (entropy_val / (n * log(k)))
  return(pmax(0, pmin(1, rel_entropy))) 
} (-sum(p * log(p + 1e-12)) / (nrow(p) * log(ncol(p))))

model_summary3 <- data.frame(Model = model_names,Rel_Entropy = sapply(model_list, entropy))
#合并BIC,CLASS%,APP,EK
model_summary <- merge(model_summary1,model_summary2,by="Model") %>%
  merge(.,model_summary3,by="Model") %>%
  dplyr::select(1,2,12,5,11,18,19) 
write.xlsx(model_summary,"model_summary.xlsx")

#基于m3_3导出分类数据
data_traj <- m3_3$pprob %>% 
  dplyr::select(1,2)
table(data_traj$class)
data <- data4 %>%
  left_join(data_traj,by="IDchild") %>%
  dplyr::mutate(IDchild=as.character(IDchild)) 
write_sav(data,"E:data_traj.sav")

# 2. 计算每组的平均轨迹
mean_traj <- data.frame(aggregate(zwhz ~ time + class, data = data, FUN = mean)) %>%
  dplyr::mutate(class=case_when(class==1~"Trajectory1",
                                class==2~"Trajectory2",
                                class==3~"Trajectory3"),
                time=case_when(time==7~36,
                               time==6~30,
                               time==5~24,
                               time==4~18,
                               time==3~12,
                               time==2~6,
                               time==1~0))

plot<- ggplot(data = mean_traj,aes(x=time,y=zwhz))+
  geom_line(aes(color=factor(class)),size=1.5,linetype=1)+
  geom_point(aes(shape = factor(class), color =factor(class)),size = 4) +
  scale_shape_manual(values = c(17,16,15))+
  scale_color_manual(values = c("#e95280","#23b1a5","#ffdd7e"))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        legend.title = element_blank())+
  scale_x_continuous(breaks = seq(0,36,by=6)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1))+
  labs(x="Months",y="WFL/Hz")+
  theme(axis.title=element_text(face="plain",colour="black",size=18))+
  theme(axis.text=element_text(face="plain",colour="black",size=18))+
  theme(legend.position = "bottom",legend.text=element_text(face = "plain",size = 18,color="black"))
plot
ggsave(plot = plot,
       filename ='E:\\0.妇幼\\0.论文\\3.肠道菌群_华大\\分析\\生长轨迹_QR\\250624\\FigureS.pdf',
       width =8,
       height =8,
)
#newdata <- as.data.frame(data4$time)
#names(newdata) <-c("time")
#plotpred <- predictY(m3_3, newdata, var.time ="time", draws = F)
#a <- plot(plotpred, lty=2,xlab="Times", ylab="WFLz", legend.loc = "topleft", cex=0.75)
