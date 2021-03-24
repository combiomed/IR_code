## assessing model performance and feature ranking

library(tidyverse)
library(dplyr)
library(purrr)
library(colorspace)
library(circlize)
library(ggthemes)
library(ggplot2)
library(pROC)

## data.files contain cleaned datasets for modeling built as described in previous steps
data.files<-c("BlMo.model.data.global.RData",
              "BlMa.model.data.global.RData",
              "BlTN.model.data.global.RData",
              "BlCM.model.data.global.RData",
              "BlEM.model.data.global.RData")
global.data<-lapply(data.files, function (x) get(load(x)))


model.data<-list()
x<-list()
y<-list()
data.down<-list()
data.down.matrix<-list()
data.down.set<-list()
data.down_train<-list()
data.down_test<-list()
data.down_trainX<-list()
data.down_testX<-list()
testy<-list()
subset<-list()
for (i in 1:length(global.data)){
  model.data[[i]]<-global.data[[i]][,c(39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,27,29,30,32,33,35,75:80,85,86,88,89,91,92,16,18:26,96)]
  x[[i]]<-model.data[[i]][,-ncol(model.data[[i]])]
  y[[i]]<-as.factor(model.data[[i]][,ncol(model.data[[i]])])
  set.seed(10857)
  data.down[[i]]<-downSample(x[[i]], y[[i]], list = FALSE, yname = "IntronType")
  data.down[[i]]<-na.omit(data.down[[i]])
  data.down.matrix[[i]]<-model.matrix( ~ -1 + ., data.down[[i]][,-ncol(model.data[[i]])])#model.matrix() allowing to automatically transform any qualitative variables (if any) into dummy variables, which is important because GLMNET/FittedModels/net() can only take numerical, quantitative inputs.
  #After creating the model matrix, we remove the intercept component at index = 1.
  data.down.set[[i]]<-data.frame(cbind(data.down.matrix[[i]], IntronType = as.factor(data.down[[i]][,ncol(model.data[[i]])])))
  data.down.set[[i]]$IntronType<-factor(data.down.set[[i]]$IntronType, levels = c(1,2), labels = c("Non.Retained", "Retained"))
  set.seed(10857)##set seed
  subset[[i]] = sample(nrow(data.down.set[[i]]), nrow(data.down.set[[i]]) * 0.80)
  data.down_train[[i]] = data.down.set[[i]][subset[[i]], ]
  data.down_test[[i]] = data.down.set[[i]][-subset[[i]], ]
  set.seed(849)
  data.down_trainX[[i]]<-as.matrix(data.down_train[[i]][,-ncol(data.down_train[[i]])])
  data.down_testX[[i]]<-data.down_test[[i]][,-ncol(data.down_test[[i]])]
  y[[i]] <- data.down_train[[i]][,ncol(data.down_train[[i]])]
  testy[[i]]<-data.down_test[[i]][,ncol(data.down_test[[i]])]
}

##loading glmnet models
glmnet.intrinall.files<-c("BlMo.glmnet.model.intrinsicall.cv10.RData",
                           "BlMa.glmnet.model.intrinsicall.cv10.RData",
                           "BlTN.glmnet.model.intrinsicall.cv10.RData",
                           "BlCM.glmnet.model.intrinsicall.cv10.RData",
                           "BlEM.glmnet.model.intrinsicall.cv10.RData")
glmnet.intrinall<-lapply(glmnet.intrinall.files, function (x) get(load(x)))

## assessing confusion matrices and roc values
pred.test<-list()
pred.obs<-list()
glmnet.result.roc<-list()
for (i in 1:length(model.data)){
  pred.test[[i]]=predict(glmnet.intrinall[[i]], data.down_testX[[i]], type = "prob") # make prediction for test set
  pred.test[[i]]$IntronType<-ifelse(pred.test[[i]]$Non.Retained>0.5, 0, 1)
  pred.obs[[i]]=data.frame(pred.test=pred.test[[i]]$IntronType, testy[[i]])
  table(pred.obs[[i]]) ##confusion matrix
  glmnet.result.roc[[i]] <- roc(response = testy[[i]], predictor = pred.test[[i]]$Retained, levels = levels(testy[[i]]) )
}

##loading conditional random forest models
cforest.intrinall.files<-c("BlMo.cforest.model.intrinsicall.cv10.RData",
                    "BlMa.cforest.model.intrinsicall.cv10.RData",
                    "BlTN.cforest.model.intrinsicall.cv10.RData",
                    "BlCM.cforest.model.intrinsicall.cv10.RData",
                    "BlEM.cforest.model.intrinsicall.cv10.RData")
cforest.intrinall<-lapply(cforest.intrinall.files, function (x) get(load(x)))
names(cforest.intrinall)<-c("BlMo", "BlMa", "BlTN", "BlCM", "BlEM")

pred.test<-list()
pred.obs<-list()
cforest.result.roc<-list()
for (i in 1:length(model.data)){
  pred.test[[i]]=predict(cforest.intrinall[[i]], data.down_testX[[i]], type = "prob") # make prediction for test set
  pred.test[[i]]$IntronType<-ifelse(pred.test[[i]]$Non.Retained>0.5, 0, 1)
  pred.obs[[i]]=data.frame(pred.test=pred.test[[i]]$IntronType, testy[[i]])
  table(pred.obs[[i]]) ##confusion matrix
  cforest.result.roc[[i]] <- roc(response = testy[[i]], predictor = pred.test[[i]]$Retained, levels = levels(testy[[i]]) )
}

## plot for Figure 2B
## data_2, data_3, data_4 and data_5 represents different cell types
data_1<-as.data.frame(glmnet.matched.result.roc[[1]]$sensitivities)
data_1$FDR<-1-glmnet.matched.result.roc[[1]]$specificities
data_1$model<-"Elastic Net"
data_1$AUC<-round(glmnet.matched.result.roc[[1]]$auc, 4)
colnames(data_1)[1]<-"Sensitivity"

data<-rbind(data_1, data_2, data_3, data_4, data_5)
data$`Model (AUC)`<-paste0(data$model, ' (', data$AUC, ')')
data$`Model (AUC)`<-factor(data$`Model (AUC)`, levels = c("Monocytes (0.92)", "Macrophages (0.87)", "T Naive Memory (0.95)", "T Central Memory (0.95)", "T Effector Memory (0.95)"))
colnames(data)[2]<-"False Positive Rate"
ggplot(data, aes(x=`False Positive Rate`, y=Sensitivity, colour=`Model (AUC)`)) +
  geom_line(size =1) +
  theme_light(base_size = 25) +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 25)) +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))+
  ggtitle('AUROC') +
  labs(x ="False Positive Rate", fill = 'Model (AUC)')+
  scale_color_manual(values = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"))+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+
  theme(legend.position=c(0.73, 0.2), legend.box="vertical", legend.title = element_blank(), legend.text = element_text(size = 14))

## plot for Figure 2C
## data_2, data_3, data_4 and data_5 represents different cell types
data_1<-as.data.frame(t(pROC::coords(glmnet.matched.result.roc[[1]], x = "all", input = "threshold", ret = c("recall", "ppv"), transpose = T)))
data_1<-data_1[order(data_1$recall),]
data_1$model<-"Elastic Net"
data_1$AUC<-round(glmnet.matched.result.roc$auc, 2)

data.ppv<-rbind(data_1, data_2, data_3, data_4, data_5)
colnames(data.ppv)[1]<-"Recall"
colnames(data.ppv)[2]<-"Precision"
data.ppv$`Model (AUC)`<-paste0(data.ppv$model, ' (', data.ppv$AUC, ')')
data.ppv$`Model (AUC)`<-factor(data.ppv$`Model (AUC)`, levels = c("Monocytes (0.92)", "Macrophages (0.87)", "T Naive Memory (0.95)", "T Central Memory (0.95)", "T Effector Memory (0.95)"))

ggplot(data.ppv, aes(x=Recall, y=Precision, colour=`Model (AUC)`)) +
  geom_line(size =1) +
  theme_light(base_size = 25) +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 25)) +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), legend.title = element_blank())+
  ggtitle("Precision Recall Curves") +
  labs(fill = "Model (AUC)") +
  ylim(0,1)+
  geom_segment(aes(x = 0, xend = 1, y = 0.5, yend = 0.5), color="darkgrey", linetype="dashed")+
  scale_color_manual(values = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E")) +
  theme(legend.position=c(0.25, 0.25), legend.box="vertical", legend.text = element_text(size = 14))


Cells<-list("1" = "Monocytes", "2"="Macrophages", "3"="T Naive", "4"="T Central Memory", "5"="T Effector Memory")
vim.crf<-list()
vim_importance<-list()
vim_importancetop5<-list()
VIMplot<-list()
for (i in 1:length(cforest.intrinall)){
  vim.crf[[i]]<-varImp(cforest.intrinall[[i]], scale = F)
  vim_importance[[i]]<-vim.crf[[i]]$importance
  vim_importance[[i]]<-vim_importance[[i]][order(-vim_importance[[i]]$Overall), , drop = F]
  vim_importancetop5[[i]]<-vim_importance[[i]][1:5, , drop = F]
  vim_importancetop5[[i]]<-tibble::rownames_to_column(vim_importancetop5[[i]], "Features")
  colnames(vim_importancetop5[[i]])[2]<-Cells[[i]]
}

## Plot for Figure 2E
all.features<- vim_importancetop5 %>% purrr::reduce(full_join, by = "Features")
all.features$Features<-factor(all.features$Features, levels = features)
all.features.long<-melt(all.features)
all.features.long$coeff<-abs(all.features.long$value)
ggplot(all.features.long, aes(x = variable, y=coeff)) +
  theme_light(base_size = 25) +
  theme(legend.position="bottom", legend.text = element_text(size = 15)) +
  geom_bar(aes(fill = Features), position = "fill", stat = "identity") +
  coord_flip() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_fill_manual(values=c("#a85400", "#701c00", "#c48c1c", "#e0c48c", "#004878", "#e0a81c","#d89060", "#f0a800","#7896b4","#96b4d2","#e37870")) #colour order for crf +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

## plots for Supplementary Figures 3C, 3D, 4A, 4C
for (i in 1:5){
  vim_importance[[i]]<-tibble::rownames_to_column(vim_importance[[i]], "Features")
  colnames(vim_importance[[i]])[2]<-Cells[[i]]
}
all.features.crf<- vim_importance %>% purrr::reduce(full_join, by = "Features")
all.features.crf$Sum<-rowSums(all.features.crf[,c(2:6)])
all.features.crf<-all.features.crf[order(-all.features.crf$Sum), , drop = F]
all.features.crf[is.na(all.features.crf)]<-0
vim_plot.heatmap<-melt(all.features.crf[,1:6])
order<-all.features.crf$Features
vim_plot.heatmap$Features <- factor(vim_plot.heatmap$Features, levels = rev(order))
color=colorRampPalette(hcl.colors(7, palette = rev("YlOrBr")))
colorpalette=color(paletteSize)

ggplot(vim_plot.heatmap, aes(x = variable, y = Features, fill = value)) +
    geom_tile() +
    theme_classic(base_size = 25)+
    labs(fill = "Mean Decrease Gini")+
    theme(legend.position = "bottom", legend.text = element_text(size = 14))+
    scale_fill_gradientn(colours = c(colorpalette[paletteSize], colorpalette[paletteSize/2], colorpalette[1]))+
    theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank()) +
    theme(axis.text.x=element_text(angle=90))
