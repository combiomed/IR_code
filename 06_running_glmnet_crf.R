library(caret)
library(party)

load("introns.ret.samp.RData")
load("introns.nonret.samp.RData")

BlMo.Hm01.retint.data<-subset(introns.ret.samp[[[[[[1]]]][,c(1:3,6:9,20,24,62:63,36,25,27,29,22,40:41,44,47:49,51,53,58:59,65:112,114,116,118,120,122,124:139)], IntronDepth>=10 & gene_FPKM>1 & Coverage >= 0.9 & trans_type == "protein_coding" & Length < 10000)
BlMo.Hm01.retint.data$IntronType<-1
BlMo.Hm01.retint.data$Sample<-"BlMo_Hm01"

BlMo.Hm01.nonretint.data<-subset(introns.nonret.samp[[[[[[1]]]][,c(1:3,6:9,20,24,62:63,36,25,27,29,22,40:41,44,47:49,51,53,58:59,65:112,114,116,118,120,122,124:139)], IntronDepth<10 & gene_FPKM>1 & trans_type == "protein_coding" & Length < 10000)
BlMo.Hm01.nonretint.data$IntronType<-0
BlMo.Hm01.nonretint.data$Sample<-"BlMo_Hm01"

BlMo.model.data.global<-rbind(BlMo.Hm01.retint.data, BlMo.Hm01.nonretint.data)
BlMo.model.data.global$Cell<-"BlMo"
BlMo.model.data.global<-BlMo.model.data.global[!is.na(BlMo.model.data.global$MaxEntScore5ss) & !is.na(BlMo.model.data.global$MaxEntScore3ss),]
BlMo.model.data.global$IntronType<-factor(BlMo.model.data.global$IntronType, levels = c(0,1), labels = c("Non.Retained", "Retained"))
BlMo.model.data.global$H3K27ac_olap_100bp5ss<-factor(BlMo.model.data.global$H3K27ac_olap_100bp5ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K27me3_olap_100bp5ss<-factor(BlMo.model.data.global$H3K27me3_olap_100bp5ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K36me3_olap_100bp5ss<-factor(BlMo.model.data.global$H3K36me3_olap_100bp5ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K4me1_olap_100bp5ss<-factor(BlMo.model.data.global$H3K4me1_olap_100bp5ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K4me3_olap_100bp5ss<-factor(BlMo.model.data.global$H3K4me3_olap_100bp5ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K9me3_olap_100bp5ss<-factor(BlMo.model.data.global$H3K9me3_olap_100bp5ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$gnomePeak_olap_100bp5ss<-factor(BlMo.model.data.global$gnomePeak_olap_100bp5ss, levels = c(0,1), labels = c("No", "Peak"))
BlMo.model.data.global$gnomePeak_olap_NFR_100bp5ss<-factor(BlMo.model.data.global$gnomePeak_olap_NFR_100bp5ss, levels = c(0,1), labels = c("No", "Peak"))
BlMo.model.data.global$H3K27ac_olap_100bp3ss<-factor(BlMo.model.data.global$H3K27ac_olap_100bp3ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K27me3_olap_100bp3ss<-factor(BlMo.model.data.global$H3K27me3_olap_100bp3ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K36me3_olap_100bp3ss<-factor(BlMo.model.data.global$H3K36me3_olap_100bp3ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K4me1_olap_100bp3ss<-factor(BlMo.model.data.global$H3K4me1_olap_100bp3ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K4me3_olap_100bp3ss<-factor(BlMo.model.data.global$H3K4me3_olap_100bp3ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K9me3_olap_100bp3ss<-factor(BlMo.model.data.global$H3K9me3_olap_100bp3ss, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$gnomePeak_olap_100bp3ss<-factor(BlMo.model.data.global$gnomePeak_olap_100bp3ss, levels = c(0,1), labels = c("No", "Peak"))
BlMo.model.data.global$gnomePeak_olap_NFR_100bp3ss<-factor(BlMo.model.data.global$gnomePeak_olap_NFR_100bp3ss, levels = c(0,1), labels = c("No", "Peak"))
BlMo.model.data.global$H3K27ac_olap_100bpmid<-factor(BlMo.model.data.global$H3K27ac_olap_100bpmid, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K27me3_olap_100bpmid<-factor(BlMo.model.data.global$H3K27me3_olap_100bpmid, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K36me3_olap_100bpmid<-factor(BlMo.model.data.global$H3K36me3_olap_100bpmid, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K4me1_olap_100bpmid<-factor(BlMo.model.data.global$H3K4me1_olap_100bpmid, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K4me3_olap_100bpmid<-factor(BlMo.model.data.global$H3K4me3_olap_100bpmid, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$H3K9me3_olap_100bpmid<-factor(BlMo.model.data.global$H3K9me3_olap_100bpmid, levels = c(0,1,2), labels = c("No", "Peak", "StrongPeak"))
BlMo.model.data.global$gnomePeak_olap_100bpmid<-factor(BlMo.model.data.global$gnomePeak_olap_100bpmid, levels = c(0,1), labels = c("No", "Peak"))
BlMo.model.data.global$gnomePeak_olap_NFR_100bpmid<-factor(BlMo.model.data.global$gnomePeak_olap_NFR_100bpmid, levels = c(0,1), labels = c("No", "Peak"))

lMa.model.data<-BlMo.model.data.global[,c(39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,96)] ##selecting varibles representing HMs only
x<-BlMo.model.data[,-ncol(BlMo.model.data)]
y<-as.factor(BlMo.model.data[,ncol(BlMo.model.data)])
set.seed(10857)
BlMo.data.down<-downSample(x, y, list = FALSE, yname = "IntronType")
BlMo.data.down<-na.omit(BlMo.data.down)
BlMo.data.down.matrix<-model.matrix( ~ -1 + ., BlMo.data.down[,-ncol(BlMo.model.data)])#model.matrix() allowing to automatically transform any qualitative variables (if any) into dummy variables, which is important because GLMNET/FittedModels/net() can only take numerical, quantitative inputs.
                                                              #After creating the model matrix, we remove the intercept component at index = 1.
BlMo.data.down.set<-data.frame(cbind(BlMo.data.down.matrix, IntronType = as.factor(BlMo.data.down[,ncol(BlMo.model.data)])))
BlMo.data.down.set$IntronType<-factor(BlMo.data.down.set$IntronType, levels = c(1,2), labels = c("Non.Retained", "Retained"))
set.seed(10857)##set seed
subset = sample(nrow(BlMo.data.down.set), nrow(BlMo.data.down.set) * 0.80)
BlMo.data.down_train = BlMo.data.down.set[subset, ]
BlMo.data.down_test = BlMo.data.down.set[-subset, ]
set.seed(849)
BlMo.data.down_trainX<-as.matrix(BlMo.data.down_train[,-ncol(BlMo.data.down_train)])
BlMo.data.down_testX<-BlMo.data.down_test[,-ncol(BlMo.data.down_test)]
y <- BlMo.data.down_train[,ncol(BlMo.data.down_train)]
testy<-BlMo.data.down_test[,ncol(BlMo.data.down_test)]

cctrl1 <- trainControl(method="cv", number=10, returnResamp="all", classProbs=TRUE, summaryFunction=twoClassSummary, verboseIter = T)
set.seed(5478)
BlMo.glmnet.model.down.cv10 <- train(BlMo.data.down_trainX, y, method = "glmnet", trControl = cctrl1, metric = "ROC", tuneLength = 10, preProc = c("center", "scale"), importance = "permutation", verbose = T)

BlMo.glmnet.model.HMs.cv10<-BlMo.glmnet.model.down.cv10
save(BlMo.glmnet.model.HMs.cv10, file = "BlMo.glmnet.model.HMs.cv10.RData")

cctrl1 <- trainControl(method="cv", number=10, returnResamp="all", classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(849)
test_class_cv_model <- train(BlMo.data.down_trainX, y,
                               method = "cforest",
                               trControl = cctrl1,
                               metric = "ROC",
                               preProc = c("center", "scale"),
                               controls = cforest_unbiased(ntree = 50, mtry = sqrt(ncol(BlMo.data.down_trainX))))

BlMo.cforest.model.HMs.cv10<-test_class_cv_model
save(BlMo.cforest.model.HMs.cv10, file = "BlMo.cforest.model.HMs.cv10.RData")

### same is repeated fro the remaining biological replicates
