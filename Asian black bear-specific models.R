#●Asian bear-specific models 250418

#最初に以下のようにデータの読み込み・出力を行うPC上のフォルダを指定してください
#ディレクトリを指定
setwd("XXX")
#指定されたかどうか確認
getwd()


MD <- read.csv("measurement_data.csv")
MD[,2] <- (MD[,2] - 36.585) / 12.342
MD[,3] <- (MD[,3] - 32.920) / 12.203
MD[,4] <- (MD[,4] - 24.984) / 9.9268
MD[,5] <- (MD[,5] - 31.146) / 11.636
write.csv(MD, "measurement_data_standardized.csv", row.names = FALSE)
MDS <- read.csv("measurement_data_standardized.csv")

ABBBS <- read.csv("Asian_black_bear_blood_standardized.csv")
MEAN <- 9.6959
SD <- 6.7553

#measurement_data.csvに測定したメチル化レベルを入力してください。
#標準化（（元の測定値-トレーニングデータの平均値）÷トレーニングデータの標準偏差）もRで行います。
#標準化されたものがmeasurement_data_standardized.csvとして出力されます。

#有効数字を考え、平均値と標準偏差は5桁にしています。トレーニングデータでの元の値は以下の通りです。
#年齢	SLC12A5-1	SLC12A5-2	SLC12A5-3	SLC12A5-4
#平均値	9.695895186	36.58480769	32.91990385	24.98384615	31.14557692
#標準偏差	6.755301022	12.34236483	12.20324742	9.926796455	11.63555204

#●Single regression (SLC12A5-4)
SRM <-lm(formula=age~SLC12A5_4_methylation_rate_ave,data=ABBBS)
predicted_age_SRM <- predict(SRM,MDS)*SD+MEAN	
predicted_age_SRM[predicted_age_SRM < 0] <- 0

#sample IDと推定年齢をR上で示す
data.frame(Sample_ID = MDS$Sample_ID, Predicted_Age = predicted_age_SRM)

#●Principal component regression (PC1)
pca_slc <- prcomp(ABBBS[, c("SLC12A5_1_methylation_rate_ave", 
                            "SLC12A5_2_methylation_rate_ave", 
                            "SLC12A5_3_methylation_rate_ave", 
                            "SLC12A5_4_methylation_rate_ave")], 
                  scale = FALSE)
ABBBS$SLC12A5_PC1 <- pca_slc$x[,1]
PCRM <- lm(formula = age ~ SLC12A5_PC1, data = ABBBS)

#MDSの4CpGの値を1つに 
MDS$SLC12A5_PC1 <- as.numeric(as.matrix(MDS[,c("SLC12A5_1_methylation_rate_ave", "SLC12A5_2_methylation_rate_ave", "SLC12A5_3_methylation_rate_ave", "SLC12A5_4_methylation_rate_ave")]) %*% pca_slc$rotation[,1])

predicted_age_PCRM <- predict(PCRM,MDS)*SD+MEAN	
predicted_age_PCRM[predicted_age_PCRM < 0] <- 0

#sample IDと推定年齢をR上で示す
data.frame(Sample_ID = MDS$Sample_ID, Predicted_Age = predicted_age_PCRM)

#●Elastic net regression (SLC12A5-1, -2, -3, -4)
library(glmnet)
ENM <- glmnet(x = cbind(ABBBS$SLC12A5_1_methylation_rate_ave,ABBBS$SLC12A5_2_methylation_rate_ave,ABBBS$SLC12A5_3_methylation_rate_ave,ABBBS$SLC12A5_4_methylation_rate_ave),
              y = ABBBS$age, family = "gaussian", lambda = 0.009484228, alpha = 0.01, standardize = FALSE)
MDS_ENM <- cbind(MDS$SLC12A5_1_methylation_rate_ave, MDS$SLC12A5_2_methylation_rate_ave, MDS$SLC12A5_3_methylation_rate_ave, MDS$SLC12A5_4_methylation_rate_ave)
predicted_age_ENM <- c(predict(ENM,MDS_ENM,s= 0.009484228)*SD+MEAN)
predicted_age_ENM[predicted_age_ENM < 0] <- 0

#sample IDと推定年齢をR上で示す
data.frame(Sample_ID = MDS$Sample_ID, Predicted_Age = predicted_age_ENM)

#●Support vector regression (SLC12A5-1, -2, -4)
library(e1071)
SVRM<-
  svm(age~SLC12A5_1_methylation_rate_ave+SLC12A5_2_methylation_rate_ave+SLC12A5_4_methylation_rate_ave,data=ABBBS, cost=10^4.9, gamma=10^-3.1, epsilon=0.1, scale = FALSE)
predicted_age_SVRM<-predict(SVRM,MDS)*SD+MEAN
predicted_age_SVRM[predicted_age_SVRM < 0] <- 0

#sample IDと推定年齢をR上で示す
data.frame(Sample_ID = MDS$Sample_ID, Predicted_Age = predicted_age_SVRM)

#●Output to a csv file
sample_ID <- MDS[, 1, drop = FALSE]
predicted_age<-cbind(sample_ID,predicted_age_SRM,predicted_age_PCRM,predicted_age_ENM,predicted_age_SVRM)
write.csv(predicted_age, "predicted_age_result.csv", row.names = FALSE)

#sample IDと4種の推定年齢の一覧がcsvファイルとして出力されます。

#この時点で、measurement_data_standardized.csv、predicted_age_result.csvの2つのファイルが出力されているはずです。
#これらをそのまま残しておくと次回以降エラーが出るので、ファイルの名前を変えるか別のフォルダに移動させてください。
