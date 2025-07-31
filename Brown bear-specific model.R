####Brown bear-specific models 250512####

#最初に以下のようにデータの読み込み・出力を行うPC上のフォルダを指定してください
#ディレクトリを指定
setwd("XXX")
#指定されたかどうか確認
getwd()


MD <- read.csv("measurement_data.csv")
MD[,2] <- (MD[,2] - 44.0627)/ 14.6684
MD[,3] <- (MD[,3] - 39.8913)/ 17.4028
MD[,4] <- (MD[,4] - 29.5218)/ 14.0691
MD[,5] <- (MD[,5] - 35.6504)/ 17.1062
write.csv(MD, "measurement_data_standardized.csv", row.names = FALSE)
MDS <- read.csv("measurement_data_standardized.csv")

BBBN<-read.csv("brown_bear_blood_standardized.csv")
MEAN <- 12.5503
SD <- 10.1130

# measurement_data.csvに測定したメチル化レベルを入力してください。
# 標準化（（元の測定値-トレーニングデータの平均値）÷トレーニングデータの標準偏差）もRで行います。
# 標準化されたものがmeasurement_data_standardized.csvとして出力されます。
 
# 有効数字を考え、平均値と標準偏差は5桁にしています。トレーニングデータでの元の値は以下の通りです。
# 年齢	SLC12A5-1	SLC12A5-2	SLC12A5-3	SLC12A5-4
# 平均値	12.55030455	44.06265306	39.89132653	29.52183673	35.65040816
# 標準偏差	10.11299266	14.66844243	17.4027674	14.06912513	17.10624161


# Support vector regression (SLC12A5-1, -2, -4)
library(e1071)
SVRMsl<-
  svm(age~SLC12A5_1_methylation_rate_ave+SLC12A5_2_methylation_rate_ave+SLC12A5_3_methylation_rate_ave+SLC12A5_4_methylation_rate_ave,data=BBBN,
      cost=12589.25   ,gamma=0.001995262   ,epsilon=0.1)
predicted_age_SVRM<-predict(SVRMsl,MDS)*SD+MEAN
predicted_age_SVRM[predicted_age_SVRM < 0] <- 0

#sample IDと推定年齢をR上で示す
data.frame(Sample_ID = MDS$Sample_ID, Predicted_Age = predicted_age_SVRM)

# Output to a csv file
sample_ID <- MDS[, 1, drop = FALSE]
predicted_age<-cbind(sample_ID,predicted_age_SVRM)
write.csv(predicted_age, "predicted_age_result.csv", row.names = FALSE)

# sample IDとSVRモデルの推定年齢の一覧がcsvファイルとして出力されます。

# この時点で、measurement_data_standardized.csv、predicted_age_result.csvの2つのファイルが出力されているはずです。
# これらをそのまま残しておくと次回以降エラーが出るので、ファイルの名前を変えるか別のフォルダに移動させてください。
