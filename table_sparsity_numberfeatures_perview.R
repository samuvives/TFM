library(dplyr)
library(gt)
library(webshot2)
path = "/home/vant/Escritorio/TFM/outputforR/tableuniquevaluessparsity.csv"
data <- read.csv(path)
colnames(data)[2] <- "Number of features"
colnames(data)[3] <- "Sparsity (%)"
data <- data %>%
    arrange(desc(`Number of features`))
head(data)
gttable <- gt(data) %>%
    tab_header(title = "Information of the views")
outputpath <- "/home/vant/Escritorio/TFM/defensaTFM/tableuniquevaluessparsity.png"
gtsave(gttable, outputpath)
