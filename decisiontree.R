table <- read.csv("all_panc_lookup_scores_2.tsv_uniprot_revel.tsv", sep = '\t', stringsAsFactors = FALSE)
table <- table[table$Chromosome != "error searching variant", ]
clinvar <- gsub(table$Clinvar2.accessions, pattern = "\\[|\\'|\\]", replacement = "")
clinvar <- strsplit(clinvar, ", ")
a <- lapply(clinvar, function(x){
  if(length(x) > 1){
    x[1]
  }
  else{
    x
  }
})
table$Clinvar2.accessions <- unlist(a)

for(col in numerictests[numerictests != "revel"]){
  cur <- gsub(table[[col]], pattern = "\\[|\\'|\\]", replacement = "")
  cur <- strsplit(cur, ", ")
  a <- lapply(cur, function(x){
    if(length(x) > 1){
      x[1]
    }
    else{
      x
    }
  })
  table[[col]] <- unlist(a)
  table[[col]] <- as.numeric(as.character(table[[col]]))
}
table$revel <- as.numeric(as.character(table$revel))

numerictests <- c("revel","Phylop.Prediction", "Dann.score", "MutationTaster.Score","Sift.Score","MutationAssessor.Score","Fathmm.mkl.Score","Fathmm.Score","Metasvm.Score","Metalr.Score","Provean.Score", "lrt.Score")
t <- table[,c("marker","result","interpretation","Clinvar2.accessions",numerictests)]
benpath <- t[t$Clinvar2.accessions %in% c("Benign","Pathogenic","Likely benign","Likely pathogenic"),]
benpath$Clinvar2.accessions[benpath$Clinvar2.accessions == "Likely benign"] <- "Benign"
benpath$Clinvar2.accessions[benpath$Clinvar2.accessions == "Likely pathogenic"] <- "Pathogenic"

library(rpart)
treeaccuracy <- c()
for(n in 1:100){
  trainrows <- sample(nrow(benpath), as.integer(.7*nrow(benpath)))
  trainbenpath <- benpath[trainrows, ]
  testbenpath <- benpath[-trainrows, ]
  rpartfit <- rpart(Clinvar2.accessions ~ ., data = trainbenpath[4:ncol(trainbenpath)], method = 'class')
  predicted <- predict(rpartfit, testbenpath[5:ncol(trainbenpath)])
  testbenpath$predicted <- predicted[,1]
  testbenpath$predicted <- ifelse(testbenpath$predicted >= .5, "Benign","Pathogenic")
  treeaccuracy <- c(treeaccuracy,sum(testbenpath$Clinvar2.accessions == testbenpath$predicted) / nrow(testbenpath))
}
mean(treeaccuracy)
sd(treeaccuracy)
library(randomForest)
forestaccuracy <- c()
for(n in 1:100){
  trainrows <- sample(nrow(benpath), as.integer(.7*nrow(benpath)))
  trainbenpath <- benpath[trainrows, ]
  testbenpath <- benpath[-trainrows, ]
  trainbenpath$Clinvar2.accessions <- as.factor(trainbenpath$Clinvar2.accessions)
  forestfit <- randomForest(Clinvar2.accessions ~ ., trainbenpath[4:ncol(trainbenpath)], ntree = 500, na.action=na.exclude)
  predicted <- predict(forestfit, testbenpath[5:ncol(trainbenpath)])
  testbenpath$predicted <- predicted
  forestaccuracy <- c(forestaccuracy,sum(testbenpath$Clinvar2.accessions == testbenpath$predicted, na.rm = TRUE) / sum(!is.na(testbenpath$predicted)))
}
mean(forestaccuracy)
sd(forestaccuracy)
library(e1071)
svmaccuracy <- c()
for(n in 1:100){
  trainrows <- sample(nrow(benpath), as.integer(.7*nrow(benpath)))
  trainbenpath <- benpath[trainrows, ]
  testbenpath <- benpath[-trainrows, ]
  svmfit <- svm(Clinvar2.accessions ~ ., data = trainbenpath[4:ncol(trainbenpath)], type = "C")
  predicted <- predict(svmfit, testbenpath[5:ncol(trainbenpath)])
  testbenpath[complete.cases(testbenpath), "predicted"] <- predicted
  svmaccuracy <- c(svmaccuracy,sum(testbenpath$Clinvar2.accessions == testbenpath$predicted, na.rm = TRUE) / sum(!is.na(testbenpath$predicted)))
}
mean(svmaccuracy)
sd(svmaccuracy)

library(caret)
gbmaccuracy <- c()
for(n in 1:100){
  trainrows <- sample(nrow(benpath), as.integer(.7*nrow(benpath)))
  trainbenpath <- benpath[trainrows, ]
  testbenpath <- benpath[-trainrows, ]
  fitControl <- trainControl( method = "repeatedcv", number = 4, repeats = 4)
  fit <- train(Clinvar2.accessions ~ ., data = trainbenpath[4:ncol(trainbenpath)], method = "gbm", na.action  = na.pass, trControl = fitControl, verbose = FALSE)
  predicted <- predict(fit, testbenpath[5:ncol(trainbenpath)], type = "prob")
  testbenpath[complete.cases(testbenpath), "predicted"] <- predicted[,1]
  testbenpath$predicted <- ifelse(testbenpath$predicted >= .5, "Benign","Pathogenic")
  gbmaccuracy <- c(gbmaccuracy,sum(testbenpath$Clinvar2.accessions == testbenpath$predicted, na.rm = TRUE) / sum(!is.na(testbenpath$predicted)))
}
mean(gbmaccuracy)
sd(gbmaccuracy)
