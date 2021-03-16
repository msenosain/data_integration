# Set up environment
library(caret)
library(dplyr)


path='/Users/senosam/Documents/Massion_lab/data_integration/models'

# read file names in data
file_ls <- list.files(file.path(path, 'data'))
# remove info files from list
file_ls <- file_ls[-grep('info', file_ls)]
# remove radiomics for now
file_ls <- file_ls[-c(3,5)]
# read train/test info file
traintest_info <- read.csv(file.path(path, 'traintest_info.csv'))

TrainModel <- function(data, traintest_info, alg, tuneLength_n=5, input_ready=FALSE, data_train, data_test){
    
    if(!input_ready){
    # Replace NA with 0
        data <- data %>%
                mutate_all(~replace(., is.na(.), 0))

        # Get train and test pt_ID
        pts_train <- as.character(traintest_info$pt_ID[which(traintest_info$ds == 'training')])
        pts_test <- as.character(traintest_info$pt_ID[which(traintest_info$ds == 'test')])

        data_train <- data[pts_train,]
        data_test <- data[pts_test,]

        data_train['SILA'] <- traintest_info$SILA[which(traintest_info$ds == 'training')]
        data_test['SILA'] <- traintest_info$SILA[which(traintest_info$ds == 'test')]
    }
    # Set seed
    set.seed(101)
    
    # Parallel runnning 
    # cluster <- parallel::makeCluster(8, setup_strategy = "sequential")
    # doParallel::registerDoParallel(cluster)

    if(alg=='PLS') {
        fitControl <- trainControl(
            ## Repeated 5–fold CV 
            method = "repeatedcv",
            number = 5,
            ## repeated 10 times
            repeats = 10,
            verboseIter = FALSE,
            allowParallel =TRUE,
            returnResamp = "all")

        model_fit <- train(
            SILA ~ .,
            data = data_train,
            method = 'pls',
            preProc = c("center", "scale"),
            tuneLength = tuneLength_n,
            trControl = fitControl
            )

        # Prediction on test dataset
        pred <- predict(model_fit, data_test)
        rmse_test <- RMSE(pred, data_test[['SILA']])
    }

    if(alg=='ENET') {
        fitControl <- trainControl(
            ## Repeated 5–fold CV 
            method = "repeatedcv",
            number = 5,
            ## repeated 10 times
            repeats = 10,
            verboseIter = FALSE,
            allowParallel =TRUE,
            returnResamp = "all")

        model_fit <- train(
            SILA ~ .,
            data = data_train,
            method = 'glmnet',
            #standardize=standardize,
            preProc = c("center", "scale"),
            tuneLength = tuneLength_n,
            trControl = fitControl
            )

        # Prediction on test dataset
        pred <- predict(model_fit, data_test)
        rmse_test <- RMSE(pred, data_test[['SILA']])
    }


    if(alg=='RF') {
        fitControl <- trainControl(
            ## Repeated 5–fold CV 
            method = "repeatedcv",
            number = 5,
            ## repeated 10 times
            repeats = 10,
            verboseIter = FALSE,
            allowParallel =TRUE,
            returnResamp = "all")

        model_fit <- train(
            SILA ~ .,
            data = data_train,
            method = 'ranger',
            preProc = c("center", "scale"),
            tuneLength = tuneLength_n,
            trControl = fitControl,
            num.trees = 700,
            importance = "permutation"
            )

        # Prediction on test dataset
        pred <- predict(model_fit, data_test)
        rmse_test <- RMSE(pred, data_test[['SILA']])
    }
    
    model_ls <- list(model_fit=model_fit, pred=pred, rmse_test=rmse_test,
            TrainSet=data_train, TestSet=data_test)
    # parallel::stopCluster(cluster)
    # foreach::registerDoSEQ()
    
    message("Model ", alg, " completed!")
    return(model_ls)

}

w_AllModels <- function(path_to_dataset, traintest_info){
    dt <- read.csv(path_to_dataset, row.names = 1)

    model_RF <- TrainModel(data=dt, traintest_info=traintest_info, alg='RF')
    model_PLS <- TrainModel(data=dt, traintest_info=traintest_info, alg='PLS')
    model_ENET <- TrainModel(data=dt, traintest_info=traintest_info, alg='ENET')

    models_ls <- list(model_RF=model_RF, model_PLS=model_PLS, model_ENET=model_ENET)

    message("Data set ", path_to_dataset, " completed!")

    return(models_ls)

}

models_by_ds <- list()
for (i in file_ls){
    models_by_ds[[gsub('.csv', '',i)]] <- w_AllModels(file.path(path, 'data', i), traintest_info)
}

save(models_by_ds, traintest_info, file_ls, file=file.path(path, 'models_res', 'models_by_ds.RData'))
# run models for concatetaned files

# concat

# Focus on:
# cytof (1 or both)
# rna reactome
# wes binary

file_ls <- file_ls[-4]

pts_train <- as.character(traintest_info$pt_ID[which(traintest_info$ds == 'training')])
pts_test <- as.character(traintest_info$pt_ID[which(traintest_info$ds == 'test')])

d_train <- data.frame('SILA'=traintest_info$SILA[which(traintest_info$ds == 'training')])
d_test <- data.frame('SILA'=traintest_info$SILA[which(traintest_info$ds == 'test')])

for (i in file_ls) {
    path_to_dataset = file.path(path, 'data', i)
    dt <- read.csv(path_to_dataset, row.names = 1)
    dt <- dt %>%
        mutate_all(~replace(., is.na(.), 0))

    data_train <- dt[pts_train,]
    data_test <- dt[pts_test,]
    colnames(data_train) <- paste0(sapply(strsplit(i, "_"), "[[", 1), '_', colnames(data_train))
    colnames(data_test) <- paste0(sapply(strsplit(i, "_"), "[[", 1), '_', colnames(data_test))

    d_train <- cbind(d_train, data_train)
    d_test <- cbind(d_test, data_test)
}


model_RF <- TrainModel(traintest_info=traintest_info, alg='RF', input_ready=TRUE, data_train=d_train, data_test=d_test)
model_PLS <- TrainModel(traintest_info=traintest_info, alg='PLS', input_ready=TRUE, data_train=d_train, data_test=d_test)
model_ENET <- TrainModel(traintest_info=traintest_info, alg='ENET', input_ready=TRUE, data_train=d_train, data_test=d_test)

models_CONCAT_ls <- list(model_RF=model_RF, model_PLS=model_PLS, model_ENET=model_ENET)

save(models_CONCAT_ls, traintest_info, file=file.path(path, 'models_res', 'models_CONCAT_ls.RData'))
















model_ls <- models_by_ds$cytof_freq$model_ENET
model_name = 'Elastic Net'
cat('RMSE test: ', model_ls$rmse_test, '\n')
model_ls$model_fit
plot(model_ls$model_fit)

coefficients <- as.data.frame(as.matrix(coef(model_ls$model_fit$finalModel, model_ls$model_fit$bestTune$lambda)))
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])


