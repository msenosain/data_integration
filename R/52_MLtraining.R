
# ftimp_rf$importance %>% 
#   dplyr::arrange(desc(.)) %>%
#   dplyr::top_n(25) %>%
#   mutate(Features = rownames(.)) %>%
#   ggplot(., aes(x=Features, y=Overall)) +
#      scale_x_discrete(limits=rev(x$Features)) +
#      geom_bar(stat="identity") +
#      coord_flip()


TrainModel <- function(data, trainIndx, alg = c('RF', 'XGB'), 
    class_col = class_col, seed = 40, allowParallel = TRUE, 
    workers = 4, save_model = T, label = label){
    library(caret)

    # Set seed
    set.seed(seed)    
    alg <- match.arg(alg, c('RF', 'XGB'))

    # Data partition
    TrainSet <- data[ trainIndx, ]
    TestSet  <- data[-trainIndx, ]

    if(alg == 'RF'){
        message('Running Random Forest')
        # Initiating parallelization
        if(allowParallel){
        cluster <- parallel::makeCluster(workers, setup_strategy = "sequential") 
        doParallel::registerDoParallel(cluster)
        }
        # Computing training control parameters
        fitControl <- trainControl(
            ## Repeated 5–fold CV 
            method = "repeatedcv",
            number = 5,
            ## repeated 10 times
            repeats = 10,
            verboseIter = TRUE,
            allowParallel =TRUE,
            returnResamp = "all")

        # Fitting the model
        model_rf <- train(TrainSet[,-1], TrainSet[,1],
                        method = 'ranger',
                        # should be set high at least p/3
                        tuneLength = 10, 
                        trControl = fitControl,
                        ## parameters passed onto the ranger function
                        # the bigger the better.
                        num.trees = 700,
                        importance = "permutation")

        message('Model training done!')

        # Assesing feature importance
        ftimp_rf <- caret::varImp(model_rf)

        # Test
        # Prediction on test dataset
        pred_rf <- predict(model_rf, TestSet)
        rmse_test <- RMSE(pred_rf, TestSet[[class_col]])

        message('Testing done!')

        # Finalizing parallelization
        if(allowParallel){
        parallel::stopCluster(cluster)
        foreach::registerDoSEQ()
        }
        message('Random Forest completed')
        
        if(save_model){
            #save RData
            save(model_rf, ftimp_rf, pred_rf, rmse_test, TrainSet, TestSet,
                file = paste0(label, '_RFmodel.RData'))            
        } else {
            x = list(model_rf=model_rf, ftimp_rf=ftimp_rf, pred_rf=pred_rf, 
                rmse_test=rmse_test, TrainSet=TrainSet, TestSet=TestSet)
            return(x)
        }
    }
    

    if(alg == 'XGB'){
        print('Running XGBoost')
        # Rearrangement
        X_train = xgboost::xgb.DMatrix(as.matrix(TrainSet[,-1]))
        y_train = TrainSet[,1]

        # Computing training control parameters
        xgb_trcontrol = caret::trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 3,
            search = 'grid',
            savePredictions = 'final',
            summaryFunction = twoClassSummary,
            allowParallel = allowParallel,
            classProbs = TRUE,
            verboseIter = FALSE,
            returnData = FALSE)

        # Specifying grid space
        xgbGrid <- expand.grid(nrounds = c(100,200),  # this is n_estimators in the python code above
                               max_depth = c(10, 15, 20, 25),
                               colsample_bytree = seq(0.5, 0.9, length.out = 5),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.1,
                               gamma=0,
                               min_child_weight = 1,
                               subsample = 1)

        # Fitting the model
        model_xgb = caret::train(
            X_train, y_train,
            trControl = xgb_trcontrol,
            tuneGrid = xgbGrid,
            method = "xgbTree",
            metric = 'ROC',
            nthread = workers)

        # Assesing feature importance
        ftimp_xgb <- caret::varImp(model_xgb)
        
        # Test
        X_test = xgboost::xgb.DMatrix(as.matrix(TestSet[,-1]))
        y_test = TestSet[,1]

        # Prediction on test dataset
        pred_xgb <- predict(model_xgb, X_test)
        # Compute confusion matrix
        conf_xgb <- caret::confusionMatrix(
            reference = as.factor(y_test),
            data = pred_xgb,
            mode = 'everything')

        # Save RData
        save(model_xgb, ftimp_xgb, pred_xgb, conf_xgb, 
            TrainSet, TestSet, file = paste0(label, '_XGBmodel.RData'))
        print('XGBoost completed')
    }

}

ftimp_gg <- function(model_ls, n_topft = 10){
    model_ls$ftimp_rf$importance %>% 
    dplyr::arrange(desc(.)) %>%
    dplyr::top_n(n_topft) %>%
    mutate(Features = rownames(.)) %>%
    ggplot(., aes(x=reorder(Features, Overall), y=Overall)) +
        geom_bar(stat="identity", fill='steelblue') +
        labs(x ="Features", y = "Importance")+
        coord_flip()
}

frac_dif <- function(model_ls){
    frac_dif <- ((model_ls$pred_rf/model_ls$TestSet$SILA_S)-1)*100
}

