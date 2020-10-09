library(randomForest)
library(caret)
library(caTools)

# Reading data
data <- read.csv('/Users/senosam/Documents/Massion_lab/data_integration/ML/m_cytof.csv')
sample = sample.split(data$SILA_S, SplitRatio = .75)
train = subset(data, sample == TRUE)
test  = subset(data, sample == FALSE)


# ML
workers = 10
cluster <- parallel::makeCluster(workers) 
doParallel::registerDoParallel(cluster)




fitControl <- trainControl(
  ## Repeated 5–fold CV 
  method = "repeatedcv",
  number = 5,
  ## repeated 10 times
  repeats = 10,
  verboseIter = TRUE,
  allowParallel =TRUE,
  returnResamp = "all")


rrfFit <- train(SILA_S ~ ., 
                 data = train,
                 method = 'ranger',
                # should be set high at least p/3
                 tuneLength = 10, 
                 trControl = fitControl,
                ## parameters passed onto the ranger function
                # the bigger the better.
                 num.trees = 700,
                 importance = "permutation")

rrfFit
trellis.par.set(caretTheme())
plot(rrfFit)

# Assesing feature importance
ftimp_rf <- caret::varImp(rrfFit)
ftimp_rf$importance %>% 
  dplyr::arrange(desc(.)) %>%
  dplyr::top_n(25) %>%
  mutate(Features = rownames(.)) %>%
  ggplot(., aes(x=Features, y=Overall)) +
     scale_x_discrete(limits=rev(x$Features)) +
     geom_bar(stat="identity") +
     coord_flip()

# Model performance
rf_pred <- predict(rrfFit, test)
RMSE(rf_pred, test$SILA_S)

