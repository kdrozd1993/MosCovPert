validationPredictions = trainedModel.predictFcn(classifierTableSig);
validationResponse = StartMosCovLabelX;
c = confusionmat(validationResponse',validationPredictions')
plotConfMat(c)
