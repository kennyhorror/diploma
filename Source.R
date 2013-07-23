library('R.matlab');
library('kernlab');
library('fastICA');
library('RSNNS');

RemoveMean <- function(data) {
  avg = mean(data);
  return (data - avg);
}

Rereference <- function(data) {
  avg = c();
  for (i in 1:dim(data)[1]) {
    avg = rbind(avg, mean(data[i,]));
  }
  avg = rep(avg, dim(data)[2]);
  return (data - avg);
}

CalculateARFeatures <- function(data) {
  return(t(arima(data, c(3, 0, 0), include.mean = FALSE)$coef));
}

CalculateFARFeatures <- function(data) {
  return(cbind(CalculateARFeatures(data), CalculateFreqFeatures(data)));
}

CalculateFreqFeatures <- function(data, Fs = 250, medicalSeparation = TRUE) {
  lims = c(1:31); # Freqz
  if (medicalSeparation) {
    lims = c(1, 4, 7, 9, 10.5, 13, 20, 30);
  }
  result = rep(0, length(lims) - 1);
  sp = fft(data);
  current = 1;
  total = 0;
  for (f in 1:length(sp)) {
    if ((f - 1) * Fs / length(data)  > lims[current]) {
      if (current > 1) {
        total = 1;
        result[current - 1] = result[current - 1] / total;
        total = 0;
      }
      current = current + 1;
    }
    if (current == length(lims) + 1) {
      break();
    }
    if (current < 2) {
      next();
    }
    result[current - 1] = result[current - 1] + abs(sp[f]);
    total = total + 1;
  }
  total = sum(result);
  result = result / total;
  return(t(result));
}

ContainsBlink <- function(eyes_channel) {
  if (max(abs(runmed(diff(runmed(eyes_channel, 11)), 11))) > 4) {
  }
  m = length(eyes_channel);
  return(max(abs(runmed(diff(runmed(eyes_channel, 11)), 11))[11:(m - 11)]) > 4);
}

PreprocessData <- function(data, label, preprocesser) {
  items = dim(data)[1];
  # Eyes channel is the last one.
  channels = dim(data)[2] - 1;
  result = c();
  for (j in 1:19) {
    # Split data on frames. 1 second each. 0.5 second overlap.
    segment = (125 * (j-1) + 1):(125 * (j + 1));
    eyes_channel = RemoveMean(data[segment, 7]);
    if (ContainsBlink(eyes_channel)) {
      next();
    }
    current_segment = Rereference(data[segment, 1:channels]);    
    tryCatch({
      current = c();
      for (i in 1:channels) {
        current = cbind(current, preprocesser(current_segment[,i]));
      }
      current = cbind(current, label);
      result = rbind(result, current);
    }, error = function(err) {
      print("Bad segment");
    })
  }
  return(result);
}

PreprocessRecords <- function(dataset, subject) {
  test_names = c('baseline', 'multiplication', 'letter-composing', 'rotation', 'counting');
  dudes = list();
  data = list();
  q = 1;
  for (i in 1:length(dataset$data)) {
    record_data = c();
    if (is.na(dataset$data[[i]][[4]][1,1])) {
      next();
    }
    if (dataset$data[[i]][1] != subject) {
      next();
    }
    current_data = PreprocessData(t(dataset$data[[i]][[4]][,]), which(dataset$data[[i]][[2]][,] == test_names), 
                                  preprocesser = CalculateFARFeatures);
    record_data = rbind(record_data, current_data);    
    if (!is.null(record_data)) {
      record_data = as.data.frame(record_data);
      record_data[,dim(record_data)[2]] = factor(record_data[,dim(record_data)[2]],
                                                 levels = c(1, 2, 3, 4, 5, '?'), ordered = TRUE);
      colnames(record_data) = c(paste("F", 1:(dim(record_data)[2] - 1), sep = ""), 'action');
    }
    data[[q]] = record_data;
    q = q + 1;
  }
  return(data);
}

SubsetOfRecords <- function(records, exclude = c()) {
  result = list();
  id = 1;
  for (i in 1:length(records)) {
    if (i %in% exclude) {
      next();
    }
    result[[id]] = records[[i]];
    id = id + 1;
  }
  return(result);
}

MergeRecords <- function(records, exclude = c()) {
  result = c();
  for (i in 1:length(records)) {
    if (i %in% exclude) {
      next();
    }
    result = rbind(result, records[[i]])
  }
  return(result);
}

BuildNeuralNetworkClassifier <- function(train_set) {
  m = dim(train_set)[2];
  train_Y = matrix(FALSE, dim(train_set)[1], 5)
  for (i in 1:5) {
    train_Y[, i] = (train_set[, m] == i);
  }
  model = mlp(train_set[, 1:(m-1)], size=c(15, 15), train_Y);
  return(model);
}

BuildSVMClassifier <- function(train_set) {
  model = ksvm(action~., train_set, type = 'spoc-svc');
  return(model);
}

ClassifyUsingNeuralNetwork <- function(test_set, model) {
  m = dim(test_set)[2];
  prediction = max.col(predict(model, test_set[, 1:(m - 1)], type = "class"));
  return(prediction);
}

ClassifyUsingSVMClassifier <- function(test_set, model) {
  prediction = predict(model, test_set);
  return(prediction);
}

BuildClassifier <- function(records, Builder, Runner, components = 15) {
  train_set = MergeRecords(records);
  m = dim(train_set)[2];
  labels = train_set[, m];
  train_set = train_set[,1:(m - 1)];
  ica = fastICA(train_set, components, verbose=FALSE, maxit = 100);
  train_set = as.data.frame(cbind(as.matrix(train_set) %*% ica$K, labels));
  colnames(train_set) = c(paste("F", 1:(dim(train_set)[2] - 1), sep = ""), 'action');
  result = list();
  result$ica = ica$K;
  result$runner = Runner;
  result$model = Builder(train_set);
  return(result);
}

ClassifyRecord <- function(test_set, model) {
  if (is.null(test_set)) {
    print("No valid semgents. Failed to classify.");
    return(-1);
  }
  m = dim(test_set)[2];
  test_set = cbind(as.matrix(test_set[, 1:(m - 1)]) %*% model$ica, test_set[,m]);
  test_set = as.data.frame(test_set);
  test_set[, dim(test_set)[2]] = factor('?', levels = c(1, 2, 3, 4, 5, '?'));
  m = dim(test_set)[2];
  colnames(test_set) = c(paste("F", 1:(m - 1), sep = ""), 'action');
  predictions = model$runner(test_set, model$model);
  return(names(which.max(table(predictions))));
}

CVOnRecords <- function(records, builder, runner, components = 11, parts = 25) {
  errors = 0;
  order = sample(1:length(records));
  parts = min(parts, length(records));
  for (i in 0:(parts - 1)) {
    positions_test = order[which((1:length(records)) %% parts == i)];
    train_set = SubsetOfRecords(records, positions_test);
    model = BuildClassifier(train_set, builder, runner, components);
    for (j in 1:length(positions_test)) {
      test_set = records[[positions_test[j]]];
      if (is.null(test_set)) {
        errors = errors + 1;
        print("No valid semgents. Failed to classify.");
        next();
      }
      correct_answer = test_set[1, dim(test_set)[2]];
      my_answer = ClassifyRecord(test_set, model);
      if (my_answer != correct_answer) {
        errors = errors + 1;
      }
    }
  }
  return(errors/length(records));
}


dataset = readMat('eegdata.mat');

results_svm = c();
for (i in 1:7) {
  train_set = PreprocessRecords(dataset, paste("subject", i));
  errors = CVOnRecords(train_set, BuildSVMClassifier, ClassifyUsingSVMClassifier, 11);
  results_svm = cbind(results_svm, errors);
  print(paste("Completed person", i, "with error rate", errors));
}
print("Results with SVM:")
print(results_svm);

results_nn = c();
for (i in 1:7) {
  train_set = PreprocessRecords(dataset, paste("subject", i));
  errors = CVOnRecords(train_set, BuildNeuralNetworkClassifier, ClassifyUsingNeuralNetwork, 15);
  results_nn = cbind(results_nn, errors);
  print(paste("Completed person", i, "with error rate", errors));
}
print("Results with Neural Network:")
print(results_nn);