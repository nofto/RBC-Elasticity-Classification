import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
import time

def model_descriptions():
    return [
        'XGB 4',
        'XGB 2',
        'RF 4',
        'RF 2',
    ]

def calculate_importances(frames, feature_indices_to_be_used, starting_step):

    start_time = time.time()
    #print("  Loading train data... ", end="")
    train_data = pd.read_csv(f'output/features_train_frames_{frames}.csv')
    end_time = time.time()
    #print(end_time - start_time, "seconds")
    
    train_data = train_data[train_data['step'] >= starting_step]


    start_time = time.time()
    #print("  Loading validation data... ", end="")
    validation_data = pd.read_csv(f'output/features_test_frames_{frames}.csv')
    end_time = time.time()
    #print(end_time - start_time, "seconds")
    
    validation_data = validation_data[validation_data['step'] >= starting_step]


    label_column = 'cell_type'

    if frames > 1:
        feature_indices_to_be_used = feature_indices_to_be_used + [x + 44 for x in feature_indices_to_be_used];
    
    feature_columns = [train_data.columns[i] for i in feature_indices_to_be_used]

    # Extract features and labels from the training data
    X_train = train_data[feature_columns]
    y_train = train_data[label_column]
    y_train_2class = [0 if label == 0 else 1 for label in y_train]

    # Extract features and labels from the validation data
    X_validation = validation_data[feature_columns]
    y_validation = validation_data[label_column]
    y_validation_2class = [0 if label == 0 else 1 for label in y_validation]


    model_xgb_4class = xgb.XGBClassifier()
    model_xgb_2class = xgb.XGBClassifier()
    model_rf_4class = RandomForestClassifier()
    model_rf_2class = RandomForestClassifier()

    start_time = time.time()
    #print("  Training 4-class XGBClassifier... ", end="")
    model_xgb_4class.fit(X_train, y_train)
    end_time = time.time()
    #print(end_time - start_time, "seconds")
    
    start_time = time.time()
    #print("  Training 2-class XGBClassifier... ", end="")
    model_xgb_2class.fit(X_train, y_train_2class)
    end_time = time.time()
    #print(end_time - start_time, "seconds")
    
    start_time = time.time()
    #print("  Training 4-class RandomForestClassifier... ", end="")
    model_rf_4class.fit(X_train, y_train)
    end_time = time.time()
    #print(end_time - start_time, "seconds")
    
    start_time = time.time()
    #print("  Training 2-class RandomForestClassifier... ", end="")
    model_rf_2class.fit(X_train, y_train_2class)
    end_time = time.time()
    #print(end_time - start_time, "seconds")

    importances = []
    for index, model in enumerate([model_xgb_4class, model_xgb_2class, model_rf_4class, model_rf_2class]):
        importances.append(model.feature_importances_)

    return importances


    accuracies = []
    start_time = time.time()
    #print("  Predicting... ", end="")
    for index, model in enumerate([model_xgb_4class, model_xgb_2class, model_rf_4class, model_rf_2class]):
        y_prediction = model.predict(X_validation)
        if index % 2 == 0:
            accuracy4 = accuracy_score(y_validation, y_prediction)
            y_prediction_2class = [0 if label == 0 else 1 for label in y_prediction]
            accuracy2 = accuracy_score(y_validation_2class, y_prediction_2class)
            accuracies.append(accuracy4)
            accuracies.append(accuracy2)
        else:
            accuracy = accuracy_score(y_validation_2class, y_prediction)
            accuracies.append(accuracy)
    end_time = time.time()
    #print(end_time - start_time, "seconds")

    return accuracies            
