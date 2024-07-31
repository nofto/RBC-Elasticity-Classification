import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import sys
from calculate_perm_imp import calculate_importances, model_descriptions

descriptions = model_descriptions()

starting_step = 100

experiments = {
    'features': list(range(6, 47)),
    #'features_all_but_deltas': list(range(6, 29)),
    'features_without_triangulation': [6, 7, 8, 20, 21, 22, 23, 24, 25, 26, 27, 28]
    #'features_xyz': [6, 7, 8, 23, 24, 25, 26, 27, 28],
    #'features_xy': [6, 7, 23, 24, 26, 27],
    #'features_span_xyz': [23, 24, 25],
    #'features_span_xy': [23, 24]
}

for experiment_desc, features in experiments.items():
    print(f'    Running experiment {experiment_desc}...')
    with open(f'permimp_{experiment_desc}.csv', 'w') as file:
        print(','.join(['frames'] + descriptions), file=file)
        for frames in [80]:
            print(f"        Frames = {frames}:")
            importances = calculate_importances(frames, features, starting_step)
            for index in range(len(importances[0])):
                print(index, file=file, end='')
                for model in range(4):
                    print(f',{importances[model][index]}', file=file, end='')
                print('', file=file)
        

