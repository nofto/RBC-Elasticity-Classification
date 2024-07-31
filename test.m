function accuracies = test(rf, frames, featureIndices, skipStartRatio)
file = ['output\features_test_frames_', num2str(frames), '.csv'];
[features, response] = readfile(file, featureIndices, skipStartRatio);

%disp('Testing on random forest ...')
tic
prediction = str2double(rf.predict(features));
elapsedTime = toc;
%disp(['                         ... done in ', num2str(elapsedTime), ' sec.'])
disp(['Testing  ', num2str(frames), ': ', num2str(elapsedTime)])
accuracy4classes = sum(prediction == response) / length(response);
accuracy2classes = sum(~xor(prediction, response)) / length(response);
accuracies = [accuracy2classes, accuracy4classes];