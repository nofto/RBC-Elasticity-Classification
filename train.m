function rf = train(frames, featureIndices, skipStartRatio)
file = ['output\features_train_frames_', num2str(frames), '.csv'];
[features, response] = readfile(file, featureIndices, skipStartRatio);

%disp('Training random forest ...')
tic
trees = 100;
rf = TreeBagger(trees, features, response, 'Method', 'classification');
elapsedTime = toc;
%disp(['                         ... done in ', num2str(elapsedTime), ' sec.'])
disp(['Training ', num2str(frames), ': ', num2str(elapsedTime)])
