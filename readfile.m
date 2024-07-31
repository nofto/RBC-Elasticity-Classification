function [features, response] = readfile(file, featureIndices, skipStartRatio)

%disp(['Reading data from ', file, '...'])
tic
data = table2array(readtable(file));
[samples, columns] = size(data);
elapsedTime = toc;
%disp(['                         ... done in ', num2str(elapsedTime), ' sec. --> ', num2str(samples), ' samples × ', num2str(columns), ' columns'])

if any(featureIndices > columns | featureIndices < 4)
    error('Feature indices out of range!')
end

firstIndex = 1 + round(skipStartRatio*samples);

features = data(firstIndex : end, featureIndices);
response = data(firstIndex : end, 3);