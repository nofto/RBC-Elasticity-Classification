function acc = test_selected(featureIndicesToBeUsed, skipStartRatio)

frames = [1, 10*2.^(0:7)];
acc = zeros(length(frames), 2);

for i = 1 : length(frames)
    fr = frames(i);
    featureIndices = featureIndicesToBeUsed;
    if fr > 1
        featureIndices = [featureIndicesToBeUsed, featureIndicesToBeUsed + 44];
    end
    acc(i, :) = test(train(fr, featureIndices, skipStartRatio), fr, featureIndices, skipStartRatio);
end
