% frames = [1, 10*2.^(0:7)];
% acc_only_means_span_stdevs = zeros(length(frames), 2);
% 
% for i = 1 : length(frames)
%     fr = frames(i);
%     featureIndices = 24 : 26;
%     if fr > 1
%         featureIndices = [24 : 26, (24 : 26) + 26];
%     end
%     acc_only_means_span_stdevs(i, :) = test(train(fr, featureIndices), fr, featureIndices);
% end

acc_only_means_span_stdevs = [
    0.755230328624824	0.378416096306005
    0.813852339181287	0.457175925925926
    0.823923152227117	0.462946647087616
    0.853074074074074	0.480617283950617
    0.870751633986928	0.516251885369533
    0.876656233698487	0.538119457485655
    0.879075014100395	0.518739424703892
    0.879680134680135	0.487508417508418
    0.880720572057206	0.487156215621562
];