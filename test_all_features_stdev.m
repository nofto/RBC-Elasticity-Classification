% frames = [1, 10*2.^(0:7)];
% acc = zeros(length(frames), 2);
% 
% for i = 1 : length(frames)
%     fr = frames(i);
%     featureIndices = 4:29;
%     if fr > 1
%         featureIndices = 4:55;
%     end
%     acc(i, :) = test(train(fr, featureIndices), fr, featureIndices);
% end
% acc_also_stdevs = acc;

%results:
acc_also_stdevs = [
    0.967173923595942	0.708011747002573
    0.970991715399610	0.758516081871345
    0.973580518844836	0.767988252569750
    0.978592592592593	0.816938271604938
    0.982780291603821	0.848227752639517
    0.987154407929056	0.881977047470005
    0.985265087422448	0.888959390862944
    0.985976430976431	0.898922558922559
    0.994774477447745	0.878960396039604
];