% frames = [1, 10*2.^(0:7)];
% acc_only_means_noabs_also_stdevs = zeros(length(frames), 2);
% 
% for i = 1 : length(frames)
%     fr = frames(i);
%     featureIndices = 7 : 29;
%     if fr > 1
%         featureIndices = [7 : 29, 33 : 55];
%     end
%     acc_only_means_noabs_also_stdevs(i, :) = test(train(fr, featureIndices), fr, featureIndices);
% end

%results:
acc_only_means_noabs_also_stdevs = [
    0.969795155574972	0.819741760108733
    0.972197855750487	0.828715886939571
    0.974461576113559	0.839243759177680
    0.980506172839506	0.853111111111111
    0.982491201608849	0.867621920563097
    0.987180490349504	0.895800730307773
    0.985941906373378	0.897236322617033
    0.984225589225589	0.901801346801347
    0.995049504950495	0.892821782178218
];