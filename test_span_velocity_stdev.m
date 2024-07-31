frames = [1, 10*2.^(0:7)];
acc_also_devs_span_velocity = zeros(length(frames), 2);

for i = 1 : length(frames)
    fr = frames(i);
    featureIndices = [7:9,24:26];
    if fr > 1
        featureIndices = [7:9,24:26, ([7:9,24:26]) + 26];
    end
    acc_also_devs_span_velocity(i, :) = test(train(fr, featureIndices), fr, featureIndices);
end

acc_also_devs_span_velocity = [
    0.783178001067909	0.408366098733071
    0.867458576998051	0.521247563352827
    0.869273127753304	0.518502202643172
    0.882283950617284	0.523617283950617
    0.890786827551534	0.546983408748115
    0.892044861763172	0.570970266040689
    0.885758601240835	0.589142695995488
    0.889444444444445	0.590639730639731
    0.891584158415842	0.645792079207921
];