clear all
%models = [1, 4]; benchmark = 0.6174; ymaxatleast = 0.65;
models = [2, 5, 3, 6]; benchmark = 0.9354; ymaxatleast = 0.95;
experiments = {'features_span_xy', 'features_xy', 'features_xyz', 'features_without_triangulation', 'features_all_but_deltas', 'features_all'};
legend_texts = {'XGB', 'XGB 4 classes', 'XGB, 2 classes', 'RF', 'RF, 4 classes', 'RF, 2 classes'};
S = [1, 10, 20, 40, 80, 160, 320, 640, 1280];
%close all
results = zeros(6,9,6);
for i = 1 : 6
    table = readtable(['acc_', experiments{i}, '.csv'], 'VariableNamingRule', 'preserve');
    results(i,:,:) = table2array(table(:,2:end));
end

for i = 1 : 9

    %f = figure;
    %f.OuterPosition = [1000, 100, 1000, 800];
    subplot(3, 3, i)
    hold off
    plot(1 : 6, squeeze(results(:, i, models)), 'x-' ,'LineWidth', 1.5)
    hold on
    plot([1,6], [benchmark, benchmark], '--k', 'LineWidth', 1.5)
    title(['{\it S} = ', num2str(S(i))])
    xlabel('Feature set')
    ylabel('Accuracy')
    %ylim([0.3, 1])
    set(gca, 'XTick', 1:6)
    set(gca, 'XTickLabel', {'1', '2', '3', '4', '5', '6'})
    if i == 1, legend(legend_texts{models}, 'location', 'best'), end
end
