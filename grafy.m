clear all
%models = [1, 4]; benchmark = 0.6174; ymaxatleast = 0.65;
models = [2, 5, 3, 6]; benchmark = 0.9354; ymaxatleast = 0.95;
experiments = {'features_span_xy', 'features_xy', 'features_xyz', 'features_without_triangulation', 'features_all_but_deltas', 'features_all'};
legend_texts = {'XGB', 'XGB, 4 classes', 'XGB, 2 classes', 'RF', 'RF, 4 classes', 'RF, 2 classes'};
S = [1, 10, 20, 40, 80, 160, 320, 640, 1280];
%close all
for i = 1 : 6
    table = readtable(['acc_', experiments{i}, '.csv'], 'VariableNamingRule', 'preserve');
    results{i} = table2array(table);
    %f = figure;
    %f.OuterPosition = [1000, 100, 1000, 800];
    subplot(3, 2, i)
    hold off
    plot(1 : 9, results{i}(:, 1 + models), 'x-' ,'LineWidth', 1.5)
    hold on
    plot([1,9], [benchmark, benchmark], '--k', 'LineWidth', 1.5)
    yl = get(gca,'ylim');
    if yl(2) < ymaxatleast, ylim([yl(1), ymaxatleast]), end
    title(['Feature set ', num2str(i)])
    xlabel('S', 'FontAngle', 'italic')
    ylabel('Accuracy')
    %ylim([0.7, 1])
    set(gca, 'XTick', 1:9)
    set(gca, 'XTickLabel', {'1', '10', '20', '40', '80', '160', '320', '640', '1280'})
    if i == 1, legend(legend_texts{models}, 'location', 'best'), end
end
