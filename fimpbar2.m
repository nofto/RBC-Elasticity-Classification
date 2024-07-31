feature_desc = ["rbc velocity x", "rbc velocity y", "rbc velocity z", "cell axis length", "equator diameter", "equator diameter min", "span x", "span y", "span z", "vel diff x", "vel diff y", "vel diff z", "rbc velocity x std", "rbc velocity y std", "rbc velocity z std", "cell axis length std", "equator diameter std", "equator diameter min std", "span x std", "span y std", "span z std", "vel diff x std", "vel diff y std", "vel diff z std"];
model_desc = ["XGB, 4 classes", "XGB, 2 classes", "RF, 4 classes", "RF, 2 classes"];
fimp = table2array(readtable("permimp_features_without_triangulation.csv", "VariableNamingRule", "preserve"));
fimp(:,1) = fimp(:,1) + 1;
for i = 1 : 4,
    fimp = sortrows(fimp, i + 1, "descend");
    subplot(2, 2, i)
    top = fimp(10:-1:1, 1);
    barh(feature_desc(top), fimp(10:-1:1,i + 1))
    title(model_desc(i))
    set(gca,'FontSize',7,'FontName','Times');
end
set(gcf,'units','centimeters','position',[1, 1, 13.9, 10]);
print('predictor_importances_2','-depsc2');

