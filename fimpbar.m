feature_desc = ["rbc velocity x", "rbc velocity y", "rbc velocity z", "volume", "surface", "edge angle mean", "edge angle deviation", "edge angle skewness", "node angle mean", "node angle deviation", "node angle skewness", "edge length mean", "edge length deviation", "edge length skewness", "cell axis length", "equator diameter", "equator diameter min", "span x", "span y", "span z", "vel diff x", "vel diff y", "vel diff z", "edge angle delta mean", "edge angle delta deviation", "edge angle delta skewness", "node angle delta mean", "node angle delta deviation", "node angle delta skewness", "edge length delta mean", "edge length delta deviation", "edge length delta skewness", "edge angle delta abs mean", "edge angle delta abs deviation", "edge angle delta abs skewness", "node angle delta abs mean", "node angle delta abs deviation", "node angle delta abs skewness", "edge length delta abs mean", "edge length delta abs deviation", "edge length delta abs skewness", "rbc velocity x std", "rbc velocity y std", "rbc velocity z std", "volume std", "surface std", "edge angle mean std", "edge angle deviation std", "edge angle skewness std", "node angle mean std", "node angle deviation std", "node angle skewness std", "edge length mean std", "edge length deviation std", "edge length skewness std", "cell axis length std", "equator diameter std", "equator diameter min std", "span x std", "span y std", "span z std", "vel diff x std", "vel diff y std", "vel diff z std", "edge angle delta mean std", "edge angle delta deviation std", "edge angle delta skewness std", "node angle delta mean std", "node angle delta deviation std", "node angle delta skewness std", "edge length delta mean std", "edge length delta deviation std", "edge length delta skewness std", "edge angle delta abs mean std", "edge angle delta abs deviation std", "edge angle delta abs skewness std", "node angle delta abs mean std", "node angle delta abs deviation std", "node angle delta abs skewness std", "edge length delta abs mean std", "edge length delta abs deviation std", "edge length delta abs skewness std"];
model_desc = ["XGB, 4 classes", "XGB, 2 classes", "RF, 4 classes", "RF, 2 classes"];
fimp = table2array(readtable("permimp_features.csv", "VariableNamingRule", "preserve"));
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
print('predictor_importances_1','-depsc2');

