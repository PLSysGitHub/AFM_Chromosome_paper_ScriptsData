close all;
clear all;


num_gauss = 5;
rng(4); % seed
options = statset('MaxIter',1000);
% load, restructure and clean data
load('masked_height_updated.mat');
data = data(:);
%data = 1e9* data(data>0); %unit is nm


model = fitgmdist(data,num_gauss,'Options',options);

x_plot = [0:0.1:150]';
figure()
histogram(data,100,'Normalization','pdf')
hold on
%[b,a] = ksdensity(data);
%plot(a,b,'linewidth',2)
y = model.pdf(x_plot);
plot(x_plot,y,'linewidth',2)
for i=[1:1:num_gauss]
    y1 = gmdistribution(model.mu(i),model.Sigma(i)).pdf(x_plot);
    y1 = y1*model.ComponentProportion(i);
    plot(x_plot,y1,'linewidth',2)
end
