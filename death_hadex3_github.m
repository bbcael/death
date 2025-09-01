%% clear workspace & load asian HadEX3 annual precip extreme metrics
clear all; close all; clc; 
r1 = ncread('HadEX3_Rx1day_ANN.nc','Rx1day'); r1 = r1(:,:,end-30:end);
r = r1; r(1:96,:,:) = r1(97:end,:,:); r(97:end,:,:) = r1(1:96,:,:); r1 = r; clear r;
a1 = r1(134:173,79:125,:); clear r1;
r5 = ncread('HadEX3_Rx5day_ANN.nc','Rx5day'); r5 = r5(:,:,end-30:end);
r = r5; r(1:96,:,:) = r5(97:end,:,:); r(97:end,:,:) = r5(1:96,:,:); r5 = r; clear r;
a5 = r5(134:173,79:125,:); clear r5;
r10 = ncread('HadEX3_R10mm_ANN.nc','R10mm'); r10 = r10(:,:,end-30:end);
r = r10; r(1:96,:,:) = r10(97:end,:,:); r(97:end,:,:) = r10(1:96,:,:); r10 = r; clear r;
a10 = r10(134:173,79:125,:); clear r10;
r20 = ncread('HadEX3_R20mm_ANN.nc','R20mm'); r20 = r20(:,:,end-30:end);
r = r20; r(1:96,:,:) = r20(97:end,:,:); r(97:end,:,:) = r20(1:96,:,:); r20 = r; clear r;
a20 = r20(134:173,79:125,:); clear r20;
r95 = ncread('HadEX3_61-90_R95p_ANN.nc','R95p'); r95 = r95(:,:,end-30:end);
r = r95; r(1:96,:,:) = r95(97:end,:,:); r(97:end,:,:) = r95(1:96,:,:); r95 = r; clear r;
a95 = r95(134:173,79:125,:); clear r95;
r99 = ncread('HadEX3_61-90_R99p_ANN.nc','R99p'); r99 = r99(:,:,end-30:end);
r = r99; r(1:96,:,:) = r99(97:end,:,:); r(97:end,:,:) = r99(1:96,:,:); r99 = r; clear r;
a99 = r99(134:173,79:125,:); clear r99;
%% look for linear trends in sum-based variables over each grid cell
for i = 1:40;
    for j = 1:47;
        y = squeeze(a10(i,j,:));
        if sum(isnan(y))<31;
            mdl = fitlm(1:31,y);
            t10(i,j) = table2array(mdl.Coefficients(2,3));
        end
    end
end
for i = 1:40;
    for j = 1:47;
        y = squeeze(a20(i,j,:));
        if sum(isnan(y))<31;
            mdl = fitlm(1:31,y);
            t20(i,j) = table2array(mdl.Coefficients(2,3));
        end
    end
end
for i = 1:40;
    for j = 1:47;
        y = squeeze(a95(i,j,:));
        if sum(isnan(y))<31;
            mdl = fitlm(1:31,y);
            t95(i,j) = table2array(mdl.Coefficients(2,3));
        end
    end
end
for i = 1:40;
    for j = 1:47;
        y = squeeze(a99(i,j,:));
        if sum(isnan(y))<31;
            mdl = fitlm(1:31,y);
            t99(i,j) = table2array(mdl.Coefficients(2,3));
        end
    end
end   
%% look for trends in daily/5-day max precip
a1(17,18,:) = NaN; a1(34,33,:) = NaN; a1(35,33,:) = NaN; a1(17,17,:) = NaN; a1(21,40,:) = NaN; a1(28,19,:) = NaN; % remove a few points that were causing convergence issues -- no idea why, but not enough to affect conclusions
for i = 1:40;
    for j = 1:47;
        y = squeeze(a1(i,j,:));
        if sum(isnan(y))==0;
            t = linspace(-1,1,31)';
            t = t(~isnan(y));
            phat = mle(y,'distribution','gev');
            [phat,pci] = mle(y,'pdf',@(x,a,b,c,d)gevpdf(x,a,abs(b),c+d.*t),'Start',[phat(1) phat(2) phat(3) 0]);
            t1(i,j) = 3.96.*phat(4)./(pci(2,4)-pci(1,4));
       end
    end
end   
a5(8,42,:) = NaN; a5(17,15,:) = NaN; a5(21,37,:) = NaN; a5(22,36,:) = NaN; a5(25,14,:) = NaN; a5(26,40,:) = NaN; a5(27,40,:) = NaN; a5(30,41,:) = NaN; a5(33,35,:) = NaN; a5(34,34,:) = NaN; a5(39,36,:) = NaN; a5(8,35,:) = NaN; % remove a few points that were causing convergence issues -- no idea why, but not enough to affect conclusions
for i = 1:40;
    for j = 1:47;
        y = squeeze(a5(i,j,:));
        if sum(isnan(y))==0;
            t = linspace(-1,1,31)';
            t = t(~isnan(y));
            phat = mle(y,'distribution','gev');
            [phat,pci] = mle(y,'pdf',@(x,a,b,c,d)gevpdf(x,a,abs(b),c+d.*t),'Start',[phat(1) phat(2) phat(3) 0]);
            t5(i,j) = 3.96.*phat(4)./(pci(2,4)-pci(1,4));
        end
    end
end   
clear t mdl i j pci phat y a*;
clc;
%% check percentages relative to null expectations
sum(sum(t10<-1.96))./sum(sum(t10~=0))./.025
sum(sum(t20<-1.96))./sum(sum(t20~=0))./.025
sum(sum(t95<-1.96))./sum(sum(t95~=0))./.025
sum(sum(t99<-1.96))./sum(sum(t99~=0))./.025
sum(sum(t1<-1.96))./sum(sum(t1~=0))./.025
sum(sum(t5<-1.96))./sum(sum(t5~=0))./.025
sum(sum(t10>1.96))./sum(sum(t10~=0))./.025
sum(sum(t20>1.96))./sum(sum(t20~=0))./.025
sum(sum(t95>1.96))./sum(sum(t95~=0))./.025
sum(sum(t99>1.96))./sum(sum(t99~=0))./.025
sum(sum(t1>1.96))./sum(sum(t1~=0))./.025
sum(sum(t5>1.96))./sum(sum(t5~=0))./.025
%% make figure
[y,x] = ecdf(t1(t1~=0));
p1 = plot(x,y,'linewidth',2,'color',[.7 .3 .1])
hold on;
p7 = plot(linspace(-3,3),normcdf(linspace(-3,3)),'--','color',[.7 .3 .1])
[y,x] = ecdf(t5(t5~=0));
p2 = plot(2+x,y,'linewidth',2,'color',[.5 .05 .5])
plot(2+linspace(-3,3),normcdf(linspace(-3,3)),'--','color',[.5 .05 .5])
[y,x] = ecdf(t10(t10~=0));
p3 = plot(4+x,y,'linewidth',2,'color',[.3 .1 .7])
plot(4+linspace(-3,3),normcdf(linspace(-3,3)),'--','color',[.3 .1 .7])
[y,x] = ecdf(t20(t20~=0));
p4 = plot(6+x,y,'linewidth',2,'color',[.05 .5 .5])
plot(6+linspace(-3,3),normcdf(linspace(-3,3)),'--','color',[.05 .5 .5])
[y,x] = ecdf(t95(t95~=0));
p5 = plot(8+x,y,'linewidth',2,'color',[.1 .7 .3])
plot(8+linspace(-3,3),normcdf(linspace(-3,3)),'--','color',[.1 .7 .3])
[y,x] = ecdf(t99(t99~=0));
p6 = plot(10+x,y,'linewidth',2,'color',[.5 .5 .05])
plot(10+linspace(-3,3),normcdf(linspace(-3,3)),'--','color',[.5 .5 .05])
clear x y;
lgnd = legend([p1 p2 p3 p4 p5 p6 p7],'Max 1-day','Max 5-day','Days $\geq$10mm','Days $\geq$20mm','\% in $\geq$95th \%ile Days','\% in $\geq$99th \%ile Days','Null Expectations');
set(lgnd,'interpreter','latex','fontsize',16,'location','northwest')
title('Annual Asian Precipitation Extreme Metrics','interpreter','latex','fontsize',16)
ylabel('Cumulative Distribution Function','interpreter','latex','fontsize',16)
xlabel('$z$-Scores of Grid Cell Trends (Shifted Left/Right for Visualization)','interpreter','latex','fontsize',16)
set(gca,'ytick',[],'xtick',[])
%% show r1 anomalies pretty much all in Irkutsk
x = ncread('HadEX3_61-90_R95p_ANN.nc','longitude');
y = ncread('HadEX3_61-90_R95p_ANN.nc','latitude');
y = y(79:125);
x = x(134:173);
x = x(1:39);
x = repmat(x,1,46);
y = y(1:46);
y = repmat(y,1,39)';
scatter(x(t10<-1.96),y(t10<-1.96));
hold on; 
scatter(x(t1<-1.96),y(t1<-1.96));
scatter(x(t5<-1.96),y(t5<-1.96));
