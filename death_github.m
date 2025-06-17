%% clear workspace, load data, process data
clear all; close all; clc;
opts = spreadsheetImportOptions("NumVariables", 32);
opts.Sheet = "EM-DAT Data";
opts.DataRange = "A1:AF26939";
opts.VariableNames = ["Var1", "Var2", "Var3", "c1", "c2", "c3", "Var7", "Var8", "Var9", "Var10", "Var11", "r2", "r1", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "y", "m", "Var28", "Var29", "Var30", "Var31", "d"];
opts.SelectedVariableNames = ["c1", "c2", "c3", "r2", "r1", "y", "m", "d"];
opts.VariableTypes = ["char", "char", "char", "categorical", "categorical", "categorical", "char", "char", "char", "char", "char", "categorical", "categorical", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "double", "char", "char", "char", "char", "double"];
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var7", "Var8", "Var9", "Var10", "Var11", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var28", "Var29", "Var30", "Var31"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "c1", "c2", "c3", "Var7", "Var8", "Var9", "Var10", "Var11", "r2", "r1", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var28", "Var29", "Var30", "Var31"], "EmptyFieldRule", "auto");
tbl = readtable("/Users/bbcael/Documents/work/dogs/mort/public_emdat_incl_hist_2025-04-15.xlsx", opts, "UseExcel", false);
clear opts;
c1 = tbl.c1; % categories
c2 = tbl.c2;
c3 = tbl.c3;
r2 = tbl.r2; % regions
r1 = tbl.r1;
y = tbl.y;
d = tbl.d;
m = tbl.m;
clear tbl;
% start in 1988: most (82%) of recorded events, also when cataloging started; can see this with ecdf(y)
r2 = r2(y>1987 & y<2025); m = m(y>1987 & y<2025); r1 = r1(y>1987 & y<2025); c2 = c2(y>1987 & y<2025); c1 = c1(y>1987 & y<2025); c3 = c3(y>1987 & y<2025); d = d(y>1987 & y<2025); y = y(y>1987 & y<2025);
y = y(~isnan(d)); m = m(~isnan(d)); r1 = r1(~isnan(d)); r2 = r2(~isnan(d)); c1 = c1(~isnan(d)); c2 = c2(~isnan(d)); c3 = c3(~isnan(d)); d = d(~isnan(d));
y = y(c1=="Natural"); m = m(c1=="Natural"); r1 = r1(c1=="Natural"); r2 = r2(c1=="Natural"); c2 = c2(c1=="Natural"); c3 = c3(c1=="Natural"); d = d(c1=="Natural"); clear c1; % remove technological disasters
y = y(c2~="Biological"); m = m(c2~="Biological"); r1 = r1(c2~="Biological"); r2 = r2(c2~="Biological"); c3 = c3(c2~="Biological"); d = d(c2~="Biological"); c2 = c2(c2~="Biological"); % remove biological disasters
y = y(c2~="Geophysical"); m = m(c2~="Geophysical"); r1 = r1(c2~="Geophysical"); r2 = r2(c2~="Geophysical"); c3 = c3(c2~="Geophysical"); d = d(c2~="Geophysical"); c2 = c2(c2~="Geophysical"); clear c2; % remove geophysical disasters
% only look at extreme temperatures, floods, & storms: most of the weather/climate deaths
y = y(c3~="Glacial lake outburst flood"); m = m(c3~="Glacial lake outburst flood"); r1 = r1(c3~="Glacial lake outburst flood"); r2 = r2(c3~="Glacial lake outburst flood"); d = d(c3~="Glacial lake outburst flood"); c3 = c3(c3~="Glacial lake outburst flood"); 
y = y(c3~="Wildfire"); m = m(c3~="Wildfire"); r1 = r1(c3~="Wildfire"); r2 = r2(c3~="Wildfire"); d = d(c3~="Wildfire"); c3 = c3(c3~="Wildfire"); 
y = y(c3~="Drought"); m = m(c3~="Drought"); r1 = r1(c3~="Drought"); r2 = r2(c3~="Drought"); d = d(c3~="Drought"); c3 = c3(c3~="Drought"); 
y = y(c3~="Mass movement (wet)"); m = m(c3~="Mass movement (wet)"); r1 = r1(c3~="Mass movement (wet)"); r2 = r2(c3~="Mass movement (wet)"); d = d(c3~="Mass movement (wet)"); c3 = c3(c3~="Mass movement (wet)"); 
c = zeros(size(c3)); c(c3=="Flood") = 1; c(c3=="Storm") = 2; clear c3;
c = c(r1~="Oceania"); m = m(r1~="Oceania"); d = d(r1~="Oceania"); y = y(r1~="Oceania"); r2 = r2(r1~="Oceania"); r1 = r1(r1~="Oceania"); 
% only look at 30+ deaths: 95% of deaths
y = y(d>29.9); m = m(d>29.9); r2 = r2(d>29.9); r1 = r1(d>29.9); c = c(d>29.9); d = d(d>29.9); 
y = y-1988;
r = zeros(size(r1)); r(r2=="Latin America and the Caribbean") = 3; r(r1=="Europe") = 2; r(r2=="Northern America") = 4; r(r1=="Asia") = 1;
clear r1 r2;
%% figure S1: data are variable, need a distribution
figure
colormap(turbo)
subplot(133)
scatter(1988+y(c==0),d(c==0),100,r(c==0),'filled','markerfacealpha',0.75);
set(gca,'yscale','log','ticklabelinterpreter','latex','fontsize',18,'ytick',[100 1000 10000 100000 1000000])
box on;
axis([1987.1 2024.9 25 200000])
title('Temperatures','interpreter','latex')
subplot(131)
scatter(1988+y(c==1),d(c==1),100,r(c==1),'filled','markerfacealpha',0.75);
set(gca,'yscale','log','ticklabelinterpreter','latex','fontsize',18,'ytick',[100 1000 10000 100000 1000000])
box on;
axis([1987.1 2024.9 25 200000])
title('Floods','interpreter','latex')
hold on;
l0 = scatter(-100,1e9,100,0,'filled','markerfacealpha',0.75);
l1 = scatter(-100,1e9,100,1,'filled','markerfacealpha',0.75);
l2 = scatter(-100,1e9,100,2,'filled','markerfacealpha',0.75);
l3 = scatter(-100,1e9,100,3,'filled','markerfacealpha',0.75);
l4 = scatter(-100,1e9,100,4,'filled','markerfacealpha',0.75);
lgnd = legend([l0 l1 l2 l3 l4],'Africa','Asia','Europe','Latin America','North America')
set(lgnd,'interpreter','latex','fontsize',16,'location','northeast')
ylabel('Deaths','interpreter','latex')
subplot(132)
scatter(1988+y(c==2),d(c==2),100,r(c==2),'filled','markerfacealpha',0.75);
set(gca,'yscale','log','ticklabelinterpreter','latex','fontsize',18,'ytick',[100 1000 10000 100000 1000000])
box on;
axis([1987.1 2024.9 25 200000])
title('Storms','interpreter','latex')
clear ans l0 l1 l2 l3 l4 lgnd;

%% figure 1: GP is a good fit

dmin = 30; cval = 1;
[phat,pci] = mle(d(c==cval & d>dmin)-dmin,'distribution','gp')
[Y,X] = ecdf(sort(d(c==cval & d>dmin))); X = X(1:end-1); 
ec = Y; Y = Y(1:end-1);
Y = gpcdf(sort(d(c==cval & d>dmin)),phat(1),phat(2),dmin); tc = unique(Y);
[~,pval] = kstest2(ec(1:end-1),tc)
Y = gpinv(unique(Y),phat(1),phat(2),dmin);
plot(1:5e5,1:5e5,'k','linewidth',2)
hold on;
p1 = scatter(X,Y,50,[.7 .35 .1],'filled');
set(gca,'xscale','log','yscale','log')
box on;
set(gca,'ticklabelinterpreter','latex','fontsize',16)
xlabel('Empirical Quantiles','interpreter','latex')
ylabel('GP Quantiles','interpreter','latex')
clear ans pci phat X Y err pval ec tc cval dmin;

dmin = 30; cval = 2;
[phat,pci] = mle(d(c==cval & d>dmin)-dmin,'distribution','gp')
[Y,X] = ecdf(sort(d(c==cval & d>dmin))); X = X(1:end-1); 
ec = Y; Y = Y(1:end-1);
Y = gpcdf(sort(d(c==cval & d>dmin)),phat(1),phat(2),dmin); tc = unique(Y);
[~,pval] = kstest2(ec(1:end-1),tc)
Y = gpinv(unique(Y),phat(1),phat(2),dmin);
plot(3.*(1:5e5),1:5e5,'k','linewidth',2)
p2 = scatter(3.*X,Y,50,[.5 .05 .5],'filled');
clear ans pci phat X Y err pval ec tc cval dmin;

dmin = 30; cval = 0;
[phat,pci] = mle(d(c==cval & d>dmin)-dmin,'distribution','gp')
[Y,X] = ecdf(sort(d(c==cval & d>dmin))); X = X(1:end-1); 
ec = Y; Y = Y(1:end-1);
Y = gpcdf(sort(d(c==cval & d>dmin)),phat(1),phat(2),dmin); tc = unique(Y);
[~,pval] = kstest2(ec(1:end-1),tc)
Y = gpinv(unique(Y),phat(1),phat(2),dmin);
plot(9.*(1:5e5),1:5e5,'k','linewidth',2)
p3 = scatter(9.*X,Y,50,[.05 .5 .5],'filled');
clear ans pci phat X Y err pval ec tc cval dmin;

axis([25 5e5 25 2e5])
lgnd = legend([p1 p2 p3],'Floods','Storms','Temperatures')
set(lgnd,'interpreter','latex','location','northwest','fontsize',16)
axis square;

clear p1 p2 p3 lgnd;

%% check which ones are changing with time

D = d(c==0 & r==0); Y = y(c==0 & r==0); % too few

D = d(c==1 & r==0); Y = y(c==1 & r==0); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==2 & r==0); Y = y(c==2 & r==0); % African storms: +scale
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==0 & r==1); Y = y(c==0 & r==1); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==1 & r==1); Y = y(c==1 & r==1); % Asian floods: -shape
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==2 & r==1); Y = y(c==2 & r==1); % Asian storms: -shape
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==0 & r==2); Y = y(c==0 & r==2); % European temps: +shape -- robust to excluding 2022+2023
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==1 & r==2); Y = y(c==1 & r==2); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==2 & r==2); Y = y(c==2 & r==2); % not enough

D = d(c==0 & r==3); Y = y(c==0 & r==3); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==1 & r==3); Y = y(c==1 & r==3); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==2 & r==3); Y = y(c==2 & r==3); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==0 & r==4); Y = y(c==0 & r==4); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==1 & r==4); Y = y(c==1 & r==4); % not enough

D = d(c==2 & r==4); Y = y(c==2 & r==4); % nope
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

% European temps: +shape: due to European warming -- fig 3
% Asian storms & floods: -shape: not climate or exposure, has to be reduced vulnerability -- fig 4
% African storms: +scale: due to storm daniel -- fig 5

clear ans D pci phat Y;

%% double-check with quantile regression

D = d(c==1 & r==1); Y = y(c==1 & r==1); % Asian floods: -shape
for i = 1:99;
    [p,stats] = quantreg(Y,D,i./100);
    P(i) = sum(stats.pboot(:,1)>0);
    i
end
P_AF = P; % 77 of the 13th-96th %iles decreasing with 95% confidence
D = d(c==2 & r==1); Y = y(c==2 & r==1); % Asian storms: -shape
for i = 1:99;
    [p,stats] = quantreg(Y,D,i./100);
    P(i) = sum(stats.pboot(:,1)>0);
    i
end
P_AS = P; % 55 of the 23-84th %iles decreasing with 95% confidence
D = d(c==0 & r==2); Y = y(c==0 & r==2); % European temps: +shape -- robust to excluding 2022+2023
for i = 1:99;
    [p,stats] = quantreg(Y,D,i./100);
    P(i) = sum(stats.pboot(:,1)>0);
    i
end
P_ET = P; % 50 of the 15-83rd %iles increasing with 95% confidence
clear i D Y P p stats ans;

%% double-check with bootstrapping

nb = 200;
for i = 1:nb;
    D = d(c==1 & r==1); Y = y(c==1 & r==1);
    ind = randi(length(D),1,length(D));
    D = D(ind); Y = Y(ind);
    phat = mle(D-29.9,'distribution','gp');
    [phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
    tAF(i) = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
    clear D ind Y phat pci;
    D = d(c==2 & r==1); Y = y(c==2 & r==1);
    ind = randi(length(D),1,length(D));
    D = D(ind); Y = Y(ind);
    phat = mle(D-29.9,'distribution','gp');
    [phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
    tAS(i) = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
    clear D ind Y phat pci;
    D = d(c==0 & r==2); Y = y(c==0 & r==2);
    ind = randi(length(D),1,length(D));
    D = D(ind); Y = Y(ind);
    phat = mle(D-29.9,'distribution','gp');
    [phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
    tET(i) = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
    clear D ind Y phat pci;
    i
end
clear i nb ans

%%

D = d(c==2 & r==1); Y = y(c==2 & r==1); % Asian storms: -shape
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))

D = d(c==0 & r==2); Y = y(c==0 & r==2); % European temps: +shape -- robust to excluding 2022+2023
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))


%% figure S2: show European Temps due to climate change

D = d(c==0 & r==2);
M = m(c==0 & r==2);
Y = y(c==0 & r==2);
C = zeros(size(M)); C(M>4 & M<9) = 1;
scatter(Y(C==0)+1988,D(C==0),100,[.05 .05 .5],'filled','markerfacealpha',0.75)
hold on;
scatter(Y(C==1)+1988,D(C==1),100,[.5 .05 .05],'filled','markerfacealpha',0.75)
set(gca,'yscale','log','ticklabelinterpreter','latex','fontsize',16)
box on;
clear D M Y C;
axis([1987.1 2024.9 25 70000])
ylabel('Deaths','interpreter','latex')
title('European Temperatures','interpreter','latex')
lgnd = legend('Oct-Feb','May-Aug');
set(lgnd,'interpreter','latex','fontsize',16,'location','northwest')
clear lgnd;
axis square;
clear pci phat ans m;

%% figure 2: show Asian storms+floods due to vulnerability

subplot(121)
D = d(c>0 & r==1);
Y = y(c>0 & r==1);
C = c(c>0 & r==1);
scatter(Y(C==1)+1988,D(C==1),100,[.7 .35 .1],'filled','markerfacealpha',0.5)
hold on;
scatter(Y(C==2)+1988,D(C==2),100,[.5 .05 .5],'filled','markerfacealpha',0.5)
set(gca,'yscale','log','ticklabelinterpreter','latex','fontsize',16)
box on;
clear D M Y C;
axis([1987.1 2024.9 25 200000])
ylabel('Deaths','interpreter','latex')
title('Asian Floods + Storms','interpreter','latex')
lgnd = legend('Floods','Storms');
set(lgnd,'interpreter','latex','fontsize',16,'location','northwest')
clear lgnd;
subplot(122)
D = d(c>0 & r==1); Y = y(c>0 & r==1); % Asian floods: -shape
[phat,pci] = mle(D-29.9,'distribution','gp')
[p,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
clear D Y pci ans;
X = 30:1000;
g0 = gppdf(X,p(1),p(3),30);
g1 = gppdf(X,p(1)+36.*p(4),p(3)+36.*p(2),30);
plot(X,1.3.*g0,'--','color',[.5 .5 .05],'linewidth',3)
hold on;
plot(X,g1,'color',[.05 .5 .5],'linewidth',3)
set(gca,'xscale','log','ytick',[],'ticklabelinterpreter','latex','fontsize',16)
xlabel('Deaths','interpreter','latex')
clear p g0 g1 ans phat X;
ylabel('Probability Density','interpreter','latex')
lgnd = legend('GP(1988)','GP(2024)')
set(lgnd,'interpreter','latex','fontsize',16','location','northeast')
axis([-Inf Inf -Inf Inf])
clear lgnd;

%% african nonstationarity driven entirely by Storm Daniel

D = d(c>0 & r==0);
Y = y(c>0 & r==0);
[~,ind] = max(D);
D(ind) = []; Y(ind) = [];
[phat,pci] = mle(D-29.9,'distribution','gp')
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
3.96.*phat(4)./(pci(2,4)-pci(1,4))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
phat = mle(D-29.9,'distribution','gp');
[phat,pci] = mle(D-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
3.96.*phat(2)./(pci(2,2)-pci(1,2))
clear D Y phat pci ind ans;

%% figure 3: show African storms due to Daniel

subplot(121)
D = d(c>0 & r==0);
Y = y(c>0 & r==0);
C = c(c>0 & r==0);
scatter(Y(C==1)+1988,D(C==1),100,[.7 .35 .1],'filled','markerfacealpha',0.5)
hold on;
scatter(Y(C==2)+1988,D(C==2),100,[.5 .05 .5],'filled','markerfacealpha',0.5)
set(gca,'yscale','log','ticklabelinterpreter','latex','fontsize',16)
box on;
clear C;
axis([1987.1 2024.9 25 20000])
ylabel('Deaths','interpreter','latex')
title('African Floods + Storms','interpreter','latex')
lgnd = legend('Floods','Storms');
set(lgnd,'interpreter','latex','fontsize',16,'location','northwest')
clear lgnd;
subplot(122)
[phat,pci] = mle(D(Y<35)-29.9,'distribution','gp')
pci = (pci(2,:)-pci(1,:))./3.96;
Dboot = [];
for i = 1:1e6;
    R = 1.30+randn(1).*.144+35.*(.0335+randn(1).*.00597);
    N = poissrnd(exp(R));
    A = phat(1)+randn(1).*pci(1);
    B = phat(2)+randn(1).*pci(2);
    dboot = gprnd(A,B,30,1,N);
    if length(dboot)>0
    Dboot(end+1) = max(dboot);
    else
    Dboot(end+1) = 0;
    end
end
clear i A B pci phat dboot R N Y Rboot ind;
1./(sum(Dboot>max(D))./length(Dboot))
[Y,X] = ksdensity(log10(Dboot(Dboot>0)));
plot(X,Y,'k','linewidth',3)
set(gca,'xtick',[2 3 4],'xticklabel',{"100","1000","10000"},'ticklabelinterpreter','latex','fontsize',16,'ytick',[])
axis([1.5 4.5 0 1.15])
xlabel('Deaths','interpreter','latex')
ylabel('Probability Density','interpreter','latex')
hold on;
plot(log10(linspace(max(D),max(D))),linspace(0,1.2),'--','linewidth',3,'color',[.5 .05 .5])
lgnd = legend('Deadliest 2023 Event','Storm Daniel');
set(lgnd,'interpreter','latex','fontsize',16,'location','northeast')
clear Dboot lgnd X Y ans D;

%% check which ones' frequencies are also changing

clear D Y;
Y = unique(y); 
for i = 1:length(Y); 
    %D(i) = sum(y==Y(i) & c>0 & r==0); 
    %D(i) = sum(y==Y(i) & c==0 & r==0); % N = 5
    %D(i) = sum(y==Y(i) & c==1 & r==0); % african floods increasing
    %D(i) = sum(y==Y(i) & c==2 & r==0); % not significant after bonferroni correction
    %D(i) = sum(y==Y(i) & c==0 & r==1); 
    %D(i) = sum(y==Y(i) & c==1 & r==1);
    %D(i) = sum(y==Y(i) & c==2 & r==1); % asian storms decreasing
    %D(i) = sum(y==Y(i) & c==0 & r==2); % driven by outlier year-pair -- not significant after bonferroni correction without these
    %D(i) = sum(y==Y(i) & c==1 & r==2); % N = 21
    %D(i) = sum(y==Y(i) & c==2 & r==2); % N = 9
    %D(i) = sum(y==Y(i) & c==0 & r==3); % N = 20
    %D(i) = sum(y==Y(i) & c==1 & r==3);
    %D(i) = sum(y==Y(i) & c==2 & r==3);
    %D(i) = sum(y==Y(i) & c==0 & r==4); % N = 16
    %D(i) = sum(y==Y(i) & c==1 & r==4); % N = 9
    %D(i) = sum(y==Y(i) & c==2 & r==4);
end
clear i;
mdl = fitglm(Y,D,'linear','Distribution','poisson')
%sum(D)
scatter(Y,D,100,'filled')
hold on;
q = table2array(mdl.Coefficients(1:2,1));
plot(Y,exp(q(1)+Y.*q(2)),'k','linewidth',3)
clear mdl D Y ans q;

%% fig. S3: show changes in frequency

Y = unique(y); 
for i = 1:length(Y); 
    D1(i) = sum(y==Y(i) & c==1 & r==0); 
    D2(i) = sum(y==Y(i) & c==2 & r==1);
end
clear i;
mdl = fitglm(Y,D1,'linear','Distribution','poisson')
subplot(131)
scatter(Y+1988,D1,100,[.7 .35 .1],'filled')
hold on;
q = table2array(mdl.Coefficients(1:2,1));
plot(Y+1988,exp(q(1)+Y.*q(2)),'k','linewidth',3)
clear mdl D1 ans q;
box on;
set(gca,'ticklabelinterpreter','latex','fontsize',16)
ylabel('Disasters per year','interpreter','latex')
axis([1987.1 2024.9 0 15])
title('African Floods','interpreter','latex')

mdl = fitglm(Y,D2,'linear','Distribution','poisson')
subplot(133)
scatter(Y+1988,D2,100,[.5 .05 .5],'filled')
hold on;
q = table2array(mdl.Coefficients(1:2,1));
plot(Y+1988,exp(q(1)+Y.*q(2)),'k','linewidth',3)
clear mdl D2 ans q;
set(gca,'ticklabelinterpreter','latex','fontsize',16)
box on;
axis([1987.1 2024.9 0 23])
title('Asian Storms','interpreter','latex')

%%

clear all; clc;
opts = spreadsheetImportOptions("NumVariables", 32);
opts.Sheet = "EM-DAT Data";
opts.DataRange = "A1:AF26939";
opts.VariableNames = ["Var1", "Var2", "Var3", "c1", "c2", "c3", "Var7", "Var8", "Var9", "Var10", "Var11", "r2", "r1", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "y", "m", "Var28", "Var29", "Var30", "Var31", "d"];
opts.SelectedVariableNames = ["c1", "c2", "c3", "r2", "r1", "y", "m", "d"];
opts.VariableTypes = ["char", "char", "char", "categorical", "categorical", "categorical", "char", "char", "char", "char", "char", "categorical", "categorical", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "double", "char", "char", "char", "char", "double"];
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var7", "Var8", "Var9", "Var10", "Var11", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var28", "Var29", "Var30", "Var31"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "c1", "c2", "c3", "Var7", "Var8", "Var9", "Var10", "Var11", "r2", "r1", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var28", "Var29", "Var30", "Var31"], "EmptyFieldRule", "auto");
tbl = readtable("/Users/bbcael/Documents/work/dogs/mort/public_emdat_incl_hist_2025-04-15.xlsx", opts, "UseExcel", false);
clear opts;
c1 = tbl.c1; % categories
c2 = tbl.c2;
c3 = tbl.c3;
r2 = tbl.r2; % regions
r1 = tbl.r1;
y = tbl.y;
d = tbl.d;
m = tbl.m;
clear tbl;
% start in 1988: most (82%) of recorded events, also when cataloging started; can see this with ecdf(y)
r2 = r2(y>1987 & y<2025); m = m(y>1987 & y<2025); r1 = r1(y>1987 & y<2025); c2 = c2(y>1987 & y<2025); c1 = c1(y>1987 & y<2025); c3 = c3(y>1987 & y<2025); d = d(y>1987 & y<2025); y = y(y>1987 & y<2025);
y = y(~isnan(d)); m = m(~isnan(d)); r1 = r1(~isnan(d)); r2 = r2(~isnan(d)); c1 = c1(~isnan(d)); c2 = c2(~isnan(d)); c3 = c3(~isnan(d)); d = d(~isnan(d));
y = y(c1=="Natural"); m = m(c1=="Natural"); r1 = r1(c1=="Natural"); r2 = r2(c1=="Natural"); c2 = c2(c1=="Natural"); c3 = c3(c1=="Natural"); d = d(c1=="Natural"); clear c1; % remove technological disasters
y = y(c2~="Biological"); m = m(c2~="Biological"); r1 = r1(c2~="Biological"); r2 = r2(c2~="Biological"); c3 = c3(c2~="Biological"); d = d(c2~="Biological"); c2 = c2(c2~="Biological"); % remove biological disasters
y = y(c2~="Geophysical"); m = m(c2~="Geophysical"); r1 = r1(c2~="Geophysical"); r2 = r2(c2~="Geophysical"); c3 = c3(c2~="Geophysical"); d = d(c2~="Geophysical"); c2 = c2(c2~="Geophysical"); clear c2; % remove geophysical disasters
% only look at extreme temperatures, floods, & storms: most of the weather/climate deaths
y = y(c3~="Glacial lake outburst flood"); m = m(c3~="Glacial lake outburst flood"); r1 = r1(c3~="Glacial lake outburst flood"); r2 = r2(c3~="Glacial lake outburst flood"); d = d(c3~="Glacial lake outburst flood"); c3 = c3(c3~="Glacial lake outburst flood"); 
y = y(c3~="Wildfire"); m = m(c3~="Wildfire"); r1 = r1(c3~="Wildfire"); r2 = r2(c3~="Wildfire"); d = d(c3~="Wildfire"); c3 = c3(c3~="Wildfire"); 
y = y(c3~="Drought"); m = m(c3~="Drought"); r1 = r1(c3~="Drought"); r2 = r2(c3~="Drought"); d = d(c3~="Drought"); c3 = c3(c3~="Drought"); 
y = y(c3~="Mass movement (wet)"); m = m(c3~="Mass movement (wet)"); r1 = r1(c3~="Mass movement (wet)"); r2 = r2(c3~="Mass movement (wet)"); d = d(c3~="Mass movement (wet)"); c3 = c3(c3~="Mass movement (wet)"); 
c = zeros(size(c3)); c(c3=="Flood") = 1; c(c3=="Storm") = 2; clear c3;
c = c(r1~="Oceania"); m = m(r1~="Oceania"); d = d(r1~="Oceania"); y = y(r1~="Oceania"); r2 = r2(r1~="Oceania"); r1 = r1(r1~="Oceania"); 
% only look at 30+ deaths: 95% of deaths
%y = y(d>29.9); m = m(d>29.9); r2 = r2(d>29.9); r1 = r1(d>29.9); c = c(d>29.9); d = d(d>29.9); 
y = y-1988;
r = zeros(size(r1)); r(r2=="Latin America and the Caribbean") = 3; r(r1=="Europe") = 2; r(r2=="Northern America") = 4; r(r1=="Asia") = 1;
clear r1 r2;
y = y(r==0 & c==1); d = d(r==0 & c==1); clear r m c;

thresh = 19.*exp(y.*.0253);
y = y(d>thresh);
d = d(d>thresh);
clear thresh;

Y = unique(y); 
for i = 1:length(Y); 
    D(i) = sum(y==Y(i));
end
mdl = fitglm(Y,D,'linear','Distribution','poisson')

subplot(132)
scatter(Y+1988,D,100,[.7 .35 .1],'filled')
%hold on;
%q = table2array(mdl.Coefficients(1:2,1));
%plot(Y+1988,exp(q(1)+Y.*q(2)),'k','linewidth',3)
clear mdl D ans q;
set(gca,'ticklabelinterpreter','latex','fontsize',16,'ytick',[0 5 10])
box on;
axis([1987.1 2024.9 0 10.9])
title('Population-Adjusted Threshold','interpreter','latex')
clear i d y Y;

%% fig. S4: lives saved calculation

clear all; close all; clc;
asiapop = [3087 3148 3210 3269 3326 3381 3435 3488 3541 3593 3645 3697 3748 3799 3848 3897 3945 3994 4043 4092 4140 4190 4240 4289 4340 4389 4437 4483 4528 4572 4614 4652 4688 4718 4748 4778 4807]./3087;
opts = spreadsheetImportOptions("NumVariables", 32);
opts.Sheet = "EM-DAT Data";
opts.DataRange = "A1:AF26939";
opts.VariableNames = ["Var1", "Var2", "Var3", "c1", "c2", "c3", "Var7", "Var8", "Var9", "Var10", "Var11", "r2", "r1", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "y", "m", "Var28", "Var29", "Var30", "Var31", "d"];
opts.SelectedVariableNames = ["c1", "c2", "c3", "r2", "r1", "y", "m", "d"];
opts.VariableTypes = ["char", "char", "char", "categorical", "categorical", "categorical", "char", "char", "char", "char", "char", "categorical", "categorical", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "double", "char", "char", "char", "char", "double"];
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var7", "Var8", "Var9", "Var10", "Var11", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var28", "Var29", "Var30", "Var31"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "c1", "c2", "c3", "Var7", "Var8", "Var9", "Var10", "Var11", "r2", "r1", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var28", "Var29", "Var30", "Var31"], "EmptyFieldRule", "auto");
tbl = readtable("/Users/bbcael/Documents/work/dogs/mort/public_emdat_incl_hist_2025-04-15.xlsx", opts, "UseExcel", false);
clear opts;
c1 = tbl.c1; % categories
c2 = tbl.c2;
c3 = tbl.c3;
r2 = tbl.r2; % regions
r1 = tbl.r1;
y = tbl.y;
d = tbl.d;
m = tbl.m;
clear tbl;
% start in 1988: most (82%) of recorded events, also when cataloging started; can see this with ecdf(y)
r2 = r2(y>1987 & y<2025); m = m(y>1987 & y<2025); r1 = r1(y>1987 & y<2025); c2 = c2(y>1987 & y<2025); c1 = c1(y>1987 & y<2025); c3 = c3(y>1987 & y<2025); d = d(y>1987 & y<2025); y = y(y>1987 & y<2025);
y = y(~isnan(d)); m = m(~isnan(d)); r1 = r1(~isnan(d)); r2 = r2(~isnan(d)); c1 = c1(~isnan(d)); c2 = c2(~isnan(d)); c3 = c3(~isnan(d)); d = d(~isnan(d));
y = y(c1=="Natural"); m = m(c1=="Natural"); r1 = r1(c1=="Natural"); r2 = r2(c1=="Natural"); c2 = c2(c1=="Natural"); c3 = c3(c1=="Natural"); d = d(c1=="Natural"); clear c1; % remove technological disasters
y = y(c2~="Biological"); m = m(c2~="Biological"); r1 = r1(c2~="Biological"); r2 = r2(c2~="Biological"); c3 = c3(c2~="Biological"); d = d(c2~="Biological"); c2 = c2(c2~="Biological"); % remove biological disasters
y = y(c2~="Geophysical"); m = m(c2~="Geophysical"); r1 = r1(c2~="Geophysical"); r2 = r2(c2~="Geophysical"); c3 = c3(c2~="Geophysical"); d = d(c2~="Geophysical"); c2 = c2(c2~="Geophysical"); clear c2; % remove geophysical disasters
% only look at extreme temperatures, floods, & storms: most of the weather/climate deaths
y = y(c3~="Glacial lake outburst flood"); m = m(c3~="Glacial lake outburst flood"); r1 = r1(c3~="Glacial lake outburst flood"); r2 = r2(c3~="Glacial lake outburst flood"); d = d(c3~="Glacial lake outburst flood"); c3 = c3(c3~="Glacial lake outburst flood"); 
y = y(c3~="Wildfire"); m = m(c3~="Wildfire"); r1 = r1(c3~="Wildfire"); r2 = r2(c3~="Wildfire"); d = d(c3~="Wildfire"); c3 = c3(c3~="Wildfire"); 
y = y(c3~="Drought"); m = m(c3~="Drought"); r1 = r1(c3~="Drought"); r2 = r2(c3~="Drought"); d = d(c3~="Drought"); c3 = c3(c3~="Drought"); 
y = y(c3~="Mass movement (wet)"); m = m(c3~="Mass movement (wet)"); r1 = r1(c3~="Mass movement (wet)"); r2 = r2(c3~="Mass movement (wet)"); d = d(c3~="Mass movement (wet)"); c3 = c3(c3~="Mass movement (wet)"); 
c = zeros(size(c3)); c(c3=="Flood") = 1; c(c3=="Storm") = 2; clear c3;
c = c(r1~="Oceania"); m = m(r1~="Oceania"); d = d(r1~="Oceania"); y = y(r1~="Oceania"); r2 = r2(r1~="Oceania"); r1 = r1(r1~="Oceania"); 
% only look at 30+ deaths: 95% of deaths
y = y-1988;
r = zeros(size(r1)); r(r2=="Latin America and the Caribbean") = 3; r(r1=="Europe") = 2; r(r2=="Northern America") = 4; r(r1=="Asia") = 1;
clear r1 r2;
y = y(r==1 & c>0); d = d(r==1 & c>0); clear r m c;
yl = y(d<30); dl = d(d<30);
y = y(d>29.9); d = d(d>29.9); pl = asiapop(yl+1)';
D = sum(d);
pop = asiapop(y+1)'; clear asiapop; % to add exposure to counterfactual
[phat,pci] = mle(d-29.9,'distribution','gp');
nb = 10000; % bootstrap iterations
for i = 1:nb;
    ib = randsample(length(y),length(y),true);
    db = d(ib);
    yb = y(ib);
    popb = pop(ib);
    [p,pci] = mle(db-29.9,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*yb)),'Start',[phat(1) 0 phat(2)]);
    vul = p(3)./(p(3)+p(2).*yb); % to add vulnerability to counterfactual
    D_alt(i) = sum(d.*pop.*vul);
    dlb = dl.*(p(3)./(p(3)+p(2).*yl)).*pl;
    Dl(i) = sum(dlb(dlb>29.9)); % deaths from events that cause 30+ deaths that wouldn't've otherwise
    i
end
D_alt = D_alt+Dl;
[Y,X] = ksdensity(D_alt,min(D_alt):1e4:max(D_alt));
p2 = plot(X,smooth(Y),'linewidth',3,'color',[.5 .05 .5]);
hold on;
p1 = plot(linspace(D,D),linspace(-1,2.*max(Y)),'--k','linewidth',3)
axis([.9.*D 1.02.*max(X) 0 1.02.*max(Y)])
box on;
set(gca,'ticklabelinterpreter','latex','fontsize',16,'ytick',[]);
xlabel('Deaths from Asian Floods and Storms','interpreter','latex')
ylabel('Probability Density','interpreter','latex')
lgnd = legend([p1 p2],'Recorded','Counterfactual')
set(lgnd,'interpreter','latex','fontsize',16,'location','northeast')