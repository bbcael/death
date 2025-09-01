%% clear workspace, load data, set sensitivity parameters, process data

clear all; close all; clc;

YEAR = 1988-1; % default is 1988
THRESH = 30-.1; % default is 30
YEAR2 = 2025; % default is 2025: can show results are robust to excluding 2019-2024

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
r2 = r2(y>YEAR & y<YEAR2); m = m(y>YEAR & y<YEAR2); r1 = r1(y>YEAR & y<YEAR2); c2 = c2(y>YEAR & y<YEAR2); c1 = c1(y>YEAR & y<YEAR2); c3 = c3(y>YEAR & y<YEAR2); d = d(y>YEAR & y<YEAR2); y = y(y>YEAR & y<YEAR2);
y = y(~isnan(d)); m = m(~isnan(d)); r1 = r1(~isnan(d)); r2 = r2(~isnan(d)); c1 = c1(~isnan(d)); c2 = c2(~isnan(d)); c3 = c3(~isnan(d)); d = d(~isnan(d));
y = y(c1=="Natural"); m = m(c1=="Natural"); r1 = r1(c1=="Natural"); r2 = r2(c1=="Natural"); c2 = c2(c1=="Natural"); c3 = c3(c1=="Natural"); d = d(c1=="Natural"); clear c1; % remove technological disasters
y = y(c2~="Biological"); m = m(c2~="Biological"); r1 = r1(c2~="Biological"); r2 = r2(c2~="Biological"); c3 = c3(c2~="Biological"); d = d(c2~="Biological"); c2 = c2(c2~="Biological"); % remove biological disasters
y = y(c2~="Geophysical"); m = m(c2~="Geophysical"); r1 = r1(c2~="Geophysical"); r2 = r2(c2~="Geophysical"); c3 = c3(c2~="Geophysical"); d = d(c2~="Geophysical"); c2 = c2(c2~="Geophysical"); clear c2; % remove geophysical disasters
y = y(c3~="Glacial lake outburst flood"); m = m(c3~="Glacial lake outburst flood"); r1 = r1(c3~="Glacial lake outburst flood"); r2 = r2(c3~="Glacial lake outburst flood"); d = d(c3~="Glacial lake outburst flood"); c3 = c3(c3~="Glacial lake outburst flood"); 
y = y(c3~="Wildfire"); m = m(c3~="Wildfire"); r1 = r1(c3~="Wildfire"); r2 = r2(c3~="Wildfire"); d = d(c3~="Wildfire"); c3 = c3(c3~="Wildfire"); 
y = y(c3~="Drought"); m = m(c3~="Drought"); r1 = r1(c3~="Drought"); r2 = r2(c3~="Drought"); d = d(c3~="Drought"); c3 = c3(c3~="Drought"); 
y = y(c3~="Mass movement (wet)"); m = m(c3~="Mass movement (wet)"); r1 = r1(c3~="Mass movement (wet)"); r2 = r2(c3~="Mass movement (wet)"); d = d(c3~="Mass movement (wet)"); c3 = c3(c3~="Mass movement (wet)"); 
c = zeros(size(c3)); c(c3=="Flood") = 1; c(c3=="Storm") = 2; clear c3;
c = c(r1~="Oceania"); m = m(r1~="Oceania"); d = d(r1~="Oceania"); y = y(r1~="Oceania"); r2 = r2(r1~="Oceania"); r1 = r1(r1~="Oceania"); 
y = y(d>THRESH); m = m(d>THRESH); r2 = r2(d>THRESH); r1 = r1(d>THRESH); c = c(d>THRESH); d = d(d>THRESH); 
y = y-(YEAR+1);
r = zeros(size(r1)); r(r2=="Latin America and the Caribbean") = 3; r(r1=="Europe") = 2; r(r2=="Northern America") = 4; r(r1=="Asia") = 1;
clear r1 r2;

%% to make sure GP is still a good fit

dmin = min(d); cval = 1;
[phat,pci] = mle(d(c==cval & d>dmin)-dmin,'distribution','gp')
[Y,X] = ecdf(sort(d(c==cval & d>dmin))); X = X(1:end-1); 
ec = Y; Y = Y(1:end-1);
Y = gpcdf(sort(d(c==cval & d>dmin)),phat(1),phat(2),dmin); tc = unique(Y);
[~,~,KS1] = kstest2(ec(1:end-1),tc)

dmin = min(d); cval = 2;
[phat,pci] = mle(d(c==cval & d>dmin)-dmin,'distribution','gp')
[Y,X] = ecdf(sort(d(c==cval & d>dmin))); X = X(1:end-1); 
ec = Y; Y = Y(1:end-1);
Y = gpcdf(sort(d(c==cval & d>dmin)),phat(1),phat(2),dmin); tc = unique(Y);
[~,~,KS2] = kstest2(ec(1:end-1),tc)

dmin = min(d); cval = 0;
[phat,pci] = mle(d(c==cval & d>dmin)-dmin,'distribution','gp')
[Y,X] = ecdf(sort(d(c==cval & d>dmin))); X = X(1:end-1); 
ec = Y; Y = Y(1:end-1);
Y = gpcdf(sort(d(c==cval & d>dmin)),phat(1),phat(2),dmin); tc = unique(Y);
[~,~,KS3] = kstest2(ec(1:end-1),tc)

clc;

clear cval dmin pci phat X Y;

%% to make sure intensity trends are still significant

D = d(c>0 & r==1); Y = y(c>0 & r==1); % Asian floods: -shape
[phat,pci] = mle(D-THRESH,'distribution','gp')
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
sfs1 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
sfs2 = 3.96.*phat(4)./(pci(2,4)-pci(1,4));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
sfs3 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
sfs4 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));

D = d(c==1 & r==1); Y = y(c==1 & r==1); % Asian floods: -shape
[phat,pci] = mle(D-THRESH,'distribution','gp')
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
sf1 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
sf2 = 3.96.*phat(4)./(pci(2,4)-pci(1,4));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
sf3 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
sf4 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));

D = d(c==2 & r==1); Y = y(c==2 & r==1); % Asian storms: -shape
[phat,pci] = mle(D-THRESH,'distribution','gp')
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
ss1 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
ss2 = 3.96.*phat(4)./(pci(2,4)-pci(1,4));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
ss3 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
ss4 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));

D = d(c==0 & r==2); Y = y(c==0 & r==2); % European temps: +shape -- robust to excluding 2022+2023
[phat,pci] = mle(D-THRESH,'distribution','gp')
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c,d)gppdf(x,a+d.*Y,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2) 0]);
et1 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
et2 = 3.96.*phat(4)./(pci(2,4)-pci(1,4));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a+b.*Y,abs(c)),'Start',[phat(1) 0 phat(2)]);
et3 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));
phat = mle(D-THRESH,'distribution','gp');
[phat,pci] = mle(D-THRESH,'pdf',@(x,a,b,c)gppdf(x,a,abs(c+b.*Y)),'Start',[phat(1) 0 phat(2)]);
et4 = 3.96.*phat(2)./(pci(2,2)-pci(1,2));

clear ans D pci phat Y;
clc;

%% to check frequency trends are still significant 

clear D Y;
Y = unique(y); 
for i = 1:length(Y); 
    Dff(i) = sum(y==Y(i) & c==1 & r==0); % african floods increasing
    Dss(i) = sum(y==Y(i) & c==2 & r==1); % asian storms decreasing
end
clear i;
mdl = fitglm(Y,Dff,'linear','Distribution','poisson');
ffr1 = table2array(mdl.Coefficients(2,3));
mdl = fitglm(Y,Dss,'linear','Distribution','poisson');
ssr1 = table2array(mdl.Coefficients(2,3));

%% print output

VALS = [max([KS1 KS2 KS3]) max(abs([sfs1 sfs2 sfs3 sfs4])) max(abs([sf1 sf2 sf3 sf4])) max(abs([ss1 ss2 ss3 ss4])) max(abs([et1 et2 et3 et4])) ffr1 ssr1]

clearvars -EXCEPT VALS