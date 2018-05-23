%This codes plots the orientations from fcc_bcc_nrss_materialia on (100) pole figures colored according to their max, min, mean, and standard deviation nRSS values
%It makes use of the function jet_white for the colormap
%Garrison Hommer, 22MAY2018

figure
plotPDF(o, Miller(1, 0, 0, cs),'Property',tauR_111_fcc_max,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize',5) %'grid','grid_res',10*degree);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0:0.1:0.5])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text

figure
plotPDF(o, Miller(1, 0, 0, cs),'Property',tauR_111_fcc_min,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize',4) %'grid','grid_res',10*degree);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0:0.1:0.5])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text

figure
plotPDF(o, Miller(1, 0, 0, cs),'Property',tauR_111_fcc_mean,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize',4) %'grid','grid_res',10*degree);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0:0.1:0.5])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text

figure
plotPDF(o, Miller(1, 0, 0, cs),'Property',tauR_111_fcc_stdev,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize',4) %'grid','grid_res',10*degree);
cMap = jet_white;
colormap(cMap)
colorbar;
%caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0:0.1:0.5])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text
