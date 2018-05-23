%This codes plots the orientations from hcp_nrss_materialia on (0001) pole figures colored according to their max, min, mean, and standard deviation nRSS values
%It makes use of the function jet_white for the colormap
%Garrison Hommer, 22MAY2018

figure;
plotPDF(o, Miller(0, 0, 0, 1, cs),'Property',tauR_basal_max,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize', 3.4);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0:0.1:0.5])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text

figure;
plotPDF(o, Miller(0, 0, 0, 1, cs),'Property',tauR_basal_min,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize', 3.4);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0 0.015 0.03 0.045 0.06 0.073])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text

figure;
plotPDF(o, Miller(0, 0, 0, 1, cs),'Property',tauR_basal_mean,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize', 3.4);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.5]);
hcb=colorbar;
%set(hcb,'YTick',[0 0.015 0.03 0.045 0.06 0.073])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text

figure;
plotPDF(o, Miller(0, 0, 0, 1, cs),'Property',tauR_pya_stdev,'points','all','xAxisDirection','east','zAxisDirection','OutPlane','MarkerSize', 3.4);
cMap = jet_white;
colormap(cMap)
colorbar;
caxis([0.0 0.07]);
hcb=colorbar;
%set(hcb,'YTick',[0 0.015 0.03 0.045 0.06 0.073])
% h = findall(gcf,'type','text'); %remove all text
% set(h,'visible','off') %remove all text