%Garrison Hommer, 1APR2021
%This code calculates hcp normalized resolved shear stresses (nRSS) using the coordinate conventions of Staroselsky PhD Thesis 1998: Crystal Plasticicty Due to Slip and Twinning
%Normalizing the resolved shear stress by the difference of the maximum and minimum stress tensor eigenvalues bounds the nRSS values by [0 0.5]
%The codes makes use of the MTEX toolbox
%The codes uses the functions dyad11 and scalar22 for computing the dyadic and scalar products between 2nd order tensors, respectivley
%Each pole orientation is rotated about the c-axis for determining min, max and standard deviation

%% Crystal Parameters

%Ti7Al lattice parameters determined by Hommer et al.
c = 4.6767; %+-0.0005
a = 2.9300; %+_0.0005
%c = 4.5; %+-0.0005
%a = 2.25; %+_0.0005


cs = crystalSymmetry('6/mmm',[2.93,2.93,4.6767],'X||b', 'Y||-a*', 'Z||c'); %symmetry used in Staroselsky PhD Thesis 1998, Ti7Al lattice parameters determined by Hommer et al.
%cs = crystalSymmetry('6/mmm',[2.25,2.25,4.5],'X||b', 'Y||-a*', 'Z||c'); %symmetry used in Staroselsky PhD Thesis 1998, Ti7Al lattice parameters determined by Hommer et al.
ss = specimenSymmetry('1'); %no specimen symmetry
num_orientations = 60000; %number of orientations to sample

%a few arbitrary stress tensors
sig_0x1y = [0 0 0;0 1 0;0 0 0];
sig_025x1y = [0.25 0 0;0 1 0;0 0 0];
sig_neg025x1y = [-0.25 0 0;0 1 0;0 0 0];
sig_neg025xneg1y = [-0.25 0 0;0 -1 0;0 0 0];
sig_05x1y = [0.5 0 0;0 1 0;0 0 0];
sig_neg05x1y = [-0.5 0 0;0 1 0;0 0 0];
sig_075x1y = [0.75 0 0;0 1 0;0 0 0];
sig_1x1y = [1 0 0;0 1 0;0 0 0];
sig_neg1x1y = [-1 0 0;0 1 0;0 0 0];

stress = sig_neg05x1y; %applied stress tensor

%%%%SLIP SYSTEMS%%% Staroselsky PhD Thesis: Crystal Plasticicty Due to Slip and Twinning
denom1 = sqrt(4*c^2 + 3*a^2);
denom2 = sqrt(c^2 + a^2);
numa = sqrt(3)*a;
numc = sqrt(3)*c;

%basal slip systems
ns(:,1) = [0 0 1]; 
ms(:,1) = [1/2 -sqrt(3)/2 0];
ns(:,2) = [ 0 0 1];
ms(:,2) = [1/2 sqrt(3)/2 0];
ns(:,3) = [0 0 1];
ms(:,3) = [-1 0 0];

%prismatic slip systems
ns(:,4) = [0 1 0]; 
ms(:,4) = [1 0 0];
ns(:,5) = [-sqrt(3)/2 1/2 0];
ms(:,5) = [1/2 sqrt(3)/2 0];
ns(:,6) = [-sqrt(3)/2 -1/2 0];
ms(:,6) = [-1/2 sqrt(3)/2 0];

%Pyramidal <a> slip systems
ns(:,7) = [0 -2*c/denom1 numa/denom1];
ms(:,7) = [1 0 0];
ns(:,8) = [numc/denom1 -c/denom1 numa/denom1];
ms(:,8) = [1/2 sqrt(3)/2 0];
ns(:,9) = [numc/denom1 c/denom1 numa/denom1];
ms(:,9) = [-1/2 sqrt(3)/2 0];
ns(:,10) = [0 2*c/denom1 numa/denom1];
ms(:,10) = [-1 0 0];
ns(:,11) = [-numc/denom1 c/denom1 numa/denom1];
ms(:,11) = [-1/2 -sqrt(3)/2 0];
ns(:,12) = [-numc/denom1 -c/denom1 numa/denom1];
ms(:,12) = [1/2 -sqrt(3)/2 0];

%1st order pyramidal <c+a> slip systems
ns(:,13) = [0 -2*c/denom1 numa/denom1];
ms(:,13) = [-a/(2*denom2) numa/(2*denom2) c/denom2];
ns(:,14) = [numc/denom1 -c/denom1 numa/denom1];
ms(:,14) = [-a/denom2 0  c/denom2];
ns(:,15) = [numc/denom1 c/denom1 numa/denom1];
ms(:,15) = [-a/(2*denom2) -numa/(2*denom2) c/denom2];
ns(:,16) = [0 2*c/denom1 numa/denom1];
ms(:,16) = [a/(2*denom2) -numa/(2*denom2) c/denom2];
ns(:,17) = [-numc/denom1 c/denom1 numa/denom1];
ms(:,17) = [a/denom2 0 c/denom2];
ns(:,18) = [-numc/denom1 -c/denom1 numa/denom1];
ms(:,18) = [a/(2*denom2) numa/(2*denom2) c/denom2];
ns(:,19) = [0 -2*c/denom1 numa/denom1];
ms(:,19) = [a/(2*denom2) numa/(2*denom2) c/denom2];
ns(:,20) = [numc/denom1 -c/denom1 numa/denom1];
ms(:,20) = [-a/(2*denom2) numa/(2*denom2) c/denom2];
ns(:,21) = [numc/denom1 c/denom1 numa/denom1];
ms(:,21) = [-a/denom2 0  c/denom2];
ns(:,22) = [0 2*c/denom1 numa/denom1];
ms(:,22) = [-a/(2*denom2) -numa/(2*denom2) c/denom2];
ns(:,23) = [-numc/denom1 c/denom1 numa/denom1];
ms(:,23) = [a/(2*denom2) -numa/(2*denom2) c/denom2];
ns(:,24) = [-numc/denom1 -c/denom1 numa/denom1];
ms(:,24) = [a/denom2 0 c/denom2];

%2nd order pyramidal <c+a> slip systems
ns(:,25) = [c/(2*denom2) -numc/(2*denom2) a/denom2]; 
ms(:,25) = [-a/(2*denom2) numa/(2*denom2) c/denom2];
ns(:,26) = [c/denom2 0 a/denom2]; 
ms(:,26) = [-a/denom2 0  c/denom2];
ns(:,27) = [c/(2*denom2) numc/(2*denom2) a/denom2]; 
ms(:,27) = [-a/(2*denom2) -numa/(2*denom2) c/denom2];
ns(:,28) = [-c/(2*denom2) numc/(2*denom2) a/denom2]; 
ms(:,28) = [a/(2*denom2) -numa/(2*denom2) c/denom2];
ns(:,29) = [-c/denom2 0 a/denom2]; 
ms(:,29) = [a/denom2 0 c/denom2];
ns(:,30) = [-c/(2*denom2) -numc/(2*denom2) a/denom2]; 
ms(:,30) = [a/(2*denom2) numa/(2*denom2) c/denom2];

%% Resolved Shear Stress Calculations
odf = uniformODF(cs,ss); %create uniform orientation distribution function
o = orientation(calcOrientations(odf,num_orientations), cs);
num_orientations = length(o);
rot_steps = [0:1:60]; %degree rotation steps about c-axis
m = Miller(0,0,0,1,cs); %pole used for pole figures

%create arrays
% tauR_basal = zeros(length(rot_steps),1);
% tauR_basal_max = zeros(num_orientations,1);
% tauR_basal_min = zeros(num_orientations,1);
% tauR_basal_mean = zeros(num_orientations,1);
% tauR_basal_stdev = zeros(num_orientations,1);
% 
% tauR_pris = zeros(length(rot_steps),1);
% tauR_pris_max = zeros(num_orientations,1);
% tauR_pris_min = zeros(num_orientations,1);
% tauR_pris_mean = zeros(num_orientations,1);
% tauR_pris_stdev = zeros(num_orientations,1);
% 
% tauR_pya = zeros(length(rot_steps),1);
% tauR_pya_max = zeros(num_orientations,1);
% tauR_pya_min = zeros(num_orientations,1);
% tauR_pya_mean = zeros(num_orientations,1);
% tauR_pya_stdev = zeros(num_orientations,1);
% 
% tauR_py1 = zeros(length(rot_steps),1);
% tauR_py1_max = zeros(num_orientations,1);
% tauR_py1_min = zeros(num_orientations,1);
% tauR_py1_mean = zeros(num_orientations,1);
% tauR_py1_stdev = zeros(num_orientations,1);
% 
% tauR_py2 = zeros(length(rot_steps),1);
% tauR_py2_max = zeros(num_orientations,1);
% tauR_py2_min = zeros(num_orientations,1);
% tauR_py2_mean = zeros(num_orientations,1);
% tauR_py2_stdev = zeros(num_orientations,1);

tauR_basal_pris = zeros(length(rot_steps),1);
tauR_basal_pris_max = zeros(num_orientations,1);
tauR_basal_pris_min = zeros(num_orientations,1);
tauR_basal_pris_mean = zeros(num_orientations,1);
tauR_basal_pris_stdev = zeros(num_orientations,1);

tauR_py_all = zeros(length(rot_steps),1);
tauR_py_all_max = zeros(num_orientations,1);
tauR_py_all_min = zeros(num_orientations,1);
tauR_py_all_mean = zeros(num_orientations,1);
tauR_py_all_stdev = zeros(num_orientations,1);

counter = 1;
for ii = 1:num_orientations
    rot_axis = o(ii)*m; %rotation axis for each orientation from uniformODF
    rot = rotation('axis',rot_axis,'angle',rot_steps*degree); %Euler angles rotation
    o_rot_mat = matrix(rot * o(ii)); %rotated orientation matrix
    for jj = 1:length(rot_steps) %loop through rotations
        nsR = o_rot_mat(:,:,jj) * ns; %rotating slip normal
        msR = o_rot_mat(:,:,jj) * ms; %rotating slip direction
        for ll = 1:30
            R(:,:,ll) = 0.5*(dyad11(nsR(:,ll),msR(:,ll)) + dyad11(nsR(:,ll),msR(:,ll))); %Schmid tensor
            tauR(:,ll) = abs(scalar22(squeeze(R(:,:,ll)),stress))/(max(eig(stress)) - min(eig(stress))); %normalized resolved shear stress
        end
%         tauR_basal(jj) = max(tauR(:,1:3)); %basal slip RSS
%         tauR_pris(jj) = max(tauR(:,4:6)); %prismatic slip RSS
%         tauR_pya(jj) = max(tauR(:,7:12)); %pyramidal <a> slip RSS
%         tauR_py1(jj) = max(tauR(:,13:24)); %pyradmidal <c+a> 1st RSS
%         tauR_py2(jj) = max(tauR(:,25:30)); %pyradmidal <c+a> 2nd RSS
        tauR_basal_pris(jj) = max(tauR(:,1:6)); %basal-prismatic composite slip RSS
        tauR_py_all(jj) = max(tauR(:,7:30)); %basal-prismatic composite slip RSS
    end
%     tauR_basal_max(ii) = max(tauR_basal);
%     tauR_basal_min(ii) = min(tauR_basal);
%     tauR_basal_mean(ii) = mean(tauR_basal);
%     tauR_basal_stdev(ii) = std(tauR_basal);
%     
%     tauR_pris_max(ii) = max(tauR_pris);
%     tauR_pris_min(ii) = min(tauR_pris);
%     tauR_pris_mean(ii) = mean(tauR_pris);
%     tauR_pris_stdev(ii) = std(tauR_pris);
% 
%     tauR_pya_max(ii) = max(tauR_pya);
%     tauR_pya_min(ii) = min(tauR_pya);
%     tauR_pya_mean(ii) = mean(tauR_pya);
%     tauR_pya_stdev(ii) = std(tauR_pya);
%     
%     tauR_py1_max(ii) = max(tauR_py1);
%     tauR_py1_min(ii) = min(tauR_py1);
%     tauR_py1_mean(ii) = mean(tauR_py1);
%     tauR_py1_stdev(ii) = std(tauR_py1);
%     
%     tauR_py2_max(ii) = max(tauR_py2);
%     tauR_py2_min(ii) = min(tauR_py2);
%     tauR_py2_mean(ii) = mean(tauR_py2);
%     tauR_py2_stdev(ii) = std(tauR_py2);
    
    tauR_basal_pris_max(ii) = max(tauR_basal_pris);
    tauR_basal_pris_min(ii) = min(tauR_basal_pris);
    tauR_basal_pris_mean(ii) = mean(tauR_basal_pris);
    tauR_basal_pris_stdev(ii) = std(tauR_basal_pris);
    
    tauR_py_all_max(ii) = max(tauR_py_all);
    tauR_py_all_min(ii) = min(tauR_py_all);
    tauR_py_all_mean(ii) = mean(tauR_py_all);
    tauR_py_all_stdev(ii) = std(tauR_py_all);
end

% tauR_basal = tauR_basal';
% tauR_pris = tauR_pris';
% tauR_pya = tauR_pya';
% tauR_py1 = tauR_py1';
% %tauR_py1_norm = tauR_py1/3; %normalized for CtauR = 3CtauR basal and prismatic
% tauR_py2 = tauR_py2';
% %tauR_py2_norm = tauR_py2/3; %normalized for CtauR = 3CtauR basal and prismatic
% tauR_basal_pris = tauR_basal_pris';