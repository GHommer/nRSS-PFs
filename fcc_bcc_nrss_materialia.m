%Garrison Hommer, 1APR2021
%This code calculates fcc and bcc normalized resolved shear stresses (nRSS)
%Normalizing the resolved shear stress by the difference of the maximum and minimum stress tensor eigenvalues bounds the nRSS values by [0 0.5]
%The codes makes use of the MTEX toolbox
%The codes uses the functions dyad11 and scalar22 for computing the dyadic and scalar products between 2nd order tensors, respectivley
%Each pole orientation is rotated about the PF plotting pole for determining min, max and standard deviation

%% Crystal Parameters

%By specifing no crystal symmetry, rotation about crystallographically equivalent poles (e.g., [100] [010] [001]) is distinct, which allows for correct PF plotting
cs = crystalSymmetry('1');
ss = specimenSymmetry('1');
num_orientations = 3000; %number of orientations to sample

%2-norm = 1 stress tensors
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

%%%%SLIP SYSTEMS%%%

%{111} <110> (<1-10> Callister) fcc slip systems 
ns(:,1) = [1 1 -1]; ms(:,1) = [0 1 1];
ns(:,2) = [1 1 -1]; ms(:,2) = [1 0 1];
ns(:,3) = [1 1 -1]; ms(:,3) = [1 -1 0];
ns(:,4) = [1 -1 -1]; ms(:,4) = [0 1 -1];
ns(:,5) = [1 -1 -1]; ms(:,5) = [1 0 1];
ns(:,6) = [1 -1 -1]; ms(:,6) = [1 1 0];
ns(:,7) = [1 -1 1]; ms(:,7) = [0 1 1];
ns(:,8) = [1 -1 1]; ms(:,8) = [1 0 -1];
ns(:,9) = [1 -1 1]; ms(:,9) = [1 1 0];
ns(:,10) = [1 1 1]; ms(:,10) = [0 1 -1];
ns(:,11) = [1 1 1]; ms(:,11) = [1 0 -1];
ns(:,12) = [1 1 1]; ms(:,12) = [1 -1 0];


%{110} <111> (<-111> Callister) bcc slip systems
ns(:,13) = [0 1 1]; ms(:,13) = [1 1 -1];
ns(:,14) = [1 0 1]; ms(:,14) = [1 1 -1];
ns(:,15) = [1 -1 0]; ms(:,15) = [1 1 -1];
ns(:,16) = [0 1 -1]; ms(:,16) = [1 -1 -1];
ns(:,17) = [1 0 1]; ms(:,17) = [1 -1 -1];
ns(:,18) = [1 1 0]; ms(:,18) = [1 -1 -1];
ns(:,19) = [0 1 1]; ms(:,19) = [1 -1 1];
ns(:,20) = [1 0 -1]; ms(:,20) = [1 -1 1];
ns(:,21) = [1 1 0]; ms(:,21) = [1 -1 1];
ns(:,22) = [0 1 -1]; ms(:,22) = [1 1 1];
ns(:,23) = [1 0 -1]; ms(:,23) = [1 1 1];
ns(:,24) = [1 -1 0]; ms(:,24) = [1 1 1];

%{211} <111> (<-111> Callister) bcc slip systems
ns(:,25) = [-2 1 1]; ms(:,25) = [1 1 1];
ns(:,26) = [1 -2 1]; ms(:,26) = [1 1 1];
ns(:,27) = [1 1 -2]; ms(:,27) = [1 1 1];
ns(:,28) = [2 1 1]; ms(:,28) = [-1 1 1];
ns(:,29) = [-1 -2 1]; ms(:,29) = [-1 1 1];
ns(:,30) = [-1 1 -2]; ms(:,30) = [-1 1 1];
ns(:,31) = [-2 -1 1]; ms(:,31) = [1 -1 1];
ns(:,32) = [1 2 1]; ms(:,32) = [1 -1 1];
ns(:,33) = [1 -1 -2]; ms(:,33) = [1 -1 1];
ns(:,34) = [2 -1 1]; ms(:,34) = [-1 -1 1];
ns(:,35) = [-1 2 1]; ms(:,35) = [-1 -1 1];
ns(:,36) = [-1 -1 -2]; ms(:,36) = [-1 -1 1];

%{321} <111> (<-111> Callister) bcc slip systems
ns(:,37) = [3 -1 -2]; ms(:,37) = [1 1 1];
ns(:,38) = [3 -2 -1]; ms(:,38) = [1 1 1];
ns(:,39) = [-1 3 -2]; ms(:,39) = [1 1 1];
ns(:,40) = [-2 3 -1]; ms(:,40) = [1 1 1];
ns(:,41) = [-1 -2 3]; ms(:,41) = [1 1 1];
ns(:,42) = [-2 -1 3]; ms(:,42) = [1 1 1];
ns(:,43) = [3 1 2]; ms(:,43) = [-1 1 1];
ns(:,44) = [3 2 1]; ms(:,44) = [-1 1 1];
ns(:,45) = [1 3 -2]; ms(:,45) = [-1 1 1];
ns(:,46) = [2 3 -1]; ms(:,46) = [-1 1 1];
ns(:,47) = [1 -2 3]; ms(:,47) = [-1 1 1];
ns(:,48) = [2 -1 3]; ms(:,48) = [-1 1 1];
ns(:,49) = [3 1 -2]; ms(:,49) = [1 -1 1];
ns(:,50) = [3 2 -1]; ms(:,50) = [1 -1 1];
ns(:,51) = [1 3 2]; ms(:,51) = [1 -1 1];
ns(:,52) = [2 3 1]; ms(:,52) = [1 -1 1];
ns(:,53) = [-1 2 3]; ms(:,53) = [1 -1 1];
ns(:,54) = [-2 1 3]; ms(:,54) = [1 -1 1];
ns(:,55) = [3 -1 2]; ms(:,55) = [-1 -1 1];
ns(:,56) = [3 -2 1]; ms(:,56) = [-1 -1 1];
ns(:,57) = [-1 3 2]; ms(:,57) = [-1 -1 1];
ns(:,58) = [-2 3 1]; ms(:,58) = [-1 -1 1];
ns(:,59) = [1 2 3]; ms(:,59) = [-1 -1 1];
ns(:,60) = [2 1 3]; ms(:,60) = [-1 -1 1];

%normalize slip systems to unit vectors
for ii = 1:60
    ns(:,ii) = ns(:,ii)/norm(ns(:,ii));
    ms(:,ii) = ms(:,ii)/norm(ms(:,ii));
end

odf = uniformODF(cs,ss); %create uniform orientation distribution function
o = orientation(calcOrientations(odf,num_orientations), cs);
num_orientations = length(o);
rot_steps = [0:1:90]; %degree rotation steps about c-axis
m100 = Miller(1,0,0,cs); %pole used for pole figures

%create arrays
tauR_111_fcc = zeros(length(rot_steps),1);
tauR_111_fcc_max = zeros(num_orientations,1);
tauR_111_fcc_min = zeros(num_orientations,1);
tauR_111_fcc_mean = zeros(num_orientations,1);
tauR_111_fcc_stdev = zeros(num_orientations,1);

tauR_110_bcc = zeros(length(rot_steps),1);
tauR_110_bcc_max = zeros(num_orientations,1);
tauR_110_bcc_min = zeros(num_orientations,1);
tauR_110_bcc_mean = zeros(num_orientations,1);
tauR_110_bcc_stdev = zeros(num_orientations,1);

tauR_211_bcc = zeros(length(rot_steps),1);
tauR_211_bcc_max = zeros(num_orientations,1);
tauR_211_bcc_min = zeros(num_orientations,1);
tauR_211_bcc_mean = zeros(num_orientations,1);
tauR_211_bcc_stdev = zeros(num_orientations,1);

tauR_321_bcc = zeros(length(rot_steps),1);
tauR_321_bcc_max = zeros(num_orientations,1);
tauR_321_bcc_min = zeros(num_orientations,1);
tauR_321_bcc_mean = zeros(num_orientations,1);
tauR_321_bcc_stdev = zeros(num_orientations,1);

for ii = 1:num_orientations
    rot_axis = o(ii)*m100; %rotation axis for each orientation from uniformODF
    rot = rotation('axis',rot_axis,'angle',rot_steps*degree); %Euler angles rotation
    o_rot_mat = matrix(rot * o(ii)); %rotated orientation matrix
    for jj = 1:length(rot_steps) %loop through rotations
        nsR100 = o_rot_mat(:,:,jj) * ns; %rotating slip normal
        msR100 = o_rot_mat(:,:,jj) * ms; %rotating slip direction
        for ll = 1:60
            R(:,:,ll) = 0.5*(dyad11(nsR100(:,ll),msR100(:,ll)) + dyad11(nsR100(:,ll),msR100(:,ll))); %Schmid tensor
            tauR(:,ll) = abs(scalar22(squeeze(R(:,:,ll)),stress))/(max(eig(stress)) - min(eig(stress))); %normalized resolved shear stress
        end
        tauR_111_fcc(jj) = max(tauR(:,1:12)); %fcc slip rss
        tauR_110_bcc(jj) = max(tauR(:,13:24)); %bcc slip rss
        tauR_211_bcc(jj) = max(tauR(:,25:36)); %bcc slip rss
        tauR_321_bcc(jj) = max(tauR(:,37:60)); %bcc slip rss
    end

    tauR_111_fcc_max(ii) = max(tauR_111_fcc);
    tauR_111_fcc_min(ii) = min(tauR_111_fcc);
    tauR_111_fcc_mean(ii) = mean(tauR_111_fcc);
    tauR_111_fcc_stdev(ii) = std(tauR_111_fcc);

    tauR_110_bcc_max(ii) = max(tauR_110_bcc);
    tauR_110_bcc_min(ii) = min(tauR_110_bcc);
    tauR_110_bcc_mean(ii) = mean(tauR_110_bcc);
    tauR_110_bcc_stdev(ii) = std(tauR_110_bcc);

    tauR_211_bcc_max(ii) = max(tauR_211_bcc);
    tauR_211_bcc_min(ii) = min(tauR_211_bcc);
    tauR_211_bcc_mean(ii) = mean(tauR_211_bcc);
    tauR_211_bcc_stdev(ii) = std(tauR_211_bcc);

    tauR_321_bcc_max(ii) = max(tauR_321_bcc);
    tauR_321_bcc_min(ii) = min(tauR_321_bcc);
    tauR_321_bcc_mean(ii) = mean(tauR_321_bcc);
    tauR_321_bcc_stdev(ii) = std(tauR_321_bcc);
end