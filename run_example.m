%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GUIDED EXAMPLE FOR HOPF MODEL WITH VS WITHOUT COMPETITIVE INTERACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB code provides an example of running the Hopf whole-brain
% computational model of brain activity to obtain the generative effective
% connectivity (GEC) from a structural connectome and functional MRI BOLD
% signals, for one Human Connectome Project example subject.
%
% The main code is fcn_Hopf_NonLinear_cooperative_competitive.m
% the /utils directory contains all required custom MATLAB functions
%
% A working MATLAB installation and valid MATLAB license is required.
% This code has been tested on MATLAB R2019a running on Windows 11 with 8GB RAM;
% and on MATLAB R2021b and R2023b, both running on macOS Sonoma 14 with 32GB RAM
% No additional software installation is required.
% No non-standard hardware is required.
%
% Required inputs/data are:
% FMRI: N-by-T matrix of BOLD signals with N regions and T timepoints;
% SC: N-by-N structural connectome;
% TR: the sampling rate of the BOLD data (in seconds);
% filter_low, filter_high: the bandpass filter frequencies (in Hz);
% Examples of each are provided with this demo in the /data directory.
%
% To run this demo, ensure that you are located in the code's directory.
% Then, simply type 'run_example' in the MATLAB command prompt, and press
% 'enter'.
%
% The total run-time for this example was ~1.5h when run on the Macbook
% devices, and ~4h when run on the Windows laptop.
% Note that run-time will increase if other processes are running at the
% same time.
% Run-time will also vary for different input subjects, and will increase
% with the number of brain regions N and with the  density and the SC matrix.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set-up
% 1) Ensure the required code is on the MATLAB path
addpath(genpath('utils'))
%sub_list = ["S02","S03","S04","S05","S06","S07","S08","S09","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20"];
sub_list = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20"];
% sub_list = ["S07","S08","S09","S10","S11","S12","S13","S14","S15"];
conditions = ["rivalry","replay"];

% Loop over subjects
for i = 1:length(sub_list)
    sub = sub_list(i)
    % Loop over conditions
    for j =1:length(conditions)
        cond = conditions(j)
        % 2) Load example SC and BOLD timeseries from one HCP subject
        % load('data/SC.mat', 'SC');
        % load('data/FMRI.mat', 'FMRI');
        sub_fmri_path = "D:/PhD/BinocularRivalry/"+sub+"/MRI/fMRI_forMatlab/"+cond+"_concat_"+sub+".mat";
        load(sub_fmri_path, 'FMRI'); % Load the current subject
        sc_path = "D:/PhD/BinocularRivalry/Struct/av_sc_100subs_mnf_HCP.mat";
        load(sc_path, 'SC'); % Load the current subject
    
    
        
        % 3) set parameters of FMRI acquisition and bandpass filtering
        % TR = 0.72; % TR of AL example
        TR = 2.2; % TR of br data
        filter_low = 0.008;
        filter_high = 0.09;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Now run the two Hopf models to simulate FC from SC 
        
        % 1) with only *positive connections*
        tic
        flag_cooperativeOnly_YN = true
        [GEC_coop, FC_corr_coop, FCsim_coop, FCemp_coop] = fcn_Hopf_NonLinear_cooperative_competitive(FMRI,...
            SC, TR, filter_low, filter_high, flag_cooperativeOnly_YN);
        toc
        % expected run-time on personal laptop: ~20-30 min
        
        
        % 2) with both cooperative *and competitive connections*
        tic
        flag_cooperativeOnly_YN = false
        [GEC_compet, FC_corr_compet, FCsim_compet, FCemp_cmpet] = fcn_Hopf_NonLinear_cooperative_competitive(FMRI,...
            SC, TR, filter_low, filter_high, flag_cooperativeOnly_YN);
        toc
        % expected run-time on personal laptop: ~1-3h
        
        
        %% Plot results
        %figure; imagesc(GEC_coop); colormap(jet); colorbar
        title({['GEC with only cooperative interactions']; ['fit with empirical FC = ', num2str(FC_corr_coop)]})
        
        %figure; imagesc(GEC_compet); colormap(jet); colorbar
        title({['GEC with competitive interactions']; ['fit with empirical FC = ', num2str(FC_corr_compet)]})
        
        filename = sprintf('output/GEC_coop_%s_sub_%s.mat', cond, sub);
        save (sprintf(filename),'GEC_coop', 'FC_corr_coop', 'FCsim_coop', 'FCemp_coop');
        filename = sprintf('output/GEC_compet_%s_sub_%s.mat', cond, sub);
        save (sprintf(filename),'GEC_compet', 'FC_corr_compet', 'FCsim_compet', 'FCemp_cmpet');
    
    
        %% Expected output: 
        % the matrix GEC_coop should have no negative entries, and
        % FC_corr_coop (correlation between real and simulated FC) should be ~0.36.
        % the matrix GEC_compet should have positive as well as negative entries, and
        % FC_corr_compet (correlation between real and simulated FC) should be ~0.68.
    end
end
    