function[GEC, FC_corr, FCsim, FCemp] = fcn_Hopf_NonLinear_cooperative_competitive(TS, SC, TR, filter_low, filter_high, flag_cooperativeOnly_YN)
% INPUTS:
% TS is N-by-T matrix of BOLD signals with N regions and T timepoints
% SC is N-by-N structural connectome
% TR is the sampling rate of the BOLD data (in seconds)
% filter_low and filter_high are the ends of the bandpass filter (in Hz)
% flag_cooperativeOnly_YN: if true, will set any negatively-valued GEC edges to zero
% (default: true, as per usual Hopf) 
%
% OUTPUTS:
% GEC (N-by-N matrix): matrix of generative effective connectivity, sparse, possibly asymmetric 
% FC_corr (scalar): correlation between real and final simulated FC
% FCsim (N-by-N matrix): final simulated FC
%
% Made by Andrea Luppi;
% Original code ofr cooperative-only Hopf model is from Gustavo Deco and Yonatan Sanz Perl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP PARAMS

if not(exist('flag_cooperativeOnly_YN', 'var'))
    disp('Using the traditional cooperative-only Hopf model')
    flag_cooperativeOnly_YN = true;
end

%Set parameters
Tau=1;
sigma=0.01;
learningRateFromFC = 0.0002;
learningRateFromFCtau=0.00004;
maxC=0.1;

% Other parameters
n_iters = 3000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, T] = size(TS);

    SC = SC;
    assert(size(SC,1)==size(SC,2))
assert(size(SC,1)==N)

%Ensure plausible SC range
SC = 0.2 .* SC./max(SC(:));

mask_triu = find(tril(ones(N),-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical info from the data

%Frequencies
[regionalFrequencies] = fcn_extract_frequencies({TS}, TR, filter_low, filter_high);

%FC and lagged FC (COVtau)
clear signal_filt;
for seed=1:N
    TS(seed,:)=detrend(TS(seed,:)-nanmean(TS(seed,:)));
end

% FC(0-lag)
FCemp=corrcoef(TS');
COVemp=cov(TS');

% COV(tau)
TS_transpose=TS';
for i=1:N
    for j=1:N
        sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
        [clag lags] = xcov(TS_transpose(:,i),TS_transpose(:,j),Tau);
        indx=find(lags==Tau);
        FClag_emp(i,j)=clag(indx)/size(TS_transpose,1);
    end
end
FClag_emp=FClag_emp.*sigratio;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise
GEC= SC;
olderror=100000;

for iter=1:n_iters

    if mod(iter, 100) == 0
        disp(['iter = ', num2str(iter)])
    end


    %% Simulated FC and lagged FC
    [sim_TS] = fcn_Hopf_simulate_BOLD_from_GEC(GEC,regionalFrequencies,sigma, T, TR);

    % FC
    FCsim=corrcoef(sim_TS');
    COVsim=cov(sim_TS');

    % lag-FC
    TS_sim_transpose=sim_TS';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            [clag lags] = xcov(TS_sim_transpose(:,i),TS_sim_transpose(:,j),Tau);
            indx=find(lags==Tau);
            FClagSim(i,j)=clag(indx)/size(TS_sim_transpose,1);
        end
    end
    FClagSim=FClagSim.*sigratio;

   if mod(iter,100)<0.1

        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((FClag_emp-FClagSim).^2));
        disp(errornow)

        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
   end
    %% Optimise the GEC to match target matrix
    for i=1:N  
        for j=1:N
            if SC(i,j)~=0

                %FC error
                GEC(i,j)=GEC(i,j)+learningRateFromFC*(FCemp(i,j)-FCsim(i,j));

                %FC-lag error
                GEC(i,j)=GEC(i,j)+learningRateFromFCtau*(FClag_emp(i,j)-FClagSim(i,j));

                %If only positives allowed, any negative is set to zero
                if flag_cooperativeOnly_YN==true
                    if GEC(i,j)<0
                        GEC(i,j)=0;
                    end
                end

            end
        end
    end
    GEC = GEC/max(max(GEC))*maxC;
end

FC_corr = corr(FCsim(mask_triu), FCemp(mask_triu));

