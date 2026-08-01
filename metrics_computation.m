clear all;

N_sub = 20;
N_sub = 10;
% sub_list = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20"];
sub_list = ["S11","S12","S13","S14","S15"];
conditions = ["rivalry", "replay"];
coop_or_comp = ["coop","compet"];


ec_folder='output/';

for n_cond=1:2
    for c_or_c=1:2
        for nsub=1:N_sub
            sub=sub_list{nsub};
            filename=ec_folder+"GEC_"+coop_or_comp(c_or_c)+"_"+conditions(n_cond)+"_sub_"+sub+".mat";
            S=load(filename);
            % Find the field that starts with 'GEC'
            vars = fieldnames(S);
            isGEC = startsWith(vars, 'GEC');
            varname = vars{isGEC};
            
            % Extract it
            Ceff = S.(varname);      % your connectivity matrix
            A=Ceff';
            d=sum(A)'; %income connections
            delta=sum(A,2); %outcome connections
            u=d+delta;
            v=d-delta;
            Lambda=diag(u)-A-A';
            Lambda(1,1)=0;
            gamma=linsolve(Lambda,v);
            gamma=gamma-min(gamma);
            hierarchicallevels(nsub,:)=gamma';
            H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
            F0=sum(sum((A.*H)))/sum(sum(A));
            trophiccoherence(nsub) = 1-F0; 
            imbalance(nsub, :) = v';
        end
        save (sprintf(strcat('output/GEC_metrics_',coop_or_comp(c_or_c),'_',conditions(n_cond),'.mat')),'trophiccoherence','hierarchicallevels','imbalance');
    end
end