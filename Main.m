% Please read the following paper for a better understanding of the code:
% "Noisy Group Testing with Side Information" available on: https://arxiv.org/abs/2202.12284
%%
clc
clear
close all
n = 500; % total number of individuals
p = 0.01; % prevalence rate at time 0
M = 200; % number of tests
tau = -8:0.2:8; % threshold range
l_tau = length(tau);
rho = 0.01; % noise parameter
q = 0.1; % contagion probability: the probability that an infected individual infects a healthy individual
theta = 0.008; % interaction probability: the probability that individual i interacts with individual j
T_BPIP = 20; % number of iterations for the Belief Propagation using Initial Prior probabilities (BPIP)
T_BPUP = 20; % number of iterations for the Belief Propagation using Updated Prior probabilities (BPUP)
T_BPCG = 40; % number of iterations for the Belief Propagation on Combined Graph (BPCG)
v = log(2); % Bernoulli design parameter
A_0 = rand(M,n) <(v/(n*p)); % measurement matrix with no knowledge on the interaction model
n_trial = 100; % number of trials
SP_BPIP = zeros(n_trial,l_tau); % success probability for BPIP
SP_BPUP = zeros(n_trial,l_tau); % success probability for BPUP
SP_BPCG = zeros(n_trial,l_tau); % success probability for BPCG
FNR_BPIP = zeros(T_BPIP,n_trial,l_tau); FPR_BPIP = zeros(T_BPIP,n_trial,l_tau); % FNR and FPR for BPIP
FNR_BPUP = zeros(T_BPUP,n_trial,l_tau); FPR_BPUP = zeros(T_BPUP,n_trial,l_tau); % FNR and FPR for BPUP
FNR_BPCG = zeros(T_BPCG,n_trial,l_tau); FPR_BPCG = zeros(T_BPCG,n_trial,l_tau); % FNR and FPR for BPCG
Def=zeros(1,n_trial); % defective items for each trial
%%
% If you have parallel computing toolbox, you can leave the code as is;
% If not, please change the following "parfor" to "for".
parfor trial = 1:n_trial
    disp(trial)
    B=rand(n,n) < theta; B=B - diag(diag(B));
    B=triu(B); B= B+B'; % A symmetric matrix that captures interactions between individuals
    P_1=zeros(1,n); % represents prior probability of infection for individuals at time 1
    for i=1:n
        P_1(i)=1-(1-p)*(1-p*q)^(sum(B(i,:))); % prior probability that individual i is infected at time 1
    end
    k_1=sum(P_1); % average number of infected individuals at time 1
    A_1=rand(M,n) < (v/k_1); % measurement matrix with the knowledge on the interaction model
    N=rand(M,1)<rho; % noise vector
    X_0=rand(1,n)< p; % represents the status of individuals at time 0
    X_1=zeros(1,n); % represents the status of individuals at time 1
    if sum(X_0)==0
        ind_1=[]; % indeces of infected individuals at time 1
        Y_0=N; % represents test results when there is no knowledge on the interaction model
        Y_1=N; % represents test results when the knowledge on the interaction model is available
    else
        for i=1:n
            if X_0(i)==1
                X_1(i)=1;
            elseif B(i,:)*X_0' > 0
                X_1(i)=rand < 1-(1-q)^(B(i,:)*X_0');
            else
                X_1(i)=0;
            end
        end
        ind_1=(find(X_1)); %index of defective items at time 1
        Y_0=xor(or(sum(A_0(:,ind_1),2),zeros(M,1)),N);
        Y_1=xor(or(sum(A_1(:,ind_1),2),zeros(M,1)),N);
    end
    Def(trial)=length(ind_1);
    [SP_BPIP(trial,:),FNR_BPIP(:,trial,:),FPR_BPIP(:,trial,:)] = BPIP(n,p,M,rho,A_0,Y_0,ind_1,T_BPIP,tau);
    [SP_BPUP(trial,:),FNR_BPUP(:,trial,:),FPR_BPUP(:,trial,:)] = BPUP(n,P_1,M,rho,A_1,Y_1,ind_1,T_BPUP,tau);
    [SP_BPCG(trial,:),FNR_BPCG(:,trial,:),FPR_BPCG(:,trial,:)] = BPCG(n,p,M,rho,q,B,A_1,Y_1,ind_1,T_BPCG,tau);
    
end
SP_BPIP = mean(SP_BPIP); SP_BPUP = mean(SP_BPUP); SP_BPCG = mean(SP_BPCG);

figure(1)
plot(tau,SP_BPIP,'-o','LineWidth',2)
hold on
plot(tau,SP_BPUP,'-*','LineWidth',2)
hold on
plot(tau,SP_BPCG,'->','LineWidth',2)
hold off
ylabel('Success Probability','fontweight','bold','FontSize',14,'interpreter','latex')
xlabel('Threshold ($\tau$)','fontweight','bold','FontSize',14,'interpreter','latex')
legend('BPIP','BPUP','BPCG')
grid on
set(gca,'FontSize',16)

NDef=n-Def;
Def=sum(Def);
NDef=sum(NDef);
% Note that unlike the success probability, FNR and FPR do not converge. 
% Instead, after a certain number of iterations, FNR and FPR oscillate around an average value. 
% Thus, for each value of ?, we compute the average FNR and FPR over a range of iterations.
FNR_BPIP=sum(mean(FNR_BPIP(10:20,:,:))); FNR_BPIP=reshape(FNR_BPIP,1,l_tau); 
FPR_BPIP=sum(mean(FPR_BPIP(10:20,:,:))); FPR_BPIP=reshape(FPR_BPIP,1,l_tau);

FNR_BPUP=sum(mean(FNR_BPUP(10:20,:,:))); FNR_BPUP=reshape(FNR_BPUP,1,l_tau); 
FPR_BPUP=sum(mean(FPR_BPUP(10:20,:,:))); FPR_BPUP=reshape(FPR_BPUP,1,l_tau);

FNR_BPCG=sum(mean(FNR_BPCG(30:40,:,:))); FNR_BPCG=reshape(FNR_BPCG,1,l_tau); 
FPR_BPCG=sum(mean(FPR_BPCG(30:40,:,:))); FPR_BPCG=reshape(FPR_BPCG,1,l_tau);

figure(2)
plot(FPR_BPIP/NDef,FNR_BPIP/Def,'-o','LineWidth',2)
hold on
plot(FPR_BPUP/NDef,FNR_BPUP/Def,'-*','LineWidth',2)
hold on
plot(FPR_BPCG/NDef,FNR_BPCG/Def,'->','LineWidth',2)
hold off
ylabel('False-Negative Rate','fontweight','bold','FontSize',14,'interpreter','latex')
xlabel('False-Positive Rate','fontweight','bold','FontSize',14,'interpreter','latex')
legend('BPIP','BPUP','BPCG')
grid on
set(gca,'FontSize',16)

