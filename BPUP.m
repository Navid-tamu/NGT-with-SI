function [SP_BPUP,FNR_BPUP,FPR_BPUP]=BPUP(n,P_1,M,rho,A_1,Y_1,ind_1,T_BPUP,tau)
l_tau = length(tau);
FNR_BPUP=zeros(T_BPUP,l_tau);
FPR_BPUP=zeros(T_BPUP,l_tau);
SP_BPUP = zeros(1,l_tau);
k = length(ind_1);
C = cell(1,n);
for i = 1:n
    C{i} = find(A_1(:,i)); % the set of tests where in individual i is involved
end
D = cell(1,M);
for i = 1:M
    D{i} = find(A_1(i,:)); % the set of individuals involved in test i
end
Mu_ind = zeros(M,n,2); % Matrix Mu_ind contains the messages from individuals to tests
for i=1:n
    Mu_ind(C{i},i,1) = 1-P_1(i);% prior distribution; prob(individual_i=0)
    Mu_ind(C{i},i,2) = P_1(i); % prior distribution; prob(individual_i=1)
end
Mu_tst = zeros(M,n,2); % Matrix H contains the messages from tests to items

for w = 1:T_BPUP
    for t = 1:M % this for loop computes the messages from tests to individuals
        l = length(D{t});
        if Y_1(t) == 0
            for i = 1:l
                temp = 1;
                xt = setdiff(D{t},D{t}(i));
                for j = 1:l-1
                    temp=temp*Mu_ind(t,xt(j),1);
                end
                norm = rho+(1-2*rho)*temp + rho;
                Mu_tst(t,D{t}(i),1) = (rho + (1 - 2*rho)*temp) / norm;
                Mu_tst(t,D{t}(i),2) = rho / norm;
            end
        else
            for i = 1:l
                temp = 1;
                xt = setdiff(D{t},D{t}(i));
                for j = 1:l-1
                    temp = temp*Mu_ind(t,xt(j),1);
                end
                norm = (1 - rho) - (1 - 2*rho)*temp + (1 - rho);
                Mu_tst(t,D{t}(i),1) = ((1 - rho) - (1 - 2*rho)*temp) / norm;
                Mu_tst(t,D{t}(i),2) = (1-rho) / norm;
            end
        end
    end
    
    for i = 1:n % this for loop computes the messages from individuals to tests
        l = length(C{i});
        for t = 1:l
            temp1 = 1;
            temp2 = 1;
            xi = setdiff(C{i},C{i}(t));
            for j = 1:l-1
                temp1 = temp1*Mu_tst(xi(j),i,1);
                temp2 = temp2*Mu_tst(xi(j),i,2);
            end
            norm = (1-P_1(i))*temp1 + P_1(i)*temp2;
            Mu_ind(C{i}(t),i,1) = ((1-P_1(i))*temp1) / norm;
            Mu_ind(C{i}(t),i,2) = P_1(i)*temp2 / norm;
        end
    end
    
    X_hat = zeros(1,n);
    for i = 1:n % this for loop finds the probability of each item being defective
        temp = log(P_1(i) / (1-P_1(i)));
        l = length(C{i});
        for j = 1:l
            temp = temp + log(Mu_tst(C{i}(j),i,2) / Mu_tst(C{i}(j),i,1));
        end
        X_hat(i) = temp;
    end
    
    j = 1;
    for thresh = tau
        Ind = find( X_hat > thresh);
        FNR_BPUP(w,j) = length(setdiff(ind_1,Ind));
        FPR_BPUP(w,j) = length(setdiff(Ind,ind_1));
        j=j+1;
    end
    
end

X_hat = zeros(1,n); % represents an estmation of the status of individuals at time 1
for i = 1:n % this for loop computes the Log-Likelihood Ratio (LLR) of the marginal probabilities
    temp = log(P_1(i) / (1-P_1(i)));
    l = length(C{i});
    for j = 1:l
        temp = temp + log(Mu_tst(C{i}(j),i,2) / Mu_tst(C{i}(j),i,1));
    end
    X_hat(i) = temp;
end

i=1;
for thresh = tau
    Ind = find( X_hat > thresh); % indeces of infected individuals
    if k == 0
        if isempty(Ind)
            SP_BPUP(i) = 1;
        else
            SP_BPUP(i) = 0;
        end
    else
        if length(Ind) == k
            if sort(ind_1) == sort(Ind)
                SP_BPUP(i) = 1;
            end
        end
    end
    i=i+1;
end


