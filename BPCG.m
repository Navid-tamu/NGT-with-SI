function [SP_BPCG,FNR_BPCG,FPR_BPCG] = BPCG(n,p,M,rho,q,B,A_1,Y_1,ind_1,T_BPCG,tau)
l_tau = length(tau);
FNR_BPCG=zeros(T_BPCG,l_tau);
FPR_BPCG=zeros(T_BPCG,l_tau);
SP_BPCG = zeros(1,l_tau);
k = length(ind_1);
E=cell(1,n);
for i=1:n
    E{i}=find(B(i,:)); % the set of factor nodes connected to individual i at time 0
end
F=cell(1,n);
for i=1:n
    F{i}=find(B(:,i)); % the set of individuals at time 0 connected to factor node i
end

C = cell(1,n);
for i = 1:n
    C{i} = find(A_1(:,i)); % the set of tests where in individual i is involved
end
D = cell(1,M);
for i = 1:M
    D{i} = find(A_1(i,:)); % the set of individuals involved in test i
end

Gamma_ind = zeros(n,n,2);% Matrix Gamma contains the messages from individuls at time 0 to factor nodes (rows:FN, columns:individuals)
Gamma_ind(:,:,1) = 1-p;
Gamma_ind(:,:,2) = p;
Delta_FN = zeros(n,1,2);% Matrix Delta_FN contains the messages from factor nodes to individuls at time 1
Mu_ind = zeros(M,n,2);% Matrix Mu_ind contains the messages from individuls at time 1 to tests
Mu_tst = zeros(M,n,2);  % Matrix Mu_tst contains the messages from tests to individuls at time 1
for i = 1:M
    Mu_tst(i,D{i},:) = 0.5;
end
Delta_ind = zeros(1,n,2);% Matrix Delta_ind contains the messages from individuls at time 1 to factor nodes
Gamma_FN = zeros(n,n,2);% Matrix W contains the messages from factor nodes to items at time 0 (FN:rows, items:column)

for w = 1:T_BPCG
    for i = 1:n % this for loop computes the messages from factor nodes to items at time 1
        l = length(F{i});
        temp = 1;
        for j = 1:l
            temp = temp*(1 - q*Gamma_ind(i,F{i}(j),2));
        end
        temp1 = 0;
        for j = 1:l
            temp2 = 0;
            jj = nchoosek(l,j);
            mtrx = nchoosek(1:l,j);
            for z = 1:jj
                temp2 = temp2 + prod(Gamma_ind(i,F{i}(mtrx(z,:)),2));
            end
            temp1 = temp1 + temp2*(-q)^j;
        end
        norm = Gamma_ind(i,i,1)*temp + Gamma_ind(i,i,2)-Gamma_ind(i,i,1)*temp1;
        Delta_FN(i,1,1) = (Gamma_ind(i,i,1)*temp) / norm;
        Delta_FN(i,1,2) = (Gamma_ind(i,i,2) - Gamma_ind(i,i,1)*temp1) / norm;
    end
    
    for i = 1:n % this for loop computes the messages from items at time 1 to tests
        l = length(C{i});
        for t = 1:l
            temp1 = 1;
            temp2 = 1;
            xi = setdiff(C{i},C{i}(t));
            for j = 1:l-1
                temp1 = temp1*Mu_tst(xi(j),i,1);
                temp2 = temp2*Mu_tst(xi(j),i,2);
            end
            norm = Delta_FN(i,1,1)*temp1 + Delta_FN(i,1,2)*temp2;
            Mu_ind(C{i}(t),i,1) = Delta_FN(i,1,1)*temp1 / norm;
            Mu_ind(C{i}(t),i,2) = Delta_FN(i,1,2)*temp2/ norm;
        end
    end
    
    for t = 1:M % this for loop computes the messages from tests to items at time 1
        l = length(D{t});
        if Y_1(t) == 0
            for i = 1:l
                temp = 1;
                xt = setdiff(D{t},D{t}(i));
                for j = 1:l-1
                    temp = temp*Mu_ind(t,xt(j),1);
                end
                norm = rho + (1 - 2*rho)*temp + rho;
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
                Mu_tst(t,D{t}(i),2) = (1 - rho) / norm;
            end
        end
    end
    
    for i = 1:n % this for loop computes the messages from items at time 1 to factor nodes
        l = length(C{i});
        temp1 = 1;
        temp2 = 1;
        for t = 1:l
            temp1 = temp1*Mu_tst(C{i}(t),i,1);
            temp2 = temp2*Mu_tst(C{i}(t),i,2);
        end
        norm = temp1 + temp2;
        Delta_ind(1,i,1) = temp1 / norm;
        Delta_ind(1,i,2) = temp2 / norm;
    end
    
    
    for i = 1:n % this for loop computes the messages from factor nodes to items at time 0
        l = length(F{i});
        for j = 1:l
            temp = Delta_ind(1,i,1)*Gamma_ind(i,i,1);
            xt = setdiff(F{i},F{i}(j));
            temp1 = 0;
            for d = 1:l-1
                temp = temp*(1 - q*Gamma_ind(i,xt(d),2));
                temp2 = 0;
                jj = nchoosek(l-1,d);
                mtrx = nchoosek(1:l-1,d);
                for z = 1:jj
                    temp2 = temp2 + prod(Gamma_ind(i,xt(mtrx(z,:)),2));
                end
                temp1 = temp1 + temp2*(-q)^d;
            end
            norm = Delta_ind(1,i,2)*Gamma_ind(i,i,2) + temp - Delta_ind(1,i,2)*Gamma_ind(i,i,1)*temp1...
                + Delta_ind(1,i,2)*Gamma_ind(i,i,2) + (1-q)*temp + Delta_ind(1,i,2)*Gamma_ind(i,i,1)*(q-(1-q)*temp1);
            Gamma_FN(i,F{i}(j),1) = (Delta_ind(1,i,2)*Gamma_ind(i,i,2) + temp - Delta_ind(1,i,2)*Gamma_ind(i,i,1)*temp1) / norm;
            Gamma_FN(i,F{i}(j),2) = (Delta_ind(1,i,2)*Gamma_ind(i,i,2) + (1-q)*temp + Delta_ind(1,i,2)*Gamma_ind(i,i,1)*(q-(1-q)*temp1)) / norm;
        end
        temp = Delta_ind(1,i,1);
        temp1 = 0;
        for j = 1:l
            temp = temp*(1 - q*Gamma_ind(i,F{i}(j),2));
            temp2 = 0;
            jj = nchoosek(l,j);
            mtrx = nchoosek(1:l,j);
            for z = 1:jj
                temp2 = temp2 + prod(Gamma_ind(i,F{i}(mtrx(z,:)),2));
            end
            temp1 =temp1 + temp2*(-q)^j;
        end
        norm = temp - Delta_ind(1,i,2)*temp1 + Delta_ind(1,i,2);
        Gamma_FN(i,i,1) = (temp - Delta_ind(1,i,2)*temp1) / norm;
        Gamma_FN(i,i,2) = Delta_ind(1,i,2) / norm;
    end
    
    for i = 1:n
        E{i} = union(E{i},i);
        l = length(E{i});
        for j = 1:l
            ft = setdiff(E{i},E{i}(j));
            norm = (1-p)*prod(Gamma_FN(ft,i,1)) + p*prod(Gamma_FN(ft,i,2));
            Gamma_ind(E{i}(j),i,1) = ((1-p)*prod(Gamma_FN(ft,i,1))) / norm;
            Gamma_ind(E{i}(j),i,2) = (p*prod(Gamma_FN(ft,i,2))) / norm;
        end
    end
    
    X_hat = zeros(1,n);
    for i = 1:n % this for loop finds the probability of each item being defective
        temp = log(Delta_FN(i,1,2) / Delta_FN(i,1,1));
        l = length(C{i});
        for j = 1:l
            temp = temp + log(Mu_tst(C{i}(j),i,2) / Mu_tst(C{i}(j),i,1));
        end
        X_hat(i) = temp;
    end
    
    j = 1;
    for thresh = tau
        Ind = find( X_hat > thresh);
        FNR_BPCG(w,j) = length(setdiff(ind_1,Ind));
        FPR_BPCG(w,j) = length(setdiff(Ind,ind_1));
        j=j+1;
    end
    
end

X_hat = zeros(1,n); % represents an estmation of the status of individuals at time 1
for i = 1:n % % this for loop computes the Log-Likelihood Ratio (LLR) of the marginal probabilities
    temp = log(Delta_FN(i,1,2) / Delta_FN(i,1,1));
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
            SP_BPCG(i) = 1;
        else
            SP_BPCG(i) = 0;
        end
    else
        if length(Ind) == k
            if sort(ind_1) == sort(Ind)
                SP_BPCG(i) = 1;
            end
        end
    end
    i=i+1;
end



