clc
clear
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
cd 'directory for the data based on season and location'
p1 = csvread('GB.results.csv',1,1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
p2 = csvread('WGB.results.csv',1,1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
T = csvread('WT.csv',1,1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
pa = csvread('Actual.csv',1,1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
cd 'directory for the temperature data'

minT = csvread('minT.csv');
maxT = csvread('maxT.csv');

copulaT = csvread('copula.csv');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
maximum.p1 = max(p1);
maximum.p2 = max(p2);
maximum.pa = max(pa);
%--------------------------------------------------------------------------
maximum.p = max(maximum.p1,maximum.p2);
maximum.p = max(maximum.p,maximum.pa)+2;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
minimum.p1 = min(p1);
minimum.p2 = min(p2);
minimum.pa = min(pa);
%--------------------------------------------------------------------------
minimum.p = min(minimum.p1,minimum.p2);
minimum.p = min(minimum.p,minimum.pa);
minimum.p = 0;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
p1 = (p1-minimum.p)/(maximum.p-minimum.p);
p2 = (p2-minimum.p)/(maximum.p-minimum.p);
pa = (pa-minimum.p)/(maximum.p-minimum.p);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plot(p1,'r');
% hold on
% plot(p2,'b');
% plot(pa,'g');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Meshsize = 100;
sigma = 5;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
r(1:size(p1)) = 0;
probability.p1(1:size(p1),1:Meshsize) = 0;
for i=1:size(p1)
    r(i) = round(p1(i)*Meshsize);
    for j=1:Meshsize
        probability.p1(i,j) = exp(-(((r(i)-j)/Meshsize*100)^2/sigma^2));
    end
end

summ.p1 = sum(probability.p1,2);

probability.p1 = probability.p1./summ.p1;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
r(1:size(p2)) = 0;
probability.p2(1:size(p2),1:Meshsize) = 0;
for i=1:size(p2)
    r(i) = round(p2(i)*Meshsize);
    for j=1:Meshsize
        probability.p2(i,j) = exp(-(((r(i)-j)/Meshsize*100)^2/sigma^2));
    end
end

summ.p1 = sum(probability.p2,2);

probability.p2 = probability.p2./summ.p1;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
T = (T-minT)/(maxT-minT);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
probability.p1updated(1:size(p1),1:Meshsize) = 0;
for i=1:size(p1)
    r(i) = round(T(i)*Meshsize);
    if (r(i)==0)
        r(i) = 1;
    end
    for j=1:Meshsize
        probability.p1updated(i,j) = probability.p1(i,j).*copulaT(j,r(i)).^4.8;
    end
    
end

summ.p1updated = sum(probability.p1updated,2);

probability.p1updated = probability.p1updated./summ.p1updated;

%probability.p1 = probability.p1updated;

surfl(probability.p2);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% probability.p1updated(1:size(p1),1:Meshsize) = 0;
% for i=1:size(p1)
%     r(i) = round(W(i)*Meshsize);
%     if (r(i)==0)
%         r(i) = 1;
%     end
%     for j=1:Meshsize
%         probability.p1updated(i,j) = probability.p1(i,j).*copula(j,r(i));
%     end
%     
% end
% 
% summ.p1updated = sum(probability.p1updated,2);
% 
% probability.p1updated = probability.p1updated./summ.p1updated;
% 
% probability.p1 = probability.p1updated;
% surfl(probability.p2);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
NP = 9; % Number of quantiles
P1.p(1:NP,1:size(p1)) = 0;
P1.done = 0;
P1.sum = 0;
for i=1:size(p1)
    P1.sum = 0;
    P1.done = 0;
    for j=1:Meshsize
        P1.sum = P1.sum+probability.p1updated(i,j);
        for k=1:NP
            if(P1.sum>(1/(NP+1))*(P1.done+1) && P1.done<k)
                P1.p(k,i) = j;
                P1.done=P1.done +1;
            end
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
NP = 9;
P1.pnot(1:NP,1:size(p1)) = 0;
P1.done = 0;
P1.sum = 0;
for i=1:size(p1)
    P1.sum = 0;
    P1.done = 0;
    for j=1:Meshsize
        P1.sum = P1.sum+probability.p1(i,j);
        for k=1:NP
            if(P1.sum>(1/(NP+1))*(P1.done+1) && P1.done<k)
                P1.pnot(k,i) = j;
                P1.done=P1.done +1;
            end
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
P2.p(1:NP,1:size(p2)) = 0;
P2.done = 0;
P2.sum = 0;
for i=1:size(p2)
    P2.sum = 0;
    P2.done = 0;
    for j=1:Meshsize
        P2.sum = P2.sum+probability.p2(i,j);
        for k=1:NP
            if(P2.sum>(1/(NP+1))*(P2.done+1) && P2.done<k)
                P2.p(k,i) = j;
                P2.done=P2.done +1;
            end
        end
    end
end
P2.p(1:NP,1:29) = 0;
P2.p(1:NP,size(p2)-26:size(p2)) = 0;
for i=1:size(p1)
    if(pa(i)==0)
        P1.p(1:NP,i)=0;
        P1.pnot(1:NP,i)=0;
        P2.p(1:NP,i)=0;
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
plot(P1.p(1,:)/100,'b');
hold on
plot(P1.pnot(1,:)/100,'r');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plot(P2.p(1,:)/100,'-r');
% plot(P2.p(NP,:)/100,'-r');
plot(pa,'.k');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Pinball.p1(1:NP,1:size(p1)) = 0;
for i=1:NP
    for j=1:size(p1)
        if(P1.p(i,j)/100>pa(j))
            Pinball.p1(i,j)=i/10*(P1.p(i,j)/100-pa(j));
        else
            Pinball.p1(i,j)=-(10-i)/10*(P1.p(i,j)/100-pa(j));
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Pinball.p1not(1:NP,1:size(p1)) = 0;
for i=1:NP
    for j=1:size(p1)
        if(P1.pnot(i,j)/100>pa(j))
            Pinball.p1not(i,j)=i/10*(P1.pnot(i,j)/100-pa(j));
        else
            Pinball.p1not(i,j)=-(10-i)/10*(P1.pnot(i,j)/100-pa(j));
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Pinball.p2(1:NP,1:size(p2)) = 0;
for i=1:NP
    for j=1:size(p2)
        if(P2.p(i,j)/100>pa(j))
            Pinball.p2(i,j)=i/10*(P2.p(i,j)/100-pa(j));
        else
            Pinball.p2(i,j)=-(10-i)/10*(P2.p(i,j)/100-pa(j));
        end
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Pinball.P1 = sum(sum(Pinball.p1))/NP/7;
Pinball.P1not = sum(sum(Pinball.p1not))/NP/7;
Pinball.P2 = sum(sum(Pinball.p2))/NP/7;


% hold off
% plot(p1,'b');
% hold on
% plot(pa,'r');
% plot(p2,'g');
   
            

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
probability.p1weighted(1:size(p1)/7,1:Meshsize) = 0;
selecttau = 1;
Pinball.P1weighted = Pinball.P1;
for tau=0.5:0.1:3
    for i=1:size(p1)/7
        r(i) = round(T(i)*Meshsize);
        if (r(i)==0)
            r(i) = 1;
        end
        for j=1:Meshsize
            probability.p1weighted(i,j) = probability.p1(i,j).*copulaT(j,r(i)).^tau;
        end

    end
    
    summ.p1weighted = sum(probability.p1weighted,2);

    probability.p1weighted = probability.p1weighted./summ.p1weighted;
    
    
    NP = 9; 
    P1.pw(1:NP,1:size(p1)/7) = 0;
    P1.done = 0;
    P1.sum = 0;
    for i=1:size(p1)/7
        P1.sum = 0;
        P1.done = 0;
        for j=1:Meshsize
            P1.sum = P1.sum+probability.p1weighted(i,j);
            for k=1:NP
                if(P1.sum>(1/(NP+1))*(P1.done+1) && P1.done<k)
                    P1.pw(k,i) = j;
                    P1.done = P1.done +1;
                end
            end
        end
    end
    
    
    for i=1:size(p1)/7
    if(pa(i)==0)
        P1.pw(1:NP,i) = 0;
    end
    end
    
    
    
    
    Pinball.p1w(1:NP,1:size(p1)/7) = 0;
    for i=1:NP
        for j=1:size(p1)/7
            if(P1.pw(i,j)/100>pa(j))
                Pinball.p1w(i,j)=i/10*(P1.pw(i,j)/100-pa(j));
            else
                Pinball.p1w(i,j)=-(10-i)/10*(P1.pw(i,j)/100-pa(j));
            end
        end
    end
    Pinball.P1w = sum(sum(Pinball.p1w))/NP/1;
    if(Pinball.P1w <= Pinball.P1weighted)
        Pinball.P1weighted = Pinball.P1w;
        selecttau = tau;
    end
end






tau = selecttau;

    for i=1:size(p1)
        r(i) = round(T(i)*Meshsize);
        if (r(i)==0)
            r(i) = 1;
        end
        for j=1:Meshsize
            probability.p1weighted(i,j) = probability.p1(i,j).*copulaT(j,r(i)).^tau;
        end

    end
    
    summ.p1weighted = sum(probability.p1weighted,2);

    probability.p1weighted = probability.p1weighted./summ.p1weighted;
    
    
    NP = 9;
    P1.pw(1:NP,1:size(p1)) = 0;
    P1.done = 0;
    P1.sum = 0;
    for i=1:size(p1)
        P1.sum = 0;
        P1.done = 0;
        for j=1:Meshsize
            P1.sum = P1.sum+probability.p1weighted(i,j);
            for k=1:NP
                if(P1.sum>(1/(NP+1))*(P1.done+1) && P1.done<k)
                    P1.pw(k,i) = j;
                    P1.done = P1.done +1;
                end
            end
        end
    end
    
    
    for i=1:size(p1)
    if(pa(i)==0)
        P1.pw(1:NP,i) = 0;
    end
    end
    
    
    
    
    Pinball.p1w(1:NP,1:size(p1)) = 0;
    for i=1:NP
        for j=1:size(p1)
            if(P1.pw(i,j)/100>pa(j))
                Pinball.p1w(i,j)=i/10*(P1.pw(i,j)/100-pa(j));
            else
                Pinball.p1w(i,j)=-(10-i)/10*(P1.pw(i,j)/100-pa(j));
            end
        end
    end
    Pinball.P1w = sum(sum(Pinball.p1w))/NP/7;
    
