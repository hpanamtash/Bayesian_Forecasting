clc
clear
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
cd 'directory for power data'
p = csvread('BayshoreElementary-all-converted.csv'); % Reading power output
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
cd 'directory for weather data'
t = csvread('W-BayshoreElementary-all-converted.csv',0,0); % Reading Weather data

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Datasize = size(p);
TP(1:2, 1:Datasize)= 0;
TP(1, 1:Datasize) = p;
TP(2, 1:Datasize) = t;
TP = TP';
TP = rmmissing(TP);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
MS=0; % Modified data size
for i=1:Datasize
    if(p(i)>0)
        MS = MS+1;
    end
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

TP2(1:MS,1:2) = 0;
j=1;
for i=1:34890
    if(TP(i,1)>0)
        TP2(j,1:2) = TP(i,1:2);
        j=j+1;
    end
end

%--------------------------------------------------------------------------
% TP2 = TP;
%--------------------------------------------------------------------------

P = TP2(:,1);
T = TP2(:,2);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
maxT = max(T);
minT = min(T);
P = (P-min(P))./(max(P)-min(P));
T = (T-min(T))./(max(T)-min(T));

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

Meshsize = 100;
copula(1:Meshsize,1:Meshsize)=0;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
sigma = 1;


for i=1:Meshsize
    for j=1:Meshsize
        for k=1:size(P)
            if(P(k)>=(i-1)/Meshsize && P(k)<i/Meshsize && T(k)>=(j-1)/Meshsize && T(k)<j/Meshsize)
                for m=i-Meshsize:i+Meshsize
                    for n=j-Meshsize:j+Meshsize
                        if(m>0 && m<=Meshsize && n>0 && n<=Meshsize)
                            copula(m,n) = copula(m,n) + exp(-((((m-i)/Meshsize*10)^2+((n-j)/Meshsize*10)^2)/sigma^2));
                        end
                    end
                    
                end
            end
        end
    end
    fprintf('Done %d \n',round(i/Meshsize*100));
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%copula = copula./sum(sum(copula));
% surfl(copula);
% colormap(pink)    % change color map
% shading interp    % interpolate colors across lines and faces
% xticks([0 Meshsize/5 Meshsize/5*2 Meshsize*3/5 Meshsize*4/5 Meshsize])
% xticklabels({'0','0.2','0.4','0.6','0.8','1'})
% yticks([0 Meshsize/5 Meshsize/5*2 Meshsize*3/5 Meshsize*4/5 Meshsize])
% yticklabels({'0','0.2','0.4','0.6','0.8','1'})

cd 'directory for the copula'
csvwrite('Copula.csv',copula);
csvwrite('maxT.csv',maxT);
csvwrite('minT.csv',minT);