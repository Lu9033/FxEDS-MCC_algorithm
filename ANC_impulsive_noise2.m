clear
clc
clf
close all

load Sw.mat
load Pw.mat
Pw=Pw(1:3:768);
Sw=[Sw(1:2:248) Sw(200:4:212)];

T = 20000; 
inde = 1; 
Len = 192;
ANR_EDS = zeros(1,T);
ANR_MCC = zeros(1,T);
ANR_RMCC = zeros(1,T);
ANR_EDSMCC = zeros(1,T);
fg=0.9999;
w=zeros(T,Len);


AE_EDS = zeros(1,T);
AD = zeros(1,T);
AE_MCC = zeros(1,T);
AE_RMCC = zeros(1,T);
AE_EDSMCC = zeros(1,T);

for ind = 1:inde
Shx = zeros(1,Len);  
Shw = zeros(1,Len);  
bulin=zeros(1,Len-size(Sw,2));
Shw = [Sw bulin];
pr=0.001;
R=binornd(1,pr,1,T);
r_EDS=randn(1,T);
    for i=1:T
        if r_EDS(i)>0
            r_EDS(i)=sqrt(2000);
        else
            r_EDS(i)=-sqrt(2000);
        end
    end
impulsive_noise=r_EDS.*R;

X_BG=awgn(impulsive_noise,1);


Y=filter(Pw,1,X_BG);

Cx = zeros(1,Len); 
%FxEDS
Cw_EDS = zeros(1,Len); 
Sx_EDS = zeros(size(Sw)); 
e_cont_EDS=zeros(1,T);

%FxMCC
Cw_MCC = zeros(1,Len); 
Sx_MCC = zeros(size(Sw)); 
e_cont_MCC = zeros(1,T);

%FxRMCC
Cw_RMCC = zeros(1,Len); 
Sx_RMCC = zeros(size(Sw)); 
e_cont_RMCC = zeros(1,T);

%FxEDS-MCC
Cw_EDSMCC = zeros(1,Len); 
Sx_EDSMCC = zeros(size(Sw)); 
e_cont_EDSMCC = zeros(1,T);

%alpha = 1.7
deita_MCC = 1;
mu_MCC = 0.0005;

%FxRMCC
forget_factor_RMCC = 0.999;
del_RMCC = 5;
p_RMCC = (1/del_RMCC)*eye(Len,Len);
xige = 1;





%FxEDS-MCC
deita_EDSMCC = 1.5;
lamda_EDSMCC = 0.9999;


Xhx=zeros(1,Len); 

%FxEDS
Q_EDS = zeros(Len,Len);
r_EDS = zeros(1,Len);
Cw_EDS = zeros(1,Len);
Sx_EDS = zeros(size(Sw));

%FxMCC
Cw_MCC = zeros(1,Len);
Sx_MCC = zeros(size(Sw));

%FxRMCC
Cw_RMCC = zeros(1,Len);
Sx_RMCC = zeros(size(Sw));

%FxRMCC_variable_width
Cw_RMCC_vw = zeros(1,Len);
Sx_RMCC_vw = zeros(size(Sw));
my =zeros(1,T);

%FXEDSMCC
Q_EDSMCC = zeros(Len,Len);
r_EDSMCC = zeros(1,Len);
Cw_EDSMCC = zeros(1,Len);
Sx_EDSMCC = zeros(size(Sw));

for p=1:Len
    Cx = [X_BG(p) Cx(1:Len-1)];
    Shx = [X_BG(p) Shx(1:Len-1)];  
    Xhx = [sum(Shx.*Shw) Xhx(1:Len-1)]; 
end


for k=Nw:T
    
    Cx = [X_BG(k) Cx(1:Len-1)]; 
    Shx = [X_BG(k) Shx(1:Len-1)];  
    Xhx = [sum(Shx.*Shw) Xhx(1:Len-1)]; 
    
    %FxEDS
    Cy_EDS = sum(Cx.*Cw_EDS);
    Sx_EDS = [Cy_EDS Sx_EDS(1:length(Sx_EDS)-1)]; 
    e_cont_EDS(k) = Y(k) - sum(Sx_EDS.*Sw);
    
    %FxMCC
    Cy_MCC = sum(Cx.*Cw_MCC);
    Sx_MCC = [Cy_MCC Sx_MCC(1:length(Sx_MCC)-1)]; 
    e_cont_MCC(k) = Y(k) - sum(Sx_MCC.*Sw);
    vv_MCC(k) = exp(-e_cont_MCC(k)^2/(2*deita_MCC^2));
    Cw_MCC = Cw_MCC + mu_MCC*vv_MCC(k)*e_cont_MCC(k)*Xhx;
    
    %FxRMCC
    Cy_RMCC = sum(Cx.*Cw_RMCC);               	
    Sx_RMCC = [Cy_RMCC Sx_RMCC(1:length(Sx_RMCC)-1)];
    e_cont_RMCC(k) = Y(k)-sum(Sx_RMCC.*Sw);  
    vv_RMCC(k) = exp(-e_cont_RMCC(k)^2/(2*xige));
    kar_RMCC = vv_RMCC(k)*p_RMCC*Xhx' / (forget_factor_RMCC + vv_RMCC(k)*Xhx*p_RMCC*Xhx');
    Cw_RMCC = Cw_RMCC + kar_RMCC'*e_cont_RMCC(k);
    p_RMCC = (1/forget_factor_RMCC)*(p_RMCC - kar_RMCC*Xhx*p_RMCC);
    
    %FxEDSMCC
    Cy_EDSMCC = sum(Cx.*Cw_EDSMCC);
    Sx_EDSMCC = [Cy_EDSMCC Sx_EDSMCC(1:length(Sx_EDSMCC)-1)]; 
    e_cont_EDSMCC(k) = Y(k) - sum(Sx_EDSMCC.*Sw);
    vv_EDSMCC(k) = exp(-e_cont_EDSMCC(k)^2/(2*deita_EDSMCC^2));
    
    %FxEDS
    Q_EDS = lamda_EDS*Q_EDS+Xhx'*Xhx;
    r_EDS = lamda_EDS*r_EDS+Y(k)*Xhx;
    
    %FxEDSMCC
    Q_EDSMCC = lamda_EDSMCC*Q_EDSMCC+vv_EDSMCC(k)*Xhx'*Xhx;
    r_EDSMCC = lamda_EDSMCC*r_EDSMCC+vv_EDSMCC(k)*Y(k)*Xhx;
    
    

        for ii=1:Len
            gg = zeros(1,Len);
            gg(ii)=1;
            %FxEDS
            test(k,:) = Cw_EDS;
            alpha_EDS = -((Cw_EDS*Q_EDS - r_EDS)*gg')/(gg*Q_EDS*gg' + 0.001);
            a_EDS(k,:) = alpha_EDS;
            Cw_EDS(ii) = Cw_EDS(ii) + alpha_EDS;
            %FxEDSMCC
            test(k,:)=Cw_EDSMCC;
            alpha_EDSMCC = -((Cw_EDSMCC*Q_EDSMCC - r_EDSMCC)*gg')/(gg*Q_EDSMCC*gg' + 0.001);
            a_EDSMCC(k,:) = alpha_EDSMCC;
            Cw_EDSMCC(ii) = Cw_EDSMCC(ii) + alpha_EDSMCC;
            
        end
        AE_EDS(k) = fg*AE_EDS(k-1) + (1 - fg)*abs(e_cont_EDS(k));
        AE_MCC(k) = fg*AE_MCC(k-1) + (1 - fg)*abs(e_cont_MCC(k));
        AE_RMCC(k) = fg*AE_RMCC(k-1) + (1 - fg)*abs(e_cont_RMCC(k));
        AE_EDSMCC(k) = fg*AE_EDSMCC(k-1) + (1 - fg)*abs(e_cont_EDSMCC(k));
        
        AD(k) = fg*AD(k-1) + (1 - fg)*abs(Y(k));
        
        ANR_EDS(k) = ANR_EDS(k) + 20*log(AE_EDS(k)/AD(k));
        ANR_MCC(k) = ANR_MCC(k) + 20*log(AE_MCC(k)/AD(k));
        ANR_RMCC(k) = ANR_RMCC(k) + 20*log(AE_RMCC(k)/AD(k));
        ANR_EDSMCC(k) = ANR_EDSMCC(k) + 20*log(AE_EDSMCC(k)/AD(k));
        count=count+1

end

end
        
        ANR_EDS = ANR_EDS/inde;
        ANR_MCC = ANR_MCC/inde;
        ANR_RMCC = ANR_RMCC/inde;
        ANR_EDSMCC = ANR_EDSMCC/inde;
        
        HL=fir1(20,0.1);%Frequencies must fall in range between 0 and 1.
        ANR_EDS = filter(HL,1,ANR_EDS);
        ANR_MCC = filter(HL,1,ANR_MCC);
        ANR_RMCC = filter(HL,1,ANR_RMCC);
        ANR_EDSMCC = filter(HL,1,ANR_EDSMCC);
        figure(7)
        plot(ANR_EDS,'c-.','LineWidth',2);
        hold on
        plot(ANR_MCC,'b-.','LineWidth',2);
        hold on
        plot(ANR_RMCC,'g:','LineWidth',2); 
        hold on
        plot(ANR_EDSMCC,'r','LineWidth',2);
        hold off
        ylabel('ANR,dB');
        xlabel('time index k')
        legend('FxEDS','FxMCC','FxRMCC','FxEDSMCC')%
 