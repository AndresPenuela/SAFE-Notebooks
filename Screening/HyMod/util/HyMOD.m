clc
clearvars -except P Qm PET

disp ('***************************************************************************')
disp ('HyMOD MATLAB code developed by M. H. Alipour (malipour@knights.ucf.edu)')
disp ('February 4, 2017')
disp ('Papers and other scientific works based in whole or in part on this code shall cite:')
disp ('Alipour, M.H., Kibler, K.M., 2018. A framework for streamflow prediction in the ') 
disp ('world’s most severely data-limited regions: Test of applicability and performance ') 
disp ('in a poorly-gauged region of China. Journal of Hydrology, 557, pp.41-54.')
disp ('***************************************************************************')

year = input('Please insert the first year when observed runoff data is available (example: 2016): ');
month = input('Please insert the first month when observed runoff data is available (example: 1 for January): ');
day = input('Please insert the first day when observed runoff data is available (example: 1 for the first day of the month): '); 
check = 0;

while check == 0
Bdis = input('Please insert number of intervals for initial discretization of parameter B (feasible range: 0.01-4): ');
Alphadis = input('Please insert number of intervals for initial discretization of parameter Alpha (feasible range: 0.01-0.99): ');
Rsdis = input('Please insert number of intervals for initial discretization of parameter Rs (feasible range: 0.0001-0.01  1/day): ');
Rqdis = input('Please insert number of intervals for initial discretization of parameter Rq (feasible range: 0.01-0.99  1/day): ');
Cmaxdis = input('Please insert number of intervals for initial discretization of parameter Cmax (feasible range: 5-1500 mm): ');
Nseq = input('Please insert number of iterations: ');
Ndiv = input('Please insert number of divisions: ');
NC = Nseq*(Bdis+1)*(Alphadis+1)*(Rsdis+1)*(Rqdis+1)*(Cmaxdis+1);
disp ('***************************************************************************')
disp (sprintf('Number of Combinations = %0.0f', NC))
check = input('Would you like to proceed (Please type 0 for No and 1 for Yes)?');
IP = input('Please insert number of time steps for initialization period: ');
clc
end

tic
syms S Ss S1 S2 S3 C Cstar ER1 ER2 e Sstar Qs Qq1 Qq2 Qq3 Q
tmax = length(Qm);
Sstar = 0;
RESNEW = 10^12;
denominator = sum((Qm(IP+1:length(Qm))-mean(Qm(IP+1:length(Qm)))).^2);

for seq=1:Nseq

if seq == 1
    Bul = 4;
    Bll = 0.01;
    Alphaul = 0.99;
    Alphall = 0.01;
    Rsul = 0.01;
    Rsll = 0.0001;
    Rqll = 0.01;
    Rqul = 0.99;
    Cmaxll = 5;
    Cmaxul = 1500;
end

if seq > 1
    if B+(Bul-Bll)/(2*Ndiv) <= Bul && B-(Bul-Bll)/(2*Ndiv) >= Bll
        BUL=Bul;
        Bul = B+(Bul-Bll)/(2*Ndiv);
        Bll = B-(BUL-Bll)/(2*Ndiv);
    else
        if B+(Bul-Bll)/(2*Ndiv) <= Bul && B-(Bul-Bll)/(2*Ndiv) < Bll
            Bul = Bll+(Bul-Bll)/Ndiv;
            Bll = Bll;
        else
            if B+(Bul-Bll)/(2*Ndiv) > Bul && B-(Bul-Bll)/(2*Ndiv) >= Bll
                Bul = Bul;
                Bll = Bul-(Bul-Bll)/Ndiv;
            end
        end
    end
    
    if Alpha+(Alphaul-Alphall)/(2*Ndiv) <= Alphaul && Alpha-(Alphaul-Alphall)/(2*Ndiv) >= Alphall
        AlphaUL=Alphaul;
        Alphaul = Alpha+(Alphaul-Alphall)/(2*Ndiv);
        Alphall = Alpha-(AlphaUL-Alphall)/(2*Ndiv);
    else
        if Alpha+(Alphaul-Alphall)/(2*Ndiv) <= Alphaul && Alpha-(Alphaul-Alphall)/(2*Ndiv) < Alphall
            Alphaul = Alphall+(Alphaul-Alphall)/Ndiv;
            Alphall = Alphall;
        else
            if Alpha+(Alphaul-Alphall)/(2*Ndiv) > Alphaul && Alpha-(Alphaul-Alphall)/(2*Ndiv) >= Alphall
                Alphaul = Alphaul;
                Alphall = Alphaul-(Alphaul-Alphall)/Ndiv;
            end
        end
    end
    
    if Rs+(Rsul-Rsll)/(2*Ndiv) <= Rsul && Rs-(Rsul-Rsll)/(2*Ndiv) >= Rsll
        RsUL=Rsul;
        Rsul = Rs+(Rsul-Rsll)/(2*Ndiv);
        Rsll = Rs-(RsUL-Rsll)/(2*Ndiv);
    else
        if Rs+(Rsul-Rsll)/(2*Ndiv) <= Rsul && Rs-(Rsul-Rsll)/(2*Ndiv) < Rsll
            Rsul = Rsll+(Rsul-Rsll)/Ndiv;
            Rsll = Rsll;
        else
            if Rs+(Rsul-Rsll)/(2*Ndiv) > Rsul && Rs-(Rsul-Rsll)/(2*Ndiv) >= Rsll
                Rsul = Rsul;
                Rsll = Rsul-(Rsul-Rsll)/Ndiv;
            end
        end
    end
    
    if Rq+(Rqul-Rqll)/(2*Ndiv) <= Rqul && Rq-(Rqul-Rqll)/(2*Ndiv) >= Rqll
        RqUL=Rqul;
        Rqul = Rq+(Rqul-Rqll)/(2*Ndiv);
        Rqll = Rq-(RqUL-Rqll)/(2*Ndiv);
    else
        if Rq+(Rqul-Rqll)/(2*Ndiv) <= Rqul && Rq-(Rqul-Rqll)/(2*Ndiv) < Rqll
            Rqul = Rqll+(Rqul-Rqll)/Ndiv;
            Rqll = Rqll;
        else
            if Rq+(Rqul-Rqll)/(2*Ndiv) > Rqul && Rq-(Rqul-Rqll)/(2*Ndiv) >= Rqll
                Rqul = Rqul;
                Rqll = Rqul-(Rqul-Rqll)/Ndiv;
            end
        end
    end
    
    if Cmax+(Cmaxul-Cmaxll)/(2*Ndiv) <= Cmaxul && Cmax-(Cmaxul-Cmaxll)/(2*Ndiv) >= Cmaxll
        CmaxUL=Cmaxul;
        Cmaxul = Cmax+(Cmaxul-Cmaxll)/(2*Ndiv);
        Cmaxll = Cmax-(CmaxUL-Cmaxll)/(2*Ndiv);
    else
        if Cmax+(Cmaxul-Cmaxll)/(2*Ndiv) <= Cmaxul && Cmax-(Cmaxul-Cmaxll)/(2*Ndiv) < Cmaxll
            Cmaxul = Cmaxll+(Cmaxul-Cmaxll)/Ndiv;
            Cmaxll = Cmaxll;
        else
            if Cmax+(Cmaxul-Cmaxll)/(2*Ndiv) > Cmaxul && Cmax-(Cmaxul-Cmaxll)/(2*Ndiv) >= Cmaxll
                Cmaxul = Cmaxul;
                Cmaxll = Cmaxul-(Cmaxul-Cmaxll)/Ndiv;
            end
        end
    end
end
    
for B = Bll:(Bul-Bll)/Bdis:Bul
    for Alpha = Alphall:(Alphaul-Alphall)/Alphadis:Alphaul
        for Rs = Rsll:(Rsul-Rsll)/Rsdis:Rsul
            for Rq = Rqll:(Rqul-Rqll)/Rqdis:Rqul
                for Cmax = Cmaxll:(Cmaxul-Cmaxll)/Cmaxdis:Cmaxul

Sstar = 0;                    
Smax = Cmax/(B+1);
Ss = 0;
S1 = 0;    
S2 = 0;
S3 = 0;
RES = 0;

for t = 1:tmax
    C = Cmax*(1-(1-((B+1)*Sstar)/(Cmax))^(1/(B+1)));
    ER1 = max(P(t)+C-Cmax,0);
    Cstar = min(P(t)+C,Cmax);
    S = (Cmax/(B+1))*(1-(1-(Cstar/Cmax))^(B+1));
    e = min(PET(t)*Cstar/Cmax,S);
    ER2 = max((Cstar-C)-(S-Sstar),0);
    Sstar = S-e;
    Ss = (1-Rs)*Ss+(1-Rs)*(1-Alpha)*ER2;
    Qs = (Rs/(1-Rs))*Ss;
    S1 = (1-Rq)*S1+(1-Rq)*(ER1+Alpha*ER2);
    Qq1 = (Rq/(1-Rq))*S1;
    S2 = (1-Rq)*S2+(1-Rq)*Qq1;
    Qq2 = (Rq/(1-Rq))*S2;
    S3 = (1-Rq)*S3+(1-Rq)*Qq2;
    Qq3 = (Rq/(1-Rq))*S3;
    Q = Qs+Qq3;
    QQQ(t) = Q;
    if S<0 | Ss<0 | S1<0 | S2<0 | S3<0 | C<0 | Cstar<0 | e<0 | ER1<0 | ER2<0 | Qs<0 | Qq1<0 | Qq2<0 | Qq3<0
        display ('infeasible')
        display ('*******************************************')
        break
    end
    if t > IP
        RES = RES + (Qm(t)-Q)^2/(denominator);
    end
end

if RES < RESNEW && t == tmax
    RESNEW = RES;
    NSE = 1-RES
    BNew = B
    AlphaNew = Alpha
    RsNew = Rs
    RqNew = Rq
    CmaxNew = Cmax
    QQ = QQQ;
end

                end
            end
        end   
    end
end
B = BNew;
Alpha = AlphaNew;
Rs = RsNew;
Rq = RqNew;
Cmax = CmaxNew;
end


clc 
B = BNew
Alpha = AlphaNew
Rs = RsNew
Rq = RqNew
Cmax = CmaxNew
NSE

disp ('***************************************************************************')
disp ('HyMOD MATLAB code developed by M. H. Alipour (February 4, 2017)')
disp ('Papers and other scientific works based in whole or in part on this code shall cite:')
disp ('Alipour, M.H., Kibler, K.M., 2018. A framework for streamflow prediction in the ') 
disp ('world’s most severely data-limited regions: Test of applicability and performance ') 
disp ('in a poorly-gauged region of China. Journal of Hydrology, 557, pp.41-54.')
disp ('***************************************************************************')

t0 = datenum(year,month,day);
dt = 1;
t = datestr(t0 + (IP:tmax-1)*dt);
tt = datenum(t);
plot(tt,Qm(IP+1:length(Qm)))
xlabel ('Date','FontName','Times New Roman','fontsize',12);
ylabel ('Q (mm/day)','FontName','Times New Roman','fontsize',12);
mid = floor((length(Qm) - IP - 1)/4);
set(gca, 'XTick', tt(1:mid:end));
dateFormat = 'dd-mmm-yyyy';
datetick('x',dateFormat,'keepticks');
hold on
plot (tt,QQ(IP+1:length(Qm)),'r');
plot (tt,QQ(IP+1:length(Qm)),'r');
h_legend = legend('Observed Streamflow', 'Modeled Streamflow', sprintf('Nash-Sutcliffe Efficiency = %0.2f', NSE));
set(h_legend,'FontName','Times New Roman','fontsize',12);
toc