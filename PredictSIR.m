function [Pred_Inc,Pred_Incb,A_r,R_c,R_e,d,z1,z2,Imax,S_Inf]=PredictSIR(m,T)
%close all

%% Initialize states
%Estimated parameters
x0(1)=m(1);
x0(2)=m(2);
x0(3)=m(3);
beta=m(4);
rho1=m(5);
rho2=0.1;

        function xdot = SIModel(~,x)
             xdot = zeros(3,1);
              S = x(1);
              I = x(2);
              R = x(3);
                          
            xdot(1) = -beta*I*S/(S + I + R);
            xdot(2) = (beta*I*S/(S + I + R))-(rho1+rho2)*I;
            xdot(3) = (rho1+rho2)*I ;
        end
tspan=(1:T);
[time, Modelsolutionx] = ode45(@SIModel,tspan,x0);
% Determination of the peak of reported cases  Peak time of reported cases
% (vector d)
xa=rho1*Modelsolutionx(:,2);
Pred_Inc =interp1(tspan,xa,time);
Prd=[(1:T)',Pred_Inc(1:T,1)/10000];
d1=sortrows(Prd,2);d=d1(T,:);
% Determination of the True peak time and size (vector z1)
xb=(beta*Modelsolutionx(:,2).*Modelsolutionx(:,1))./(Modelsolutionx(:,1)+Modelsolutionx(:,2)+Modelsolutionx(:,3));
Pred_Incb =interp1(tspan,xb,time);
Prdb=[(1:T)',Pred_Incb(1:T,1)/10000];
d2=sortrows(Prdb,2);z1=d2(T,:);
% Minimum number of Susceptible and time at which it is reached
M=[(1:T)',Modelsolutionx(1:T,1)];d3=sortrows(M,2);z2=d3(1,:);
% Control reproduction number
R_c=beta/(rho1+rho2);
% Maximum number of active cases (Imax)
Imax=(1-(1/R_c)*(1+log(R_c)))*(m(1)+m(2)+m(3))*(1/10000);
% Effective reroduction number
R_e=(R_c*Modelsolutionx(1:T,1))./(Modelsolutionx(1:T,1)+Modelsolutionx(1:T,2)+Modelsolutionx(1:T,3));
N=m(1)+m(2)+m(3);
% Attack ratio (A_r)
A_r=(N-Modelsolutionx(1:T,1))/N;S0=m(1);int=log((1:S0)'/N);
% Final epidemic size
SInf=[(1:S0)',abs((R_c*(((1:S0)'/N)-1)-int))];
clear Modelsolutionx M d1 d2 d3 Prdb Prd int S0
vc=sortrows(SInf,2);S_Inf=vc(1,1)*(1/10000);
clear d1 vc d1 Prdb xa xb 
end