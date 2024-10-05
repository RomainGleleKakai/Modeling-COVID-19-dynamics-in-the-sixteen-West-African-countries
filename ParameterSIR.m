function [namesP,LBParam,m,UBParam]=ParameterSIR(k,Country)

tstart=cputime;
% INPUTS
% k=number of replications considered to construct 95%CI around estimated
% parameters values (k=1000)
% Country='Benin','Burkina_Faso','Cape_Verde','Cote_Ivoire','Gambia','Guinea',
% 'Ghana','Guinea_Bissau','Liberia','Mali','Mauritania','Niger','Nigeria','Senegal',
% 'Sierra_Leone','Togo','West_Africa'

% OUTPUTS
% LBParam=Lower Bound of the 16 estimated parameters
% m= mean estimated values for the parameters considered
% UBParam=Upper Bound of the 16 estimated parameters

% m=[S0/1000000, I0, R0, beta, rho1, RMSE1, Rc, RPeakT, RPeakS/10000, PeakT, PeakS/10000, TS_Inf, S_Inf/1000000, Imax/10000, S_InfC/10000, RMSE2];
% For West Africa as a whole, I0=I0/1000 and R0=R0/1000;

% beta=transmission rate; rho1=detection rate of infected individuals;
% RMSE1=RMSE computed on the 90% of the observations; Rc=Basic reproduction
% number (since control measures were being implemented at this time, this
% parameter is more linked to a control reproduction number)
% RPeakT=Peak time of reported cases;RPeakS=peak size of reported cases
% PeakT=True peak time; PeakS=True peak size
% TS_Inf=Time at which a minimum number of Susceptible is reached
% S_Inf= Minimum number of Susceptible (over the course of the epidemic)
% Imax= Maximum number of active cases; S_InfC=Final epidemic size
% S_InfC=Final epidemic size
% RMSE2=RMSE computed on the remaining observations (10%).

% Author: Prof. Romain Glèlè Kakaï, LABEF, University of Abomey-Calavi
% (Benin)


rho2=0.1;
if (strcmp(Country, 'West_Africa') == 1)
Prop=0.95;%proportion of the training observations in a cross validation process
else
Prop=0.9;%proportion of the training observations in a cross validation process
end

T=300;%number of timepoints considered for the estimation of the peak
xc=0.5;%Proportion of country's population representing upper bound of S0

if ischar(Country)
switch Country
case 'Benin' 
data = readmatrix('Benin');
u1=(0.01*12123200:xc*12123200)';u2=(1:100)';
%                S0      R0  beta    rho_1 
LB = [0.01*12123200      0   0.1    0.005];
UB = [xc*12123200      100   0.5    0.10];
case 'Burkina_Faso'
data = readmatrix('Burkina_Faso');
%% initialize parameter space (a guess is fine)
u1=(0.01*20946992:xc*20946992)';u2=(0:200)';
%% Bounds for fitted parameters
%                S0       R0  beta    rho_1 
LB = [0.01*20946992       0   0.1    0.005];
UB = [xc*20946992       200   1.0    0.10];
case 'Cape_Verde'
data = readmatrix('Cape_Verde');
u1=(0.01*556498:xc*556498)';u2=(0:200)';
%              S0       R0  beta    rho_1 
LB = [0.01*556498       0   0.1    0.01];
UB = [xc*556498       100   1.0    0.10];
case 'Cote_Ivoire'
data = readmatrix('Cote_Ivoire');
u1=(0.01*26428999:xc*26428999)';u2=(0:200)';
%              S0         R0  beta    rho_1
LB = [0.01*26428999       0   0.1    0.0001];
UB = [xc*26428999       100   1.0    0.15];
case 'Gambia'
data = readmatrix('Gambia');
u1=(0.01*2421823:xc*2421823)';u2=(0:100)';
%               S0        R0  beta    rho_1 
LB = [0.01*2421823        0   0.1    0.001];
UB = [xc*2421823        100   1.0    0.15];
case 'Guinea'
data = readmatrix('Guinea');
u1=(0.01*13160021:xc*13160021)';u2=(0:200)';
%                S0       R0  beta    rho_1 
LB = [0.01*13160021       0   0.1    0.001];
UB = [xc*13160021       200   0.1    0.10];
case 'Ghana'
data = readmatrix('Ghana');
%% initialize parameter space (a guess is fine)
u1=(0.01*31072945:xc*31072945)';u2=(0:100)';
%% Bounds for fitted parameters
%                S0       R0  beta    rho_1 
LB = [0.01*31072945       0   0.1    0.001];
UB = [xc*31072945       100   1.0    0.15];
case 'Guinea_Bissau'
data = readmatrix('Guinea_Bissau');
%% initialize parameter space (a guess is fine)
u1=(0.01*1971640:xc*1971640)';u2=(0:200)';
%% Bounds for fitted parameters
%               S0       R0  beta    rho_1 
LB = [0.01*1971640       0   0.1    0.001];
UB = [xc*1971640       100   1.0    0.15];
case 'Liberia'
data = readmatrix('Liberia');
%% initialize parameter space (a guess is fine)
u1=(0.01*5066990:xc*5066990)';u2=(0:300)';
%% Bounds for fitted parameters
%               S0       R0  beta    rho_1 
LB = [0.01*5066990       0   0.1    0.001];
UB = [xc*5066990       100   1.0    0.15];
case 'Mali'
data = readmatrix('Mali');
%% initialize parameter space (a guess is fine)
u1=(0.01*20294900:xc*20294900)';u2=(0:200)';
%% Bounds for fitted parameters
%               S0      R0  beta    rho_1 
LB = [0.01*20294900     0   0.1    0.001];
UB = [xc*20294900     200   1.0    0.15];
case 'Mauritania'
data = readmatrix('Mauritania');
%% initialize parameter space (a guess is fine)
u1=(0.01*4659052:xc*4659052)';u2=(0:100)';
%% Bounds for fitted parameters
%               S0       R0  beta    rho_1 
LB = [0.01*4659052       0   0.1    0.001];
UB = [xc*4659052       100   0.5    0.20];
case 'Niger'
data = readmatrix('Niger');
%% initialize parameter space (a guess is fine)
u1=(0.01*24269389:xc*24269389)';u2=(0:400)';
%% Bounds for fitted parameters
%               S0       R0  beta    rho_1 
LB = [0.01*24269389      0   0.1    0.0008];
UB = [xc*24269389      200   1.0    0.20];
case 'Nigeria'
data = readmatrix('Nigeria');
%% initialize parameter space (a guess is fine)
u1=(0.01*206522290:xc*206522290)';u2=(0:500)';
%% Bounds for fitted parameters
%                 S0       R0  beta    rho_1 
LB = [0.01*206522290       0   0.1    0.003];
UB = [xc*206522290       500   1.0    0.10];
case 'Senegal'
data = readmatrix('Senegal');
u1=(0.01*16776618:xc*16776618)';u2=(0:100)';
%                S0      R0  beta    rho_1 
LB = [0.01*16776618      0   0.1    0.015];
UB = [xc*16776618      100   1.0    0.15];
case 'Sierra_Leone'
data = readmatrix('Sierra_Leone');
u1=(0.01*7989949:xc*7989949)';u2=(0:50)';
%                S0       R0  beta    rho_1 
LB = [0.01*7989949        0   0.1    0.001];
UB = [xc*7989949         50   0.5    0.10];
case 'Togo'
data = readmatrix('Togo');
u1=(0.01*8293924:xc*8293924)';u2=(0:100)';
%               S0       R0  beta    rho_1 
LB = [0.01*8293924       0   0.1    0.0001];
UB = [xc*8293924       100   1.0    0.15];
case 'West_Africa'
data = readmatrix('West_Africa');
u1=(0.01*402555230:xc*402555230)';u2=(1:5000)';
%                 S0      R0  beta    rho_1 
LB = [0.01*402555230      0   0.1    0.005];
UB = [xc*402555230     5000   1.5    0.15];
end
end
n=height(data);
function d=CodeFitDataWAentire(K)
tdata = data(1:round(n*K),1);%training dataset

%% Initialize states
%% initialize parameter space 
parameter(1) = datasample(u1,1);               %S0
parameter(2) = datasample(u2,1);               %R0 
parameter(3:4) = rand(1,2);
x0(1) =parameter(1);                           %x0(1)=S(0)  
x0(2) =data(1,3)/parameter(4);                 %x0(2)=I(0)
x0(3) =parameter(2);                           %x0(3)=R(0)


options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',2000,'TolX',1E-4,'TolFun',1E-4);
[Parameter_min, RMSE_min] = fminsearchbnd(@SumOfSquaresOfErrors,parameter,LB,UB,options);
Parameter_min=[Parameter_min(1),data(1,3)/Parameter_min(4),Parameter_min(2:4)];
[~,~,~,R_c,~,U1,U2,U3,Imax,S_Inf]=PredictSIR(Parameter_min,T);
d=[Parameter_min,RMSE_min/10000,R_c,U1,U2,U3,Imax,S_Inf];

clear Parameter_min Parameter_mina R_c U1 U2 u3 Imax S_Inf

 function RMSE1 = SumOfSquaresOfErrors(parameter)
              
        x0(1)=parameter(1);   
        x0(2)=data(1,3)/parameter(4);
        x0(3)=parameter(2);
        beta=parameter(3);
        rho1=parameter(4);
                           
        function xdot = SIModel(~,x)
              xdot = zeros(3,1);
              S = x(1);
              I = x(2);
              R = x(3);
                          
            xdot(1) = -beta*I*S/(S + I + R);
            xdot(2) = (beta*I*S/(S + I + R))-(rho1+rho2)*I;
            xdot(3) = (rho1+rho2)*I ;
        end
%  Prepare to solve ode
 tspan = 1:1:n; % set time points
 [time,Modelsolution] = ode45(@SIModel,tspan,x0); % solve ode
 PredictIncidence = rho1*Modelsolution(:,2);
 Pred_Inc =(interp1(time,PredictIncidence,tdata));
 M1=(Pred_Inc(1:round(n*K),1)-data(1:round(n*K),3)).^2;
 RMSE1=sqrt(sum(M1)/round(n*K));
 clear ModelSolution Pred_Inc PredictIncidence M1 
end
end
da=zeros(k,15);
for e=1:k
da(e,:)=CodeFitDataWAentire(Prop);
end
za=sortrows(da,6);
RMSE=zeros(k,1);
clear da
for w=1:k
[Pred_Inci,~,~,~,~,~,~,~,~,~]=PredictSIR(za(w,1:5),T);
Pred_Inc1=Pred_Inci;
IncidenceData1=data(1:n,3);
RM1=(Pred_Inc1(round(n*Prop)+1:n,1)-IncidenceData1(round(n*Prop)+1:n,1)).^2;
RMSE(w,1)=sqrt(sum(RM1)/(n-round(n*Prop)));
end
clear Pred_Inci Pred_Inc1 RM1
if (strcmp(Country, 'West_Africa') == 1)
du1=[za,RMSE/1000];
du=sortrows(du1,16);du(:,14)=du(:,14)/1000000;du(:,13)=du(:,13)/1000000;
du(:,15)=du(:,15)/1000000;du(:,1)=du(:,1)/1000000;du(:,2)=du(:,2)/1000;du(:,3)=du(:,3)/1000;
else
du1=[za,RMSE/1000];
du=sortrows(du1,16);du(:,14)=du(:,14)/1000000;du(:,13)=du(:,13)/1000000;
du(:,15)=du(:,15)/1000000;du(:,1)=du(:,1)/1000000;
end
clear za RMSE
namesP=['  S0 ','      I0 ','     R0 ','       beta ','     rho1 ','     RMSE1','      Rc','     RPeakT','     RPeakS','   PeakT','      PeakS','   TS_Inf','     S_Inf','      Imax','    S_InfC','     RMSE2']; 
[LBParam,m,UBParam]=Curves(T,k,du,Country);
 for i=1:16
    if LBParam(i)<0
       LBParam(i)=0;
    end
end
time=(cputime-tstart)/60;
display(time)
end