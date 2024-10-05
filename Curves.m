function [LBParam,m,UBParam]=Curves(T,k,du,Country)

if ischar(Country)
switch Country
case 'Benin' 
data = readmatrix('Benin');
case 'Burkina_Faso'
data = readmatrix('Burkina_Faso');
case 'Cape_Verde'
data = readmatrix('Cape_Verde');
case 'Cote_Ivoire'
data = readmatrix('Cote_Ivoire');
case 'Gambia'
data = readmatrix('Gambia');
case 'Guinea'
data = readmatrix('Guinea');
case 'Ghana'
data = readmatrix('Ghana');
case 'Guinea_Bissau'
data = readmatrix('Guinea_Bissau');
case 'Liberia'
data = readmatrix('Liberia');
case 'Mali'
data = readmatrix('Mali');
case 'Mauritania'
data = readmatrix('Mauritania');
case 'Niger'
data = readmatrix('Niger');
case 'Nigeria'
data = readmatrix('Nigeria');
case 'Senegal'
data = readmatrix('Senegal');
case 'Sierra_Leone'
data = readmatrix('Sierra_Leone');
case 'Togo'
data = readmatrix('Togo');
case 'West_Africa'
data = readmatrix('West_Africa');
end
end
Day=data(:,5);n=height(data);
m=du(1,:);mb=m;

if (strcmp(Country, 'West_Africa') == 1)
mb(1)=mb(1)*1000000;mb(2)=mb(2)*1000;mb(3)=mb(3)*1000;
else
mb(1)=mb(1)*1000000;   
end

s=std(du);LBParam=m-((2.04*s)/sqrt(k));UBParam=m+((2.04*s)/sqrt(k));
[Pa1,Pa2,A_r,~,R_e,~,~,~,~,~]=PredictSIR(mb,T);

% Daily reported over time
set(0,'DefaultAxesFontName', 'Times New Roman','defaultaxesfontsize',15,'defaulttextinterpreter','latex');
X=Day(1):Day(1)+n-1;
date = (datetime(X,'ConvertFrom','excel'))';
Xf=Day(1):Day(1)+T-1;
datef = (datetime(Xf,'ConvertFrom','excel'))';
figure(1)
plot(date(1:n,1),data(1:n,3),'or','MarkerSize',5);
hold on
plot(datef(1:T,1),Pa1,'r-','LineWidth',2);
plot(datef(1:T,1),Pa2,'m-','LineWidth',2);
hold off
ylabel('Daily cases')
legend('Reported cases', 'Predicted Reported cases', 'Predicted Total cases')
ax = gca;
ax.XTickLabelRotation = 45;
figure(2)
plot(datef(1:T,1),A_r,'r','LineWidth',2);
ylabel('Attack Ratio')
figure(3)
plot(datef(1:T,1),R_e,'r','LineWidth',2);
ylabel('Effective reproduction number')
ax = gca;
ax.XTickLabelRotation = 45;
end