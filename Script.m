% Beleg 2 - Autokovarianzfunktion und Energiespektrum

clear; close all; clc; format long;
%*******aus Beleg 1***********

% Einlesen Zeitreihe
y = readmatrix("aufgabe9.txt");
y = y(1:1000,1:3);
t = y(:,2);
an = [88.692361725097570;-42.535052398364826;40.250629677348464;23.135242413340411];
bn = [100.1013697891288;-33.8270656240537;15.4510264918595;-33.2975572844197];
omega = [28.98*pi/180;15.12*pi/180;13.86*pi/180;30.06*pi/180];

% Bereinigung von Trend und Offset
detrended = detrend(y(:,3));
trend_fkt = polyfit(y(:,2),y(:,3),1);
trend = trend_fkt(1,1)*y(:,2)+trend_fkt(1,2);


%*****************************

% Aufgabe 1

% Ausgleichung

I = length(detrended);
I2 = I/2;

L = y(:,3);



% Anzahl der Beobachtungen:
n = 1000;
% Anzahl der Unbekannten:
u = 10;
% Freiheitsgrade:
f = n-u;

% Naeherungsvektor der Unbekannten:
m_ = trend_fkt(1);
q_ = trend_fkt(2);
a1_ = an(1);
a2_ = an(2);
a3_ = an(3);
a4_ = an(4);
b1_ = bn(1);
b2_ = bn(2);
b3_ = bn(3);
b4_ = bn(4);

X0 = [a1_, b1_, a2_, b2_, a3_, b3_, a4_, b4_, m_, q_];

A = [];
    A(:,1) = cos(omega(1).*t);
    A(:,2) = sin(omega(1).*t);
    A(:,3) = cos(omega(2).*t);
    A(:,4) = sin(omega(2).*t);
    A(:,5) = cos(omega(3).*t);
    A(:,6) = sin(omega(3).*t);
    A(:,7) = cos(omega(4).*t);
    A(:,8) = sin(omega(4).*t);
    A(:,9) = t(:);
    A(:,10) = 1;

Li = ones(1000,1)*9999;
Ldach = ones(1000,1)*-9999;
while round(Li,6) ~= round(Ldach,6)

L0 = [];
L0(:,1) = m_ .*t + q_ + a1_*cos(omega(1).*t) + b1_*sin(omega(1).*t) + a2_*cos(omega(2).*t) + b2_*sin(omega(2).*t) + a3_*cos(omega(3).*t) + b3_*sin(omega(3).*t) + a4_*cos(omega(4).*t) + b4_*sin(omega(4).*t);

% gekuerzter Beobachtungsvektor:
l = L - L0;



% Stochastisches Modell
% ---------------------
% P = Einheitsmatrix


% Ausgleichung
% ------------


[N,Qxx,xdach,v,W,var0,Sigma_xx] = ausgleichung(A,l);



% Ausgeglichene Beobachtungen
Ldach = L + v;

% Hauptrechenprobe
Xdach = [ X0 ]' + xdach; % Ausgeglichene Parameter


a1_= Xdach(1);
b1_= Xdach(2);
a2_= Xdach(3);
b2_= Xdach(4);
a3_= Xdach(5);
b3_= Xdach(6);
a4_= Xdach(7);
b4_= Xdach(8);
m_= Xdach(9);
q_ = Xdach(10);


Li = [];
Li(:,1) = m_ .*t + q_ + a1_*cos(omega(1).*t) + b1_*sin(omega(1).*t) + a2_*cos(omega(2).*t) + b2_*sin(omega(2).*t) + a3_*cos(omega(3).*t) + b3_*sin(omega(3).*t) + a4_*cos(omega(4).*t) + b4_*sin(omega(4).*t);

if round(Li,6) == round(Ldach,6)
    fprintf('Hauptrechenprobe war erfolgreich.\n')
else
    fprintf('Hauptrechenprobe ist fehlgeschlagen.\n')
end




end

sigma_xx = diag(Sigma_xx);

% clearvars -except Xdach sigma_xx t L detrended omega

%
%%

a1= Xdach(1);
b1= Xdach(2);
a2= Xdach(3);
b2= Xdach(4);
a3= Xdach(5);
b3= Xdach(6);
a4= Xdach(7);
b4= Xdach(8);
m= Xdach(9);
q = Xdach(10);

stdabw = sqrt(sigma_xx);

% Rekonstruktion

ydach = [];
T = 2000;

% ydach(:,1) = m.*t + q + a1*cos(1*2*pi.*t/T) + b1*sin(1*2*pi.*t/T) + a2*cos(2*2*pi.*t/T) + b2*sin(2*2*pi.*t/T) + a3*cos(3*2*pi.*t/T) + b3*sin(3*2*pi.*t/T) + a4*cos(4*2*pi.*t/T) + b4*sin(4*2*pi.*t/T);
ydach(:,1) = m.*t + q + a1*sin(omega(1).*t) + b1*cos(omega(1).*t) + a2*sin(omega(2).*t) + b2*cos(omega(2).*t) + a3*sin(omega(3).*t) + b3*cos(omega(3).*t) + a4*sin(omega(4).*t) + b4*cos(omega(4).*t);



figure(1)
hold on
plot(t,ydach, '-r')
plot(t,L, 'b.')
legend('reduzierte Zeitreihe','Ausgangszeitreihe')
hold off

saveas(1,'expdavid/1.png')

figure(2)
hold on
plot(t,L-ydach,'b.')
legend('Differenz zwischen reduzierter und Ausgangszeitreihe')
hold off

saveas(2,'expdavid/2.png')


%
%%
% Mittelwert und Varianz
close all

dL = L-ydach;

v = var(dL);
mu = mean(dL);
sd = sqrt(v);

%
%%
% Häufigkeitsverteilung
% clearvars -except dL
close all
lauf = 17;
z = -400:lauf:400;
b = hist(dL,z);

figure(3)
hold on
bhist = histogram(dL,z);

saveas(3,'expdavid/3_histogram.png')

% bhist

haeuf = b/1000;
r = haeuf/lauf;
x = linspace(-400,400,50);
f = 1/(sqrt(2*pi)*sd)*exp((-1/2)*(x-mu).^2/sd^2);
f1 = @(x) 1/(sqrt(2*pi)*sd)*exp((-1/2)*(x-mu).^2/sd^2);
integral(f1,-400,400)
plot(x,f)


hold off

figure(4)
histfit(dL,48,'normal')

saveas(4,'expdavid/4_histfit.png')

% syms n
% limit(haeuf(:)/n, n, inf)

figure(5)
hold on
bar(z,r,'BarWidth', 1)
plot(x,f,'LineWidth',2)
hold off

saveas(5,'expdavid/5_angep.png')


%
%%
% empirische Schätzung AKF
dt = t(2)-t(1);
tau_max = T/50;
N = tau_max/dt;
r = 0:N;
n = length(dL);

Cemp = []; 
for i = 1:length(r)

    C1 = 1/(n-r(i)-1);
    Cemp(i,1) = 0;
    for j = 1:(n-r(i))
        Cemp(i,1) = Cemp(i,1) + (dL(j)-mu)*(dL(j+r(i))-mu);
    end
    Cemp(i,1) = Cemp(i,1) * C1;
end

figure(6)
hold on
plot(r,Cemp,'b-')
hold off

saveas(6,'expdavid/6_empAKF.png')


%
%%
% Anpassung empirische AKF an analytische AKF
clearvars -except dL mu sd Cemp C1 t T dt r tau_max
close all

alpha_ = 1;
Canal = [];
Canal(:,1) = (sd^2) * exp(-alpha_*abs(r.*dt));

% z = polyfit(r, Cemp,6);


P = ones(1,21);
P(1,1:2) = 1000;
P(1,3:21) = 1/1000;
P = diag(P);


Li = ones(length(Cemp),1)*9999;
Ldach = ones(length(Cemp),1)*-9999;
k = 0;
% while round(Li,6) ~= round(Ldach,6)
while k ~= 10

    A = [];
    A(:,1) = (sd^2) .* abs(r.*dt) .* (-exp(-alpha_.*abs(r.*dt)));

    L0(:,1) = (sd^2) * exp(-alpha_*abs(r.*dt));

    l = Cemp - L0;

    [N,Qxx,xdach,v,W,var0,Sigma_xx] = ausgleichung2(A,P,l,1);

    % Ausgeglichene Beobachtungen
    Ldach = Canal - v;
    
    % Hauptrechenprobe
    alpha_ = alpha_ + xdach; % Ausgeglichene Parameter
    
    
    
    Li = [];
    Li(:,1) = (sd^2) * exp(-alpha_*abs(r.*dt));
    
    if round(Li,6) == round(Ldach,6)
        fprintf('Hauptrechenprobe war erfolgreich.\n')
    else
        fprintf('Hauptrechenprobe ist fehlgeschlagen.\n')
    end

    k = k+1;

end


figure(7)
hold on
plot(r,Cemp,'b-')
plot(r, Ldach,'r-')
hold off

saveas(7,'expdavid/7.png')

Canal = [];
r1=0:0.01:20;
Canal(:,1) = (sd^2) * exp(-alpha_*abs(r1.*dt));

figure(8)
hold on
plot(r,Cemp,'b-')
plot(r1,Canal,'r-')
hold off

saveas(8,'expdavid/8_AKFangep.png')

%
%%
% Spektraldichte

omega_ny = pi/((r(2)*dt)-(r(1)*dt));
domega = 2*pi/tau_max;
omega_line = [];
n = 1:round(omega_ny/domega);
omega_line(:,1) = n .* domega; 

s = [];
s(:,1) = 2*alpha_./(omega_line.^2+alpha_^2);

figure(9)
plot(omega_line, s)

saveas(9,'expdavid/9_Spektraldichteexp.png')

s2(:,1) = (sqrt(pi)/alpha_) * exp(-(omega_line/(2*alpha_)).^2);

figure(10)
plot(omega_line, s2)

saveas(10,'expdavid/10_Spektraldichtegauss.png')
