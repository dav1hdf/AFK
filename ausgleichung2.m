function [N,Qxx,xdach,v,W,var0,Sigma_xx] = ausgleichung(A,P,l,theo_var)

N = A'*P*A; % Normalgleichungsmatrix
n = A'*P*l; % Absolutglied
Qxx = inv(N); % Inversion Normalgleichungsmatrix
xdach = Qxx*n; % Gekürzter Parametervektor
v = A*xdach - l; % Verbesserungsvektor
W = v'*P*v; % Verbesserungsquadratsumme
if round(W,4) ~= round(l'*P*l-n'*xdach,4) %Rechenprobe
    fprintf('Rechenprobe ist falsch!\n')
end
z = size(A); n = z(1); u = z(2); % Freiheitsgrade
var0 = (v'*P*v)/(n-u); % Varianz der Gewichtseinheit
Sigma_xx = var0 * Qxx; % Kovarianzmatrix der Unbekannten
if theo_var == 0 
   fprintf('Sigma0 ist ungueltig.');  % Testgröße
else
    TG = (n-u)*var0/theo_var;
end
TS = chi2inv(0.95,(n-u)); % Testschranke
if TG <= TS %Globaltest
    fprintf('H0 wird angenommen.\n')
else
    fprintf('H0 wird verworfen.\n')
end
end