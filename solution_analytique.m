% Résolution par un schéma au différences finies explicite, du problème suivant :
%
%  Ks*d/dz(Kr*d/dz(Phi)) = d/dt(theta)
%
% z = déplacement 1D (vertical)
% h = pression hydraulique
% Phi = h+z ; Charge hydraulique
% Ks = cst ; conductivité hydrolique d'un sol saturé
% Kr(h) = exp(alpha*h) ; conductivité hydrolique d'un sol non saturé
% alpha = cst ;  paramètre dépendant du type de sol
% Theta(h) = thetad + (thetas - thetad)Kr(h) ; teneur en eau
% thetad = cst ;  teneur en eau pour un sol sec
% thetas = cst ;  teneur en eau pour un sol saturé

clear; close all;

%% Données
Tf = 1; % temps final en jour
N = 100; % Nombre de subdivision de l'intervalle de temps [0;TF]
L = 50; % Limite de l'espace considéré en m
J = 200; % Nombre de subdivision de l'intervalle en espace [0;L]

Ks = 0.1 ; % m/Jour
alpha = 0.1 ; % m^-1
Kr = @(h) exp(alpha*h) ;
thetad = 0.15 ;
thetas = 0.45 ;
Theta = @(h) thetad + (thetas - thetad)*Kr(h);
dTheta = @(h) alpha*(thetas - thetad)*Kr(h); % Dérivée de theta en fonction de h

%% Résolution en h par une méthode semi-implicite
h = zeros(J+1,N+1); % Pression hydraulique
hd = -20; % Pression hydraulique d'un sol sec
h(:,1)= hd*ones(J+1,1); % Condition initiale
h_0=hd*ones(1,N); % Condition au bord h(0,t) avec t>0
h(1,2:N+1)=h_0;
h_L=zeros(1,N); % Condition au bord h(L,t) avec t>0
h(J+1,2:N+1)=h_L;

temps = cputime;
% En argument, h contient déjà les conditions initiales à l'intérieure
[h]=Semi_implicite(J,L,N,Tf,h,Ks,Kr,dTheta); %

disp(["Temps d'execution :",num2str(cputime-temps)]);

%% Affichage de la solution
dz=L/J; % Pas d'espace
tz=0:dz:L;
dt=Tf/N; % Pas de temps
t=0:dt:Tf;
%[T,Z]=meshgrid(t,tz);

figure(1);
plot(tz,h(:,1),"Linewidth",2,tz,h(:,round(0.1*N)),"Linewidth",2,tz,h(:,round(0.3*N)),"Linewidth",2,...
      tz,h(:,round(0.5*N)),"Linewidth",2,tz,h(:,round(0.8*N)),"Linewidth",2,tz,h(:,N),"Linewidth",2);
xlabel("Hauteur z");
ylabel("Presion hydraulique h");
legend(["h(z) au temps t=",num2str(0),"jours"],["h(z) au temps t=",num2str(round(Tf)/10),"jours"],...
       ["h(z) au temps t=",num2str(round(3*Tf)/10),"jours"],["h(z) au temps t=",num2str(round(5*Tf)/10),"jours"],...
       ["h(z) au temps t=",num2str(round(8*Tf)/10),"jours"],["h(z) au temps t=",num2str(Tf),"jours"],'Location','northwest');


%% Solution analytique

% Pour la série, le 50 ième termes est de l'ordre de 1e-20 alors
% que les précédents convergent vers un ordre de 1e-5

% Calcul de la solution exacte
nbTermes = 50; % Nombre de termes pour la série
tt = [0.1,0.3,0.5,0.8,1]*Tf;
hex = zeros(J+1,length(tt)); % solution exacte

epsilon = exp(alpha*hd);
c = alpha*(thetas - thetad)/Ks;
lambda_k = @(k) pi/L*k;
mu_k = @(k) 1/c*(alpha^2/4 + lambda_k(k)^2);

for j=1:length(tt);
  t = tt(j);
  for i=1:length(tz)
    serie=0; z=tz(i);
    for k=1:nbTermes
      serie = serie + (-1)^k*lambda_k(k)/mu_k(k)*sin(lambda_k(k)*z)*exp(-mu_k(k)*t);
    end
    hbar = (1-epsilon)*exp(alpha/2*(L-z))*(sinh(alpha/2*z)/sinh(alpha/2*L) + 2/(L*c)*serie);
    hex(i,j) = 1/alpha*log(hbar+epsilon);
  end
end

% Erreur et visualisation de la solution
erreur(1) = norm(hex(:,1) - h(:,round(0.1*N)))/norm(hex(:,3));
erreur(2) = norm(hex(:,2) - h(:,round(0.3*N)))/norm(hex(:,3));
erreur(3) = norm(hex(:,3) - h(:,round(0.5*N)))/norm(hex(:,3));
erreur(4) = norm(hex(:,4) - h(:,round(0.8*N)))/norm(hex(:,3));
erreur(5) = norm(hex(:,5) - h(:,round(1*N)))/norm(hex(:,3));
disp(["Erreur relative au temps t = ", num2str(tt)]); disp(erreur);

figure(2); hold on;
plot(tz,hex(:,3));
plot(tz,h(:,round(0.5*N)),'o');
hold off;
xlabel("Hauteur z");
ylabel("Presion hydraulique h");
title("Différence entre solution analytique et solution trouvée");
legend("Solution exacte","Solution numérique");

figure(3); hold on;
% plot(tz,Theta(h(:,round(0.5*N))),"Linewidth",2);
plot(tz,Theta(hex(:,1)));
plot(tz,Theta(hex(:,2)));
plot(tz,Theta(hex(:,3)));
plot(tz,Theta(hex(:,4)));
xlabel("Hauteur z");
ylabel("Theta(h)");
legend(["t=",num2str(round(1*Tf)/10)],["t=",num2str(round(3*Tf)/10)],["t=",num2str(round(5*Tf)/10)],["t=",num2str(round(8*Tf)/10)]);
title("Teneur en eau en fonction du temps");

