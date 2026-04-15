% Résolution par un schéma au différences finies explicite, du problème suivant :
%
%  Ks*d/dz(Kr*d/dz(Phi)) = d/dt(Teta)
%
% z = déplacement 1D (vertical)
% h = pression hydraulique
% Phi = h+z ; Charge hydraulique
% Ks = cst ; conductivité hydrolique d'un sol saturé
% Kr(h) = exp(alpha*h) ; conductivité hydrolique d'un sol non saturé
% alpha = cst ;  paramètre dépendant du type de sol
% Teta(h) = tetad + (tetas - tetad)Kr(h) ; teneur en eau
% tetad = cst ;  teneur en eau pour un sol sec
% tetas = cst ;  teneur en eau pour un sol saturé

clear; close all;


Tf = 1; % temps final en jour
N = 100; % Nombre de subdivision de l'intervalle de temps [0;TF]
L = 50; % Limite de l'espace considéré en m
J = 200; % Nombre de subdivision de l'intervalle en espace [0;L]

dz=L/J; % Pas d'espace
z=dz*[1:J-1];

dt=Tf/N;
t=dt*[1:N-1];


h = zeros(J+1,N+1); % Pression hydraulique
hd = -20; % Pression hydraulique d'un sol sec
h(:,1)= hd*ones(J+1,1); % Condition initiale
h_0=hd*ones(1,N); % Condition au bord h(0,t) avec t>0
h(1,2:N+1)=h_0;
h_L=zeros(1,N); % Condition au bord h(L,t) avec t>0
h(J+1,2:N+1)=h_L;

Phi = zeros(J+1,N+1);
Phi(:,1)= h(:,1) + transpose([0,z,L]);
Phi_0=hd*ones(1,N);
Phi(1,2:N+1)=Phi_0;
Phi_L=zeros(1,N) + L*ones(1,N);
Phi(J+1,2:N+1)=Phi_L;


Ks = 0.1 ; % m/Jour
alpha = 0.1 ; % m^-1
Kr = @(h) exp(alpha*h) ;
tetad = 0.15 ;
tetas = 0.45 ;
Teta = @(h) tetad + (tetas - tetad)*Kr(h);
dTeta = @(h) alpha*(tetas - tetad)*Kr(h); % Dérivée de Teta en fonction de h

beta = (Ks*dt)/(dz^2); % Pour alléger l'écriture



temps = cputime;
for n=1:N

    h_n = h(:,n); % Vecteur de la pression en tout point de l'espace au temps n

    % Calcul des vecteurs K_r+, K_r- et dTeta à chaque étape temporelle
    hm_p = 0.5*(h_n(2:J)+h_n(3:J+1)); % 0.5(h_j + h_j+1)
    hm_m = 0.5*(h_n(2:J)+h_n(1:J-1)); % 0.5(h_j + h_j+1)

    Kr_p = Kr(hm_p);
    Kr_m = Kr(hm_m);

    cTeta = 1./dTeta(h_n(2:J));

    % Construction de la matrice tridiagonale du système linéaire
    a = -beta * cTeta .* Kr_m;
    b = 1 + beta * cTeta .* (Kr_m + Kr_p);
    c = -beta * cTeta .* Kr_p;
    A = spdiags([[a(2:end);0], b, [0;c(1:end-1)]],[-1,0,1],J-1,J-1);

    % Construction du second membre du système linéaire (avec conditions aux bords)
    B = Phi(2:J,n);
    B(1) = B(1) + beta*cTeta(1)*Kr_m(1)*Phi(1,n+1);
    B(J-1) = B(J-1) + beta*cTeta(J-1)*Kr_p(J-1)*Phi(J+1,n+1);

    % Résolution du système linéaire par l'algorithme de Thomas
    Phi(2:J,n+1) = A\B;

    % Calcul de h à l'étape n+1
    h(2:J,n+1) = Phi(2:J,n+1) - z';

end

disp(["Temps d'execution :",num2str(cputime-temps)]);

tz=[0,z,L];
t=[0,t,Tf];
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
% Calcul de la solution exacte
nbTermes = 40;
tt = [0.1,0.3,0.5,0.8,1]*Tf;
hex = zeros(J+1,length(tt)); % solution exacte

epsilon = exp(alpha*hd);
c = alpha*(tetas - tetad)/Ks;
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
erreur = norm(hex(:,3) - h(:,round(0.5*N)))/norm(hex(:,3));
disp(["Erreur relative au temps t = ", num2str(t)]); disp(erreur);

figure(2); hold on;
plot(tz,hex(:,3));
plot(tz,h(:,round(0.5*N)));
hold off;
xlabel("Hauteur z");
ylabel("Presion hydraulique h");
title("Différence entre solution analytique et solution trouvée");
legend("Solution exacte","Solution numérique");

figure(3); hold on;
% plot(tz,Teta(h(:,round(0.5*N))),"Linewidth",2);
plot(tz,hex(:,1));
plot(tz,hex(:,2));
plot(tz,hex(:,3));
plot(tz,hex(:,4));
xlabel("Hauteur z");
ylabel("Theta(h)");
legend(["t=",num2str(round(1*Tf)/10)],["t=",num2str(round(3*Tf)/10)],["t=",num2str(round(5*Tf)/10)],["t=",num2str(round(8*Tf)/10)]);
title("Teneur en eau en fonction du temps");

