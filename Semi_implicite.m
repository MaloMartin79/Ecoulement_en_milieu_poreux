function [h]=Semi_implicite(J,L,N,Tf,h,Ks,Kr,dTheta)
% Résolution par un schéma volumes finies implicite, du problème suivant :
%               Ks*d/dz(Kr*d/dz(Phi)) = d/dt(Teta)
%
% z = déplacement 1D (vertical)
% h = pression hydraulique
% Phi = h+z ; Charge hydraulique
% Ks = cst ; conductivité hydrolique d'un sol saturé
% Kr(h) = exp(alpha*h) ; conductivité hydrolique d'un sol non saturé
% alpha = cst ;  paramètre dépendant du type de sol
% Theta(h) = tetad + (tetas - tetad)Kr(h) ; teneur en eau
% thetad = cst ;  teneur en eau pour un sol sec
% thetas = cst ;  teneur en eau pour un sol saturé
%
% On part du principe que le h passé en argument contient déjà les conditions initiales à l'intérieure
% et est présenté sous la forme h(z,t) avec z la variable d'espace et t la variable de temps
%

dz=L/J; % Pas d'espace
z=dz*[1:J-1];
dt=Tf/N; % Pas de temps
t=dt*[1:N-1];

beta = (Ks*dt)/(dz^2); % Pour alléger l'écriture

Phi= h + transpose([0,z,L]);

for n=1:N

    h_n = h(:,n); % Vecteur de la pression en tout point de l'espace au temps n

    % Calcul des vecteurs K_r+, K_r- et dTeta à chaque étape temporelle
    hm_p = 0.5*(h_n(2:J)+h_n(3:J+1)); % 0.5(h_j + h_j+1)
    hm_m = 0.5*(h_n(2:J)+h_n(1:J-1)); % 0.5(h_j + h_j+1)

    Kr_p = Kr(hm_p);
    Kr_m = Kr(hm_m);

    cTheta = 1./dTheta(h_n(2:J));

    % Construction de la matrice tridiagonale du système linéaire
    a = -beta * cTheta .* Kr_m;
    b = 1 + beta * cTheta .* (Kr_m + Kr_p);
    c = -beta * cTheta .* Kr_p;
    A = spdiags([[a(2:end);0], b, [0;c(1:end-1)]],[-1,0,1],J-1,J-1);

    % Construction du second membre du système linéaire (avec conditions aux bords)
    B = Phi(2:J,n);
    B(1) = B(1) + beta*cTheta(1)*Kr_m(1)*Phi(1,n+1);
    B(J-1) = B(J-1) + beta*cTheta(J-1)*Kr_p(J-1)*Phi(J+1,n+1);

    % Résolution du système linéaire par l'algorithme de Thomas
    Phi(2:J,n+1) = A\B;

    % Calcul de h à l'étape n+1
    h(2:J,n+1) = Phi(2:J,n+1) - z';
end

