function [h]=Semi_implicite2D(I,a,J,L,N,Tf,h,Ks,Kr,dTheta)
% Résolution par un schéma volumes finies implicite, du problème suivant :
%               Ks*d/dz(Kr*d/dz(Phi)) = d/dt(Teta)
%
% z = déplacement vertical
% x = déplacement horizontal
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


dx=a/I; % Pas d'espace horizontal
dz=L/J; % Pas d'espace vertical

if dx ~= dz
  disp("Le pas d'epace doit être homogène suivant les deux dimmensions");
  return
end

z=dz*[1:J-1];
dt=Tf/N; % Pas de temps
t=dt*[1:N-1];

beta = (Ks*dt)/(dz^2); % Pour aléger l'écriture

Phi = zeros((J+1)*(I+1),N+1);
for i=1:I+1
  Phi(i:I+1:J*(I+1)+i,:)= h(i:I+1:J*(I+1)+i,:) + transpose([0,z,L]);
end


for n=1:N

    % Matrice (grille) de la pression en tout point de l'espace au temps n
    H_n = reshape(h(:,n), I+1, J+1);


    % Calcul des matrices K_r+, K_r- selon x à chaque étape temporelle
    hm_px = 0.5*(H_n(2:I,2:J) + H_n(3:I+1,2:J));
    hm_mx = 0.5*(H_n(2:I,2:J) + H_n(1:I-1,2:J));

    Kr_px = Kr(hm_px);
    Kr_mx = Kr(hm_mx);

    % Calcul des matrices K_r+, K_r- selon z à chaque étape temporelle
    hm_pz = 0.5*(H_n(2:I,2:J) + H_n(2:I,3:J+1));
    hm_mz = 0.5*(H_n(2:I,2:J) + H_n(2:I,1:J-1));

    Kr_pz = Kr(hm_pz);
    Kr_mz = Kr(hm_mz);

    % Calcul du vecteur des coéficients suivant dTeta à chaque étape temporelle
    cTheta = 1./dTheta(H_n(2:I,2:J));

    % Passage en vecteur
    Kr_px = Kr_px(:);
    Kr_mx = Kr_mx(:);
    Kr_pz = Kr_pz(:);
    Kr_mz = Kr_mz(:);
    cTheta = cTheta(:);


    % Construction de la matrice A tridiagonale par bloc du système linéaire
    % La matrice A est de taille (I-1)*(J-1) avec des blocs de taille I-1

    % Construction des blocs sous et sur diagonale
    ep = -beta * cTheta .* Kr_pz;
    em = -beta * cTheta .* Kr_mz;

    % Construction du bloc central
    b = 1 + beta * cTheta .* (Kr_mx + Kr_px + Kr_mz + Kr_pz); % Diagonale
    a = -beta * cTheta .* Kr_mx; % Sous diagonale
    c = -beta * cTheta .* Kr_px; % Sur diagonale




    emA = [em(I:end); zeros(I-1,1)], [a(2:end); 0];
    aA = [a(2:end); 0];
    cA = [0; c(1:end-1)];
    epA = [zeros(I-1,1); ep(1:end-I+1)];

    % Coupures entre chaques blocs
    for j=1:J-2
      s = j*(I-1);
      aA(s) = 0;
      cA(s+1) = 0;
    end

    A = spdiags([ emA, aA, b, cA, epA ], [-(I-1), -1, 0, 1, (I-1)], (J-1)*(I-1),(J-1)*(I-1));
    spy(A);


    % Construction du second membre du système linéaire (avec conditions aux bords)

    Phi_n = reshape(Phi(:,n), I+1, J+1);
    Bmat = Phi_n(2:I,2:J);

    % Forme matricielle
    amat = reshape(a, I-1, J-1);
    cmat = reshape(c, I-1, J-1);
    emmat = reshape(em, I-1, J-1);
    epmat = reshape(ep, I-1, J-1);

    % Ajout des conditions du bord inférieur
    Bmat(:,1) = Bmat(:,1) - emmat(:,1) .* Phi_n(2:I,1);

    % Ajout des conditions du bord supérieur
    Bmat(:,end) = Bmat(:,end) - epmat(:,end) .* Phi_n(2:I,J+1);

    % Ajout des conditions du bord gauche
    Bmat(1,:) = Bmat(1,:) - amat(1,:) .* Phi_n(1,2:J);

    % Ajout des conditions du bord droit
    Bmat(end,:) = Bmat(end,:) - cmat(end,:) .* Phi_n(I+1,2:J);

    % Vectorisation
    B = Bmat(:);


    % Résolution du système linéaire
    U = A\B;
    U = reshape(U, I-1, J-1);

    Phi_np = reshape(Phi(:,n+1), I+1, J+1);
    Phi_np(2:I,2:J) = U;

    % Mise à jour de h
    H_np = Phi_np;

    for j = 1:J+1
       H_np(:,j) = Phi_np(:,j) - (j-1)*dz;
    end

    % Revectorisation
    Phi(:,n+1) = Phi_np(:);
    h(:,n+1) = H_np(:);

end

