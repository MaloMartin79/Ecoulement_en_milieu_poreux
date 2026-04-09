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

dz=L/J % Pas d'espace
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
dh_Teta = @(h) alpha*(tetas - tetad)*Kr(h); % Dérivée de Teta en fonction de h



beta = (Ks*dt)/(dz^2); % Pour alléger l'écriture

I = speye(J-1,J-1);
Ap = spdiags([ones(J-1,1),-ones(J-1,1)],[0,1],J-1,J-1);
Am = spdiags([-ones(J-1,1),ones(J-1,1)],[-1,0],J-1,J-1);
F = sparse(J-1,1);


%A = eye(J-1) + alpha*(ht/hx^2)*(diag(-2*ones(1,J-1))+diag(ones(1,J-2),1)+diag(ones(1,J-2),-1));
%A=sparse(A);

for n=1:5

    h_n = h(:,n); % Vecteur de la pression en tout point de l'espace au temps n


    % Calcul des vecteurs K_r+, K_r- et dTeta à chaque étape temporelle
    disp(["Iteration :",num2str(n)]);
    hm_p = 0.5*(h_n(2:J)+h_n(3:J+1)); % 0.5(h_j + h_j+1)
    hm_m = 0.5*(h_n(2:J)+h_n(1:J-1)); % 0.5(h_j + h_j+1)

    Kr_p = Kr(hm_p);
    Kr_p = spdiags(Kr_p,0,J-1,J-1);
    Kr_m = Kr(hm_m);
    Kr_m = spdiags(Kr_m,0,J-1,J-1);

    coef_Teta = 1./dh_Teta(h_n);
    coef_Teta = spdiags(coef_Teta,0,J-1,J-1);


    F(1) = beta*coef_Teta(1)*Kr_m(1)*Phi(1,n);
    F(J-1) = beta*coef_Teta(J+1)*Kr_p(J-1)*Phi(J+1,n);


    Phi(2:J,n+1) = (I - beta*coef_Teta*(Ap*Kr_p + Am*Kr_m))*Phi(2:J,n) + F;

    h(2:J,n+1) = Phi(2:J,n+1) - transpose(z);
end

z=[0,z,L];
t=[0,t,Tf];
[T,Z]=meshgrid(t,z);

size(T)
size(Z)
size(h)
h;
figure(1);
meshc(T,Z,h);
%view([0 90])
%colorbar
xlabel("Temps en secondes");
ylabel("Pression hydraulique");
