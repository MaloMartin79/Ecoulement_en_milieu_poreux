% Résolution par un schéma au différences finies explicite, du problème suivant :
%
%  Ks*d/dz(Kr*d/dz(Phi)) = d/dt(theta)
%
% x = déplacement horizontal
% z = déplacement vertical
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
a = 30; % Largeur de l'espace considéré en m
I = 20; % Nombre de subdivision de l'intervalle en espace horizontal [0;a]
L = 30; % Hauteur de l'espace considéré en m
J = 20; % Nombre de subdivision de l'intervalle en espace vertical [0;L]
nbn = (J+1)*(I+1); % Nb total de noeuds du maillage en espace

Ks = 0.1 ; % m/Jour
alpha = 0.1 ; % m^-1
Kr = @(h) exp(alpha*h) ; %ones(length(h),1);
thetad = 0.15 ;
thetas = 0.45 ;
Theta = @(h) thetad + (thetas - thetad)*Kr(h);
dTheta = @(h) alpha*(thetas - thetad)*Kr(h); % Dérivée de theta en fonction de h
hd = -20; % Pression hydraulique d'un sol sec
eps = exp(alpha*hd);


%% Résolution en h par une méthode semi-implicite
% On numérote les noeuds du maillage de l'espace de gauche à droite puis bas vers le haut
% i.e x croisant puis z croissant
h = zeros(nbn,N+1); % Pression hydraulique
h(:,1)= hd*ones(nbn,1); % Condition initiale h(x,z,0)
h(1:I+1,2:N+1)=hd*ones(I+1,N); % Condition au bord h(x,0,t) avec t>0
h(1:I+1:(I+1)*J+1,2:N+1)=hd*ones(I+1,N); % Condition au bord h(0,z,t) avec t>0
h(I+1:I+1:nbn,2:N+1)=hd*ones(I+1,N); % Condition au bord h(a,z,t) avec t>0
dx = a/I;
x = dx*[1:I-1]';
h_L=log(eps + (1-eps)*(0.75*sin(pi*x/a) - 0.25*sin(3*pi*x/a)))/alpha;
h((I+1)*J+2:nbn-1,2:N+1)=h_L .* ones(length(h_L),N); % Condition au bord h(x,L,t) avec t>0


%% Affichage de la solution
x = x';
tx=[0,x,a]; % Espace horizontal
dz=L/J; % Pas d'espace doit être égal pour les deux dimensions en espace
tz=0:dz:L; % Discrétisation de l'espace vertical
dt=Tf/N; % Pas de temps
tt=0:dt:Tf;
%[T,Z]=meshgrid(t,tz);


h_mid1 = h(:,end);
[X,Z]=meshgrid(tx,tz);
h_mid1=reshape(h_mid1,I+1,J+1);
figure(1);
meshc(X,Z,h_mid1');
xlabel("x")
ylabel("z")
zlabel("h")
title("Pression hydraulique au bout d'une journée")


temps = cputime;
% En argument, h contient déjà les conditions initiales à l'intérieure
[h]=Semi_implicite2D3(I,a,J,L,N,Tf,h,Ks,Kr,dTheta);

disp(["Temps d'execution :",num2str(cputime-temps)]);

%% Affichage de la solution


h_mid = h(:,end);
[X,Z]=meshgrid(tx,tz);
h_mid=reshape(h_mid,I+1,J+1);
figure(2);
meshc(X,Z,h_mid');
xlabel("x")
ylabel("z")
zlabel("h")
title("Pression hydraulique à la mi-journée")





