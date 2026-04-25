clear; close all;

%% Données
Tf = 5; % temps final en jour
N = 100; % Nombre de subdivision de l'intervalle de temps [0;TF]
L = 30; % Limite de l'espace considéré en m
J = 1000; % Nombre de subdivision de l'intervalle en espace [0;L]

% Sable
Ks = 8.25*1e-5 ; % m/Jour
alpha = 0.145*1e2 ; % m^-1
m = 0.627;
n = 1/(1-m);
thetad = 0.045 ;
thetas = 0.430 ;
%

% Yolo light clay
Ks = 9.22*1e-5 ; % m/Jour
alpha = 0.0335*1e2 ; % m^-1
m = 0.5;
n = 1/(1-m);
thetad = 0.102 ;
thetas = 0.368 ;
%

Se = @(h) (1+(alpha*h).^n).^(-m);
Kr = @(h) sqrt(Se(h)).*(1-(alpha*h).^(n-1).*Se(h)).^2 ;
Theta = @(h) thetad + (thetas - thetad)*Se(h);
dTheta = @(h) - (thetas - thetad)*alpha*m*n*(alpha*h).^(n-1).*(1+(alpha*h).^n).^(-m-1); % Dérivée de theta en fonction de h

%% Résolution en h par une méthode semi-implicite
h = zeros(J+1,N+1); % Pression hydraulique
h(:,1)= -100 * ones(J+1,1); % Condition initiale
h_0= -7.5 * ones(1,N); % Condition au bord h(0,t) avec t>0
h(1,2:N+1)=h_0;
h_L= -100 * ones(1,N); % Condition au bord h(L,t) avec t>0
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

figure(1); hold on;
plot(h(:,10),tz); % Temps : 0.5 jours
plot(h(:,20),tz); % Temps : 1 jours
plot(h(:,40),tz); % Temps : 2 jours
plot(h(:,end),tz); % Temps : 5 jours
axis('ij');
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
ylabel("altitude z");
xlabel("Pression hydraulique h");
legend("0.5 jour","1 jour","2 jour","5 jour",'location','southeast');
