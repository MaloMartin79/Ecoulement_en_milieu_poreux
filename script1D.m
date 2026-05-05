clear; close all;

%% Données
Tf = 3600; % temps final en jour
N = 3600; % Nombre de subdivision de l'intervalle de temps [0;TF]
L = 30; % Limite de l'espace considéré en m
J = 500; % Nombre de subdivision de l'intervalle en espace [0;L]

% Argile limoneuse
Ks = 5.50*1e-8; % m/s
alpha = 0.005 ; % cm^-1
m = 0.083;
n = 1/(1-m);
thetad = 0.070;
thetas = 0.360;
%

% Sable
Ks = 8.25*1e-5 ; % m/s
alpha = 0.145 ; % cm^-1
m = 0.627;
n = 1/(1-m);
thetad = 0.045 ;
thetas = 0.430 ;
%

% Yolo light clay
Ks = 9.22*1e-5; % m/s
alpha = 0.0335 ; % cm^-1
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
hd = 100;
h(:,1)= -hd * ones(J+1,1); % Condition initiale
h_0= -hd * ones(1,N); % Condition au bord h(0,t) avec t>0
h(1,2:N+1)=h_0;
h_L= -7.5 * ones(1,N); % Condition au bord h(L,t) avec t>0
h(J+1,2:N+1)=h_L;

temps = cputime;
% En argument, h contient déjà les conditions initiales à l'intérieure
[h,Stab]=Semi_implicite(J,L,N,Tf,h,Ks,Kr,dTheta); %

disp(["Temps d'execution :",num2str(cputime-temps)]);

%% Affichage de la solution
dz=L/J; % Pas d'espace
tz=0:dz:L;
dt=Tf/N; % Pas de temps
t=0:dt:Tf;

figure(1); hold on
plot(h(:,2),tz); % Temps : 0.5 jours
plot(h(:,3),tz); % Temps : 0.5 jours
plot(h(:,round(N/10)),tz); % Temps : 0.5 jours
plot(h(:,round(N/4)),tz); % Temps : 1 jours
plot(h(:,round(N/2)),tz); % Temps : 2 jours
plot(h(:,round(3*N/4)),tz); % Temps : 1 jours
plot(h(:,end),tz); % Temps : 5 jours
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
ylabel("altitude z");
xlabel("Pression hydraulique h");
legend("0.5 jour","1 jour","2 jour","5 jour",'location','southeast');

figure(2); hold on
plot(Theta(h(:,10)),tz); % Temps : 0.5 jours
plot(Theta(h(:,30)),tz); % Temps : 1 jours
plot(Theta(h(:,round(N/2))),tz); % Temps : 2 jours
plot(Theta(h(:,end)),tz); % Temps : 5 jours
plot(Theta(h(:,2)),tz);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
ylabel("altitude z");
xlabel("Teneur en eau");
legend("0.5 jour","1 jour","2 jour","5 jour",'location','southeast');
