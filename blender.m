%% ---------------------------Set up variables-----------------------------
clear all
tic;

%---------------------------Set constants----------------------------------
Re = 1000; % Reynolds number
gamma = 1.5; % Set height of cylinder
radius = 1; % Set radius of cylinder

%---------------Set # grid points and interior grid points-----------------
N = 40; % Set number of points in r axis
M = 1.5*N; % Number of points along z axis, z/r aspect ratio = 1.5.

%-------------------Set grid spacing/timestep------------------------------
dz = 1/150;
dr = 1/100;

%dt = 10^-5; % for Re = 1
%dt = 10^-4; % for Re = 100
dt = 10^-3; % for Re = 1000

%% -------------------------Set up matrices--------------------------------
%-------------------------1/r and 1/r^2 matrices---------------------------
Rad = sparse(eye(N,N));
for i = 1:N
    Rad(i,i) = 1/(i*dr);
end
Rad2 = Rad*Rad;

%-----------------------Poisson r difference matrix------------------------
for i = 1:N % Set coefficients for tridiagonal entries
    a(i) = -(2+1/i^2) ;  % Main diag
    j(i) = 1+1/(2*i); % Upper diag
    k(i) = 1-1/(2*i); % Lower diag
end

kk(1:N-1) = k(2:N); % Have to shift k/j vectors due to spdiags indexing
kk(N) = 0;
jj(1) = 0;
jj(2:N) = j(1:N-1);
G = [kk; a; jj]';

% Poisson 2nd-order r difference matrix
PoissonR2 = (1/dr^2)*spdiags(G, -1:1, N, N);

%-------------------General r and z difference matrices--------------------
e_n = sparse(ones(N,1));
% Set up general 1st-order and 2nd-order difference matrices for r 
R1 = (0.5/dr)*spdiags([-e_n zeros(N,1) e_n], -1:1, N, N); 
R2 = (1/dr^2)*spdiags([e_n -2*e_n e_n], -1:1, N, N);

e_m = sparse(ones(M,1));
% Set up general 1st-order and 2nd-order difference matrices for z
Z1 = (0.5/dz)*spdiags([-e_m zeros(M,1) e_m], -1:1, M, M);
Z2 = (1/dz^2)*spdiags([e_m -2*e_m e_m], -1:1, M, M);

%% ----------------Pre-loop steps/Initial conditions-----------------------
% Find eigenvector matrix Z, its inverse Zi, and eigenvalues E of matrix
% PoissonR2 for use inside loop
[Z,e] = eig(full(PoissonR2));
E = eig(full(PoissonR2));
Zi = sparse(inv(Z));

T = cell(1,N); % Create cells to store N arrays for use inside loop
I = sparse(eye(M));
for i = 1:N
    T{1,i} = sparse(Z2+E(i)*I);
end

%---------------------V zero with floor forcing----------------------------
V = sparse(zeros(N,M));
v0 = linspace(0,1,N);
V0(:,1) = v0;

%--------------------------Phi and Eta = 0---------------------------------
Psi = sparse(zeros(N,M));
Eta = sparse(zeros(N,M));

%% ------------------------------LOOP--------------------------------------
%----------------Set up V/Eta for predictor/corrector----------------------
V_k = V;
Eta_k = Eta;

V2 = sparse(zeros(N+2,M+2));
P2 = sparse(zeros(N+2,M+2));
E2 = sparse(zeros(N+2,M+2));

err = 1;
Eta_temp = 0*Eta_k;
V_temp = 0*V_k;

U = sparse(zeros(N,M));

t = 1; % Initiate time step
w = 0; % Initiate counter for pngs

while err >= 10^(-16)
%------------------------------V advance-----------------------------------
    % Predictor/corrector
    % For NON-LINEAR SOLVE  RHS = ...
    RHS = (Rad*Psi*Z1).*(R1*V_k + Rad*V_k) - (Rad*R1*Psi).*(V_k*Z1) + ...
         (1/Re)*(R2*V_k + Rad*R1*V_k - Rad2*V_k + V_k*Z2); 
    % For LINEAR SOLVE RHS = ...
%      RHS = (1/Re)*(R2*V_k + Rad*R1*V_k - Rad2*V_k + V_k*Z2);
     
    V_p = V_k + dt*RHS; %Prediction
    
    % For NON-LINEAR SOLVE RHSp = ...
    RHSp = (Rad*Psi*Z1).*(R1*V_p + Rad*V_p) - (Rad*R1*Psi).*(V_p*Z1) + ...
         (1/Re)*(R2*V_p + Rad*R1*V_p - Rad2*V_p + V_p*Z2);
    % For LINEAR SOLVE RHSp = ...
%      RHSp = (1/Re)*(R2*V_p + Rad*R1*V_p - Rad2*V_p + V_p*Z2);
    
    V_k = V_k + 0.5*dt*(RHS + RHSp); %Correction
    
    V_k(:,1) = v0; %Reset boundary conditions
    V_k(:,end) = 0;
    V_k(1,:) = 0;
    V_k(end,:) = 0;

%----------------------------Eta advance-----------------------------------
    RHS = (Rad*Psi*Z1).*(R1*Eta_k) - (Rad*R1*Psi).*(Eta_k*Z1) - ...
        (Rad2*Eta_k).*(Psi*Z1) + (2*Rad*V_k).*(V_k*Z1) + ...
        (1/Re)*(R2*Eta_k + Rad*R1*Eta_k - Rad2*Eta_k + Eta_k*Z2);
    Eta_p = Eta_k + dt*RHS;
    
    RHSp = (Rad*Psi*Z1).*(R1*Eta_p) - (Rad*R1*Psi).*(Eta_p*Z1) - ...
        (Rad2*Eta_p).*(Psi*Z1) + (2*Rad*V_k).*(V_k*Z1) + ...
        (1/Re)*(R2*Eta_p + Rad*R1*Eta_p - Rad2*Eta_p + Eta_p*Z2);
    Eta_k = Eta_k + 0.5*dt*(RHS + RHSp);

% %----------------------With new Eta, update Psi----------------------------
    F = Eta_k; % Set up forcing matrix
    for i = 1:N
        F(i,:) = -i*dr*F(i,:);
    end

%     % Right multiply transpose(F) by transpose(inv(Z)) to get H on r.h.s
    H = sparse(F'*Zi');
%     
%     % --> (Z2+eI)u = h, where e is  the ith eigenvalue of PoissonR2, u is 
%     % the ith row of U, and h is the ith column of H. Solve for u and 
%     % combine into U.
%     
    for i = 1:N
        U(i,:) = T{1,i}\H(:,i);
    end

%     % Finalize solution Psi = ZU.
    Psi = sparse(Z*U);
% 
% %----------------With updated Psi, update Eta walls------------------------
    E2(2:N+1,2:M+1) = Eta_k; % SWITCH TO FULL DOMAIN, CALCULATE BOUNDARIES
    P2(2:N+1,2:M+1) = Psi;
    
    for j = 1:M+2
        E2(N+2,j) = (1/(2*N*dr^3))*(P2(N,j)-8*P2(N+1,j));
    end
    for i = 1:N+2
       E2(i,M+2) = (1/(2*i*dr*dz^2))*(P2(i,M)-8*P2(i,M+1));
       E2(i,1) = (1/(2*i*dr*dz^2))*(P2(i,3)-8*P2(i,2));
    end
    
   Eta_k = E2(2:N+1,2:M+1); % SWITCH BACK TO INTERIOR
 
   % Saving pngs for movie
%    if mod(t,50)==1
%        
%         w = w + 1;
%         
%         V2(2:N+1,1:M) = V_k;
% %         P2(2:N+1,1:M) = Phi; 
% %         E2(2:N+1,1:M) = Eta_k;
% 
%         image = contourf(V2');
%         fName = ['v_transient_Re1000_',num2str(w,'%06.f'),'.png'];
%         saveas(gcf,fName);
% 
%    end

   if t > 1 %Recalculate error every time until it's small enough
       
       err(t) = abs((normest(Eta_k,2)-normest(Eta_temp,2))/normest(Eta_temp,2)); %
       % Error in terms of Eta

       % err(t) =
       % abs((normest(V_k,2)-normest(V_temp,2))/normest(V_temp,2)); %
       % Error in terms of V
   end
   
   Eta_temp = Eta_k;
   V_temp = V_k;
   
   disp(['timestep =  ' num2str(t) ', err = ' num2str(err(t))]) %Display current timestep and error
   t = t + 1; %Advance timestep
    
end

%%-------------------------Plotting---------------------------------------
%----------------Set domain with V,Eta,Psi interiors-----------------------
V2(2:N+1,1:M) = V_k;
% P2(2:N+1,1:M) = Phi; 
% E2(2:N+1,1:M) = Eta_k;

%--------------------------Plot contour------------------------------------
figure(1)
contourf(V2')
title(['Numerical approximation to non-linear BVP (Re = ' num2str(Re) ')'])
xlabel('r = [0,1]'), ylabel('z = [0,\Gamma]')
colorbar
set(gca,'DataAspectRatio', [1 1.5 1]);
toc
