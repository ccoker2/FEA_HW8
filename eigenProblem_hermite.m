%{
Christian Coker
27 March 2018
%}
clc, clear, close all


%% Specify Solver Parameters
EI = 20.25;
rhoA = 4.05;
L = 0.3;
f = 0;
beta = 1/4;
gamma = 1/2;

k = EI;


m = 16;
n = 1*m+1;                       % number of nodes
h = L/m;                         % element length

plotEig = [1,2];

%% Specify the Problem
w0t = 0;
wLt = 0;




%% Set Up Problem Domain
xNm = linspace(0,L,n);
U = zeros(n,1);                  % pre-allocate solution
xAn = linspace(0,L,100);         % analytical plotting domain (arbitrary)


%% Element Matrices
% Hermite Shape Functions

% Element Stiffness Matrix
kuu = k * (2/h)^3 * [3/2, -3/2; -3/2, 3/2];                             
kut = k * (2/h)^3 * [3*h/4, 3*h/4; -3*h/4, -3*h/4];                             
ktu = k * (2/h)^3 * [3*h/4, -3*h/4; 3*h/4, -3*h/4];                                                                                    
ktt = k * (2/h)^3 * [h^2/2, h^2/4; h^2/4, h^2/2];

% Element Mass Matrix
muu = 2*k*h   * [26/35, 9/35; 9/35, 26/35];
mut = 2*k*h^2 * [11/105, -13/210; 13/210, -11/105];
mtu = 2*k*h^2 * [11/105, 13/210; -13/210, -11/105];
mtt = 2*k*h^3 * [2/105, -1/70; -1/70, 2/105];
                       
% % Element Force Vector                       
% fu = f * (h/2) *    [  1;    1];                        
% ft = f * (h/2) * h* [1/6; -1/6]; 


%% Global Assembly

kuuGlobal = kAssemble(m, kuu);
kutGlobal = kAssemble(m, kut);
ktuGlobal = kAssemble(m, ktu);
kttGlobal = kAssemble(m, ktt);

muuGlobal = kAssemble(m, muu);
mutGlobal = kAssemble(m, mut);
mtuGlobal = kAssemble(m, mtu);
mttGlobal = kAssemble(m, mtt);

% fuGlobal = fAssemble(m, fu);
% ftGlobal = fAssemble(m, ft);


%% Apply BCs
j = 1;

uIndices = 1:2*n;
dIndices = uIndices;


%% Displacement Conditions
if exist('w0t','var') == 1
    uIndices(1) = 0;
    uKnown(j) = w0t;
    j = j + 1;
end

if exist('wLt','var') == 1
    uIndices(n) = 0;
    uKnown(j) = wLt;
    j = j + 1;
end


uIndices(uIndices==0) = [];

%% Solve
kGlobal = [kuuGlobal, ktuGlobal; kutGlobal, kttGlobal];
mGlobal = [muuGlobal, mtuGlobal; mutGlobal, mttGlobal];
% fGlobal = [fuGlobal; ftGlobal];

dIndices = nonzeros(~ismember(dIndices,uIndices).*dIndices)';


% extract submatrices
uKnown = uKnown';

kUnknown = kGlobal(uIndices, uIndices);
kKnown   = kGlobal(uIndices, dIndices);

mUnknown = mGlobal(uIndices, uIndices);
mKnown   = mGlobal(uIndices, dIndices);

% fUnknown = fGlobal(uIndices,:);

% solve EVP
[vSub,dSub] = eigs(kUnknown,mUnknown,length(plotEig),'smallestabs');

% rebuild
v = zeros(n,length(plotEig));  d = zeros(n,length(plotEig));
    

% re-assemble submatrices
v(uIndices,1:length(plotEig)) = vSub;
d(1:length(plotEig),1:length(plotEig)) = dSub;



u = U(1:n);
theta = U(n+1:end);

xNm = linspace(0,L,n);

disp(num2str(d(d~=0)))

%% Plot the eigenvector

figure(1);
hold on
grid on
title(['Eigenvector for n = 1']);
xlabel('x');
ylabel('u(x)');
plot(xNm, v(1:n,1), '-');

figure(2);
hold on
grid on
title(['Eigenvector for n = 2']);
xlabel('x');
ylabel('u(x)');
plot(xNm, v(1:n,2), '-');
