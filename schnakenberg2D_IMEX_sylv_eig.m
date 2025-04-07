clear all
close all

addpath(genpath(fileparts(which('schnakenberg2D.m'))))

T = 0.25; 
m = 6000;

% experiment 1: T = 0.25, m = [3000,4000,5000,6000]

t = linspace(0,T,m+1);

du = 1;
dv = 10;
rho = 1000;
au = 0.1;
av = 0.9;
ue = au+av;
ve = av/(au+av)^2;
gu = @(t,u,v) rho*(au-u+(u.^2).*v);
gv = @(t,u,v) rho*(av-(u.^2).*v);

dim = 2;

a = [0,0];
b = [1,1];
n = [150,150];

for mu = 1:dim
    h{mu} = (b(mu)-a(mu))/(n(mu)-1);
    x{mu} = linspace(a(mu),b(mu),n(mu));

    D2{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/(h{mu}^2),-1:1,n(mu),n(mu));
    D2{mu}(1,2) = 2/h{mu}^2; 
    D2{mu}(n(mu),n(mu)-1) = 2/h{mu}^2;

    Aus{mu} = du*D2{mu};
    Avs{mu} = dv*D2{mu};

    Au{mu} = full(Aus{mu});
    Av{mu} = full(Avs{mu});
end
[X{1:dim}] = ndgrid(x{1:dim});

rng(0)
U0 = ue+1e-5*rand(n);
V0 = ve+1e-5*rand(n);

tic

k = 0:n-1;
j = (0:n-1)';
tau = T/m;
for mu = 1:dim
    eigu{mu} = -du*4/h{mu}^2*sin(pi*k/(2*(n(mu)-1))).^2;
    Pu{mu} = cos(pi*j*k/(n(mu)-1));
    Pu{mu} = Pu{mu}./vecnorm(Pu{mu});
    iPu{mu} = Pu{mu}\eye(n(mu));

    eigv{mu} = -dv*4/h{mu}^2*sin(pi*k/(2*(n(mu)-1))).^2;
    Pv{mu} = cos(pi*j*k/(n(mu)-1));
    Pv{mu} = Pv{mu}./vecnorm(Pv{mu});
    iPv{mu} = Pv{mu}\eye(n(mu));
end
dMu{1} = 1-tau/2*eigu{1};
dMu{2} = -tau/2*eigu{2};
DMu = dMu{1}(:)+dMu{2};

dMv{1} = 1-tau/2*eigv{1};
dMv{2} = -tau/2*eigv{2};
DMv = dMv{1}(:)+dMv{2};

t = 0;

U = U0;
V = V0;

for n = 1:m
    Ui = Pu{1}*((iPu{1}*(U+tau/2*gu(t,U,V))*iPu{2}')./DMu)*Pu{2}';
    Vi = Pv{1}*((iPv{1}*(V+tau/2*gv(t,U,V))*iPv{2}')./DMv)*Pv{2}';
    U = Pu{1}*((iPu{1}*(U+tau/2*(Au{1}*U+U*Au{2}')+tau*gu(t+tau/2,Ui,Vi))*iPu{2}')./DMu)*Pu{2}';
    V = Pv{1}*((iPv{1}*(V+tau/2*(Av{1}*V+V*Av{2}')+tau*gv(t+tau/2,Ui,Vi))*iPv{2}')./DMv)*Pv{2}';

    t = t+tau;
end

toc

load('schnakenberg_2D_Uref.mat')
Vref = Uref{2};
Uref = Uref{1};
normUref = norm(Uref,'fro');
normVref = norm(Vref,'fro');
err = norm([norm(U-Uref,'fro')/normUref;norm(V-Vref,'fro')/normVref]);
fprintf('Error: %.3e\n',err)