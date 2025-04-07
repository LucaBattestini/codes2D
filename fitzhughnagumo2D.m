clear all
close all

addpath(genpath(fileparts(which('fitzhughnagumo2D.m'))))

method = 'ode23';
error = true;
save_res = true;

T = 10; 
m = 27500;
tol = 2e-7;
tol_phiks = 1e-3;

% experiment 1: T = 10, 
% ERK2ds, m = [20000,22500,25000,27500]
% ERK2, m = [20000,22500,25000,27500]
% Lawson2, m = [600000,675000,750000,825000] 
% IMEX2_LU, m = [20000,22500,25000,27500]
% IMEX2_Sylvester, m = [20000,22500,25000,27500]
% ode23tb, tol = [5e-7,1e-7,8e-8,6e-8]
% ode23, tol = [8e-7,6e-7,4e-7,2e-7]
%
% experiment 2: T = 50, ERK2ds, m = 30000

options.RelTol = tol;
options.AbsTol = tol;
par = m;
t = linspace(0,T,m+1);

du = 1;
dv = 42.1887;
rho = 65.731;
av1 = 11;
av2 = 0.1;
gu = @(t,u,v) rho*(-u.*(u.^2-1)-v);
gv = @(t,u,v) rho*av1*(u-av2*v);

dim = 2;

a = [0,0];
b = [pi,pi];
n = [100,100];

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
U0 = 1e-3*rand(n);
V0 = 1e-3*rand(n);

N = prod(n);
dgudu = @(t,w) rho*(-3*w(1:N).^2+1);
dgudv = @(t,w) -rho*ones(N,1);
dgvdu = @(t,w) rho*av1*ones(N,1);
dgvdv = @(t,w) -rho*av1*av2*ones(N,1);
g = @(t,w) [gu(t,w(1:N),w(N+1:2*N));gv(t,w(1:N),w(N+1:2*N))];
Ku = kronsum(Au);
Kv = kronsum(Av);
odefun = @(t,w) [Ku*w(1:N);Kv*w(N+1:2*N)]+g(t,w);
options.OutputFcn = @(t,u,flag) myoutfcn(t,u,flag,T);
options.Jacobian = @(t,w) [Ku+spdiags(dgudu(t,w),0,N,N),spdiags(dgudv(t,w),0,N,N);...
    spdiags(dgvdu(t,w),0,N,N),Kv+spdiags(dgvdv(t,w),0,N,N)];

disp(['Method: ',method])

switch method
    case 'ERK2ds'
        tic
        [U,V,Umean,Udiff] = erk2ds(U0,V0,T/m,m,Au,Av,gu,gv,x);
        time = toc;
        name = ['result_',method,'_',num2str(m),'.mat'];
    case 'ERK2'
        tic
        [U,V] = erk2(U0,V0,T/m,m,Au,Av,gu,gv,tol_phiks);
        time = toc;
        name = ['result_',method,'_',num2str(m),'.mat'];
    case 'Lawson2'
        tic
        [U,V] = lawson(U0,V0,T/m,m,Au,Av,gu,gv);
        time = toc;
        name = ['result_',method,'_',num2str(m),'.mat'];
    case 'IMEX2_LU'
        tic
        [u,v] = IMEX2_2D_lu(U0(:),V0(:),T/m,m,Ku,Kv,gu,gv);
        time = toc;
        U = reshape(u,n);
        V = reshape(v,n);
        name = ['result_',method,'_',num2str(m),'.mat'];
    case 'IMEX2_Sylvester'
        tic
        [U,V] = IMEX2_2D_sylvester(U0,V0,T/m,m,Au,Av,gu,gv);
        time = toc;
        name = ['result_',method,'_',num2str(m),'.mat'];
    case 'ode23tb'
        tic
        ode23tb(odefun,[0,T],[U0(:);V0(:)],options);
        time = toc;
        load('sol')
        U = reshape(app(1:N),n);
        V = reshape(app(N+1:2*N),n);
        name = ['result_',method,'_',num2str(tol),'.mat'];
        par = tol;
    case 'ode23'
        tic
        ode23(odefun,[0,T],[U0(:);V0(:)],options);
        time = toc;
        load('sol')
        U = reshape(app(1:N),n);
        V = reshape(app(N+1:2*N),n);
        name = ['result_',method,'_',num2str(tol),'.mat'];
        par = tol;
end
fprintf('Time: %.2f s\n',time)

if error && T == 10
    load('fitzhughnagumo_2D_Uref.mat')
    Vref = Uref{2};
    Uref = Uref{1};
    normUref = norm(Uref,'fro');
    normVref = norm(Vref,'fro');
    err = norm([norm(U-Uref,'fro')/normUref;norm(V-Vref,'fro')/normVref]);
    fprintf('Error: %.3e\n',err)
    if save_res
        save(['files/results/fitzhughnagumo2Dresults/',method,'/',name],'par','time','err')
    end
else
    if strcmp(method,'ERK2ds')
        figure
        plot(t,Umean)
        figure
        semilogy(t(2:m+1),Udiff)
    end
    figure
    pl_erk2ds = surf(X{:},U);
    view(2)
    set(pl_erk2ds,'EdgeColor','none')
    colorbar
    xlim([a(1),b(1)])
    ylim([a(2),b(2)])
end