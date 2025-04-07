clear all
close all

addpath(genpath(fileparts(which('DIB2D.m'))))

method = 'ERK2ds';
error = true;
save_res = false;

T = 2.5;
m = 1250;
tol = 5e-5;
tol_phiks = 1e-3;

% experiment 1: T = 2.5, 
% ERK2ds, m = [1250,1500,1750,2000]
% ERK2, m = [1250,1500,1750,2000]
% Lawson2, m = [3750,4500,5250,6000] 
% IMEX2_LU, m = [1250,1500,1750,2000]
% IMEX2_Sylvester, m = [1250,1500,1750,2000]
% IMEX2_Sylvester_eig, m = [1250,1500,1750,2000]
% ode23tb, tol = [5e-5,3e-5,2e-5,1e-5]
% ode23, tol = [9e-4,8e-4,7e-4,6e-4]
%
% experiment 2: T = 100, ERK2ds, m = 50000

options.RelTol = tol;
options.AbsTol = tol;
par = m;
t = linspace(0,T,m+1);

du = 1;
dv = 20;
rho = 25/4;
au1 = 10;
au2 = 1;
au3 = 66;
au4 = 0.5;
av1 = 3;
av2 = 2.5;
av3 = 0.2;
av4 = (av1*(1-au4)*(1-av3+av3*au4))/(au4*(1+av3*au4));
av5 = 1.5;
ue = 0;
ve = au4;
gu = @(t,u,v) rho*(au1*(1-v).*u-au2*u.^3-au3*(v-au4));
gv = @(t,u,v) rho*(av1*(1+av2*u).*(1-v).*(1-av3*(1-v))-av4*v.*(1+av5*u).*(1+av3*v));

dim = 2;

a = [0,0];
b = [20,20];
n = [200,200];

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

rng(123)
U0 = ue+1e-5*rand(n);
V0 = ve+1e-5*rand(n);

N = prod(n);
dgudu = @(t,w) rho*(au1*(1-w(N+1:2*N))-3*au2*w(1:N).^2);
dgudv = @(t,w) rho*(-au1*w(1:N)-au3);
dgvdu = @(t,w) rho*(av1*av2*(1-w(N+1:2*N)).*(1-av3*(1-w(N+1:2*N)))-av4*av5*w(N+1:2*N).*(1+av3*w(N+1:2*N)));
dgvdv = @(t,w) rho*(av1*(1+av2*w(1:N)).*(2*av3-1-2*av3*w(N+1:2*N))-av4*(1+av5*w(1:N)).*(1+2*av3*w(N+1:2*N)));
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
    case 'IMEX2_Sylvester_eig'
        tic
        [U,V] = IMEX2_2D_sylvester_eig(U0,V0,T/m,m,Au,Av,du,dv,gu,gv,n,h);
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

if error && T == 2.5
    load('DIB_2D_Uref.mat')
    Vref = Uref{2};
    Uref = Uref{1};
    normUref = norm(Uref,'fro');
    normVref = norm(Vref,'fro');
    err = norm([norm(U-Uref,'fro')/normUref;norm(V-Vref,'fro')/normVref]);
    fprintf('Error: %.3e\n',err)
    if save_res
        save(['files/results/DIB2Dresults/',method,'/',name],'par','time','err')
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