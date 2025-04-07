function [U,V] = IMEX2_2D_sylvester_eig(U0,V0,tau,m,Au,Av,du,dv,gu,gv,n,h)

for mu = 1:2
    k = 0:n(mu)-1;
    j = (0:n(mu)-1)';

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