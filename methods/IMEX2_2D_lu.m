function [u,v] = IMEX2_2D_lu(u0,v0,tau,m,Ku,Kv,gu,gv)

N = length(u0);

Mu = speye(N)-tau/2*Ku;
[Lu,Uu,pu,qu] = lu(Mu,'vector');
iqu(qu) = 1:N;
Mv = speye(N)-tau/2*Kv;
[Lv,Uv,pv,qv] = lu(Mv,'vector');
iqv(qv) = 1:N;

t = 0;

u = u0;
v = v0;

for n = 1:m
    ui = systemsol(Lu,Uu,pu,iqu,u+tau/2*gu(t,u,v));
    vi = systemsol(Lv,Uv,pv,iqv,v+tau/2*gv(t,u,v));
    u = systemsol(Lu,Uu,pu,iqu,u+tau/2*Ku*u+tau*gu(t+tau/2,ui,vi));
    v = systemsol(Lv,Uv,pv,iqv,v+tau/2*Kv*v+tau*gv(t+tau/2,ui,vi));

    t = t+tau;
end
end

function x = systemsol(L,U,p,iq,b)
x = U\(L\b(p));
x = x(iq);
end