function [U,V] = lawson(U0,V0,tau,m,Au,Av,gu,gv)

d = length(size(U0));

for mu = 1:d
    Eu{mu} = expm(tau*Au{mu});
    Ev{mu} = expm(tau*Av{mu});
end
U = U0;
V = V0;

t = 0;
for n = 1:m
    guev = gu(t,U,V);
    gvev = gv(t,U,V);
    Un2 = tucker(U+tau*guev,Eu);
    Vn2 = tucker(V+tau*gvev,Ev);
    U = tucker(U+tau/2*guev,Eu)+tau/2*gu(t,Un2,Vn2);
    V = tucker(V+tau/2*gvev,Ev)+tau/2*gv(t,Un2,Vn2);

    t = t+tau;
end