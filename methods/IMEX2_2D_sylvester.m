function [U,V] = IMEX2_2D_sylvester(U0,V0,tau,m,Au,Av,gu,gv)

Mu{1} = full(speye(length(Au{1}))-tau/2*Au{1});
Mu{2} = full(-tau/2*Au{2});
Mv{1} = full(speye(length(Av{1}))-tau/2*Av{1});
Mv{2} = full(-tau/2*Av{2});

for mu = 1:2
    Au{mu} = full(Au{mu});
    Av{mu} = full(Av{mu});
end

t = 0;

U = U0;
V = V0;

for n = 1:m
    Ui = sylvester(Mu{1},Mu{2}',U+tau/2*gu(t,U,V));
    Vi = sylvester(Mv{1},Mv{2}',V+tau/2*gv(t,U,V));
    U = sylvester(Mu{1},Mu{2}',U+tau/2*(Au{1}*U+U*Au{2}')+tau*gu(t+tau/2,Ui,Vi));
    V = sylvester(Mv{1},Mv{2}',V+tau/2*(Av{1}*V+V*Av{2}')+tau*gv(t+tau/2,Ui,Vi));

    t = t+tau;
end