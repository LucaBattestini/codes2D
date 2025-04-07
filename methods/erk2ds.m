function [U,V,Umean,Udiff] = erk2ds(U0,V0,tau,m,Au,Av,gu,gv,x)

d = length(size(U0));

Umean = NaN(m+1,1);
Udiff = NaN(m,1);

for mu = 1:d
    [P2u{mu},P1u{mu}] = phi2m(tau*Au{mu});
    [P2v{mu},P1v{mu}] = phi2m(tau*Av{mu});
end

U = U0;
V = V0;

mean = trapz(x{1},U,1)/(x{1}(end)-x{1}(1));
for mu = 2:d
    mean = trapz(x{mu},mean,mu)/(x{mu}(end)-x{mu}(1));
end
Umean(1) = mean;

t = 0;
for n = 1:m
    Un2 = U+tau*tucker(kronsumv(U,Au)+gu(t,U,V),P1u);
    Vn2 = V+tau*tucker(kronsumv(V,Av)+gv(t,U,V),P1v);
    Unew = Un2+(2^(d-1)*tau)*tucker(gu(t,Un2,Vn2)-gu(t,U,V),P2u);
    V = Vn2+(2^(d-1)*tau)*tucker(gv(t,Un2,Vn2)-gv(t,U,V),P2v);
    
    mean = trapz(x{1},Unew,1)/(x{1}(end)-x{1}(1));
    for mu = 2:d
        mean = trapz(x{mu},mean,mu)/(x{mu}(end)-x{mu}(1));
    end
    Umean(n+1) = mean;
    Udiff(n) = norm(Unew-U,'fro'); 

    U = Unew;

    t = t+tau;
end