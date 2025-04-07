function [U,V] = erk2(U0,V0,tau,m,Au,Av,gu,gv,tol_phiks)

U = U0;
V = V0;

Fu = @(t,U,V) kronsumv(U,Au)+gu(t,U,V);
Fv = @(t,U,V) kronsumv(V,Av)+gv(t,U,V);

t = 0;
for n = 1:m
    Uk1 = phiks(tau,Au,Fu(t,U,V),1,tol_phiks);
    Vk1 = phiks(tau,Av,Fv(t,U,V),1,tol_phiks);
    Un2 = U+tau*Uk1{2};
    Vn2 = V+tau*Vk1{2};
    
    Uk2 = phiks(tau,Au,gu(t,Un2,Vn2)-gu(t,U,V),2,tol_phiks);
    Vk2 = phiks(tau,Av,gv(t,Un2,Vn2)-gv(t,U,V),2,tol_phiks);
    U = Un2+tau*Uk2{3};
    V = Vn2+tau*Vk2{3};

    t = t+tau;
end