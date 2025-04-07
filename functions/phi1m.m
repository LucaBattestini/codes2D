function [N,PHI0] = phi1m(A)
%
% function PHI1 = phi1m(A)
% function [PHI1,PHI0] = phi1m(A)
%
% PHI1 is the approximation of phi_1(A) by the Pade' approximation (9,9)
% and PHI0 is the approximation of expm(A)

% Scale A by power of 2 so that its norm is < 2 .
[f,e] = log2(norm(A,1));
s = min(max(0,e-1),1023);
A = A/2^s;
% Pade approximation for phi1(z)
ID = eye(size(A));
a = [1,1/38,2/57,7/7752,7/25840,1/155040,1/1627920,1/84651840,...
     1/3047466240,1/335221286400];
b = [1,-9/19,2/19,-14/969,7/5168,-7/77520,1/232560,-1/7054320,...
     1/338607360,-1/33522128640];
% Paterson-Stockmeyer
A2 = A*A;
A3 = A2*A;
B0 = a(3)*A2+a(2)*A+a(1)*ID;
B1 = a(6)*A2+a(5)*A+a(4)*ID;
B2 = a(9)*A2+a(8)*A+a(7)*ID;
B3 = a(10)*ID;
N = ((B3*A3+B2)*A3+B1)*A3+B0;
B0 = b(3)*A2+b(2)*A+b(1)*ID;
B1 = b(6)*A2+b(5)*A+b(4)*ID;
B2 = b(9)*A2+b(8)*A+b(7)*ID;
B3 = b(10)*ID;
D = ((B3*A3+B2)*A3+B1)*A3+B0;
N = full(D\N);
% Undo scaling by repeated squaring
PHI0 = A*N+ID;
for i = 1:s
  N = (PHI0+ID)*N/2;
  PHI0 = PHI0*PHI0;
end
%!test
%! A = randn(4)+1i*randn(4);
%! [P1,P0] = phi1m(A);
%! E = expm(A);
%! assert(P1,A\(E-eye(4)),256*eps)
%! assert(P0,E,256*eps)
