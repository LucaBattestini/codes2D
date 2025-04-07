function [N,PHI1,PHI0] = phi2m(A)
%
% function PHI2 = phi2m(A)
% function [PHI2,PHI1] = phi2m(A)
% function [PHI2,PHI1,PHI0] = phi2m(A)
%
% PHI2 is the approximation of phi_2(A) by the Pade' approximation (9,9),
% PHI1 is the approximation of phi_1(A), and PHI0 is the approximation
% of expm(A)

% Scale A by power of 2 so that its norm is < 2 .
[f,e] = log2(norm(A,1));
s = min(max(0,e-1),1023);
A = A/2^s;
% Pade approximation for phi2(z)
ID = eye(size(A));
a = [1/2,-7/120,4/285,-7/9120,21/258400,-3/1447040,1/8139600,-1/1172102400,...
     1/30474662400,1/6704425728000];
b = [1,-9/20,9/95,-7/570,7/6460,-7/103360,7/2325600,-1/10852800,...
     1/564345600,-1/60949324800];
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
PHI1 = A*N+ID;
PHI0 = A*PHI1+ID;
for i = 1:s
  N = ((PHI0+ID)*N+PHI1)/4;
  PHI1 = (PHI0+ID)*PHI1/2;
  PHI0 = PHI0*PHI0;
end
%!test
%! A = randn(4)+1i*randn(4);
%! [P2,P1,P0] = phi2m(A);
%! E = expm(A);
%! assert(P2,A\(A\(E-eye(4)-A)),256*eps)
%! assert(P1,A\(E-eye(4)),256*eps)
%! assert(P0,E,256*eps)
