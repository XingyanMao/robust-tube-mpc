function theta = tightenconstraints( A, B, K, W, N)
% theta = tightenconstraints(Z, A, B, K, W, N)
%
% Tightens state and input constraints for a box W.
%
% Z should be a struct with fields G, H, and psi to define the feasible set Z as
%
%    Gx + Hu <= psi
%
% A, B, and K should be the system model and controller gain. W should be a
% vector defining the maximum absolute values for w. The system evolves as
%
%    x^+ = Ax + Bu + w,     u = v + Kx,     -W <= w <= W
%
% This means that W must be a symmetric box in R^n (we make this restriction
% because optimization becomes particularly easy).
%
% Returns theta such that the tighter constraints
%
%    Gx + H(v + Kx) <= psi - theta
%
% are satisfied for any N realizations of W.

G=[1,0;-1,0;0,1;0,-1;0,0;0,0];
H=[0;0;0;0;1;-1];
C = G + H*K;
A = A + B*K;
theta = zeros(size(size(C, 1)));
Ak = eye(size(A));
W = abs(W);
for i = 0:N
    theta = theta + abs(C*Ak)*W;
    Ak = A*Ak;
end

end%function
