F = [-4 2 7; 3 -3 -7; -12 2 15]
[P D] = eig(F)
D = [2 0 0; 0 1 0; 0 0 1]
G = P*D*inv(P)
eig(G)

A = [2 3; 3 2]
[P D] = eig(A)