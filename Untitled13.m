A = [-1 -2 3 -4 -5; 3 6 -1 4 2; -2 -4 0 -2 0; -2 -4 1 -3 1]
rref(A)

B = [-1, 3, -5; 3, -1, 2; -2 0 0; -2 1 1]
rref([B, eye(4)])

C =  null(A)
D = [-1 -2; 0 1; 1 0; 1 0; 0 0]
b = [-7 4 -1 -1 0]
rref([D b'])

A = [2 1 1; 0 1 1; 0 0 1]
F = A'*A
rref(A)
eig(A'*A)

[P, D] = eig(A'*A)

P =[1 3 3; 4 5 6; 7 8 9]
D = [2 0 0; 0 1 0; 0 0 1]
B = P*D*inv(P)
[P, D] = eig(B)

A = [6 8 -2; -3 -4 1; 9 12 -3]
rref(A)

rref(F)

A = [1/sqrt(5), 0;1/sqrt(2), 0; 2/sqrt(5), 1/sqrt(2)]; 
b =[1, 0; 0, 1; 2, 1]
rref([A b])