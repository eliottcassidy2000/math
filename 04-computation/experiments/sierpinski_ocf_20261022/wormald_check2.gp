\\ Wormald's non-circulant regular 7-tournament (p. 6): cyclic triples abc, a'b'c'; arcs aa', bb', cc', ad, bd, cd; rest forced by regularity.
\\ vertices: a=1 b=2 c=3 a'=4 b'=5 c'=6 d=7.  Expected sqrt(d/H) = 3088/1287.
arcs = [[1,2],[2,3],[3,1],[4,5],[5,6],[6,4],[1,4],[2,5],[3,6],[1,7],[2,7],[3,7],[7,4],[7,5],[7,6],[4,2],[4,3],[5,1],[5,3],[6,1],[6,2]];
A = matrix(7,7); for(k=1,#arcs, A[arcs[k][1],arcs[k][2]] = 1);
print("outdegrees ", vector(7,i,sum(j=1,7,A[i,j])), "  tournament: ", A + A~ == matrix(7,7,i,j,i!=j));
S = A - A~; P = charpoly(S, t);
print("rho = chi_S(d)/|chi_S(i d)| = ", subst(P,t,7)/abs(subst(P,t,7*I)), "   Wormald 3088/1287 = ", 3088/1287*1.);
print("normal matrix? A*A~ == A~*A : ", A*A~ == A~*A);
quit
