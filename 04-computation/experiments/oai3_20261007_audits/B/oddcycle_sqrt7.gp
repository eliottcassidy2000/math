read("oddcycle_lib.gp");
w7 = quadgen(28);
V = List();
addv(A, C) = {my(den=A^2+C^2, x, y); if(den != 0, x=(A^2-C^2)/den; y=2*A*C/den; listput(V, [real(x), imag(x), real(y), imag(y)]~));}
{for(a=-3,3, for(b=-1,1, for(c=-3,3, for(d=-1,1, addv(a+b*w7, c+d*w7)))));}
V = Vec(Set(Vec(V))); K = #V;
print("Q(sqrt7)^2: unit vectors generated: ", K);
ok = 1; for(j=1,K, my(v=V[j], x=v[1]+v[2]*w7, y=v[3]+v[4]*w7); if(x^2+y^2 != 1, ok=0)); print("all unit: ", ok);
r = oddtest(V); print("no bipartition functional (odd closed walk exists): ", r[1], "   rank ", r[2]);
\\ small explicit odd relation using few vectors: try subsets of the 40 lowest-height vectors
ht(v) = vecmax(apply(t->max(abs(numerator(t)),denominator(t)), Vec(v)));
W = vecsort(V, (u,v)->ht(u)-ht(v));
for(n=8, 60, my(S = W[1..n]); if(oddtest(S)[1], my(o = oddrel(S)); print("odd relation among the ", n, " lowest-height vectors; closed walk length ", o[2]); for(i=1,n, if(o[1][i], print("   ", o[1][i], " x (", S[i][1], " + ", S[i][2], " r7, ", S[i][3], " + ", S[i][4], " r7)"))); print("   sum check: ", matconcat(S)*o[1]); break));
