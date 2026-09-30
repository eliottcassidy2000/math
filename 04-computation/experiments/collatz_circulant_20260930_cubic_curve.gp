E = ellfromeqn(x^3 + y^3 + 1 - x*y);
print("coeffs: ", E);
E = ellinit(E); print("j: ", E.j);
G = ellglobalred(E); print("conductor: ", G[1]);
M = ellminimalmodel(E); print("minimal: [", M.a1, ",", M.a2, ",", M.a3, ",", M.a4, ",", M.a6, "]");
print("torsion: ", elltors(E)[1]);
default(timer,0); print("2-Selmer / rank bounds: ", ellrank(E, 0));
