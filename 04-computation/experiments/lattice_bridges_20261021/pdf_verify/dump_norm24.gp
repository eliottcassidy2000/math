read("Lraw_cand2.gp");
U = qflllgram(G); R = U~*G*U;
i = 1; while(R[i,i] != 24, i++);
v = U[,i]; print("norm = ", v~*G*v);
write("Lraw_cand2_norm24_vector.txt", v~);
quit
