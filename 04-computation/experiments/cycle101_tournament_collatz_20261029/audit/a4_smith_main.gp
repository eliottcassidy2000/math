\\ Smith forms of the audit's own run-observer matrices (written by a3_runobs.py into audit/mats/).
\\ For each matrix: the elementary divisors, then for each prime p dividing det, the multiplicity of each
\\ exponent e (Z/p^e summands).  Usage: gp -q a4_smith.gp
default(parisizemax, 4000000000)
localprof(D) = {
  my(F, P = Map(), out = List(), pr);
  for (i = 1, #D, if (D[i] != 1 && D[i] != 0,
      F = factor(D[i], 10^8);
      for (j = 1, #F~, pr = F[j,1];
        my(key = [pr, F[j,2]], c = 0);
        if (mapisdefined(P, key, &c), mapput(P, key, c + 1), mapput(P, key, 1)))));
  Mat(P)
};
runone(name) = {
  my(D, t);
  read(Str("mats/", name, ".gp"));
  t = getabstime();
  D = matsnf(M);
  print("== ", name, " (", #M, "x", #M, "), matsnf in ", getabstime() - t, " ms");
  print("   det = ", factor(matdet(M), 10^8));
  print("   number of non-unit elementary divisors: ", sum(i = 1, #D, D[i] != 1));
  print("   [prime, exponent] -> multiplicity:");
  my(P = localprof(D));
  for (i = 1, #P~, print("     ", P[i,1], " : ", P[i,2]));
};
foreach(["P7_k2_M", "P7_k2_B0", "NQR7_k2_M", "R7_k2_M", "R7_k2_B0", "P11_k2_M", "P11_k2_B0", "P7_k3_B0"], n, runone(n));
quit;
