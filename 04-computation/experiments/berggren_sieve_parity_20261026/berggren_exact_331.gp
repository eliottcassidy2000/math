act(q, t, g) = {
  if (g == 3, return(if (t == q, q, (t + 2) % q)));
  if (t == q, return(2 % q));
  if (t == 0, return(q));
  my(i = lift(Mod(t, q)^(-1)));
  if (g == 1, (2 - i) % q, (2 + i) % q);
}
realroots(q, sg) = {
  my(N = q + 1, M = matrix(N, N), F, out = List());
  for (t = 0, q, for (g = 1, 3, my(s = act(q, t, g)); M[t + 1, s + 1] += sg[g]));
  F = factor(charpoly(M));
  for (k = 1, #F~, my(f = F[k, 1]); if (poldegree(f) > 1, my(rr = polroots(f)); for (j = 1, #rr, if (abs(abs(rr[j]) - sqrt(3)) > 1e-20, listput(out, [f, rr[j]])))));
  out;
}
print("q=331 B: ", realroots(331, [1, 1, 1]));
print("q=251 B_chi: ", realroots(251, [1, -1, 1]));
