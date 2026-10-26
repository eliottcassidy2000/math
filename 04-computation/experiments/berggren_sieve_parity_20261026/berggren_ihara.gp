\\ Exact spectra of the Berggren digraph on P1(F_q): B = P_U + P_A + P_D (and signed B_chi = P_U - P_A + P_D),
\\ U: t -> 2 - 1/t, A: t -> 2 + 1/t, D: t -> t + 2, points 0..q-1 and infinity (index q).
\\ For each q: factor charpoly over Q; every irreducible factor other than x-3, x-1, x+1 must have all roots of
\\ modulus sqrt(3) (Ihara/Ramanujan-type); report the multiplicities of x-1, x+1, x-3 and the worst modulus deviation.
act(q, t, g) = {
  if (g == 3, return(if (t == q, q, (t + 2) % q)));
  if (t == q, return(2 % q));
  if (t == 0, return(q));
  my(i = lift(Mod(t, q)^(-1)));
  if (g == 1, (2 - i) % q, (2 + i) % q);
}
check(q, sg) = {
  my(N = q + 1, M = matrix(N, N), cp, F, dev = 0., m1 = 0, p1 = 0, m3 = 0, nf = 0);
  for (t = 0, q, for (g = 1, 3, my(s = act(q, t, g)); M[t + 1, s + 1] += sg[g]));
  cp = charpoly(M); F = factor(cp);
  for (k = 1, #F~,
    my(f = F[k, 1], e = F[k, 2]);
    if (f == x - 1, p1 += e; next);
    if (f == x + 1, m1 += e; next);
    if (f == x - 3, m3 += e; next);
    nf += e;
    my(r = polroots(f)); for (j = 1, #r, dev = max(dev, abs(abs(r[j]) - sqrt(3)))));
  [q, p1, m1, m3, nf, dev];
}
{
  forprime (q = 3, 211,
    my(a = check(q, [1, 1, 1]), b = check(q, [1, -1, 1]));
    printf("q=%d  B: mult(x-1)=%d mult(x+1)=%d mult(x-3)=%d other-degree=%d maxdev=%.2e | B_chi: mult(x-1)=%d mult(x+1)=%d other-degree=%d maxdev=%.2e\n",
           q, a[2], a[3], a[4], a[5], a[6], b[2], b[3], b[5], b[6]));
}
