\\ Galois-orbit traces at primes p = 3 (mod 4), p <= 820, of newforms of levels 2^a * m (m = 3, 9, 5, 7), weights 2, 4, 6,
\\ trivial character.  Output: label p Tr(a_p).  mac-mini-2026-10-25-pyth.
default(parisize, 2000000000);
P = 820;
plist = select(p -> p % 4 == 3, primes([3, P]));
levels = [3,6,12,24,48,96, 9,18,36,72,144, 5,10,20,40,80, 7,14,28,56,112, 15,30,60,21,42,84,45,63,90,126];
{
for (a = 1, #levels, for (kk = 1, 3,
  k = 2 * kk; N = levels[a];
  mf = mfinit([N, k, 1], 0);
  if (mfdim(mf) == 0, next);
  L = mfeigenbasis(mf); F = mffields(mf);
  for (j = 1, #L,
    c = mfcoefs(L[j], P); pol = F[j]; R = if (poldegree(pol) > 0, polroots(pol), [0]);
    lab = Str("N", N, "k", k, "c1o", j, "d", poldegree(pol));
    for (t = 1, #plist,
      p = plist[t]; aa = c[p + 1];
      v = if (type(aa) == "t_POLMOD", lift(aa), aa);
      tr = round(real(sum(i = 1, #R, if (type(v) == "t_POL", subst(v, variable(v), R[i]), v))));
      print(lab, " ", p, " ", tr, " 0")))));
}
quit;
