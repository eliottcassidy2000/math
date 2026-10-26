\\ Galois-orbit traces of newform coefficients at primes p = 3 (mod 4), p <= 820, for 2-power levels:
\\ trivial character, weights 2, 4, 6; chi_8, weights 2, 4, 6; chi_-4 and chi_-8, weights 3, 5.
\\ Output lines: label p Tr(a_p) Tr(|a_p|^2) (the second via complex embeddings, rounded).  mac-mini-2026-10-25-pyth.
default(parisize, 1000000000);
P = 820;
plist = select(p -> p % 4 == 3, primes([3, P]));
specs = [[16,6,1],[32,6,1],[64,6,1],[128,6,1], [16,4,1],[32,4,1],[64,4,1],[128,4,1],[256,4,1], [32,2,1],[64,2,1],[128,2,1],[256,2,1], [8,2,8],[16,2,8],[32,2,8],[64,2,8],[8,4,8],[16,4,8],[32,4,8],[64,4,8],[8,6,8],[16,6,8],[32,6,8], [8,3,-4],[16,3,-4],[32,3,-4],[64,3,-4],[128,3,-4],[8,3,-8],[16,3,-8],[32,3,-8],[64,3,-8],[8,5,-4],[16,5,-4],[32,5,-4],[64,5,-4],[8,5,-8],[16,5,-8],[32,5,-8],[64,5,-8]];
{
for (s = 1, #specs,
  sp = specs[s]; mf = mfinit(sp, 0);
  if (mfdim(mf) == 0, next);
  L = mfeigenbasis(mf); F = mffields(mf);
  for (j = 1, #L,
    c = mfcoefs(L[j], P);
    pol = F[j]; R = if (poldegree(pol) > 0, polroots(pol), [0]);
    lab = Str("N", sp[1], "k", sp[2], "c", sp[3], "o", j, "d", poldegree(pol));
    for (t = 1, #plist,
      p = plist[t]; a = c[p + 1];
      v = if (type(a) == "t_POLMOD", lift(a), a);
      emb = vector(#R, i, if (type(v) == "t_POL", subst(v, variable(v), R[i]), v));
      tr = round(real(sum(i = 1, #emb, emb[i])));
      nr = round(sum(i = 1, #emb, norml2(emb[i])));
      print(lab, " ", p, " ", tr, " ", nr)));
);
}
quit;
