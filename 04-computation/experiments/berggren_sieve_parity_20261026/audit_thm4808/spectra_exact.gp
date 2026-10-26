\\ Independent EXACT spectral audit of THM-4808 (G) and HYP-9418 (4)-(5) (own code; no lane script used).
\\ Points of P^1(F_q): 0..q-1 and q = infinity.  V = {1, iota, S, S.iota}, chi = (1,-1,1,-1) on V.
\\ A[O,O']    = #{v in V : D(v s_O) in O'}                         (weighted adjacency of G_q, row sums 4)
\\ Achi[O,O'] = sum_{v : D(v s_O) = u s_O'} chi(v) chi(u)          (twisted adjacency, chi-equivariant functions)
\\ By the (proved) Ihara-Bass identity (checked exactly in g_structure.py), every eigenvalue of B other than
\\ +-1 is a root of l^2 - mu l + 3 with mu in spec(A), mu = 4 giving l in {3, 1}; |l| = sqrt 3 <=> |mu| <= 2 sqrt 3.
\\ Exact test: F1 = charpoly(A)/(x-4); R(y) with R(x^2) = +-F1(x)F1(-x); #{mu : |mu| >= 2 sqrt 3} = #roots of R in [12, oo).
inv1(q,t) = if(t==q, 0, if(t==0, q, lift(Mod(t,q)^-1)));
neg1(q,t) = if(t==q, q, (-t)%q);
Dact(q,t) = if(t==q, q, (t+2)%q);
Sact(q,t) = neg1(q, inv1(q,t));
Vact(q,t,k) = if(k==1, t, if(k==2, neg1(q,t), if(k==3, Sact(q,t), Sact(q, neg1(q,t)))));
CHI = [1,-1,1,-1];

build(q) = {
  my(oid = vector(q+1), reps = List(), nO = 0, A, Ac, alive, pos, nA);
  for (t = 0, q, if (oid[t+1] == 0, nO++; my(o = vector(4, k, Vact(q,t,k))); for (k = 1, 4, oid[o[k]+1] = nO);
       listput(reps, vecmin(o))));
  A = matrix(nO, nO);
  for (a = 1, nO, my(s = reps[a]); for (k = 1, 4, my(p = Dact(q, Vact(q, s, k))); A[a, oid[p+1]] += 1));
  \\ twisted
  alive = List(); pos = vector(nO);
  for (a = 1, nO, my(s = reps[a], good = 1); for (k = 1, 4, if (Vact(q,s,k) == s && CHI[k] == -1, good = 0));
       if (good, listput(alive, a); pos[a] = #alive));
  nA = #alive; Ac = matrix(nA, nA);
  for (i = 1, nA, my(a = alive[i], s = reps[a]);
    for (k = 1, 4, my(p = Dact(q, Vact(q, s, k)), b = oid[p+1]);
      if (pos[b], my(u = 0); for (j = 1, 4, if (Vact(q, reps[b], j) == p, u = j));
        Ac[i, pos[b]] += CHI[k] * CHI[u])));
  [A, Ac];
}

\\ number of roots mu of F with |mu| >= 2 sqrt 3 (distinct values), and whether some |mu| = 2 sqrt 3 exactly
outside(F) = {
  my(G = F * subst(F, x, -x), R, Rs);
  R = substpol(G, x^2, y); R = subst(R, y, x);
  Rs = R / gcd(R, R');
  [polsturm(Rs, [12, +oo]), subst(R, x, 12) == 0];
}

check(q) = {
  my(M = build(q), A = M[1], Ac = M[2], F, F1, FA, oA, oC, r4);
  F = charpoly(A, x, 4);
  r4 = valuation(F, x - 4);                \\ multiplicity of the trivial eigenvalue
  F1 = F / (x - 4)^r4;
  oA = outside(F1);
  FA = if (matsize(Ac)[1], charpoly(Ac, x, 4), 1);
  oC = if (poldegree(FA) > 0, outside(FA), [0, 0]);
  [q, matsize(A)[1], r4, subst(F, x, -4) == 0, oA, matsize(Ac)[1], subst(FA, x, 4) == 0 || subst(FA, x, -4) == 0, oC];
}

\\ the factor carrying the exceptional root, its degree, and the B-side root lambda = (mu +- sqrt(mu^2-12))/2
exception(q, which) = {
  my(M = build(q), A = M[which], F = charpoly(A, x), fa, out = List());
  if (which == 1, F = F / (x - 4));
  fa = factor(F);
  for (k = 1, #fa~, my(f = fa[k,1]);
    if (outside(f)[1] > 0,
      my(rr = polrootsreal(f), big = [r | r <- rr, abs(r) > 2*sqrt(3)]);
      for (j = 1, #big, my(mu = big[j], l1 = (mu + sign(mu)*sqrt(mu^2 - 12))/2, l2 = 3/l1);
        listput(out, [poldegree(f), fa[k,2], mu, l1, l2]))));
  out;
}
