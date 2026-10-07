\\ generic: V = vector of column vectors (rational) = unit vectors; decide whether the subgroup they generate
\\ admits phi: Gamma -> Z/2 with phi(v)=1 for all v (bipartite) ; if not, return an odd relation.
oddtest(V) = {
  my(K = #V, M = matconcat(V), Dn = denominator(M), H, C, s, ones);
  H = mathnf(Dn*M);
  C = matrix(K, matsize(H)[2]);
  for(j=1,K, my(c = matsolve(H, Dn*V[j])); if(denominator(c)!=1, error("not integral")); for(k=1,#c, C[j,k] = c[k]));
  ones = vector(K, j, 1)~;
  s = matsolvemod(C, 2, ones);
  [type(s)=="t_INT", matsize(H)[2]]
};
\\ find a short explicit odd relation via matkerint on a subset
oddrel(V) = {
  my(K = #V, M = matconcat(V), Dn = denominator(M), Kr, best = 0, bl = 0);
  Kr = matkerint(Dn*M);
  for(j=1, matsize(Kr)[2], my(r = Kr[,j], t = sum(i=1,K,r[i])); if(t % 2, my(L = sum(i=1,K,abs(r[i]))); if(!bl || L < bl, best = r; bl = L)));
  [best, bl]
};
