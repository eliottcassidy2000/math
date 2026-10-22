\\ AUDIT cross-check: PARI qfbclassno(-q) and quadclassunit(-q) against the reduced-form counts (a1_hlist_gp.txt, written
\\ from a1_class_number.py), and the Dirichlet sum sum_{a<q/2} (a/q) = (2 - (2/q)) h, for every prime q = 3 mod 4, 7 <= q < 20000.
{
L = readvec("a1_hlist_gp.txt");
bad = 0; cnt = 0;
for(i = 1, #L, q = L[i][1]; h = L[i][2]; if(q < 7, next);
  cnt++;
  h1 = qfbclassno(-q); h2 = quadclassunit(-q).no;
  s = sum(a = 1, (q-1)/2, kronecker(a, q));
  if(h1 != h || h2 != h || s != (2 - kronecker(2, q)) * h, bad++; print("MISMATCH q=", q, " forms=", h, " qfbclassno=", h1, " quadclassunit=", h2, " sum=", s)));
print("primes checked: ", cnt, "; mismatches (reduced forms vs qfbclassno vs quadclassunit vs Dirichlet sum): ", bad);
}
\q
