\\ Paley vs rotational: B_q = #{a in [1,(q-1)/2] : a non-residue mod q} equals ((q-1)/2 - (2 - (2/q)) h(-q))/2 for primes q = 3 mod 4, q > 3
\\ (Dirichlet's class number formula).  Checks all such q < 20000; lists B_q <= 1.  mac-mini-2026-10-21-fanopsp.
{
bad = 0; cnt = 0; small = List();
forprime(q = 7, 20000, if(q % 4 == 3,
  m = (q-1)/2; B = sum(a = 1, m, kronecker(a, q) == -1);
  h = qfbclassno(-q);
  pred = (m - (2 - kronecker(2, q)) * h) / 2;
  cnt++; if(B != pred, bad++; print("MISMATCH q=", q, " B=", B, " pred=", pred));
  if(B <= 1, listput(small, [q, B, h]))));
print("primes q = 3 mod 4, 7 <= q < 20000: ", cnt, "; formula mismatches: ", bad);
print("B_q <= 1 at: ", small);
print("Heegner primes q = 3 mod 4: ", [[q, sum(a=1,(q-1)/2, kronecker(a,q)==-1), qfbclassno(-q)] | q <- [7,11,19,43,67,163]]);
}
\q
