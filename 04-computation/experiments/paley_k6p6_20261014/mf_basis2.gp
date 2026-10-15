\\ Newforms of 2-power level N = 4..256, quadratic 2-power characters, weights 2..6.
\\ For each Galois orbit: the trace of a_p over the coefficient field, for all odd primes p < B.
B = 1000;
P = primes([3, B]);
{
for(k = 2, 6,
  foreach([4, 8, 16, 32, 64, 128, 256], N,
    foreach([1, -4, 8, -8], D,
      if(N % abs(D) != 0 && D != 1, next);
      if((D < 0) != (k % 2 == 1), next);
      iferr(
        mf = mfinit([N, k, D], 0);
        if(mfdim(mf) == 0, next);
        L = mfeigenbasis(mf);
        for(i = 1, #L,
          co = mfcoefs(L[i], B);
          v = vector(#P, j, my(c = co[P[j] + 1]); if(type(c) == "t_POLMOD", trace(c), c));
          d = if(type(co[2]) == "t_POLMOD", poldegree(co[2].mod), 1);
          print(k, " ", N, " ", D, " ", d, " : ", v)
        ),
      E, next)
    )
  )
);
}
quit
