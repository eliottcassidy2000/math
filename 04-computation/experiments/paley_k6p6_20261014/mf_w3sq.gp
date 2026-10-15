\\ Weight-3 newforms of 2-power level with odd quadratic character: for each Galois orbit print dim, and for odd primes
\\ p < B the trace of a_p^2 over the coefficient field (and the trace of a_p).
B = 1000;
P = primes([3, B]);
{
foreach([4, 8, 16, 32, 64, 128, 256, 512], N,
  foreach([-4, -8], D,
    if(N % abs(D) != 0, next);
    iferr(
      mf = mfinit([N, 3, D], 0);
      if(mfdim(mf) == 0, next);
      L = mfeigenbasis(mf);
      for(i = 1, #L,
        co = mfcoefs(L[i], B);
        d = if(type(co[2]) == "t_POLMOD", poldegree(co[2].mod), 1);
        v2 = vector(#P, j, my(c = co[P[j] + 1]); if(type(c) == "t_POLMOD", trace(c^2), c^2));
        print(N, " ", D, " ", d, " : ", v2)
      ),
    E, print("ERR ", N, " ", D, " ", E))
  )
);
}
quit
