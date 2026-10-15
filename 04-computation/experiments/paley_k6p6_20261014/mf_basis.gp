\\ Rational traces of newform Galois orbits a_p, for primes p = 3 mod 4, p < B.
\\ Output: one line per orbit: "k N D dim : a_p list"
B = 1500;
P = select(p -> p % 4 == 3, primes([3, B]));
{
for(k = 2, 6,
  forstep(N = 4, 128, 4,
    foreach([1, -3, -4, -8, 5, 8, 12, -20, -24, 13, -7], D,
      if(N % abs(D) != 0 && D != 1, next);
      if((D < 0) != (k % 2 == 1), next);      \\ odd weight needs odd character
      iferr(
        mf = mfinit([N, k, D], 0);
        if(mfdim(mf) == 0, next);
        L = mfeigenbasis(mf);
        for(i = 1, #L,
          f = L[i];
          co = mfcoefs(f, B);
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
