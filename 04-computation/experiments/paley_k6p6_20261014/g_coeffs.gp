\\ The newform orbit of S_3^new(Gamma_0(32), (-4/.)) used in HYP-9276: dimension, CM test, and a_p for p < 100
\\ (written with the exact Gaussian-integer value, y = i).
mf = mfinit([32, 3, -4], 0);
L = mfeigenbasis(mf);
print("dim S_3^new(32, (-4/.)) = ", mfdim(mf), "; Galois orbits: ", #L);
mfall = mfinit([32, 3, -4], 1);
print("dim S_3(32, (-4/.)) (full cuspidal) = ", mfdim(mfall));
f = L[1];
print("CM? ", mfisCM(f), "; field ", mffields(mf));
co = mfcoefs(f, 100);
{
for(p = 3, 100, if(isprime(p),
  c = co[p + 1]; c = if(type(c) == "t_POLMOD", lift(c), c);
  print(p, " ", polcoef(c, 0, y), " ", polcoef(c, 1, y))));
}
quit
