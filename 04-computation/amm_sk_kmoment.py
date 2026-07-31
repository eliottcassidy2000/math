from mpmath import (mp, mpf, hyp3f2, quad, pi, sqrt, log, atan, ellipk, nstr,
                    pslq, catalan, gamma, asin, atanh)
mp.dps = 120
def S(k): return hyp3f2(mpf(1)/4, mpf(3)/4, mpf(1)/k, 1, 1+mpf(1)/k, 1)

# verify the K-moment formula for several k
print("K-moment formula check  S(k) = 16/(k pi sqrt2) int_0^1 mu^{4/k-1}(2-mu^2)^{-2/k-1/2}K(mu)dmu")
for k in [1,2,3,4,5,6]:
    val = 16/(k*pi*sqrt(2))*quad(lambda m: m**(mpf(4)/k-1)*(2-m**2)**(-mpf(2)/k-mpf(1)/2)*ellipk(m**2), [0,1])
    print(f"  k={k}: formula {nstr(val,25)}  3F2 {nstr(S(k),25)}  "
          f"{'MATCH' if abs(val-S(k))<mpf(10)**-20 else 'MISMATCH'}")

M = quad(lambda m: ellipk(m**2)/(2-m**2), [0,1])
print("\nM = int_0^1 K(mu)/(2-mu^2) dmu =", nstr(M, 40))
print("  check pi*S(4) = 2 sqrt2 M :", nstr(2*sqrt(2)*M, 30), "vs", nstr(pi*S(4), 30))

s2 = sqrt(2); L2 = log(1+s2)
cands = {
 "pi^2": pi**2, "L2^2": L2**2, "pi*L2": pi*L2, "G": catalan, "log2^2": log(2)**2,
 "pi*log2": pi*log(2), "L2*log2": L2*log(2), "pi^2/8": pi**2/8, "1": mpf(1),
}
names = list(cands); vals = [cands[n] for n in names]
r = pslq([M]+vals, tol=mpf(10)**-90, maxcoeff=10**8, maxsteps=200000)
print("\nPSLQ M against", names, "->", r)

# targeted: is M a combination of pi^2 and L2^2 only?
for nm, b in [("{pi^2,L2^2}", [pi**2, L2**2]),
              ("{pi^2,L2^2,pi L2}", [pi**2, L2**2, pi*L2]),
              ("{pi^2,L2^2,G}", [pi**2, L2**2, catalan]),
              ("{pi^2,L2^2,pi L2,G,log2^2}", [pi**2,L2**2,pi*L2,catalan,log(2)**2])]:
    r = pslq([M]+b, tol=mpf(10)**-90, maxcoeff=10**7, maxsteps=200000)
    print(f"  M vs {nm} -> {r}")
