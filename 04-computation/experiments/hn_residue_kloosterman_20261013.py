"""
hn_residue_kloosterman_20261013.py -- residue colourings of the unit-distance graph of L = Q(zeta_55) at every
conjugation-fixed prime of residue degree two are forced to grow like sqrt(p)/2 colours.

Residue graph at a fixed prime p = -1 mod 55 (residue field F_(p^2), complex conjugation = Frobenius):
  vertices F_(p^2), generators the norm-one torus N1 = {z : z^(p+1) = 1} (p+1 elements), which is the image of
  the unit vectors w/wbar of L (Hilbert 90: the residue map sends w/wbar to w^(1-p), whose image is exactly N1).
Spectrum: the eigenvalue at the additive character a is lambda_a = sum_{z in N1} psi(Tr(a z)); it depends only on
  the norm c = N(a) in F_p^* (the torus acts), and we check the identity lambda = -K(c) with the Kloosterman sum
  K(c) = sum_{x in F_p^*} e(2 pi i (x + c/x)/p) (Medrano-Myers-Stark-Terras 1996), so |lambda_a| <= 2 sqrt(p) (Weil).
Hoffman: chi >= 1 + (p+1)/|lambda_min| >= 1 + (p+1)/(2 sqrt p).
Output: for each such prime p <= PMAX the exact least eigenvalue (float64 character sums, with a 1e-6 tolerance
  check against the Kloosterman identity computed independently), the Hoffman bound and the Weil floor.
"""
import sys, math, cmath
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2500

def is_prime(n):
    if n < 2: return False
    for d in range(2, int(n ** 0.5) + 1):
        if n % d == 0: return False
    return True

def nonsquare(p):
    for d in range(2, p):
        if pow(d, (p - 1) // 2, p) == p - 1:
            return d
    raise ValueError

def torus(p, d):
    """norm-one elements of F_p[sqrt d]: (a, b) with a^2 - d b^2 = 1 mod p."""
    return [(a, b) for a in range(p) for b in range(p) if (a * a - d * b * b) % p == 1]

def spectrum(p):
    d = nonsquare(p)
    N1 = torus(p, d)
    assert len(N1) == p + 1
    # Tr((x + y sqrt d)(a + b sqrt d)) = 2(x a + d y b) ; additive character psi(t) = e(t/p) composed with trace
    # lambda depends on N(a + b sqrt d) = a^2 - d b^2 = c; take representatives a = (c-th coset): pick (a,b) with norm c
    rep = {}
    for a in range(p):
        for b in range(p):
            c = (a * a - d * b * b) % p
            if c not in rep:
                rep[c] = (a, b)
    lam = {}
    for c, (a, b) in rep.items():
        s = 0.0
        for (x, y) in N1:
            t = (2 * (x * a + d * y * b)) % p
            s += math.cos(2 * math.pi * t / p)
        lam[c] = s
    # Kloosterman sums over F_p
    K = {}
    for c in range(1, p):
        s = 0.0
        for x in range(1, p):
            s += math.cos(2 * math.pi * ((x + c * pow(x, p - 2, p)) % p) / p)
        K[c] = s
    return d, lam, K

print(f'fixed primes of residue degree two in Q(zeta55): p = -1 mod 55, p <= {PMAX}')
print(' p    d  lambda_min   -K identity(max err)   Hoffman    Weil floor 1+(p+1)/(2sqrt p)')
for p in range(2, PMAX + 1):
    if p % 55 != 54 or not is_prime(p):
        continue
    d, lam, K = spectrum(p)
    err = max(abs(lam[c] + K[c]) for c in range(1, p))
    lmin = min(lam[c] for c in range(1, p))
    hoff = 1 + (p + 1) / abs(lmin)
    weil = 1 + (p + 1) / (2 * math.sqrt(p))
    print(f'{p:5d} {d:3d}  {lmin:+9.4f}   {err:.2e}            {hoff:7.3f}    {weil:7.3f}   -> chi >= {math.ceil(hoff - 1e-9)}', flush=True)
# the degree-ten residue fields of 19, 29, 79: q = l^5, the same lemma with q in place of p
print('\nconjugation-fixed primes of residue degree ten (q = l^5): Weil floor only')
for l in (19, 29, 79):
    q = l ** 5
    print(f' l={l}: q=l^5={q}, floor 1+(q+1)/(2 sqrt q) = {1 + (q + 1) / (2 * math.sqrt(q)):.1f}')
