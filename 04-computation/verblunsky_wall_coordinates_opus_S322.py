# opus-2026-07-16-S322 -- HYP-6985 / THM-883: THE VERBLUNSKY WALL COORDINATES.
# (1) EXACT: tight-locus measure (uniform on primitive 14th roots): moments =
#     Ramanujan sums c_14(k)/6; rational Gram-Schmidt OPUC; predictions:
#     Phi_6 = 14th cyclotomic, alpha_5 = -1 (termination = tightness).
# (2) EXACT: the half-turn law: primitive-7th-roots measure -> alpha_n flips
#     sign as (-1)^(n+1) under z -> -z (THM-882's transport in coordinates).
# (3) FLOW: alpha_n(delta) of the AP {1..13} safe measure, delta -> 1/14:
#     the blow-up |alpha_5| -> 1 and its rate.
# (4) WALL COORDINATES: loose sets' alpha-vectors at delta = 1/14.
from fractions import Fraction
from math import gcd, pi, cos, sin
import cmath

# ---------- exact rational OPUC from rational real moments ----------
def opuc_exact(moments, N):
    # moments: dict k -> Fraction for |k| <= N (real symmetric: c_{-k} = c_k)
    def c(k): return moments[abs(k)]
    # Gram-Schmidt on 1, z, ..., z^N with <z^i, z^j> = c(i - j)
    # (monic orthogonal polys; rational arithmetic)
    polys = []   # list of coefficient lists (ascending)
    alphas = []
    for n in range(N + 1):
        coef = [Fraction(0)]*(n+1)
        coef[n] = Fraction(1)   # z^n
        for m in range(n):
            pm = polys[m]
            # <z^n, Phi_m> / <Phi_m, Phi_m>
            num = sum(pm[j] * c(n - j) for j in range(len(pm)))
            den = sum(pm[i] * pm[j] * c(i - j)
                      for i in range(len(pm)) for j in range(len(pm)))
            if den == 0: return polys, alphas, m   # rank hit
            r = num / den
            for j in range(len(pm)): coef[j] -= r * pm[j]
        polys.append(coef)
        if n >= 1:
            alphas.append(-coef[0] if True else None)  # alpha_{n-1} = -Phi_n(0) (real)
    return polys, alphas, None

def ramanujan(q, k):
    # c_q(k) = sum over units p mod q of e(pk/q): integer via Mobius/gcd formula
    from math import gcd as g
    d = g(k, q)
    # c_q(k) = mu(q/d) * phi(q) / phi(q/d)
    def mu(n):
        r, p, m = n, 2, 1
        while p*p <= r:
            if r % p == 0:
                r //= p
                if r % p == 0: return 0
                m = -m
            p += 1
        if r > 1: m = -m
        return m
    def phi(n):
        r, p, res = n, 2, n
        while p*p <= r:
            if r % p == 0:
                while r % p == 0: r //= p
                res -= res // p
            p += 1
        if r > 1: res -= res // r
        return res
    return mu(q//d) * phi(q) // phi(q//d)

print("(1) THE TIGHT-LOCUS OPUC (exact):")
mom14 = {k: Fraction(ramanujan(14, k), 6) for k in range(0, 8)}
print(f"   moments c_k = c_14(k)/6: {[str(mom14[k]) for k in range(8)]}")
polys, alphas, rank = opuc_exact(mom14, 6)
print(f"   Phi_6 coefficients (ascending): {[str(x) for x in polys[6]]}")
cyc14 = [1, -1, 1, -1, 1, -1, 1]
print(f"   == 14th cyclotomic z^6-z^5+z^4-z^3+z^2-z+1: "
      f"{[Fraction(c) for c in cyc14] == polys[6]}")
print(f"   Verblunsky alphas (alpha_0..alpha_5): {[str(a) for a in alphas]}")
print(f"   termination |alpha_5| = 1: {abs(alphas[5]) == 1}")

print("\n(2) THE HALF-TURN LAW (primitive 7th roots vs primitive 14th):")
mom7 = {k: Fraction(ramanujan(7, k), 6) for k in range(0, 8)}
p7, a7, _ = opuc_exact(mom7, 6)
flip = [((-1)**(n+1)) * a7[n] for n in range(6)]
print(f"   heptagon alphas:  {[str(a) for a in a7]}")
print(f"   sign-flipped:     {[str(a) for a in flip]}")
print(f"   == tight-locus alphas: {flip == alphas}")

print("\n(3) THE CRITICALITY FLOW alpha(delta) for the AP {1..13}:")
def safe_arcs(S, delta):
    # exact arcs of {t in (0,1): ||v t|| > delta for all v in S}
    ivs = [(Fraction(0), Fraction(1))]
    for v in S:
        bands = []
        for k in range(v):
            lo = (Fraction(k) + delta)/v
            hi = (Fraction(k+1) - delta)/v
            if lo < hi: bands.append((lo, hi))
        new = []
        for (a, b) in ivs:
            for (cc, d) in bands:
                lo, hi = max(a, cc), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def moments_arcs(arcs, K):
    tot = float(sum(b - a for a, b in arcs))
    ms = {}
    for k in range(K + 1):
        if k == 0: ms[0] = 1.0 + 0j; continue
        z = 0j
        for (a, b) in arcs:
            z += (cmath.exp(-2j*pi*k*float(b)) - cmath.exp(-2j*pi*k*float(a)))/(-2j*pi*k)
        ms[k] = z / tot
    return ms

def opuc_float(ms, N):
    polys, alphas = [], []
    def c(k): return ms[abs(k)] if k >= 0 else ms[abs(k)].conjugate()
    for n in range(N + 1):
        coef = [0j]*(n+1); coef[n] = 1.0+0j
        for m in range(n):
            pm = polys[m]
            num = sum(pm[j].conjugate()*0 + pm[j]*0 for j in [0])  # placeholder
            num = sum(pm[j].conjugate() * c(n - j) for j in range(len(pm)))
            den = sum(pm[i].conjugate() * pm[j] * c(i - j)
                      for i in range(len(pm)) for j in range(len(pm)))
            r = num / den
            for j in range(len(pm)): coef[j] -= r * pm[j]
        polys.append(coef)
        if n >= 1: alphas.append(-coef[0].conjugate())
    return alphas

S = list(range(1, 14))
print("   delta        1-|a5|        a5           a0")
for num, den in ((0, 1), (1, 28), (1, 20), (1, 16), (5, 72), (1, 14.0001), (0, 0)):
    if den == 0:
        delta = Fraction(1, 14) - Fraction(1, 100000)
    elif isinstance(den, float):
        delta = Fraction(1, 14) - Fraction(1, 10000)
    else:
        delta = Fraction(num, den)
    arcs = safe_arcs(S, delta)
    if not arcs: continue
    ms = moments_arcs(arcs, 7)
    al = opuc_float(ms, 6)
    print(f"   {float(delta):.6f}  {1-abs(al[5]):.6e}  {al[5]:.4f}  {al[0]:.4f}")

print("\n(4) WALL COORDINATES of loose sets at delta = 1/14:")
for name, S2 in (('GAP {1..7}u{10..16}', list(range(1,8))+list(range(10,17))),
                 ('near-AP {21..32}', list(range(21,33))),
                 ('floor {1..13}-{6}', [v for v in range(1,14) if v != 6])):
    arcs = safe_arcs(S2, Fraction(1, 14))
    if not arcs:
        print(f"   {name}: safe set EMPTY at 1/14")
        continue
    ms = moments_arcs(arcs, 9)
    al = opuc_float(ms, 8)
    print(f"   {name}: {len(arcs)} arcs, mu = {float(sum(b-a for a,b in arcs)):.5f}")
    print(f"      alpha = {[f'{a:.3f}' for a in al]}")
