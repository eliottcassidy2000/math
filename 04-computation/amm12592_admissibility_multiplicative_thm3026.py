"""AMM 12592: ADMISSIBILITY IS MULTIPLICATIVE.

The delta's are coefficients in the basis  B_{d,k}(x) = x^{d-k} (1-x)^k, and a
block is ADMISSIBLE at degree d when
    |delta_k| <= binom(d,k)   and   delta_k = binom(d,k)  mod 2.

LEMMA (M).  If (delta) is admissible at d and (eps) at e, then the product
polynomial has basis coefficients (delta * eps)_kappa = sum_{k+k'=kappa}
delta_k eps_k', and this is ADMISSIBLE AT d+e:
    B_{d,k} B_{e,k'} = x^{d+e-(k+k')} (1-x)^{k+k'} = B_{d+e, k+k'}    (exact)
    |sum| <= sum binom(d,k)binom(e,k') = binom(d+e, kappa)            (Vandermonde)
    sum = sum binom(d,k)binom(e,k') = binom(d+e,kappa) mod 2          (same identity mod 2)
LEMMA (L).  Admissibility LIFTS to any larger degree d' >= d, because
x + (1-x) = 1 gives B_{d,k} = sum_t binom(d'-d,t) B_{d', k+d'-d-t}, and
Vandermonde again returns exactly binom(d',.).
LEMMA (D).  For gamma = 3/5, d_i + d_i' <= d'_{i+i'} ALWAYS, since
floor(a)+floor(b) <= floor(a+b) and 3(R+i)/5 + 3(R+i')/5 = 3(2R+i+i')/5.

So each product block Delta_i Delta_i' IS admissible in the doubled epoch at
its own row.  What the naive squaring q^{2R-1} = (1-x)(sum x^i Delta_i)^2
breaks is only the ROW DISTRIBUTION: row j receives ~j pairs (i,i'), and the
capacities add.  The obstruction to a doubling induction is therefore purely
COMBINATORIAL (redistribution across rows), not archimedean -- the per-block
capacity is exactly preserved.
"""
from math import comb
import amm12592_gamma35_beam_deathstar as beam

def conv(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                r[i + j] += u * v
    return r

def admissible(delta, d):
    if len(delta) - 1 > d: return False
    for k, v in enumerate(delta):
        if abs(v) > comb(d, k): return False
        if (v - comb(d, k)) % 2: return False
    return True

def basis(d, k):
    return beam.basis_poly(d, k)

def as_poly(delta, d):
    out = [0] * (d + 1)
    for k, v in enumerate(delta):
        if v:
            for t, c in enumerate(basis(d, k)): out[t] += v * c
    return out

print("LEMMA (M): product of admissible blocks is admissible at the SUM of degrees.")
import random
random.seed(5)
ok_m = True
for trial in range(400):
    d = random.randint(1, 9); e = random.randint(1, 9)
    dl = [random.choice([v for v in range(-comb(d, k), comb(d, k) + 1)
                         if (v - comb(d, k)) % 2 == 0]) for k in range(d + 1)]
    ep = [random.choice([v for v in range(-comb(e, k), comb(e, k) + 1)
                         if (v - comb(e, k)) % 2 == 0]) for k in range(e + 1)]
    c = conv(dl, ep)
    if not admissible(c, d + e): ok_m = False; print("   FAIL", d, e); break
    # and the basis identity is exact at the polynomial level
    if conv(as_poly(dl, d), as_poly(ep, e)) != as_poly(c, d + e):
        ok_m = False; print("   POLY FAIL", d, e); break
print(f"   400 random admissible pairs (d,e <= 9): all products admissible at d+e "
      f"AND the polynomial identity is exact: {ok_m}")

print("\nLEMMA (L): admissibility lifts to any larger degree.")
ok_l = True
for trial in range(400):
    d = random.randint(1, 8); dp = d + random.randint(0, 5)
    dl = [random.choice([v for v in range(-comb(d, k), comb(d, k) + 1)
                         if (v - comb(d, k)) % 2 == 0]) for k in range(d + 1)]
    # lift: multiply by the admissible block for degree dp-d with delta_k = binom(dp-d,k)
    one = [comb(dp - d, k) for k in range(dp - d + 1)]
    if not admissible(one, dp - d): ok_l = False; break
    c = conv(dl, one)
    if not admissible(c, dp): ok_l = False; break
    def _tr(a):
        a = list(a)
        while a and a[-1] == 0: a.pop()
        return a
    if _tr(as_poly(c, dp)) != _tr(as_poly(dl, d)):
        ok_l = False; break
print(f"   400 random lifts d -> d' (the lifting block is delta_k = binom(d'-d,k), "
      f"which is (x+(1-x))^(d'-d) = 1): {ok_l}")

print("\nLEMMA (D): for gamma=3/5, d_i + d_i' <= d'_{i+i'} for all i,i' and R.")
bad = []
for R in (8, 16, 32, 64, 128):
    d = [(3 * (R + i)) // 5 for i in range(R)]
    dp = [(3 * (2 * R + j)) // 5 for j in range(2 * R)]
    for i in range(R):
        for i2 in range(R):
            if i + i2 < 2 * R and d[i] + d[i2] > dp[i + i2]:
                bad.append((R, i, i2))
print(f"   violations over R = 8..128: {bad if bad else 'NONE'}")
print("\n=> the per-block capacity is EXACTLY preserved by doubling (Vandermonde);")
print("   the only obstruction to a doubling induction is redistributing the")
print("   ~j pairs (i,i') that land in row j.  That is combinatorial, not archimedean.")
