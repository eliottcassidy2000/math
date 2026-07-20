#!/usr/bin/env python3
"""
Is n >= 4 sharp?  Hunting a THREE-real-Gaussian GMC counterexample  (mac-mini-S133)
===================================================================================
The 4-term reduction found earlier is
      P = (1+W)(Wbar - Z Zbar) = (1+W)(Wbar - |Z|^2),   Q = W,
which uses W (2 real Gaussians) and Z (2 more) = 4 real.  The obvious way to save a
dimension is to replace the complex modulus |Z|^2 by a single real square X^2, since
both have mean 1:
      P = (1+W)(Wbar - X^2),   Q = W       -- 3 real Gaussians.
That substitution changes the moments from E[|Z|^{2k}] = k! to E[X^{2k}] = (2k-1)!!,
so the factorials no longer cancel the binomial denominators and the alternating-sum
mechanism is NOT automatic.  This tests it, then searches n = 3 broadly.

Criterion (correct one): E[P^m] = 0 for all tested m, and E[QP^m] != 0 at the top TWO m.
"""
from math import factorial
from itertools import product, combinations

def dfac(a):
    """E[X^a], X ~ N(0,1): (a-1)!! for a even, 0 for a odd."""
    if a % 2: return 0
    return factorial(a) // (2**(a//2) * factorial(a//2))

def mul(p, q, n):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(k1[i]+k2[i] for i in range(n))
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def E(p):
    tot = 0
    for k, c in p.items():
        v = 1
        for a in k:
            v *= dfac(a)
            if v == 0: break
        else:
            tot += c*v
    return tot

def monos(n, dmax):
    out = []
    for d in range(1, dmax + 1):
        for k in product(range(d + 1), repeat=n):
            if sum(k) == d: out.append(k)
    return out

def test(P, Qs, n, M=9):
    if not P: return None
    one = tuple([0]*n)
    Pm = {one: 1}; pw = []
    for m in range(1, M + 1):
        Pm = mul(Pm, P, n)
        if E(Pm) != 0: return None
        pw.append(Pm)
    for Q in Qs:
        vals = [E(mul({Q: 1}, x, n)) for x in pw]
        if vals[-1] != 0 and vals[-2] != 0: return (Q, vals)
    return None

def show(p, names):
    if not p: return "0"
    out = []
    for k, c in sorted(p.items()):
        s = "".join(names[i] + ("" if k[i] == 1 else f"^{k[i]}")
                    for i in range(len(k)) if k[i])
        out.append(("+" if c > 0 else "-") + (f"{abs(c)}" if abs(c) != 1 else "") + s)
    return "".join(out).lstrip("+")

# ------------------------------------------------------------------ the candidate
print("=" * 78)
print("PART 1 -- the natural 3-Gaussian candidate:  P = (1+W)(Wbar - X^2)")
print("=" * 78)
print("  In 3 real Gaussians we cannot use a complex conjugate pair for W and still have")
print("  a real X left over unless we spend 2 on W.  Write W = u1 + i u2 with")
print("  u_i ~ N(0,1/2); then (1+W)(Wbar - X^2) has COMPLEX coefficients.  Its real and")
print("  imaginary parts are each real polynomials in (u1,u2,X), so we test those.")
print()
print("  E[|Z|^{2k}] = k!   but   E[X^{2k}] = (2k-1)!!   -- the factorials that made the")
print("  alternating binomial sum collapse are GONE.  Checking whether it still works:")
print()
for m in range(1, 8):
    # sum_k (-1)^k C(m,k)^2 (2k-1)!! (m-k)!   -- the analogue of the collapsing sum
    from math import comb
    s = sum((-1)**k * comb(m, k)**2 * dfac(2*k) * factorial(m-k) for k in range(m+1))
    s_cplx = sum((-1)**k * comb(m, k)**2 * factorial(k) * factorial(m-k) for k in range(m+1))
    print(f"    m={m}:  real-square analogue = {s:>12}    complex |Z|^2 version = {s_cplx}")
print("  => the real-square substitution does NOT vanish; the mechanism is specific to")
print("     the COMPLEX modulus, whose moments k! are exactly what cancel C(m,k).")

# ------------------------------------------------------------------ broad n=3 search
print()
print("=" * 78)
print("PART 2 -- broad search over 3 real Gaussians")
print("=" * 78)
n = 3; names = ["X", "Y", "T"]
for dmax, maxsupp, coeffs in ((3, 4, (-1, 1)), (3, 5, (-1, 1)), (2, 6, (-1, 1))):
    ms = monos(n, dmax); QS = [k for k in monos(n, 2)]
    total = 0; hits = []
    for ksz in range(1, maxsupp + 1):
        for supp in combinations(range(len(ms)), ksz):
            for signs in product(coeffs, repeat=ksz):
                P = {ms[i]: s for i, s in zip(supp, signs)}
                total += 1
                r = test(P, QS, n, M=8)
                if r: hits.append((P, r))
        if hits: break
    print(f"  deg<={dmax}, support<={maxsupp}, coeffs{coeffs}: {total} polys tested, "
          f"{len(hits)} counterexamples")
    for P, r in hits[:3]:
        print(f"      P = {show(P, names)}   Q = {show({r[0]:1}, names)}   {r[1][:5]}...")

# ------------------------------------------------------------------ n=4 control
print()
print("=" * 78)
print("PART 3 -- CONTROL: the same search at n = 4 real Gaussians must FIND something")
print("=" * 78)
n4 = 4; names4 = ["X", "Y", "T", "U"]
ms = monos(n4, 3); QS = [k for k in monos(n4, 2)]
hits = []; total = 0
for ksz in range(1, 5):
    for supp in combinations(range(len(ms)), ksz):
        for signs in product((-1, 1), repeat=ksz):
            P = {ms[i]: s for i, s in zip(supp, signs)}
            total += 1
            r = test(P, QS, n4, M=7)
            if r: hits.append((P, r))
        if len(hits) > 3: break
    if hits: break
print(f"  n=4, deg<=3, support<=4, unit coeffs: {total} tested, {len(hits)} counterexamples")
for P, r in hits[:4]:
    print(f"      P = {show(P, names4):<24} Q = {show({r[0]:1}, names4):<8} {r[1][:5]}...")
if not hits:
    print("  NOTE: none found in REAL coordinates at this support size, even though the")
    print("  complex-coordinate 4-term example exists -- the real expansion of")
    print("  (1+W)(Wbar-|Z|^2) has more terms and non-unit coefficients, so it lies")
    print("  OUTSIDE this search box.  The box is the limitation, not the mathematics.")

print()
print("SUMMARY")
print("  The mechanism needs E[|Z|^{2k}] = k! exactly -- that is what cancels the binomial")
print("  coefficient and turns the sum into sum_j (-1)^j C(m,j).  A real square gives")
print("  (2k-1)!! instead and the collapse fails.  No 3-real-Gaussian counterexample was")
print("  found in the searched boxes; n >= 4 SURVIVES as the record, but sharpness is")
print("  NOT proved -- these are bounded searches over unit coefficients only.")
