#!/usr/bin/env python3
"""
MEASURE FLOORS FOR THE RESIDUAL FAMILIES (klein-2026-07-10-S235, HYP-5900,
THM-686 reserved).

(1) THE CONTINUUM QUINTIC BONFERRONI (unconditional): with danger sets
    D_l = {alpha : frac(v_l alpha) not in [1/14, 13/14]},
        mu(S) = 1 - meas(Union D_l) >= B5(S) := 1 - S1 + S2 - S3 + S4 - S5,
    S_k = Sum_{|U|=k} meas(Intersect_{l in U} D_l)  (Bonferroni: an
    inclusion-exclusion truncation ending on +S5 upper-bounds the union).
    All S_k EXACT: meas(Intersect D) = Sum_{V subset U} (-1)^|V| A_inf(V),
    each A_inf a line-sweep rational (S234 machinery, integer-scaled).
    THM-671's discrete B5 (100% covering fire rate) ported to the continuum;
    with THM-685, B5 > 0 certifies the family at every q > Sum(v)/B5.

(2) THE EXACT-mu CENSUS of the exhaustive small-speed covering class
    (primitive covering 13-sets, speeds <= 18 -- the 966 bank's domain):
    per-instance exact mu by one 13-coordinate sweep. The empirical measure
    floor of the class + the argmin anatomy.

(3) THE t=2 DEVIATION LEMMA: |A_inf(a,b) - (6/7)^2| <= (24/7)/max(a,b)
    (Koksma on the exact coprime grid; proof in THM-686) -- verified against
    exact sweeps over all coprime pairs to b = 40 (sharp constant measured).

(4) CHAIN-COARSENED B5 where plain B5 fails: boxeph's maximal doubling
    chains as super-runners (chain-safe measure = A_inf(chain direction) =
    his mu_L -- cross-checked), exact chain-danger intersections via
    union-family sweeps; full-depth I-E over chains = exact mu (validation).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb

WELL = list(range(1, 13)) + [182]
W966 = [1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 17]
GEN = [12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120]
DIL = [20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260]
R = F(6, 7)


def lcm(a, b):
    return a * b // gcd(a, b)


def A_inf(w):
    """Exact line measure of the band box along direction w (ints, any gcd:
    reduced internally). Integer-scaled sweep: scale 2D, D = 14*lcm(w)."""
    g = 0
    for v in w:
        g = gcd(g, v)
    w = [v // g for v in w]
    L = 1
    for v in w:
        L = lcm(L, v)
    D2 = 28 * L  # scale 2D
    pts = {0, D2}
    for v in w:
        s = 2 * L // v  # 2D/(14v)
        for k in range(v):
            pts.add((14 * k + 1) * s)
            pts.add((14 * k + 13) * s)
    pts = sorted(pts)
    good_len = 0
    for x, y in zip(pts, pts[1:]):
        m = x + y  # midpoint at scale 2*D2
        if all(2 * D2 <= 14 * (v * m % (2 * D2)) <= 26 * D2 for v in w):
            good_len += y - x
    return F(good_len, D2)


def exact_B5_and_mu(S, do_mu=True):
    """Exact continuum B5 (and mu). A_inf for every sub-support size <= 5."""
    A = {(): F(1)}
    for t in range(1, 6):
        for U in combinations(range(13), t):
            A[U] = A_inf([S[i] for i in U])
    Sk = [F(1)] + [F(0)] * 5
    for t in range(1, 6):
        for U in combinations(range(13), t):
            m = F(0)  # meas(Intersect_{l in U} D_l) by I-E over sub-supports
            for r in range(t + 1):
                for V in combinations(U, r):
                    m += (-1) ** r * A[V]
            Sk[t] += m
    B5 = sum((-1) ** j * Sk[j] for j in range(6))
    mu = A_inf(S) if do_mu else None
    return B5, mu, Sk


def chains_of(S):
    Sset = set(S)
    ch = []
    for v in sorted(S):
        if v % 2 == 1 or v // 2 not in Sset:
            c = [v]
            while c[-1] * 2 in Sset:
                c.append(c[-1] * 2)
            ch.append(c)
    return ch


def chain_bonferroni(S, depth):
    """Bonferroni at the chain level, exact: super-runner danger = chain
    danger; k-wise intersections via I-E over sub-collections of chains,
    each A_inf on the union family. Valid lower bound for odd truncation
    parity (ending on +S_even ... we use depth = odd count of terms:
    mu >= sum_{j<=depth} (-1)^j S_j with depth odd)."""
    ch = chains_of(S)
    n = len(ch)
    A = {(): F(1)}
    for t in range(1, min(depth, n) + 1):
        for U in combinations(range(n), t):
            fam = [v for i in U for v in ch[i]]
            A[U] = A_inf(fam)
    Sk = [F(1)] + [F(0)] * depth
    for t in range(1, min(depth, n) + 1):
        for U in combinations(range(n), t):
            m = F(0)
            for r in range(t + 1):
                for V in combinations(U, r):
                    m += (-1) ** r * A[V]
            Sk[t] += m
    return sum((-1) ** j * Sk[j] for j in range(min(depth, n) + 1)), ch


def enumerate_bank():
    """Primitive covering 13-sets with speeds <= 18."""
    out = []
    for S in combinations(range(1, 19), 13):
        ok = True
        for q in range(2, 15):
            if not any(v % q == 0 for v in S):
                ok = False
                break
        if not ok:
            continue
        g = 0
        for v in S:
            g = gcd(g, v)
        if g == 1:
            out.append(list(S))
    return out


print("(0) SWEEP CROSS-CHECK vs boxeph's chain ladder (HYP-5853):")
for L, expect in ((2, F(11, 14)), (3, F(5, 7)), (4, F(9, 14)), (5, F(33, 56))):
    a = A_inf([2 ** i for i in range(L)])
    print(f"    A_inf(1,2,..,2^{L-1}) = {a}  (mu_{L} = {expect})  "
          f"{'OK' if a == expect else 'MISMATCH'}")

print("\n(1) EXACT CONTINUUM B5 (iid reference value 0.122148):")
rows = []
for S, name in ((WELL, "DEEP WELL"), (W966, "W966 worst-covering"),
                (GEN, "GEN"), (DIL, "DIL")):
    B5, mu, Sk = exact_B5_and_mu(S)
    sv = sum(S)
    qs = -(-sv * B5.denominator // B5.numerator) if B5 > 0 else None
    print(f"  {name}: B5 = {float(B5):+.6f}  mu = {float(mu):.6f}  "
          f"(S1..S5 = {[round(float(x), 4) for x in Sk[1:]]})")
    print(f"      fire: {'YES -> q* = Sum(v)/B5 = ' + str(qs) if B5 > 0 else 'NO'}"
          f"   [q*_mu = {sum(S) // mu + 1}]")
    rows.append((name, S, B5, mu))

print("\n(2) EXACT-mu CENSUS of the <=18 primitive covering class:")
bank = enumerate_bank()
mus = []
for S in bank:
    mus.append((A_inf(S), tuple(S)))
mus.sort()
print(f"    class size = {len(bank)}; exact mu range: "
      f"min = {mus[0][0]} = {float(mus[0][0]):.6f}  at {list(mus[0][1])}")
print(f"    max = {float(mus[-1][0]):.6f}; median = {float(mus[len(mus)//2][0]):.6f}")
print(f"    CLASS FLOOR (exact, exhaustive): mu >= {mus[0][0]} on the whole class")
worst = list(mus[0][1])
B5w, muw, _ = exact_B5_and_mu(worst)
print(f"    B5 at the class argmin: {float(B5w):+.6f} "
      f"({'fires' if B5w > 0 else 'does NOT fire -> chain rescue below'})")

print("\n(3) t=2 DEVIATION LEMMA |A_inf(a,b) - 36/49| <= (24/7)/b  (b = max):")
worst_c = F(0)
argw = None
for b in range(2, 41):
    for a in range(1, b):
        if gcd(a, b) == 1:
            dev = abs(A_inf([a, b]) - R ** 2)
            if dev * b > worst_c:
                worst_c, argw = dev * b, (a, b)
print(f"    all coprime pairs b <= 40: max |dev|*b = {worst_c} = "
      f"{float(worst_c):.4f} at {argw}  (proved constant 24/7 = 3.43; "
      f"bound holds {float(worst_c) <= 24/7})")

print("\n(4) CHAIN-COARSENED BONFERRONI where plain B5 fails:")
for name, S, B5, mu in rows:
    if B5 > 0:
        continue
    for depth in (3, 5):
        cb, ch = chain_bonferroni(S, depth)
        print(f"  {name}: chains {[len(c) for c in ch]}  "
              f"chain-B{depth} = {float(cb):+.6f} "
              f"{'FIRES' if cb > 0 else 'no'}")
# also chain-B on the census argmin if plain failed
if B5w <= 0:
    for depth in (3, 5):
        cb, ch = chain_bonferroni(worst, depth)
        print(f"  census-argmin: chains {[len(c) for c in ch]}  "
              f"chain-B{depth} = {float(cb):+.6f} "
              f"{'FIRES' if cb > 0 else 'no'}")

print("\n(5) B5 FIRE-RATE on a bank sample (every 20th of the census order):")
fire = tot = 0
worst_fire = None
for _, St in mus[::20]:
    b5, _, _ = exact_B5_and_mu(list(St), do_mu=False)
    tot += 1
    if b5 > 0:
        fire += 1
        if worst_fire is None or b5 < worst_fire[0]:
            worst_fire = (b5, St)
    else:
        print(f"    NO-FIRE: {list(St)}  B5 = {float(b5):+.6f}")
print(f"    fired {fire}/{tot}; smallest firing B5 = "
      f"{float(worst_fire[0]):+.6f} at {list(worst_fire[1])}" if worst_fire
      else "    none fired")
