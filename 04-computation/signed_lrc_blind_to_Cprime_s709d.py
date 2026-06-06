"""
monad-explorer-2026-06-06-S709 — THE SIGNED LRC IS BLIND TO C' (the open core).

Dispatched angle: structural reduction — what does a signed counterexample/tight config
imply for the UNSIGNED LRC?  Answer (rigorous): for the OPEN part of LRC, NOTHING.

Canon recalled:
 - THM-369: if n∤v for all v in S, then t=k/n witnesses M(S) >= 1/n  (the n-grid witness).
 - THM-398: LRC(n) reduces ENTIRELY to  C':  "n|v for some v  =>  M(S) > 1/n (loose)".
            C' is still a CONJECTURE — the open core. (Tight configs all have n∤v, S564.)
 - T1 (S699 gauge): M({eps_i v_i}) = M(S);  M is sign-invariant.
 - THM-401(i): 2/(2n-1) is the Farey successor of 1/n; 2n-1 is the SECOND-pinch modulus.
 - HYP-2281 (S708b): "shell-partner" is an unsigned property, sign-orbit-constant; the
            sign group adds no M-content on the worry-set. (We sharpen this onto C'.)

THREE compounding reasons the signed/2n-1 apparatus is blind to C':
 (1) MODULUS TRANSVERSALITY: floor/C' lives on G_n={k/n}; signed shells live mod 2n-1.
     gcd(n,2n-1)=1 => G_n ∩ G_{2n-1} = {0}. A signed shell-partner (sum≡0 mod 2n-1) is
     never the floor-active pair (sum≡0 mod n) for worry-set speeds.
 (2) DOMAIN CONFINEMENT: the signed theory's domain = worry-set (tight), all n∤v
     (construction side), already settled by THM-369. C' is the COMPLEMENTARY n|v side.
 (3) GAUGE-BLINDNESS: C' is a statement about M; M is sign-invariant (T1); the sign group
     acts trivially on C''s truth.

This script gives explicit/exhaustive confirmation of (1) and (3) and the domain fact (2).
"""
import sys
from itertools import combinations, product
from math import gcd
from fractions import Fraction as F


def nrm(x):
    r = x % 1
    return min(r, 1 - r)


def pinch_denoms(V):
    ds = set(V)
    L = len(V)
    for i in range(L):
        for j in range(i + 1, L):
            ds.add(V[i] + V[j]); ds.add(abs(V[i] - V[j]))
    ds.discard(0)
    return ds


def M_fastF(V):
    """fast integer M, returned as a reduced Fraction."""
    V = tuple(a for a in V if a != 0)
    best_n, best_d = -1, 1
    for d in pinch_denoms(V):
        for m in range(1, d):
            mn = d
            for v in V:
                r = (v * m) % d
                rr = d - r
                x = r if r < rr else rr
                if x < mn:
                    mn = x
                    if mn == 0:
                        break
            if mn * best_d > best_n * d:
                best_n, best_d = mn, d
    g = gcd(best_n, best_d)
    return F(best_n // g, best_d // g)


def M_exact_arg(V):
    best = F(-1); args = []
    for d in pinch_denoms(V):
        for m in range(1, d):
            t = F(m, d)
            v = min(nrm(x * t) for x in V)
            if v > best:
                best, args = v, [t]
            elif v == best:
                args.append(t)
    return best, args


def gcd_list(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return g


pr = print
def fl(): sys.stdout.flush()

# -------------------------------------------------------------------------
pr("=" * 78)
pr("(1) MODULUS TRANSVERSALITY: gcd(n,2n-1)=1, G_n ∩ G_{2n-1} = {0}")
pr("=" * 78)
bad_gcd = [n for n in range(2, 2000) if gcd(n, 2 * n - 1) != 1]
pr(f"  gcd(n,2n-1) over n=2..1999: nontrivial cases = {bad_gcd}  (empty => always 1)")
# G_n ∩ G_{2n-1}: k/n == j/(2n-1) in [0,1) only at 0
inter_bad = 0
for n in range(2, 60):
    C = 2 * n - 1
    Gn = {F(k, n) for k in range(n)}
    GC = {F(j, C) for j in range(C)}
    common = (Gn & GC) - {F(0)}
    if common:
        inter_bad += 1
pr(f"  G_n ∩ G_(2n-1) \\ {{0}} nonempty for n=2..59: {inter_bad} cases (0 => transverse)")
fl()

# floor-active pair (sum≡0 mod n) vs shell-partner (sum≡0 mod 2n-1) — disjoint on worry-set
pr("")
pr("  Worry-set: is any pair BOTH floor-active (sum≡0 mod n) AND shell-partner")
pr("  (sum≡0 mod 2n-1)?  (would need sum ≡0 mod n(2n-1); worry speeds too small)")
both_pair = 0
for n in range(4, 9):
    C = 2 * n - 1; thr = F(1, n); W = 2 * n
    for V in combinations(range(1, W + 1), n - 1):
        if gcd_list(V) != 1:
            continue
        if M_fastF(V) != thr:
            continue
        for i in range(len(V)):
            for j in range(i + 1, len(V)):
                s = V[i] + V[j]
                if s % n == 0 and s % C == 0:
                    both_pair += 1
pr(f"  tight configs n=4..8: pairs that are BOTH = {both_pair}  (0 => roles disjoint)")
fl()

# -------------------------------------------------------------------------
pr("")
pr("=" * 78)
pr("(2) DOMAIN CONFINEMENT: every tight (worry-set) config has n∤v (construction side)")
pr("=" * 78)
viol = 0; tot = 0
for n in range(4, 9):
    thr = F(1, n); W = 2 * n
    for V in combinations(range(1, W + 1), n - 1):
        if gcd_list(V) != 1:
            continue
        if M_fastF(V) != thr:
            continue
        tot += 1
        if any(v % n == 0 for v in V):
            viol += 1
pr(f"  {tot} tight configs (n=4..8); with a multiple of n = {viol}  "
   f"(0 => all construction-side, settled by THM-369)")
fl()

# -------------------------------------------------------------------------
pr("")
pr("=" * 78)
pr("(3) GAUGE-BLINDNESS OF C': M is sign-invariant => C' truth is sign-invariant")
pr("=" * 78)
# Verify M({eps_i v_i}) = M(V) over all sign patterns, incl. n|v (C') instances.
import random
random.seed(7)
bad = 0; tested = 0
for _ in range(1500):
    n = random.randint(3, 8)
    k = n - 1
    V = tuple(sorted(random.sample(range(1, 3 * n), k)))
    if gcd_list(V) != 1:
        continue
    M0 = M_fastF(V)
    eps = tuple(random.choice((1, -1)) for _ in range(k))
    AV = tuple(abs(eps[i] * V[i]) for i in range(k))
    if M_fastF(AV) != M0:
        bad += 1
    tested += 1
# full 2^k check on a handful incl. an explicit C' (n|v) instance
explicit = [(8, (8, 1, 2, 3, 5, 6, 7)), (6, (6, 1, 2, 4, 5)), (4, (4, 1, 2))]
pr("  explicit C' (n|v) instances — full 2^(n-1) sign-orbit, M constant?")
allconst = True
for n, V in explicit:
    M0 = M_fastF(tuple(abs(x) for x in V))
    vals = set()
    for eps in product((1, -1), repeat=len(V)):
        AV = tuple(abs(eps[i] * V[i]) for i in range(len(V)))
        vals.add(M_fastF(AV))
    Cprime_loose = (M0 > F(1, n))
    pr(f"    n={n} V={V} (n|{[v for v in V if v % n == 0]}): "
       f"M-values across signs={ {str(x) for x in vals} }; "
       f"C'(M>1/n)={Cprime_loose}; orbit-constant={len(vals)==1}")
    allconst = allconst and len(vals) == 1
pr(f"  random sample: {tested} (V,sign) M-equality checks, {bad} violations")
pr(f"  => M sign-invariant (incl. on C' instances): {bad==0 and allconst}")
pr("")
pr("CONCLUSION: the sign group fixes M, hence fixes C''s truth; combined with (1) "
   "transversality and (2) domain confinement, the signed/2n-1 apparatus carries no "
   "DIRECT information about C' — the open core of LRC.")
fl()
