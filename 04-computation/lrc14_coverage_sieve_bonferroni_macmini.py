"""
REFRAMING C: Coverage sieve over N=7 sectors + Bonferroni/Brun truncation.

p0(E) = meas{ x in [0,1) : orbit {frac(e*x): e in E} hits ALL 7 sectors }.

COVERAGE SIEVE (over the 7 sectors):
  Let A_j = {x : sector j is NEVER hit by the orbit} = {x : no e in E has frac(e*x) in sector j}.
  "hits all 7" = complement of union_j A_j.
  p0(E) = 1 - meas(union_j A_j) = sum_{S subset {0..6}} (-1)^|S| m_S(E),
  where m_S(E) = meas{ x : every sector in S is missed } = meas(intersect_{j in S} A_j).
  (m_empty = 1.)

BONFERRONI / BRUN INEQUALITIES on union_j A_j:
  Let U = meas(union A_j) = sum_{k>=1} (-1)^{k+1} Sigma_k,  Sigma_k = sum_{|S|=k} m_S.
  Truncating: partial sums of U alternate as bounds.
    U <= Sigma_1                       (Boole, level 1)
    U >= Sigma_1 - Sigma_2             (level 2)
    U <= Sigma_1 - Sigma_2 + Sigma_3   (level 3)
    ...
  Hence on p0 = 1 - U:
    p0 >= 1 - Sigma_1
    p0 <= 1 - Sigma_1 + Sigma_2
    p0 >= 1 - Sigma_1 + Sigma_2 - Sigma_3
    ...
  So EVEN cumulative levels of the p0-sieve are UPPER bounds on p0; ODD are LOWER.
  To show p0 <= cap we want an EVEN-level upper bound <= cap.

KEY VANISHING (which m_S vanish?):
  The orbit of |E| runners is a set of |E| points on the circle (for generic x).
  It can MISS at most 7 - 1 = 6 sectors only if it lands in 1 sector... but more useful:
  with |E| points you occupy at most |E| sectors, so you MISS at least 7-|E| sectors generically.
  m_S = meas(miss all of S). For m_S>0 we need the orbit to avoid |S| sectors on a
  positive-measure set. Combined with APEX truncation: p0(E)=0 for |E|<7. Here we probe |E|>=7.
"""
from fractions import Fraction as F
from itertools import combinations

def sector_of(p):
    return int((p % 1) * 7)

def breakpoints(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7*e + 1):
            bps.add(F(m, 7*e))
    return sorted(bps)

def cells(E):
    """Return list of (length, set_of_sectors_hit) for each maximal constant cell."""
    E = sorted(set(E))
    bps = breakpoints(E)
    out = []
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0:
            continue
        mid = (x0 + x1) / 2
        hit = frozenset(sector_of(e * mid) for e in E)
        out.append((x1 - x0, hit))
    return out

def p0_exact(E):
    tot = F(0)
    for L, hit in cells(E):
        if len(hit) == 7:
            tot += L
    return tot

def m_S(E, Sset, _cells=None):
    """meas{x : every sector in Sset is MISSED} = sum of cell lengths whose hit-set is disjoint from Sset."""
    Sset = frozenset(Sset)
    cl = _cells if _cells is not None else cells(E)
    tot = F(0)
    for L, hit in cl:
        if hit.isdisjoint(Sset):
            tot += L
    return tot

def sigma_k(E, k, _cells=None):
    cl = _cells if _cells is not None else cells(E)
    tot = F(0)
    for S in combinations(range(7), k):
        tot += m_S(E, S, cl)
    return tot

def sieve_report(E, cap=None, label=""):
    cl = cells(E)
    p0 = sum(L for L, h in cl if len(h) == 7)
    # Sigma_k for k=1..7
    Sig = [None]  # 1-indexed
    for k in range(1, 8):
        Sig.append(sigma_k(E, k, cl))
    # p0 cumulative sieve: p0 = 1 - U,  U = sum (-1)^{k+1} Sig_k
    # p0_partial(L) = 1 - sum_{k=1}^{L} (-1)^{k+1} Sig_k
    print(f"--- {label}  E={E}  (|E|={len(set(E))}) ---")
    print(f"  p0(exact) = {p0} = {float(p0):.6f}")
    if cap is not None:
        print(f"  cap        = {cap} = {float(cap):.6f}   (need p0 <= cap; margin {float(cap-p0):+.6f})")
    U = F(0)
    print(f"  Bonferroni cumulative bounds on p0:")
    for L in range(1, 8):
        U += (-1)**(L+1) * Sig[L]
        p0_partial = 1 - U
        kind = "UPPER" if L % 2 == 0 else "LOWER"
        tight = ""
        if cap is not None and kind == "UPPER":
            tight = "  <= cap? " + ("YES" if p0_partial <= cap else "NO")
        print(f"    level {L} ({kind:5s}): Sig_{L}={float(Sig[L]):.5f}  p0_bound={float(p0_partial):+.6f}{tight}")
    return p0, Sig


if __name__ == "__main__":
    caps = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(4, 7)}

    print("=" * 70)
    print("CONSECUTIVE co-offset clusters E = {0,1,...,k-1}, k = 8,9,10")
    print("=" * 70)
    for k in (8, 9, 10):
        E = list(range(k))
        sieve_report(E, cap=caps.get(k), label=f"consec k={k}")
        print()

    print("=" * 70)
    print("Smaller diagnostic: k=7 (apex boundary), consecutive")
    print("=" * 70)
    sieve_report(list(range(7)), label="consec k=7")


def deep_dive(E, cap, label):
    """Examine the m_S distribution by size, and test a 'gap-aware' truncation.
       Key idea: with |E| runners the orbit occupies <= |E| sectors but for the
       DANGEROUS region we care about x where it occupies exactly 6 (misses 1).
       Test: how much of Sigma_2..Sigma_4 mass sits on 'large gaps' (|S| consecutive)?
    """
    from itertools import combinations
    cl = cells(E)
    print(f"=== DEEP DIVE {label}  E={E} ===")
    # distribution of cells by number of sectors hit
    from collections import defaultdict
    by_hit = defaultdict(F)
    for L, h in cl:
        by_hit[len(h)] += L
    print("  measure of x by #sectors-hit:")
    for nh in sorted(by_hit):
        print(f"    hits {nh} sectors: meas {float(by_hit[nh]):.5f}")
    # so p0 = meas(hits 7). The 'miss exactly 1' band = hits 6.
    # m_S for |S|=1: the single-sector miss measures
    print("  m_{j} (single sector missed), j=0..6:")
    singles = []
    for j in range(7):
        v = m_S(E, {j}, cl)
        singles.append(v)
        print(f"    j={j}: {float(v):.5f}")
    # Are consecutive S the heavy ones? Compare m_S for consecutive vs spread doubletons
    print("  m_S for |S|=2  (consecutive vs antipodal):")
    print(f"    consec {{0,1}}: {float(m_S(E,{0,1},cl)):.5f}")
    print(f"    skip   {{0,2}}: {float(m_S(E,{0,2},cl)):.5f}")
    print(f"    far    {{0,3}}: {float(m_S(E,{0,3},cl)):.5f}")
    print()


if __name__ == "__main__" and False:
    pass
