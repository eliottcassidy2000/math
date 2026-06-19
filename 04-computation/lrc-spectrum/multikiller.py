"""
Multi-killer / lcm-pattern construction at a single k, exact M.
Try to push level a above the single-killer a=2 and the a~3 of primorial-2-3 k.

Constructions tested (all exact M):
 (1) block {1..k-1} with j gaps + j killers, killers = safe central speeds at t*=m/q.
 (2) speeds = j * lcm(1..m) patterns + small generators.
 (3) "stacked killers": several multiples of (k+1) and of 2(k+1)+1 etc.
We brute small parameter grids and report the smallest exact M.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations
sys.path.insert(0, ".")
from fast_M import M_exact_fast
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def Mval(S):
    return M_exact_fast(sorted(set(S)))[0]


def central_safe(m, q, a, hi):
    """safe speeds sorted by how central their residue is (closest to q/2 first)."""
    cand = []
    for v in range(1, hi + 1):
        r = (v * m) % q
        gg = r if r <= q - r else q - r
        if gg >= a:
            cand.append((abs(r - q / 2), v))
    cand.sort()
    return [v for _, v in cand]


def construct1(k, a, ngap, verbose=False):
    """block {1..k-1+? } with ngap interior gaps replaced by ngap central killers."""
    q = a * (k + 1) - 1
    floor = Fraction(1, k + 1)
    best = None
    ms = [m for m in range(1, q) if gcd(m, q) == 1]
    if len(ms) > 20:
        ms = ms[:: max(1, len(ms) // 20)][:20]
    for m in ms:
        # small block: the smallest safe speeds
        small = []
        for v in range(1, 4 * k):
            r = (v * m) % q
            gg = r if r <= q - r else q - r
            if gg >= a:
                small.append(v)
        if len(small) < k:
            continue
        block = small[: k - ngap]
        # killers: central safe speeds not in block
        cen = [v for v in central_safe(m, q, a, 3 * q) if v not in set(block)]
        if len(cen) < ngap:
            continue
        # try the top few central killers in combinations
        cen = cen[:12]
        for kills in combinations(cen, ngap):
            S = block + list(kills)
            if len(set(S)) != k or setgcd(S) != 1:
                continue
            M = Mval(S)
            if M <= floor:
                continue
            if best is None or M < best[0]:
                best = (M, sorted(S), (a, q, m, ngap))
                if verbose:
                    print(f"   a={a} q={q} m={m} ngap={ngap}: M={M} ~{float(M):.7f} "
                          f"a_act={float(level_a(M,k)):.3f}", flush=True)
    return best


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 31
    floor = Fraction(1, k + 1)
    overall = None
    for a in range(2, int(k**0.5) + 4):
        for ngap in range(1, 4):
            b = construct1(k, a, ngap, verbose=True)
            if b and (overall is None or b[0] < overall[0]):
                overall = b
    if overall:
        M, S, info = overall
        print(f"\nBEST k={k}: M={M} (~{float(M):.10f}) g*k^2={float((M-floor)*k*k):.6f} "
              f"a={float(level_a(M,k)):.4f} info={info}")
        print(f"  S={S}", flush=True)
