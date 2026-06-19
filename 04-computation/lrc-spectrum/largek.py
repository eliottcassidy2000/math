"""
Large-k targeted search. Same minimizer shape (near-prefix with <=G gaps + r killers),
but with killers drawn from a FOCUSED set: speeds that are safe at a target t*=m/q,
q in {a(k+1)-1 : a small}, so we directly probe the level-a family.

For each target level a (q=a(k+1)-1) and coprime m:
  - block = {1,...,k-r} possibly with up to G interior gaps (chosen to be the residues
    that are NOT safe, i.e. forcing the optimum to t*),
  - killers = r safe speeds in a window [q-W, q+W] (near the modulus, the natural killers).
Compute EXACT M, keep smallest.
Also always include the plain single-killer baseline.
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


def is_safe(v, m, q, a):
    r = (v * m) % q
    gg = r if r <= q - r else q - r
    return gg >= a


def baseline_single_killer(k):
    floor = Fraction(1, k + 1)
    base = list(range(1, k))
    best = None
    for X in range(k, 5 * k):
        S = base + [X]
        if setgcd(S) != 1:
            continue
        M = Mval(S)
        if M <= floor:
            continue
        if best is None or M < best[0]:
            best = (M, sorted(S))
    return best


def level_a_probe(k, a, G=1, r=1, W=None):
    q = a * (k + 1) - 1
    if W is None:
        W = 2 * k
    floor = Fraction(1, k + 1)
    best = None
    ms = [m for m in range(1, q) if gcd(m, q) == 1]
    if len(ms) > 16:
        step = max(1, len(ms) // 16)
        ms = ms[::step][:16]
    for m in ms:
        # killers: safe speeds in [k+something, q+W]
        killer_cands = [v for v in range(k, q + W + 1) if is_safe(v, m, q, a)]
        # block region: small speeds; identify which small speeds are safe
        small_safe = [v for v in range(1, k + 8) if is_safe(v, m, q, a)]
        if len(small_safe) < k - r:
            continue
        # block = first (k-r) safe small speeds; also try dropping up to G of them (gaps)
        base_block = small_safe[:k - r + G]
        # choose which (k-r) to keep among first (k-r+G): i.e. drop G
        from itertools import combinations as comb
        for drop in comb(range(len(base_block)), G):
            block = [base_block[i] for i in range(len(base_block)) if i not in drop]
            if len(block) != k - r:
                continue
            # killers: choose r from killer_cands restricted to a window to keep small
            kc = [v for v in killer_cands if v > max(block)]
            if len(kc) < r:
                continue
            kc = kc[:25]
            for killers in comb(kc, r):
                S = block + list(killers)
                if len(set(S)) != k or setgcd(S) != 1:
                    continue
                M = Mval(S)
                if M <= floor:
                    continue
                if best is None or M < best[0]:
                    best = (M, sorted(S), (a, q, m))
    return best


if __name__ == "__main__":
    import sympy
    ks = [int(x) for x in sys.argv[1:]] or [30, 31, 42, 43, 60, 61]
    for k in ks:
        floor = Fraction(1, k + 1)
        om = len(sympy.factorint(k - 1))
        bestM, bestS, why = None, None, None
        b0 = baseline_single_killer(k)
        if b0:
            bestM, bestS, why = b0[0], b0[1], "single-killer(a=2)"
        amax = int(k**0.5) + 3
        for a in range(2, amax + 1):
            for r in (1, 2):
                for G in (0, 1, 2):
                    pr = level_a_probe(k, a, G=G, r=r)
                    if pr and (bestM is None or pr[0] < bestM):
                        bestM, bestS, why = pr[0], pr[1], f"level-a a={a} q={pr[2][1]} m={pr[2][2]} r={r} G={G}"
        g = bestM - floor
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} omega(k-1)={om} ====")
        print(f"  BEST M={bestM} (~{float(bestM):.10f}) g*k^2={float(g*k*k):.6f} "
              f"a={float(level_a(bestM,k)):.4f} via {why}")
        print(f"  S={bestS}", flush=True)
