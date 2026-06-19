"""
Promote r elements of the AP {1..k} to large 'killer' speeds.
Take base {1..k}, choose r positions to REPLACE by killers = chosen multiples that
are safe at the target t* = m/q, q = a(k+1)-1, level a/q. Search over which positions
to promote and which killer residues, evaluate EXACT M, maximize level a.

This is the structure of the observed minimizers: k=7 {1,2,3,4,5,7,18} = base {1..7}
with 6 promoted-out (gap at 6) ... actually it's {1..5,7} + killer 18 (dropped 6, kept 7,
added 18). So: a (k-1)-subset of small integers (with a few gaps) + 1 killer.

We search: choose t*=m/q (q=a(k+1)-1). Safe residues mod q. Build S greedily/exhaustively
from safe speeds but allow the search to PICK which safe speeds (small subset search via
itertools over a limited safe pool), to find a set whose true M == a/q (level a realized).
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


def safe_pool(m, q, a, hi):
    out = []
    for v in range(1, hi + 1):
        r = (v * m) % q
        gg = r if r <= q - r else q - r
        if gg >= a:
            out.append(v)
    return out


def try_level(k, a, hi_mult=2, pool_cap=None, killer_count=2):
    """For target level a (q=a(k+1)-1), find a k-subset of the safe pool whose exact
    M equals a/q. Strategy: take the k-killer_count SMALLEST safe speeds as the 'block',
    then try all combinations of killer_count additional safe speeds from a window near q.
    Return best (smallest M) realized."""
    q = a * (k + 1) - 1
    target = Fraction(a, q)
    floor = Fraction(1, k + 1)
    best = None
    for m in range(1, q // 2 + 1):
        if gcd(m, q) != 1:
            continue
        pool = safe_pool(m, q, a, hi_mult * q)
        if len(pool) < k:
            continue
        block = pool[:k - killer_count]
        # killer candidates: safe speeds beyond the block, prefer those near multiples of q
        killer_pool = [v for v in pool if v not in block]
        if len(killer_pool) < killer_count:
            continue
        # limit killer_pool size for combinatorial feasibility
        killer_pool = killer_pool[:30]
        for kill in combinations(killer_pool, killer_count):
            S = sorted(set(block) | set(kill))
            if len(S) != k:
                continue
            if min(S) <= 0 or setgcd(S) != 1:
                continue
            M, t = M_exact_fast(S)
            if M <= floor:
                continue
            if best is None or M < best['M']:
                g = M - floor
                best = dict(S=S, M=M, t=t, g=g, gk2=g * k * k,
                            a=level_a(M, k), m=m, q=q, target=target,
                            hit=(M == target))
    return best


if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [12, 16, 20]
    for k in ks:
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} ====", flush=True)
        overall = None
        for a in range(2, int(k**0.5) + 4):
            b = try_level(k, a, killer_count=2)
            if b and (overall is None or b['gk2'] < overall['gk2']):
                overall = b
            if b:
                print(f"   target a={a} q={b['q']}: realized M={b['M']} (~{float(b['M']):.8f}) "
                      f"a_act={float(b['a']):.3f} gk2={float(b['gk2']):.5f} hit_target={b['hit']} S={b['S']}", flush=True)
        if overall:
            print(f"  >>> BEST k={k}: a={float(overall['a']):.4f} gk2={float(overall['gk2']):.5f} M={overall['M']}", flush=True)
