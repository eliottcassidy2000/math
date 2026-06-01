from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce

def frac_dist(x):
    """||x|| = distance to nearest integer, for Fraction x."""
    f = x - int(x)        # fractional part in (-1,1)
    if f < 0: f += 1
    return min(f, 1 - f)

def D(t, speeds):
    """min_i ||v_i t||"""
    return min(frac_dist(v * t) for v in speeds)

def witness_times(speeds):
    """Distinct candidate optimal times t = k/(v_i+v_j) in (0,1), with provenance."""
    times = {}   # t -> set of pairs (i,j) that generated it
    for i, j in combinations(range(len(speeds)), 2):
        s = speeds[i] + speeds[j]
        for k in range(1, s):
            t = Fraction(k, s)
            times.setdefault(t, set()).add((i, j))
    return times

def M_witness(speeds):
    """M(v) using witness-time encoding. Returns (M, best_t, generating_pairs)."""
    times = witness_times(speeds)
    bestM = Fraction(0)
    best_t = None
    best_pairs = None
    for t, pairs in times.items():
        d = D(t, speeds)
        if d > bestM:
            bestM = d
            best_t = t
            best_pairs = pairs
    return bestM, best_t, best_pairs

def M_brute(speeds, denom_cap=None):
    """Brute force M over all rationals t = a/q, q <= Q. Q = lcm-ish bound.
    For verification only. Uses denominator = product-ish; we use all q up to cap."""
    if denom_cap is None:
        denom_cap = 2 * sum(speeds)  # plenty: optimum has denom v_i+v_j <= 2max < this
    bestM = Fraction(0)
    best_t = None
    for q in range(1, denom_cap + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            t = Fraction(a, q)
            d = D(t, speeds)
            if d > bestM:
                bestM = d
                best_t = t
    return bestM, best_t

def n_of(speeds):
    return len(speeds) + 1

def binding_pair(speeds, best_t):
    """Among runners, which two are equidistant on opposite sides at best_t and are the min?
    Return the pair (i,j) achieving the min distance with opposite-side structure."""
    n = len(speeds)
    dmin = D(best_t, speeds)
    # runners achieving the min
    achievers = [i for i in range(n) if frac_dist(speeds[i]*best_t) == dmin]
    # find a pair on opposite sides
    sides = {}
    for i in achievers:
        f = speeds[i]*best_t
        ff = f - int(f)
        if ff < 0: ff += 1
        # side: <0.5 means just above integer (positive side), >0.5 means below (negative)
        s = '+' if ff <= Fraction(1,2) else '-'
        sides[i] = s
    pos = [i for i in achievers if sides[i]=='+']
    neg = [i for i in achievers if sides[i]=='-']
    if pos and neg:
        return (min(pos), min(neg)) if min(pos)<min(neg) else (min(neg),min(pos))
    # fallback: just the two smallest achievers
    if len(achievers)>=2:
        return tuple(sorted(achievers[:2]))
    return (achievers[0], achievers[0])

if __name__ == "__main__":
    # sanity checks against known small LRC results
    tests = [
        [1,2,3],      # m=3,n=4
        [1,2,3,4],    # m=4,n=5
        [1,3,4,7],    # known tight-ish set
    ]
    for sp in tests:
        Mw, tw, pw = M_witness(sp)
        Mb, tb = M_brute(sp)
        n = n_of(sp)
        print(sp, "n=",n, "M_wit=",Mw, "M*n=",float(Mw*n), "agree=", Mw==Mb, "t=",tw)
