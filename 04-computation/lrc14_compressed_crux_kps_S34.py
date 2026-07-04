"""kps-S34: RESOLVE opus/mac-mini disagreement -- does the census denominator q grow for
COMPRESSED (ratio<=2, no far runner) covering band-blockers, or stay bounded (17-19)?
Strong compressed adversary: 13 speeds in [N,2N] divisibility-blocking {2..Q}, Q maximal."""
from math import gcd
def lcm(a,b): return a*b//gcd(a,b)

def is_lonely_mod(speeds, a, q):
    """t=a/q lonely (margin 1/14): every runner keeps dist(v a/q) >= 1/14."""
    for v in speeds:
        r = (v % q) * a % q
        m = min(r, q-r)
        if 14*m < q:              # min(r,q-r)/q < 1/14  => danger
            return False
    return True

def q_min(speeds, Qmax=400):
    """smallest denominator q with a lonely a/q, or None if > Qmax."""
    for q in range(2, Qmax+1):
        for a in range(1, q):
            if is_lonely_mod(speeds, a, q):
                return q, a
    return None, None

def first_free_modulus(speeds, Qmax=400):
    """first q dividing no runner (q | v => residue-0 doom for all a)."""
    for q in range(2, Qmax+1):
        if all(v % q != 0 for v in speeds):
            return q
    return None

def build_compressed_blocker(N, Qtry):
    """13 numbers in [N,2N], each = lcm(bin) * cofactor, divisibility-blocking {2..Qtry}.
    Greedy bin-pack {2..Qtry} into bins with lcm<=2N; need <=13 bins."""
    bins = []                              # each: [lcm_so_far, set_of_q]
    for q in range(2, Qtry+1):
        placed = False
        for b in bins:
            l = lcm(b[0], q)
            if l <= 2*N:
                b[0] = l; b[1].add(q); placed = True; break
        if not placed:
            if q <= 2*N:
                bins.append([q, {q}])
            else:
                return None
    if len(bins) > 13:
        return None
    runners = []
    for b in bins:
        l = b[0]
        c = max(1, (N + l - 1)//l)         # smallest multiple of l that is >= N
        r = l*c
        if r > 2*N:                        # l itself already in (N,2N]?
            if N <= l <= 2*N: r = l
            else: return None
        runners.append(r)
    # pad to 13 with fillers in [N,2N] that are multiples of an already-blocked small modulus
    filler_base = bins[0][0] if bins else 2
    k = 1
    while len(runners) < 13:
        r = filler_base * max(1, (N + filler_base*k - 1)//(filler_base*k))  # rough
        r = filler_base * ((N)//filler_base + k)
        if r < N: r = N + k
        if r > 2*N: r = N + k
        runners.append(r); k += 1
    return runners[:13]

def max_blockable_Q(N):
    """largest Q such that {2..Q} is divisibility-blockable by 13 compressed runners."""
    best = None
    Q = 14
    while True:
        fam = build_compressed_blocker(N, Q)
        if fam is None:
            break
        best = (Q, fam)
        Q += 1
        if Q > 300: break
    return best

print(f"{'N':>14} {'maxQ block':>11} {'q_min':>7} {'first_free':>11} {'ratio':>7} {'log2(N)':>8}")
import math
for e in range(3, 16):
    N = 10**e
    res = max_blockable_Q(N)
    if res is None:
        print(f"{N:>14,}  (no block)"); continue
    Q, fam = res
    fam = sorted(set(fam))
    ratio = max(fam)/min(fam)
    qm, am = q_min(fam, 400)
    ff = first_free_modulus(fam, 400)
    print(f"{N:>14,} {Q:>11} {str(qm):>7} {str(ff):>11} {ratio:>7.2f} {math.log2(N):>8.1f}")

print()
print("KEY: if q_min GROWS with log(N) => mac-mini right (no bounded census, crux=CRT-capacity).")
print("     if q_min BOUNDED (17-19)  => opus right (13 runners bound the coverable band).")
print("Note: ratio<=2 confirms COMPRESSED (no far runner to peel).")
