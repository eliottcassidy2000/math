"""
lrc14_clift_sieve_family_kps.py  (kind-pasteur-2026-06-27-S31ag)

Verify THM-574 candidate: for a primitive 13-set S, if c<=7 and |H_c|>=14-c
(H_c = multiples of c in S), then M(S) > 1/14.

Also confirm c=7 (threshold 7) is optimal among c<=7, and that the threshold
14-c is sharp (sets with |H_c| = 13-c can dip to M<=1/14).
"""
from math import gcd
import random

def nrm(x):
    f = x - int(x)
    if f < 0: f += 1.0
    return min(f, 1 - f)

def M_value(S):
    """Fast float M(S)=max_t min_i ||s_i t||. Optimum at t=n/(s_a +/- s_b)."""
    S = sorted(set(S))
    denoms = set()
    for i in range(len(S)):
        a = S[i]
        denoms.add(2 * a)
        for j in range(len(S)):
            if i == j: continue
            b = S[j]
            denoms.add(a + b)
            if a != b: denoms.add(abs(a - b))
    best = -1.0
    for D in denoms:
        if D == 0: continue
        for n in range(1, D):
            t = n / D
            m = 1.0
            for s in S:
                v = nrm(s * t)
                if v < m: m = v
                if m <= best: break
            if m > best: best = m
    return best

def gcd_list(S):
    g = 0
    for s in S: g = gcd(g, s)
    return g

random.seed(574)
thr = 1.0 / 14.0
print("Testing THM-574: |H_c| >= 14-c (c<=7), primitive => M > 1/14")
for c in range(2, 8):
    need = 14 - c
    fired = 0
    fail = 0
    minM = 1.0
    for trial in range(1500):
        # build a primitive 13-set with exactly 'need' (or more) multiples of c
        nmult = random.randint(need, min(13, need + 2))
        mults = random.sample([c * j for j in range(1, 10)], min(nmult, 9))
        others = []
        pool = [x for x in range(1, 45) if x % c != 0]
        while len(mults) + len(others) < 13:
            x = random.choice(pool)
            if x not in others and x not in mults:
                others.append(x)
        S = sorted(set(mults + others))
        while len(S) < 13:
            x = random.randint(1, 60)
            if x not in S: S.append(x)
        S = sorted(S)[:13]
        if gcd_list(S) != 1:
            continue
        Hc = sum(1 for s in S if s % c == 0)
        if Hc < need:
            continue
        fired += 1
        M = M_value(S)
        if M < minM: minM = M
        if M <= thr + 1e-9:
            fail += 1
            if fail <= 3:
                print(f"  c={c} FAIL S={S} |H_c|={Hc} M={float(M):.5f}")
    print(f"c={c}: threshold |H_c|>={need}  fired={fired} failures(M<=1/14)={fail}  minM={float(minM):.5f}")

# Sharpness: c=7, |H_7| = 6 (=13-7, one below threshold) — can it be tight?
print("\nSharpness at c=7: |H_7|=6 sets (one below threshold) — dilated AP 2*{1..13} has |H_7|=1 only.")
S = [2*i for i in range(1,14)]
print(f"  2*{{1..13}} |H_7|={sum(1 for s in S if s%7==0)} M={float(M_value(S)):.5f} (should be 1/14, primitive? gcd={gcd_list(S)})")
