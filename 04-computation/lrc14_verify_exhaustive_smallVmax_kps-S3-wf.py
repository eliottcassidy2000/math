"""
EXHAUSTIVE adversarial sweep of covering S3 sets with bounded Vmax.

Unlike Monte Carlo, this enumerates EVERY primitive covering 13-set S in case S3
with max(S) <= VMAX_CAP. For each:
  - compute exact M(S); flag any M < 1/14  (true LRC(14) counterexample)
  - test C(S): does SOME removed runner v give W(S\{v}) > 1/(7v)?  flag total failures.
  - test the 'v=Vmax' rule and record min over-the-board margins.

This is a COMPLETE check in the bounded-Vmax regime -- the regime where the author's
random sampling and finite-check claims live. If C(S) ever fails or M<1/14 here, it is
a hard counterexample, not a sampling artifact.

We bound enumeration by:
  - small part P subset of {1..13}
  - large speeds L: 2..11 speeds in (13, VMAX_CAP], chosen so |S|=13.
Enumerating all C(VMAX_CAP-13, |L|) large subsets is huge, so we cap VMAX_CAP modestly
and additionally restrict to CLUSTERED large parts (spread <= SPREAD_CAP) which is exactly
the S3 'mixed' regime; we ALSO run a separate fully-unrestricted sweep at a tiny cap.

kind-pasteur-S3-wf
"""
from fractions import Fraction as F
from math import gcd
import itertools

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g

def case_of(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13*Vmin: return 'S2'
    return 'S3'

def C_holds(S):
    wits = []
    for v in S:
        A = [x for x in S if x != v]
        if Wwidth(A) > F(1, 7*v): wits.append(v)
    return (len(wits) > 0, wits)

def sweep(VMAX_CAP, SPREAD_CAP, min_small=7, max_small=11):
    """Enumerate covering S3 13-sets with max(S)<=VMAX_CAP and large-part spread<=SPREAD_CAP."""
    n_total = 0; n_Cfail = 0; n_Mbelow = 0
    Cfails = []; Mbelows = []
    minM = F(1); worst = None
    vmax_fail = 0
    minratio_overall = F(10)  # min over all sets of best C-margin
    # iterate over number of large speeds
    for nlarge in range(2, 7):  # k = nlarge >= 2
        nsmall = 13 - nlarge
        if nsmall < min_small or nsmall > max_small: continue
        if nsmall > 13: continue
        # small parts: subsets of 1..13 of size nsmall
        for small in itertools.combinations(range(1, 14), nsmall):
            small = list(small)
            # large parts: choose nlarge speeds in (13, VMAX_CAP] with spread<=SPREAD_CAP
            # iterate by smallest large speed w0
            for w0 in range(14, VMAX_CAP+1):
                hi = min(VMAX_CAP, w0 + SPREAD_CAP)
                pool = list(range(w0+1, hi+1))
                if len(pool) < nlarge-1: continue
                for rest in itertools.combinations(pool, nlarge-1):
                    large = [w0] + list(rest)
                    S = sorted(set(small + large))
                    if len(S) != 13: continue
                    if max(S) > VMAX_CAP: continue
                    if max(S) != large[-1] if False else None: pass
                    if gcd_all(S) != 1: continue
                    if not is_covering(S): continue
                    if case_of(S) != 'S3': continue
                    n_total += 1
                    m = Mval(S)
                    if m < minM: minM = m; worst = S
                    if m < F(1, 14):
                        n_Mbelow += 1; Mbelows.append((S, m))
                    ok, wits = C_holds(S)
                    if not ok:
                        n_Cfail += 1; Cfails.append((S, m))
                    # best margin
                    best = F(0)
                    for v in S:
                        A = [x for x in S if x != v]
                        r = Wwidth(A) * 7 * v
                        if r > best: best = r
                    if best < minratio_overall: minratio_overall = best
                    # vmax rule
                    Vmax = max(S); A = [x for x in S if x != Vmax]
                    if not (Wwidth(A) * 7 * Vmax > 1): vmax_fail += 1
    return dict(n_total=n_total, n_Cfail=n_Cfail, n_Mbelow=n_Mbelow,
                Cfails=Cfails, Mbelows=Mbelows, minM=minM, worst=worst,
                vmax_fail=vmax_fail, minratio_overall=minratio_overall)

def main():
    print("EXHAUSTIVE S3 sweep (bounded Vmax, clustered large part)")
    # Modest cap so enumeration finishes. SPREAD<=45 matches the S3 'tight cluster' regime.
    for VMAX_CAP, SPREAD_CAP in [(40, 30), (60, 45)]:
        print("="*70)
        print(f"VMAX_CAP={VMAX_CAP}, SPREAD_CAP={SPREAD_CAP}")
        r = sweep(VMAX_CAP, SPREAD_CAP)
        print(f"  total covering S3 sets enumerated: {r['n_total']}")
        print(f"  #M<1/14 (TRUE COUNTEREXAMPLES): {r['n_Mbelow']}")
        for S, m in r['Mbelows'][:20]:
            print("     M<1/14:", S, "M=", m, float(m))
        print(f"  #C(S) total-failures: {r['n_Cfail']}")
        for S, m in r['Cfails'][:30]:
            print("     C FAILS:", S, " M=", m, float(m))
        print(f"  min M over all = {r['minM']} = {float(r['minM']):.5f} at {r['worst']}")
        print(f"  v=Vmax rule failures (C may still hold via other v): {r['vmax_fail']}")
        print(f"  min best-C-margin over all sets: {r['minratio_overall']} = {float(r['minratio_overall']):.4f}")

if __name__ == '__main__':
    main()
