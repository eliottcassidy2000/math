"""
ROUTE A, multi-anchor TAIL: does the few-anchor tail P(max_{a in A} gap_a > 1/7)
discharge the WHOLE (A') ledger?  boxeph-2026-07-07-S1, HYP-4760.

Key: maxgap >= max_{a in A} gap_a pointwise, so
    mu_{1/7}(E) = P(maxgap>1/7) >= P(max_{a in A} gap_a > 1/7)  =: PA_{|A|}(E).
Klein-S153's anchor floor was on the MEAN (knife-edge ~0.147); the TAIL is what the
skeleton actually needs (mu_{1/7} = a tail), and it is comfortable (Route A single
anchor already ~0.28).  Question: what anchor count discharges every leg
inf_k PA >= T_k (union-bound thresholds)?  1-anchor cleared k=9..13, missed k=8.

Test A={0}, {0,1/2}, {0,1/3,2/3} per k.  Each PA is STILL a lower bound on mu, and
uses only |A| gaps (not the full max) -- a bounded, tractable object.
"""
import random
def gap_at(pts_sorted, a):
    # gap of the circular config (sorted floats in [0,1)) containing anchor a
    rn=None
    for p in pts_sorted:
        if p>a: rn=p; break
    if rn is None: rn=pts_sorted[0]+1.0
    ln=None
    for p in pts_sorted:
        if p<=a: ln=p
    if ln is None: ln=pts_sorted[-1]-1.0
    return rn-ln

def PA_fast(E, anchors, G=12000, thr=1/7):
    c=0
    for aa in range(G):
        x=(aa+0.5)/G
        pts=sorted((e*x)%1.0 for e in E)
        best=0.0
        for a in anchors:
            g=gap_at(pts,a)
            if g>best: best=g
            if best>thr: break
        if best>thr: c+=1
    return c/G

Tk = {8:0.6185, 9:0.5057, 10:0.3956, 11:0.2747, 12:0.1429, 13:0.0565}
anchor_sets = {
    "{0}":       [0.0],
    "{0,1/2}":   [0.0, 0.5],
    "{0,1/3,2/3}":[0.0, 1/3, 2/3],
}
random.seed(99)
# adversarial families per k (spread APs dominate the minimizer; add random)
def cands_for(k):
    # the P(gap0>1/7) minimizer is always a spread AP {a+dj}; focus there + a
    # modest random probe. (leaner than S1 first pass so it finishes fast.)
    cs=[list(range(1,k+1))]
    for d in range(2,13):
        for a0 in range(0,d):
            E=sorted(set(e for e in (a0+d*j for j in range(k)) if e!=0))
            if len(E)==k: cs.append(E)
    for _ in range(80):
        cs.append(sorted(random.sample(range(1,70),k)))
    return cs

print("inf_k P(max_{a in A} gap_a > 1/7) vs union-bound thresholds T_k")
print(f"{'k':>3} | {'T_k':>7} | " + " | ".join(f"{lbl:>12}" for lbl in anchor_sets))
print("-"*70)
discharged={lbl:[] for lbl in anchor_sets}
for k in range(8,14):
    cs=cands_for(k)
    row=f"{k:>3} | {Tk[k]:>7.4f} |"
    for lbl,A in anchor_sets.items():
        inf=min(PA_fast(E,A) for E in cs)
        ok = inf>=Tk[k]
        if ok: discharged[lbl].append(k)
        row += f" {inf:>7.4f}{'  Y' if ok else '  .'} |"
    print(row)

print("\nDISCHARGED legs by anchor set:")
for lbl in anchor_sets:
    d=discharged[lbl]
    print(f"  A={lbl:12s}: k = {d}  {'<== ALL k=8..13' if set(d)==set(range(8,14)) else ''}")
print("""
Interpretation: if a fixed 2- or 3-anchor TAIL discharges every leg, the whole
(A') load-bearing lemma reduces from 'full maxgap AP-minimality' to a
'few-anchor avoidance tail' -- a bounded functional (|A| gaps), and each anchor's
gap is a translated-Steinhaus three-distance object (opus-S134 Farey roof extends
to a shifted origin).  That is a materially simpler proof target than the full max.
""")
