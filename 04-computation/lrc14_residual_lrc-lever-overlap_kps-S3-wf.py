#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_kps-S3-wf

ANGLE: "lrc-lever-overlap".
Goal: close the ONLY open case S3 of LRC(14):
  S3 = covering 13-sets with k>=2 large runners and Vmax >= 13*Vmin.

Structure of an S3 set:  S = P  u  L, where
  - L = a TIGHT large CLUSTER, all speeds in a window [V0, V0+s], internal spread s small.
  - P = the rest (small runners and possibly other groups).
We want to prove the criterion C(S):  exists v in S with W(S\{v}) > 1/(7v).
Empirically v = max(L) works.  So set v=Vmax, A = S\{v} = P u (L\{v}).

KEY IDEA (lever-overlap):
  G_A = G_P  cap  G_{L'}   (L' = L\{v}).
  - G_P: safe set of the small part. WIDE arcs (it's a small lonely-runner system).
  - G_{L'}: safe set of the cluster minus top. The cluster is "almost periodic":
      each u in L' is u = V0 + d_u (d_u = u - V0, the OFFSET, a BOUNDED integer 0..s).
      As V0 -> infty the cluster collapses; G_{L'} has a fine almost-periodic structure
      whose PATTERN depends only on the offset multiset {d_u}, not on V0.

This script:
  (0) tools (exact Fraction).
  (1) S3 generator + decomposition into (P, L).
  (2) MEASURE meas(G_L) and widest arc of G_L exactly; relate to offsets {d_u}.
  (3) The OVERLAP experiment: measure widest arc of G_P cap G_{L'} and compare to 1/(7 Vmax).
  (4) The ALMOST-PERIODIC structure: show G_{L'} restricted to a unit V0-cell depends only
      on offsets; quantify the "cluster band" width and how it sweeps tau in [0,1).
  (5) Honest residual: where (if anywhere) the closed-form overlap bound fails.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()

C = F(1, 14)   # the LRC(14) level

# ---------- exact tools (verbatim from prompt) ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r

def safe_components(A, h=C):
    """Return list of maximal safe arcs (lo,hi) of G_A = {tau: ||a tau||>=h all a in A}."""
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe

def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append((sc[0][1])+(1-sc[-1][0]))
    return max(ws)

def measure_safe(A):
    """Lebesgue measure of G_A = sum of safe-arc widths."""
    return sum(b-a for a,b in safe_components(A))

def cand(S):
    S=sorted(set(S)); Cs=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cs.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cs.add(F(k,d)); k+=1
    Cs.add(F(1,2)); return Cs

def Mval(S):
    b=F(0); at=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; at=t
    return b, at

def is_covering(S):
    return all(any(v%q==0 for v in S) for q in range(2,15))

def gcd_all(S): return reduce(gcd, S)

# ---------- interval intersection (exact) ----------
def intersect_arclists(L1, L2):
    """Intersect two lists of disjoint sorted (lo,hi) arcs on [0,1)."""
    out=[]; i=j=0
    while i<len(L1) and j<len(L2):
        a,b=L1[i]; c,d=L2[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: out.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    return out

def widest_arc(arclist):
    """Widest arc with circular wrap on [0,1)."""
    if not arclist: return F(0)
    ws=[b-a for a,b in arclist]
    if arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        ws.append(arclist[0][1] + (1-arclist[-1][0]))
    return max(ws)

# ---------- S3 set generator ----------
def make_S3_set(rng, max_attempts=20000):
    """Generate covering primitive 13-sets in case S3:
       k>=2 large (>13) AND Vmax >= 13*Vmin.
       Bias toward a tight large cluster + small part."""
    for _ in range(max_attempts):
        # small part: pick a subset of 1..13 covering as much as possible
        small = set(rng.sample(range(1,14), rng.randint(4,9)))
        # cluster: pick a base V0 and tight offsets
        V0 = rng.randint(20, 400)
        s = rng.randint(8, 45)  # internal spread of cluster
        clsize = rng.randint(2, 13 - len(small))
        cluster = set()
        tries=0
        while len(cluster) < clsize and tries < 200:
            cluster.add(V0 + rng.randint(0, s)); tries+=1
        S = small | cluster
        # top up to 13 if short, with more large runners
        while len(S) < 13:
            S.add(rng.randint(V0, V0+s) if rng.random()<0.5 else rng.randint(1,13))
        if len(S) != 13: continue
        S = set(list(S))
        if len(S)!=13: continue
        if gcd_all(S)!=1: continue
        if not is_covering(S): continue
        Vmin=min(S); Vmax=max(S)
        k = sum(1 for v in S if v>13)
        if k>=2 and Vmax >= 13*Vmin:
            return S
    return None

def decompose(S, cluster_gap=14):
    """Split S into small part P and large cluster(s) L. Return (P, L) where
       L is the TOP cluster (around Vmax). We define L = the maximal run of large
       speeds within a window of width < 13*Vmin-ish near the top."""
    S=sorted(S)
    Vmin=S[0]
    larges=[v for v in S if v>13]
    if not larges: return set(S), set()
    # cluster = the topmost contiguous group within spread bound
    larges_sorted=sorted(larges)
    # group consecutive larges where gap < some threshold relative to value
    top=larges_sorted[-1]
    L=[top]
    i=len(larges_sorted)-2
    while i>=0 and top - larges_sorted[i] < min(larges_sorted)*0 + 14*3:
        # within a tight band of the top
        if top - larges_sorted[i] <= 50:
            L.append(larges_sorted[i]); i-=1
        else: break
    L=set(L)
    P=set(S)-L
    return P, L

# =========================================================
# EXPERIMENT 1: verify criterion C(S) holds via v=Vmax on S3 sets
#               and record the OVERLAP structure.
# =========================================================
def run_experiment(n_sets=300, seed=20260617):
    rng=random.Random(seed)
    flush("="*70)
    flush("EXPERIMENT 1: criterion C(S) via v=Vmax, overlap structure of G_P cap G_{L'}")
    flush("="*70)
    fails_C=0
    margins=[]
    overlap_data=[]
    got=0
    attempts=0
    while got<n_sets and attempts<n_sets*40:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S)
        A=set(S)-{Vmax}
        W=Wwidth(A)
        thresh=F(1,7*Vmax)
        margin = W/thresh   # >1 means C fires via Vmax
        margins.append(margin)
        if W <= thresh:
            fails_C+=1
            flush(f"  C-FAIL via Vmax: S={sorted(S)} W={W} thresh={thresh} margin={float(margin):.3f}")
            # check if ANOTHER v works
            anyv=False
            for v in S:
                if Wwidth(set(S)-{v}) > F(1,7*v): anyv=True; break
            if not anyv:
                flush(f"    *** TRUE C-FAILURE (no v works): {sorted(S)} ***")
        # overlap decomposition
        P,L=decompose(S)
        if L and Vmax in L:
            Lp = L-{Vmax}
            GP = safe_components(P) if P else [(F(0),F(1))]
            GLp = safe_components(Lp) if Lp else [(F(0),F(1))]
            inter = intersect_arclists(GP, GLp)
            wover = widest_arc(inter)
            overlap_data.append((len(P),len(L),Vmax,wover,thresh, float(wover/thresh) if thresh>0 else 0))
    flush(f"\n  S3 sets tested: {got}")
    flush(f"  C(S) via Vmax FAILURES: {fails_C}")
    if margins:
        mn=min(margins); mx=max(margins)
        flush(f"  margin = W(S\\Vmax) / (1/(7 Vmax)):  min={float(mn):.4f}  max={float(mx):.4f}")
        flush(f"           exact min margin = {mn}")
    # overlap stats
    if overlap_data:
        ratios=[d[5] for d in overlap_data]
        flush(f"\n  OVERLAP widest-arc(G_P cap G_L') / (1/(7Vmax)) : min={min(ratios):.4f} max={max(ratios):.4f}")
        worst=min(overlap_data, key=lambda d:d[5])
        flush(f"  worst overlap config: |P|={worst[0]} |L|={worst[1]} Vmax={worst[2]} wover={worst[3]} ratio={worst[5]:.4f}")
    return margins

if __name__=="__main__":
    run_experiment(n_sets=300)
