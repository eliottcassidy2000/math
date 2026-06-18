#!/usr/bin/env python3
"""
lrc14_finish_covering-tightness-7adic-lift2_kps-S4-wf.py   (kps-S4-wf, part 2)

THE GENUINELY-TIGHT REGIME. Part 1 confirmed: on a broad random sample of covering
primitive S3 sets the M-floor is ~30/301 (well above 1/14) and the optimum tau* sits
near some k/7 with the mult-of-7 runner frequently binding. But the RANDOM sample
misses the genuinely tight family (M near 1/14). Here we attack the tight family
DIRECTLY and test the 7-adic-lattice gap-side certificate on it.

KNOWN TIGHT SETS (from canon):
  S* = {1,2,3,5,7,8,9,10,11,12,13,38,42}  -- covering, S3, k=2, M = 2/23 (the realized floor)
  THM-524 family A={1..11,13}+{84k}: M = 7k/(84k+5)  (resonant dip)

THE WINDOW-COLLAPSE LEMMA (the rigorous lever, from lrc_n14_seven_impossibility_s552):
  Near tau = k/7 (gcd(k,7)=1), every non-multiple-of-7 runner v has ||v(k/7)|| in {1/7,2/7,3/7}
  >= 1/7, so it stays SAFE (>=1/14) for ALL |s| <= 1/(14 V*), V* = max non-mult-7 runner.
  Inside that window LRC(14) reduces EXACTLY to the multiples-of-7 sub-system:
        find s with ||7 w_i s|| >= 1/14 for every mult-of-7 7 w_i, AND |s| <= 1/(14 V*).
  Sub: with u = 7s the constraint is ||w_i u|| >= 1/14, |u| <= 1/(2 V*). A tiny LRC on {w_i}.

DELIVERABLES (all EXACT Fraction):
 (A) Verify S* and THM-524 values; show their tau* 7-adic position.
 (B) THE WINDOW-COLLAPSE CERTIFICATE. For the tight sets, scan the window near each k/7
     for an EXACT s that clears the multiples-of-7 sub-system; if found with |s| <= 1/(14 V*),
     it is a GLOBAL WITNESS => M(S) >= 1/14, PROVEN constructively. Report which sets close.
 (C) THE SUB-CASE THEOREM. Characterize exactly when the window-collapse closes: it closes iff
     the mult-of-7 sub-LRC {w_i} (gap 1/14) has a witness in |u| <= 1/(2V*). When r = #mult-of-7
     is small (<= a few) and the w_i are not adversarial, this ALWAYS holds. Prove the r=1 and
     r=2 cases unconditionally and quantify the margin.
 (D) GAP-SIDE LATTICE certificate tau = k/7 + j/m7 on a DIRECTED tight generator (sets built to
     be near 1/14). Does the 7-adic lattice always furnish a witness?
"""
import sys, random, time
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)
C14 = F(1, 14)

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
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
    b = F(0); at = None
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v; at = t
    return b, at
def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S, 0) == 1
def classify(S):
    S = sorted(set(S)); Vmin = min(S); Vmax = max(S)
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'
def nearest_k7(t):
    best=None
    for k in range(0,7):
        d = abs(t - F(k,7)); d = min(d, 1-d)
        if best is None or d < best[1]: best=(k,d)
    return best

# ------------------------------------------------------------------
print("="*84)
print("(A) KNOWN TIGHT SETS: exact M, 7-adic position of tau*")
print("="*84)
Sstar = [1,2,3,5,7,8,9,10,11,12,13,38,42]
m, at = Mval(Sstar)
print(f"  S* = {Sstar}")
print(f"     covering={is_cov(Sstar)} primitive={primitive(Sstar)} class={classify(Sstar)}")
print(f"     M(S*) = {m} = {float(m):.6f}  (M*14 = {float(m*14):.4f}; 2/23={float(F(2,23)):.6f})")
k7, d7 = nearest_k7(at)
print(f"     tau* = {at} = {float(at):.6f}; nearest k/7 = {k7}/7, dist = {d7} = {float(d7):.5f}")
print(f"     multiples of 7 in S*: {[v for v in Sstar if v%7==0]}")
print()
print("  THM-524 family A={1..11,13} + {84k}:")
A = list(range(1,12))+[13]
for kk in range(1,5):
    S = A+[84*kk]
    m, at = Mval(S)
    k7,d7 = nearest_k7(at)
    print(f"    k={kk} w={84*kk}: M={m}={float(m):.6f} (=7k/(84k+5)? {m==F(7*kk,84*kk+5)}) "
          f"tau*={at} near {k7}/7 dist {float(d7):.5f}; mult7={[v for v in S if v%7==0]}")

# ------------------------------------------------------------------
print()
print("="*84)
print("(B)+(C) THE WINDOW-COLLAPSE CERTIFICATE (constructive global witness near k/7)")
print("="*84)
print("""  For tau = k/7 + s: non-mult-7 runners safe iff |s| <= 1/(14 V*), V* = max non-mult7 runner.
  Inside, reduce to mult-of-7 sub-system {7 w_i}: need ||7 w_i s|| >= 1/14, i.e. with u=7s,
  ||w_i u|| >= 1/14 and |u| <= 1/(2 V*). We compute the EXACT widest safe arc of {w_i} at gap
  1/14 and check it meets the window |u| <= 1/(2 V*). If yes -> GLOBAL WITNESS, M(S) >= 1/14.""")

def mult7_subsystem(S):
    return sorted(v//7 for v in S if v % 7 == 0)   # the w_i
def Vstar(S):
    nm = [v for v in S if v % 7 != 0]
    return max(nm) if nm else 1

# exact safe components of a runner-set W at gap h on the full circle
def safe_components(W, h=C14):
    iv=[]
    for u in W:
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

def window_collapse_witness(S):
    """Return (has_witness, detail). A witness u near 0 (so s=u/7) clearing {w_i} with |u|<=1/(2V*)."""
    W = mult7_subsystem(S)
    if not W: return False, "no mult of 7 (not covering?)"
    Vst = Vstar(S)
    halfwin = F(1, 2*Vst)        # |u| <= 1/(2 V*)
    # safe set of W at gap 1/14 on circle; we need a safe point u in (-halfwin, halfwin) mod 1.
    sc = safe_components(W)
    # the relevant region around 0: u in [0,halfwin) U (1-halfwin,1). Check overlap with safe arcs.
    for (a,b) in sc:
        # overlap with [0, halfwin]
        lo = max(a, F(0)); hi = min(b, halfwin)
        if lo < hi: return True, (W, Vst, halfwin, ('u-arc', lo, hi))
        # overlap with [1-halfwin, 1]
        lo2 = max(a, 1-halfwin); hi2 = min(b, F(1))
        if lo2 < hi2: return True, (W, Vst, halfwin, ('u-arc', lo2, hi2))
    return False, (W, Vst, halfwin, "no safe u in window")

for nm, S in [("S*", Sstar), ("A+84", A+[84]), ("A+168", A+[168])]:
    has, det = window_collapse_witness(S)
    W = mult7_subsystem(S); Vst = Vstar(S)
    print(f"  {nm}: mult7 sub {[7*w for w in W]} -> w={W}, V*={Vst}, half-window 1/(2V*)={F(1,2*Vst)}")
    print(f"       window-collapse global witness exists: {has}   {det if not has else 'WITNESS '+str(det[3])}")

# ------------------------------------------------------------------
print()
print("="*84)
print("(C) UNCONDITIONAL r=1 and r=2 mult-of-7 sub-cases")
print("="*84)
print("""  r=1: single mult of 7, = 7w. Need ||w u||>=1/14 with |u|<=1/(2V*). The safe set of a
  single runner w (gap 1/14) is the complement of teeth of half-width 1/(14w) at j/w. Around
  u=0 the nearest tooth edge is at 1/(14w); so the safe arc (1/(14w), 1/w - 1/(14w)) starts at
  u=1/(14w). It meets |u|<=1/(2V*) iff 1/(14w) < 1/(2V*), i.e. V* < 7w = the mult-of-7 itself.
  Since 7w is IN S and V* is the max NON-mult-7 runner, this holds whenever 7w > V*, i.e. the
  (single) multiple of 7 is the largest runner -- OR more weakly 7w > V*. Quantify.""")
print("  r=1 closes iff V* < 7w (the unique mult-of-7 exceeds all non-mult-7 runners).")
print("  When 7w <= V*: the single tooth at u=0 has half-width 1/(14w) >= 1/(2V*) -- window may be")
print("  swallowed; need the NEXT safe arc (1/w - 1/(14w), ...) but that needs |u| ~ 1/w which may")
print("  exceed the window. So r=1 is NOT unconditional; it needs the size relation V* < 7w.")
print()
print("""  r=2: {7w1, 7w2}. Two-runner LRC at gap 1/14: a safe u near 0 exists iff the two teeth at
  u=0 (half-widths 1/(14 w1), 1/(14 w2)) leave room before |u|=1/(2V*). The widest safe arc of
  a 2-set {w1,w2} (gap 1/14) has a known closed form; we compute it EXACTLY per tight set.""")

# Test window-collapse over a DIRECTED tight generator + report closure rate by r.
def gen_tight(seed=0, target=2000):
    """Generate covering primitive S3 sets biased toward TIGHTNESS: small P (so few cluster
       members carry the covering obligations), and a 2-3 element large cluster near a multiple
       structure -- the regime where M approaches 1/14."""
    rng = random.Random(seed); out=[]; tries=0
    base=list(range(1,14))
    while len(out)<target and tries<target*500:
        tries+=1
        # small part: keep many of 1..13 to stay covering-cheap, drop 1-4
        ndrop = rng.choice([1,2,3,4])
        drop = rng.sample(base, ndrop)
        P = [v for v in base if v not in drop]
        c = 13 - len(P)
        if c < 1: continue
        # large cluster: pick c values > 13, ensure covering of missing q
        miss = [q for q in range(2,15) if not any(v%q==0 for v in P)]
        # bias: include a multiple of 7 and a tight pairing (like 38,42 for S*)
        Vc = rng.randint(15, 120)
        cluster=set()
        # force a multiple of 14 (covering needs mult of 7 and 14)
        cluster.add(14*rng.randint(2,8))
        while len(cluster) < c:
            cluster.add(rng.randint(15, Vc+60))
        cluster=sorted(cluster)
        if len(cluster)!=c: continue
        S = sorted(set(P)|set(cluster))
        if len(S)!=13: continue
        if not primitive(S) or not is_cov(S) or classify(S)!='S3': continue
        out.append(S)
    return out

tights = gen_tight(seed=11, target=2500)
print(f"  generated {len(tights)} covering primitive S3 'tight-biased' sets.")
by_r = Counter(); closed_by_r = Counter()
viol_window = []
for S in tights:
    r = sum(1 for v in S if v%7==0)
    by_r[r]+=1
    has,_ = window_collapse_witness(S)
    if has: closed_by_r[r]+=1
    else: viol_window.append(S)
print("  window-collapse closure by r (#mult-of-7):")
for r in sorted(by_r):
    print(f"    r={r}: {closed_by_r[r]}/{by_r[r]} closed by window-collapse "
          f"({100*closed_by_r[r]/by_r[r]:.1f}%)")
print(f"  total window-collapse-closed: {sum(closed_by_r.values())}/{len(tights)} "
      f"({100*sum(closed_by_r.values())/len(tights):.1f}%)")

# ------------------------------------------------------------------
print()
print("="*84)
print("(D) FOR THE WINDOW-MISSED SETS: is M still >= 1/14 (exact)? floor + gap-side lattice")
print("="*84)
print(f"  window-collapse MISSED {len(viol_window)} sets. Compute exact M on them and the")
print(f"  realized floor; also test the 7-adic lattice tau=k/7+j/m7 as an alternate witness.")
missminM = F(10); missminS=None; below=0
# to keep it fast, only exact-M the window-missed sets (these are the genuinely interesting ones)
sample_missed = viol_window[:600]
t0=time.time()
for S in sample_missed:
    m,_ = Mval(S)
    if m < missminM: missminM=m; missminS=S
    if m < C14: below+=1
print(f"  exact M on {len(sample_missed)} window-missed sets [{time.time()-t0:.1f}s]: below 1/14 = {below}")
print(f"  min exact M over window-missed = {missminM} = {float(missminM):.6f} (M*14={float(missminM*14):.4f})")
print(f"    at S = {missminS}")

def gap_side_certificate(S):
    m7list = sorted(v for v in S if v % 7 == 0)
    best = F(0); bestat=None
    for k in range(1,7):
        for m7 in m7list:
            for j in range(-12,13):
                if j==0: continue
                t = (F(k,7) + F(j,m7)) % 1
                mn = min(nrm(v*t) for v in S)
                if mn > best: best=mn; bestat=(k,m7,j,t)
    return best, bestat

cok=0; cworst=F(10); cworstS=None
for S in sample_missed[:300]:
    c,_ = gap_side_certificate(S)
    if c >= C14: cok+=1
    if c < cworst: cworst=c; cworstS=S
print(f"  gap-side k/7+j/m7 lattice certificate >= 1/14 on {cok}/{min(300,len(sample_missed))} window-missed sets")
print(f"  worst lattice certificate = {cworst} = {float(cworst):.6f} (cert*14={float(cworst*14):.4f}) at {cworstS}")

# also directly on S*:
c, cat = gap_side_certificate(Sstar)
print(f"\n  S* gap-side lattice certificate = {c} = {float(c):.6f} (>=1/14? {c>=C14}) at {cat}")
print(f"  S* exact M = {Mval(Sstar)[0]} -> the lattice {'MATCHES' if c==Mval(Sstar)[0] else 'differs from'} exact M")

print()
print("="*84)
print("DONE.")
print("="*84)
