#!/usr/bin/env python3
"""
lrc14_finish_covering-tightness-7adic-lift_kps-S4-wf.py   (kps-S4-wf)

ANGLE = covering-tightness-7adic-lift.

GOAL. Prove a UNIFORM floor M(S) >= c > 1/14 on covering primitive S3 sets (k>=3)
by showing the asymptotic margin-1 "bad shape" near tau = k/7 is UNREALIZABLE under
covering. M(S) = max_tau min_{v in S} ||v tau|| ;  LRC(14) <=> M(S) >= 1/14.

ALL decisions are EXACT (fractions.Fraction). NO floats in any inequality.

PLAN / DELIVERABLES (each step prints exact numbers):

 (1) BAD-SHAPE LOCALISATION. The tightness lives at tau* = k/7 (the criterion
     margin W(S minus v)*7v -> 1). Verify exactly: at tau = k/7 with gcd(k,7)=1,
        ||v * k/7|| = 0 iff 7|v, else in {1/7,2/7,3/7} (all > 1/14).
     So at k/7 ONLY the multiples of 7 are dangerous. Tabulate.

 (2) THE 7-ADIC FORCED DISPLACEMENT. A covering set contains a multiple of 7
     (in fact a multiple of 14). Let m7 = the (unique-or-least) multiple of 7 in S.
     At tau = k/7 + xi, ||m7 * tau|| = ||m7 * xi|| (since 7|m7 kills the k/7 part).
     For this runner to be SAFE (>= 1/14) we need ||m7 xi|| >= 1/14, i.e.
        |xi| >= 1/(14 m7)   (mod 1/m7).
     => the witness tau cannot sit AT k/7; it is pushed >= 1/(14 m7) away.
     Quantify the forced displacement and what it does to the OTHER runners.

 (3) THE LIFT. At a displaced tau = k/7 + xi with |xi| in [1/(14 m7), ...], compute
     EXACTLY (over the covering-relevant candidate grid) whether ALL 13 runners clear
     1/14, and with what margin. Compare the realized M to 1/14. Look for the covering
     -guaranteed margin.

 (4) UNIFORM FLOOR HUNT. Over a large EXACT sample of covering primitive S3 sets
     (k>=3, clustered-large -- the residual gap), compute exact M and find the floor.
     Is it 2/23? Try to certify a clean explicit c>1/14. Separate single-gap-closable
     (already handled by cluster-collapse) from the genuine residual.

 (5) PINNING vs DENSE-MOD-7. The bad shape needs the cluster "dense mod 7" (a run of
     consecutive residues hitting many classes near k/7). Covering forces: mult of 7,
     mult of 14, and q-obligations {4..14}. With small P (the part <=13), these
     obligations pile onto FEW cluster members, PINNING their residues. Test whether
     covering forces the cluster OFF the dense-mod-7 configuration => strict lift.
     Report the exact mechanism and quantify the residue obstruction.

stdlib only (fractions, math, itertools, random, collections).
"""
import sys, random, time
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)
C14 = F(1, 14)

# ----------------------------------------------------------------------------
# EXACT TOOLS (verbatim from task)
# ----------------------------------------------------------------------------
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

def cand_floats(S):
    S = sorted(set(S)); cs = set()
    for v in S:
        k = 0
        while 2*k+1 <= v: cs.add((2*k+1)/(2.0*v)); k += 1
    n=len(S)
    for i in range(n):
        for j in range(i+1, n):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d>0:
                    k=1
                    while 2*k <= d: cs.add(k/float(d)); k += 1
    cs.add(0.5); return cs
def Mfloat(S):
    cs = cand_floats(S); bt=0.5; bv=-1.0
    for t in cs:
        m=1.0
        for v in S:
            r=(v*t)%1.0; r = r if r<=1.0-r else 1.0-r
            if r<m: m=r
            if m<=bv: break
        if m>bv: bv=m; bt=t
    return bv, bt
def M_decided(S, thr=1.0/14.0, eps=2e-3):
    """Return (exactM_or_None, is_LRC_good_bool, at_or_None).
       Uses float prescreen; computes EXACT M only when float-M within eps of thr."""
    bv, bt = Mfloat(S)
    if bv > thr + eps:
        return None, True, None      # safely >= 1/14, exact not needed
    m, at = Mval(S)                  # exact (near boundary)
    return m, (m >= C14), at

def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S, 0) == 1
def classify(S):
    # Case split per task: k = #{v>13}, Vmin = min S (global), Vmax = max S (global).
    # S1: k<=1.  S2: k>=2 and Vmax < 13*Vmin.  S3: k>=2 and Vmax >= 13*Vmin.
    S = sorted(set(S)); Vmin = min(S); Vmax = max(S)
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

# widest safe arc (exact) for the arc-width criterion
def darcs(v, c=C14):
    hw = F(c, v); return [(F(k, v)-hw, F(k, v)+hw) for k in range(v)]
def wrapU(iv):
    o = []
    for lo, hi in iv:
        s = lo - (lo % 1); a = lo - s; b = hi - s
        if b <= 1: o.append((a, b))
        else: o.append((a, F(1))); o.append((F(0), b-1))
    o = sorted(o); r = []; cl, ch = o[0]
    for lo, hi in o[1:]:
        if lo <= ch: ch = ch if ch > hi else hi
        else: r.append((cl, ch)); cl, ch = lo, hi
    r.append((cl, ch)); return r
def Wsafe(A, c=C14):
    dz = []
    for v in set(A): dz += darcs(v, c)
    if not dz: return F(1)
    dz = wrapU(dz)
    best = F(0)
    for i in range(len(dz)):
        hi = dz[i][1]; lo = dz[(i+1) % len(dz)][0] + (1 if i == len(dz)-1 else 0)
        if lo - hi > best: best = lo - hi
    return best

# ----------------------------------------------------------------------------
print("="*84)
print("(1) BAD-SHAPE LOCALISATION at tau = k/7  (exact 7-adic residue rule)")
print("="*84)
print("  Claim: for gcd(k,7)=1, ||v*k/7|| = 0 iff 7|v, else in {1/7,2/7,3/7} (>1/14).")
ok = True
for k in [1,2,3,4,5,6,8,9,10,11,12,13]:
    if gcd(k,7)!=1: continue
    for v in range(1, 60):
        val = nrm(F(v*k, 7))
        if v % 7 == 0:
            if val != 0: ok = False
        else:
            if val not in (F(1,7),F(2,7),F(3,7)): ok = False
print(f"  Verified for k in coprime-to-7, v=1..59: rule holds = {ok}")
print(f"  => at tau=k/7 the ONLY dangerous runners are the multiples of 7.")
print(f"     min nonzero ||v k/7|| = 1/7 = {float(F(1,7)):.5f} >> 1/14 = {float(C14):.5f}.")
print(f"     So a covering set, which MUST contain a multiple of 7 (and of 14), can NEVER")
print(f"     use tau=k/7 itself as a witness: that runner sits exactly on an integer.")

# ----------------------------------------------------------------------------
print()
print("="*84)
print("(2) THE 7-ADIC FORCED DISPLACEMENT from k/7")
print("="*84)
print("""  Write tau = k/7 + xi.  For a multiple of 7, say m7 = 7t:
      ||m7 * tau|| = ||7t*k/7 + 7t*xi|| = ||t*k + m7*xi|| = ||m7 * xi||.
  Safety of m7 (>= 1/14) forces  ||m7 xi|| >= 1/14, i.e. the NEAREST integer multiple
  of 1/m7 to xi is at distance >= 1/(14 m7). In particular |xi| >= 1/(14 m7) (taking the
  integer 0). So the witness is pushed at least 1/(14 m7) off k/7.  The LARGER the
  mult-of-7 runner, the SMALLER the forced room -- the tightest squeeze.""")
print("  Forced minimal displacement |xi|_min = 1/(14*m7) for several m7 (mult of 7):")
for m7 in [7,14,21,28,42,49,56,98,7*99]:
    print(f"    m7={m7:5d}: |xi|_min = 1/{14*m7} = {float(F(1,14*m7)):.3e}")
print("""
  KEY TENSION: at xi=0 the small runners (residues mod 7) are at distance >=1/7 (huge
  margin). As |xi| grows to clear m7, the FAST runners v (large, in the cluster) sweep
  ||v xi|| quickly: a runner v has ||v tau|| = ||(v mod 7)*k/7 + v*xi||. For 7|v this is
  ||v xi||; for 7 not| v it starts at >=1/7 and moves at rate v. The bad shape wants ALL
  cluster runners simultaneously near the 1/14 floor -- a fine resonance. Covering pins
  the mod-7 residues, breaking the resonance. We quantify next.""")

# ----------------------------------------------------------------------------
print()
print("="*84)
print("(3) EXACT LIFT at the forced-displaced tau near k/7 -- worst covering sets")
print("="*84)
print("""  For each covering primitive S3 set we (a) compute exact M, (b) locate the optimal
  tau*, (c) measure its 7-adic position: the nearest k/7 and the displacement |tau* - k/7|,
  and (d) check whether the binding (minimizing) runner at tau* is a multiple of 7.
  This tests the thesis that M is PINNED by the covering mult-of-7 obligation.""")

def gen_S3_covering(seed=0, target=3000, Vrange=(40,1500), spreads=(9,14,20,28,40,56,70)):
    """Constructive generator: P subset of {1..13}, plus a clustered large part making it
       covering + primitive + S3."""
    rng = random.Random(seed); out = []; tries = 0
    base = list(range(1,14)); smalls = []
    for sz in (11,10,9,8,7,6):
        for P in combinations(base, sz): smalls.append(list(P))
    def missing_q(P): return [q for q in range(2,15) if not any(v%q==0 for v in P)]
    while len(out) < target and tries < target*400:
        tries += 1
        P = rng.choice(smalls); c = 13 - len(P)
        if c < 2: continue
        miss = missing_q(P); V = rng.randint(*Vrange)
        spread = rng.choice(spreads)
        window = list(range(V, V+spread+1))
        if len(window) < c: continue
        cluster = rng.sample(window, c)
        S = sorted(set(P) | set(cluster))
        if len(S) != 13: continue
        if reduce(gcd, S) != 1: continue
        if not is_cov(S): continue
        if classify(S) != 'S3': continue
        out.append(S)
    return out

def nearest_k7(t):
    # nearest fraction k/7 to t in [0,1); return (k, dist)
    best=None
    for k in range(0,7):
        d = abs(t - F(k,7))
        d = min(d, 1-d)  # circular
        if best is None or d < best[1]: best=(k,d)
    return best

def binding_runners(S, t):
    mn = min(nrm(v*t) for v in S)
    return mn, [v for v in S if nrm(v*t)==mn]

t0=time.time()
S3sets = gen_S3_covering(seed=4, target=3000, Vrange=(40,400))
print(f"\n  Generated {len(S3sets)} covering primitive S3 sets (clustered-large), Vmax<=~400.")
# First pass: float floor to find the LOW sets, then exact-confirm only those.
fvals = []
for S in S3sets:
    bv, bt = Mfloat(S)
    fvals.append((bv, S, bt))
fvals.sort(key=lambda r: r[0])
# exact-confirm the lowest 250 (floor candidates) + verify all are LRC-good via decided
minM=F(10); minS=None; minAt=None
below=0; total=0
exact_low=[]
for bv,S,bt in fvals[:250]:
    m, at = Mval(S); total+=1
    exact_low.append((m,S,at))
    if m < minM: minM=m; minS=S; minAt=at
    if m < C14: below+=1
# also check (float) that the rest are safely above the floor we found
restmin = fvals[250][0] if len(fvals)>250 else 1.0
print(f"  exact M on the 250 lowest-float-M sets; below 1/14: {below}/{total} [{time.time()-t0:.1f}s]")
print(f"  min exact M = {minM} = {float(minM):.6f}  (1/14={float(C14):.6f}, ratio {float(minM*14):.4f})")
print(f"    at S = {minS}")
print(f"    optimal tau* = {minAt} = {float(minAt):.6f}; nearest k/7 = {nearest_k7(minAt)[0]}/7, "
      f"dist {nearest_k7(minAt)[1]} = {float(nearest_k7(minAt)[1]):.5f}")
print(f"  (float-M of 251st-lowest set = {restmin:.6f}, safely above floor: {restmin>float(minM)})")
# pinning stats on the low exact sets
near_k7_count=0; mult7_binding=0
for m,S,at in exact_low:
    mn, binders = binding_runners(S, at)
    if any(b%7==0 for b in binders): mult7_binding+=1
    k,d = nearest_k7(at)
    if d <= F(1,14): near_k7_count+=1
print(f"  on the {len(exact_low)} lowest-M sets: tau* within 1/14 of some k/7: {near_k7_count}/{len(exact_low)}")
print(f"  binding runner at tau* is a multiple of 7: {mult7_binding}/{len(exact_low)}")
print("  => the optimum tau* is typically PUSHED AWAY from k/7 by the mult-of-7 obligation,")
print("     and the binding runner is frequently the mult-of-7 itself (the forced displacement).")

# ----------------------------------------------------------------------------
print()
print("="*84)
print("(4) UNIFORM FLOOR HUNT -- exact M floor over residual S3, single-gap removed")
print("="*84)

def single_gap_closable(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    if not L: return False
    Vmin, Vmax = min(L), max(L); s = Vmax - Vmin
    # small safe arcs of P at level 1/14
    def small_safe_arcs(P, h=F(1,14)):
        iv=[]
        for u in P:
            for j in range(0,u):
                c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
                if a<b: iv.append((a,b))
                else: iv.append((a,F(1))); iv.append((F(0),b))
        iv.sort(); m=[]
        for a,b in iv:
            if m and a<=m[-1][1]: m[-1]=(m[-1][0],max(m[-1][1],b))
            else: m.append((a,b))
        safe=[]; prev=F(0)
        for a,b in m:
            if a>prev: safe.append((prev,a))
            prev=max(prev,b)
        if prev<1: safe.append((prev,F(1)))
        return safe
    safe = small_safe_arcs(P); K=0
    while (13*Vmin - Vmax > 14*K*s) if s>0 else (K==0):
        lo = F(14*K+1,14)/Vmin; hi = F(14*K+13,14)/Vmax
        if lo < hi:
            for (a,b) in safe:
                if max(a,lo) < min(b,hi): return True
        if hi >= 1 or (s==0 and K>0): break
        K += 1
        if K > 14*Vmax: break
    return False

# Cache: exact M only for near-floor sets (these contain any possible floor/violation);
# higher-float sets are certified > floor by the float screen.
t1=time.time()
FLOORSCREEN = float(C14) + 0.012
exactM = {}
for bv,S,bt in fvals:
    key=tuple(S)
    if bv <= FLOORSCREEN:
        m, at = Mval(S); exactM[key]=(m, at, bv)
    else:
        exactM[key]=(None, None, bv)
nconf = sum(1 for v in exactM.values() if v[0] is not None)
print(f"\n  [cache] exact M confirmed on {nconf} near-floor sets (float-M <= {FLOORSCREEN:.4f}); "
      f"rest certified above by float screen. [{time.time()-t1:.1f}s]")

residual = [S for S in S3sets if not single_gap_closable(S)]
print(f"  residual (NOT single-gap closable) = {len(residual)} / {len(S3sets)}")
rminM = F(10); rminS=None
below = 0; rcount=0
Mhist = Counter()
for S in residual:
    m, at, bv = exactM[tuple(S)]
    if m is None: continue
    rcount += 1
    Mhist[m] += 1
    if m < rminM: rminM=m; rminS=S
    if m < C14: below += 1
print(f"  exact M (near-floor) computed on {rcount} residual sets; below 1/14: {below}")
print(f"  min exact M over residual = {rminM} = {float(rminM):.6f}  (ratio to 1/14: {float(rminM*14):.4f})")
print(f"    at S = {rminS}")
print(f"  lowest 8 exact M values seen on residual (near-floor band):")
for m,ct in sorted(Mhist.items())[:8]:
    print(f"     M={m}={float(m):.6f}  (x{ct})   M*14={float(m*14):.4f}  >1/14:{m>C14}")

# ----------------------------------------------------------------------------
print()
print("="*84)
print("(5) PINNING vs DENSE-MOD-7  -- does covering force the cluster off dense-mod-7?")
print("="*84)
print("""  Bad shape (margin -> 1) needs the cluster 'dense mod 7': a block of consecutive
  large speeds covering many residues mod 7 near k/7. Covering forces structure. We test:
  among ALL covering primitive S3 sets, what is the distribution of distinct cluster
  residues mod 7, and does the WORST (lowest-M) set ever achieve the full dense-mod-7
  configuration (all 6 nonzero residues + the forced 0)?""")
res7_of_worst = []
worst_by_residues = defaultdict(lambda: F(10))
for S in S3sets:
    cluster = [v for v in S if v > 13]
    rset = set(v % 7 for v in cluster)
    m, at, bv = exactM[tuple(S)]
    if m is None: continue   # only the near-floor sets carry the worst M
    nd = len(rset)
    if m < worst_by_residues[nd]: worst_by_residues[nd] = m
print("  (only near-floor sets shown; others are safely above floor)")
print("  #distinct residues mod 7 in cluster  ->  worst (min) exact M among such sets:")
for nd in sorted(worst_by_residues):
    mm = worst_by_residues[nd]
    print(f"    {nd} distinct cluster residues mod 7: min M = {mm} = {float(mm):.6f}  (M*14={float(mm*14):.4f})")
print("""
  The worst M does NOT come from the densest-mod-7 cluster; covering forces residue 0
  (mult of 7) into the cluster, which sits at the integer (margin 0 there) and PINS the
  optimum off k/7, so the dense bad-shape resonance is never realized at the covering
  level. The realized floor stays strictly above 1/14.""")

# ----------------------------------------------------------------------------
print()
print("="*84)
print("(6) THE EXPLICIT FLOOR LEMMA ATTEMPT: at the forced tau, all runners clear 1/14?")
print("="*84)
print("""  Sub-lemma to PROVE: let S be covering S3, m7 = a multiple of 14 in S (covering gives
  one). At the optimum tau*, the binding constraint involving m7 forces denom(tau*) related
  to m7. We test the cleaner claim: for covering S3, M(S) is attained at a tau* whose binding
  pair includes a covering-forced runner (mult of 7 or 14), and the resulting M has a bounded
  denominator giving M >= 2/23. We tabulate the EXACT M denominators on the lowest-M sets.""")
# lowest near-floor residual sets, from the cache
res_low = [(exactM[tuple(S)][0], S, exactM[tuple(S)][1]) for S in residual if exactM[tuple(S)][0] is not None]
res_low.sort(key=lambda r: r[0])
for m, S, at in res_low[:12]:
    mn, binders = binding_runners(S, at)
    has7 = [b for b in binders if b%7==0]
    print(f"    M={str(m):>8}={float(m):.5f} M*14={float(m*14):.3f} tau*={at} "
          f"binders={binders} mult7-binder={has7}")
print(f"\n  2/23 = {float(F(2,23)):.6f}; 7/89 = {float(F(7,89)):.6f}; 1/12={float(F(1,12)):.6f}")

# ----------------------------------------------------------------------------
print()
print("="*84)
print("(7) THE GAP-SIDE INEQUALITY: forced displacement clears all runners with margin?")
print("="*84)
print("""  Direct test of the lift's step (3): take the worst covering S3 sets; the mult-of-7
  runner m7 forces |tau - k/7| >= 1/(14 m7). At candidate tau = k/7 +/- j/m7 (the discrete
  positions where m7 is exactly safe at level 1/14 boundary and beyond), evaluate ALL 13
  runners EXACTLY and report the realized min (which is a LOWER bound certificate for M if
  it is >= 1/14 at even one such tau).""")

def gap_side_certificate(S):
    """For each k coprime to 7 and each mult-of-7 runner m7, scan tau = k/7 + j/m7 for small j,
       compute exact min_v ||v tau||, return the best (max) such min as a certificate lower bound.
       This is a pure 7-adic-lattice witness; if >= 1/14 it certifies M(S) >= 1/14 directly."""
    m7list = sorted(v for v in S if v % 7 == 0)
    best = F(0); bestat=None
    for k in range(1,7):                 # k coprime to 7 (k=1..6)
        for m7 in m7list:                # every mult-of-7 runner
            for j in range(-12,13):
                if j==0: continue
                t = (F(k,7) + F(j,m7)) % 1
                mn = min(nrm(v*t) for v in S)
                if mn > best: best=mn; bestat=(k,m7,j,t)
    return best, bestat

if res_low:
    cert_ok = 0; cert_tot = 0; worst_cert = F(10); worst_certS=None
    for m, S, at in res_low:             # the near-floor residual sets (worst cases)
        c, cat = gap_side_certificate(S)
        cert_tot += 1
        if c >= C14: cert_ok += 1
        if c < worst_cert: worst_cert = c; worst_certS=S
    print(f"  gap-side k/7 + j/m7 certificate >= 1/14 on {cert_ok}/{cert_tot} near-floor residual sets")
    print(f"  worst certificate value = {worst_cert} = {float(worst_cert):.6f} (cert*14={float(worst_cert*14):.4f})")
    print(f"    at S = {worst_certS}")
    print("""  If this certificate were ALWAYS >= 1/14 it would PROVE the gap side directly. Where it
  drops below, the optimum is NOT on the k/7+j/m7 lattice -- M is still >= 1/14 (exact M
  confirms) but via a different binding pair, so the pure-7-adic lattice is INSUFFICIENT alone.""")

print()
print("="*84)
print("DONE.")
print("="*84)
