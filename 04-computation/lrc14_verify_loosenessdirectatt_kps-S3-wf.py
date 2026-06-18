#!/usr/bin/env python3
"""
lrc14_verify_loosenessdirectatt_kps-S3-wf.py  (LEAN, flushing)

ADVERSARIAL VERIFICATION of the "looseness-direct" claimed advance on LRC(14) case S3.

Claims under test:
  (T1) TWO-GAP BAND-FIT LEMMA: P (small) union L (cluster). If exists tau*,int m>=1 with
       (A) frac(u tau*) in [1/14,13/14] for all u in P, and
       (B) [V tau*,(V+s) tau*] subset (m+1/14,m+13/14),
       then tau* safe for all S => M(S)>=1/14.  Suff pair (B1) s*tau*<=6/7 + (B2) gap m fits.
  (T2) GAP-VECTOR IDENTITY: exists nonempty A(k)  <=>  M(S)>1/14.
  (T3) FINITE FLOOR: every S3 set in {1..22} has M>=1/12, equality at {1,2,3,4,10..18}.
  (T4)/(T5): hunt ANY covering primitive S3 13-set with M<1/14.  ONE => LRC(14) refuted.
"""
from fractions import Fraction as F
from math import gcd, floor
from itertools import combinations
import random, sys

def P(*a): print(*a, flush=True)

# ---------------- core tools ----------------
def nrm(x):
    r = x - int(x); r = r+1 if r<0 else r
    return r if r <= F(1,2) else 1-r

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2,15))

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2):
                        C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C

def Mval_full(S):
    """returns (M, argmax tau)."""
    b = F(0); bt=None
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v; bt=t
    return b, bt

def Mval(S):
    return Mval_full(S)[0]

def primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1

def is_S3(S):
    large = [v for v in S if v > 13]
    if len(large) < 2: return False
    return max(S) >= 13*min(S)

def split_PL(S):
    return sorted(v for v in S if v <= 13), sorted(v for v in S if v > 13)

H = F(1,14)

# Independent exhaustive-breakpoint M (superset of cand) for AUDIT on small-speed sets.
def Mval_exhaustive(S):
    S = sorted(set(S)); C = set([F(1,2)])
    for v in S:
        for k in range(v):
            p = F(2*k+1, 2*v)
            if p <= F(1,2): C.add(p)
    for i in range(len(S)):
        for j in range(len(S)):
            if i==j: continue
            a,b = S[i],S[j]
            for d in (a+b, abs(a-b), a, b):
                if d>0:
                    k=1
                    while F(k,d) < 1:
                        C.add(F(k,d)); k+=1
    best=F(0)
    for t in C:
        t=t%1
        v=min(nrm(x*t) for x in S)
        if v>best: best=v
    return best

# True M via uniform grid (valid only as a LOWER bound; used to catch cand() missing the max).
def Mval_grid(S, denom):
    best=F(0)
    for i in range(denom):
        t=F(i,denom)
        v=min(nrm(x*t) for x in S)
        if v>best: best=v
    return best

# ====================================================================
P("="*70); P("TEST A: audit cand() Mval against exhaustive breakpoints + grid (small sets)")
P("="*70)
random.seed(12345)
disc=0; gridbeat=0; nA=0
for trial in range(400):
    n=random.randint(3,7)
    S=sorted(random.sample(range(1,26),n))
    if not primitive(S): continue
    nA+=1
    m1=Mval(S); m2=Mval_exhaustive(S)
    if m1!=m2:
        disc+=1; P(f"  DISCREPANCY cand vs exhaustive: S={S} cand={m1} exh={m2}")
    # grid lower-bound check: grid must NOT exceed cand (else cand missed the true max)
    if max(S)<=26 and n<=6:
        dens=1
        for v in S: dens=dens*v//gcd(dens,v)
        if dens<=50000:
            mg=Mval_grid(S, dens)
            if mg>m1:
                gridbeat+=1; P(f"  GRID BEATS cand: S={S} cand={m1} grid={mg}")
P(f"  sets audited: {nA}")
P(f"  cand-vs-exhaustive discrepancies: {disc}  (expect 0 => cand() is the correct M)")
P(f"  grid-beats-cand events: {gridbeat}  (expect 0 => cand() never undercounts)")

# ====================================================================
P(""); P("="*70); P("TEST B: two-gap lemma soundness + (B1)+(B2) -> M>=1/14 hunt (bounded speed)")
P("="*70)
def gap_safe(u,tau,h=H):
    f=(u*tau)%1; return h<=f<=1-h
def lemma_safe(S,tau,h=H):
    return all(nrm(v*tau)>=h for v in S)

def b1b2_witness(Pp,L,h=H,max_m=80):
    V=min(L); Vs=max(L); s=Vs-V
    for m in range(1,max_m+1):
        lo=F(m)+h; hi=F(m)+1-h
        tlo=lo/V; thi=hi/Vs
        if tlo>thi: continue
        if s>0: thi=min(thi, F(6,7)/s)
        if tlo>thi: continue
        cands=set([tlo,thi])
        for u in Pp:
            jlo=int(tlo*u)-1; jhi=int(thi*u)+1
            if jhi-jlo>200:  # safety cap; window is tiny so this never triggers legitimately
                jhi=jlo+200
            for j in range(max(0,jlo),jhi+1):
                for pt in (F(j)+h,F(j)+1-h):
                    t=pt/u
                    if tlo<=t<=thi: cands.add(t)
        cs=sorted(cands); more=set()
        for a,b in zip(cs,cs[1:]): more.add((a+b)/2)
        for t in sorted(cands|more):
            if all(gap_safe(u,t,h) for u in Pp):
                if F(m)+h<=V*t and Vs*t<=F(m)+1-h and s*t<=F(6,7):
                    return (t,m)
    return None

def gen_S3(num, seed, maxV):
    random.seed(seed); out=[]; tries=0
    while len(out)<num and tries<num*300:
        tries+=1
        ks=random.randint(3,11)
        Pp=sorted(random.sample(range(1,14),ks))
        kL=13-len(Pp)
        if kL<2: continue
        V0=random.randint(14,maxV); sp=random.randint(14,45)
        L=set()
        while len(L)<kL: L.add(V0+random.randint(0,sp))
        S=sorted(set(Pp)|L)
        if len(S)!=13: continue
        if not primitive(S): continue
        if not is_covering(S): continue
        if not is_S3(S): continue
        out.append(S)
    return out

S3B=gen_S3(500, 999, 220)   # V0<=220 keeps Mval fast; structure/minimizers live here
P(f"  generated {len(S3B)} covering primitive S3 sets (V0<=220)")
lemv=0; b1b2hold=0; minM=None; argm=None
for idx,S in enumerate(S3B):
    if idx%100==0: P(f"   ...TEST B progress {idx}/{len(S3B)} running min M={minM}")
    Pp,L=split_PL(S)
    if len(L)>=2:
        w=b1b2_witness(Pp,L)
        if w:
            b1b2hold+=1
            tau,m=w
            if not lemma_safe(S,tau):
                lemv+=1; P(f"  !!! T1 REFUTED: S={S} tau={tau} not safe (min={min(nrm(v*tau) for v in S)})")
    Ms=Mval(S)
    if minM is None or Ms<minM: minM=Ms; argm=S
    if Ms<H: P(f"  !!! S3 COUNTEREXAMPLE: S={S} M={Ms}<1/14")
P(f"  (B1)+(B2)+(A) witnesses found: {b1b2hold}/{len(S3B)}")
P(f"  T1 logical violations (witness not actually safe): {lemv}  (expect 0)")
P(f"  min M over these S3 sets: {minM}={float(minM):.5f} at {argm}  (margin x{float(minM/H):.3f})")

# ====================================================================
P(""); P("="*70); P("TEST C: gap-vector identity T2 (exists A(k) nonempty  <=>  M>1/14)")
P("="*70)
random.seed(77); t2fail=0; nC=0
for trial in range(500):
    n=random.randint(3,8)
    S=sorted(random.sample(range(1,55),n))
    if not primitive(S): continue
    nC+=1
    M,bt=Mval_full(S)
    # exists nonempty A(k) (open) <=> M>H. Construct from argmax when M>H.
    if M>H:
        ok=all(H<((u*bt)%1)<1-H for u in S)
        if not ok:
            # argmax may sit on closed boundary even when M>H; nudge search over cand
            ok=False
            for t in cand(S):
                if all(H<((u*t)%1)<1-H for u in S): ok=True; break
        exists=ok
    else:
        # M<=H: claim no nonempty A(k). Double-check no cand tau is strictly interior.
        exists=any(all(H<((u*t)%1)<1-H for u in S) for t in cand(S))
    want=(M>H)
    if exists!=want:
        t2fail+=1; P(f"  T2 MISMATCH: S={S} M={M} exists={exists} want={want}")
P(f"  sets tested: {nC};  T2 mismatches: {t2fail}  (expect 0)")

# ====================================================================
P(""); P("="*70); P("TEST D: exhaustive S3 floor in window {1..22}")
P("="*70)
floor_min=None; floor_arg=None; below1_12=0; below1_14=0; checked=0
universe=range(1,23)
for S in combinations(universe,13):
    S=list(S)
    if not is_covering(S): continue
    if not primitive(S): continue
    if not is_S3(S): continue
    checked+=1
    M=Mval(S)
    if floor_min is None or M<floor_min: floor_min=M; floor_arg=S
    if M<F(1,12): below1_12+=1
    if M<H: below1_14+=1; P(f"  !!! BELOW 1/14: S={S} M={M}")
    if checked%2000==0: P(f"   ...checked {checked} S3 sets, running floor={floor_min}")
P(f"  S3 sets checked: {checked}")
P(f"  floor (min M): {floor_min}={float(floor_min):.6f}  argmin={floor_arg}")
P(f"  floor == 1/12 ? {floor_min==F(1,12)}")
P(f"  count below 1/12: {below1_12};  below 1/14: {below1_14}")
cm=[1,2,3,4,10,11,12,13,14,15,16,17,18]
P(f"  claimed minimizer {cm}: cov={is_covering(cm)} prim={primitive(cm)} S3={is_S3(cm)} M={Mval(cm)}")

# ====================================================================
P(""); P("="*70); P("TEST E: broad adversarial hunt for S3 with M<1/14 (bounded speed, completes)")
P("="*70)
random.seed(20260618); gmin=None; garg=None; brk=False; nE=0
def gen_adv():
    # speeds bounded (V0<=260) so Mval stays fast; the documented minimizers (156,182,
    # 10..18, 280-family) all live in this range, and the claim asserts scale-invariance.
    s=random.randint(0,4)
    if s==0:
        Pp=[1,2,3,4]; c=random.randint(10,250); L=[c+i for i in range(9)]; S=sorted(set(Pp)|set(L))
    elif s==1:
        Pp=list(range(1,12)); a=random.randint(14,250); b=a+random.randint(1,80); S=sorted(set(Pp)|{a,b})
    elif s==2:
        ks=random.randint(4,10); Pp=sorted(random.sample(range(1,14),ks)); kL=13-len(Pp)
        if kL<2: return None
        V0=random.randint(14,260); sp=random.randint(14,60); L=set()
        while len(L)<kL: L.add(V0+random.randint(0,sp))
        S=sorted(set(Pp)|L)
    elif s==3:
        Pp=sorted(random.sample(range(1,14),random.randint(5,9))); kL=13-len(Pp)
        if kL<2: return None
        V0=random.choice([156,182,140,168,210,280,154,196,120,260]); L=set()
        while len(L)<kL: L.add(V0+random.choice([0,14,28,2,4,7,21,35,42,3,6]))
        S=sorted(set(Pp)|L)
    else:
        ks=random.randint(3,9); Pp=sorted(random.sample(range(1,14),ks)); kL=13-len(Pp)
        if kL<2: return None
        V0=random.randint(50,260); sp=random.randint(14,50); L=set()
        while len(L)<kL: L.add(V0+random.randint(0,sp))
        S=sorted(set(Pp)|L)
    if len(S)!=13: return None
    if not primitive(S): return None
    if not is_covering(S): return None
    if not is_S3(S): return None
    return S
target=4000; tries=0
while nE<target and tries<target*200:
    tries+=1
    S=gen_adv()
    if S is None: continue
    nE+=1
    M=Mval(S)
    if gmin is None or M<gmin: gmin=M; garg=S
    if M<H:
        brk=True; P(f"  !!!! LRC(14) COUNTEREXAMPLE: S={S} M={M}={float(M):.6f}<1/14 !!!!")
    if nE%500==0: P(f"   ...tested {nE}, running min M={gmin}={float(gmin):.5f}")
P(f"  adversarial S3 sets tested: {nE}")
P(f"  min M: {gmin}={float(gmin):.6f} at {garg}  (margin x{float(gmin/H):.4f})")
P(f"  any break below 1/14: {brk}")

# ---- separate HIGH-SPEED scale-invariance probe (small count, big V0) ----
P(""); P("  -- high-speed probe (V0 in [800,5000], 120 sets, slow Mval) --")
random.seed(31337); hmin=None; harg=None; nH=0; hbrk=False
while nH<120:
    ks=random.randint(4,9); Pp=sorted(random.sample(range(1,14),ks)); kL=13-len(Pp)
    if kL<2: continue
    V0=random.randint(800,5000); sp=random.randint(14,50); L=set()
    while len(L)<kL: L.add(V0+random.randint(0,sp))
    S=sorted(set(Pp)|L)
    if len(S)!=13 or not primitive(S) or not is_covering(S) or not is_S3(S): continue
    nH+=1
    M=Mval(S)
    if hmin is None or M<hmin: hmin=M; harg=S
    if M<H: hbrk=True; P(f"  !!!! HIGH-SPEED COUNTEREXAMPLE: S={S} M={M}<1/14")
    if nH%30==0: P(f"   ...high-speed {nH}/120 running min={hmin}={float(hmin):.5f}")
P(f"  high-speed sets tested: {nH}; min M={hmin}={float(hmin):.6f} at {harg} (>=1/14? {hmin>=H})")

# claimed specific minimizer {1..11,156,182}
chk=sorted(set(range(1,12))|{156,182})
Mc,tc=Mval_full(chk)
P(f"  {{1..11,156,182}}: cov={is_covering(chk)} prim={primitive(chk)} S3={is_S3(chk)} M={Mc}={float(Mc):.5f} (=5/61? {Mc==F(5,61)}) tau={tc}")

P(""); P("="*70); P("SUMMARY")
P("="*70)
P(f"  T1 logical violations: {lemv}")
P(f"  T2 mismatches: {t2fail}")
P(f"  T3 floor {{1..22}}: {floor_min}={float(floor_min):.6f}  (>=1/14? {floor_min>=H}; ==1/12? {floor_min==F(1,12)})")
P(f"  T5 adversarial min M: {gmin}={float(gmin):.6f}  (>=1/14? {gmin>=H})")
P(f"  high-speed probe min M: {hmin}={float(hmin):.6f}  (>=1/14? {hmin>=H})")
P(f"  NO break below 1/14 anywhere: {(not brk) and (not hbrk) and (minM>=H) and (floor_min>=H) and (gmin>=H) and (hmin>=H)}")
