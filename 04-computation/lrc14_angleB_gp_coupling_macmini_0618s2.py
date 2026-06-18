#!/usr/bin/env python3
"""
lrc14_angleB_gp_coupling_macmini_0618s2.py   (mac-mini-2026-06-18-S2)

ANGLE B — the G_P coupling for THM-527 (LRC(14) case S3).

GOAL: rho*(P,E) = meas{ x in G_P : maxgap{frac(e_i x)} > 2/7 } >= c0 > 0,
given Good_E = {x: maxgap>2/7} has meas >= c(k)>0 (Angle A) and meas(G_P) >= 7/858.
Danger: Good_E and G_P anti-correlated (Good concentrated where G_P empty).

This script proves/verifies the SCALE-DECOUPLING that resolves the coupling:

PROVED LEMMA (four-window / q-grid + Lipschitz).
  For ANY offset set E (0 in E, spread = maxE = max E), Good_E contains an interval
  around each of the four low-denominator rationals {0, 1/2, 1/3, 2/3}:
      near 0   : (0, 5/(7 maxE))                                [one-sided; reflection gives near 1]
      near 1/2 : (1/2 - 3/(28 maxE), 1/2 + 3/(28 maxE))
      near 1/3 : (1/3 - 1/(42 maxE), 1/3 + 1/(42 maxE))
      near 2/3 : (2/3 - 1/(42 maxE), 2/3 + 1/(42 maxE))
  PROOF: at x = a/q (q in {1,2,3}) the points {frac(e_i a/q)} sit on the q-grid
  {0,1/q,...,(q-1)/q}, so the circular max-gap is >= 1/q >= 1/3 > 2/7. Each point
  moves at speed <= maxE in x, so a base gap G0 (>= 1/q) realised between two consecutive
  occupied slots stays >= G0 - 2 maxE |x-a/q| (its two endpoints each drift <= maxE|delta|,
  and no third point can enter while maxE|delta| < 1/q). Hence > 2/7 for
  |x-a/q| < (1/q - 2/7)/(2 maxE), which gives the stated half-widths. (q=1: gap = 1,
  one point cluster, half-width 5/(7 maxE).)  [verified exhaustively below]

COUPLING REDUCTION (this is the whole point of Angle B).
  rho*(P,E)  >=  SUM_{(a,q) in {0,1/2,1/3,2/3}} meas( G_P  ∩  window_{a/q}(maxE) ).
  The RHS depends on E ONLY through the single scalar maxE, and on P through the FIXED
  small-speed safe set. The "different scales" decoupling is literal: G_P is built from
  low frequencies p <= 13 (coarse comb), the windows sit at the four coarsest rationals,
  and Good_E's fine structure (the cluster orbit frac(e_i x)) is irrelevant beyond maxE.
  ==> the coupling crux is a PURE statement about (P subset {1..13}, scalar maxE).

RESULTS (all exact rationals):
 (I)   verify the four-window lemma never fails (exhaustive small E + random bounded-spread).
 (II)  PROVABLE conservative floor: min over (P, consecutive E) of the four-window LB.
       Closes ALL consecutive cases EXCEPT a finite explicit list of 4 (P,k) pairs.
 (III) the 4 exceptions: exhibit certified rational safe+good sub-arcs (rho* > 0 exactly).
       (these survive on the SHARP near-1/2 or near-1/3 window, wider than the conservative one).
 (IV)  EXACT true floor (consecutive) and the correlation ratio R = rho*/(meas G_P * meas Good),
       bounded below => no destructive anti-correlation (quasi-independence holds with R >= ~0.068).
 (V)   broad bounded-spread scan: min EXACT rho* over perforated / spread shapes (no zeros).

HONEST STATUS:
 PROVED: the four-window lemma; rho* >= (four-window LB); the LB is positive for every
   (P, consecutive E) outside the 4 explicit exceptions; the 4 exceptions have rho* > 0 by
   exact certificate. => rho*(P, consecutive E) > 0 for ALL admissible P, UNCONDITIONALLY.
 CONJECTURE/OPEN (the residual crux, unchanged from THM-527-G): a UNIFORM c0 over the FULL
   bounded-spread shape space (not just consecutive), and the integer- vs real-offset passage.
   Part (V) is strong computational evidence (no zeros, R bounded below) but not a theorem.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(6182)
H=F(1,14); TWO7=F(2,7)

# ---------- exact small-speed safe set G_P ----------
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A,h=H):
    if not A: return [(F(0),F(1))]
    dz=merge([iv for u in A for iv in danger_arcs(u,h)]); safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def meas(arcs): return sum((b-a for a,b in arcs), F(0))
def intersect(A,B):
    out=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: out.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return out
def wrap_window(c, hw):
    lo=(c-hw)%1; hi=(c+hw)%1
    if lo<hi: return [(lo,hi)]
    return [(lo,F(1)),(F(0),hi)]

# ---------- exact good set of E: {x: maxgap{frac(e x)} > 2/7} ----------
def good_set_exact(E):
    E=sorted(set(E)); k=len(E)
    diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1)
    good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1,E[t]) for t in range(k))
        order=[e for _,e in pts]; floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx]; f_cur=floors[idx]
            if idx<k-1: e_nx=order[idx+1]; f_nx=floors[idx+1]; wrap=F(0)
            else: e_nx=order[0]; f_nx=floors[0]; wrap=F(1)
            A=F(e_nx-e_cur); Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>TWO7: good.append((x0,x1))
                continue
            xb=(TWO7-Cc)/A
            if A>0: lo=max(x0,xb); hi=x1
            else: lo=x0; hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return merge(good)

def maxgap_at(E,x):
    pts=sorted(set((F(e)*x)%1 for e in E)); g=F(0)
    for i in range(len(pts)-1): g=max(g,pts[i+1]-pts[i])
    g=max(g,(pts[0]+1)-pts[-1]); return g

# conservative four-window measures
WINDOWS=[(F(0),lambda m:F(5,7*m)),(F(1,2),lambda m:F(3,28*m)),
         (F(1,3),lambda m:F(1,42*m)),(F(2,3),lambda m:F(1,42*m))]
def four_window_LB(gp, maxE):
    tot=F(0)
    for c,hf in WINDOWS:
        tot+=meas(intersect(gp, wrap_window(c, hf(maxE))))
    return tot

# ============================================================================
print("="*80)
print("(I) PROVED-LEMMA CHECK: Good_E contains the four conservative windows (never fails)")
print("="*80)
allok=True; tested=0
# exhaustive small E
small_E=[]
for k in range(3,8):
    for rest in itertools.combinations(range(1,9),k-1):
        small_E.append([0]+list(rest))
# random bounded-spread E
for _ in range(4000):
    k=random.randint(3,13); spread=random.randint(k-1,40)
    small_E.append([0]+sorted(random.sample(range(1,spread+1),min(k-1,spread))))
for E in small_E:
    mE=max(E); tested+=1
    for c,hf in WINDOWS:
        hw=hf(mE)
        for fr in (F(-99,100),F(-1,2),F(1,2),F(99,100)):
            x=(c+fr*hw)%1
            if maxgap_at(E,x)<=TWO7:
                allok=False
                print(f"   FAIL E={E} center={c} hw={hw} fr={fr}")
print(f"   tested {tested} offset sets (exhaustive |E|<=7 spread<=8 + 4000 random spread<=40)")
print(f"   four-window lemma holds: {allok}")

# ============================================================================
print("\n"+"="*80)
print("(II) PROVABLE conservative floor over (P, consecutive E={0..k-1}), maxE=k-1")
print("     LB(P,k) = sum_q meas(G_P ∩ window_q(k-1))  <=  rho*(P, consecutive)")
print("="*80)
LB_zero=[]; glob_LB=(F(10),None,None)
for k in range(5,14):
    psz=13-k; mlb=(F(10),None)
    for P in itertools.combinations(range(1,14),psz):
        gp=safe_set(list(P)); lb=four_window_LB(gp,k-1)
        if lb<mlb[0]: mlb=(lb,P)
        if lb==0: LB_zero.append((k,P))
    print(f"   k={k:2d}: min provable LB = {mlb[0]} = {float(mlb[0]):.6f}  at P={mlb[1]}")
    if mlb[0]<glob_LB[0]: glob_LB=(mlb[0],mlb[1],k)
print(f"\n   GLOBAL provable LB floor = {glob_LB[0]} = {float(glob_LB[0]):.6f} (at k={glob_LB[2]} P={glob_LB[1]})")
print(f"   #(P,k) with conservative LB = 0 (the finite exception list): {len(LB_zero)}")
for k,P in LB_zero: print(f"        EXCEPTION: k={k} P={P}")

# ============================================================================
print("\n"+"="*80)
print("(III) THE EXCEPTIONS: exact rho* > 0 with certified rational safe+good sub-arcs")
print("="*80)
for k,P in LB_zero:
    gp=safe_set(list(P)); Gk=good_set_exact(list(range(k)))
    inter=intersect(Gk,gp)
    print(f"   k={k} P={list(P)}: EXACT rho* = {meas(inter)} = {float(meas(inter)):.6f} > 0")
    for a,b in inter[:4]:
        print(f"        certificate sub-arc ({a},{b}) center~{float((a+b)/2):.4f} width {float(b-a):.5f}")

# ============================================================================
print("\n"+"="*80)
print("(IV) EXACT true floor (consecutive) + correlation ratio R = rho*/(measG_P * measGood)")
print("     R bounded below => quasi-independence (no destructive anti-correlation)")
print("="*80)
glob_rho=(F(10),None,None); glob_R=(F(10),None,None)
for k in range(3,14):
    Gk=good_set_exact(list(range(k))); mG=meas(Gk); psz=13-k
    mr=(F(10),None); mR=(F(10),None)
    for P in itertools.combinations(range(1,14),psz):
        gp=safe_set(list(P)); mP=meas(gp)
        ex=meas(intersect(Gk,gp))
        R = ex/(mP*mG) if mP*mG>0 else F(1)
        if ex<mr[0]: mr=(ex,P)
        if R<mR[0]: mR=(R,P)
    print(f"   k={k:2d}: min rho*={float(mr[0]):.6f}({mr[0]}) P={mr[1]} | min R={float(mR[0]):.5f} P={mR[1]}")
    if mr[0]<glob_rho[0]: glob_rho=(mr[0],mr[1],k)
    if mR[0]<glob_R[0]: glob_R=(mR[0],mR[1],k)
print(f"\n   EXACT consecutive floor rho* = {glob_rho[0]} = {float(glob_rho[0]):.6f} (k={glob_rho[2]} P={glob_rho[1]})")
print(f"   MIN correlation ratio R      = {glob_R[0]} = {float(glob_R[0]):.5f} (k={glob_R[2]} P={glob_R[1]})")
print(f"   => R >= {float(glob_R[0]):.5f} > 0 on all consecutive cases: quasi-independence holds.")

# ============================================================================
print("\n"+"="*80)
print("(V) BROAD bounded-spread scan: min EXACT rho* over perforated / spread shapes (zeros?)")
print("="*80)
best=(F(10),None,None); zeros=0; ntest=0
for _ in range(3000):
    k=random.randint(5,11); psz=13-k
    P=sorted(random.sample(range(1,14),psz))
    spread=random.choice([k-1,k,k+1,k+3,2*k,3*k])
    rest=sorted(random.sample(range(1,spread+1), min(k-1,spread)))
    E=[0]+rest
    if len(set(E))<3: continue
    ntest+=1
    gp=safe_set(P); Gk=good_set_exact(E)
    ex=meas(intersect(Gk,gp))
    if ex<best[0]: best=(ex,P,E)
    if ex==0: zeros+=1
print(f"   tested {ntest} (P, bounded-spread E) pairs")
print(f"   min EXACT rho* found = {best[0]} = {float(best[0]):.6f}  at P={best[1]} E={best[2]}")
print(f"   #(rho* == 0): {zeros}")
print("\nDONE.")
