#!/usr/bin/env python3
"""ROUTE E / THREAD E -- the CUT/CYCLE split of measS7 (opus-2026-06-20-S2).

GOAL (THREAD E deliverable). THM-559 split c3 into a CUT-space 2-body Ising energy
(regular = ground state, provably extremal) and a CYCLE-space many-body remainder.
measS7 (the LRC(14) Z/7 cover) is the analogous object on the OFFSET side. We give an
EXACT, computable CUT/CYCLE split of measS7 and ask whether the CUT part is consec-
extremal by a single-particle / 2-body argument, with the CYCLE correction bounded.

THE SPLIT (exact, via the proved Fourier-on-Z/7 identity of mac-mini-0620).

  measS7(E) = sum_{a in (Z/7)^k} Khat(a) J(E,a),
     J(E,a)  = int_0^1 prod_e omega^{a_e c_e(x)} dx        (joint clock moment)
     Khat(a) = (1/7^k) sum_{c:{c_e}=Z/7} prod_e omega^{-a_e c_e}   (cover indicator coeff)
     c_e(x)  = floor(7 frac(e x)) in Z/7,   omega = e^{2 pi i/7}.

  DEFINE the per-clock MARGINAL moment  m(a) := int_0^1 omega^{a c_1(x)} dx
     where c_1(x)=floor(7x) on [0,1].  By a 1-line computation m(0)=1 and
     m(a)=0 for a != 0 mod 7 (the clock is exactly uniform on Z/7 -- the
     APEX-PRIME vanishing).  More generally for clock e, int_0^1 omega^{a c_e(x)}dx
     = m(a) too (each clock is individually uniform), so the *single-clock* moment
     is m(a) for EVERY e.

  CUT (single-particle / decorrelated) part:
     J_cut(E,a) := prod_e m(a_e)   = [a = 0]   (1 if a==0, else 0).
     measS7_cut(E) := sum_a Khat(a) J_cut(E,a) = Khat(0).
     Khat(0) = (1/7^k) #{c in (Z/7)^k : {c_e}=Z/7} = 7! S(k,7)/7^k = iid_k.
     => THE CUT PART IS EXACTLY THE IID SURJECTION PROBABILITY iid_k,
        which is INDEPENDENT OF E (depends only on k).  Trivially consec-extremal:
        it is constant over all k-shapes; consec attains it (no shape beats it).

  CYCLE (correlation) part:
     J_cyc(E,a) := J(E,a) - J_cut(E,a) = J(E,a) - [a=0].
     measS7_cyc(E) := measS7(E) - iid_k = corr(E).
     This is the GENUINELY E-DEPENDENT, relation-lattice-carried part:
        corr(E) = sum_{a != 0} Khat(a) J(E,a),  supported on a with 7 | sum a_e e.

  So measS7(E) = iid_k + corr(E)  EXACTLY, and the brief's wide-bound mechanism
  is precisely this CUT(=iid) + CYCLE(=corr) split.

WHAT THREAD E ESTABLISHES / TESTS HERE.
  (1) The split is an EXACT identity: measS7 = iid_k + corr, iid_k = Khat(0) = J_cut sum.
  (2) The CUT part is consec-extremal in the STRONGEST possible (degenerate) sense:
      it is constant = iid_k, the same for every shape.  So "the cut-part is consec-
      extremal" is TRUE but VACUOUS -- all extremality lives in the CYCLE part.
      THIS IS THE THM-559 ANALOGY *FAILING* in a precise, instructive way:
      for c3 the cut part (score variance) is where ALL the action is; for measS7
      the cut part is FLAT and the action is entirely in the cycle correction.
  (3) The 2-BODY (pair) layer of the cycle part:  corr(E) = sum_{t>=1} corr_t(E)
      where corr_t collects a of Hamming weight t.  We ask: is the 2-body layer
      corr_2 alone consec-extremal (the would-be Ising term)?  Test exactly.
  (4) Sign / magnitude budget: |corr(E)| and the per-weight |corr_t| vs the wide
      budget cap_k - iid_k, to see which layers a wide bound must control.

All exact: Fractions for measS7/iid; Fourier diagnostics cross-checked to the exact
geometric value.  stdlib only.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

w = cmath.exp(2j*math.pi/7)

# ---------------- exact geometric measS7 ----------------
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi<=lo: continue
        mid=(lo+hi)/2
        if len(set(int(((e*mid)%1)*7) for e in E))==7: tot+=hi-lo
    return tot

def stirling2(n,k): return sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1))//math.factorial(k)
def iid_k(k): return F(math.factorial(7)*stirling2(k,7), 7**k)

# ---------------- Fourier pieces ----------------
def J_exact(E, avec):
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=0j
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; val=1+0j
        for e,a in zip(E,avec):
            if a==0: continue
            val *= w**(a*int((e*xm % 1)*7))
        total += val*float(x1-x0)
    return total

def cM(M,a): return sum(w**(-a*j) for j in range(7) if j not in M)/7
_kc={}
def Khat(avec):
    key=tuple(sorted(avec))
    if key in _kc: return _kc[key]
    tot=0j
    for r in range(8):
        for Mset in itertools.combinations(range(7),r):
            p=1+0j
            for a in avec:
                p*=cM(Mset,a)
                if abs(p)<1e-18: break
            tot+=((-1)**r)*p
    _kc[key]=tot
    return tot

def corr_low_weight(E, wmax):
    """corr_t(E)=sum_{a: weight t, 7|sum a_e e} Khat(a) J(E,a), real, for t=1..wmax.
    Khat(a)=0 unless multiset {a_e} hits >=6 distinct nonzero values won't happen at
    low weight, so low-weight corr_t are exactly computable but Khat may still be tiny.
    We compute exactly. Heavy weights are NOT enumerated (left as a residual)."""
    k=len(E); byw=defaultdict(float)
    for r in range(1,wmax+1):
        for S in itertools.combinations(range(k), r):
            es=[E[i] for i in S]
            for vals in itertools.product(range(1,7), repeat=r):
                if sum(v*e for v,e in zip(vals,es))%7!=0: continue
                a=[0]*k
                for idx,v in zip(S,vals): a[idx]=v
                kv=Khat(a)
                if abs(kv)<1e-15: continue
                c=(kv*J_exact(E,a)).real
                if abs(c)<1e-15: continue
                byw[r]+=c
    return byw

# ---------------- single-clock marginal check (apex-prime vanishing) ----------------
def marginal_moment(e,a):
    """int_0^1 omega^{a*floor(7 e x)} dx, exact via breakpoints; should be m(a)=[a=0]."""
    bps=set([F(0),F(1)])
    for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=0j
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        total += w**(a*int((e*xm%1)*7))*float(x1-x0)
    return total

CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}

if __name__=="__main__":
    print("#"*86)
    print("# THREAD E: CUT/CYCLE split of measS7 = iid_k (cut) + corr(E) (cycle)")
    print("#"*86)

    print("\n[A] CUT part is the single-clock-product moment; verify J_cut(E,a)=[a=0].")
    print("    single-clock marginal m(e,a)=int omega^{a c_e}dx should be [a=0] (apex-prime).")
    ok=True
    for e in [1,2,3,5,8,13]:
        row=" ".join(f"a={a}:{abs(marginal_moment(e,a)):.3e}" for a in range(7))
        m0=marginal_moment(e,0)
        nz=max(abs(marginal_moment(e,a)) for a in range(1,7))
        if abs(m0-1)>1e-9 or nz>1e-9: ok=False
        print(f"    e={e:>2}: m(e,0)={m0.real:.6f}  max_{{a!=0}}|m|={nz:.2e}")
    print(f"    => each clock individually uniform on Z/7: m(a)=[a=0]  (verified={ok})")
    print("    => J_cut(E,a)=prod_e m(a_e)=[a=0]; measS7_cut = Khat(0) = iid_k (E-independent).")

    print("\n[B] EXACT IDENTITY measS7(E) = iid_k + corr(E)  [corr := measS7 - iid_k, FREE].")
    print("    + low-weight (2-body, 3-body) cycle layers via the Fourier-on-Z/7 split.")
    shapes={
        8:[("consec8",list(range(8))),("drop1_8",[0,1,2,3,4,5,6,8]),
           ("spread8",[0,1,2,3,4,5,6,9]),("twoblk8",[0,1,2,3,7,8,9,10]),
           ("iidlike8",[0,1,7,8,21,22,49,50])],
    }
    WMAX=3
    for k,shs in shapes.items():
        ik=iid_k(k)
        print(f"  k={k}: iid_k = {ik} = {float(ik):.6f}   cap={float(CAP[k]):.6f}  "
              f"wide-budget cap-iid={float(CAP[k]-ik):.6f}")
        for lab,E in shs:
            true=measS7(E)
            corr=true-ik                          # EXACT, free
            byw=corr_low_weight(E,WMAX)           # low-weight layers only
            low=sum(byw.values())
            prof=" ".join(f"c{t}:{byw.get(t,0.0):+.5f}" for t in sorted(byw))
            print(f"    {lab:<9} measS7={float(true):.6f} = iid {float(ik):.4f} + corr {float(corr):+.6f}")
            print(f"        cycle layers (t<= {WMAX}): [{prof}]  sum_low={low:+.6f}  "
                  f"high-residual(corr-low)={float(corr)-low:+.6f}")

    print("\n[C] LAYER EXTREMALITY: is consec the MAX of each cut/cycle layer separately?")
    print("    cut layer = iid_k (constant -> consec ties, never beaten: vacuously extremal).")
    print("    2-body cycle layer corr_2: does consec maximize it over all bounded shapes?")
    import itertools as it
    for k in [8]:
        cons=list(range(k)); cons_c2=corr_low_weight(cons,2).get(2,0.0)
        cons_full=float(measS7(cons))
        print(f"  k={k}: consec corr_2={cons_c2:+.6f}  consec measS7={cons_full:.6f}")
        worst_c2=(-9.9,None); worst_full=(-9.9,None); n_c2_beat=0; total=0
        for combo in it.combinations(range(1,14),k-1):
            E=[0]+list(combo); total+=1
            c2=corr_low_weight(E,2).get(2,0.0)
            full=float(measS7(E))
            if c2>worst_c2[0]: worst_c2=(c2,E)
            if full>worst_full[0]: worst_full=(full,E)
            if c2>cons_c2+1e-9: n_c2_beat+=1
        print(f"    over {total} bounded shapes (span<=13):")
        print(f"      MAX corr_2 = {worst_c2[0]:+.6f} at {worst_c2[1]}  (consec corr_2={cons_c2:+.6f}; "
              f"#shapes beating consec in corr_2 = {n_c2_beat})")
        print(f"      MAX measS7 = {worst_full[0]:.6f} at {worst_full[1]}  consec is max: "
              f"{worst_full[1]==cons}")
        print(f"    VERDICT: 2-body layer consec-extremal? {n_c2_beat==0}")
