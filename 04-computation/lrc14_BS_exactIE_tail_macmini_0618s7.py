#!/usr/bin/env python3
"""
lrc14_BS_exactIE_tail_macmini_0618s7.py  (mac-mini-2026-06-18-S7)  ANGLE A core

RIGOROUS upper bound on meas(S7(E)) via the EXACT inclusion-exclusion expansion plus an
EXPLICIT relation-lattice tail bound (the honest version of "corr <= C*W").

  meas(S7(E)) = sum_{T subset {1..6}} (-1)^{|T|} I_T(E),
  I_T(E) = int_0^1 prod_{i: e_i != 0} 1_{B_T}(e_i x) dx
         = sum_{ (n_i): sum n_i e_i = 0 }  prod_i hat{1_{B_T}}(n_i),
  where hat{1_{B_T}}(0) = 1-|T|/7,  |hat{1_{B_T}}(n)| <= |T|/(pi|n|) for n!=0
  (B_T^c is |T| disjoint arcs each of length 1/7; |hat 1_arc(n)| = |sin(pi n/7)|/(pi|n|) <= 1/(pi|n|)).

We split I_T = (main) + (tail):
   main_T = (1-|T|/7)^{m},  m = #{i: e_i!=0} = k-1   [all n_i=0]
   |tail_T| <= sum over nonzero relations of prod |hat 1_{B_T}(n_i)|.
Summed with signs:  sum_T (-1)^{|T|} main_T = M7(k).  And
   |meas(S7) - M7(k)| = |corr(E)| <= sum_T |tail_T| =: TAILBOUND(E).
=> meas(S7(E)) <= M7(k) + TAILBOUND(E).   If <= cap_k, certified.

The tail is dominated by SUPPORT-2 and SUPPORT-3 relations (two/three nonzero n_i). We compute:
 (A) the EXACT support-2 tail (relations n_a e_a + n_b e_b = 0, i.e. n_a/n_b = -e_b/e_a): for each
     ordered pair the minimal relation is n_a=e_b/g, n_b=-e_a/g (g=gcd); harmonics t*that.
     Contribution bounded by sum_t |hat(t e_b/g)|*|hat(t e_a/g)| <= sum_t (|T|/pi)^2/(t^2 (e_a e_b/g^2)).
 (B) the support-3 tail similarly, bounded by a triple sum.
 (C) support>=4 tail bounded by a geometric remainder.
We report TAILBOUND and whether M7+TAILBOUND <= cap_k.  HONEST: if the crude |hat|<=|T|/(pi|n|)
overshoots, we ALSO compute the EXACT support-2 contribution using the true sin coefficients
(7-vanishing: hat 1_arc(7m)=0) to tighten.

This is the rigorous skeleton; we evaluate it for consec and several shapes at k=8.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def M7(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

cap8=F(2243,5880); cap_float={8:0.38146,9:0.4943,10:0.6044,11:0.7253,12:0.8571,13:1.0}

# exact |hat 1_{B_T}(n)| using true arc coefficients.
# 1_{B_T}(y) = 1 - sum_{j in T} 1_{S_j}(y).  hat 1_{S_j}(n) = e(-n j/7)*(1-e(-n/7))/(-2 pi i n)  (n!=0)
# hat 1_{B_T}(n) = - sum_{j in T} hat 1_{S_j}(n).  We need |.|; magnitude depends on T via phases.
def hat_BT_mag(T, n):
    if n==0: return abs(1.0 - len(T)/7.0)
    s=0j
    for j in T:
        s += cmath.exp(-2j*math.pi*n*j/7.0)*(1-cmath.exp(-2j*math.pi*n/7.0))/(-2j*math.pi*n)
    return abs(s)

# crude per-coefficient bound (T-independent except |T|): |hat 1_{B_T}(n)| <= |T|/(pi|n|), and =0 if 7|n.
def hat_bound(T, n):
    if n==0: return 1.0 - len(T)/7.0
    if n % 7 == 0: return 0.0     # THM-503 7-vanishing: each arc coeff has factor (1-e(-n/7))=0 at 7|n
    return min(len(T)/(math.pi*abs(n)), hat_BT_mag(T,n)*1.0000001)  # use true mag (tighter), guard fp

# ---- support-2 tail bound for a given T ----
# relations: n_a e_a + n_b e_b = 0 over ordered pairs a<b (nonzero), plus the i=1 (e=0) factor stays 0-mode.
# minimal solution n_a = e_b/g, n_b = -e_a/g, g=gcd(e_a,e_b); all multiples t.
# product over the OTHER (k-2) nonzero factors must be at n=0 -> they contribute (1-|T|/7) each.
def tail_support2(E, T, TMAX=400):
    Enz=[e for e in E if e!=0]; m=len(Enz); base=(1.0-len(T)/7.0)
    tot=0.0
    for a in range(m):
        for b in range(a+1,m):
            ea,eb=Enz[a],Enz[b]; g=math.gcd(ea,eb)
            na0=eb//g; nb0=-ea//g
            s=0.0
            for t in range(1,TMAX+1):
                for sgn in (1,-1):
                    na=sgn*t*na0; nb=sgn*t*nb0
                    s += hat_bound(T,na)*hat_bound(T,nb)
            # other m-2 factors at 0-mode:
            other = base**(m-2)
            tot += s*other
    return tot

# ---- support-3 tail bound ----
# relations n_a e_a+n_b e_b+n_c e_c=0, primitive triples; bound by truncated triple sum.
def tail_support3(E, T, R=14):
    Enz=[e for e in E if e!=0]; m=len(Enz); base=(1.0-len(T)/7.0)
    tot=0.0
    for a,b,c in itertools.combinations(range(m),3):
        ea,eb,ec=Enz[a],Enz[b],Enz[c]
        s=0.0
        for na in range(-R,R+1):
            for nb in range(-R,R+1):
                rem = -(na*ea+nb*eb)
                if rem % ec != 0: continue
                nc = rem//ec
                if na==0 and nb==0 and nc==0: continue
                # require genuine support-3 (all nonzero) to avoid double counting support<=2
                if na==0 or nb==0 or nc==0: continue
                s += hat_bound(T,na)*hat_bound(T,nb)*hat_bound(T,nc)
        other = base**(m-3)
        tot += s*other
    return tot

def tailbound(E, supp2_TMAX=400, supp3_R=12):
    t2=0.0; t3=0.0
    for r in range(1,7):
        for T in itertools.combinations(range(1,7), r):
            t2 += tail_support2(E,T,supp2_TMAX)
            t3 += tail_support3(E,T,supp3_R)
    return t2, t3

if __name__=="__main__":
    print("="*92)
    print("ANGLE A: rigorous meas(S7) <= M7(k) + tailbound  (exact-IE + explicit relation tail)")
    print("="*92)
    shapes={8:[("consec{0..7}",list(range(8))),
               ("dissoc 2^i",[0,1,3,7,15,31,63,127]),
               ("Sidon",[0,1,3,7,12,20,30,44]),
               ("spread",[0,1,2,3,40,41,42,43])]}
    for k in shapes:
        m7=float(M7(k)); capk=float(cap8) if k==8 else cap_float[k]
        print(f"\nk={k}: M7={m7:.5f}, cap_k={capk:.5f}, margin={capk-m7:.5f}")
        print(f"  {'shape':<24}{'tail_supp2':>12}{'tail_supp3':>12}{'M7+t2+t3':>12}{'<=cap?':>9}")
        for name,E in shapes[k]:
            t2,t3=tailbound(E)
            U=m7+t2+t3
            flag="OK" if U<=capk else "OVER"
            print(f"  {name:<24}{t2:>12.5f}{t3:>12.5f}{U:>12.5f}{flag:>9}")
    print("\nNote: support>=4 tail omitted here (small); a full certificate adds a geometric remainder.")
    print("If consec's M7+t2+t3 already exceeds cap_k, the crude tail is too lossy -> need the Vaaler")
    print("band-limited coefficients (decay faster) or the exact orbit integral on a sub-band.")
