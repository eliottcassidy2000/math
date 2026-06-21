#!/usr/bin/env python3
"""
lrc_corr_support_ladder_macmini.py
THREAD A deliverable: decompose corr(E) by the SUPPORT SIZE of the offset relation (Lambda(E)),
to test whether the SUPPORT-3 additive-triangle (Schur-triple) layer is the LEADING term, and
whether CROSS-BLOCK relations vanish as blocks separate.

This is the LRC analog of the OCF cycle-space layer decomposition (THM-560 degree ladder:
deg_b(c_{2k+1})=2k; the support/length of a relation = the analog of cycle length).

corr(E) = sum_{0 != n in Lambda(E)} W(n),   Lambda(E)={n in Z^k: sum n_e e =0},
   W(n) = sum_{T subset Z/7} (-1)^|T| prod_e a_{n_e}(T).

We enumerate relations n with bounded entries |n_e|<=B, GROUP THEM BY |support(n)| = #nonzero,
and sum W(n) within each support layer.  Compare partial sums to the exact corr.

KEY TESTS:
 (1) corr(E) ~ sum over support layers; does it converge with the truncation B and reproduce
     exact corr?  (validates lattice picture at k=7,8 where corr != 0)
 (2) Is the SUPPORT-3 layer (Schur triples a+b=c) the dominant CROSS-BLOCK contribution?
 (3) For two-block E, partition relations into WITHIN-BLOCK (support inside one block) and
     CROSS-BLOCK (support touches >=2 blocks). Does the cross-block contribution -> 0 as M->inf,
     governed by the support-3 (Schur) cross-block relations?  --> the multi-block lever.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): bps.add(F(mm,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            if e==0: continue
            secs.add(int(((e*xm)%1)*7))
        if 0 in E: secs.add(0)
        if len(secs)==7: total+=x1-x0
    return total

from math import comb
def iid_k(k):
    return sum((-1)**j*comb(7,j)*F(7-j,7)**k for j in range(8))

# Fourier coeff a_m(T)
_acache={}
def a_m(m, T):
    key=(m,T)
    if key in _acache: return _acache[key]
    Tset=set(T)
    if m==0:
        v=complex((7-len(Tset))/7.0,0.0)
    else:
        s=0j
        for j in range(7):
            if j in Tset: continue
            hi=cmath.exp(-2j*math.pi*m*(j+1)/7); lo=cmath.exp(-2j*math.pi*m*j/7)
            s+=(hi-lo)/(-2j*math.pi*m)
        v=s
    _acache[key]=v
    return v

# Precompute T-subset list with signs
_allT=[]
for r in range(0,8):
    for T in itertools.combinations(range(7), r):
        _allT.append((T,(-1)**r))

def W(n):
    """relation weight: sum_T (-1)^|T| prod a_{n_e}(T)."""
    w=0j
    for T,sgn in _allT:
        prod=1.0+0j
        for ne in n:
            if ne==0:
                prod*=a_m(0,T)
            else:
                prod*=a_m(ne,T)
            if abs(prod)<1e-300: break
        w+=sgn*prod
    return w

def corr_by_support(E, B, blocks=None):
    """Enumerate relations |n_e|<=B; group W(n) by support size. Optionally split within/cross block.
       blocks: list of index-sets partitioning range(k); a relation is CROSS if support meets >=2 blocks."""
    E=sorted(set(E)); k=len(E)
    layer=defaultdict(complex)        # support-size -> sum W
    within=defaultdict(complex)
    cross=defaultdict(complex)
    block_of=None
    if blocks is not None:
        block_of=[None]*k
        for bi,bs in enumerate(blocks):
            for i in bs: block_of[i]=bi
    for n in itertools.product(range(-B,B+1), repeat=k):
        if all(v==0 for v in n): continue
        if sum(n[i]*E[i] for i in range(k))!=0: continue
        supp=[i for i in range(k) if n[i]!=0]
        s=len(supp)
        wv=W(n)
        layer[s]+=wv
        if blocks is not None:
            bset={block_of[i] for i in supp}
            if len(bset)>=2: cross[s]+=wv
            else: within[s]+=wv
    return layer, within, cross

if __name__=="__main__":
    print("="*92)
    print("(1) corr by SUPPORT layer, validate vs exact (k=7,8)")
    print("="*92)
    tests=[("consec7",list(range(7)),3),
           ("consec8",list(range(8)),2),
           ("[0,1,2,3,4,5,7]",[0,1,2,3,4,5,7],3)]
    for name,E,B in tests:
        exact=float(measS7(E)-iid_k(len(set(E))))
        layer,_,_=corr_by_support(E,B)
        tot=sum(layer.values()).real
        print(f"  {name:<18} B={B}  exact_corr={exact:+.5f}  lattice(trunc)={tot:+.5f}  diff={tot-exact:+.5f}")
        for s in sorted(layer):
            print(f"      support {s}: sum W = {layer[s].real:+.6f}  (cumulative {sum(layer[t].real for t in sorted(layer) if t<=s):+.6f})")
        print()

    print("="*92)
    print("(2/3) TWO-BLOCK: within vs cross-block contribution by support, as M grows (k=7)")
    print("="*92)
    half=3
    base=list(range(half))            # [0,1,2]
    B=3
    for M in (10, 30, 100, 300):
        E=base+[M+i for i in range(4)]  # k=7: block0={0,1,2}, block1 far
        E=sorted(set(E))
        if len(E)!=7: continue
        blocks=[[i for i,e in enumerate(E) if e<M],[i for i,e in enumerate(E) if e>=M]]
        exact=float(measS7(E)-iid_k(7))
        layer,within,cross=corr_by_support(E,B,blocks=blocks)
        tot=sum(layer.values()).real; wtot=sum(within.values()).real; ctot=sum(cross.values()).real
        print(f"  M={M:>4}  E={E}")
        print(f"      exact_corr={exact:+.5f}  lattice={tot:+.5f}  within={wtot:+.5f}  cross={ctot:+.5f}")
        for s in sorted(set(within)|set(cross)):
            print(f"        support {s}: within={within[s].real:+.6f}  cross={cross[s].real:+.6f}")
        print()
    print("  LEVER TEST: if 'cross' total -> 0 as M->inf, dominated by the SUPPORT-3 (Schur) cross")
    print("  layer, then |corr - within| <= (Schur-cross bound) closes the multi-block case.")
