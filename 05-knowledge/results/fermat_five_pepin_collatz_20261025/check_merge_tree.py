# Exact 2-adic merge tree for the Fermat line pairs with L = 1 mod 4.
# State at Terras time 2L+1: v = p = 3^(L+1)+2, u = p-1, chain (k,e) = (0,-1).  The next N parity bits of v depend only on
# p mod 2^N, hence on r = (L+1) mod 2^(N-2).  Absorption (0,0) at post-prefix time n' <= N is decided exactly per class.
from fractions import Fraction as Fr
import json, sys
def step(k,e,beta):
    f=e*Fr(3)**max(0,-k); assert f.denominator==1
    sigma=f.numerator%2
    if sigma==0 and beta==0: return (k, e/2)
    if sigma==0 and beta==1: return (k, (3*e+1-Fr(3)**k)/2)
    if sigma==1 and beta==0: return (k+1, (3*e+1)/2)
    return (k-1, (e-Fr(3)**(k-1))/2)
N=int(sys.argv[1]) if len(sys.argv)>1 else 20
M=2**N; cls_mod=2**(N-2)
absorb={}   # r -> (n', odd steps)
for r in range(0, cls_mod):              # all classes; L = r-1 mod 2^(N-2)
    p=(pow(3,r,M)+2)%M
    v=p; k,e=0,Fr(-1); odd=0; t=None
    for n in range(1,N+1):
        beta=v%2
        if beta: odd+=1
        k,e=step(k,e,beta)
        v=(3*v+1)//2 % M if beta else v//2      # residue mod 2^N suffices for parity of the next N-n steps
        if (k,e)==(0,0): t=(n,odd); break
    absorb[r]=t
dist={}
for r,t in absorb.items():
    key=t[0] if t else None
    dist[key]=dist.get(key,0)+1
total=len(absorb)
# controls by type of L = r-1: even L (r odd) must absorb at n'=3 with 1 odd step; L = 3 mod 4 (r = 0 mod 4) at n'=4 with 2 odd steps
ctrl_even=all(absorb[r]==(3,1) for r in range(1,cls_mod,2))
ctrl_3=all(absorb[r]==(4,2) for r in range(0,cls_mod,4))
ones={r:t for r,t in absorb.items() if r%4==2}
dist={}
for r,t in ones.items():
    key=t[0] if t else None
    dist[key]=dist.get(key,0)+1
total=len(ones)
absorb=ones
# minimal classes: for each absorption time n', the class is determined mod 2^(n'-2); list residues of L (= r-1) mod 2^(n'-2)
minimal={}
for r,t in absorb.items():
    if t:
        n,odd=t; mod=2**(n-2)
        minimal.setdefault(n,set()).add(((r-1)%mod, odd))
out={'N':N,'control_even_L_absorb_at_3_with_1_odd':ctrl_even,'control_L3mod4_absorb_at_4_with_2_odd':ctrl_3,'classes_L_1mod4':total,'absorbed_by_time':{str(k):v for k,v in sorted(dist.items(), key=lambda kv: (kv[0] is None, kv[0]))},
     'fraction_unabsorbed_at_N':round(dist.get(None,0)/total,4),
     'minimal_classes_L_mod_2^(n-2)_with_odd_steps':{str(n):sorted(v) for n,v in sorted(minimal.items()) if n<=14}}
# cross-check with direct data for L <= 1500 (from oddL_depth.py logic)
def T(x): return x>>1 if x%2==0 else (3*x+1)>>1
bad=0; checked=0
for L in range(1,1501,4):
    x=2**(2*L+1)+1; y=2**(2*L+2)+1
    u=T(y); v=x; n=0; direct=None
    while v>2:
        if u==v: direct=n; break
        v=T(v); u=T(u); n+=1
    r=(L+1)%cls_mod; pred=absorb[r]
    post = None if direct is None else direct-(2*L+1)
    checked+=1
    if pred is None:
        if post is not None and post<=N: bad+=1
    else:
        if post!=pred[0]: bad+=1
out['direct_crosscheck']={'checked':checked,'mismatches':bad}
print(json.dumps(out,indent=0))
