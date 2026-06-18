"""
DECISIVE: explicitly exhibit ALL-MULT7-LARGE covering primitive S3 sets for which
the Window-Collapse witness construction FAILS (mult-of-7 subsystem teeth cover the
ENTIRE window (0,1/(2V*)]), so NO global witness tau=k/7+u/7 exists.

Show: (a) exact window coverage = 100% (no gap); (b) the set is genuinely in scope;
(c) M(S)>=1/14 nevertheless (so closure CONCLUSION may survive, but the stated proof
mechanism / exhaustive-sweep claim does not). Find the MINIMAL V* such a config exists.
kind-pasteur S4-wf.
"""
from fractions import Fraction as F
from math import gcd
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def window_gaps(W,Vstar):
    """Return list of safe gaps in (0,R], R=1/(2V*), for subsystem W."""
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            a=max(F(0),lo); b=min(R,c+h/F(w))
            if a<=b: forb.append((a,b))
            j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0); gaps=[]
    for a,b in merged:
        if a>prev: gaps.append((prev,a))
        prev=max(prev,b)
    if prev<R: gaps.append((prev,R))
    return [(lo,hi) for lo,hi in gaps if hi>lo and hi>0], R

print("="*70)
print("CLEANEST EXPLICIT COUNTEREXAMPLE to witness-existence (mechanism failure)")
print("="*70)
S=[1, 6, 8, 10, 15, 16, 18, 22, 24, 26, 28, 378, 581]
m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]; Vstar=max(nm7)
W=[v//7 for v in m7]
gaps,R=window_gaps(W,Vstar)
print("S =",S)
print("  covering:",is_cov(S)," primitive:",is_prim(S)," |S|=",len(S))
Vmin=min(S);Vmax=max(S);k=sum(1 for v in S if v>13)
print("  S3? (k>=2,Vmax>=13Vmin):",k>=2 and Vmax>=13*Vmin,"  (k=",k,"Vmin=",Vmin,"Vmax=",Vmax,")")
print("  mult-of-7:",m7," V*=",Vstar," ALL-MULT7-LARGE (all>V*):",all(x>Vstar for x in m7))
print("  subsystem W=v/7:",W," window R=1/(2V*)=",R)
print("  SAFE GAPS in window:",gaps,"  => witness exists?",bool(gaps))
print("  Mval(S) =",Mval(S),"~",float(Mval(S))," >=1/14?",Mval(S)>=F(1,14))
print("  CONCLUSION: ALL-MULT7-LARGE & in scope, but NO tau=k/7+u/7 witness; M>=1/14 via other tau.")

print()
print("="*70)
print("MINIMAL-V* search: smallest V* with an ALL-MULT7-LARGE no-witness subsystem")
print("="*70)
# For each V*, look for W (size<=4), all w>V*/7, that fully cover (0,R].
# Strategy: wmin gives wide first tooth; add larger w to fill remaining gaps.
import itertools
def covers_window(W,Vstar):
    g,_=window_gaps(W,Vstar); return len(g)==0
found=[]
for Vstar in range(14,140):
    wmin=Vstar//7+1
    R=F(1,2*Vstar)
    # w candidates: all w with 7w>V*, i.e. w>=wmin, up to a bound where teeth still
    # reach the window: need some j with j/w in (0,R], i.e. w>=1/R=2V* gives center
    # at 1/w<=R. Actually any w has its j=0 tooth at 0; far teeth need w<= ... we just
    # allow w up to 2*Vstar+ a margin (teeth of period 1/w; to land inside (0,R] need
    # 1/w <= R => w>=2V*; but w slightly below also contributes via j=0 right half).
    wpool=list(range(wmin, 2*Vstar+2))
    # too large to brute force C(.,4); use greedy: pick wmin, then search to fill gap.
    base=[wmin]
    g,_=window_gaps(base,Vstar)
    # try to cover the (single) remaining gap with up to 3 more w's
    ok=False; chosen=None
    if not g:
        ok=True; chosen=tuple(base)
    else:
        # greedy fill: repeatedly find a w whose tooth covers the left end of current first gap
        def try_fill(cur, depth):
            gg,_=window_gaps(cur,Vstar)
            if not gg: return tuple(sorted(set(cur)))
            if len(set(cur))>=4: return None
            lo,hi=gg[0]
            # need a w (>=wmin) with a tooth-center j/w in [lo-1/(14w), ...] covering lo.
            # i.e. exists j>=1 with |j/w - lo| <= 1/(14 w)  => |j - w lo|<=1/14.
            for w in wpool:
                if w in cur: continue
                # nearest j to w*lo
                wl=F(w)*lo
                jj=round(wl)
                for j in (int(wl), int(wl)+1, jj):
                    if j<1: continue
                    c=F(j,w)
                    if abs(c-lo) <= F(1,14)/w + F(0):  # tooth covers lo
                        res=try_fill(cur+[w], depth+1)
                        if res: return res
            return None
        res=try_fill(base,1)
        if res: ok=True; chosen=res
    if ok and chosen and len(set(chosen))>=2 and all(7*w>Vstar for w in chosen):
        found.append((Vstar,chosen))
        if len(found)<=6:
            g2,_=window_gaps(chosen,Vstar)
            print("  V*=",Vstar," W=",chosen," 7W=",tuple(7*w for w in chosen),
                  " gaps=",g2," (fully blocked:",len(g2)==0,")")
        if len(found)>=6: break
if found:
    print("  MINIMAL V* with no-witness subsystem:",found[0][0]," W=",found[0][1])
else:
    print("  none found in V*<140 by greedy (try the bulk-hunt examples instead)")
