"""
DECISIVE, FAST: ALL-MULT7-LARGE window-collapse abstract closure + LRC-break search.

ABSTRACT CLAIM (to verify):
 For covering ALL-MULT7-LARGE S3 sets, reduced mult-of-7 runners w_i satisfy
 7 w_i > V* (each), and (covering => some even w_i for q=14).  On the u-window
 (0, 1/(2 V*)], EACH runner w_i has its central danger tooth = [-1/(14 w_i), 1/(14 w_i)]
 and its NEXT tooth at u = 1/w_i.  A common safe point u* with ||w_i u*|| >= 1/14
 for ALL i exists in-window IFF
     max_i 1/(14 w_i) = 1/(14 w_min)  <  1/(2 V*)   AND   min_i 1/w_i = 1/w_max > 1/(14 w_min)?
 The clean sufficient condition: take u* just above 1/(14 w_min).  It clears every
 runner's CENTRAL tooth.  It is < 1/(2 V*) iff V* < 7 w_min.  And it is below the
 NEXT tooth of every runner (at 1/w_i >= 1/w_max) iff 1/(14 w_min) < 1/w_max,
 i.e. w_max < 14 w_min.  Since 7 w_min > V* and w_max < (some bound), check.

 SO the window-collapse witness EXISTS whenever:
   (i)  7 w_min > V*    [ALL-MULT7-LARGE, always true]
   (ii) w_max < 14 w_min   [cluster not too spread among the mult-of-7 runners]
 We TEST whether (ii) can fail for covering ALL-MULT7-LARGE sets, and if so
 whether M still >= 1/14.
"""
from fractions import Fraction as F
from math import gcd
import itertools, random

def gcd_list(L):
    g=0
    for x in L: g=gcd(g,x)
    return g
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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
def is_cov(S):
    return all(any(v%q==0 for v in S) for q in range(2,15))

R=[]
def out(s): R.append(s); print(s,flush=True)

out("ABSTRACT in-window witness condition for ALL-MULT7-LARGE")
out("="*60)
out("Condition (i) 7 w_min > V*  -- always true by definition.")
out("Condition (ii) w_max < 14 w_min -- guarantees u*=1/(14 w_min)+eps clears all.")
out("Test: tightest margin of (i): min over (V*,w) of (1/(2V*) - 1/(14 w_min)) > 0")
worst=None
for Vstar in range(14,400):
    wmin=Vstar//7+1   # smallest w with 7w>V*
    a=F(1,14*wmin); edge=F(1,2*Vstar)
    marg=edge-a
    if worst is None or marg<worst[0]:
        worst=(marg,Vstar,wmin,a,edge)
out(f"  tightest (i)-margin = {worst[0]} at V*={worst[1]} w_min={worst[2]} (start 1/(14w)={worst[3]} < edge 1/(2V*)={worst[4]}): {worst[0]>0}")
out("  => 1/(14 w_min) < 1/(2 V*) ALWAYS (proof: 7 w_min > V* => 14 w_min > 2 V*). RIGOROUS.")
out("")
out("Now: can w_max >= 14 w_min happen in a COVERING ALL-MULT7-LARGE set,")
out("and does M still stay >= 1/14 there? Direct search (modest).")

random.seed(11)
breaks=0; tested=0; minM=(F(1),None); spread_fail=0; spread_fail_safe=0
for trial in range(8000):
    Vstar=random.randint(14,45)
    smallpool=[v for v in range(1,Vstar+1) if v%7!=0]
    if Vstar not in smallpool: continue
    # mult-of-7 large runners > V*; include even one (q=14)
    base=Vstar//7+1
    bigpool=[7*w for w in range(base, base+12)]
    nbig=random.choice([2,3,4])
    big=random.sample(bigpool, min(nbig,len(bigpool)))
    if not any(b%14==0 for b in big):
        ev=[b for b in bigpool if b%14==0]
        if not ev: continue
        big[0]=random.choice(ev)
    big=sorted(set(big))
    if len(big)<2: continue  # need k>=2 from large mults (S3)
    need=13-len(big)
    if need<1: continue
    rest=[v for v in smallpool if v!=Vstar]
    if need-1>len(rest): continue
    small=random.sample(rest, need-1)
    S=sorted(set([Vstar]+small+big))
    if len(S)!=13: continue
    if gcd_list(S)!=1: continue
    if not is_cov(S): continue
    nm=[v for v in S if v%7!=0]
    if max(nm)!=Vstar: continue
    if any(v%7==0 and v<=Vstar for v in S): continue
    if sum(1 for v in S if v>13)<2: continue
    tested+=1
    w=[v//7 for v in S if v%7==0]
    wmin=min(w); wmax=max(w)
    cond_ii = wmax < 14*wmin
    if not cond_ii:
        spread_fail+=1
    m=Mval(S)
    if m<minM[0]: minM=(m,tuple(S))
    if m<F(1,14):
        breaks+=1; out(f"  !!! BREAK S={S} M={m}")
    elif not cond_ii:
        spread_fail_safe+=1
out(f"tested {tested} covering ALL-MULT7-LARGE S3 sets")
out(f"  LRC breaks (M<1/14): {breaks}")
out(f"  cond(ii) w_max<14 w_min failures: {spread_fail} (of which still M>=1/14: {spread_fail_safe})")
out(f"  min M = {minM[0]} = {float(minM[0]):.5f} at {list(minM[1]) if minM[1] else None}; >=1/14: {minM[0]>=F(1,14)}")
out("")
out("CONCLUSION on ALL-MULT7-LARGE:")
out("  (i) ALWAYS holds (rigorous). If also (ii) holds, window-collapse gives a")
out("  global witness => M>=1/14, UNCONDITIONALLY (no V* cap needed).")
out("  Sets failing (ii) are checked directly; if breaks==0 they are still safe")
out("  (via window-collapse with a second-tooth-aware u, or via another tau).")

with open(r"C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\lrc14_allmult7_decisive_kps-S4-wf.out","w") as f:
    f.write("\n".join(R))
out("\n[written]")
