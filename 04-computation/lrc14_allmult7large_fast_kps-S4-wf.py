"""
FAST decisive probe: does the window-collapse witness ALWAYS exist for covering
ALL-MULT7-LARGE S3 sets?  And independently: does ANY covering ALL-MULT7-LARGE
set break LRC (M<1/14)?  Capped reduced-w-sets r<=2, plus a direct random/struct
search over actual sets for an LRC break.

Key structural fact under test:
 In ALL-MULT7-LARGE, the multiples of 7 in S are the LARGE runners 7 w_i with
 7 w_i > V*.  Covering forces a multiple of 14 => some 7 w_i even => some w_i even,
 with 7 w_i > V* => w_i >= ceil((V*+1)/7).  The window is u in (0, 1/(2 V*)].
 The {w_i}-safe set near 0: nearest tooth of runner w_i is at u=1/(2 w_i) (its
 half-period) -- NO, teeth (danger) are centered at j/w_i width 1/(7 w_i); the
 safe region around 0 extends to u = (1/2 - 1/14)/w_min_effect... we just compute.
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
def safe_components(A,h=F(1,14)):
    iv=[]
    for u in A:
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

def window_witness_exists(wlist, Vstar):
    bound=F(1,2*Vstar)
    sc=safe_components(sorted(set(wlist)))
    for (a,b) in sc:
        lo=max(a,F(0)); hi=min(b,bound)
        if hi>lo and hi>F(0): return True
        lo2=max(a,F(1)-bound); hi2=min(b,F(1))
        if hi2>lo2: return True
    return False

R=[]
def out(s): R.append(s); print(s, flush=True)

out("FAST ALL-MULT7-LARGE probe (r<=2 reduced w-sets, V*<=200)")
out("="*72)
blocked=[]
checked=0
for Vstar in range(14,201):
    wmin=Vstar//7+1
    wmax=(6*Vstar)//7
    wpool=list(range(wmin,wmax+1))
    if not wpool: continue
    # r=1: single even mult-of-7 large runner
    for w in wpool:
        if w%2==0:
            checked+=1
            if not window_witness_exists([w],Vstar): blocked.append((Vstar,(w,)))
    # r=2: two large mult-of-7, at least one even
    for combo in itertools.combinations(wpool,2):
        if not any(w%2==0 for w in combo): continue
        checked+=1
        if not window_witness_exists(combo,Vstar): blocked.append((Vstar,combo))
out(f"checked {checked} reduced w-sets (r<=2); blocked windows: {len(blocked)}")
if blocked:
    out("  sample blocked: "+str(blocked[:15]))
else:
    out("  => NO blocked window for any covering-forced (even-containing) w-set r<=2.")
    out("     The even w_i is the binding one; its first danger tooth is at u=1/(2w),")
    out("     leaving a safe gap (0, something) that always overlaps (0,1/(2V*)] since")
    out("     V* >= 7 w/... let us print the tightest margin.")

# tightest margin diagnostic for r=1 even-w case: safe gap right of 0 ends where?
out("")
out("r=1 even-w margin: for 7w>V* (so w>=ceil((V*+1)/7)), is (0,1/(2V*)] inside safe(0)?")
worst=None
for Vstar in range(14,201):
    wmin=Vstar//7+1
    for w in range(wmin, (6*Vstar)//7+1):
        if w%2: continue
        if 7*w<=Vstar: continue
        sc=safe_components([w])
        # find safe component containing points just >0
        # safe set of single runner w: gaps (danger) around j/w of half-width 1/(14w).
        # around 0 the danger is [-1/(14w),1/(14w)]; safe starts at 1/(14w).
        # first safe interval to the right of 0 is (1/(14w), 1/w - 1/(14w)).
        a=F(1,14*w); b=F(1,w)-F(1,14*w)
        bound=F(1,2*Vstar)
        # need a point in (0,bound] that is safe: safe right of 0 starts at a=1/(14w).
        # is a < bound?  i.e. 1/(14w) < 1/(2 V*) <=> 2 V* < 14 w <=> V* < 7 w. TRUE (7w>V*).
        ok = a < bound
        marg = bound - a
        if worst is None or marg<worst[0]:
            worst=(marg, Vstar, w, a, bound, ok)
out(f"  tightest (smallest) margin bound-1/(14w): {worst[0]} at V*={worst[1]} w={worst[2]}")
out(f"    safe-start 1/(14w)={worst[3]}, window-edge 1/(2V*)={worst[4]}, witness-room>0: {worst[5]}")
out("  PROOF NOTE: for any even w with 7w>V*, 1/(14w) < 1/(2V*) <=> V* < 7w (TRUE),")
out("  so the point u just above 1/(14w) is in-window AND w-safe.  But OTHER w_i may")
out("  re-block; r<=2 sweep shows they never do for V*<=200.  A second tooth of the")
out("  even w_i sits at 1/w > 1/(2V*) (since w<=6V*/7<V* ... check): need 1/w>1/(2V*)")
out("  i.e. w<2V*, true since w<=6V*/7. So only the central tooth is in-window.")

# Direct LRC-break search over actual covering ALL-MULT7-LARGE sets (structured)
out("")
out("DIRECT: search actual covering ALL-MULT7-LARGE S3 sets for M<1/14")
out("="*72)
random.seed(7)
breaks=0; tested=0; minM=(F(1),None)
for trial in range(60000):
    Vstar=random.randint(14,60)
    # small runners: subset of non-mult-of-7 in 1..V*, must include V*
    smallpool=[v for v in range(1,Vstar+1) if v%7!=0]
    if Vstar not in smallpool: continue
    # large mult-of-7: pick 2-3 multiples of 7 above V*, include an even one
    bigpool=[7*w for w in range(Vstar//7+1, Vstar//7+1+10)]
    nbig=random.choice([2,3])
    big=random.sample(bigpool, min(nbig,len(bigpool)))
    if not any(b%14==0 for b in big):
        # force an even mult of 7 for q=14
        evenbig=[b for b in bigpool if b%14==0]
        if not evenbig: continue
        big[0]=random.choice(evenbig)
    big=sorted(set(big))
    need=13-len(big)
    if need<1 or need-1>len(smallpool)-1: continue
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
    m=Mval(S)
    tested+=1
    if m<minM[0]: minM=(m,tuple(S))
    if m<F(1,14): breaks+=1; out(f"  !!! BREAK S={S} M={m}")
out(f"tested {tested} covering ALL-MULT7-LARGE S3 sets; LRC breaks: {breaks}")
out(f"min M = {minM[0]} = {float(minM[0]):.5f} at S={list(minM[1]) if minM[1] else None}")
out(f"min M >= 1/14 ? {minM[0]>=F(1,14)}")

with open(r"C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\lrc14_allmult7large_fast_kps-S4-wf.out","w") as f:
    f.write("\n".join(R))
out("\n[written]")
