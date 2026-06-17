"""
DISPROVE-side: pin the LOGIC of a genuine LRC(14) counterexample and search for it.

LRC(14) (in our normalization, 13 speeds, gap 1/14) asks: for every 13-set S of
distinct positive integers, does there exist tau with min_v ||v tau|| >= 1/14
(equivalently MM(S) >= 1/14)? A COUNTEREXAMPLE is S with MM(S) < 1/14.

Relation to L (lonely measure, OPEN inequality ||v tau||>1/14):
  - L(S) > 0  => exists tau with all ||v tau|| > 1/14 => MM(S) > 1/14 > can't be CE.
  - L(S) = 0  => the open lonely set has measure 0. Two sub-cases:
       (a) MM(S) = 1/14 : LRC holds with EQUALITY (closed arcs cover, a single
           witness tau* attains 1/14). NOT a counterexample.
       (b) MM(S) < 1/14 : closed danger arcs (radius 1/14) cover [0,1] with margin.
           This IS a counterexample.
  So: counterexample <=> MM(S) < 1/14 <=> the CLOSED radius-1/14 arcs cover [0,1].

We therefore directly test coverage by CLOSED arcs of radius 1/14 with a STRICT-gap
detector: a counterexample exists iff the closed arcs cover [0,1] AND there is no
single touch point -- i.e. MM<1/14 strictly. We hunt MM < 1/14 over:
  (1) random 13-sets with entries up to N (broad net),
  (2) hill-climb to MINIMIZE MM starting from random and from near-AP seeds,
  (3) structured 'covering-greedy' configs: pick speeds whose danger arcs (radius
      1/14) tile [0,1] as evenly as possible (the only way to beat 1/14 is to cover
      with margin everywhere, which needs the 13 speeds' 1/7-each coverage, total
      13/7 ~ 1.857, to overlap-cover [0,1] WITHOUT leaving a >1/14-deep gap).
"""
from fractions import Fraction as Fr
import random

def MM_float(S):
    """fast float max-min via dense+refine; conservative for screening."""
    # f(tau)=min_v ||v tau||; sample then golden-refine local maxima.
    import math
    def f(t):
        m=1.0
        for v in S:
            x=v*t; fr=x-math.floor(x); d=fr if fr<=0.5 else 1-fr
            if d<m: m=d
        return m
    N=4000
    best=0.0; bt=0.0
    prev=f(0.0)
    for i in range(1,N+1):
        t=i/N; cur=f(t)
        if cur>best: best=cur; bt=t
    # refine around bt
    lo=max(0.0,bt-1.0/N); hi=min(1.0,bt+1.0/N)
    for _ in range(60):
        m1=lo+(hi-lo)/3; m2=hi-(hi-lo)/3
        if f(m1)<f(m2): lo=m1
        else: hi=m2
    c=(lo+hi)/2
    return max(best,f(c))

THRESHf=1.0/14

print("=== Logic check: does ANY config approach MM < 1/14? ===")
print(f"   1/14 = {THRESHf:.8f}")

# (1) Broad random net
random.seed(12345)
best=1.0; bestS=None; hits=0
NTRIAL=25000
for trial in range(NTRIAL):
    N=random.choice([14,20,30,50,80,150])
    S=random.sample(range(1,N+1),13)
    mm=MM_float(S)
    if mm<best:
        best=mm; bestS=sorted(S)
    if mm<THRESHf-1e-7:
        hits+=1
        if hits<=10:
            print(f"   float MM<1/14 candidate: {sorted(S)} MM~{mm:.8f}")
print(f"   random ({NTRIAL}): min float MM = {best:.8f} at {bestS}; #(float<1/14)={hits}")

# (2) hill-climb to MINIMIZE MM
def neighbors(S):
    S=list(S)
    out=[]
    for i in range(len(S)):
        for d in (-1,1,-2,2,-3,3):
            T=S[:]; T[i]+=d
            if T[i]>=1 and len(set(T))==13:
                out.append(tuple(sorted(T)))
    return out
climb_best=1.0; climb_S=None
for start in range(60):
    N=random.choice([14,20,40,80])
    cur=tuple(sorted(random.sample(range(1,N+1),13)))
    curv=MM_float(cur)
    for _ in range(200):
        improved=False
        for nb in neighbors(cur):
            v=MM_float(nb)
            if v<curv-1e-12:
                cur,curv=nb,v; improved=True; break
        if not improved: break
    if curv<climb_best:
        climb_best=curv; climb_S=cur
print(f"   hill-climb-MIN: min float MM = {climb_best:.8f} at {sorted(climb_S)}")

print()
print("=== Exact verification of the lowest float candidates ===")
def MM_exact(S):
    cand=set()
    for v in S:
        for k in range(v): cand.add(Fr(2*k+1,2*v))
    Sl=sorted(S); n=len(Sl)
    for a in range(n):
        u=Sl[a]
        for b in range(a+1,n):
            w=Sl[b]
            for i in range(u+1):
                for j in range(w+1):
                    if u!=w:
                        t=Fr(i-j,u-w)
                        if 0<=t<=1: cand.add(t)
                    t2=Fr(i+j,u+w)
                    if 0<=t2<=1: cand.add(t2)
    best=Fr(0)
    for t in cand:
        m=None
        for v in S:
            x=v*t; fl=x.numerator//x.denominator; fr=x-fl
            d=fr if fr<=Fr(1,2) else 1-fr
            if m is None or d<m: m=d
        if m>best: best=m
    return best

for S in [bestS, sorted(climb_S)]:
    mm=MM_exact(S)
    print(f"   {S}: exact MM = {mm} = {float(mm):.8f}  ({'COUNTEREXAMPLE!' if mm<Fr(1,14) else 'not (>= 1/14)'})")
