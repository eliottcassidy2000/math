# opus-2026-07-17-S373 -- THE CLEANER REFORMULATION: THE LARGEST UNCOVERED GAP.
#
# min-denominator is awkward (not dilation-invariant, and my searches kept
# finding worse: 25 random -> 32 -> 39).  The quantity that CONTROLS it is
# cleaner: if the uncovered set contains an interval of length L, that interval
# contains a rational of denominator ~1/L.  So
#     "minimal lonely denominator is bounded"  <==  "largest uncovered gap is
#      bounded BELOW by an absolute constant L0"
# and the latter is a natural strengthening of LRC(14): not merely that the
# uncovered set is nonempty, but that it contains a gap of definite size.
from fractions import Fraction as F
from functools import reduce
from math import gcd
import random
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def gaps(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    u=union(allv); g=[]; prev=F(0)
    for a,b in u:
        if a>prev: g.append(a-prev)
        prev=b
    if prev<1: g.append(1-prev)
    return g
def largest_gap(V):
    g=gaps(V)
    return max(g) if g else F(0)

print("(9) LARGEST UNCOVERED GAP on the adversarial witnesses and named families")
named=[("S373 hunt witness (min-den 39)",[31,32,36,74,145,210,231,260,304,459,500,552,616]),
       ("S373 witness (min-den 32)",[93,147,200,212,215,216,243,247,253,267,280,286,289]),
       ("THM-1055 primitive failure",[27,36,46,70,101,114,117,121,140,160,194,277,293]),
       ("{1,...,13} (tight)",list(range(1,14))),
       ("AP d=8",[1+8*i for i in range(13)]),
       ("odd {1,3,...,25}",[2*i+1 for i in range(13)]),
       ("THM-1060 L=31",[248+8*i for i in range(13)])]
for nm,V in named:
    L=largest_gap(V); tot=sum(gaps(V))
    print(f"    {nm:32s} largest gap = {float(L):.6f}   total uncovered = {float(tot):.6f}   1/L = {1/float(L) if L>0 else float('inf'):.1f}")

print()
print("(10) DISTRIBUTION OF THE LARGEST GAP over primitive families")
random.seed(3733)
vals=[]
for _ in range(150):
    V=sorted(random.sample(range(1,500),13))
    if reduce(gcd,V)!=1: continue
    vals.append(float(largest_gap(V)))
vals.sort()
print(f"    n={len(vals)}  min {vals[0]:.6f}   median {vals[len(vals)//2]:.6f}   max {vals[-1]:.6f}")
print(f"    smallest largest-gap found: {vals[0]:.6f}  =>  needs denominator ~{1/vals[0]:.0f}")
print()
print("  THE QUESTION THE ROUTE REDUCES TO: is the largest uncovered gap bounded")
print("  below by an absolute L0 > 0 over PRIMITIVE 13-families?  If yes with an")
print("  explicit L0, LRC(14) becomes a FINITE check on residues mod lcm(1..ceil(1/L0)).")
