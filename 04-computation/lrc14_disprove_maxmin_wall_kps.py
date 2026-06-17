"""
DISPROVE-side: the max-min WALL and the near-counterexample hunt.

A genuine LRC(14) counterexample = a 13-set S with
    MM(S) := max_tau min_{v in S} ||v tau|| < 1/14.

EXACT max-min: f(tau)=min_v ||v tau||. ||v tau|| as a function of tau is a sawtooth
with peaks (value reaching up toward 1/2) at tau=(2k+1)/(2v) and zeros at tau=k/v.
f(tau) is piecewise linear; its breakpoints are the union over v of {k/v} (zeros,
where that v's contribution is 0 and rising) and {(2k+1)/(2v)} (that v's local peaks).
Between consecutive *global* breakpoints every ||v tau|| is linear, so f (a min of
linear pieces) is concave and its max on each cell is at a cell endpoint OR at an
interior intersection of two of the linear pieces. We therefore enumerate:
  - all breakpoints B = {k/v} U {(2k+1)/(2v)} for v in S,
  - all pairwise intersections of the active linear pieces within each cell.
Evaluating f exactly at the finite candidate set and taking the max gives MM exactly.

Simpler exact equivalent we use: MM = sup{r : radius-r arcs don't cover [0,1]}.
For each maximal GAP (uncovered interval) of the radius-(1/14) arcs we'd get nothing
(if covered). So instead we directly compute candidate values: for every pair of
speeds (u,w) and integers (i,j) the point tau where ||u tau|| = ||w tau|| as the
binding constraints, the common value is a candidate for MM. We bound i,j by the
speeds. We ALSO include single-speed peaks (2k+1)/(2v) -> value 1/2 clipped by others.

To stay exact and tractable we use: candidate tau set =
   { (i*w + j*u) / (2*u*w) : signs } intersections of the two sawtooth lines,
for all pairs u<w in S and relevant i,j. We evaluate f exactly at each candidate in
[0,1] and take the max. Verified against known values (AP -> 1/14, sporadic -> ...).
"""
from fractions import Fraction as Fr
from math import gcd, lcm

def f_exact(S, tau):
    """min_v ||v tau||, exact."""
    best=None
    for v in S:
        x=v*tau
        fl=x.numerator//x.denominator
        fr=x-fl
        d=fr if fr<=Fr(1,2) else 1-fr
        if best is None or d<best: best=d
    return best

def maxmin_exact(S):
    cand=set()
    # single-speed peaks (clipped) and zeros
    for v in S:
        for k in range(v):
            cand.add(Fr(2*k+1,2*v))  # peak of v's sawtooth
    # pairwise intersections of the two sawtooth linear pieces:
    # ||u tau|| rising = u*tau - i ; ||w tau|| rising/falling = +-(w*tau - j)
    # Solve u*tau - i = w*tau - j  -> tau=(i-j)/(u-w); and u*tau-i = -(w*tau-j) -> tau=(i+j)/(u+w)
    Sl=sorted(S)
    n=len(Sl)
    for a in range(n):
        u=Sl[a]
        for b in range(a+1,n):
            w=Sl[b]
            # i in 0..u, j in 0..w
            for i in range(0,u+1):
                for j in range(0,w+1):
                    if u!=w:
                        t=Fr(i-j,u-w)
                        if 0<=t<=1: cand.add(t)
                    t2=Fr(i+j,u+w)
                    if 0<=t2<=1: cand.add(t2)
    best=Fr(0); bt=None
    for t in cand:
        val=f_exact(S,t)
        if val>best: best=val; bt=t
    return best,bt

THRESH=Fr(1,14)
AP=list(range(1,14))

print("=== max-min EXACT for key families ===")
tests={
 "AP {1..13}": AP,
 "sporadic {1..11,13,24}": [1,2,3,4,5,6,7,8,9,10,11,13,24],
 "champ {1..11,13,36} L=1/1260": [1,2,3,4,5,6,7,8,9,10,11,13,36],
 "drop10 ->20  L=1/980": [1,2,3,4,5,6,7,8,9,11,12,13,20],
 "drop6 w=69   L=19/10626": [1,2,3,4,5,7,8,9,10,11,12,13,69],
}
for name,S in tests.items():
    mm,t=maxmin_exact(S)
    rel = '<' if mm<THRESH else ('=' if mm==THRESH else '>')
    print(f"  {name:34s} MM = {mm!s:>10s} = {float(mm):.6f}  {rel} 1/14   at tau={t}")
print(f"  1/14 = {float(THRESH):.6f}")
