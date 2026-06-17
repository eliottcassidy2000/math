from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import sys

# Generalized danger arcs at gap 1/D (D=14 our case; D=13 the LRC(12) critical gap)
def danger_arcs_D(v,D):
    w=Fr(1,D*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A
def L_D(S,D):
    A=[]
    for v in S: A.extend(danger_arcs_D(v,D))
    A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

# PROVE-side lever: a 12-subset is a 12-speed lonely-runner instance.
# LRC(12) (gap 1/13) => safe set NONEMPTY for every 12-subset.
# Check: over 12-subsets of {1..14}, is any tight (meas 0) at gap 1/14? And at the critical gap 1/13?
print("12-subsets of {1..14}: lonely measure at gap 1/14 vs touching (meas 0) at gap 1/13")
sub14_tight14=0; sub14_tight13=0
mn14=None
both=[]
for C in combinations(range(1,15),12):
    m14=L_D(C,14)
    m13=L_D(C,13)
    if m14==0: sub14_tight14+=1
    if m13==0: sub14_tight13+=1
    if mn14 is None or m14<mn14[0]: mn14=(m14,C)
print(f"  # tight at gap 1/14 (meas=0): {sub14_tight14}")
print(f"  # tight at gap 1/13 (meas=0): {sub14_tight13}")
print(f"  min lonely measure at gap 1/14: {mn14[0]}={float(mn14[0]):.6f} at {mn14[1]}")
print(f"    == 7/858 ? {mn14[0]==Fr(7,858)}")
sys.stdout.flush()
