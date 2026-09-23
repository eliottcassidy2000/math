# Direct test of candidate hostile rational points of the forward E-game (2-adic):
# exhaustive DFS over E-paths from x0 (exact rational arithmetic), up to B halvings,
# pruning when cost - (remaining halvings)*log2 >= 0.  Reports whether any prefix has 3^a < 2^b.
from fractions import Fraction as F
import math, sys
L2,L3=math.log(2),math.log(3)
def two_adic_even(x):  # x = p/q with q odd
    return x.numerator % 2 == 0
def descends(x0, BMAX=40, nodelim=2_000_000):
    stack=[(x0,0,0)]; nodes=0
    while stack:
        x,a,b=stack.pop(); nodes+=1
        if nodes>nodelim: return None
        cost=a*L3-b*L2
        if a+b>0 and cost< -1e-12: return (a,b)
        if b>=BMAX: continue
        if cost-(BMAX-b)*L2>=-1e-12: continue
        if not two_adic_even(x):
            stack.append((3*x+1,a+1,b))
        else:
            stack.append((3*x+1,a+1,b))
            stack.append((x/2,a,b+1))
    return False
cands=["-1","-13/9","-35/27","-97/81","-113/81","-275/243","-307/243","-355/243","-371/243","-793/729","-857/729","-953/729","-985/729","-1049/729","-1081/729","-2315/2187","-2443/2187","-6817/6561",
       "-5/3","-7/9","-11/9","-17/9","-19/9","-23/9","-25/27","-29/27","-31/27","-37/27","-41/27","-43/27","-1/3","-2/3","-4/3","-7/3","-3/5","-7/5","-9/7","13/9","1/3","-3/2" ]
for c in cands:
    x=F(c)
    if x.denominator%2==0: print(c,"not 2-adic integer"); continue
    r=descends(x)
    print(f"{c:>12} = {float(x):8.4f}: ", "NO descent within 40 halvings (hostile candidate)" if r is False else ("node limit" if r is None else f"descends at (a,b)={r}"))
