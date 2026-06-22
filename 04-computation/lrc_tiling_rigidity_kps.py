"""
The complete-residue-system tiling rigidity (kind-pasteur S31n).
Safe(S) = meas{t in [0,1): ||s t|| >= 1/14  for all s in S}.  mac-mini: =0 iff S=d*{1..13}.
Exact via breakpoints of U_s = {||st||<1/14} = union_m ((m-1/14)/s,(m+1/14)/s).
"""
from fractions import Fraction as F

def safe_measure(S):
    S=[s for s in S if s!=0]
    bps=set([F(0),F(1)])
    for s in S:
        a=abs(s)
        for m in range(0,a+1):
            for sign in (-1,1):
                x=F(14*m+sign,14*a)
                if 0<=x<=1: bps.add(x)
    B=sorted(bps); safe=F(0)
    for lo,hi in zip(B,B[1:]):
        mid=(lo+hi)/2
        covered=False
        for s in S:
            r=(abs(s)*mid)%1
            if r<F(1,14) or r>F(13,14): covered=True; break
        if not covered: safe+=hi-lo
    return safe

if __name__=="__main__":
    AP=list(range(1,14))
    print("THE AP {1..13} and dilates (should be exactly 0 = exact tiling):")
    for d in [1,2,3,5]:
        S=[d*x for x in AP]; print(f"  {d}*{{1..13}}: safe = {float(safe_measure(S)):.6f}")
    print("\nPERTURBATIONS (should be > 0 = NOT tiling => the rigidity):")
    for name,S in [("{1..12,14}",list(range(1,13))+[14]),
                   ("{2..14}",list(range(2,15))),
                   ("{1..11,13,14}",list(range(1,12))+[13,14]),
                   ("{1..13} swap 7->15",[x for x in AP if x!=7]+[15]),
                   ("{1..13} swap 1->14",[x for x in AP if x!=1]+[14])]:
        print(f"  {name:22s}: safe = {float(safe_measure(S)):.6f}")
    print("\nGoddyn-Wong-type {1..11,13,24} (other tight family, expect ~0):")
    print(f"  {{1..11,13,24}}: safe = {float(safe_measure(list(range(1,12))+[13,24])):.6f}")
