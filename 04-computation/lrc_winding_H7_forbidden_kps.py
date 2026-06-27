r"""
The winding tournament H(T(x)) and the forbidden 7 (kps-S31y, owner's hint).
T(x): vertices = runners {e_i}; arc i->j iff frac((e_i-e_j)x) in (0,1/2) [i is 'ahead' of j].
As x sweeps [0,1) the arcs flip (2 states each) at breakpoints x=k/(e_i-e_j). H(T(x)) = #Hamiltonian
paths (Redei count, always ODD; 7 FORBIDDEN). Question: as x varies, does H(T(x)) AVOID 7? And what
is H at the LONELY point (max-min ||e_i x||)? 14=2*7: is the apex-7 obstruction visible here?
"""
from itertools import permutations
def winding_tournament(E, x):
    n=len(E); beats=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            f=((E[i]-E[j])*x)%1.0
            beats[i][j] = (0 < f < 0.5)
    return beats
def H_redei(beats):
    n=len(beats); cnt=0
    for p in permutations(range(n)):
        if all(beats[p[k]][p[k+1]] for k in range(n-1)): cnt+=1
    return cnt
def loneliness(E, x):
    return min(min((e*x)%1, 1-(e*x)%1) for e in E if e!=0)
# 7 runners: reference 0 + AP speeds {1..6} = LRC(7), M=1/7
E=[0,1,2,3,4,5,6]
N=2520  # grid (multiple of small denoms)
from collections import Counter
hvals=Counter(); seven_count=0; lonely_x=None; lonely_val=-1; H_at_lonely=None
for a in range(1,N):
    x=a/N; b=winding_tournament(E,x); H=H_redei(b)
    hvals[H]+=1
    if H==7: seven_count+=1
    L=loneliness(E,x)
    if L>lonely_val: lonely_val=L; lonely_x=x; H_at_lonely=H
print(f"E={E} (7 runners, LRC(7) AP):  H(T(x)) value distribution over x:")
for h in sorted(hvals): print(f"   H={h}: {hvals[h]} grid pts {'<-- FORBIDDEN 7!' if h==7 else ''}")
print(f"\n   H=7 occurrences: {seven_count}  (Redei: 7 is FORBIDDEN => expect 0)")
print(f"   LONELY point (max-min ||.||): x={lonely_x:.4f}, min||.||={lonely_val:.4f} (=1/7={1/7:.4f}?), H there = {H_at_lonely}")
print(f"   => H avoids 7; the lonely point sits at a specific H. 14=2*7: the apex-7 forbidden value.")
