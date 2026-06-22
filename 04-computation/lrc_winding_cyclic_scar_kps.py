"""
LRC = the parity-dual clique scar's runner-shadow (kind-pasteur S31i).

Creative application of the parity-dual clique-scar structure to LRC(14):
  14 = 2 * 7 = (parity/even factor) * (the unique forbidden odd-clique value 7 = I(K_3,2)).
The LRC complement-halving (x -> -x, Z/2) makes the 7 inner sectors -- the SAME Z/2 as the
even-graph "forget orientation".  The 7 sectors = the apex prime = the odd-clique-scar value.

GROUNDING TEST (HYP-2605, the winding tournament): LRC loneliness at phase x = the winding
tournament T(x) has a scale-1/7 LOCAL SINK (a >1/7 empty arc).  A sink = HIERARCHICAL
(transitive-ish, FEW 3-cycles).  NON-loneliness = no sink = CYCLIC (MANY 3-cycles), the
direction of the forbidden purely-cyclic K_3 apex.  So: lonely <-> low cyclic content;
non-lonely <-> high cyclic content (toward the K_3 scar).  Test on a small AP cluster.
"""
import math
from itertools import combinations

def winding_tournament(e, x):
    k = len(e)
    beats = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i==j: continue
            f = ((e[i]-e[j])*x) % 1.0
            beats[i][j] = (0 < f < 0.5)
    return beats

def c3(beats):  # number of directed 3-cycles
    k = len(beats); cnt=0
    for a,b,c in combinations(range(k),3):
        # count cyclic triangles among the 3
        for (x,y,z) in [(a,b,c)]:
            if beats[x][y] and beats[y][z] and beats[z][x]: cnt+=1
            if beats[x][z] and beats[z][y] and beats[y][x]: cnt+=1
    return cnt

def maxgap(e, x):
    pts = sorted(((ei*x)%1.0) for ei in e)
    if not pts: return 1.0
    gaps = [pts[(i+1)%len(pts)]-pts[i] for i in range(len(pts)-1)]
    gaps.append(1.0 - pts[-1] + pts[0])
    return max(gaps)

def run(e, N=20000):
    k=len(e)
    lonely_c3=[]; nonlonely_c3=[]
    for n in range(1,N):
        x = n/N
        g = maxgap(e,x)
        c = c3(winding_tournament(e,x))
        if g > 1/7: lonely_c3.append(c)
        else: nonlonely_c3.append(c)
    import statistics as st
    print(f"cluster e={e}  (k={k}):")
    print(f"  LONELY (maxgap>1/7, has 1/7 sink): frac={len(lonely_c3)/(N-1):.3f}, "
          f"mean #3-cycles={st.mean(lonely_c3) if lonely_c3 else 0:.3f}")
    print(f"  NON-lonely (no 1/7 sink, cyclic):  frac={len(nonlonely_c3)/(N-1):.3f}, "
          f"mean #3-cycles={st.mean(nonlonely_c3) if nonlonely_c3 else 0:.3f}")
    if lonely_c3 and nonlonely_c3:
        print(f"  => non-lonely has {st.mean(nonlonely_c3)/max(st.mean(lonely_c3),1e-9):.2f}x the "
              f"cyclic content of lonely  (cyclic = toward the K_3 scar)")

if __name__=="__main__":
    for e in (list(range(8)), list(range(10)), list(range(12))):
        run(e)
        print()
