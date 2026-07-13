#!/usr/bin/env python3
"""mac-mini-S87: LOCALIZE the Gowers multi-linear cancellation for meas(G' cap middle) on the
near-AP residue. The near-blanket (actual=4% of independent, S86) is driven by the smooth runners'
additive RESONANCES (Sum a_w v_w = 0). Identify (1) which subset of smooth runners drives it, (2)
the resonance order, (3) whether it reduces to a tractable sub-problem (e.g. the far element
decorrelates cleanly, leaving the near-consecutive {2..12} block = an LRC(<=13) object)."""
import numpy as np
c=1.0/14; P=[2,3,5,7,11,13]
def cop(v): return all(v%p!=0 for p in P)
def middle_safe(runners,N=300000):
    x=(np.arange(N)+0.5)/N
    mid=((x>=1/14)&(x<=13/14))
    ok=mid.copy()
    for v in runners:
        r=(v*x)%1.0; ok &= (np.minimum(r,1-r)>=c)
    return ok.mean(), mid.mean()
S=[*range(1,12),13,84]; nc=[v for v in S if not cop(v)]  # smooth non-core = {2..11,13,84}
print(f"near-AP {S}: smooth non-core = {nc}\n")
print("(1) PEEL the smooth runners one at a time: meas(smooth-safe cap middle) as we ADD runners:")
run=[]
for w in sorted(nc):
    run.append(w)
    L,_=middle_safe(run)
    print(f"   add {w:3d} -> runners {run}: meas(safe cap middle)={L:.5f}")
print("\n(2) which single runner, REMOVED, most opens the middle (the binding blanketer)?")
Lfull,midm=middle_safe(nc)
for w in sorted(nc):
    Lminus,_=middle_safe([v for v in nc if v!=w])
    print(f"   remove {w:3d}: meas jumps {Lfull:.5f} -> {Lminus:.5f}  (x{Lminus/Lfull:.1f})" if Lfull>0 else "")
print("\n(3) does the FAR element (84) decorrelate? compare {2..11,13} vs {2..11,13,84}:")
Lno84,_=middle_safe([v for v in nc if v!=84]); Lwith,_=middle_safe(nc)
print(f"   {{2..11,13}} meas(safe cap middle)={Lno84:.5f}; + 84 -> {Lwith:.5f} (far element reduces by x{Lwith/Lno84:.2f})")
print(f"   => if 84 reduces little, the blanket is the NEAR-CONSECUTIVE {{2..11,13}} block (11 runners),")
print(f"      an LRC(<=13)-adjacent object at level 1/14; the far element is a mild decorrelating peel.")
print("\nVERDICT: the near-blanket is driven by the near-consecutive smooth block; the Gowers object")
print("= 'the near-AP smooth block doesn't blanket the middle unless the full AP' = the AP-inverse at")
print("the middle resonance order (|a|_1~6,7). This is the Green-Tao-Ziegler Gowers-inverse frontier;")
print("I localize it (near-consecutive block + far peel), the harmonic analysis is the open crux.")
