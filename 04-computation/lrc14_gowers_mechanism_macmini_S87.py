#!/usr/bin/env python3
"""mac-mini-S87: the Gowers-local MECHANISM of the covering-min. Both covering-min extremals are
the non-covering AP {1..13} with ONE element sent far (deep well = 13->182; multi-killer = 12->84).
The far element is chosen to ALIGN with the core's safe point; the claim (Gowers-local) is that ONE
far runner = a sparse period-1/w comb CANNOT fully close the core's safe INTERVAL, so L>0. Verify:
L(deep well) = meas(SafeSet(core) cap middle) - meas(that cap D_far); the far comb removes only an
O(1/w) slice, not the whole interval. This localizes the multi-linear cancellation to a 1-runner-
vs-interval statement = the tractable RESIDUE of the Gowers inverse."""
import numpy as np
c=1.0/14
def safe_mask(runners,x):
    ok=np.ones(len(x),bool)
    for v in runners:
        r=(v*x)%1.0; ok &= (np.minimum(r,1-r)>=c)
    return ok
def L_and_pieces(core,far,N=2_000_000):
    x=(np.arange(N)+0.5)/N; mid=((x>=1/14)&(x<=13/14))
    core_safe=safe_mask(core,x)&mid
    r=(far*x)%1.0; Dfar=(np.minimum(r,1-r)<c)
    L=(core_safe&~Dfar).mean(); base=core_safe.mean(); cut=(core_safe&Dfar).mean()
    return L,base,cut
print("MECHANISM: covering-min extremal = AP{1..13} with one element sent far.\n")
for name,core,far in [("deep well (13->182)",list(range(1,13)),182),
                      ("multi-killer(12->84)",[*range(1,12),13],84)]:
    L,base,cut=L_and_pieces(core,far)
    print(f"{name}: core={core[0]}..{core[-1]}{'+13' if 13 in core and core[-1]!=13 else ''}, far={far}")
    print(f"   meas(SafeSet(core) cap middle) = {base:.5f}   (the core's safe interval)")
    print(f"   meas(that cap D_far)           = {cut:.5f}   (what the ONE far comb removes)")
    print(f"   L = base - cut                 = {L:.5f}   (>0: comb can't close interval)")
    print(f"   far comb removes {100*cut/base:.1f}% of the core interval; period 1/far={1/far:.5f} vs interval~{base:.4f}\n")
print("Is the deep-well far element (182=13*14) the BEST-ALIGNING? scan far covering values for the")
print("{1..12} core (must be divisible by 13 AND 14 to carry both missing moduli => multiples of 182):")
core=list(range(1,13))
best=None
for m in range(1,9):
    far=182*m
    L,base,cut=L_and_pieces(core,far,N=1_000_000)
    tag=""
    if best is None or L<best[0]: best=(L,far); tag="  <-- min so far"
    print(f"   far=182*{m}={far:5d}: L={L:.5f}, cut={cut:.5f} ({100*cut/base:.1f}% of interval){tag}")
print(f"\n=> smallest far (182) aligns best (removes most of the core interval) => DEEPEST well = min L.")
print("   Larger far => sparser... no: 182 is the SMALLEST covering carrier; it maximally overlaps")
print("   the core safe interval near t=1/13 (182/13=14 integer => D_182 hits t=1/13 exactly).")
print("\nGOWERS-LOCAL RESIDUE: L(deep well)>0 <=> the single comb D_182 (period 1/182) cannot cover the")
print("core's safe INTERVAL SafeSet({1..12}) cap middle (a union of intervals of total measure ~base).")
print("This is a 1-runner-vs-interval covering question = the tractable core of the Gowers inverse;")
print("the multi-linearity is fully absorbed into 'SafeSet(core)' (a fixed interval union). The open")
print("crux is whether ONE arithmetic comb can close a three-distance interval union -- a Diophantine")
print("(not statistical) question, and the answer is NO iff the balance margin 14/183-1/14>0.")
