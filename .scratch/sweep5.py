"""
Part 4: distance-graph / circular chromatic reformulation.

Classical LRC view: LRC for speeds v_1..v_m holds iff there is t with ||v_i t|| >= 1/n for all i,
where ||x|| = distance to nearest integer. Equivalently (Cantor/Wills) iff the point
(v_1 t, ..., v_m t) mod 1 visits the open box (1/n, 1-1/n)^m.

Distance-graph angle: Build a graph on Z/N where N = n * lcm(v_i). Map time t = a/N for
a in Z/N. We want a "coloring" c: Z/N -> {colors} ... The cleaner classical object is:
the circulant/Cayley graph X(Z_N, S) with connection set S = "bad differences", and ask for an
independent set hitting the lonely condition. We instead directly compute, on the time-lattice
Z/N (N=n*lcm), the indicator L(a) = 1 if observer lonely at t=a/N (i.e. ||v_i a/N||>=1/n for all i).

We also build the circular-chromatic reformulation: consider the graph H whose vertices are
the m runners and connect... that doesn't capture it. The known clean statement:

  LRC(speeds) TRUE  <=>  the set T_lonely = { t in [0,1) : N(t)=0 } is NONEMPTY.

We computed N(t) cell counts already. Here we tie it to a CIRCULAR COLORING bound:
the time-circle map t -> (point configuration) and the chromatic number of the danger graph
chi(G(t)) <= n is AUTOMATIC (n vertices). The nontrivial direction is min_t chi. We test:

  CLAIM A: observer lonely at t  <=>  the observer vertex is isolated in G(t)
           <=>  chi(G(t)) = chi(G(t) - observer).
  CLAIM B (circular chromatic for distance graph): For the *single-speed* danger structure of
  runner i, the bad set on Z/N is an arc; the union over i gives a circular-arc cover. LRC asks:
  is there a point of Z/N (a time) covered by NO runner-bad-arc? i.e. do the m arcs (each of
  measure 2/n) FAIL to cover the circle. Total measure = m*(2/n) = 2m/n = 2(n-1)/n = 2 - 2/n.
  For n>=2 that's >1, so measure does not preclude covering. The question is purely combinatorial
  arc-cover -> NOT a clean chromatic bound by measure alone. We quantify the cover.
"""
from lrc_coloring_core import *
from fractions import Fraction as F
from math import lcm
from functools import reduce

def lonely_set_cells(speeds, n):
    """Return list of cells (open intervals) where observer is lonely, and total measure."""
    crit = critical_times(speeds, n)
    mids = cell_midpoints(crit)
    cl = crit + [F(1)]
    intervals = list(zip(cl, cl[1:]))
    lonely = []
    for (a,b),tm in zip(intervals, mids):
        if observer_danger_count(speeds,n,tm)==0:
            lonely.append((a,b))
    meas = sum(b-a for a,b in lonely)
    return lonely, meas

# Circular-arc COVER reformulation: runner i is "dangerous to observer" when v_i t mod 1 in
# (-1/n, 1/n) i.e. t mod 1 in union of arcs around k/v_i. Observer lonely <=> t avoids all.
# We compute the cover measure and whether uncovered (lonely) region exists -> LRC for that set.
def cover_analysis(speeds, n):
    lonely, meas = lonely_set_cells(speeds, n)
    return len(lonely)>0, meas, len(lonely)

print("Part 4 / LRC arc-cover reformulation (true LRC instances should ALL be lonely-nonempty):")
print(f"{'speeds':28} {'n':>3} {'LRC(lonely?)':>12} {'lonely_measure':>16} {'#lonely cells':>14} {'2-2/n bound'}")
for speeds in [(1,2,3,4),(1,2,3,5),(1,3,4,7),(2,3,5,7),(1,2,4,8),
               (1,2,3,4,5),(1,2,3,4,6),(1,3,4,5,9),(2,3,4,5,6),
               (1,2,3,4,5,6),(1,2,3,4,5,7),(1,3,5,7,9,11),(1,2,4,8,16,32)]:
    n=len(speeds)+1
    lrc, meas, nc = cover_analysis(speeds,n)
    print(f"{str(speeds):28} {n:>3} {str(lrc):>12} {str(meas):>16} {nc:>14}   cover<= {2*(n-1)}/{n}")
