# Exact strong span-2 ribbon floor and cycle ceiling

**PROVED (symbolic integer certificate).** Let \(T\) be a strong
tournament that admits a span-2 ribbon presentation: in the presentation
order, every backwards arc has span at most two. For \(n\ge3\) let
\(m(n)\) be the minimum Hamiltonian-path count \(H(T)\) in this class,
and let \(C(T)\) count directed Hamiltonian cycles up to cyclic rotation.
Then
\[
\begin{aligned}
m(n)&=f(n) &&(3\le n\le15),\\
m(16+2k)&=3375\,3^k &&(k\ge0),\\
m(17+2k)&=5625\,3^k &&(k\ge0),\\
\min_T (H(T)-2C(T))&=m(n)-2 &&(n\ge3).
\end{aligned}
\]
Here \(f\) denotes the [Busch floor for strong tournaments](https://doi.org/10.37236/1141),
\(f(3)=3,f(4)=5,f(5)=9,f(n+3)=5f(n)\). In particular
[HYP-9457(a)](../../hypotheses/HYP-9457-single-cycle-floor-band-and-lorentzian-join-crevasses.md)
holds at **every order in span 2**, including its equality clause:
\(H(T)\ge f(n)+2(C(T)-1)\), with equality if and only if \(H(T)=f(n)\).
The inequality is strict throughout span 2 at \(n=16\) and \(n\ge18\);
the other orders have at least one equality ribbon.

The stronger exact floor \(m(n)\) eventually grows faster than \(f(n)\):
their ratio multiplies by \(27/25\) under \(n\mapsto n+6\).
Thus at a fixed large order it excludes an entire near-floor band
below the ribbon minimum, complementing the existing absolute ribbon
obstruction at 47.
This does **not** exclude these integers at a different order or in
general tournaments, prove the span-3 version, or prove Chan–Pak.

## 1. Positive transfers and the two-rotation subtraction

A letter is a two-bit mask on the backward arcs from the newest vertex
at distances 1 and 2. The
[existing ribbon transfer](../ribbon_tournaments_bounded_span_transfer_20261104.md)
builds a path by joining fragments. After reading a prefix, restrict
the prospective path to its old vertices and retain its fragments
as pairs (head age, tail age), where age 0 is newest. An endpoint
outside the two-vertex window is recorded as -1. A head outside the
window cannot receive an arc from a future vertex, so two such heads
rule out eventual completion to one path. An old tail can attach to
every future vertex. A new vertex can be a singleton, prepend,
append, or join two distinct fragments; the same-fragment join would
close a cycle and is not used. Forget endpoints once they age out.
Thus the transfer is exact and nonnegative, and its reachable
27-coordinate profile vector counts histories with multiplicity.
The terminal row for (H) sums one-fragment states.

For (H-2C), do not subtract two separate counts in the transfer.
Instead, mark precisely those paths which start at either of the
last two vertices of the given presentation and close by the arc
from their tail back to their head. Each directed Hamiltonian
cycle gives exactly two such paths by rotation, and each such path
closes one cycle. Therefore (H-2C) has a **nonnegative terminal
row** on the same profile, indexed by the last two letters. Apply
the ordinary fragment transitions for these two letters but
increase endpoint retention to ages 0..2 on the first and ages 0..3
on the second. At termination, a head of age 0 or 1 is in the
last two vertices, and the closing arc can be computed from the
two letters and the recorded tail age. If the tail is old, the
closing arc is always forward. Discard exactly the marked paths.
This terminal construction works at any width (w), retaining
(w+1,w+2) vertices at the final two steps. More generally, selecting
the (j) last possible starting vertices yields a nonnegative
terminal row for (H-jC) with a (j)-letter ending, for \(j\le n\).

For strong connectivity, carry the existing condensation automaton
alongside the profile. Each state's SCC list is in its ordinal
order, recording the ages of its window vertices and whether it has
any older vertex. When the new vertex beats any component vertex, let
(i) be the first component it hits. It merges with the interval of
components from (i) through the last component that can send an arc
to the new vertex; a terminal suffix entirely beaten by the new
vertex is left alone. Indeed every component earlier than (i)
dominates the newcomer; a directed cycle through the newcomer and a
hit at (i) includes the entire intervening interval, and no beaten
terminal component can enter that cycle. If the merge interval is
empty, insert a singleton at (i); if there is no hit, append a
singleton. A first component with no window vertex will never merge
again and can be killed, since all future arcs at that distance
point away from it. At termination one component is needed. This is
an exact finite state test; notably, it does not equate strength with
just checking the original order cuts.

## 2. The induction certificate: affine in (3^k), not scaling by 3

For each live condensation state (s), keep only the componentwise
minimal reachable profile vectors at each length. Nonnegative
transfers preserve domination, and both terminal rows are
nonnegative. Thus this exact pruning never loses a smaller value
of either objective. Every retained vector is physically attained.

From lengths 16,18,20 and 17,19,21 respectively, fit a finite
family of symbolic templates in each condensation state, of form
\[
 L_{v,d}(k)=v+\frac{3^k-1}{2}d,\qquad k\ge0,
\]
where \(v\) is the initial vector and \(d\) is its change in the next
two steps. Fit by exact pairing from \(v,w,4w-3v\) in the three
frontiers; this is not assumed to continue. Pure multiplicative
scaling of whole profiles would fail because some coordinates
stay equal to one. The affine term preserves them.

For each feasible two-letter product \(A\) from \(s\) to \(t\) and each
template \((v,d)\) at \(s\), the replay verifies that some template
\((z,e)\) at \(t\) satisfies
\[
 Av\ge z+e,\qquad Ad\ge3e
\]
componentwise. These imply \(A L_{v,d}(k)\ge L_{z,e}(k+1)\) for every \(k\ge0\). Conversely, for every target template there exists
a predecessor and a two-letter product for which **both**
coefficient relations are exact equalities. This proves physical
attainment at every level, not just a lower envelope. It also checks
that each ordered pair of different templates in the same state has
a coordinate where the first initial value is strictly larger and
its slope is at least as large. Thus every template stays on the
Pareto antichain. The finite inequalities extend the **whole**
antichain, in both parities, to every order.

The even bank has 77 templates and 980 checked transitions; the odd
bank has 83 and 1058. For the (H) row over the corresponding
accepting templates, the least intercept and least slope bound the
same line, with ((m,2m)) for (m=3375,5625) respectively.
For the (H-2C) row applied after any admissible two-letter ending,
the corresponding attained pair is ((3m-2,6m)).
Both statements are checked as domination of intercepts and slopes,
not numerical fitting. This gives exactly the claimed minima. The
initial rows through 18 are checked directly on the finite fronts.

## 3. Equality and independent controls

Compare (m(n)) and (f(n)) on \(16\le n\le21\), then use the exact
ratio factor (27/25>1) for each advance of six. Only (n=17)
can be an equality order above 15. At the finite equality orders,
keeping only minimum vectors would not in itself prove the
rigidity clause: a discarded vector can tie a terminal row with
zero coordinates. Instead, propagate Pareto-minimal (H-2C) output
rows backwards over all accepting endings, then traverse the
forward word tree. A prefix vector is pruned only if its dot
product with **every** minimal remaining row exceeds (f(n)-2).
This preserves all tied words. Every retained word at the
threshold indeed has (H=f(n)). Conversely if a strong ribbon has
(H=f(n)), it has (Cge1), and the proved inequality gives
(C=1) and equality.

[Verifier](verify_near_floor.py): run with standard-library
Python 3.8+; checks raise exceptions even in optimized mode.
It builds the matrices and the condensation transitions from scratch,
independently checks (H,C) and strong connectivity on all 680
distinct labelled ribbon graphs of orders 3–6 by permutation counting
and reachability, regenerates the banks, and checks the symbolic
coefficient relations. Its listed equality multiplicities count
presentations with the ignored leading letter bits, not
isomorphism classes.

An unrestricted proof would still need a theorem passing from
general strong tournaments to an appropriate bounded interface (or
a new argument for unbounded interface). A finite prefix of an
unbounded-width enumeration is not a substitute for that step.
