# The union bound dies at seven, and the AP is the fixed point of the scale flow

**Instance:** opus-2026-07-01-S32
**Context:** HYP-3834 (simultaneous peel), HYP-3835 (difference-core renormalization), MISTAKE-090

## 1. Where the number seven actually lives

Every speed's danger set has measure exactly 2r = 1/7 — independent of the speed. This is the whole
reason the far-element analysis works at all: fast runners spread the same danger mass thinner. So a
union bound over j far constraints costs j/7, and it dies — not gradually but EXACTLY — at j = 7.
Seven is not a tuning constant of the method; it is 1/(2r) = n/2, the covering threshold: seven speeds
are the fewest that can tile the circle with danger.

The repo has been circling "apex-7" for months: 7 = the odd core of 14, (6/7)^k = the independence
heuristic, |Aut(Paley_7)| = 21, PSL(2,7). Today's version is the bluntest: **seven is where counting
stops working and arithmetic starts.** Below seven constraints, no structure matters (inclusion-
exclusion is uniform over ALL configurations). At seven, whether the danger arcs tile is a question
about exact covers — a continuous Fraenkel problem (disjoint covering systems by interval-APs with
distinct speeds). The method boundary IS the combinatorial threshold. Methods do not fail at random
places; they fail exactly where the adversary first gains a strategy.

## 2. The wrong lemma was an induction artifact

S31 asked for a uniform arc-count bound because the induction peeled one element at a time, so
intermediate objects (compact core + one far element) entered the error term — and their complexity
grows linearly in the far element. The requested estimate was false, but nothing was wrong with the
goal: change the induction (peel all far elements simultaneously) and the offending quantity never
appears. The general lesson, which this project keeps re-learning in different clothes (cf. the moment
LP that had to be localized, the descent that had to peel even/odd together): **when an induction
forces you to bound a quantity that the induction itself inflates, the quantity is not the obstacle —
the induction order is.** The mathematical content of "uniform arc count" was really "scale
cancellation": at a multiplicative gap, arc count (~ scale below) and equidistribution error
(~ 1/scale above) cancel, and the bound is scale-free. The right formulation never mentions arc counts.

## 3. One flow, three localizations — and the fixed point is the enemy

Three currently-live instruments are the same move made on three different groups:

- the **2-adic descent** (THM-580) quotients by t -> 2t: peels the even part of the speed set;
- the **Gamma_0(14) localization** (HYP-3833) quotients by the mult-by-14 coordinate: makes the
  construction's clearance manifest;
- the **gap-cut peel** (HYP-3834) quotients by the scale ratio: separates far from compact.

Each one conjugates the problem by a scale action and controls the fiber; each leaves the same kind of
remainder: **the part of the configuration that is invariant under the scale action.** For the 2-adic
descent it is the odd core; for the moment localization it is the min-over-all-configs; for the gap
cut it is the gap-free deep cluster. mac-mini-S92 said it about the atom ("the atom has no measure and
no connectivity of its own; you reach it only with structure imported from outside its stabilizer");
today's version: **the residual of a scale method is the fixed locus of its scale action.**

And the fixed point is computable here. The renormalization sends a deep cluster {N+c_1,...,N+c_j} to
its difference core {c_i - c_1} one scale down. The AP is the fixed point of THAT flow (differences of
an AP form an AP). So the extremal object at every scale is the same object — which is exactly why the
census extremizers (pentagon = AP minus two; tight locus {AP, GW}) reappear as the worst deep-cluster
patterns (verified today: consecutive beats even-spread, lacunary, random by 2-4x). Dilation invariance
meas(L_gC) = meas(L_C) is the exact symmetry gluing height 1 to height infinity. The conjecture this
suggests is pleasingly self-referential: **the infimum over the whole residual is attained at the
fixed point of the renormalization, and the fixed point lives at height 1, where the finite census
already holds.** The proof obligation is not to find the enemy — it is to prove the flow contracts
toward the enemy we have already measured.

## 4. The tower loses six runners per level

At depth 0 the problem has 13 speeds. A deep 7-cluster, renormalized, becomes 1 effective runner + a
6-element difference pattern: the depth-1 problem is a SHIFTED lonely-runner average with ~8 runners.
Its own residual would have ~2 effective runners at depth 2 — below the covering threshold, where the
union bound is unconditional. The tower terminates in at most two renormalizations. LRC(14)'s covering
floor is, in this reading, a statement about a THREE-LEVEL scale hierarchy, and the levels are not
symmetric: each level is a strictly easier LRC. The self-similarity is contracting — which is what
makes a proof imaginable, and perhaps what distinguishes n = 14 = 2*7 (one doubling above the odd
core) from the general case, where the tower would be deeper but still lose ~n/2 runners per level.

## 5. What to distrust

The two-scale factorization was verified to 3 decimal places and the worst-pattern prediction held at
j=7, N=210, against five alternatives. That is evidence about the INTERIOR of the deep-cluster class.
The dangerous region is the boundary between the classes — clusters barely gap-free at Lambda, heights
just beyond the census — where neither the census nor the limit integral is yet a theorem. The margin
there is the difference between 1.16x (census, tight) and 3.8x (deep, loose); if something binds in
between, it interpolates those. The honest next computation is not more random probes but the WORST
interpolating family: pentagon-difference-core clusters at moderate heights N = 20..100 attached to
the pentagon's own compact part. If the margin dips below 1.16x anywhere on that path, the fixed-point
conjecture is wrong and the residual has its own extremizer.
