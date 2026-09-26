# Small arithmetic graphs, families beyond 27, and a natural-density minimum

**Status: PROVED + independently audited, with separate FINITE-EXACT controls;
ordinary Collatz OPEN.** Session 2026-09-26, kind-pasteur. Two canonical
records were reserved in commit `6d5865d4e` and promoted after audit:
[THM-4500, pairing Bellman/natural density](../../01-canon/theorems/THM-4500-pairing-bellman-contraction-and-natural-density.md)
and [THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md).

## 1. Outcome and inheritance

The strongest advance is a complete natural-density optimization result
for the repository's modified two-step pairing model. A contracting
Bellman operator gives the exact limiting constant and certified bounds;
a slow variation of policy along logarithmic phase constructs one global
pairing **attaining** that natural density. This closes the natural
attainment boundary left open by the previous session. It is not a proof
about the unmodified Collatz map.

For the user's family question, we distinguish three objects: numbers
that actually reach 27, numbers sharing its later tail, and numbers with
arbitrarily long initial growth. The last group admits an explicit
recursive construction carrying a small prime graph, with a proved
occurrence rate. The fractal language of all expanding parity prefixes
has a different rate and a different limiting set.

Inheritance pass:

- Closest proved mechanisms: THM-4496's harmonic pairing optimum and
  phase repairs; THM-4497's weighted prime-graph kernel; THM-4495's sharp
  no-descent word count; the corrected square-sum endpoint census.
- Canonical hostiles: 27 has no odd predecessor; a selected square-sum
  path omits graph edges; arithmetic progressions lose their odd CRT
  factor upon dyadic closure; a Haar equivalence class does not specify
  its values at ordinary integers.
- Corrected near misses: long total stopping time versus long first
  descent; positive finite realizations versus one infinite realization;
  isolated optimal policies versus a globally compatible limiting policy.
- Least-used coordinates: ancestor resets, logarithmic phase, odd-prime
  residue clocks, unused square edges, and the ambient completion.

Portfolio: the anchor was the family's rate and small arithmetic carrier;
the niche was square-sum path rigidity and local exchanges; the wildcard
was the signed harmonic root-cost difference. The wildcard produced the
strongest new theorem and became the principal target.

Incoming main was integrated through `e9aef5afe` before promotion.
THM-4498's independent audit repaired an endpoint inequality and a floating
tie test; THM-4499 was independently audited SOUND. This session's fractal
corollary uses the already audited THM-4495, not a new inference from the
thin-orbit theorem. The latter concerns visits of a single injective orbit,
whereas our occurrence counts concern many different starting integers.

## 2. A family carrying a ten-node arithmetic graph

[The orbit note](crossroads_family_20260926_orbits.md) proves:

1. The positive ancestors of 27 are exactly `27*2^r`.
2. Its next odd node 41 has the direct predecessor comb
   `27,109,437,1749,...`, with recursion `n ->4n+1`. These share a long
   tail, but all except 27 descend within two shortcut steps.
3. A stronger construction uses base 4347, whose first ten odd values
   contain the pair-exclusive prime triangle `7,23,29`. Let

       M=19124224,
       F_k={2^k u-1: u>=4348, u=4348*3^(-k) mod M}.

   Every member has at least `k+11` initial shortcut steps above its
   starting value, with the designated prime motif in the block beginning
   after k steps. The disjoint
   tail union has the exact frequency

       #((union_(k>=K)F_k) intersect [1,X])
         =2^(1-K)X/19124224+O(log X).

   The phase clock has period 59136, giving exact set recursion
   `F_(k+59136)=2^59136(F_k+1)-1`.

The small graph preserves prime-incidence constraints, while the phase
clock preserves congruence and an unbounded scale supplies height. These
are distinct coordinates; the ten vertices alone do not store the entire
arithmetic family.

The tail closures in Z_2 have Haar measure **4669 times** their natural
density, because the odd congruence factor disappears. Every closed tail
set has dimension 1, yet their intersection is only `{-1}`. The positive
integer sets have no common member. This is an exact hostile to inferring
one divergent integer from increasingly long finite examples.

## 3. What the aggregate occurrence rate really says about fractal recursion

[The fractal note](crossroads_family_20260926_fractal.md) combines exact
parity coding with THM-4495. If W_k counts words whose multiplier exceeds
one at every nonempty prefix, then

    W_k=Theta(2^(hk)k^(-3/2)),
    h=H_2(log_3 2)=0.949955527188... .

Thus the fixed-horizon density is
`Theta(2^(-(1-h)k)k^(-3/2))`. The corresponding infinite dyadic set E has
Hausdorff and box dimension h, but its h-dimensional Hausdorff measure is
zero. Dimension h was already present in the repository; the exact cover
order and critical-measure corollary use the newer sharp count.

A finite arithmetic graph genuinely approximates E: use one vertex with
one labelled loop for each allowed b-bit block, carrying its exact affine
inverse. Arbitrary loop concatenations remain expanding; their dimensions
`log_2 W_b/b` tend to h. Increasing b increases the graph's alphabet, so
this is not a claim of a bounded-state exact recognizer for E.

The constructed motif family has frequency of order `2^(-k)`, far smaller
than the full long-growth language. It exposes a rigorous recursive piece,
not the dominant source of that language's entropy. Moreover known
convergent integers occur in every finite parity cylinder: the basin of
41 is dense in Z_2. Finite parity realizability therefore cannot alone
distinguish eventual convergence from a putative exception.

One explicit two-loop hostile makes the ambient-space issue tangible:
`(2y-1)/3,(2y-2)/3` have real attractor `[-2,-1]` and dyadic attractor
`Z_2`. The same rational prefix sequence can converge to a negative real
and to the positive dyadic integer 1. These are control maps, not Collatz.

## 4. The density result that emerged from the graph approach

[The flow note](crossroads_family_20260926_flow.md) proves that the global
P2 pairing constant B_* obeys

    0.3166003458 < B_* < 0.3166306675.

The exact rational endpoints come from 2^24 residual states, with a
uniform proof of the remaining error. The result is stronger than a
numerical search: a Haar L1 contraction computes B_*, complete finite
policy LPs approach it at a known rate, and one explicit global pairing
has count

    A_F(X)=B_*X+O(X/(log X)^(1/3)).

The construction chooses increasingly precise policies slowly along
logarithmic phase. A parent roughly multiplies an index by 2/3, so it
subtracts one from log_(3/2) of that index and nearly preserves its
fractional part. Outside a vanishing set of phase boundaries, several
ancestor steps use the same policy; resets erase earlier mismatches.
This retains a single globally consistent positive-index construction.

The cheap hostile was decisive: directly patching the stationary bit
sequences fails the complement constraint at indices **26 and 39**.
Recomputing bits recursively with the phase-selected free choices gives
bits `(0,1)` there and satisfies all constraints. A proof of small average
cost was insufficient until that legality issue was repaired.

This is the attained minimum natural density and minimum upper natural
density for the modified model. The smaller minimum lower natural density
from THM-4491 remains a different quantity. Finite-period attainment and
rationality of B_* remain open.

## 5. Why the square-sum list begins at 15, and where local flexibility begins

[The square note](crossroads_family_20260926_squares.md) independently
recovers the unique nontrivial first path

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

The full graph Q_15 is a five-cycle with two pendant paths attached at
adjacent cycle vertices. Deleting their common cycle edge leaves this
unique spanning path. The leaves force endpoints 8 and 9. This also
explains the adjacent successes at 16 and 17 and failure at 18.

The first square-sum four-cycle occurs at **46**, forced by the first
relevant square-label collision `4+81=36+49`. An explicit Hamiltonian
cycle at 46 can be changed by a two-edge switch to another Hamiltonian
cycle, so 46 is also the sharp first switch threshold. Numerous cycles
below 46 can exist while remaining isolated under this particular move.

The exact bridge to the prime graph is the signless operator. Under
`x'_v=q^2 x_v+c_v`, retaining an edge's designated scaled square root
requires `c_u+c_v=0`. A bipartite component has an alternating offset;
an odd cycle forces that offset to zero. This is the equal-exponent
specialization of THM-4497's valuation equations and the same linear
operator in THM-2521's LRC potential note. Its variables have different
meanings in those settings; no trajectory-preserving map follows.

For Q_15 the selected path lifts, but its unused odd-cycle edge obstructs
lifting the whole graph. Restoring discarded equations is the useful
operation to transfer to Collatz investigations.

## 6. Live board after the synthesis

| Concept | Preserved object / operation | Obstruction and strongest next test |
|---|---|---|
| Bellman density | Legal pairing tree; exact harmonic telescope and resets | Extend to three-step descent: first determine whether the constraint carrier still has a reset or bounded-repair property |
| Motif occurrence | Actual orbit prefixes; CRT lift and odd-step prepend | Measure motif incidence across increasing horizons; keep one source, all supports, carries, and height |
| Entropy fractal | All finite expanding words; block concatenation | Separate ordinary integer intersection from topological closure; basin41 is the mandatory hostile |
| Square-sum lifts | Designated edge roots; signless kernel | Solve root-correction compatibility and then full interval coverage; sparse embedded paths alone are insufficient |
| Finite-state approximation | Policy LP and residue phase clock | Test stabilization of optimal policies; finite precision cannot assert rationality or finite-state attainment |

Every lane changed another: sharp word counts gave a precise fractal cover;
the explicit motif exposed closure's lost odd factors; square unused edges
suggested the legality test that direct policy patching failed; ancestor
phase and resets then repaired that gluing in the pairing model. The
generalization suggested by this last proof requires a contracting parent
scale, a uniform reset bound, and uniform policy discrepancy. Those
hypotheses have not been established for longer-horizon pairing models or
for ordinary Collatz, so no blanket transfer is claimed.

## 7. Reproduction and status boundary

Run `python -B 04-computation/experiments/crossroads_family_20260926_audit.py`.
It checks all four lane programs in normal and optimized Python against
their retained outputs. NumPy is required only for guarded deep integer
arrays in the flow lane. The [audit record](crossroads_family_20260926_audit.md)
gives universes, hostile controls, and independent proof reviews. The
manifest hashes committed Git blobs, including text normalization.

These results supply exact constructions, a solved optimization boundary
in a modified model, and several falsifiable next targets. They do not
settle ordinary Collatz, LRC(14), or an unproved global square-sum claim.
