# Focused review of the third-wave interfaces

**Verdict: PASS for the new proved statements reviewed below. No material
mathematical correction required.** This is a focused internal review, not an
external-person audit, a literature-priority review, or a replay of the large
shell censuses. Inherited results were accepted at their stated scopes.

Reviewed sources:

- `../bridges/DOMINATION_TRANSFER.md`
- `../constructors/TWO_EDGE_REPAIR.md`
- `../arithmetic/QUASIPERFECT_NORMS.md`

## 1. Domination adapter

The tournament side is exact: deletion of the apex is transitive, so every
cycle uses the apex; the allowed odd cycles correspond to even-cardinality
increasing subsequences beginning at A and ending at B. The exceptional
adjacent-endpoint coefficient is correctly kept as2. All longer endpoint
coefficients are2^(distance−1) after the OCF weight.

The four-way domination recurrence is valid, including the boundary cases:
with the first A selected and the last B omitted, choosing the last selected
B gives the powers2^(j−2), while no B selected forces precisely the full A
shore. With both endpoints omitted, the removed A-only/B-only sets exist
exactly under the two stated boundary-letter conditions. The middle word is
nonempty for n≥3, so the two excluded shore sets cannot coincide. Isolated
leading B and trailing A vertices may be removed without changing the count.

This proves equality of counts. The note correctly does not claim that its
inductive proof identifies individual dominating sets with Hamilton paths.
The word/interleaving is retained as the required constructive sidecar.

The all-order29 exclusion is sound. The extreme A/B pair forces
1+2^(n−2)≤29, hence core length≤6. The31 binary cores of lengths2 through6
are the complete reduced universe. Prime29 cannot appear by multiplying
counts of several chain-graph components unless one component already has
count29. The same logic applies to13 with its smaller core bound. The given
tree's only nontrivial automorphism swaps two leaves and reverses no edge;
the count29/evenness example therefore lies outside the adapter as claimed.

The two-register update uses the full binary prefix x, rather than h alone;
the parity correction is necessary and correct. The final arbitrary-bipartite
inequality is explicitly only finite evidence/open, not promoted to a theorem.

## 2. All-apex gap and rank-template repair

The adjacent01→10 sorting move increases the count strictly: the local pair
adds2 and all other affected weights move to a weakly larger distance. At a
fixed number of ones, the unique sorted word is therefore the maximum.
Every unsorted sorting path has the same unique penultimate boundary-swap
word, proving the asserted second level at fixed weight. This argument does
not assume that nonadjacent arbitrary flips are monotone.

For N=2a+1, the two balanced weights a,a+1 have maximum n_a+3. Their best
unsorted words lose3·2^(a−1)−2. The next sorted weights lose3·2^(a−1),
which is larger, and all further sorted weights are smaller. Thus the global
top-two-value claim and all-a absence of n_a+1 among apex reorientations are
valid, including a=2.

The five-fixed-vertex/two-repeatable-module adjacency used for the four-edge
template agrees with the actual flips. In particular the flipped Alast/Blast
arc and the three apex reversals are all retained in the quotient. The P/S
symbols form a DAG, so the compressed-word language is genuinely finite;
it does not silently truncate unbounded graph paths. Surjection/Stirling
expansion gives the displayed integer exponential formula. The a=3 empty-P
boundary is treated separately instead of deleting its zero-base terms.

The subtraction from the target and the positive lower bound for a≥5 are
correct. The a=4 hit is consequently isolated within this specified rank
pattern. The note carefully distinguishes complete two-edge shells, a chosen
three-edge branch, and fourth-edge continuations of selected near misses;
no complete distance-four exclusion or all-a two-edge exclusion is claimed.

## 3. Quasiperfect/norm interface

The graph identity has a direct path bijection. When the singleton z is
interior, a Hamilton path has exactly the ordered block form
A_left B_left z A_right B_right, with B_left and A_right nonempty. The two
independent choices are Q2(A)−H(A) and Q2(B)−H(B). Placing z at either end
contributes the two H(A)H(B) terms. Thus the diagonal quadratic form is an
actual tournament construction, not just an arithmetic resemblance.

The (n+1)H injection into ordered two-path covers is valid and retains both
the original path and its cut. Its size bound supplies a complete finite
universe for each candidate norm profile. The automorphism divisibility
statement uses ordered vertex lists, so the group action is free as asserted.

The classical ring step is sufficiently supplied in the note itself:
nearest-integer approximation in Z[sqrt(−2)] gives remainder norm at most
3/4 of the divisor norm. Norm-Euclideanity therefore gives unique
factorization; for prime norm P, every norm-P element is associated to one of
the two conjugate factors. Units are±1. The resulting positive-coordinate
uniqueness, combined with x≥ny, proves the restricted prime-target barrier.
It does not obstruct a hypothetical quasiperfect source, whose divisor sum
is composite as explained. Odd-square quasiperfect necessity remains an
explicit inherited assumption, not a new existence or nonexistence result.

The factor-three switch requires the physical Q2 coordinate. Its ratio bound
reduces every H≥3 case to n≤5, while H=1 is handled separately by transitivity.
Thus the small labelled census is a justified complete base, not an
experimental cutoff. The displayed three profiles and targets3,51,99 have
the stated scope. The finite diagonal-target compiler and the Pell first
failure at Q2=42 also preserve the graph/arithmetic distinction correctly.

## 4. Cheap decisive controls

`third_wave_check.py` is a new independent small control, importing none of
the author scripts. It checks:

- all31 chain cores through length6 by direct dominating-subset enumeration,
  including the boundary recurrence and exact finite image;
- all32 apex cuts with five old vertices by literal Hamilton-path permutation
  counting, giving top levels23,19;
- the rank-template quotient's edge interface for a=3,4,5,6;
- independent subset-DP counts9,465,3639 for a=3,4,5.

All pass; see `third_wave_receipt.json`. No edge-shell census or inherited
mathematical theorem was rerun. The norm arguments were checked directly at
the proof/interface level, with the author's new small certificate inspected
rather than duplicating its complete through-order-five enumeration.

Reproduce from this directory:

    python3 third_wave_check.py
