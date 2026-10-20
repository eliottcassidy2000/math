# Regular endcap detours, a deficiency-three recurrence, and Fibonacci matchings

Research continuation, October 20, 2026. Initial base: `7c6b843800`;
incoming work subsequently integrated through `7441e88d7c`.
The results below distinguish all-order proofs, cited inputs, finite controls,
and the global problems that remain open. No literature novelty is claimed.

## Main advances

**PROVED, with a CITED input and small FINITE-EXACT boundary certificates:**
every blue line with regular outer strong endcaps of order at least five
has a black-only bypass of length at most three, retaining its marked target
and arbitrary middle. [Regular endcaps](regular.md) proves the required
Hamiltonian connectivity for all regular orders at least seven using the
cut bound, Thomassen's theorem, and independently complete endpoint checks
on the three order-seven and fifteen order-nine classes. The pentagon uses
the inherited explicit identity. [Rotation presentations](black.md) gives
a separate elementary three-move construction for every marked circulant
endcap, including the threshold connection sets where block extraction fails.
Its proper-interval rotation lemma handles composite orders too.

**PROVED:** the entire feedback-deficiency-three owner layer has a bounded
external-port kernel sum. Writing N3(k) for its capped weighted-core count,

    N3(k)/(2k-3)!! = P(k-2),   k>=2,
    deg P=7,   leading coefficient8/3.

Thus its eighth forward difference vanishes, giving an exact polynomial-
coefficient recurrence. [Owners](owners.md) derives this from kernels with
at most twelve endpoints, retaining every external matching partner and
the exact clone-pair filter. Arbitrarily many remaining matching edges are
summed, not enumerated. The finite coefficient polynomials C2,C3 are computed;
C4,C5,C6 are specified but uncomputed. This is a whole-stratum recurrence,
not an unrestricted scalar recurrence for A000568.

**PROVED elementary:** the ordered SCC language on `{1,C3,R4}` is in explicit
bijection with matchings on two parity paths. [Arithmetic](arithmetic.md)
identifies the odd-order counts as

    F_m F_(m+1) = 1,2,6,15,40,104,273,... = sum_(j=1)^m F_j^2.

At order13, retaining the second edge of the even path partitions these
273 actual tournament types into63 containing it and210 omitting it.
Triangle count is the matching size. Hamiltonian-path weight needs the full
block decoding: the R4 block has weight5, whereas two independent C3 blocks
have weight9. The arithmetic half-plane and braid/Farey comparisons likewise
retain the operation address and central integer height that their quotients
would otherwise erase.

**OPEN:** the ordinary Hamiltonian-path spectrum equals all positive odd
integers except7 and21. [Spectrum](spectrum.md) pins the primary-source status,
proves the two omissions with a small cutoff, and derives the reduced-witness
bound `n<=floor(3 log_3 h)`. Degree-preserving directed-cycle flips connect
each labelled score fibre, but the regular seven-vertex reversal175->171
refutes automatic two-unit interpolation. These cycle flips are different
from Hamiltonian-path-complement metagraph lines.

## Incoming proofs that meet at the cutoff

The incoming [effective black-defect bound](../effective_black_defect_order_bound_20261020/PROOF.txt)
gives `n<=32` for a strong tournament with zero black path orbits. The
inherited odd-circulant reduction makes `n<=31`; the independent certificate
bank covers every skew connection set through31. [The composed proof and
replay](purity_closure.md) therefore confirm the complete strong pure-blue
alphabet `{1,C3,C5}` and the all-order class series, including the empty object,

    P(z)=(1+z+z^3+z^5)/(1-z^2-z^6).

This is an independent proof route, not the first classification: the
repository already contains an audited
[uniform prefix-rerouting proof](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md).
The new connection removes the residual33..131 range from the separate
finite-atoms packet by composing its certificates with the stronger incoming
bound. It gives an exact recurrence for the pure-blue part of the metagraph,
not for all tournament classes.

The incoming [relocation theorem](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt)
also proves that all labelled Hamiltonian words of any tournament are connected
under legal single-vertex relocation. These moves remain inside one class
fibre; they are not black complement lines. Consequently, if L_n has the
`2^binom(n-1,2)` fixed-path tilings as vertices and legal normalized
relocations as edges, then `A000568(n)=number of connected components of L_n`
for n>=3. Every tournament is represented by Redei's Hamiltonian-path theorem;
relocation cannot change its class, and the new theorem connects its whole
fibre. Adding reflection instead counts converse-merged classes. This is an
exact structural generation rule without Burnside, but generating that entire
cube is still exponential in n squared; it is not the sought small recurrence.
The incoming
[bridge-hub and double-descent theorem](../bridge_hubs_and_double_descent_20261020/BOUNDARY_HUB_AND_DOUBLE_DESCENT.md)
rules out nonleaf paired-black bridges at its rigidly guarded sources and
forces any terminal-free bridge branch through a self-converse single-descent
bottleneck. It leaves bridges between terminal basins and general connectivity
open. Together these results narrow the unresolved objects without conflating
three different graphs.

## Inheritance and concept board

The closest proved mechanisms are the incoming
[SCC-basin endpoint criterion](../scc_basins_20261020d/moves.md), the
[fixed-deficiency local owner](../deficiency_three_20261020/ownership.md),
and [THM-4094 / two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
The canonical hostiles are endpoint-changing R4 presentations, local gaps
with different outside partners, and a response gcd of two with holes.
The corrected near misses are the false replacement fit in MISTAKE-055
and treating the two proved Hamiltonian omissions as exhaustive. The
least-used sidecars are ordered endpoint matrices, intrinsic heavy-gap
addresses, and the central exponent in the braid/modular map.

| Lane | Live concept | Preserved data and present outcome |
|---|---|---|
| Anchor | Marked black detours | Exact labelled target, middle arcs and endpoint guards; all regular endcaps closed |
| Anchor | Deficiency-three ownership | Full exceptional ports and local owner order; unbounded matching bulk summed |
| Niche | Ordinary H spectrum | Path multiplicities and ordered endpoints; exclusions proved, completeness open |
| Niche | SCC matching language | Exact word inverse and triangle count; Fibonacci product and marked-edge partition |
| Wildcard | Sum/product half-plane | Both parents and fixed diagonal; output-only shadows lose construction data |
| Wildcard | Farey and braid lifts | Operation word and central exponent; matrix alone loses full twists |

Each comparison is given a source, target, explicit map, preserved predicate,
lost coordinate, and small hostile in its detailed note. The matching map is
an actual bijection. Its scalar H weighting is not independent between rows.
The half-plane and braid comparisons motivate retained state; they do not
prove a tournament or Collatz claim by analogy.

## What this leaves

1. **Connectivity:** regular outer endcaps can no longer supply an unresolved
   mixed-to-mixed blue bridge. The remaining general boundary includes
   nonregular outer SCCs without the endpoint certificate, and the global
   problem of joining the resulting black basins. A detour theorem alone
   is not a proof that all active nodes lie in one component.
2. **Counting:** evaluate C4,C5,C6 for a fully numeric D3 formula; then test
   whether ownership remains stable under analogous bulk deletion at higher
   deficiency. The present lexicographic proof uses `q<=1` in an essential
   step, so its unbounded-deficiency extension is not assumed. A000568 still
   requires all deficiencies and remains outside this scalar recurrence.
3. **Spectrum:** construct uniform strong witnesses in the two prime lanes
   of THM-4094, or connected exact-response states with enough overlapping
   two-unit intervals. Score-fibre connectivity by itself supplies neither.

Independent agents audited the black macros and the matching bijection;
the owner proof received a separate logical audit of its bulk, lexicographic,
large-light and leading-coefficient steps. Parent review checked the cutoffs,
count cancellation and stated boundaries. [Validation](VALIDATION.md) gives
the exact finite universes and reproduction. The small exact computations
support the proofs; no large-order census or fitted extrapolation is used.
