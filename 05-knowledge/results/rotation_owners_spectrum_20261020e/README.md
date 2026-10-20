# Regular endcap detours, a deficiency-three recurrence, and Fibonacci matchings

Research continuation, October 20, 2026. Integrated base: `7c6b843800`.
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
