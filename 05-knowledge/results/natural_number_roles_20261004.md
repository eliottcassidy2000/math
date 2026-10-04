# Natural-number roles across the repository

**Status: curated index of exact roles, inherited scoped results, and explicitly
typed proposed connections; not an exhaustive occurrence census.** Updated
2026-10-04 by `difference-families-oct04`.

This companion fills the small-integer gap deliberately left by the earlier
[9437-file constants atlas](constants_atlas_20260926_recurrent_numbers.md).
That atlas remains the broad lexical census. The present index asks what a
number **does**: dimension, norm, clock, carry, orbit count, obstruction, or
label. Repetition becomes a mathematical connection only after specifying a
map. Read the proof/status in the linked source before importing a claim.

The strongest new connected chain is

    2,3 -> 3^3-2^3=19 -> norm-19 Eisenstein ideal
        -> three 19-vertex triangular tori -> Paley residue action
        -> Collatz exponent mod 6 / retained exponent mod 18.

The field tower connects to the same chain through its missing order-19
phase and the factorization `Phi_18(2)=57=3*19`. The precise construction,
proof boundaries and executable checks are in
[Checked route switches and the missing 19-phase](checked_switch_phase19_20261004.md).

## 1. Index

“Map” means an actual construction or identity with its specified predicate.
“Echo” means separate exact facts whose cross-domain significance remains
unproved. “Proposal” names a construction to investigate, not a theorem.

| Number | Roles recovered | Relation and source |
|---|---|---|
| **1** | odd positive Collatz root; unit for multiplication; the diagonal excluded by distinct-summand pairs | Root and unit are different types. A certificate must name its terminal basin. [Sibling grammar](creative_sibling_20260925.md), [braid](arithmetic_braids_20260917_collatz.md). |
| **2** | binary valuation; first-reset boundary; two proper subfield degrees begin with 2; two golden ideal branches | Map: exponents control the real denominator and the binary guard. A reset equal to 2 creates the debt state (4), rather than the legal switch (3). [New construction](checked_switch_phase19_20261004.md). |
| **3** | Collatz multiplier; ternary inverse integrality; cubic extension; three torus layers at 19; lost ternary digit | Maps exist, but the spaces differ. `3/2 mod 19` has order 3; `F64* -> F4* x F8*` has kernel of order 3. [Phase geometry](checked_switch_phase19_20261004.md), [subfield norms](sextic_subfield_clock_20261004.md). |
| **4** | ordinary positive Collatz cycle includes 4; sibling map `4n+1`; four isomorphism classes of 4-tournaments; gap of both `(7,11)` and `(8,12)`; `76=4*19` | Sibling relation is a map. A count of four does not canonically identify tournaments with the four denominators `{1,2,11,76}`. [Difference operations](difference_family_operations_20261004.md), [golden phase recursion](difference_diagonals_golden_phase_recursion_20261004.md). |
| **5** | spherical `(2,3,5)`; golden discriminant; norm of `1+3phi` is -5; negative Collatz basin at -5 | Exact algebra gives `(1+3phi)^2=5phi^4`. Geometry and the signed basin have no supplied dynamics-preserving map. [Sixth clock](sixth_clock_branches_20261004.md), [geometry](camion_busch_gaps_polyhedra_collatz_20261002.md). |
| **6** | flat `(2,3,6)`; binary dimension of F64; order of 2 mod 9; six units mod 14; six free edges after fixing a Hamiltonian path on five vertices; exceptional Mersenne index | Strong maps: cyclotomic degree, inverse valuation clock, and the phase action table. Six free bits represent 64 encodings, not six tournament isomorphism classes. [Subfield clock](sextic_subfield_clock_20261004.md), [THM-873 μ6 spine](../../01-canon/theorems/THM-873-mu6-spine-tight-locus.md), [new note](checked_switch_phase19_20261004.md). |
| **7** | hyperbolic `(2,3,7)`; Paley K7 torus; `F8*` order; Lucas `L4`; minus basin pair 5,7; first uncertified source in the old sibling grammar | K7 -> torus is a map. The Lucas identity `phi^8+1=7phi^4` is Cayley–Hamilton. Source 7 needs a different certificate rule, not a geometrical label. [Camion–Busch](camion_busch_gaps_polyhedra_collatz_20261002.md), [unit-minus-one clocks](collatz_lucas_monotile_discrepancy_20260930.md), [sibling](creative_sibling_20260925.md). |
| **8** | F8 subfield; `2^3`; Sylow-7 normalizer index `168/21`; 8 and 12 have gap4 | The subfield sits in F64; the group index belongs to PSL2(7). A cardinality match alone does not identify the subfield with a homogeneous space. [Subfield clock](sextic_subfield_clock_20261004.md), [Ogg audit](ogg_triangular_chowla_20260927_audit.md). |
| **9** | order of the sextic root modulo 2; `Phi9(X)=X^6+X^3+1`; inverse residues mod 9; nine quadratic residues mod 19 | Map: the nine directions partition into three cube-root cosets. Root order 9 differs from generator order 63. [Subfield clock](sextic_subfield_clock_20261004.md), [phase geometry](checked_switch_phase19_20261004.md). |
| **10** | edges of K5; triangular number T4; proper nontrivial divisors of `p^2qr`; `10=7+3` | The divisor classification uses nested subsets: seven squarefree divisors include the three prime divisors. It is not a disjoint 7+3 partition. [Divisor balance](seam_prime_balance_20260925.md). |
| **11** | `L5`; 11 halvings in the -17 cycle; golden denominator11; prime triple `{2,3,11}`; binary coefficient label of `zeta^3+zeta+1` | The first three roles have their own exact clocks. The last role depends on a chosen polynomial basis and does not explain prime 11. [Golden recursion](difference_diagonals_golden_phase_recursion_20261004.md), [unit clocks](collatz_lucas_monotile_discrepancy_20260930.md), [new phase seed](checked_switch_phase19_20261004.md). |
| **12** | `tau(p^2qr)=3*2*2`; twelve classes of 5-tournaments; 12 and 8 share gap4 | The first 12 is an exponent-box size; the second is a quotient count. Their identification is an echo. [Divisor balance](seam_prime_balance_20260925.md), [difference operations](difference_family_operations_20261004.md). |
| **13** | `2^8-3^5`; LRC(14)'s thirteen nonstationary speeds; exponent13 gives a translation on the 19-phase; first bad reset from 7 lands at 13 | New map: `2^13=3 mod 19`, hence `C13(x)=x+13`. The Pillai gap13 is a separate integer equation. [Phase geometry](checked_switch_phase19_20261004.md), [constants atlas](constants_atlas_20260926_recurrent_numbers.md). |
| **14** | LRC anchor; K7 triangular faces; Heawood vertices; 14 edges in the finite flat-star refinement graph | K7 duality gives Heawood and its 14 vertices. The flat-star graph and the LRC problem require separate predicates. [Camion–Busch](camion_busch_gaps_polyhedra_collatz_20261002.md), [flat stars](forest_20260926_tilings.md). |
| **15** | T5 and K6 edge count; `6+9` enlarged triangular array; first square-sum Hamiltonian threshold in the audited lane | Array coordinates can store information; no Collatz certificate follows from their number. [Square-sum audit](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md), [decoder seams](seam_threshold_20260925.md). |
| **16** | AMM superblock uses16 terms; `2^4`; threshold for nonplanarity of the integer graph with +2 and doubling edges; denominator of -29/16 | These are separate mechanisms. The last is a quadratic-dynamical parameter, not a Collatz multiplier. [THM-4468](../../01-canon/theorems/THM-4468-golden-zero-superblocks-beat-golden-amm12592.md), [graph minors](decoder_minors_20260925.md), [difference norms](difference_visibility_norms_20261004.md). |
| **17** | minimum absolute value in one known negative Collatz cycle; 17 unordered locally flat polygon stars | Equal numbers, different predicates. The cycle has seven odd steps and eleven halvings; the star count is an Egyptian-fraction classification. [Golden carriers](collatz_golden_carriers_20261004.md), [flat stars](forest_20260926_tilings.md). |
| **18** | second cubic field dimension; multiplicative order of 2 modulo 19; golden period associated with denominator76 | First two are directly linked: a primitive order-19 root first lives in F_(2^18). The golden clock needs its separate denominator/phase map. [Depth tower](cyclotomic_depth_towers_20261004.md), [golden recursion](difference_diagonals_golden_phase_recursion_20261004.md). |
| **19** | missing multiplicative phase; `3^3-2^3`; Eisenstein norm19; torus vertices; `76/4`; golden norm `L3^2+3` | Strongest joined hub: explicit norm factorization, quotient lattice, and guarded phase semiconjugacy. The scalar golden norm does not alone identify golden ideals with Eisenstein ideals. [New proof](checked_switch_phase19_20261004.md), [depth tower](cyclotomic_depth_towers_20261004.md). |
| **20** | proper divisors of 60 include20; one nonextendable flat star is 5.5.10, another4.5.20 | Distinct boundary mechanisms. Flat local angle sums do not guarantee global tiling. [Flat stars](forest_20260926_tilings.md), [divisors](seam_prime_balance_20260925.md). |
| **21** | edges of K7; order of its Paley affine automorphism group; nonzero rank-one F4-tensor-F8 elements; cyclic flat-star types | Paley's group acts on its21 arcs. Tensor rank-one count is `(4-1)(8-1)=21`. No isomorphism is claimed between that cyclic multiplicative subgroup and the nonabelian Paley group. [Camion–Busch](camion_busch_gaps_polyhedra_collatz_20261002.md), [tensor clock](sextic_subfield_clock_20261004.md). |
| **24** | labelled transitive 4-tournaments; labelled strong 4-tournaments; Hamiltonian cycles of P7 in the inherited convention; tetrahedron flags | Tetrahedron flags map to transitive orientations; equality with strong orientations uses the cited duality. The cycle count is separate. [Polyhedral tournament study](camion_busch_gaps_polyhedra_collatz_20261002.md). |
| **27** | order ofzeta in dimension 18; `3^3`; least positive source whose reset-two debt in this note reaches31 | Field order and inverse integrality share powers of 3; the hard source27 remains an arithmetic orbit. [New construction](checked_switch_phase19_20261004.md). |
| **29** | `|N(phi^7-1)|`; numerator of quadratic parameter -29/16 | Exact norm identity explains that arithmetic occurrence. A shared parameter does not transfer the Collatz cycle graph. [Difference norms](difference_visibility_norms_20261004.md). |
| **31** | Mersenne `2^5-1`; smaller source in 63=>31; Paley31 with some connected genus63 triangular maps; debt state from 27 | The first two have the checked reset map. The geometry uses a different genus formula. [Route switch](checked_switch_phase19_20261004.md), [Camion–Busch](camion_busch_gaps_polyhedra_collatz_20261002.md). |
| **32** | F4 x F8 as an additive Cartesian product; number of children per clock under a binary precision lift; pure-power endpoint in the 3 mod 6 row | The tensor product has 64 elements, so Cartesian and tensor products must be distinguished. [Subfields](sextic_subfield_clock_20261004.md), [sixth clock](sixth_clock_branches_20261004.md), [braid](arithmetic_braids_20260917_collatz.md). |
| **38** | triangular faces in each19-vertex torus | Exact map: every triangular edge is incident to two faces, so `F=2E/3=38`. [Phase geometry](checked_switch_phase19_20261004.md). |
| **42** | inverse curvature denominator for(2,3,7); K7 torus full affine symmetry count in the inherited map; number of binary2x3 matrices of rank2; largest polygon in the flat-star census | Each has a counting mechanism. There is no automatic identification of all 42-element objects. [Geometry](camion_busch_gaps_polyhedra_collatz_20261002.md), [tensor clock](sextic_subfield_clock_20261004.md), [flat stars](forest_20260926_tilings.md). |
| **54** | period of 2 modulo 81; phase increment in the four-step inverse-ray family | Map: the guard for four odd steps lives modulo `3^4`; its unit clock has order 54. [Boundary/codec join](collatz_boundary_compiler_20261004.md). |
| **57** | `Phi18(2)=3*19`; edges of each19-vertex torus; order of the translation-plus-cube-root subgroup of its affine symmetries | The norm identity and torus edge count are both exact; the subgroup has a concrete action. [Phase geometry](checked_switch_phase19_20261004.md). |
| **60** | least mixed balanced divisor example `2^2*3*5`; spherical(2,3,5) rotation group order | The value60 is shared, but no transport between divisor predicates and rotations has been supplied. [Divisor balance](seam_prime_balance_20260925.md), [geometry audit](ogg_triangular_chowla_20260927_audit.md). |
| **63** | `2^6-1=3^2*7`, no new prime; F64 multiplicative generator order; `63=>31` checked Collatz switch; genus of connected triangularK31 maps | First two are exactly the same field-order formula. The route is an LTE/carry identity. The genus is `(31-3)(31-4)/12`, a separate mechanism. [Sixth clock](sixth_clock_branches_20261004.md), [new note](checked_switch_phase19_20261004.md), [Camion–Busch](camion_busch_gaps_polyhedra_collatz_20261002.md). |
| **64** | F64; all labelled orientations of K4; six-bit vector space | Choosing an edge order and field basis gives a bijection. It does not make field multiplication preserve tournament isomorphism or Hamiltonian structure. [Tensor clock](sextic_subfield_clock_20261004.md), [4-tournament geometry](camion_busch_gaps_polyhedra_collatz_20261002.md). |
| **73** | `Phi9(2)=73`; factor of the dimension 18 generator order 13797 | It occurs in the half-field factor `2^9-1=7*73`, unlike the missing 19 from `2^9+1=27*19`. [Depth tower](cyclotomic_depth_towers_20261004.md). |
| **76** | `4*19`; `L9`; golden scalar denominator for the -17 cycle; `phi^18-1=76phi^9` | Exact Lucas/golden arithmetic. The finite gate at denominator76 has 240 golden cycles; integer realization selects a much smaller subset. [Golden recursion](difference_diagonals_golden_phase_recursion_20261004.md). |
| **81** | `3^4`; coefficient of a four-odd-step word; modulus whose unit period is 54 | The same guarded affine map forces all three roles. [Boundary compiler](collatz_boundary_compiler_20261004.md). |
| **100** | entries of a10-by-10 multiplication table; `10^2` | Useful proposed array, but no map from the table to universal Collatz certificates is supplied. Reconstructing a symmetric table needs boundary/axis conventions. |
| **128** | coarse-cylinder period for word(1,1,2,3); period of the exact79=>39 switch family | Both are dyadic guards, with different words. `7+128k` has the coarse endpoint oddpart(5+81k); `79+128k` has exact endpoint101+162k. [Boundary compiler](collatz_boundary_compiler_20261004.md), [new switch](checked_switch_phase19_20261004.md). |
| **139** | `3^7-2^11`; cycle clock for -17; factor of 1807 | Ordered cycle carry is 2363=`17*139`. The factor in Sylvester's sequence does not by itself predict the cycle. [Constants atlas](constants_atlas_20260926_recurrent_numbers.md), [unit clocks](collatz_lucas_monotile_discrepancy_20260930.md). |
| **168** | order of PSL2(7); Klein quartic automorphisms in the geometry lane | Actual map: Paley affine group of order 21 is a Sylow-7 normalizer, index 8. It is not the full symmetry group of the Paley tournament. [Ogg audit](ogg_triangular_chowla_20260927_audit.md). |
| **171** | edges of K19; order of its translation-times-square affine subgroup | A direct action is sharply transitive on directed Paley arcs; both counts are `19*9`. The three torus layers each contain57 edges. [Phase geometry](checked_switch_phase19_20261004.md). |
| **183** | `Phi6(14)=14^2-14+1=3*61`; LRC deep-well denominator; Eisenstein norm | Same norm form as57, with argument14 instead of 8. It links exact arithmetic form, not an established Collatz-to-LRC theorem. [THM-873 μ6 spine](../../01-canon/theorems/THM-873-mu6-spine-tight-locus.md), [constants index](../../00-navigation/CONSTANTS-INDEX.md). |
| **189** | Hamiltonian paths in P7; denominator of a strict LRC reserve; factor of 378 | Paley path count transports through its actual Sylvester-core identification. The LRC reserve needs its own inequality. [Constants atlas](constants_atlas_20260926_recurrent_numbers.md), [LRC constants](../../00-navigation/CONSTANTS-INDEX.md). |
| **323** | order of `<2,3>/<2>` modulo 87211; `17*19` | New exact map: odd-step count mod 323 is necessary to recognize unit-slope macros at the next binary phase prime. It is not a cycle-length restriction or a -17 basin label. [Next-depth probe](checked_switch_phase19_20261004.md). |
| **511,513** | `2^9-1=7*73`, `2^9+1=27*19` | Exact split of dimension 18 group order. The cubic-root-plus-one element uses the first factor and the 3-primary part of the second, leaving19. [Depth tower](cyclotomic_depth_towers_20261004.md). |
| **5779** | next golden norm after 19 in `L_(3^(k-1))^2+3` | The exact recurrence is `N_next=N^3-3N^2+3`. It differs from the binary next primitive factor 87211, so the towers must not be identified termwise. [Depth tower](cyclotomic_depth_towers_20261004.md). |
| **13797,262143** | deficient and full multiplicative orders in dimension 18 | New explicit phase restoration changes the generator from order 13797 to 262143 by adjoining an order 19 factor. [Exact field check](checked_switch_phase19_20261004.md). |
| **87211** | next binary primitive factor; next Paley phase prime with 14535 torus layers | Its branch layer clock is 18 and Frobenius layer period 9; the missing unit-slope register has 323 states. The single-branch translation seen at 19 disappears. [Next-depth probe](checked_switch_phase19_20261004.md). |

## 2. Four underlying mechanisms worth keeping separate

**Incoming work incorporated:** the concurrent
[5–6–7 surface and observation study](geometry_567_collatz_observers_20261004.md)
adds two useful cross-field routes. The same Heawood graph has a torus map
with seven hexagons or a genus-three map with three fourteen-gons, depending
on its local rotations. Separately, triangle-to-star flow transport has
kernel `A^7` and each nowhere-zero star has `|A|-3` nowhere-zero triangle
lifts. These give actual connections to coloring and the flow problems;
they require rotation, circulation and homology coordinates. The general
Paley orientation law in that work is inherited here; our new layer and
within-layer coordinates refine its two-parity-bit observation quotient.

The next incoming checkpoint adds
[actual completed routes in every 19-adic residue](mod19_route_lifts_20261004.md)
of sharp uniform odd depth four, and
[virtual contraction families](virtual_contraction_ladders_20261004.md).
The latter shares the source cylinder `79+128t` with our reset switch but
provides a different smaller dependency. This makes the retained set of
certificate alternatives a concrete interface between the two constructions.

**Unit-minus-one / cokernel clocks.** The older
[Lucas and monotile lane](collatz_lucas_monotile_discrepancy_20260930.md)
recognized Mersenne orders, golden torsion orders and Collatz cycle clocks
as instances of norm/resultant/cokernel calculations. The robust common
object is a linear or affine return map plus its **carry class**. The size of
a cokernel gives possible phase addresses; the actual carry chooses whether
a specified word realizes an integer cycle. This is much stronger than
matching numbers and much weaker than classifying all Collatz cycles.
Its elementary identities can be used independently of that note's unreviewed
claims about an attached monotile manuscript.

**Quotients losing one ternary digit.** The F64 pair of subfield norms has
kernel3. The three Paley19 layers only see a Collatz exponent modulo 6;
exponents1,7,13 require the next ternary digit to distinguish affine orders
3,3,19. These are two precise information-loss mechanisms. There is not yet
a natural isomorphism identifying their ternary digits. Proposal: use a
uniform record type for “quotient state plus lift digit,” then demand a
commuting diagram for each intended operation.

**Moving proof boundaries.** AMM12592 superblocks and the present route
compiler both improve a method by retaining information across a local
boundary. The AMM proof permits signed polynomial deficits; the Collatz
proof uses two real guarded trajectories and a smaller certified source.
They share a proof-design move, not a conserved quantity. The old
`log_2(3)` numerical lead is SUPERSEDED by
[THM-4494](../../01-canon/theorems/THM-4494-amm12592-exact-ratio-bottom-regime-c-below-log2-3.md).

**Local geometry versus global realization.** A flat polygon star need not
tile; a pair of Steiner systems need not have a single vertex-link cycle;
a compatible finite valuation word need not give a terminating infinite
Collatz orbit. The common research move is to retain the gluing coordinate
and test global realization. The respective predicates and obstructions
remain different. LRC(14) receives a useful source of phase coordinates
from its AP's six equality times, but no new LRC bound is claimed here.

## 3. The sixth root, coprimality, and a disciplined speculative lane

The requested scalar

    lambda=(5*phi^4)^(-1/6)=(1+3*phi)^(-1/3)

has an exact algebraic role. Writing `t=1/lambda` gives
`t^6-5t^3-5=0`; reducing that polynomial modulo 2 leads to the sextic field
and its63-clock. This is an algebraic map, not evidence that the scalar
predicts a physical Higgs mass. Separately,

    6/pi^2 = product over primes p of (1-1/p^2)

is the coprimality density. Both appearances of 6 invite a comparison of
local factors and global information, but no equation in this work connects
that density to the sixth-root scalar. [Earlier norm/visibility work](difference_visibility_norms_20261004.md)
is the place to retain that distinction.

One concrete speculative object is an **integer with phase and proof
coordinates**: an actual integer label, its word guards, ternary inverse
addresses, alternative common-future certificates, and compatible cyclic
phases. At deeper field levels the missing factors provide additional
coordinates. These phases may make a discarded relation visible, as19
has done here; universal termination still requires a well-founded rule for
the certificate dependencies. Finite phase cycles are expected and cannot
be treated as failures of the integer Collatz conjecture.

Arithmetic correction retained from the earlier prompts: `63 != 2*19+6`
(the right side is 44). Useful true neighbors are `63=3^2*7=3*19+6`,
`57=3*19`, and `76=4*19`. The new norm factorization explains why57 and 19
belong together; changing an incorrect equality is more productive than
building an analogy upon it.

## 4. Next tests selected from the index

1. **Anchor:** classify the 239 residual seed obligations of the current
   finite compiler by the reset-two debt boundary. Try to replace a whole
   family by one exact decreasing join; retain all alternative certificates.
2. **Niche:** the first next-depth probe is now complete: at 87211 the
   six-phase law becomes eighteen-phase and the one-step translation
   disappears; the quotient of order 323 repairs the missing word address.
   Test later prime factors with the same subgroup-quotient procedure.
3. **Wildcard:** align the cokernel carry, field phase and torus lift digit
   as typed coordinates. A successful alignment must commute with the
   guarded branch and distinguish the 1,7,13 exponent example; a cardinality
   match or an unguarded modulo 6 rule does not pass that test.

These are active routes with concrete hostile examples. Neither a complete
positive Collatz proof nor a proof that the known negative cycles exhaust
all negative basins is asserted.
