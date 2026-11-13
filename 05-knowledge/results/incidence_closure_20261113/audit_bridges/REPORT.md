# Independent audit of the completion, sheaf, and holonomy bridges

**Audit verdict: PASS for the exact reviewed bytes below.** No producer file was
edited. This is an independent mathematical and implementation review, not a
status promotion: B1–B6 remain **RESERVED** in the producer note until the root
integrates the review. The finite controls below are **FINITE-EXACT**. The cited
classical results and inherited local stress calculation are not new claims of
this packet. No theorem IDs or global files were changed.

Reviewed in the isolated `incidence-closure-20261113` worktree at base commit
`a8b6d8dbf0c6b601813785efb0b5904b46c23d43`. Producer files are uncommitted packet
artifacts at review time; the SHA256 values bind their actual contents. This
report does not certify any subsequent edit.

| Reviewed path, relative to this packet | SHA256 |
|---|---|
| `bridges/PROOF.md` | `e0e6c6428614e9ebdd0be5417d5373f9b509b43c3d57eb0081809d6aabd860ce` |
| `bridges/check_bridges.py` | `32eea5653748102011bded91edcd1c846ebd8aeac7b190253e2cd6e33436759d` |
| `bridges/check.normal.json` | `46bc599c3848ae69e04b71248c2bc96c4d474daa5fca54c391dfa6eb34b06d15` |
| `bridges/check.optimized.json` | `46bc599c3848ae69e04b71248c2bc96c4d474daa5fca54c391dfa6eb34b06d15` |

## Claim-by-claim review

**B1 — PASS.** For each connected fibre, edge equality is necessary and
sufficient for a scalar observer to descend; surjectivity of the quotient map
is what makes its descended value unique. For a connected graph with a chosen
root and bijective inverse edge transports, evaluation identifies sections with
the common fixed set of the based loop action. Backtracking cancellation is
legitimate because inverses are specified. No relations from unfilled triangles
are imposed. The proof also handles empty state sets; a chosen root presupposes
a nonempty graph. A nonidentity holonomy may still fix some states, so the
stronger conclusion “all root states extend” requires every loop to act
identically. Noninvertible edge maps need a diagram limit instead.

The independent triangle control enumerates every assignment of three
permutations of a three-element state set: among 216 assignments, 72 have no
section, 108 have one, and 36 have three. This tests the distinction between
nontrivial holonomy and an empty fixed set, beyond the regular two-state torsor.

**B2 — PASS.** The actual six raw color choices with fixed edge `03` have the
four displayed tree completions. The other two choices are disconnected
triangles plus an isolated vertex. Changing at most the two free coordinates
makes the completion graph exactly `K4`; the retained edge never moves. The
transport bits are expressly extra data, not a canonical invariant of those
trees. For any finite graph, the vertex-to-edge coboundary has kernel the
componentwise constants, hence rank `v-c`, quotient dimension `e-v+c`, and
exactly `2^c` solutions for a balanced assignment. Loops and multiple edges,
if allowed, simply contribute the corresponding cycle constraints; the exact
finite checker uses simple graphs only.

For `K4`, the eight balanced assignments have two torsor sections each; the
other 56 have none. This must not be restated for linear parallel sections:
a vector space always contains its zero section. Over a field of characteristic
not two, the associated sign line has dimension one in the balanced case and
zero otherwise. In characteristic two, every scalar sign becomes the identity.
The independent rank and enumeration controls retain these distinctions.

The signed-adjacency comparison is correctly restricted to guarded ratios for
a tested base permutation. The zero-support, absolute-value, and diagonal
conditions precede the coboundary test. A raw zero or nonunit signed matrix
entry is not an invertible transport. Nothing identifies the completion graph
with the tournament signed-support graph.

**B3 — PASS.** The bit-labelled graph cover is a genuine two-sheet covering of
the one-dimensional realization. Continuous local sections form a sheaf because
compatible continuous functions glue, and the right-inverse condition is local.
A section over an edge cannot change its discrete sheet in the interior, so
its endpoints satisfy exactly B2. Conversely compatible endpoint choices glue
the chosen lifted edges continuously. Thus a nontrivial cover can have an
already-sheafified section object and no global section. Effective gluing
applies to a matching family that exists; it does not produce one by fiat.
`F_s^+ = F_s` is understood through the canonical isomorphism, not literal
identity of implementation data structures.

**B4 — PASS, inherited mechanism.** On a genuine site, the first plus is
separated and the second is a sheaf. The empty covering of the empty open is
essential here. In the two-point model the original two labels on the empty
open force equality of the singleton labels in the first pass. That pass
replaces the empty value by a singleton, permitting four pairs in the second.
Normalizing the empty value in advance yields the four pairs in one pass.
For a one-point space one pass already suffices; the claim is an upper bound,
not exact two-step depth for every space or presheaf. The generic deterministic
observer refinement really has `m-1` strict rounds, so it supplies no universal
two-pass identification with plus. This finite example is already covered by
THM-4892(B), file `01-canon/theorems/THM-4892-tournament-sites-plus-construction.md`.

**B5 — PASS after the producer's pre-audit forgetful-map repair.** A trivial
finite-dimensional tensor factor gives a finite direct sum of the same linear
compatibility equations, so the displayed natural map on parallel sections is
an isomorphism. A scalar `-1` loop retains zero fixed-vector dimension after
such amplification in characteristic not two. The zero-dimensional tensor
factor is harmless; no assertion that a nonzero vector survives is made.

This is a statement about the specified flat transport, not its underlying
bundle. The complex section `exp(pi*i*t)` satisfies the endpoint identification
and is nowhere zero; hence the sign-holonomy complex line is topologically
trivial. The eight real rotation blocks similarly give an explicit frame of
the sixteen-dimensional underlying real bundle. Neither trivialization is
parallel. The reviewed text now expressly excludes any inference that a
nontrivial KO/KU class is preserved by these operations. Classical Bott and
Clifford periodicity are confined to their own categories and require an
actual module/transport map before transferring any predicate.

**B6 — PASS.** Every homomorphism to an abelian group kills the commutator
subgroup of a perfect normal subgroup. Hence a sign character factors through
the perfect quotient; surjectivity ensures a nontrivial character stays
nontrivial. Quillen's map has exactly that quotient on fundamental groups, so
this is a precise obstruction-preservation statement, independent of the
homology part of the plus theorem. It does not identify Quillen plus with
presheaf plus. The independent `S5/A5` control verifies normality, perfection
(by literal commutators), and the surviving two-element sign quotient.

## Inherited scope and primary-source check

The fixed-tail connectivity and diameter are correctly inherited from
[polynomial staircase diameter, Theorem 3 and Corollary 4](../../graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md).
The hypotheses include an acyclic selection of every upper-tail difference
color, and moves can change up to two free coordinates. The bound is
`r²-3r+3` for `r>=3`; the inherited small cases are separate. This does not
preserve the abstract tree type and does not prove the Graceful Tree Conjecture.

The local stress formulas agree with
[the inherited stress derivation](../../five_stress_connection_atlas_20261112/STRESS.md).
They assume homogeneous isotropic plane elasticity, no body force, a specified
holomorphic coordinate with nonzero derivative, and the stated treatment of
isotropic stress points. With `phi''=-F`, `psi'=-G`, continuation of `phi'`
first gives a rotation increment; after it vanishes, continuation of `phi` and
`psi` gives a separate translation increment. The conjugation of the contour
factor `2*pi*i` is essential: the residue translation condition is
`kappa Res(zF) - conjugate(Res G)=0`. For `Z=z+1/z`, the independent Laurent
calculation confirms the displayed residues and dimensions five for stress,
two for global displacement at `kappa=3`. Additional holes at zeros of `Z'`
contribute no residues because these rational functions are regular there.
This checks the narrow inherited calculation, not a broader elasticity theorem.

The following primary passages were independently reopened on 2026-11-13:

- [Stacks, Theorem 7.10.10, tag 00WB](https://stacks.math.columbia.edu/tag/00WB)
  explicitly gives separatedness, one further plus for a separated presheaf,
  and the general double-plus conclusion. These are cited classical inputs.
- [Bott, *The Stable Homotopy of the Classical Groups*, pp. 314–315](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf)
  gives the stable unitary period two and orthogonal period eight. It does not
  classify parallel sections of a prescribed connection.
- [Atiyah–Bott–Shapiro, *Clifford Modules*, §§2, 4, 6](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf)
  uses negative-square real generators; §4 gives the real matrix amplification
  with period eight, and Proposition 6.8 gives the degree-eight isomorphism of
  the specified module quotient. The grading convention is retained.
- [Quillen, *Higher K-theory for categories with exact sequences*, §1, pp. 95–96](https://aareyanmanzoor.github.io/assets/articles/quillen-higher-K-exact-seq.pdf)
  states the perfect normal kernel and acyclic-map construction, with local
  coefficient systems taken on the target. B6 needs only its group quotient.
- [Sáenz's publisher abstract](https://link.springer.com/article/10.1007/s00707-008-0118-3)
  confirms a broader simply connected 2–5 parameter result and links an erratum.
  Its full proof, Neményi's 1933 full article, and the erratum were not reviewed
  here. No claim about their unexamined contents is imported.

The original correction lineage matters: THM-4892 and its audit already retract
an unsupported Hamiltonian-path restriction/global-section claim. This bridge
uses actual restriction maps for continuous cover sections and finite
presheaves; it does not revive that claim.

## Implementation and independent controls

The producer's `sieves` routine enumerates every downward-closed sieve and
checks covering by union, including both covering sieves of the empty object.
`matching` enforces naturality along every subset inclusion. The colimit joins
families along every covering refinement, and restriction pulls a sieve back
by intersection with the subobject's downset. Representative independence and
presheaf identities/composition are checked. The stage tests correctly handle
the one-point exception and normalized-empty cases. Comparing the number of
blocks suffices to stop the refinement loop because each signature contains
the previous color, so a step cannot merge old blocks. No removable `assert`
statement supplies a correctness gate.

Both normal and optimized producer runs were independently replayed and match
the frozen receipts byte for byte, each reporting 5,432 gates. That replay is
execution verification, not independent mathematics. The separate
[check_independent.py](check_independent.py) imports no producer code and tests:

- every simple graph on zero through four labelled vertices: 76 graphs and
  761 transport assignments, with direct gauge enumeration versus binary rank;
- the actual six graceful choices, and all 216 three-state triangle transports;
- all 64 `K4` sign systems over `F3`, exact scalar-section enumeration versus
  rank, tensor dimensions zero/two, characteristic-two degeneration, and the
  rational rank-sixteen `-I` obstruction;
- plus via the independently justified **terminal covering sieve**, rather than
  the producer's all-cover colimit, on 16 constant-presheaf cases through four
  points and 482 further two-point presheaves with arbitrary empty-overlap maps;
- the literal `S5/A5` quotient and 81 exact Gaussian-integer stress residue cases.

The 3,504 independent gates pass in normal and optimized Python. The receipts
are byte-identical. Finite tests support their declared universes; the general
claims rely on the audited proofs. These controls do not audit a large
incidence census, prove a new sheafification theorem, or construct a canonical
transport on graceful completions.

```sh
python3 05-knowledge/results/incidence_closure_20261113/audit_bridges/check_independent.py > 05-knowledge/results/incidence_closure_20261113/audit_bridges/normal.json
python3 -O 05-knowledge/results/incidence_closure_20261113/audit_bridges/check_independent.py > 05-knowledge/results/incidence_closure_20261113/audit_bridges/optimized.json
```

| Independent artifact | SHA256 |
|---|---|
| `check_independent.py` | `e852ea253345ac03fd9f8e5e913bce4381c64c2bf8739c255e5f4aa9dff0ba5f` |
| `normal.json` | `f2e6983a07f9a61e5a1e3b7e90592343b2a6c036c0d92731d50b92907313a3a1` |
| `optimized.json` | `f2e6983a07f9a61e5a1e3b7e90592343b2a6c036c0d92731d50b92907313a3a1` |

For cross-lane reproducibility, the signed-lifting note's separately requested
stochastic-basis clarification was frozen at SHA256
`1710c4730bb2512e8574984ddb4a5a398fe46c699824571bd309b40210407ac7`.
It distinguishes the natural fibre-constant stochastic block from its
orthonormal conjugate. That correction changes neither B1–B6 nor their producer
checker and is not a status promotion of the signed-lifting lane.

## Status-only promotion rebind, 2026-11-13

The preceding report records the prepromotion review snapshot. Root has now
promoted B1–B6 to **PROVED IN NOTES**. This audit independently reversed all
eight declared bridge edits in `../PROMOTIONS.json`, including correction of
the stale RESERVED sentence at the end. The inverse recovers exactly the
reviewed SHA256 `e0e6c6428614e9ebdd0be5417d5373f9b509b43c3d57eb0081809d6aabd860ce`.
No mathematical text changed. The final promoted bridge proof is bound by
`271ad8a5a9ffe408e62b2148b88d26ce9fd255b595267523740da0811923fb2c`.
The prior statements in this report that the note remained RESERVED describe
the original audit stage, not its current promoted status.

The exact manifest used has SHA256
`dbccf7ced4fe064dbe908d196ebdf507542018b2c8acf74ff57fdcc84db1a898`.
The independent inverse-verification record is
[promotion_verification.json](promotion_verification.json), SHA256
`d69e60772899562bb8c53f3ea2f5789077d375952e36d3cd493fea6097ac703c`.
It also verifies the two root-note promotion inverses.

The original receipts are preserved byte for byte as
`normal.prepromotion.json` and `optimized.prepromotion.json`, both retaining
SHA256 `f2e6983a07f9a61e5a1e3b7e90592343b2a6c036c0d92731d50b92907313a3a1`.
The earlier audit-artifact table refers to that preserved snapshot. Only the
expected bridge-proof hash was rebound in `check_independent.py`; its checks
and mathematical controls did not change. Normal and optimized reruns again
pass all 3,504 gates and agree byte for byte. After deleting the two provenance
fields `reviewed_sha256` and `source_sha256`, their JSON payload is identical to
the archived prepromotion result.

| Current rebound artifact | SHA256 |
|---|---|
| `check_independent.py` | `bbb151d8970345ca26e2095397138c7fb3151525495890be6f57a6163f262972` |
| `normal.json` | `5539744400a9606d5d56b6dc0d20bf35fbe57d163015025a1d4063d6363c33b2` |
| `optimized.json` | `5539744400a9606d5d56b6dc0d20bf35fbe57d163015025a1d4063d6363c33b2` |
