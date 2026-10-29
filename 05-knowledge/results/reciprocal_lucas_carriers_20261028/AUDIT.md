# Scope, evidence and reproduction audit

2026-10-28. Same-agent proof and alternative-computation review. No independent
agent audit, Lean formalization or external priority review is claimed.

## Proof obligations

| Claim | Mechanism and quantifiers | Failure boundary checked |
|---|---|---|
| P1 | Factor (u-v)(uv-1); nondegenerate reciprocal pairs in characteristic !=2; tree reconstruction modulo common inversion. | A retained relation such as v=u^2 removes choices. An oriented root observation removes the common gauge. |
| P2 | Exact order in the split or inert quadratic algebra; phi/psi has norm one. | Primes 2,5 separate; splitting is not exact order; the F11 involutions have different fixed sets. |
| P2b | Future increments distinguish all phases; joint periods combine by LCM. | Named vector, not its OR; no claim of primality equivalence. |
| P2c | First-hit indices have odd quotients; positive-index gcd follows with valuations. | The admissible index set is not increment-closed; eventual silence does not prove an infinitude of composites. |
| P3 | Native forbidden-neighbor Gram splits into singleton and larger-leaf blocks. | Minimal polynomial degree is not rank or a general circuit lower bound. |
| P4 | Empty independent set is intrinsically fixed; distance is weight; converse sign is edge parity. | No implication from even automorphism order to sign-oddness, or from a positive embedding to a positive sum. |
| P5 | Closed transports act through the actual order-three centralizer. | Exact moved ports need three values; unobserved phases may be quotiented, named ones may not. |
| P6 | Freely chosen group frames modulo a common left action. | Stabilizers and non-flat natural connections are not covered by the free flat-frame model. |
| M1–M3 | Integer trace recurrences; explicit conjugates; quadratic order discriminants; g18's degree obstruction. | No inference from surd length, class number alone, or equal class counts to a complete modular carrier. |

The Lucas graph source orientation is intrinsic. No residue arrows with
ties have been forced into a tournament. The tournament controls use the
actual directed adjacency and actual permutation transporters. The word
"minimal" is always relative to the stated observer and operation menu.

## Exact finite universes and independent routes

[check.py](check.py) uses the Python standard library and explicit `need`
checks that remain active under optimization. [receipt.json](receipt.json)
contains the full outputs, not only summary counts.

- Reciprocal fibers: every one of the 2^k sign assignments for k=1..8 in
  path and heap charts; 49 rational pair identities; all vertex-sign cut
  codewords on complete graphs of sizes 2..6; four Gaussian allocations at 65.
- Prime clocks: every prime at most 1,000 other than 2,5 and every index
  0..300. Three realizations agree: recurrence residues, projective matrix
  powers, and arithmetic in the full quadratic algebra. There are 166
  primes and 460 divisibility incidences. Order exactness is checked by
  every prime divisor of the claimed order. All factors of L_2..L_60 are
  found by complete trial division, with no probable-prime assumption.
- Valuations: first-hit r<=40, multipliers 1,3,5,7,9,l; the 2-adic formula
  through index 120; all 7,260 unordered positive-index gcd pairs through 120.
  All distinct admissible indices through 300 are checked for coprimality.
- Joint banks: full output periods for the four declared banks, plus every
  divisor of the period as a candidate smaller period. The separate mathematical
  proof, not only the enumeration, supplies the minimum at arbitrary bank size.
- Nine-vertex holonomy: all 81 displayed automorphisms, 27 anti-involutions,
  and 729 based triangles; the inherited isomorphism oracle separately
  recovers the full group. Frequencies are 297,216,216. The actual vertex
  orbits have sizes 1,1,1,3,3.
- Stars: construct the native flip Gram for m=1..9 and check its polynomial
  on every coordinate vector. Modular elimination supplies an independent
  rank lower bound for m<=6. The block proof supplies the upper bound and
  the exact spectrum at every size. Toggle graphs are fully constructed
  for n=3..13 to check vertex/edge counts and the intrinsic empty vertex.
- CM: enumerate all reduced primitive positive forms at the twelve
  discriminants in the table and at -72. The ambiguous-class criterion
  identifies the doubled class groups; numerical surd comparison is not
  used to determine a class group.

[replay_inherited.py](replay_inherited.py) runs the earlier actual atom controls
without rewriting their receipts. Runtime sources include modules loaded
outside `sys.modules`; these are disclosed and pinned too. This is an inherited
implementation replay, not an independent reimplementation of 81/27.
It retains attachment orbits 176/512, transporters 27/0, and allowed normal-basis
shifts {0,2,4} at F64 and {0} at F256. The latter is the declared tournament
orientation guard, not a claim that F256 has trivial Galois group.

## Numerical versus exact evidence

The q-series coefficients are computed exactly by truncated integer formal
series. Decimal calculations use 120-digit working precision, a Machin
formula for pi and 80 product factors. Four special values are compared
independently with their defining products to within 10^-100. These are
**VERIFIED numerical comparisons**, not certified interval evaluations or
proofs of Ramanujan's evaluations. The special values are **CITED** from
his primary paper. The integer trace targets and the g18 degree obstruction
have separate exact proofs. No rational minimal polynomial is inferred
merely from matching a decimal value.

## Literature and recovered provenance

- [Ramanujan, Modular equations and approximations to pi (1914), sections 1,4–6](https://ramanujan.sirinudi.org/Volumes/published/ram06.html):
  the defining products, displayed evaluations and near-integer argument;
  [original PDF](https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf).
- [Schertz, Weber's class invariants revisited (2002)](https://numdam.org/item/JTNB_2002__14_1_325_0.pdf):
  normalization/congruence/eta-multiplier requirements, particularly Theorems
  1–2 and the multiplier discussion. No unverified transcription of a missing
  OCR formula is imported.
- [Bilu–Hanrot–Voutier, Existence of primitive divisors of Lucas and Lehmer numbers, Theorem 1.4](https://www.researchgate.net/profile/Paul_Voutier/publication/2397582_Existence_of_Primitive_Divisors_of_Lucas_and_Lehmer_Numbers/links/559999c308ae99aa62cc6983/Existence-of-Primitive-Divisors-of-Lucas-and-Lehmer-Numbers.pdf?origin=publication_detail):
  used for the U-sequence F_(2n), giving a new divisor of L_n for n>15.
  It is not a prime-value theorem for the sequence L_n.
- [Matsusaka, The Fourier coefficients of the McKay–Thompson series and the traces of CM values (2017)](https://link.springer.com/article/10.1007/s40993-017-0090-x):
  the introduction's normalized level-two eta formula. Our M1a spells out
  the different nome convention when evaluating Ramanujan's G and g.

Classical trace, finite-torus, gauge and gcd mechanisms are explicitly
attributed as such; their assembly and controls are the local contribution.
The [connection map](CONNECTIONS.md) gives the selected historical sources,
lost coordinates and decisive tests. Search matches and older log claims
are not used as proof authority.

## Repairs made during recovery

Three old outputs contained demonstrated errors and now carry prominent
historical correction notices, with lineage in `01-canon/MISTAKES.md`:

1. `fibonacci_lfunctions.out`: a table labelled F_p was F_(p+1), and the
   statement p|F_p for every prime is false already at p=3. The correct
   congruence is F_p=(5/p) modulo odd p. The prime 2 is inert in Z[phi].
2. `dihedral_fibonacci_bridge.out`: forgetting orientation turns the union
   of a tournament and its converse into K_p, with automorphism group S_p.
   The displayed dihedral action concerns an oriented pair, not that union.
3. `monstrous_moonshine_s134.out`: Ramanujan tau is attached to the weight-12
   form Delta, not a claim that the weight-zero j function has those Hecke
   eigenvalues. A transitive tournament has one Hamilton path, not n!.

Other speculative statements in those historical outputs are not promoted.
The present g18 and 130 hostiles are new warnings in this packet, not
allegations that Ramanujan claimed the false generalizations they refute.

## Reproduction and sealing

From the repository root:

```powershell
python -X utf8 -B 05-knowledge/results/reciprocal_lucas_carriers_20261028/check.py
python -X utf8 -B 05-knowledge/results/reciprocal_lucas_carriers_20261028/replay_inherited.py
python -X utf8 -B 05-knowledge/results/reciprocal_lucas_carriers_20261028/verify.py
python -X utf8 -B -O 05-knowledge/results/reciprocal_lucas_carriers_20261028/verify.py
```

The verifier checks pinned sources and packet hashes, recomputes both
receipts in memory, compares complete mathematical payloads, and checks local
document links. It does not reseal or overwrite saved evidence. The explicit
maintainer `seal.py` is separate. Normal and optimized reports differ only
in their optimization flag. `recover.py` can refresh discovery separately;
its results naturally change with the live repository and are not mathematical
test expectations. No theorem namespace is reserved or promoted by this packet.
