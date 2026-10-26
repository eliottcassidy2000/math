# Independent audit of the fibre and arithmetic packet

2026-10-26. **VERIFIED** proof audit; **FINITE-EXACT** reproduction checks.
The two collaborators read each other's proofs independently. This is an
audit of the stated arguments and finite controls, not a formal proof-assistant
certificate or a claim of historical novelty.

## Fibre proof audit by arithmetic_species

The independent reviewer checked [RESULT.md](RESULT.md), in particular:

- F2: if n<4p, every p-subgroup of the common size-preserving quotient group
  fixes the fewer than p children whose ordinary isomorphism colour can vary.
  It therefore lies in both refined colour stabilizers. The induction and
  sharp 4p construction are valid.
- F3: below order12 no odd valuation can vary. At12 the only possible bad
  prime node has three children of size4, and variation requires all three
  to lie in the unique four-vertex A/B fibre. The uniqueness conclusion is
  for the underlying unlabelled cyclic-triple hypergraph, exactly as stated.
- F4: the ordered-cover bijections preserve the vertex set of every path.
  Linear-child reorder and quotient dualization lift through an internal
  module because entering and exiting arcs are uniform on that module.
- The cited primary PDF's Theorem10, Theorem24, Proposition51 and Corollary57
  supply the common strong-module tree and independent quotient choices.
  No automorphism-order invariance is attributed to that literature.

One terminology issue was found and repaired: Aut(Q)/Stab(0) is a coset
space of size p; the stabilizer need not be normal. There was no change to
the proof or conclusion. The full hypergraph action and direct tournament
automorphism search in the checker are separate from the quotient formula.

The normal and optimized fibre receipts are byte-identical, SHA256:

    9e001620bdfa0f97c9f01ae30bd91c5bcbb13d893794d84e354d83824cf4d537

## Arithmetic proof audit by cyclic_fibre_audit

The independent reviewer read [arithmetic/RESULT.md](../arithmetic/RESULT.md)
and its checker in full. The following implications were checked separately:

- The automorphism formula is valid for the canonical invariant partition.
  Its image preserves complete oriented child types, and chosen child
  identifications give the claimed noncanonical splitting.
- The colour group acts freely on the ordered partitions with prescribed
  sizes. This proves multinomial divisibility and nonnegative local deficits.
- A fixed ambient Sylow subgroup forces uniform cross arcs between its
  independent transitive orbit factors, and recursively forces precisely
  the displayed circulant towers. Both directions of the saturation
  classification and the fixed-labelled-subgroup dimension are valid.
- Digit sum2 is exactly the obstruction to strong connectivity among the
  saturated presentations. Splitting one nontrivial digit block into p
  blocks loses exactly one p-factor and allows a strong quotient on p+1
  vertices. Digit sum1 and every digit sum at least3 attain the unrestricted
  bound. The restrictions n>=3 and odd p exclude the small exceptional
  nonexistence cases.
- The all-order maximum and strong maximum are separate extrema for each
  prime; neither proof asserts simultaneous saturation of different primes.

No theorem-changing issue was found. The full arithmetic checker was rerun
from the repository root:

```text
python -B 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/check_arithmetic.py --out 05-knowledge/results/tournament_composition_classification_20261026/fibre/arithmetic-crossaudit.json
```

It passed 160 signed-pair action cases, 8,930 exhaustively generated fixed
tournaments, the all-tournament maxima through order5, the strong attaining
constructions through order40 at four primes, the 400 carry identities,
and all eight explicit full-group computations. The independently produced
`arithmetic-crossaudit.json` is byte-identical to the author's
`arithmetic/controls-normal.json`, SHA256:

    88bb4789ddd095ef04df81d1afc831330454d2203e4637d9fbb3ad15d05806b7

The author also supplied an optimized-mode receipt. Those finite tests are
controls for the proofs, not evidence by extrapolation for the all-order
claims.
