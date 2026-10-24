# A sharp 18-clock floor for every existing fixed-exterior triple crossing

Status: **PROVED with a finite-exact original native exclusion.** This small
supplement explains the lowest value in the3015 boundary histograms. It does not
modify the frozen first audit.

Use the original sixteen cheap sources A from the native19-slot receiver, and
triple T={1457,1841,3015}. The actual triple label is1263 and its switching cost is
2+3+3=8 clocks. Every A source switches at exactly2 clocks.

## Claim and proof

Every existing crossing from the UUU branch to the VVV branch at the SAME receiver
grade and charge, with exterior records fixed, changes at least8 slots and costs
at least18 clocks. The bound is attained at grade8 and charge1.

Write u,v for the selected A-subsets in the UUU and VVV branches. Equality of
receiver populations gives |u|-|v|=3. Equality of receiver charges gives
lambda(u)-lambda(v)=1263. If h=|u symmetric-difference v|, then h is odd and h>=3.

If h=3, necessarily v is a subset of u and u-v is a three-element A-subset of
charge1263. Independent enumeration of ALL 560 three-element subsets shows that
none has this charge. Hence h>=5. The crossing changes h+3>=8 actual slots and
costs2h+8>=18 clocks. This proves the lower bound for every crossing that exists,
at any grade; it makes no assertion that a crossing exists at every grade.

For attainment, in increasing A-source order use base masks11079 and2826. The
first has eight selected A sources and triple UUU. The second has five selected A
sources and triple VVV. Both complete receiver states have V-population8 and
charge1. Their base disagreement is5, so they differ on8 native slots and cost18
clocks exactly. `receipt.json` retains their actual full paths and target records.

This is a sharp boundary-edge floor, not a global connectivity threshold for all
charge fibres: internal branch connectivity and the costs of reaching a minimizing
boundary witness are separate questions. If exterior records are allowed to
change during the crossing, the fixed-grade/fixed-receiver-charge argument no
longer applies without further conditioning.

## Authentication and reproducibility

The checker uses three byte-pinned original inputs: core full records, electrical
potentials, and canonical tournament nodes. It reconstructs A by literal path
costs and recomputes every label. It lists all560 tested triples and their charges.
It independently exhausts all19receiver guarded Hamiltonian menus against actual
canonical U/V target classes, checks the selected kernel representatives, inverse
partners, converse paths and38pairwise-disjoint receiver alternatives, and directly
verifies preservation of every original X/Y source/target row for the attaining
pair, including target diagonal terms and paired weight two.

The global82-slot resource supply and all-dyadic saturation are inherited from
the frozen first native audit, whose manifest hash is
`c77ee9067796f528f9bb45d3c619d9cdad21f19bcbd70c79a055296ad284dca6`.
Its write-up and manifest are byte-preserved reference archives, not substitutes
for re-running that original full audit. No new source, port, kernel capacity,
indexed-owner equality, all-order supply or microscopic transition contract is
introduced.

Run `python3 check.py`, `python3 -O check.py`, and
`python3 verify_manifest.py --relocate`. The last command makes a genuine fresh
manifest-only copy and executes both checker modes there. All checks survive -O.
