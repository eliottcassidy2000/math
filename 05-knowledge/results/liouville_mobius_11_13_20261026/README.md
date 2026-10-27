# Liouville and Möbius in tournament composition

2026-10-26. Liouville records prime multiplicity modulo two; Möbius records
the same sign only when no prime repeats. Their relation is an exact square
transform, not interchangeable notions of arithmetic sign. This distinction
connects the earlier H=11/13 metagraph example to the 81/27 symmetry example
and gives a new sharp symmetry threshold.

| Result | Status and location |
|---|---|
| `lambda = mu * 1_squares`, `mu = lambda * (mu supported on squares)`; the corresponding Dirichlet, Lambert and ordinary Euler-product identities | **CITED classical identities, PROVED derivations**; [arithmetic A1–A3](ARITHMETIC.md) |
| Under all multiplication contexts, lambda has two response classes; mu distinguishes every squarefree integer, while all nonsquarefree integers have zero response | **PROVED**; [A4](ARITHMETIC.md#a4-what-multiplication-contexts-recover) |
| On r squarefree prime coordinates, the mu composition matrix has integral inverse and rank `2^r`; the lambda matrix has rank one | **PROVED**, exact controls r=1..5; [A4](ARITHMETIC.md#a4-what-multiplication-contexts-recover) |
| H=11/13 are the same two five-tournament classes already in THM-466; four cyclic triples, but one versus two directed five-cycles | **FINITE-EXACT recovery**, inherited mechanism; [recovery](RECOVERY.md) |
| If one realization has squarefree automorphism order, the p-part is invariant in its C3 fibre below **5p**, sharp for every odd prime p | **PROVED** from cited C3 realization theory; [S1](SYMMETRY.md#s1-the-squarefree-threshold-is-5p) |
| `lambda(|Aut|)` first varies within a C3 fibre at 12 vertices; `mu(|Aut|)` first varies at **15**, in one hypergraph type with 1024 labelled lifts and four classes | **PROVED**, tiny classification and displayed counts **FINITE-EXACT**; [S2–S3](SYMMETRY.md#s2-the-first-mobius-failure-is-fifteen) |
| Replacing mu by lambda in primitive-period inversion adds all square-divisor period layers; at periods 121 and 169 the false primitive count has remainder 2 modulo the period | **PROVED classical convolution consequence**, exact controls; [A5](ARITHMETIC.md#a5-primitive-periods-and-the-frobenius-connection) |

For the actual H=11/13 tournaments, order-join makes H multiplicative.
Both arithmetic signs are initially -1. After adjoining the H=11 class,
Möbius gives `mu(121)=0` and `mu(143)=1`, whereas Liouville gives +1 to
both. Every order-join context leaves their Liouville responses equal.
Conversely, `lambda(81)=1` and `lambda(27)=-1`, while Möbius vanishes on
both orders and on all their multiples. The observers detect different
information; neither refines the other.

The new fifteen-vertex fibre replaces the nonrigid four-block `C3 -> 1`
by the rigid five-block `D=S4 -> 1`, where S4 is the unique strong
four-tournament. Its four ordinary types are
`C3[D,D,D]`, `C3[D,D,D*]`, `C3[D,D*,D*]`, `C3[D*,D*,D*]`, with automorphism
orders **3,1,1,3**. All share the same 131 labelled cyclic triples and
Hamilton-path count **152377413**. The common C3 hypergraph determines
path-cover counts but loses the oriented child colours.

## Inheritance and scope

The closest proved mechanism is the decorated strong-module tree and its
colour-stabilizer formula in [the earlier classification packet](../tournament_composition_classification_20261026/README.md).
The canonical hostile is its 81/27 fibre. The corrected near miss is
automorphism-order invariance below 12, not at every order. The least-used
recovered sidecar is the [Euler-product ghost sequence](../../../07-reflections/euler-product-ghosts-and-zeta-side-cancellation-gates.md).
The [recovery map](RECOVERY.md) separates H-values, vertex counts, periods,
prime coordinates and residue classes when 11 or 13 occurs.

Anchor: arithmetic observers of tournament composition. Niche: divisor
inversion and Euler ghosts. Wildcard: primitive periods and the first rigid
symmetry witness. The live concept board is prime multiplicity / divisor
inversion / composition context / directed cycles and ports / primitive
period. [RECOVERY.md](RECOVERY.md) records how the calculations connect them.

These results classify specified observer responses and the first symmetry
failures. They do not classify all prime tournaments or make the modular
tree closed under Hamilton-path complement moves. The actual metagraph
still needs its path and port data. LRC(14) remains **OPEN**.

## Reproduction

From the repository root:

```text
python3 -B 05-knowledge/results/liouville_mobius_11_13_20261026/verify.py --replay
python3 agents/check_docs.py
```

[AUDIT.md](AUDIT.md) states the finite universes and proof checks. Normal
and optimized runs must match frozen receipts. The independent implementation
uses a sieve, colour-refined permutation enumeration and a backward path
recurrence; it imports none of the producer. This is a same-agent proof
review with independent computational paths, not an independent-agent audit.
