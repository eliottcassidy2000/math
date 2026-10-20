# Exact correction candidates for the inherited prime-circulant claims

These are narrow corrections, supported by the new proof and exact checker. They do not retract the type-count identity or the existence of an affine-equivariant forward map inducing a type bijection.

## 1. Replace a labelled bijection with an orbit-bijective equivariant map

Affected current files:

- `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`: frontmatter title/status and Corollary C.
- `05-knowledge/results/even_graphs_20261019_sylow_refinement_components_walsh_duality.md`: status table C and section 3.
- `05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md`: frontmatter (b), and section (b), item 2.

Replace claims of an `AGL(1,p)`-equivariant **bijection between labelled circulant tournaments and labelled even circulant graphs** with:

> An `AGL(1,p)`-equivariant map from labelled circulant tournaments onto labelled sign-even circulant graphs exists and induces a bijection on isomorphism types.

Where “bijection” was intended to refer only to types, state that explicitly. The existence result through THM-4531 Lemma N1 is correct in precisely this form, because that lemma uses stabilizer containment, not stabilizer equality. The new derivative-path algorithm constructs such a map uniformly for every odd prime.

Minimal witness: p=3 has two labelled translation-invariant tournaments and one labelled sign-even circulant graph (the empty graph). The affine stabilizers have orders 3 and 6, respectively. At p=5 the labelled counts are 4 and 1. A forward map can collapse an orbit while inducing a bijection on orbit sets. It cannot be a labelled bijection.

## 2. Narrow the maximal odd subgroup forcing argument

Affected claims in the same full note and HYP-9340 say that all tournaments invariant under K_max, the maximal odd multiplier subgroup, must map to the empty graph for every p. This is false without a dyadic restriction.

Write p−1=2^r q, with q odd. The K_max-fixed tournament sector has exactly

    2^(2^(r−1)) / 2^r

multiplier types. This equals one only for r=1 or r=2. At p=17, r=4 and q=1, K_max is trivial and the sector is all 16 tournament types. They cannot all map to the empty graph.

Safe replacement:

> If v2(p−1)≤2, the unique tournament type invariant under the maximal odd multiplier subgroup must correspond to the empty graph under an equivariant type bijection. For larger dyadic part the fixed sector can contain several types, and that forcing argument no longer applies.

The finite experimental statement that each listed fixed integer multiplier rule fails somewhere can remain if its original computation supports it. The faulty universal forcing explanation should not be retained as its proof.

## 3. Update the open constructive target without overclaiming

The all-prime derivative-path packet gives an explicit, bounded, invertible-on-types rule, equivariant in the forward labelled realization. Thus the corresponding “uniform rule for prime circulants” target can be marked solved in this precise algorithmic sense, after independent proof audit. Questions about generator-independent, first-order ring-definable, or fixed-multiplier formulas remain separate. No general all-order S_n-equivariant statement follows.

The new source is `prime_circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md` (publisher may choose a different final directory). This candidate path should be replaced by the verified final path before linking it.
