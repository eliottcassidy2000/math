# Odd primes, symmetry descent, and the boundary of natural bijections

20 October 2026. A focused continuation on HYP-9340, prime-circulant constructions, and HYP-9172.

## Results at a glance

- **HYP-9172(a): PROVED for every n≥5**, in both the tournament and switching-class formulations. A uniform elementary obstruction uses only the primes 3 and 5. No additive-prime conjecture or group classification is needed.
- **HYP-9172(b), numerical assertion: PROVED.** The odd and level Burnside branches are integers for every n≥3. The n=2 half-integer exception is sharp. A structural decorated pairing remains open; the naive pairing on all even-automorphism types is refuted.
- **HYP-9340: extended beyond p-groups in a precise direction.** An odd pronormal base, followed by one-prime symmetry refinement in its normalizer quotient, yields mixed-prime equalities and simultaneous preservation of the associated extension predicates. The unrestricted odd-group conjecture remains refuted.
- **Prime circulants: a sharp naturality distinction.** A generator-free projection works for p≡3 mod4. For p≡1 mod4, p>5, an abstract multiplicative-group automorphism obstructs every natural type injection. This is not an obstruction to every field-based algorithm.
- **Prime powers and definability:** even-degree finite-field extensions give a genuine Frobenius obstruction in the stated multiplier-type model. A separate circuit argument excludes a first-order definition of the particular parity projection over prime fields, without excluding every other bijection.

Historical novelty beyond the repository is not asserted. The papers and proofs below state their definitions and dependencies explicitly.

## 1. The strongest closure: a uniform 3/5 obstruction

Read [Uniform equivariance obstruction](../uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md), with its [independent audit](../uniform_equivariance_obstruction_20261020/review/AUDIT.md).

The question is whether a map on labelled objects can respect every vertex relabelling and induce a bijection of isomorphism types. The two cases are:

1. tournaments → sign-even graphs;
2. tournament switching classes → graphs that are both sign-even and Eulerian.

“Sign-even” means trivial edge-reversal sign on the graph automorphism group. It is different from having even degrees.

The proof has three short steps.

**A conservation law becomes rigidity.** The signed Walsh basis is invertible over Z[1/2]. Its weighted orthogonality implies that the product of source stabilizer orders and the product of target stabilizer orders have equal valuation at every odd prime. Equivariance gives stabilizer containment at every matched pair. Every odd-prime valuation can only increase, but its total cannot increase, so it must be preserved at every pair.

**Local 3-cycles force a global block system.** Source stabilizers cannot contain a double transposition when n≥5. For switching classes this follows by taking the unique representative with a fixed vertex as a source. Two single 3-cycles with overlapping, unequal supports generate a forbidden double transposition. Consequently all their supports are disjoint triples. A full ambient Sylow-3 subgroup supplies floor(n/3) such triples, exhausting all but at most two vertices.

**The required 5-part no longer fits.** A stabilizer preserving that triple partition has 5-adic order at most v₅(floor(n/3)!), strictly less than v₅(n!) for n≥5. But the empty graph's preimage would have to preserve the full ambient 3-part and 5-part simultaneously. Contradiction.

This also gives an explicit Hall-deficient family: source types carrying a full ambient Sylow-3 or a full ambient Sylow-5 subgroup. The two source subsets are disjoint, while the corresponding target subsets intersect at least in the empty graph. Their separate cardinalities agree, so their union is too small to receive distinct images of all these source types.

The [prime-power growth addendum](../uniform_equivariance_obstruction_20261020/PRIME_POWER_HALL_GROWTH.md) makes the obstruction quantitative. Put n=3^a5^b−1 and M=min(3^a,5^b), with a,b≥1. The M graphs K_(k,n−k), 0≤k<M, are distinct sign-even targets with both full ambient Sylows: binomial(n,k) is prime to 3 and 5 because each factor (n+1−j)/j has cancelling valuations. Even k give (M+1)/2 Eulerian targets. Thus the maximum matching misses at least M tournament types, and at least (M+1)/2 switching-class types. Balancing the two prime powers gives an explicit infinite subsequence with deficits bounded below by a constant times √n. This is a lower bound along that subsequence, not an exact matching count or an all-n growth rate.

The threshold is sharp: both containment-matching problems have solutions through n=4. The theorem does not prohibit arbitrary bijections of unlabelled types. In particular, the [Royle–Praeger–Glasby–Freedman–Devillers paper](https://arxiv.org/pdf/2204.01947), v2, asks for a natural bijection of unlabelled types without imposing this precise equivariant-map definition. The present closure is of HYP-9172(a)'s explicit formulation, not every possible interpretation of that broader request.

## 2. HYP-9340: what survives beyond one prime

Read [Pronormal bases and quotient refinement](../pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md), with its [independent audit](../pronormal_quotient_bijections_20261020/review/REVIEW.md).

The six-vertex C3 case is already covered by the inherited prime-power theorem. A direct independent certificate now makes the two embeddings explicit:

| Symmetry on six vertices | Tournament types | Sign-even graph types |
|---|---:|---:|
| One 3-cycle, three fixed vertices | 11 | 11 |
| Two synchronized 3-cycles | 3 | 3 |
| Independent C3×C3 | 1 | 1 |

The first two predicates intersect in exactly the third type. Normalizer-only counting introduces extra odd graphs and extra source presentations; these must not be mistaken for full isomorphism counts.

A subgroup K is pronormal if K and each ambient conjugate are conjugate inside the subgroup they generate. For odd pronormal K, two crucial facts hold:

- Full types admitting K are exactly normalizer orbits of K-fixed objects.
- Sign-evenness tested in the normalizer is equivalent to full sign-evenness.

The K-fixed tournaments form a Boolean cube. Its Fourier dual is exactly the K-invariant graph cube because all K-orbit sizes on pairs are odd. For any odd prime q, a q-unit Walsh matching in N(K)/K preserves the entire quotient Sylow-q action. It therefore preserves every symmetry predicate for

    K ≤ H ≤ N(K), with H/K a q-group,

simultaneously on the sector admitting K. These groups H can have mixed prime divisors and need not themselves be pronormal.

Examples include four distinct C7-based symmetries at n=10, with equal counts 12, 4, 2 and 2. The diagonal F21 action on 7+3 vertices is nonpronormal, yet lies in the theorem's scope. At n=15, C5 wr C3 has exactly one type on each side.

The old C15 boundary now has a structural explanation. At n=8 and 9 its C5 base is an ambient Sylow subgroup, so the quotient theorem proves equality. At n=10 that base is no longer pronormal: the two disjoint C5 factors cannot be conjugated inside their abelian join. The K5,5 part-swap sign can escape the chosen normalizer, producing the known 32-versus-31 counterexample. Abstract normality of C5 inside C15 does not repair this ambient failure.

## 3. A generator-free rule and two different obstructions

Read [Generator-free circulant naturality](../generator_free_circulants_20261020/GENERATOR_FREE_NATURALITY.md), with its [independent audit](../generator_free_circulants_20261020/review/INDEPENDENT_AUDIT.md).

For p≡3 mod4, the nonzero squares H form a group of odd order. Record a tournament's skew connection set on H as a binary word s. Define

    w(h) = s(h) + sum_(u in H) s(u), over F₂.

Extend w symmetrically to ±H. This is the complementary augmentation projector to odd-order averaging. It is generator-free, exactly two-to-one on labelled connection sets, and bijective on circulant isomorphism types.

For an abstract cyclic group of order 2m with m even and m≥4, the exponent automorphism i↦(m+1)i fixes every symmetric graph connection set. It moves a tournament type: the skew word 0^m1^m has two cyclic transitions, while its image has 2m−2. Hence no type injection can commute with all abstract cyclic-group automorphisms. The one-type m=2 case, corresponding to p=5, is exceptional.

This obstruction forgets field addition. It therefore does not rule out a generator chosen canonically using arithmetic labels in a prime field, or every parameter-free ring construction there.

Over F_(p^r) with r even, a different obstruction is genuinely field-theoretic. A Frobenius-fixed tournament multiplier type would yield a semilinear tournament automorphism of even order, impossible. The empty target type is Frobenius-fixed. Thus a Galois-equivariant multiplier-type bijection cannot exist. At q=9 the fixed-type counts are 0 versus 2.

The prime-power target here is the even-primitive-root necklace model, not automatically the fully sign-even graph class. The q=9 rook graph explicitly demonstrates the distinction. Additive finite fields are also not cyclic rings Z/p^r.

Finally, the displayed generator-free projection computes parity of an unbounded number of input bits. The cited classical AC0 parity lower bound excludes defining this particular map by one first-order ring-plus-input formula across all p≡3 mod4. Naturality and first-order definability are therefore genuinely different requirements. Definability of some other prime-field bijection remains open.

## 4. The branch integrality problem and its correct remaining boundary

Read [Branch integrality](../switching_branch_integrality_20261020/BRANCH_INTEGRALITY.md), [pairing boundary](../switching_branch_integrality_20261020/PAIRING_BOUNDARY.md), and [independent audit](../switching_branch_integrality_20261020/review/REVIEW.md).

The numerical proof combines two complementary local facts:

1. On each stabilizer H, Frobenius's root-count theorem makes the proportion of odd-order elements a dyadic rational. Summing over switching-class types puts N_odd in Z[1/2].
2. For every odd partition except (1,1), the fixed-class cube has enough powers of two to cancel the 2-part of its centralizer denominator. Every partition summand is therefore 2-integral.

The intersection of dyadic rationals and rationals with odd denominator is Z. Thus both branches are integers for every n≥3; n=2 gives 1/2+1/2.

This proves numerical integrality, including every positive order divisible by four. It does not produce the requested structural pairing. Indeed, the naive extension of the old matching-flip fails already at n=4. At n=8 there are 27 even-automorphism switching-class types, so no fixed-point-free involution can pair that entire unadorned carrier; N_lev=15 is not half of 27. Any successful pairing explanation must change or decorate the carrier.

The general reusable principle works at any prime p: the p-regular restricted Burnside average lies in Z[1/p]. A separate p-integral conjugacy-class formula forces it to be an integer. Here p=2 and the extra input is an explicit affine fixed-space dimension bound.

## 5. Mechanisms and next questions

The common algebraic engine is dyadic Fourier linearization. The Walsh transform is invertible after dividing by powers of two, so every odd prime supports its own symmetry-preserving matching. Those matchings need not fit together as one relabelling-equivariant map. Odd-order averaging also works in characteristic two, while even-order field automorphisms produce genuine obstructions. The branch-integrality argument completes the same local picture from the opposite side, by controlling the remaining power-of-two denominator.

This produces several precise mechanisms, rather than a universal symmetry claim:

- Monotone stabilizer containment plus a conserved product forces pointwise equality.
- Odd orbit sizes make the restricted Fourier pairing nondegenerate.
- Pronormality controls both fusion of subgroup copies and escape of an odd sign.
- Odd-order averaging supplies a canonical projection; the dyadic factor obstructs its naive continuation.
- Frobenius detects descent obstructions unavailable over prime fields.
- Complementary local denominator bounds turn a rational Burnside contribution into an integer.

The next precise questions are: weaken pronormality to an object-relative fusion/sign condition; characterize first-order-definable prime-field bijections other than parity projection; formulate genuine prime-power circulants with full sign-evenness; and find a decorated combinatorial carrier for the level branch.

One concrete next prime-power test is the affine group F₂₅⁺ ⋊ C₃, where the order-three multiplier acts irreducibly over F₅. Translations reduce unordered pairs to nonzero differences modulo sign; the order-three multiplier and negation have four orbits on the 24 nonzero differences. Thus there are only 16 invariant tournaments and 16 invariant graphs before quotienting. Its regular translation group is not an ambient Sylow subgroup of S₂₅. This makes it a small fixed-cube test of global fusion/sign descent beyond the current ambient-Sylow route; no global refined count is asserted here.

All proof packets include exact scope, explicit counterexamples to tempting stronger claims, and reproducible finite controls. Publication uses ordinary commits directly to main.
