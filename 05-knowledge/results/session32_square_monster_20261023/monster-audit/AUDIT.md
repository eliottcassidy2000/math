# Session 32 independent audit: Monster prime partition and valuation cycle

**Verdict: PASS.** No mathematical blocker remains in the audited new results. This audit made no Git mutations and did not edit the author's sources. The author corrected one literal implementation comment and clarified two moduli-language points after feedback. All calculations below are exact.

## Audited scope and reproduction

Audited files:

- `../session32-monster/MONSTER_BABY_PRIME_PARTITION.md`
- `../session32-monster/audit_level_two.py`
- `../session32-monster/independent_field_audit.py`
- `../session32-monster/level_two_exact.json`
- `../session32-monster/independent_field_exact.json`
- `../session32-root/VALUATION_CYCLE_AND_QUOTIENT.md`
- `../session32-root/check_valuation_cycle.py`

`check_independent.py` is a new implementation, importing neither author script. It exhaustively enumerates roots of the Legendre Hasse polynomial in a directly implemented quadratic finite field, maps those roots to j-values, enumerates the u-coordinate fibers, then forms their Fricke quotient s-values. It does not assume at the outset that supersingular j-values are rational. Run:

```
python3 check_independent.py > independent_exact.json
```

Its only nonstandard input is the author's JSON, used for comparison after constructing the point sets. The independent script itself uses the Python standard library, not SymPy. Both author scripts were also replayed. Their outputs were byte-identical to the corresponding supplied JSON after the endpoint refinement. Environment: Python 3.12.14, SymPy 1.14.0. The only replay warning is a SymPy modular-integer ordering deprecation, not a failed assertion.

## 1. Group-order arithmetic: exact and complete

The numerical orders and factorizations agree with the [ATLAS Monster page](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/) and [ATLAS Baby Monster page](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/). No historical maximal-subgroup completeness assertion is imported from those pages.

Independent integer arithmetic verifies:

- M/B = 194,478,922,284,018,372,000.
- M/(2B) = 97,239,461,142,009,186,000.
- N = 39,046,339,332,000 = 2^5 3^7 5^3 7^4 11 13^2.
- 29·41·59·71 = 4,980,721.
- M/B = N·4,980,721.

The fifteen Monster primes partition exactly as claimed: equal positive valuations {17,19,23,31,47}; Monster-only support {29,41,59,71}; increased common valuations {2,3,5,7,11,13}. This is exhaustive for this pair of group orders, not an empirical pattern over an unspecified universe. The latter six happen also to be the primes whose Monster valuation exceeds one; that is not a general statement about arbitrary pairs of groups.

[Cummins–Matias, introduction p.2](https://arxiv.org/pdf/1710.01071) explicitly identifies 2.B as the centralizer of a Monster 2A involution. Thus the involution class has size M/(2B), not M/B. The displayed ratios are integers, not quotient groups. The Sylow statement in the note is valid: for each unchanged prime, a Sylow subgroup of the centralizer is also a Sylow subgroup of the Monster, by equality of p-valuations.

## 2. Cited theorem hypotheses and their legitimate use

[Aricheta, Theorems 1.1/3.1 and 1.3/3.3](https://arxiv.org/html/1809.07421) gives Q_p(X)/2 = g(X^p)−g(X) for p not dividing the level, where Q_p counts nonrational supersingular geometric points. For X=X0+(2), the correct target is X0(2p)/<w2,wp>, the full Atkin–Lehner quotient X0*(2p), not the Fricke-only quotient X0+(2p). Table 1/Appendix A identifies its rationality primes as 3,5,7,11,13,17,19,23,31,47. The centralizer criterion requires a Fricke element, not merely an arbitrary element of the indicated prime order. These are cited global results, not proved by the finite computations.

Here p=2 is excluded by the hypothesis. This does not assert that an abstract coarse projective line has no good model at 2. The author's p>=5 count excludes p=3 for its trace-zero test, whereas the theorem includes it. The new independent Hasse enumeration also verifies p=3 directly: lambda=2, j=0, u=2, s=1, with exactly one rational quotient point.

## 3. Exact Frobenius defects and endpoint refinement

The independent calculation agrees with all fifteen author rows: every Monster prime at least 5, plus 37 and 43. It additionally checks p=3.

| p | Supersingular j points | Supersingular u points | Supersingular s points | Nonrational s points | Unique quadratic factor |
|---:|---:|---:|---:|---:|---|
| 29 | 3 | 7 | 4 | 2 | s²−9s−13 |
| 41 | 4 | 10 | 5 | 2 | s²+8s−8 |
| 59 | 6 | 15 | 8 | 2 | s²−10s−6 |
| 71 | 7 | 18 | 10 | 2 | s²−22s−18 |

The four discriminants are 17,14,6,59 in the respective prime fields, each a nonsquare. Hence every row has exactly one nontrivial Frobenius orbit of length two. All other s-points in each row are rational. Applying the cited formula with g(X0+(2))=0 gives genus one for X0*(58), X0*(82), X0*(118), and X0*(142).

Both conjugate s-values have the same unordered pair of rational j-endpoints:

- p=29: {2,25}.
- p=41: {3,32}.
- p=59: {28,48}.
- p=71: {41,66}.

The code independently evaluates j(u) for both u-values in each nonrational s-fiber, asserts all four pairs, and verifies that the endpoints avoid 0 and 1728. The two s-values are therefore two parallel edges in the coarse unordered supersingular 2-isogeny multigraph. Frobenius fixes their vertices and exchanges the edges. This interpretation is independently verified, not merely suggested by the numerical factors.

The non-Monster examples 37 and 43 also have exactly two nonrational s-points. Therefore “one quadratic defect at level two” alone does not characterize the four missing Baby primes over all primes. The Monster/Ogg condition is essential. The prime partition is globally exhaustive because of the group-order support, while this defect statement is a finite exact conclusion at explicitly named primes.

## 4. Completeness, exceptional fibers, and coarse points

The fresh implementation finds exactly (p−1)/2 distinct Hasse roots, equal to the polynomial degree, so no Hasse root is omitted beyond Fp². Their j-images provide the complete supersingular j-set. For p>3, completeness of the u-enumeration is independently checked against

|U| = 3|J| − 2·1_{0 in J} − 1_{1728 in J}.

This follows directly from the cubic j-map: a generic j-fiber has three distinct points; j=0 has one; j=1728 has two. Those are the ramification corrections that would be lost by naively multiplying by three. The p=3 case is handled explicitly rather than applying this formula with colliding exceptional values.

All nonzero u-values are enumerated. The cusps u=0 and infinity map to j=infinity and are not supersingular elliptic points. The involution u↦4096/u preserves the enumerated set. Every s-fiber has two points except the fixed u-values ±64, which have one. Thus the quotient is formed as a set of geometric points, with fixed points counted once. This agrees with the squarefree resultants and does not substitute stack masses, multiplicities, or Fp-isomorphism classes of models for coarse geometric-point counts.

A rational s-value need not have rational u-values: Frobenius can exchange the two orientations of the isogeny. The author's method correctly forms s from the entire Fp² u-set; it does not discard this phenomenon by counting only rational u-values. The revised note states this distinction explicitly.

The author's resultant routine computes an integer resultant and then reduces modulo p, rather than computing directly over Fp(j). This is valid because resultant formation commutes with coefficient reduction here and the relevant leading coefficients do not vanish. Its originally inaccurate comment has been corrected.

## 5. Valuation-cycle note: PASS

Independent integer and graph checks verify the exponent chain 13→2→5→3→7→4 and the extra arrow 11→1. Identifying the graph vertex 4 with 13 makes a five-state cycle, not six states. With N/2 the exponent at 2 is 4 and the quotient cycle is [13,2], of length two, not three.

A ring quotient identifying 13 with 4 must have modulus dividing 9. Modulo 9 the inputs 11 and 2 agree but their images 1 and 5 disagree. Modulo 3 the inputs 2 and 5 agree but their images 5 and 3 disagree. Only the trivial modulus remains. The same modulo-3 witness rules out realizing even the selected arrows by an integer polynomial modulo 9. These are elementary proofs, not numerical conjectures.

## 6. Boundary of conclusions

The positive bridge is genuine: group-order valuations, the 2A centralizer, and supersingular Frobenius/quotient geometry are connected by specified arithmetic and cited theorems. None of these constructions supplies a map to Hamiltonian counts, tournament isomorphism classes, or the forbidden values {7,21}. Shared labels and the valuation-cycle graph do not create such a map. The author's OPEN/UNSUPPLIED boundary is appropriate.

This review does not independently reprove the historical repository's Leech/VOA theorems or every assertion in older research notes. Their mention as possible routes is not counted as evidence for the newly audited finite conclusions.


## 7. Bounded addendum: Sylow fixed-point partition — PASS

Independently reviewed `../session32-root/SYLOW_FIXED_POINT_PARTITION.md` and `../session32-root/sylow_partition.json`. No mathematical correction is needed.

For the transitive action on G/H, the stabilizer of gH is gHg^-1. Hence a subgroup K has a common fixed point precisely when K is contained in one such conjugate. The word “common” matters: the statement concerns one point fixed by every element of K, not possibly different fixed points for different elements.

If b=a, a Sylow subgroup of H is Sylow in G, and Sylow conjugacy makes every Sylow subgroup of G fix some point. If 0<b<a, Cauchy's theorem supplies an order-p subgroup in H, while orders exclude every full Sylow subgroup from all point stabilizers. If b=0, every order-p subgroup is excluded. These quantifiers match the note exactly; the middle regime does not assert that every order-p subgroup fixes a point. The maximal exponent k of a fixed-point subgroup of order p^k is b, with k=0 allowed; a Sylow subgroup of H realizes that bound.

For the Monster conjugation action on its 2A class, the point stabilizer is the actual subgroup 2.B. Independent valuation extraction from the full decimal orders agrees with every JSON entry. In particular, v_2(|2.B|)=42, v_2(|M|)=46, and the deficit is 4, not 5. Thus all three sets and all six nontrivial common-prime deficits are correct. The geometric comparison remains restricted to odd good-level primes, whereas the fixed-point lemma also covers p=2.

This is an action-theoretic interpretation derived from standard finite-group facts and the cited stabilizer, not a new classification theorem or a tournament correspondence. Exact checks and source hashes are in `sylow_addendum_exact.json`.

Audited source SHA-256:

- `session32-root/SYLOW_FIXED_POINT_PARTITION.md`: `3dc2a7db7869b4cd71d876554f6eb1592e2f26404762f675d1928aba2d274788`.
- `session32-root/sylow_partition.json`: `4072fdfebeda1a0e726cefea6c1d086c721ca32d3305fba267dc0fb45df9f346`.


## 8. Source-only addendum: Fischer support shells — PASS

Independently opened the [official ATLAS Fi22 page](https://brauer.maths.qmul.ac.uk/Atlas/v3/group/Fi22) on 2026-10-23. Its order line states |Fi22|=64,561,751,654,400=2^17·3^9·5²·7·11·13. Independent multiplication reproduces the decimal integer. Set subtraction against the previously audited Baby Monster and Monster supports gives exactly {2,3,5,7,11,13}, {17,19,23,31,47}, and {29,41,59,71} as the three consecutive prime-support shells.

`../session32-root/FISCHER_SUPPORT_SHELLS.md` and `../session32-root/fischer_support_check.json` pass this bounded review. Only the cited order and elementary support identities are certified. Neither support nesting nor divisibility alone proves subgroup inclusions, identifies element conjugacy classes, or supplies a modular level-six theorem. The concurrent PRE-AUDIT packet was not used as theorem evidence. No source edits or Git mutations were performed.

Exact check and source-qualified pin: `fischer_addendum_exact.json`.

Audited source SHA-256:

- `session32-root/FISCHER_SUPPORT_SHELLS.md`: `dc07a08481befedda81558d5c32896783241b6000aee028c847bd0b82a055d48`.
- `session32-root/fischer_support_check.json`: `910e66a6c9008271824b5686ce249e3411e36f3dd2b6b94da2c30f8cabc4c56c`.
