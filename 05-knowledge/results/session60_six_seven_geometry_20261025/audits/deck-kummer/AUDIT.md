# Session60 independent deck, Kummer, and Hodge audit

**Gate: PASS after one verification-only correction.** No theorem correction is
required in the audited frozen packets. The inherited Session59 root-cover
surface theorem is an explicit input, not silently counted as independently
reproved here. No historical priority assertion or later 63-involution
classification is in scope.

Audit input directories:

- `/tmp/session60-deck-k3`
- `/tmp/session60-hourglass-symmetries`

## 1. Found defect and authorized narrow correction

The original deck checker had 28 bare `assert` sites. Its normal run passed,
but `python -O` removed these assertions; matching JSON output did **not**
validate the skipped tests. The source has now been corrected, with parent
authorization, to explicit `require`/`RuntimeError` gates. No tested expected
value or mathematical statement changed. There are **91,064 executed gates**
in both modes. Normal, optimized, and `EXACT_CHECKS.json` are byte-identical.
A deliberately false genus expectation is rejected with exit code 1 in
**both** modes. The original checker is retained verbatim here as
`original_check_deck_characters.py`; it is provenance, not a current checker.

- Original source SHA256: `b1c17c4475c6d12b7ab4dfbd2ea65488a842eae726e599ed93bedc4a5a188031`.
- Corrected source SHA256: `0dd019b0173baa60ea10772d8afe15342f37a7ee870329d329c895b8accce689`.
- Full provenance: `SOURCE_PROVENANCE.json`.
- Corrected deck certificate: `../../deck-k3/VERIFICATION_CERTIFICATE.json`.
- Corrected source/output pins: `../../deck-k3/SHA256SUMS`.

The original `deck.normal.json` and `deck.optimized.json` in this audit
folder are retained initial-run evidence. Use `corrected-deck.*.json` for
the genuine optimization-safe validation.

The two symmetry checkers also pass in both modes, with byte-identical
outputs: 116 source checks and 372 independent branch-configuration checks.

## 2. Independent exact controls

`independent_audit.py` is a separate implementation using explicit exceptions,
not source imports. It passes **825 gates** under normal and optimized Python;
`audit.normal.json` and `audit.optimized.json` are byte-identical. It checks:

- monic-basis character indexing and genus against independently assembled
  complete-intersection Hilbert coefficients for m=5,...,9;
- character orthogonality, the seven deletion invariant lines, and the deck
  determinant twist;
- all hourglass triple incidences, four-subsets, six-root deletion pairs,
  and resolved eight-fixed-point totals;
- all actual D12 matrix symmetries and the determinant ratio on their
  explicitly computed four-dimensional quadratic-relation spaces;
- an explicit original-normal-to-Veronese matrix and all Gale relations;
- the scalar 27 by a symbolic resultant, independently of the recorded
  evaluation-product calculation;
- the sixth-root branch transformation by reduction modulo d²+3, both
  elliptic-map equations/differentials, and the CM divisor determinant.

One extra reproducibility detail: for s=1/(t-3), the conic column is

```
T*(1,s,s²)^t,  T=[[1,5,7],[1,6,8],[1,4,3]],
```

and at the six packet parameters its scaling relative to the original
normal is respectively `3,1/9,3/16,1/4,1,3/25`. This explicitly verifies
the omitted-but-standard projective column normalization in the Gale proof.

## 3. Mathematical review

### Deck eigensheaves and canonical characters: PASS

Even subsets index projective characters because the common sign is trivial.
On L_j nonzero, deleting j from an even subset is a bijection to all monic
basis supports; the eigensheaf therefore is O(-|I|/2). Duality and crepancy
give H20 character component H0(O(|I|/2-3)). At seven roots, exactly the
seven weight-six words occur. Their only F2 dependence is their total sum;
they are not a seven-element Fano character group.

The residue generator contributes the determinant sign. The independent
matrix calculation confirms the ambient determinant divided by the relation
space determinant equals det(base A) for every D12 element. This is not the
untwisted root permutation representation.

### Resolutions, symplectic group, and lower bound: PASS

The six-root symplectic deck subgroup is precisely total-parity kernel
(C2)^4, with fifteen nonidentity two-flips. Ordinary intersections contribute
eight smooth fixed points. At a triple, four nodes each contribute two
exceptional fixed directions; the complementary sign eigenspace would need
four concurrent lines and is excluded. Distinct nonidentity stabilizers
cannot share a fixed point in the smooth symplectic surface: an order-two
nonidentity element of SL2(C) is -I. Hence there are fifteen quotient nodes.

Lefschetz averaging gives invariant H2 rank seven, and a rank-fifteen
algebraic coinvariant part. For each branch triple, the sum of its four
exceptional classes is invariant, has square -8, and is orthogonal to the
hyperplane class and the other triple sums. These t+1 independent invariant
classes yield rho >=16+t, only a lower bound.

The normal double plane is **not** the quotient of the already resolved K3
as a singular model when triples occur: the former has D4 points, whereas
the latter has fifteen A1 points. The packet preserves that distinction.

For deletion Y7 -> S_j, a target node not involving j has an étale local
pair upstairs, so its blowup lifts. A source node involving j maps to a
smooth target point and its exceptional curve is contracted. Thus the
resolved map is a proper generically finite degree-two morphism, not
necessarily finite, and its quotient need not equal the minimal K3 model.

### Rational Hodge decomposition and Picard arithmetic: PASS

On transcendental H2, quotient transfer is valid after resolving the graph;
exceptional/divisor-supported terms vanish. Pullback scales the pairing by
the group order. The character projectors are rational. Every character
without H20 or H02 is entirely rational (1,1), so contributes only to NS.
The jth deletion isolates precisely c_j. Distinct character spaces are
orthogonal. Consequently

```
T(Y7)_Q = orthogonal direct sum f_j^* T(S_j)_Q,
(f_j^*a,f_j^*b)=2(a,b),
rho(Y7)=sum_j rho(S_j)-76.
```

This does **not** assert an integral lattice decomposition or an isometry
without scaling. The numerical inputs give b2(Y7)=78 and b2(S_j)=22.
Central rho20 and outer transitivity yield rho(Y7)=6rho_outer-56, with
current possibilities 52,58,64; the outer value remains undetermined here.

### Shioda scalar, CM, and original Kummer property: PASS

The independent resultant gives q2*Res_t(q,f)=27B, so the genus-two model
supplies W²=27B over Q, not literally w²=B. W=3sqrt(3)w identifies them over
Q(sqrt(3)). Neither Q-isomorphism nor Q-nonisomorphism is proved; twisting
the genus-two equation by a scalar only squares that scalar in the product
construction. The packet's arithmetic boundary is correct.

The branch Möbius transformation identifies the genus-two curve over C with
v²=u6-1. Both displayed elliptic quotients have j=0 and independent pulled
back differentials. Thus the Jacobian is isogenous, as an abelian variety,
to E0×E0; it is **not** thereby identified with the product as a principally
polarized abelian surface. The determinant -3 of the four divisor classes
proves rho=4, and Kummer resolution adds sixteen algebraic classes, giving
rho=20. Rational Hodge transfer then gives rho(S0)=20.

Crucially, being Kummer is **not** inferred from having a Kummer quotient.
The Gale matrix gives an invertible complex coordinate change to the
classical diagonal three-quadric Kummer model itself. The condition that
six normals lie on a smooth conic is used as a sufficient criterion here,
not asserted to be necessary for every possible Kummer presentation.
Root rescaling square roots are retained as a complex-field boundary.

## 4. Primary-source hypothesis audit

[Garbagnati–Sarti](https://sarti.pages.math.cnrs.fr/home/GarbagnatiSarti_RMJ.pdf),
§5.1 equation (2), printed p.22, gives the diagonal Kummer equations. Although
that section starts in a general Picard-17 setting, §10.2, printed pp.54–55,
explicitly gives the same Kummer subfamily for pairwise distinct complex
parameters. Thus the identification used here is not a generic-Picard-rank
assertion. The surrounding smooth/general-position quotient description
must not replace the packet's independent triple-point resolution argument.

[Clingher–Malmendier](https://arxiv.org/pdf/1809.02003), §4.1, printed pp.15–16,
starts with a smooth genus-two curve, equivalently six distinct branch
points. Lemma 4.2 gives birational identification of the Shioda sextic with
its Jacobian Kummer; Remark 4.3 gives the six-tangent-line form and converse.
No generic Picard assumption is imposed there. The conclusion about smooth
minimal K3 resolutions follows from birationality of smooth minimal K3s,
not an assertion that the unresolved singular models are identical.

Both primary PDFs were independently inspected on 2026-10-25. All remaining
claims above are checked from the packet's explicit arguments and local
calculations. No arithmetic square construction, full automorphism-group
classification, principally polarized product isomorphism, or integral
transcendental lattice identification is inferred.

## 5. Reproduce

```
python3 independent_audit.py > audit.normal.json
python3 -O independent_audit.py > audit.optimized.json
cmp audit.normal.json audit.optimized.json
```

Run the corrected deck and two symmetry checkers similarly. The saved
negative-control traces document failure of an intentionally incorrect
expected genus in both optimization modes. `SHA256SUMS` pins this audit;
`FROZEN_SOURCE_SHA256SUMS` pins the corrected audited packets.
