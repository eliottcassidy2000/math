# Hadwiger–Nelson recovery: results, repairs, and surviving connections

Recovered from repository commit `2d8f820db8716b168a4cb956224ffcb59c37b998`
without modifying Git. Thirteen small, targeted source files were read and
pinned by Git blob identity and SHA-256 in `SOURCE_PINS.json`. The original
four requested theorem cards are included with nine targeted continuations;
we did not hydrate the central mistakes ledger or broadly crawl repository
blobs. Source copies total 142,941 bytes.

**Scope.** This is a recovery of seven distinct methods within the HN lane,
not a census of every repository idea, not seven new results this session,
and not a priority claim. The new-to-session synthesis is the precise
pseudoforest/Euclidean/lattice split and the cycle-lift certificate in
`PROOF.md`. Statements labelled inherited computer-assisted or cited are not
silently upgraded to independently reproved.

## 1. Seven methods actually present in the source texts

| Method | Actual recovered result | Source / status / boundary |
|---|---|---|
| Logical compactness and algebraic transfer | The real-algebraic unit-distance graph has an arithmetical chromatic invariant; finite real unit-distance configurations have real-algebraic realizations by real-closed-field transfer; under ZFC the full plane and real-algebraic plane have the same chromatic number. | **THM-4622**, known in substance, elementary arguments written there. The countable compactness argument needs no choice. Its old 5 lower bound is historical; the new external six-bound status is separately recorded below. No exact plane value follows. |
| All-lattice coloring cap | Every positive norm layer of every planar rank-two lattice has a finite-period coloring with at most three colors. Nonarithmetic real Gram forms have at most two unoriented achieved-norm directions, hence bipartite layers. | **THM-419**, corrected Oct 25, and Session61 `lattice-cap/LATTICE_CAP_REPAIR.md`, with independent audit. The complete proof embeds an arithmetic normalized layer in an imaginary quadratic field and uses a ramified-prime coloring. This replaces the invalid historical recursion. |
| Arithmetic field coloring and spectral residue obstructions | Conjugation-fixed local reduction gives upper-bound colorings of a field plane; for Q(ζ55) every conjugation-fixed prime quotient needs at least seven colors, via the norm-one torus, Kloosterman eigenvalues, Hoffman, and Weil. | **THM-4558**, Lemma R known; **THM-4628**, elementary torus identity and exact residue classification, classical spectral bounds cited. This excludes a method of producing ≤6-colorings, not ≤6-colorings of the field or plane by other means. Numerical stronger floors stay numerical. |
| Global arithmetic configuration obstruction | The fish genus-one curve has no point over the real cyclic quartic F4. Fish, EI21, and the nine-vertex subconfiguration Γ9 cannot be placed in Q(ζ55) or the specified odd-degree conjugation-stable carriers. | **THM-4662**, inherited exact computer-assisted two-descent, corrected square oracle and independent audit; **not rerun in full here**. This is nonexistence in named fields, not a planar nonexistence statement or a new chromatic bound. |
| Spectral moments and elliptic configuration curves | Closed walks in finite norm-one graphs are Kloosterman moments; W2,W3,W4 are elementary; three-step open-walk counts equal point counts of a quartic genus-one curve. Higher polygon and linkage counts connect to modular forms and periods. | **THM-4663**, mixed status: elementary identities proved, exact finite moment ranges, cited modularity, numerical or conjectural Bessel claims. No graph-chromatic improvement is claimed by the source itself. |
| Perfect codes, graph quotients, and incidence | The two ideals (2−ω), (3+ω) give seven-color triangular-lattice quotients and perfect one-codes; T² has χ=7, while T has χ=3; oriented triangular-face quotients give two Fano designs. | Session61 `hexagonal/RESULTS.md`. Classical construction, explicit coordinate recovery and finite classification. T² has three lengths; quotient K7 is not a Euclidean unit-distance witness. Fano incidence comes from quotient faces, not seven Euclidean points. |
| Metric cell lifting | At center spacing h, regular hexagonal cells have diameter 2h/√3 and same-color separation h√(7/3). Thus √(3/7)<h<√3/2 gives a proper plane seven-coloring, including boundaries. | Session61 `hexagonal/RESULTS.md`, classical construction rederived with exact metric control. h=3/4 has squared diameter3/4 and squared separation21/16. The metric inequalities, not merely seven labels, do the work. |

## 2. Corrections that must travel with the results

### 2.1 All lattices: completion, not the old mod-three argument

THM-418 correctly proves the square and triangular special cases. Its broader
wording and physical reading must not replace the actual proof. THM-419's
current header directs readers to the completed Session61 proof. The old
body is historical and contains several false steps:

* A nonrational Gram form can have **two** unoriented equal-norm directions,
  not only one: `x²+y²+(√2/2)xy` has ±e1,±e2 at norm one.
* A degenerate mod-three kernel sublattice need not be similar to the
  original lattice: `x²+6y²` restricts to `9u²+6v²`, minimum six, not the
  minimum three required by similarity with index three.
* The old census category totals sum to 14,486, not 16,486, and the old
  census was not used as proof in the repair.
* Irrational rotation angles do not alone prove lattice escape: square
  lattice vectors (1,2),(2,1) have cosine4/5.

The repaired argument uses three equal-norm directions to force rational
normalized Gram coefficients, then a normalized imaginary-quadratic field
embedding and a ramified-prime quotient. It retains a finite-period
coefficient/denominator certificate. Off a rank-two lattice does **not** mean
off algebraic number fields, and no formula χ=3+rotation-rank follows.

### 2.2 Q(ζ55): even residue degree is not enough

THM-4628 and `incoming_cylinder_prime_repair_20261013.md` give exactly

```
conjugation-fixed unramified class mod55: 19,24,29,39,54
residue degree:                         10,10,10,10,2.
```

The requirement is `−1∈<ℓ>`, not merely even order. In particular
`2^10=34 mod55`, not54, so the claimed degree-twenty fixed example at2 was
false. Neither ramified prime5 nor11 is fixed. The least q in the inert
quadratic residue description is109; the degree-ten possibilities have
`q=ℓ^5≥19^5`. The elementary Hoffman–Weil floor still rules out ≤6-color
factor-through-residue colorings. No smaller full-plane coloring is thereby
excluded: the graph-map inequality has the opposite direction.

This session independently exhausts all forty unit classes modulo55.

### 2.3 Fish descent: retain the square-oracle repair and exact scope

The curve is

```
11z0²=D(1−D),  11z1²=D(4−D),  11z2²=(1−D)(D+3).
```

It is the two-covering with class (33,33) of
`E:y²=x(x−9)(x+7)`. Over `K0=Q(√5)`, the source establishes rank E=0 and
rank of its δ-twist=2, `δ=(55+11√5)/2`. Over `F4=K0(√δ)` the four-dimensional
descent image is spanned by

```
(−7,−1), (−3,−3), (1,11), (−√5 φ,−(4−√5)), φ=(1+√5)/2,
```

and excludes (33,33). The source reports Selmer group orders4 and16.
The old square oracle incorrectly rejected `5=(√5)²`; the corrected focused
replay and independent F2 calculation recover the same target exclusion.
That correction is substantive even though the conclusion survives.

The saved replay also contains an explicit placement of Γ9 in Q(√7,i).
This session independently checks all36 pair distances exactly: precisely
14 prescribed unit edges, no vertex collisions, and no extra unit edges.
Its stored three-coloring is proper and `h-a-w_ab-w_cd-d-h` is a five-cycle,
so Γ9 itself is exactly three-chromatic. **Arithmetic placement obstruction
and high chromatic obstruction are therefore different predicates.**
This exact control does not independently reproduce the Selmer calculations.

### 2.4 Polygon ladder: keep proof, finite, cited, and numerical types apart

For `G_q=Cay(F_(q²),μ_(q+1))`, THM-4663 gives

```
W_k(q)=q^−2[(q+1)^k+(q+1)Σ_{a≠0}(−K_q(a))^k],
W2=q+1, W3=(q+1)(1−(q/3)), W4=3q(q+1).
```

The equality of torus eigenvalues with negative Kloosterman sums follows
from the exact opposite trace counts `1−η(t²−4c)` and `1+η(t²−4c)`.
For odd primes and `c∉{0,1,9}`, the open-three-chain count is

```
R3(c)=# [y²=(t²−4c)(t−c−1)(t−c+3)](F_p)
```

on the smooth projective quartic (two infinity points). The source's explicit
torsion yields divisibility by six; the configuration curve is related to
Γ1(6). The hexagon modular form is the weight-four level-six eta product.
The source gives fifth/eighth-moment exact computations through2999 and cites
modularity; its real Bessel assertions include numerical comparisons and
explicitly labelled conjectural identities. Do not relabel the whole card
as a single unconditional new theorem. Its audit corrects the K3 sign and
R3(1) arithmetic, among other citation details.

This session independently replays the trace-count identity for every
nonzero norm and trace at p=5,7,11,13, the W2–W4 formulas, and all24 smooth
quartic open-walk counts at those primes. These are bounded exact controls
of the elementary formulas, not a replay of higher modularity or periods.

## 3. Why these strands connect without solving one another

The surviving connection is **predicate preservation with explicit lifts**:

* countable compactness retains every finite coloring constraint;
* lattice normalization retains exact norm-one differences;
* reduction retains unit-edge incidence but can lose exact cycle closure;
* square-class descent retains a global arithmetic obstruction invisible
  to a bare graph coloring;
* polygon moments retain closed finite-field walk counts, not automatically
  simple cycles, injective characteristic-zero configurations, or colorings;
* a seven-label quotient needs a metric cell-size/separation certificate
  before it becomes a plane coloring;
* a functional graph can be realized by free unit-circle placements, but
  exponential inverse growth forbids uniform discrete lattice separation.

Theorem6 in `PROOF.md` makes the middle item exact: fixed lifted edge labels
extend to vertex coordinates iff their sums vanish on a fundamental cycle
basis. Its quotient-triangle control demonstrates failure even for a locally
bijective graph covering. The map's fibers and cycle displacement, rather
than matching numerical labels, are the reusable data.

## 4. Current external HN statement: cited, not independently audited

Primary web sources checked on 2026-11-02:

1. OpenAI's [commit-pinned manuscript page, *The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf), dated September23,2026, is the external six-bound source retained by Session61. The page was accessible; raw-PDF extraction was not supported by the web tool. We do not claim a fresh full-PDF read this session.
2. The [same commit's Lean root](https://raw.githubusercontent.com/openai/math/adc7f1241b42e322a6451854ab7e4b4c146bf78a/lean/OAI/Geometry/PlaneColoring/Five.lean), read directly, ends in `no_proper_five_coloring`, using proper-to-Borel-weak transfer and a weak-five-coloring contradiction. We read the root statement and proof text, not its full transitive dependency proof or an executable build.
3. [DeepMind's Formal Conjectures Erdős508 catalogue](https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/ErdosProblems/508.lean), read directly, records six as the lower bound and seven as the upper bound and links the same source. Its catalogue theorem placeholders are not the external proof.
4. For historical separation, [de Grey's 2018 primary paper](https://arxiv.org/abs/1804.02385) proves the finite five-chromatic construction. The September28,2026 [Alon–Bucić–Davies typical-norm paper](https://arxiv.org/html/2609.36220v1) still describes the older Euclidean five-bound in its introduction, illustrating that surrounding source summaries are not synchronized. Its typical-norm theorem is not a theorem about the Euclidean norm.

Accordingly this packet treats **6≤χ(R²)≤7 as an externally cited claim/theorem
with the above precise source, not a locally verified theorem**. No Lean
build, transitive axiom audit, or full independent argument review was made.
The exact choice between six and seven is not resolved by this packet.
Every new deduction here is independent of adopting that external six bound.
