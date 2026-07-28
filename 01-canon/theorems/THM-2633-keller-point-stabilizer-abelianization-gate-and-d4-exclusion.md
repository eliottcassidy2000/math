---
id: THM-2633
title: "Affine Keller point-stabilizer abelianization gate and D4 exclusion"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  Let G be the
  geometric monodromy group of a dominant complex affine-space Keller map
  and H a point stabilizer in its generic-fibre action.  Then H surjects
  onto G^ab; equivalently, restriction Hom(G,C_l) -> Hom(H,C_l) is injective
  for every prime l.  The proof uses two affine facts: every target divisor
  has a generic affine inverse branch because the Keller map is open, and a
  polynomial on affine space cannot be a nonconstant unit.  Consequently
  every Jelonek component has at least one surviving sheet.  At degree four
  the groups C4, V4, and D4 fail the stabilizer gate, leaving only A4 and S4.
  Thus D4 monodromy is excluded in every dimension.  No A4/S4 exclusion,
  JC(2), G1, DC(2), or full Jacobian-conjecture closure follows.
source: root-long-frontiers-2026-07-28-point-stabilizer-gate
depends_on:
  - THM-2627-d4-jelonek-quadratic-character-rank-and-component-gate
related:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-2612-d4-deck-pole-tax-and-depressed-resolvent-gcd-gate
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-2628-d4-opposite-pair-escape-and-deck-pole-census
  - MISTAKE-297
script: 04-computation/jacobian_keller_point_stabilizer_abelianization_thm2633.py
output: 05-knowledge/results/jacobian_keller_point_stabilizer_abelianization_thm2633.out
script_sha256: 53acfb367d6fa2089492dbaeb8cacd205fc2b761c6b60ac95b7343317c8ee05b
output_sha256: 57efaa6ae3c67e70a951a18af9e5cc1647df13d5bb1937ebe88acdece019d0e3
hash_basis: LF-normalized bytes
---

# THM-2633 -- a Keller point stabilizer generates the monodromy abelianization

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**

The earlier `D4` boundary program separated three ledgers: Kummer character
parity, local inertia, and which inverse branches remain in the affine
source.  The missing global observation is that an affine-space Keller map
cannot omit an entire target divisor.  Once this is restored, every boundary
inertia group fixes an affine sheet.  A character invisible to the point
stabilizer is then invisible on every boundary component, contradicting its
nonzero Kummer class.

This yields a group-theoretic gate valid in every dimension and excludes the
entire `D4` degree-four lane.

## 1. Every target divisor has a generic affine inverse branch

Let

```text
F:X=A^n_C -> Y=A^n_C                                    (1)
```

be a dominant polynomial Keller map of generic degree `d`.  Thus `F` is
etale and in particular quasi-finite and open.  Let

```text
D=V(f) subset Y                                          (2)
```

be any irreducible target divisor.  Its inverse image is

```text
X_D=X times_Y D=V(f o F).                                (3)
```

This inverse image is nonempty.  Otherwise `f o F` would be a unit of
`C[x_1,...,x_n]`, hence a nonzero constant `c`.  Dominance makes

```text
F^*:C[Y] -> C[X]                                         (4)
```

injective, so `F^*(f-c)=0` would force `f=c`, contradicting that `D` is a
divisor.

The projection `X_D -> D` is a base change of `F`, hence is etale and open.
Its nonempty image is therefore a nonempty open subset of the irreducible
scheme `D`, so it contains the generic point.  After geometric base change,

```text
the generic fibre over every target divisor contains an affine point.     (5)
```

The openness input is standard; see
[Stacks Project, Tag 02GH, Lemma 29.37.13](https://stacks.math.columbia.edu/tag/02GH).

Now let `D` be an irreducible component of the Jelonek set and let `k_D` be
the generic geometric cardinality of the affine fibre along `D`.  Equation
(5) gives the new lower bound

```text
k_D>=1.                                                   (6)
```

The usual finite Zariski-main normalization gives the complementary strict
upper bound.  Factor `F` as an open immersion `X -> Xbar` followed by the
finite normalization `Xbar -> Y` in `C(X)`.  If all `d` geometric points of
`Xbar` over the generic point of `D` lay in `X`, then the finite image of the
closed complement `Xbar minus X` would miss that generic point.  After
shrinking `Y` around it, the complement disappears and `F` is finite, hence
proper, over that neighborhood.  This contradicts that `D` is a component of
the nonproperness set.  Hence

```text
1<=k_D<=d-1.                                              (7)
```

In particular, a codimension-one component can never be a total-loss
component.  Empty fibres may still occur in codimension at least two.

## 2. The point-stabilizer abelianization gate

Let `K=C(Y)`, let `Omega/K` be the Galois closure of `C(X)/K`, and write

```text
G=Gal(Omega/K),          C(X)=Omega^H,                    (8)
```

where `H` is a point stabilizer in the transitive action on the `d` generic
inverse sheets.  Then

```text
H -> G^ab is surjective.                                  (9)
```

Equivalently, for every prime `l`, restriction is injective:

```text
res_H^G:Hom(G,C_l) -> Hom(H,C_l).                         (10)
```

### Proof

If `d=1`, the assertion is immediate.  For `d>1`, write the reduced Jelonek
hypersurface as

```text
A_F=V(f_1...f_c),          U=Y minus A_F.                (11)
```

Over `U`, the Galois closure is a connected finite-etale `G`-cover.  At the
geometric generic point of each component `D_i`, pass to a transverse strict
henselian DVR and choose a tame inertia generator `g_i`.  By (5), at least one
affine inverse branch extends across that DVR.  Since `F` is etale there, the
branch is fixed pointwise by inertia.  Thus `g_i` fixes a sheet in `G/H`, or
equivalently

```text
g_i is conjugate into H.                                  (12)
```

Suppose that a nonzero character

```text
chi:G -> C_l                                               (13)
```

restricts trivially to `H`.  Characters are conjugacy invariant, so (12)
gives `chi(g_i)=0` for every component.  On the other hand, quotienting the
connected `G`-cover by `ker(chi)` gives a nontrivial connected `C_l`-cover of
`U`, hence a nonzero Kummer class.  Since

```text
O(U)=C[y_1,...,y_n,1/(f_1...f_c)]
```

is a localization of a UFD, `Pic(U)=0`; its units are a nonzero constant
times a Laurent monomial in the pairwise nonassociate irreducibles `f_i`.
As `C^*` is `l`-divisible, the Kummer sequence therefore gives

```text
H^1_et(U,mu_l)
 =O(U)^*/O(U)^{*l}
 =F_l^c.                                                  (14)
```

The coordinates in (14) are the divisor valuations modulo `l`, equivalently
the component-inertia values `chi(g_i)` after choosing generators.  They all
vanish, so the Kummer class is zero, contradicting (13).  This proves (10).

Finally, let `B` be the image of `H` in the finite abelian group `G^ab`.  If
`B` were proper, the nonzero finite abelian quotient `G^ab/B` would have a
quotient `C_l` for some prime `l`.  The resulting character of `G` would be
nonzero and trivial on `H`, contradicting (10).  This proves (9).

## 3. Exact degree-four exclusion

The transitive subgroups of `S_4` have the following point-stabilizer gate:

| `G` | point stabilizer `H` | `G^ab` | gate |
|:--|:--|:--|:--:|
| `C4` | `1` | `C4` | fail |
| `V4` | `1` | `V4` | fail |
| `D4` | `C2` | `C2 x C2` | fail |
| `A4` | `C3` | `C3` | pass |
| `S4` | `S3` | `C2` | pass |

Therefore a degree-four affine-space Keller map has geometric monodromy only

```text
A4 or S4.                                                 (15)
```

In particular, `D4` is impossible in every dimension.  Explicitly, its deck
character

```text
chi_deck:D4 -> C2                                         (16)
```

is nonzero but trivial on a vertex stabilizer `H=<s>`.  Its four odd elements
are the two edge reflections and two four-cycles, all derangements.  Kummer
nontriviality would force an odd Jelonek component, while (5) forces a fixed
affine sheet there.  This is the concrete `D4` contradiction.

The elimination of `C4,V4` agrees with the older polynomial-Galois gate, but
the `D4` exclusion is new: its source extension is not Galois and its generic
deck involution need not extend polynomially.

## 4. Why this does not repeat MISTAKE-297

MISTAKE-297 incorrectly extended a deck action from the finite-etale open
source across the Jelonek boundary and then applied Smith theory to the
completed affine source.  Nothing above extends a deck transformation or
uses a fixed-point theorem.  The only open morphism is the original Keller
map `F`, which is etale on all of `X`.  The proof then uses:

```text
affine-space units + dominance + openness of F
  -> one affine branch over every divisor;

one affine branch + local etaleness
  -> boundary inertia fixes a sheet;

Kummer support
  -> no character invisible to the point stabilizer.     (17)
```

The pole theorem THM-2612 and the opposite-pair census THM-2628 remain valid
conditional structure theorems, but their hypothetical `D4` polynomial
Keller branch is now empty.

## 5. Consequences for the inverse-quartic ledger

THM-2621's valuation law was proved for abstract local resultant controls with
`k_D=0,1,2,3`.  Equation (7) removes the zero row for an actual polynomial
degree-four Keller realization:

```text
k_D in {1,2,3},

the full pole is in a_1, a_2, or a_3, never a_0.          (18)
```

Likewise, THM-2628's `k_D=0` inertia rows remain correct local
normalization/group-theory possibilities but are not affine-space Keller
Jelonek components.

The previously intended two-component `D4` normal form is internally
consistent as a conditional parity table: one component is a split
diagonal-reflection singleton with a full `a_1` pole, while the other is
deck-odd, has `k_D=0`, and has a full `a_0` pole.  Equation (6) shows that the
second component cannot occur.  The normal form is therefore a vacuous
conditional boundary, not a remaining monodromy lane.

For G1 and the planar field-degree-four problem, all `D4`-specific quadratic
intermediate-field cases disappear.  The live groups `A4,S4` have maximal
point stabilizers and no proper intermediate field in the source quartic.
This sharpens those frontiers but does not exclude either surviving group.

## 6. Hostile and positive controls

Three controls locate every load-bearing hypothesis.

1. The non-etale dominant map

   ```text
   (x,y) -> (u,v)=(x,xy)                                  (19)
   ```

   meets `u=0` only at `(0,0)` in the target, not generically.  Its Jacobian
   is `x`, so the failed openness step is visible.

2. The open immersion

   ```text
   G_m x A^(n-1) -> A^n                                  (20)
   ```

   omits a divisor because its source ring has a nonconstant unit.  The
   affine-space unit hypothesis in Section 1 is therefore essential.

3. THM-2627's dominant `D4` hostile has an irreducible nonproper divisor, but
   its Jacobian exposes additional finite critical-value divisors.  Its
   Kummer support is not confined to the Jelonek boundary, so it does not
   satisfy the Keller character argument.

As a positive group control, the in-repo degree-three `S3` Keller anatomy has
point stabilizer `C2` mapping onto `S3^ab=C2`.  THM-2473/2546 give
`k_D=1` generically on its Jelonek hypersurface; its empty fibres occur only
on a codimension-two curve.  This is exactly the distinction predicted by
(6).

## 7. Exact evidence and scope

Run

```text
python3 04-computation/jacobian_keller_point_stabilizer_abelianization_thm2633.py
python3 -O 04-computation/jacobian_keller_point_stabilizer_abelianization_thm2633.py
```

The companion verifies `H[G,G]=G` or its failure for `C3,S3` and every
transitive degree-four group, checks the complete `D4` deck-odd derangement
class and fixed-sheet histogram, and retains all guards under optimized
execution.  Normal and optimized modes byte-match the stored transcript.

Two independent hostile audits checked the affine unit argument, base-change
openness, strict-henselian fixed-sheet step, Kummer-coordinate identification,
and equivalence of (9) and (10).  The exact companion is a finite-group audit;
the geometric proof is Sections 1--2, not an extrapolation from enumeration.

No `A4` or `S4` exclusion, polynomial inverse, JC(2), G1, DC(2), or general
Jacobian-conjecture closure follows.

QED.
