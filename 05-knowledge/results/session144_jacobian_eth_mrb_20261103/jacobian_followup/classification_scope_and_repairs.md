# Session 144: what is classified, and repairs to the equivariant entry point

**Status:** inheritance audit and elementary repairs; the independent Jacobian auditor signed off on the corrected weight-monomial form and gcd power law. All citations below refer to full theorem filenames in `01-canon/theorems/`, read from `origin/main = 5e1cedf452c7cc1ec447205939a9d6cf4c0e8743`. No changes were made to the source checkout/index, and the mistakes ledger was not read or edited.

## 1. The actual all-dimensional classification already in canon

`THM-3438-weighted-lift-keller-degree-spectrum.md` is the decisive prior result. For every ambient dimension `m>=3`, it proves the exact generic-degree value spectrum

`KDeg(m)={1,3,4,5,...}`.

It constructs a map in every degree `n>=3` already in dimension three; padding preserves the degree. Every member of its normalized weighted family has geometric monodromy `S_n` and is a composition atom. In the same ambient dimension:

- every degree `n>=3` contains a map atom;
- exactly the degrees `ab`, with `a,b>=3`, contain a composite of two nonunits;
- the first mixed degree is nine;
- all maps in degrees `3,...,8` are atoms, but not necessarily equivalent to one another;
- degree nine contains both an `S_9` atom and the composite of two degree-three maps.

The monodromy-to-atom implication is one-way: a primitive field extension obstructs decomposition. An arbitrary map atom need not have primitive monodromy, because a proper intermediate field need not be realizable by an intermediate affine-space Keller map.

These classify **values and existence spectra**, not maps up to equivalence. `THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule.md` explicitly warns against the latter interpretation. Its universal anatomy is necessary data; the inverse-Jelonek realization and gluing problem remains open.

The two-sided equivalence is `F~A composed with F composed with B`, with source and target polynomial automorphisms `A,B`. This differs from conjugacy, equivariant equivalence, and stable equivalence after identity padding. Any classification request must name which is intended.

## 2. Family and scope atlas

| Entry | Exact useful content | Boundary |
|---|---|---|
| `THM-1300-jacobian-counterexample-dixmier-A3-explicit.md` | Explicit determinant/collision witness in dimension three; padding in higher dimensions | One witness does not classify all maps |
| `THM-1310-jacobian-counterexample-fiber-geometry-S3-resolvent-jelonek.md` | Fixed-map degree-three elimination, discriminant, monodromy and nonproperness geometry | Particular elimination charts need boundary care; this is not all degree-three maps |
| `THM-1305-jacobian-equivariant-anatomy-rigidity-decode.md` | Useful determinant calculus in a restricted z-affine chart; finite capped tangent calculation | Literal full normal form and power-law claims need the repairs below; capped tangent rigidity is not global classification |
| `THM-1340-engine-trichotomy-zaffine-keller.md` | Exact determinant split for z-affine maps; rationality of the projective leading-coefficient curve | Engine-degree = map-degree bridge is labelled DERIVED, not a general proved classification |
| `THM-1375-reduced-jacobian-lattice-reflection-driven.md` | Degree-three counterexamples have S_3 monodromy; weight ordering is not an invariant | Corrected deck-extension warning; group tables are constraints, not realizability |
| `THM-3438-weighted-lift-keller-degree-spectrum.md` | All geometric degree values, symmetric-monodromy atoms, atom/composite grade spectra | No classification of arbitrary maps in a fixed degree |
| `THM-3448-weighted-keller-cyclic-jelonek-inertia-family.md` | Component-generic inertia and escape counts for its weighted subfamily; cyclic inertia of every length | Assumptions on roots and reconstruction matter; not arbitrary Keller boundaries |
| `THM-4562-slice-lemma-and-the-jacobian-counterexample-is-fully-non-slice.md` | Slice reduction; the three listed dual fields/components of the fixed map are non-slice/noncoordinate | Coordinatewise conclusions do not automatically survive arbitrary changes of output coordinates |
| `THM-4567-the-jacobian-counterexample-at-the-prime-two.md` | Arithmetic of the fixed integral model at 2 and integer collisions | Integral/2-adic invariants are not automatically invariants of complex Aut×Aut equivalence |
| `THM-4625-tangent-sweep-keller-maps-fibrewise-punctured-lines-hyperelliptic-fibres-and-symmetric-monodromy.md` | Tangent-sweep inverse polynomials and generic fibrations, subject to polynomiality side conditions | Its status explicitly retracts two invalid example seeds; choosing arbitrary p is insufficient |
| `THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md` | Correct quotient Jacobian mechanism; affine-scale primitive family | A reduction is not a classification; composites already have nonlinear scales |

The composition tower in the companion note closes the literal higher-degree-scale existence gap without claiming a new primitive species. More cheaply, the equivariant shear `T_k=(x,y+x^(2k-1)z^k,z)` changes the primitive witness scale to `2-3v-t-3t^k`, of arbitrary degree k, without changing its equivariant equivalence class. Thus even primitive nonlinear-scale existence is not the open question; classification must be modulo equivariant source/target automorphisms. See `equivariant_gauge_hostile.md`.

A further boundary correction matters: the total number of finite points over a special coordinate hyperplane need not equal the generic geometric degree. The weighted quartic in THM-3438 has degree four but only three finite points over a generic point of its output hyperplane `C=0`: one on `x=0`, two on `gamma=0`, and one sheet escapes. Therefore THM-4727's generic fold count `d minus degree on x=0` requires the target hyperplane not to lie in the nonproperness set. The independent audit supplies the exact three-point fiber at `(0,1,0)`. THM-3448 already records the necessary boundary ledger.

Do not call these maps global topological coverings or assume surjectivity. They are etale polynomial endomorphisms; after restricting away from nonproper values they become finite unramified coverings. For the original witness the missing target curve is nonempty, as the independent fixed-map fiber audit proves.

## 3. Correct complete weight-monomial form

Let the source weights be `(1,-1,-k)`, `k>=2`, and target weights `(-k,-1,1)`. Put `t=xy`, `s=x^k z`. The exact **unique** form of every equivariant polynomial map is

```
F_1=z A(t,s)+y^k B(t),
F_2=y C(t)+x^(k-1) z D(t,s),
F_3=x E(t,s),
```

with `A,D,E` bivariate and `B,C` univariate. This repairs THM-1305's assertion with univariate `A,D`.

**Proof.** A monomial `x^a y^b z^c` of weight `-k` satisfies `a-b-kc=-k`. If `c=0`, it is `y^k t^a`; if `c>=1`, it is `z t^b s^(c-1)`. For weight `-1`, the same split gives `y t^a` or `x^(k-1)z t^b s^(c-1)`. For weight `+1`, it is `x t^b s^c`. The split by `c=0` versus `c>=1` is disjoint, proving uniqueness. ∎

For example `x^k z^2=z s` has weight `-k` and cannot occur in THM-1305's univariate-coefficient form. The older form is a restricted chart with the first two rows affine in `z`, not every equivariant map. This explains why higher-z-degree weighted lifts and iterates are not ruled out by chart rigidity.

A module formulation allowing all coefficients in `K[t,s]` is also possible but is not unique: the generators satisfy `s y^k=t^k z` and `s y=t x^(k-1)z`. The unique form above is preferable for a computational chart.

## 4. Correct power law for the quadratic determinant coefficient

In the restricted chart with constant coefficient of `s` in `E`, suppose the displayed determinant coefficient is

`c_2=-2 A' D+(k+1) A D'`.

For **nonzero** `A,D in K[t]` over characteristic zero, set `g=gcd(k+1,2)`. Then

```
c_2=0 iff
A=a V^((k+1)/g), D=d V^(2/g)
```

for nonzero constants `a,d` and a polynomial `V` (including constants).

**Proof.** At every irreducible factor, logarithmic derivatives give `(k+1) ord(D)=2 ord(A)`. The primitive nonnegative integer solution is `ord(A)=((k+1)/g)r`, `ord(D)=(2/g)r`. Multiplying the factors with exponents `r` gives `V`; the only untracked factors are units. Conversely substitution verifies the identity. ∎

For even `k`, the exponents are `(k+1,2)`, as in the old formula. For odd `k`, they are `((k+1)/2,1)`. The smallest hostile is `k=3,A=t^2,D=t`, which satisfies `c_2=0` but is not `A=V^4,D=delta V^2` with polynomial `V`. Over a non-algebraically-closed field, scalar `a` also cannot always be absorbed into `V`. The degenerate cases `A=0` or `D=0` satisfy this one coefficient equation automatically and must be handled separately.

This does **not** construct a new Keller map: the remaining determinant equations and nonzero constant condition must still be solved. It does invalidate an exclusion search that silently imposed the stronger powers at odd `k`.

## 5. Weight sums are weaker than the tangent representation gate

For an equivariant Keller map between diagonal torus representations with no zero weights, the **entire weight multisets agree**, not merely their sums.

Indeed `F(0)=0`, since it must be fixed by every torus element and the target has no zero weight. Differentiating equivariance at the origin makes `JF(0)` an intertwiner of the two representations. It is invertible, so each weight multiplicity agrees.

For example `(1,-1,-3)` and `(1,-2,-2)` have the same weight sum and the same `W=4`, but no equivariant Keller map can connect them: the derivative has no match for either negative source weight. This sharpens the necessary gate in THM-4727, without changing its quotient determinant formula. Also, its coordinate hyperplanes are invariant, not pointwise fixed; with no zero weights the fixed locus is only the origin.

## 6. Primary-source checks and verification limits

The following accessible primary pages were actually opened on 2026-11-03:

- [AFP verification entry](https://isa-afp.org/entries/Jacobian_Counterexample.html): identifies the formal three-dimensional witness and identity padding. This worker did not execute the Isabelle session.
- [Gao, arXiv:2608.00222v1](https://arxiv.org/html/2608.00222v1): tangent sweeps, polynomiality conditions, higher-dimensional direction-field constructions, and explicit examples. It is a construction framework, not a claim that every Keller counterexample is in the family.
- [Gallagher weighted-lift source](https://github.com/algal/jacobianfun/blob/main/RESEARCH.md): the endpoint/integral seed conditions and exact degree-four construction underlying THM-3438.
- [van den Essen, arXiv:2609.17795v1](https://arxiv.org/html/2609.17795v1): an elementary reconstruction equivalent up to linear coordinate changes to the existing witness, rather than a classification of all counterexamples.

External formal-build claims remain CITED. The determinant and collision identities were recomputed here independently, and the new composition and repair statements are proved above and checked by the companion. No implications for ETH, MRB, or the planar conjecture follow from matching language or the presence of the same small integers.
