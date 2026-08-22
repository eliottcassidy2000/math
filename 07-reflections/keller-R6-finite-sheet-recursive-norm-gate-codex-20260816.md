# A fifth norm orbit factor closes the next fixed Keller rung

**Date:** 2026-08-16
**Status:** research reflection; theorem authority is
[THM-3523](../01-canon/theorems/THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing.md)

## Outcome first

The canonical finite old-`L` sheet remains a unit one rung beyond THM-3521:

```text
q=(2,5/6,-7/8),
R_6=L^1699N(R_5),
R_6(q)!=0,
v_L(N(R_6))=-10663.
```

Therefore `R_7=L^10663N(R_6)` is polynomial and `L`-coprime.  Proved
THM-3522 can then, and only then, propagate the complete packet to

```text
A(10663,3867) -> A(66907,24255).
```

Nothing here advances image equations, irreducibility, the degree-`243`
separability gate, or a general Jacobian-conjecture statement.

## Inheritance pass

- **Closest proved mechanism:** THM-3521 splits the old-`L` valuation into
  two face-controlled divergent sheets and one specialization-controlled
  finite sheet.
- **Canonical hostile:** finite-sheet vanishing.  Even a complete exposed
  face fixes only the two divergent valuations; a zero on the finite branch
  raises the norm valuation and destroys coprimality after nominal clearing.
- **Corrected near miss:** MISTAKE-415 and THM-3522 separate packet renewal
  from polynomiality.  A face recurrence may be sound while the next norm
  still has an unknown denominator.
- **Least-used sidecar:** the unrolled norm orbit.  At this rung it adds only
  `N^5(L)` and one power-of-three increase in the leaf multiplicity of `H`.

## Concept board

| object | representation | preserved predicate | lost information | decisive control |
|---|---|---|---|---|
| `R_6(q)` | recursive norm tower | exact scalar value | global support and factors | one lawful nonzero reduction |
| `N^5(L)` | `243`-sheet finite algebra | complete leaf multiplicity | branch labels | all unit/discriminant gates |
| bottom `H` | frozen 361-term polynomial on 81 sheets | same terminal norm | decomposition into `L,N(L)` | flat `81x81` determinant |
| outer norm | three split branches mod 71 | all outer labels and product | a single nested algebra | `49*22*60` versus unsplit value |
| normalization | `64` at every `H` leaf | nonmonic scalar | all other geometry | omit-`64` factor `64^-81` |
| old-`L` valuation | two Puiseux sheets plus finite sheet | exact denominator exponent | other boundary divisors | complete face plus `R_6(q)!=0` |
| `R_7` packet | THM-3522 after polynomiality | five fixed-chart faces | image/factorization data | theorem hypothesis order |

## The algebraic compression

Writing

```text
P_0=L, P_1=H, P_2=J, P_3=G, P_4=R_5, P_5=R_6,
```

the fixed transitions are

```text
P_(r+1)=c_r L^e_r N(P_r),
(c_r,e_r)=(2^6,1),(2^35,7),(1,43),(1,271),(1,1699).
```

Norm multiplicativity gives the localized identity

```text
R_6=2^1431 L^1699 N(L)^271 N^2(L)^43
     N^3(L)^7 N^4(L) N^5(L).
```

This is the essential compression.  The global polynomial `R_6` need never
be expanded.  The new work relative to THM-3521 is exactly one cubic
extension and one terminal norm factor.  The scalar exponent is forced by
the actual norm depths:

```text
1431=35*27+6*81.
```

The three complete norm-orbit rows modulo `101,103,107` are

```text
(16,12,72,9,49,97),
(12,53,22,85,76,94),
(38,45,28,3,17,17).
```

Their lawful products give `R_6(q)=26,70,69`.  No factor vanishes.

## Why the two alternate representations matter

The first alternate route stops at `H` in dimension `81`.  It evaluates the
frozen polynomial directly, descends its norm through four cubic layers,
and compares that answer with the determinant of the literal `81 by 81`
multiplication matrix.  This challenges both the bottom polynomial and the
norm aggregation.  Omitting `64` changes the final value by `64^-81`, which
checks the nonmonic normalization and the exact number of leaves at once.

The second alternate route changes the outer representation.  Modulo `71`,
the top inverse cubic splits at `w=10,23,38`.  Direct inverse formulas give
three scalar points, and complete `81`-sheet towers above them give

```text
R_5=49,22,60.
```

With `L(q)=53`, the branch product is

```text
53^1699*49*22*60=9 mod 71,
```

matching the unsplit `243`-sheet algebra.  Deleting the last branch changes
the value to `25`.  Thus the same complete universe is realized once as a
single nested algebra and once as three named outer branches.

## Connection contract

```text
source: recursive fixed-map norm definitions through R_6
target: one scalar R_6(q) in a good finite field
map: universal inverse graph -> complete cubic tower -> transitive norm
preserved: all 3^5 sheets, transition scalars, inverse equations, unit gates
destroyed: support, factorization, image multiplicity, other divisors
sidecars: L/S/discriminant units, frozen-H determinant, split outer branches
cheapest decisive test: any nonzero lawful residue, with omitted-scalar and omitted-branch hostiles
```

The lossy map is exactly matched to the question.  Nonvanishing at one
specialization is enough to prove that the regular finite branch is not
identically zero.  It is not enough to detect a new image component, and the
reflection should not pretend otherwise.

## Why “243 sheets” is not the degree-243 gate

The primary evaluator has vector-space dimension `243` because it realizes
five successive cubic inverse choices.  This certifies a norm value in a
finite algebra.  The open degree-`243` problem asks for a full-degree,
squarefree/separable eliminant or fibre suitable for the image-prime program.
Those predicates are discarded by norm aggregation.  The shared integer is
a cardinality, not an implication.

## Valuation and theorem-order audit

THM-3522 gives the complete face

```text
x^10663(3xz-2y)^3867.
```

On each divergent old-`L` sheet, `x` has valuation `-1/2` and
`3xz-2y` is a unit, so the value is exactly `-10663/2`.  Completeness rules
out equal-weight cancellation.  The finite specialization proves valuation
zero on the third sheet.  Hence the norm valuation is exactly `-10663`.

Finite étaleness over `Q[a,b,c,L^-1]` places the norm in that localization.
UFD and irreducibility of `L` then turn the valuation into polynomiality and
coprimality of `R_7`.  Only after this step does THM-3522 apply again.  This
ordering prevents a recurrent support law from silently becoming an
unproved denominator law.

## Concurrent-signal audit

Incoming LRC work during the computation emphasized complete labelled
universes before quotienting and showed that a plausible transporter can
lose exactly the needed channel.  The methodological analogy reinforced the
split-branch hostile here: retain every sheet before taking the product.
There is no mathematical LRC/Keller transfer in this result, and none is
recorded in the theorem graph.

## Next exact gates

The fixed tower now presents two sharply different continuations:

1. the next local gate is `R_7(q)!=0`; the same recursion would introduce
   `N^6(L)`, a `729`-sheet bottom `L`, and a `243`-sheet bottom `H`;
2. the independent geometric gate remains the degree-`243` separability and
   image-prime program for the earlier candidate.

The first is a local denominator question; the second is global geometry.
Neither substitutes for the other, and no all-level claim should be made
from two consecutive positive finite-sheet tests.

## Method lesson

When one more norm rung looks globally impossible, first ask whether the
question is only a value.  If so, unroll the norm orbit, represent the full
finite sheet algebra, and challenge it by changing where the tower is cut
and how the outer sheets are aggregated.  This method is powerful for local
nonvanishing and exact divisor valuations.  It is a counterindication for
support, irreducibility, and image claims.
