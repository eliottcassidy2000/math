# The first gate after the fixed-G renewal is finite and recursive

**Date:** 2026-08-16
**Status:** research reflection; theorem authority is
[THM-3521](../01-canon/theorems/THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing.md)

## Outcome first

The next old-`L` finite sheet is a unit.  For the fixed sporadic Keller map,

```text
R_5=L^271 N(G),
q=(2,5/6,-7/8),
R_5(q)!=0,
v_L(N(R_5))=-1699.
```

Thus `R_6=L^1699N(R_5)` is polynomial and `L`-coprime.  The top-only part
of the proved norm-face transform then gives exposed pair `(10663,3867)`.
The renewal faces of `R_5`, fifth image prime, degree `243`, and all-level
law remain open.

## Inheritance pass

- **Closest proved mechanism:** THM-3498 separates two divergent sheets,
  controlled by a complete Newton face, from one finite sheet, controlled by
  one exact nonzero specialization.
- **Canonical hostile:** MISTAKE-415.  A visible face pair is not a closed
  recurrence; complete residual faces and the finite sheet are independent
  obligations.
- **Corrected near miss:** THM-3506 replaced `(259,87)` by `(271,99)` and
  isolated the renewal and finite-sheet hypotheses instead of folding them
  into a scalar extrapolation.
- **Least-used sidecar:** the exact normalization chain
  `H=64LN(L)`, `J=2^35L^7N(H)`, `G=L^43N(J)`.  It turns an apparently huge
  polynomial evaluation into a short tower computation.

## Concept board

| object | predicate | operation | information at risk | cheapest test |
|---|---|---|---|---|
| top `lambda` face of `R_5` | divergent valuation | substitute old-`L` Puiseux sheets | finite branch | complete face plus one finite value |
| finite inverse point `q` | `R_5(q)!=0` | good reduction | bad denominators/ramification | invert `L,S,disc` at every tower level |
| norm tower | evaluate without expanding `J,G,R_5` | recursive cubic norm | normalization constants | explicit-`H` route and omit-`64` hostile |
| global `J` | independent audit | split outer cubic | shared tower implementation | branchwise product mod `71` |
| `R_6` | polynomial and `L`-coprime | UFD denominator clearing | other boundary primes and faces | exact old-`L` valuation only |
| fifth-image program | image prime / degree `243` | finite image and eliminant | multiplicity, separability, old factors | separate slice and squarefree-fibre gates |

## Why the cheap evaluator works

Define `P_0=L`, `P_1=H`, `P_2=J`, `P_3=G`, and `P_4=R_5`.  Each transition
has the form

```text
P_(r+1)=c_r L^e_r N(P_r),
(c_r,e_r)=(64,1),(2^35,7),(1,43),(1,271).
```

Multiplicativity unrolls this exactly on the finite-etale locus:

```text
R_5=2^477 L^271 N(L)^43 N^2(L)^7 N^3(L) N^4(L).
```

This is a value identity, not the forbidden scalar extrapolation of the
Newton packet.  It says the finite-sheet gate is only a five-factor
nonvanishing problem once the norm orbit is represented lawfully.

To evaluate `P_4(q)` modulo a good prime, there is no reason to construct
`P_4`, `P_3`, or `P_2`.  Form the universal inverse point in a cubic finite
etale algebra, recurse, and take the cubic norm on return.  The dimensions
grow only as

```text
1,3,9,27,81,
```

while the bottom polynomial has five terms.  Matrix inversion certifies
that every apparent denominator is a unit.  Direct substitution certifies
that every universal point really lies on the inverse graph.

The three norm-orbit ledgers are `(16,12,72,9,49)`,
`(12,53,22,85,76)`, and `(38,45,28,3,17)` modulo `101,103,107`; every
factor is nonzero and the unrolled product recovers `74,36,88`.

The second route stops at dimension `27` and evaluates the frozen `361`-term
`H`.  It agrees with the `L` route modulo `101,103,107`.  A literal `27 by
27` determinant agrees with iterated norm descent.  Omitting `64` changes
the final answer by `64^-27`, so the normalization is observed rather than
silently canceled.

The representation-disjoint route reconstructs all `66,146` terms of `J`.
At `p=71`, the outer inverse cubic splits into three points and direct
branch multiplication gives `R_5(q)=43`.  This agrees with the recursive
route but shares neither its bottom polynomial nor its treatment of the
outer norm.

## Connection contract

```text
source: the recursive definitions of H,J,G,R_5 in the localized cubic cover
target: one scalar R_5(q) in F_p
map: universal inverse point -> recursive cubic norms -> base-field value
preserved: exact normalization, multiplicity of all 3^4 leaves, inverse graph
lost: global support, factors, singularities, image multiplicity, other divisors
sidecar: unit/discriminant ledger at every intermediate algebra
decisive test: one nonzero lawful good reduction, challenged by a global-J split route
```

The map is intentionally lossy.  It is faithful for nonvanishing at one
point and therefore for the finite-sheet generic-unit question.  It cannot
answer irreducibility or image geometry because those are precisely among
the discarded coordinates.

## What changed on every lane

- **Divergent face:** unchanged; `(1699,615)` already supplied exact pole
  order `1699/2` twice.
- **Finite sheet:** closed positively; no hidden zero raises the valuation.
- **Polynomial tower:** the next cleared norm `R_6` now exists and is
  `L`-coprime.
- **Face tower:** only the next top face advances, to `(10663,3867)`.
  Renewal remains a two-face obligation.
- **Geometry:** unchanged.  Nothing here proves `R_5` is the fifth image
  prime or even irreducible.
- **Eliminants:** unchanged.  Degree `243` full-degree and squarefreeness
  still require a separate good-reduction fibre.

## Next exact gates

The remaining obligations should stay separate:

1. derive the `z`-top and minimum-`gamma` faces of `R_5`, or produce a sharp
   cancellation witness;
2. find a degree-`243` full-degree squarefree fibre;
3. find a bounded-degree slice of `R_5` that is squarefree and old-factor
   coprime, then use the finite-image multiplicity argument;
4. only after those steps ask whether the five-face packet renews again.

The order matters conceptually, not logically: a positive answer to any one
does not substitute for the others.

## Method lesson

When a norm tower is defined recursively but its global polynomials explode,
evaluate the recursion in the universal finite algebras instead of expanding
the intermediate norms.  Pair this with a representation-disjoint route and
print every unit gate.  This is effective exactly for local values and
finite-sheet questions; it is a counterindication for global support,
factorization, or image claims.
