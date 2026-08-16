---
id: THM-3504
title: "Level-four sporadic Keller image prime and four-component nonproperness"
status: >
  PROVED + VERIFIED-EXACT.  For the fixed sporadic Keller map of THM-2473,
  the polynomial G=L^43 N(J) from THM-3498 is absolutely irreducible and
  coprime to L,H,J.  The restriction V(J)->V(G) has generic degree one,
  closure(F(V(J)))=V(G), and S_(F^4)=V(LHJ G) has exactly four irreducible
  components.  The multiplicity and distinctness gates come from an exact
  squarefree mod-1009 slice of degree 542, not a global expansion of N(J).
  This is a theorem about one fixed map, not a Jacobian-conjecture or
  all-level newest-factor result.
source: codex/level4-norm-J/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
related:
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_level_four_norm_J_mod1009_slice_probe_20260816.py
output: 05-knowledge/results/keller_level_four_norm_J_mod1009_slice_probe_20260816.out
script_sha256: d5d6a45cbbb6fe78fc572563f5d95394df79e8c640bd7640ac4763a7b3643410
output_sha256: 84920c4dbb254aa9a9d3490ca6f5767d23d29a17d137f9e950d62cf49b5f6e69
coefficient_ledger_sha256: 47fba77866ee50d00fcae28b834e8a0b0c18a4cf52c2cd1b9c05155410c91d00
hash_basis: raw LF bytes for files; ascending coefficient ledger of G(A,1,1) in F_1009[A]
---

# THM-3504 -- the fourth image prime of the fixed sporadic Keller map

**PROVED + VERIFIED-EXACT.**

Let `F:C^3->C^3` be the fixed sporadic Keller map of THM-2473.  Retain the
normalizations of THM-3495 and THM-3498:

```text
N(H)=J/(2^35 L^7),
G=L^43 N(J) in Q[a,b,c],
gcd(G,L)=1.                                               (1)
```

Here `J` is primitive and absolutely irreducible, and `N` is the cubic
function-field norm induced by `F`.

## 1. The theorem

The fourth cleared norm is absolutely irreducible, up to its fixed nonzero
rational scalar, and is new:

```text
G is absolutely irreducible,       gcd(G,LHJ)=1.          (2)
```

Moreover,

```text
closure(F(V(J)))=V(G),                                  (3)
```

and the dominant restriction

```text
V(J) -> V(G)                                             (4)
```

has generic degree one.  Consequently the exact composition law of
THM-2576 gives

```text
S_(F^4)=V(LHJ G),                                       (5)
```

with exactly four irreducible components.  Together with THM-3498, the
degree-81 eliminant has square class

```text
[Delta_4]=[2G].                                         (6)
```

The new input is a squarefree, old-factor-coprime specialization of `G`.
No global expansion of `N(J)` is formed.

## 2. The image of V(J) is a single hypersurface

Put `U=A^3\V(L)`.  THM-2473 proves that

```text
F^-1(U) -> U                                             (7)
```

is finite etale of degree three.  The intersection of `V(J)` with the
source of (7) is nonempty.  Indeed, THM-3495 supplies

```text
q=(3,-1,0),              H(q)=0.
```

Direct evaluation of the fixed map gives

```text
p=F(q)=(10,-46,33),      L(p)=-504,
F(p)=(-1854753363,121225664,-19180),
L(F(p))=-69753247104.                                   (8)
```

Because `L(p)!=0`, the norm identity for `H` and the vanishing `H(q)=0`
give `J(p)=0`.  The last value in (8) says that `p` lies in
`V(J) intersect F^-1(U)`.  Since `V(J)` is irreducible, this is a dense
open subset of `V(J)`.  The companion also substitutes directly into the
frozen `H` and `J` ledgers and verifies both displayed zeroes independently
of this norm-support inference.

Its image under the finite map (7) is an irreducible closed surface in `U`.
Let `P` be an irreducible equation of its closure in `A^3`.  The zero support
of a finite norm is the finite image of the zero support of its argument.
THM-3498 has already removed the only possible boundary denominator and
proved that the cleared numerator has no `L` factor.  Unique factorization
therefore gives

```text
G=u P^e,       u in Q^*,       1<=e<=3,                 (9)
```

Here `J` has divisor order one at its generic point and (7) is etale, so the
norm-divisor exponent `e` is exactly the generic degree of (4).  Thus the
global problem is not an unrestricted factorization problem: only the
exponent `e` remains.

## 3. A sharp degree bound on the slice b=c=1

Write the remaining target coordinate as `A`.  On this slice,

```text
L=A(27A-2),       S=27A-1,       E(w)=Lw^3+w-2.          (10)
```

At `A=infinity`, with `t=1/A`, the Newton polygon of `E` gives

```text
v_t(w)=2/3
```

on all three Puiseux branches.  THM-3495's reduced inverse formulas then
give the bounds

```text
v_t(q_x)=2/3,       v_t(q_y)>=-2/3,       v_t(q_z)>=-2. (11)
```

For a monomial `x^i y^j z^k` of `J`, its pole order is therefore at most

```text
(2/3)(-i+j+3k).
```

Exact inspection of THM-3495's frozen 66,146-term ledger gives

```text
max_(supp J)(-i+j+3k)=228.                              (12)
```

Hence each conjugate value `J(q)` has pole order at most `152`, its cubic
norm has pole order at most `456`, and multiplication by `L^43` adds at most
`86`.  Therefore

```text
deg_A G(A,1,1) <= 456+86=542.                           (13)
```

This bound is the anti-aliasing sidecar for the finite-field computation.
It is proved before interpolation and is strictly below `1009`.

## 4. The exact F_1009 specialization

Reduce the slice modulo `1009`.  The reduced inverse coordinates have only
powers of `2`, `L`, and the primitive linear polynomial `S` in their displayed
denominators.  Since (1) is polynomial over `Q`, Gauss's lemma cancels the
`S` powers without introducing any denominator except a power of `2`.
Thus reduction modulo `1009` is lawful.

There are `1006` residues at which `L`, `S`, and the cubic discriminant are
nonzero.  At `543` of them the companion:

1. constructs `F_1009[w]/(Lw^3+w-2)`;
2. verifies the inverse graph `F(q)= (A,1,1)`;
3. evaluates the full polynomial `J(q)`;
4. takes its cubic regular-representation norm; and
5. multiplies by `L^43`.

The bound (13) makes those `543` values determining.  Twelve further regular
residues are held out and agree with fresh direct norm evaluations.  The
interpolated polynomial, denoted `G_1009`, has

```text
deg G_1009=542,
gcd(G_1009,G_1009')=1,                                 (14)
```

so (13) is attained.  FLINT factorization gives the degree/exponent ledger

```text
(1,1),(2,1),(2,1),(4,1),(12,1),(21,1),(500,1).          (15)
```

In particular every exponent is one.  The ascending coefficient ledger has
SHA256

```text
47fba77866ee50d00fcae28b834e8a0b0c18a4cf52c2cd1b9c05155410c91d00. (16)
```

The old specializations have degrees `2`, `14`, and `86`, and direct gcds
give

```text
deg gcd(G_1009,L_1009)=0,
deg gcd(G_1009,H_1009)=0,
deg gcd(G_1009,J_1009)=0.                               (17)
```

## 5. Independent routes and hostile controls

The primary evaluator groups the `J` ledger by `(i,k)` and sums along the
`y` exponent, leaving `4,160` algebra-valued groups.  Every held-out point is
also evaluated after grouping by `(i,j)` and summing along `z`, leaving
`5,657` groups.  The two aggregation orders agree.

After either evaluation, the cubic norm is computed in two representations:

```text
closed norm formula = determinant of the literal 3 by 3 multiplication map.
```

They agree at every interpolation and held-out point.  A polynomial built
from only `542` nodes fails at the omitted node, confirming that the top
degree in (14) is real rather than padding.  Two hostile controls also fire:

- replacing `G_1009` by its square makes the derivative gcd have degree
  `542`;
- injecting the old factor `H_1009` makes the old-factor gcd have degree
  `14`.

These controls distinguish the two exact failure modes needed below.

## 6. Multiplicity one and the new image prime

Specialize (9) at `b=c=1` and reduce modulo `1009`.  The result is nonconstant
of degree `542`, so the reduction of `P` is nonconstant and the scalar is
nonzero.  If `e` were `2` or `3`, then `G_1009` would have a repeated factor;
characteristic `1009` divides neither exponent.  Equation (14) rules this
out.  Therefore

```text
e=1.                                                     (18)
```

Thus `G=uP`, proving (3) and the generic-degree-one assertion (4).  The image
of the irreducible complex hypersurface `V(J)` is geometrically irreducible,
so `P`, and hence `G`, is absolutely irreducible.

THM-3498 already gives `gcd(G,L)=1`.  If the new irreducible `G` were
associated to either irreducible old factor `H` or `J`, their lawful
specializations would have a common nonconstant factor.  Equation (17)
excludes both possibilities.  This proves (2).

## 7. Four nonproperness components

THM-2576 proves the exact set law

```text
S_(F^m)=union_(r=0)^(m-1) F^r(S_F).                     (19)
```

The left side is closed.  THM-2576 and THM-3495 identify the first three
successive closures as `V(L)`, `V(H)`, and `V(J)`; equation (3) identifies
the fourth as `V(G)`.  Taking the closure of the finite union in (19) gives
(5).  Pairwise coprimality and absolute irreducibility give exactly four
irreducible components.  Equation (6) is the already-audited discriminant
consequence of THM-3498.

## 8. Exact scope

This theorem never expands the global polynomial `G`.  It does not determine
an integral primitive normalization, global term count, multidegree, or the
geometry of singularities of `V(G)`.  It proves neither a depth-five step nor
an all-level induction.  It concerns one fixed map in dimension three and
adds no Jacobian-conjecture, Dixmier-conjecture, arbitrary-family, or Lonely
Runner conclusion.

Reproduce the exact gate with

```text
python 04-computation/keller_level_four_norm_J_mod1009_slice_probe_20260816.py
python -O 04-computation/keller_level_four_norm_J_mod1009_slice_probe_20260816.py
```

**QED.**
