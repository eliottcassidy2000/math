---
id: THM-3444
title: "Commuting smooth Hensel vector fields give a free lattice action"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Pointwise independence
  of the first-carry vector fields
  of r commuting congruence-to-identity automorphisms gives a free
  (Z/p^(a-c)Z)^r action on every smooth depth-c Hensel fibre, with an exact
  orbit-bank invoice and first-lift converse.
source: root-smooth-hensel-lattice-action-2026-08-15
audit: independent mixed-valuation proof reconstruction; 729-exponent nonlinear stress; dependent and two-adic hostile audits; normal/optimized/stored replay; dependency, hash, AST, ID, routing, and documentation gates clean
depends_on:
  - THM-3442-smooth-hensel-fibre-vector-field-orbit-law
related:
  - THM-3439-near-identity-grassmannian-hensel-orbit-law
script: 04-computation/commuting_smooth_hensel_lattice_action_thm3444.py
output: 05-knowledge/results/commuting_smooth_hensel_lattice_action_thm3444.out
script_sha256: 7bf18fd6329ce80801f0712d3df37636dfec53154520b04824ecd2342407491c
output_sha256: 7f8587e4bf37ba932ecf6b5f0748284088cb58a7ccac07df8d1feaea3ae197b7
semantic_sha256: a3393fafed345ae96bea50a5dcc4dfb7dadb56320e1d7cefd57a800d3ed615c1
hash_basis: LF-normalized bytes
---

# THM-3444 -- commuting smooth Hensel vector fields give a free lattice action

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The proof and exact companion have passed an independent derivation and
replay audit.

## 1. Exact statement

Let `p` be prime, let `c>=1`, and assume

```text
p is odd, or c>=2.                                      (1)
```

Let `X` be a smooth finite-type `Z_p`-scheme of pure relative dimension
`d>=1`, with `X(F_p)` nonempty.  Let

```text
g_1,...,g_r in Aut_(Z_p)(X),       r>=1,                (2)
```

be pairwise commuting automorphisms, each scheme-theoretically equal to the
identity modulo `p^c`.  Write

```text
delta_i in H^0(X_(F_p),T_X)                             (3)
```

for their first-carry vector fields from THM-3442.  Assume that for every
`xbar in X(F_p)` the vectors

```text
delta_1(xbar),...,delta_r(xbar)                         (4)
```

are linearly independent over `F_p`.  In particular `r<=d`.

For every `a>=c`, put `n=a-c`.  The commuting action on a depth-`c` fibre
factors through

```text
G_(a,c)=(Z/p^n Z)^r.                                    (5)
```

On every nonempty fibre of

```text
X(Z/p^a Z) -> X(Z/p^c Z),                               (6)
```

this action is free.  Consequently

```text
fibre size       =p^(d n),
every orbit size =p^(r n),
number of orbits =p^((d-r)n).                           (7)
```

For `a>c`, the action is transitive on every fibre exactly when `r=d`.

The independence condition is exact at the first lift.  The following are
equivalent:

```text
the fields in (4) are independent at every F_p-point;
(F_p)^r acts freely on every depth-c fibre at level c+1;
G_(a,c) acts freely on every depth-c fibre for every a>c.                 (8)
```

More locally, if `(4)` is dependent at one `xbar`, choose a nonzero
`q=(q_i) in F_p^r` with `sum q_i delta_i(xbar)=0`.  For every lift `x_c` of
`xbar`, the nonidentity element `q` fixes all `p^d` first lifts above `x_c`.

The same statement holds for a smooth `p`-adic formal scheme locally
topologically of finite presentation and of pure relative dimension `d`,
with finite residue-point set.

## 2. Inheritance and connection ledger

THM-3442 is the one-generator law.  The new point is not an application of
that theorem to each generator separately: freeness of the product action
requires every nonzero linear combination of the first carries to remain
nonzero.  Pointwise linear independence is exactly that missing coordinate.

The canonical hostile is a redundant translation pair.  Two individually
free generators can have the same tangent direction, so a nontrivial group
combination fixes an entire fibre.  The corrected near miss is therefore
“each `delta_i` is nowhere zero”; this is necessary but not sufficient.  The
least-used sidecar is the `d-r`-dimensional orbit-bank label.

| field | exact content |
|---|---|
| source | `r` commuting smooth Hensel translations at one common depth |
| target | a free finite `p`-adic lattice action on each lifting fibre |
| map | exponent vector `(m_i)` to `product_i g_i^(m_i)` |
| preserved | base point modulo `p^c` and the full exponent lattice |
| destroyed by one orbit | `d-r` transverse tangent coordinates |
| required sidecar | one of `p^((d-r)(a-c))` orbit-bank labels |
| cheapest positive | independent coordinate translations on `A^d` |
| cheapest hostile | two translations along the same coordinate |

## 3. The finite action

THM-3442 gives, for every `i`,

```text
g_i^(p^n)=id mod p^(c+n).                               (9)
```

Thus each generator has exponent dividing `p^n` on level `a=c+n`.  Pairwise
commutation makes `(5)` a well-defined action, and every generator preserves
each fibre of `(6)`.  Smoothness gives exactly `p^(dn)` points in that fibre.

## 4. No nontrivial stabilizer

Take a nonzero exponent vector

```text
m=(m_1,...,m_r) in (Z/p^n Z)^r.                         (10)
```

Choose representatives and let

```text
t=min_i v_p(m_i)<n,       m_i=p^t q_i,                  (11)
```

where at least one `q_i` is a unit modulo `p`.  Put

```text
h=product_i g_i^(m_i).                                  (12)
```

The first-carry iteration in THM-3442, followed by multiplication of the
commuting factors, gives on a local formal chart

```text
h^*=I+p^(c+t)D_m mod p^(c+t+1),
D_m mod p=sum_i (q_i mod p) delta_i.                    (13)
```

All cross terms contain at least `p^(2(c+t))` and lie beyond the displayed
layer under `(1)`.  The right side of `(13)` is nonzero at every special
point by `(4)`.

If a point at level `a` were fixed by `h`, evaluation of `(13)` on local
functions and cancellation of `p^(c+t)` would force the vector on the right
to vanish at its reduction.  This is impossible.  Hence no nonzero element
of `(5)` fixes any point, proving freeness.  Dividing the smooth fibre count
by `|G_(a,c)|=p^(rn)` proves `(7)` and the transitivity boundary.

## 5. The exact first-lift converse

Fix `xbar` and a lift `x_c`.  The first lifting fibre is the affine torsor

```text
x_c+p^c T_(xbar)X.                                      (14)
```

By THM-3442, the element `q=(q_i) in (F_p)^r` translates `(14)` by

```text
sum_i q_i delta_i(xbar).                                (15)
```

An action of an elementary abelian `p`-group by translations is free exactly
when its translation map is injective.  Thus `(4)` is equivalent to
first-level freeness.  If dependence occurs, a nonzero kernel vector makes
`(15)` zero and fixes the whole torsor.  Together with Section 4, this proves
all directions of `(8)`.

## 6. Equality and failure boundaries

1. **Affine equality family.**  On `A^d`, translations by `p^c` in the first
   `r` coordinate directions attain every formula in `(7)`.  The action is
   regular exactly for `r=d`.
2. **Nonlinear equality family.**  Conjugate two coordinate translations on
   `A^2` by `H(x,y)=(x,y+x^2)`.  The resulting commuting maps have first
   carries `(1,2x)` and `(0,1)`, determinant one at every special point, and
   act regularly on every depth-`c` fibre.
3. **Redundant hostile.**  The translations

   ```text
   (x,y)->(x+p^c,y),       (x,y)->(x+2p^c,y)
   ```

   are individually free, but their first carries are dependent.  A
   nonzero exponent relation fixes every point.
4. **Two-adic boundary.**  At `p=2,c=1`, take the THM-3442 irreducible
   projective first carry on `P^1` and an independent affine translation.
   The two fields are independent on `P^1 x A^1`, but the projective
   generator has order two rather than four at level three, so the expected
   rank-two action is not free.  At `c=2` the full `4x4` action is regular.
5. **Commutation.**  Without pairwise commutation, exponent vectors do not
   define the abelian group `(5)` and commutator carries can enter later
   layers.  No nilpotent-group extension is asserted.
6. **No physical transport.**  Orbit lattices are not LRC owners, currents,
   or safe-time intervals.  No LRC(14), `JC(2)`, or boundary-response
   consequence follows.

## 7. Exact companion

The standard-library companion exhausts finite fibres for affine lattice
actions with `r<d` and `r=d`, all nine base fibres of the nonlinear conjugate
translation at `p=3`, a dependent first-lift hostile, and the product
projective/affine `p=2` boundary at depths one and two.  It checks point,
group, orbit, stabilizer, and transitivity invoices using integer arithmetic.

Reproduce with

```text
python3 -B 04-computation/commuting_smooth_hensel_lattice_action_thm3444.py
python3 -B -O 04-computation/commuting_smooth_hensel_lattice_action_thm3444.py
```

QED.
