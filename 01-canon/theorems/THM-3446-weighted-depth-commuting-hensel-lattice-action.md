---
id: THM-3446
title: "Weighted-depth commuting Hensel lattice action"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Commuting generators with
  unequal carry depths and pointwise-independent first carries act freely
  through their weighted direct product; the exact orbit-bank exponent
  includes the persistent depth spread, and dependence creates a stabilizer
  by the sharp universal level M+1, though possibly earlier.
source: root-weighted-depth-hensel-lattice-2026-08-15
audit: >
  independent minimum-effective-depth proof reconstruction; exact-depth,
  a=M, tied-carry, orbit-bank, transitivity, delayed-bound, p=2, and
  scheme/formal typing audits; MISTAKE-398 correction-lineage audit;
  normal/optimized/stored and clean-room semantic replay; dependency, hash,
  AST/security, ID, routing, and documentation gates clean
depends_on:
  - THM-3444-commuting-smooth-hensel-vector-field-lattice-action
related:
  - THM-3442-smooth-hensel-fibre-vector-field-orbit-law
script: 04-computation/weighted_depth_commuting_hensel_lattice_thm3446.py
output: 05-knowledge/results/weighted_depth_commuting_hensel_lattice_thm3446.out
script_sha256: 117f9a3a27aacdbc7f3576fb8edd53b35c12e39c35fe5b1d46577b8a708280b2
output_sha256: cec4e30c279d284fc657ca72b781c7163f911947dfef0c9188307d644bb2d2f3
semantic_sha256: 317fbb4d6c074e76a3886a8653f4b6bb6423f97c135cff44c2df212722e85bb2
hash_basis: LF-normalized bytes
---

# THM-3446 -- weighted-depth commuting Hensel lattice action

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The proof and exact companion passed an independent derivation, hostile, and
replay audit after MISTAKE-398 repaired the group and alignment typing.

## 1. Exact statement

Let `p` be prime.  Let `X` be a smooth finite-type `Z_p`-scheme of pure
relative dimension `d>=1`, with `X(F_p)` nonempty.  Let

```text
g_1,...,g_r in Aut_(Z_p)(X),       r>=1,                  (1)
```

be pairwise commuting.  Give generator `g_i` an integer carry depth
`c_i>=1`, meaning that it is scheme-theoretically the identity modulo
`p^(c_i)`.  Assume

```text
p is odd, or min_i c_i>=2.                               (2)
```

Let

```text
delta_i=((g_i^*-I)/p^(c_i)) mod p                       (3)
```

be its first-carry vector field.  It is allowed that `g_i` is actually
deeper, but then `delta_i=0` and the independence hypothesis below cannot
hold.  Put

```text
C=min_i c_i,       M=max_i c_i.                         (4)
```

Assume that at every `xbar in X(F_p)` the vectors

```text
delta_1(xbar),...,delta_r(xbar)                         (5)
```

are linearly independent over `F_p`.

This independence hypothesis itself forces `r<=d`; no separate rank bound is
needed in the ambient statement.  Keeping arbitrary `r` also types the
dependent converse and its cheapest hostile inside the same universe.

For every `a>=M`, the commuting action on each depth-`C` fibre factors
through the weighted exponent group

```text
G_a=product_(i=1)^r Z/p^(a-c_i)Z.                       (6)
```

The action of `G_a` is free.  Hence every nonempty fibre of

```text
X(Z/p^a Z) -> X(Z/p^C Z)                                (7)
```

has

```text
fibre size       =p^(d(a-C)),
every orbit size =p^(sum_i(a-c_i)),
number of orbits =p^B(a),                               (8)

B(a)=d(a-C)-sum_i(a-c_i)
    =(d-r)(a-C)+sum_i(c_i-C).                           (9)
```

For `a>C`, the action is transitive exactly when

```text
r=d and c_1=...=c_r=C.                                  (10)
```

Thus unequal depths create a persistent orbit-bank tariff even when `r=d`.

The independence condition is also necessary for freeness at all levels.
If `(5)` is dependent at some `xbar`, choose nonzero
`q=(q_i) in F_p^r` with `sum_i q_i delta_i(xbar)=0`.  At level

```text
a=M+1,       m_i=p^(M-c_i) q_i                          (11)
```

defines a nonzero element of `G_(M+1)`.  For every lift `x_M` of `xbar`,
this element fixes all `p^d` points over `x_M` at level `M+1`.  Therefore

```text
(5) holds at every F_p-point
iff G_a acts freely on every depth-C fibre for every a>=M.                (12)
```

When all depths agree, `(8)--(12)` reduce to THM-3444 and the dependence is
visible already on the first lift.  At unequal depths it can remain hidden
until the carries are aligned at depth `M`.

The theorem also holds for a smooth `p`-adic formal scheme locally
topologically of finite presentation and of pure relative dimension `d`,
with finite residue-point set.

## 2. Inheritance and connection ledger

THM-3444 treats one common carry depth.  Its proof suggests rescaling every
generator to the deepest common layer, but that operation destroys the
actual exponent orders and the orbit-bank count.  The weighted product `(6)`
retains exactly that missing scale coordinate.

The canonical hostile is a shallow and a deep translation along the same
coordinate.  The deep generator is invisible at the shallow first lift, but
at depth `M+1` a powered shallow translation cancels it.  The corrected near
miss is therefore to audit only the minimum-depth fields.  All fields must be
retained and aligned.  The least-used sidecar is the depth vector
`(c_i-C)_i`.

| field | exact content |
|---|---|
| source | commuting smooth automorphisms at unequal carry depths |
| target | a direct product of cyclic exponent channels acting freely |
| map | `(m_i)` to `product_i g_i^(m_i)` |
| preserved | each generator's individual exponent depth |
| destroyed by flattening depths | the persistent tariff `sum(c_i-C)` |
| required sidecar | the full depth vector, not only `C` or `M` |
| cheapest positive | coordinate translations by `p^(c_i)` |
| cheapest hostile | translations by `p^C` and `p^M` in one direction |

## 3. Factor-through and the weighted invoice

Apply THM-3442 to generator `i` at its own depth `c_i`.  At level `a`,

```text
g_i^(p^(a-c_i))=id mod p^a.                             (13)
```

Pairwise commutation therefore gives the action `(6)`.  Every generator is
the identity modulo `p^C`, so it preserves `(7)`.  Smoothness gives the fibre
cardinality in `(8)`.  Formula `(9)` is the quotient of that cardinality by
`|G_a|`, once freeness is proved.

## 4. Mixed-depth stabilizer exclusion

Take a nonzero `m=(m_i) in G_a`, using representatives

```text
0<=m_i<p^(a-c_i).                                       (14)
```

With `v_p(0)=infinity`, put

```text
L=min_i(c_i+v_p(m_i))<a,
I_L={i:c_i+v_p(m_i)=L},
m_i=p^(L-c_i)q_i       for i in I_L.                    (15)
```

At least one `q_i` is a unit.  The individual first-carry iteration from
THM-3442 and multiplication of the commuting factors give

```text
(product_i g_i^(m_i))^*
  =I+p^L D_m mod p^(L+1),
D_m mod p=sum_(i in I_L)(q_i mod p) delta_i.             (16)
```

Components of larger effective depth vanish from this layer, while every
cross term is divisible by `p^(2L)`.  Condition `(2)` is exactly what makes
the individual iteration valid; `2L>=L+1` handles the product.

By `(5)`, the vector on the right of `(16)` is nonzero at every special
point.  A fixed point of the group element at level `a` would, after
evaluation on local functions and cancellation of `p^L`, force this vector
to vanish.  Thus no nonzero element has a fixed point.  This proves freeness,
then `(8)--(9)`.  Both summands in the last form of `(9)` are nonnegative, so
its equality condition is exactly `(10)`.

## 5. Universal delayed-dependence bound and sharpness

Suppose `(5)` fails at `xbar` and choose `q` as above.  The exponent vector
in `(11)` is nonzero because each component with `q_i!=0` has valuation
exactly one below its modulus at level `M+1`.  Every active component now has
effective depth

```text
c_i+v_p(m_i)=M.                                         (17)
```

Thus its product is the identity modulo `p^M`, and its first carry at depth
`M` is `sum_i q_i delta_i`, which vanishes at `xbar`.  The first-lift
translation law of THM-3442 says that it fixes the entire tangent torsor over
every lift `x_M`.  Smoothness supplies those lifts.  This proves the reverse
direction of `(12)` and the universal alignment bound `M+1`.  This bound is
sharp, but it need not be the first failure level for every dependent packet:
a relation supported on shallower equal-depth generators can appear earlier.

## 6. Equality and failure boundaries

1. **Weighted affine equality.**  On `A^d`, translate coordinate `i` by
   `p^(c_i)` for `i<=r`.  Every formula `(8)--(9)` is attained.
2. **Persistent bank at full rank.**  For `d=r=2` and depths `(1,2)`, the
   depth-one fibre at level three has `p^4` points and the exponent group has
   size `p^3`, leaving exactly `p` orbits.  Full tangent rank does not imply
   transitivity after the scale vector is discarded.
3. **Delayed hostile.**  On `A^1`, translations by `p` and `p^2` are
   dependent.  The shallow action is free at level two, where the second
   factor is trivial.  At level three, `g_1^p g_2^(-1)` fixes every point.
   This attains the universal `M+1` bound.  Here `r=2>d=1`, which is allowed
   in the ambient dependent-converse universe; the positive independence
   hypothesis automatically rules it out.
4. **Nonlinear equality.**  Conjugating unequal-depth coordinate
   translations on `A^2` by `(x,y)->(x,y+x^2)` preserves commutation and
   independence; the exact weighted bank remains.
5. **Two-adic range.**  All depths at least two are admitted.  A depth-one
   factor at `p=2` is excluded for the same quadratic-carry reason as in
   THM-3442/3444.
6. **No untyped scale merge.**  Replacing all `c_i` by `C` changes both the
   abstract group and the bank exponent.  It is not a harmless normalization.
7. **No physical transport.**  The weighted depth vector is not an LRC
   denominator profile or a boundary-response DeathBar.  No LRC(14),
   `JC(2)`, or physical consequence follows.

## 7. Exact companion

The standard-library companion exhausts affine fibres with equal and unequal
depths, ranks below and equal to dimension, a three-generator scale vector,
all nine nonlinear conjugate-translation base fibres, the delayed dependence
hostile before and at its exact failure level, and the surviving
`p=2,c_i>=2` range.  It checks group sizes, orbit sizes/counts, exponent
orders, persistent tariffs, and transitivity using integer arithmetic only.

Reproduce with

```text
python3 -B 04-computation/weighted_depth_commuting_hensel_lattice_thm3446.py
python3 -B -O 04-computation/weighted_depth_commuting_hensel_lattice_thm3446.py
```

QED.
