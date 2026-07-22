---
id: THM-2074
title: "Strict LRC(14) holds on a density-one set of thirteen-speed rows"
status: >
  PROVED. THM-2051 confines every thirteen-speed row without a
  positive-measure strict LRC(14) witness to one of exactly
  25173854387233097811887443361297472 primitive support-three-to-five
  relation hyperplanes. Their union has at most R B^12 points in the ordered
  height-B box, while increasing primitive distinct rows have asymptotic
  B^13/(13! zeta(13)). Strict LRC(14) therefore has density one. Moreover,
  after adjoining the 78 equality hyperplanes, every prescribed prime-power
  tower contains whole congruence packets all of whose positive lifts are
  certified. This is an almost-everywhere theorem, not universal LRC(14).
source: codex-2026-07-21-LRC-density-one-hyperplane-sieve
depends_on:
  - THM-2051
related:
  - THM-934
  - THM-2065
script: 04-computation/lrc14_density_one_relation_hyperplanes_codex_20260721.py
result: 05-knowledge/results/lrc14_density_one_relation_hyperplanes_codex_20260721.out
script_sha256: 71e01a8f2e3e6c36148e247a9afb81cc8624af0acedfb1cff12b45a91070c388
result_sha256: 14369f4098b384bfe31b63a424a52767eeecdcfe1c3ea4f8b9fe003a18402fdb
hash_basis: normalized repository blobs (LF)
---

# THM-2074 -- Density-one strict LRC(14)

Put

```text
n=13,                    H=2^20.                         (1)
```

A speed row is an ordered tuple `v=(v_1,...,v_13)` of distinct positive
integers. It is primitive when `gcd(v_1,...,v_13)=1`. Scaling does not change
the Lonely Runner predicate, but keeping primitivity makes the density
normalization canonical.

## 1. The finite relation ledger

Let `A_H` contain one representative modulo global sign of every primitive
integer normal

```text
a=(a_1,...,a_13),
3<=|supp(a)|<=5,             0<|a_i|<=H on its support. (2)
```

Write `R=|A_H|`. THM-2051 Theorem B gives the exact implication

```text
no positive-measure strict LRC(14) witness for v
    => a.v=0 for some a in A_H.                          (3)
```

Dividing the relation coefficients by their gcd preserves the height bound,
so (2) loses no relation. Distinct primitive integer normals define the same
rational hyperplane only when they differ by sign. Thus `A_H` is the exact
hyperplane ledger, not merely a coefficient-vector upper bound.

For `s>=1`, define

```text
P_s(H)=sum_(d=1)^H mu(d) floor(H/d)^s.                  (4)
```

Mobius inversion says that `P_s(H)` counts positive `s`-tuples of coefficient
magnitudes with gcd one. Choosing the support and then a sign pattern modulo
global reversal gives

```text
R=sum_(s=3)^5 C(13,s) 2^(s-1) P_s(H).                  (5)
```

At `H=2^20`, the exact values are

```text
P_3=959124025074311215,
P_4=1116973047989955380768527,
P_5=1222506215916189106034284205191,                   (6)

R=25173854387233097811887443361297472.                 (7)
```

For comparison, retaining every signed coefficient vector before primitive
normalization gives the valid crude bound

```text
sum_(s=3)^5 C(13,s)(2H)^s
 =52206936149913413947000523914739712.                 (8)
```

## 2. Density-one theorem

For any nonzero normal `a`, the number of points `v in {1,...,B}^13` with
`a.v=0` is at most `B^12`: fix the other twelve coordinates and solve for one
coordinate whose coefficient is nonzero. The union bound and (3) therefore
give the completely explicit estimate

```text
#{ordered height-B rows with no positive-measure strict witness}
    <=R B^12.                                           (9)
```

This remains true after restricting to distinct or primitive rows.

On the other hand, Mobius inversion gives

```text
#{primitive positive ordered rows in {1,...,B}^13}
    =B^13/zeta(13)+O(B^12).                            (10)
```

The union of the `C(13,2)=78` equality hyperplanes has `O(B^12)` points, so
distinctness does not change the leading term. Every distinct ordered row has
exactly one of its `13!` permutations increasing. Hence

```text
#{1<=v_1<...<v_13<=B, gcd(v)=1}
    =B^13/(13! zeta(13))+O(B^12).                      (11)
```

Dividing (9) by (11) proves:

> Among increasing primitive thirteen-speed rows of height at most `B`, the
> proportion having a positive-measure strict LRC(14) witness tends to one as
> `B->infinity`.

Equivalently, the exceptional relative density is `O(R/B)`. This upgrades
THM-934's `200/200` random experiment to an asymptotic theorem. It does not
classify the relation-rich zero-density set containing AP and covering
extremals.

## 3. Whole congruence packets on any fixed prime tower

Adjoin the 78 primitive equality normals `e_i-e_j` to `A_H`, and put

```text
H_tot=R+78
     =25173854387233097811887443361297550.              (12)
```

Fix **any** prime `ell`, any `k>=1`, and `q=ell^k`. A primitive integer normal
has at least one coefficient not divisible by `ell`, hence a unit coefficient
modulo `q`. Its zero set in `(Z/qZ)^13` therefore has exactly `q^12` points:
choose the other twelve coordinates arbitrarily and solve uniquely for the
unit coordinate. The union bound gives at least

```text
q^13-H_tot q^12=q^12(q-H_tot)                           (13)
```

residue rows avoiding every relation and every equality. In particular, the
packet is nonempty whenever `ell^k>H_tot`.

Let `bar v` be any good residue row and let `v` be any positive integer lift.
Distinct residues imply distinct integer entries. If `v` satisfied a relation
allowed by THM-2051, primitive normalization would put its normal in `A_H`,
contradicting goodness modulo `q`. Hence every positive lift has a
positive-measure strict witness by THM-2051.

The lift need not itself be primitive. Dividing by its common gcd preserves
positivity, distinctness, relation-freeness, and the LRC predicate. Thus one
good residue class supplies an infinite certified family.

This is a literal no-new-prime statement: the modulus may remain on any
preassigned prime-power tower. The first guaranteed levels for four small
primes are

| `ell` | first `k` with `ell^k>H_tot` | `ell^k` |
|---:|---:|---:|
| 2 | 115 | 41538374868278621028243970633760768 |
| 3 | 73 | 67585198634817523235520443624317923 |
| 5 | 50 | 88817841970012523233890533447265625 |
| 7 | 41 | 44567640326363195900190045974568007 |

This is characteristic-zero congruence selection for integer speeds. It is
not an `F_ell` torsion theorem: in characteristic `ell`, `X^(ell^k)-1` is
`(X-1)^(ell^k)` and the reduced torsion collapses to the Hasse-jet regime of
THM-2041/2043.

## 4. Exact-period witnesses after a strict exit

There is a complementary phase-space corollary. If a fixed row already has a
strict witness, continuity supplies an open interval of strict witnesses. For
any preassigned prime `ell`, primitive points `a/ell^k` become dense as
`k->infinity`, so all sufficiently large levels contain an exact-period
`ell^k` strict witness. Again no new denominator prime is needed **after** a
strict exit exists.

The order of quantifiers is essential. This density argument cannot create a
witness on a tight boundary. For example, at threshold `1/4`, the safe set of
`{1,2,3}` consists only of `1/4,3/4`; generic boundary avoidance deletes both.
Likewise THM-2072's common-multiple core defeats every fixed finite LRC clock
bank even though an adaptive continuous certificate exists.

## 5. Relation to the live LRC residual

The theorem cleanly separates bulk from structure:

```text
generic thirteen-speed rows
    -> avoid a finite relation arrangement
    -> positive-measure strict exit (density one);

structured residual
    -> lies in the arrangement
    -> needs owner/deck/phase-height analysis.            (14)
```

On a two-anchor plane, THM-2065 pulls each nonpersistent relation hyperplane
back to one projective ray; only persistent marked coefficient-row circuits
remain. THM-2069 then filters their deletion primitivity by a code/cogirth
wheel. On the dyadic lane, THM-2073--2079 retain safe-child addresses and odd
tail owners. MISTAKE-230 is decisive: an empty full-row safe set does not
descend after those tails are deleted, so the terminal core alone is not the
missing universal theorem.

Thus THM-2074 is real LRC progress but cannot finish LRC(14): all hypothetical
counterexamples already live in its zero-density exceptional arrangement.
The exact remaining work is classification/discharge there, with the marked
sidecars retained.

## 6. Cross-domain lesson and limits

The mechanism is a finite algebraic genericity principle. On any fixed
multivariate Gaussian support whose charge convex hull contains zero, one
balanced occupation contributes a coefficient monomial with a nonzero Wick
factor. The corresponding moment is therefore a nonzero polynomial in the
independent coefficients, so generic coefficient choices are outside the
moment nullcone in every dimension.

That observation does not prove arbitrary-coefficient GMC. A tuned point may
lie on the moment hypersurface; THM-2070 embeds the hostile Laurent polynomial
`u^2+u+u^-1-u^-2`, whose return support is cofinite while every odd constant
term vanishes. Physical angular rotations act trivially on balanced channels,
and a generic coefficient-phase deformation changes the polynomial. THM-2022's
whole-face Frobenius theorem, not hyperplane avoidance, is what controls the
original arbitrary coefficients in dimension two.

The same warning applies to finite reconstruction: a generic scalar functional
can separate a finite atlas of distinct signatures, but an invariant scalar
cannot distinguish equal signatures or recover labels it discarded.

## Assumption challenge and tournament analysis

Runners are not the faithful vertices of this proof. The faithful vertices are
primitive relation normals and equality obligations; their pairwise observable
is intersection codimension, which is symmetric. There is no intrinsic switch
turning it into a directed relation and no meaningful tie Hamiltonian path.
Lexicographically orienting the hyperplanes would encode only file order and
destroy the union-bound predicate. The correct carrier is an undirected
intersection lattice plus coefficient-height and support labels.

## Computational audit

The companion uses runtime checks, not Python assertions. It verifies the
primitive-normal formula on nine small universes, directly enumerates 1,376
small relation hyperplanes and their box union, checks the ordered/increasing
`5!` quotient, and performs 414 exact prime-power hyperplane-size checks on
nine moduli. A linear Mobius sieve then reproduces (6)--(8), (12), and the four
prime-power thresholds. Normal and `python -O` outputs byte-match the frozen
result and end in `PASS`.
