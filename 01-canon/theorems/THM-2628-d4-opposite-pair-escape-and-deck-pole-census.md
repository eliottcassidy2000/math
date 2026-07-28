---
id: THM-2628
title: "D4 opposite-pair escape and deck-pole census"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At the geometric
  generic point of a Jelonek component of a hypothetical degree-four planar
  Keller map with D4 monodromy, let S be the finite affine inverse branches
  and let tau be the canonical opposite-sheet deck involution.  Every branch
  in S is fixed pointwise by local inertia.  The geometric affine branches
  carrying a pole of tau are exactly P=S minus tau(S).  The complete census
  is: k=1 gives one pole branch and identity or diagonal-reflection inertia;
  an adjacent k=2 pair gives two pole branches and identity inertia; an
  opposite k=2 pair gives no pole branch and identity or diagonal-reflection
  inertia; k=3 gives one pole branch and identity inertia.  THM-2612 therefore
  forces at least one Jelonek component of one of the three pole-positive
  types.  THM-2621 identifies the full coefficient pole a_k there, but k and
  its pole order do not determine the opposite-pair ownership.  Counts are
  geometric branch counts, not numbers or orders of irreducible pole divisors.
  This theorem alone gives no D4 exclusion; THM-2633 later excludes the whole
  D4 Keller lane and makes the k=0 rows nonrealizable there.  No component
  count, JC(2), or DC(2) follows.
source: root-long-frontiers-2026-07-28-d4-opposite-pairs
depends_on:
  - THM-2612-d4-deck-pole-tax-and-depressed-resolvent-gcd-gate
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2627-d4-jelonek-quadratic-character-rank-and-component-gate
  - THM-2633-keller-point-stabilizer-abelianization-gate-and-d4-exclusion
script: 04-computation/jacobian_d4_opposite_pair_escape_thm2628.py
output: 05-knowledge/results/jacobian_d4_opposite_pair_escape_thm2628.out
script_sha256: 43bdd6faeeebf52d7298cfe18f5891c348a124ab2f2a4825ad6e71a6312e9ce9
output_sha256: 31bbf8e936a06d3013057302255e96067568e8363fcbb06b7a3775f3ab5bef40
hash_basis: working-tree bytes (LF)
---

# THM-2628 -- a D4 deck pole is a mixed opposite-pair boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2612 proves that the unique `D_4` source-deck involution must acquire a
pole somewhere over the Jelonek boundary.  THM-2621 separately records how
many inverse sheets remain finite at a generic boundary component and which
normalized inverse coefficient has the resulting pole.  Neither result says
which surviving branch is paired with which escaping branch.

That missing coordinate is finite and exact.  The four inverse sheets form
two canonical opposite pairs under the deck involution.  A deck pole occurs
precisely when one member of such a pair remains in the affine source and the
other goes to infinity.  This theorem classifies every possible mixture and
shows why the specialized degree and coefficient pole order cannot recover
it.

## 1. Transverse branches at a Jelonek component

Let

```text
F=(P,Q):X=A^2_C -> Y=A^2_C
```

be a hypothetical Keller map of generic degree four and geometric monodromy
`D_4`.  Use THM-2612 notation

```text
A_F=the Jelonek set,             V=Y minus A_F,

U=F^(-1)(V),                    F|_U:U -> V finite etale.   (1)
```

Fix an irreducible component `D` of `A_F` and pass to a geometric generic
smooth point, a transverse strict henselian DVR, and then its punctured
fraction field.  Label the four geometric inverse branches by

```text
Omega={0,1,2,3}.                                           (2)
```

Let `I_D` be the cyclic local inertia subgroup in the transitive square
action of `D_4`.  Define

```text
S_D={branches having a finite centre in X over D},

k_D=|S_D|.                                                (3)
```

The equality in (3) agrees with THM-2621's generic surviving-fibre
cardinality.  A finite branch converges to a point of `X`; the inverse
function theorem supplies a unique local inverse through that point because
`F` is Keller.  Distinct branches therefore have distinct finite limits, and
each extends across the transverse disk.

This proves the pointwise inertia constraint

```text
S_D subseteq Fix(I_D).                                    (4)
```

The converse is false.  An inertia-fixed branch can be a single-valued
meromorphic branch which still escapes to infinity.  This is the
Zariski-main ownership sidecar missed by monodromy alone.  In particular,
one must not replace (4) by equality.

Finally

```text
0<=k_D<=3.                                                (5)
```

The zero row is a valid local normalization possibility, but THM-2633 proves
that it cannot occur on a target divisor of an affine-space Keller map.  In
the polynomial Keller scope one therefore has `1<=k_D<=3`.

If all four branches remained finite generically, the four local inverses
would account for the full degree in a neighborhood of `D`, leaving no
nonproper branch; `D` would not be a Jelonek component.

## 2. The opposite-pair deck involution

In the square action, THM-2612's nontrivial source-deck transformation is the
central opposite pairing

```text
tau=(0 2)(1 3).                                           (6)
```

For a surviving branch `s in S_D`, let `E_s` denote its geometric affine
prime over `D`.  The pullbacks `tau^*x,tau^*y` along `E_s` are exactly the
coordinate expansions of the opposite branch `tau(s)`.  Hence

```text
s has no deck pole
 iff tau(s) in S_D,

s has a deck pole
 iff tau(s) notin S_D
 iff min(v_(E_s)(tau^*x),v_(E_s)(tau^*y))<0.              (7)
```

In THM-2621's monic transverse frame, bounded `x`, bounded target
coordinates, and monicity of `P` in `y` force `y` bounded.  Thus (7) can be
tested there by the primitive coordinate alone:

```text
s has a deck pole iff v_(E_s)(tau^*x)<0.                  (8)
```

Define the geometric pole-branch set and count

```text
P_D=S_D minus tau(S_D),                 p_D=|P_D|.         (9)
```

Equation (7) proves that this is not merely a combinatorial proxy: it is the
exact branchwise support of the affine pole divisor over `D`.

## 3. Complete survivor/inertia/pole census

Since `tau` partitions `Omega` into two opposite pairs, (4) and (9) give the
complete table.

| `k_D` | shape of `S_D` | `p_D` | possible local inertia |
|---:|---|---:|---|
| `0` | empty | `0` | any cyclic `D_4` inertia |
| `1` | singleton | `1` | identity or the diagonal reflection fixing its opposite pair |
| `2` | adjacent pair | `2` | identity only |
| `2` | opposite pair | `0` | identity or the diagonal reflection fixing that pair |
| `3` | complement of a singleton | `1` | identity only |

Indeed, the identity fixes all subsets.  Every nonidentity `D_4` element
with a fixed vertex is a diagonal reflection, and its fixed set is exactly
one opposite pair.  Because finite local inverses are fixed **pointwise**, an
adjacent pair or a three-point set permits no nonidentity inertia.  Equations
(6) and (9) then give the displayed `p_D` values.

There is one descent guardrail.  The geometric sets `S_D,P_D` must be stable
under the decomposition/residue action.  Thus `p_D` counts geometric pole
branches, equivalently the sum of residue degrees of affine pole primes.  It
is not in general the number of irreducible pole divisors, and it says
nothing about their pole orders.  For diagonal-reflection inertia,

```text
N_(D4)(I_D)/I_D=C_2                                      (10)
```

swaps the two fixed sheets.  A singleton `S_D` therefore requires trivial
residual pairing; a nontrivial residual action forces either both fixed
sheets or neither.

## 4. Global D4 consequence

THM-2612 forces at least one affine prime on which `tau` has a pole.  Its
image lies in some Jelonek component `D`, and (7)--(9) give `p_D>0`.
Therefore every hypothetical `D_4` Keller map has at least one boundary
component of exactly one of the following three types:

```text
k_D=1, singleton,                    p_D=1,

k_D=2, adjacent pair,                p_D=2,

k_D=3, complement of singleton,      p_D=1.               (11)
```

In particular, the finite branch set cannot be a union of complete opposite
pairs at every Jelonek component.  A rotational inertia element, an edge
reflection, or the central half-turn has no finite fixed branch and cannot
carry the affine deck-pole prime forced by THM-2612.  If the pole component
has nontrivial inertia, the table sharpens it uniquely to

```text
k_D=1 with diagonal-reflection inertia.                   (12)
```

This is a boundary-ownership gate, not a `D_4` exclusion: all three types in
(11) are locally consistent.

## 5. Interface with the inverse-coefficient pole ledger

Write THM-2621's transverse resultant as

```text
R=cT^4+r_3T^3+r_2T^2+r_1T+r_0,            f=R/c.         (13)
```

If `e_D=nu_D(c)`, that theorem gives

```text
nu_D(a_(k_D))=-e_D.                                      (14)
```

Combining (14) with the census gives

| survivor type | forced normalized coefficient pole | deck-pole branches |
|---|---|---:|
| singleton `k=1` | `a_1` | `1` |
| adjacent `k=2` | `a_2` | `2` |
| opposite `k=2` | `a_2` | `0` |
| triple `k=3` | `a_3` | `1` |

Thus odd `k_D` detects one mixed opposite pair, but the even value `k_D=2`
is ambiguous.  The same specialized degree, the same leading-coefficient
order, and the same forced pole in `a_2` can describe either two deck-pole
branches or none.  Opposite-pair ownership is an independent sidecar, not a
function of the resultant coefficient ledger.

The local inertia cycle type in this theorem is likewise distinct from
THM-2612's monic integral depressed-order tax.  Importing that gcd gate at a
Jelonek prime requires a separate integral model; the normalized inverse
quartic here already has coefficient poles.  No gcd support is inferred from
the present table alone.

## 6. Exact local hostiles

Three strict-henselian normalization/resultant models show every delicate
boundary.  They are local algebra controls, not global Keller maps.

First,

```text
R_1(T)=t^2 T(T-t^(-1))(T^2-t^(-1))
      =t^2T^4-tT^3-tT^2+T.                              (15)
```

Its roots are

```text
0, t^(-1), +t^(-1/2), -t^(-1/2).                        (16)
```

Inertia swaps the square-root pair and fixes the first two opposite sheets.
Take `tau` to pair `0` with `t^(-1)` and the two square roots with each other.
Only `0` is finite, so

```text
k_D=1,              p_D=1,              R_1 mod t=T.     (17)
```

This is the diagonal-reflection singleton allowed in the table.  It also
exhibits explicitly why `S_D` can be a proper subset of `Fix(I_D)`.

Second,

```text
R_2(T)=(tT-1)(tT-2)T(T-1)                               (18)
```

has identity inertia, two finite roots `0,1`, two escaping roots
`t^(-1),2t^(-1)`, specialized degree two, and leading order two.  Pair the
two finite roots together and the two escaping roots together to get
`p_D=0`; pair across the boundary twice to get `p_D=2`.  The polynomial
`R_2`, `k_D=2`, and `e_D=2` are identical in both cases.  This is the sharp
hostile to recovering deck ownership from (13)--(14).

Third,

```text
R_3(T)=(tT-1)T(T-1)(T-2)                                (19)
```

has identity inertia and three finite roots.  Pair `t^(-1)` with `0` and
pair `1` with `2`; then

```text
k_D=3,              p_D=1,

R_3 mod t=-T(T-1)(T-2).                                  (20)
```

This realizes the triple row of the census.

## 7. Exact verification and scope

Run

```text
python3 04-computation/jacobian_d4_opposite_pair_escape_thm2628.py
python3 -O 04-computation/jacobian_d4_opposite_pair_escape_thm2628.py
```

Both executions byte-match the stored transcript.  The companion:

1. enumerates all eight square permutations in `D_4` and checks that `tau`
   is the central opposite pairing;
2. exhausts all sixteen survivor subsets, their pole counts, and every
   pointwise-fixing inertia element;
3. checks the exact three pole-positive shapes and the full labelled census;
4. expands (15), (18), and (19), verifies specialized degrees `1,2,3`, and
   checks the full normalized coefficient pole orders; and
5. retains explicit optimized-mode guards for every truth-bearing check.

The independent hostile audit rederived the transverse local-inverse
argument, the valuation equivalence (7), every group-theoretic row, and the
decomposition/residue guardrail.  It supplied (15) and (19) as positive local
models and (18) as the same-resultant ownership hostile.

The theorem itself proves no lower bound on the number of Jelonek components,
pole order beyond THM-2621, integral gcd gate, global compactification model,
JC(2), or DC(2).  THM-2633 independently excludes `D_4`, so this census is now
a conditional/local boundary atlas rather than a live Keller lane.  It
isolates the exact missing coordinate:

```text
surviving-sheet degree + coefficient pole
       + opposite-pair boundary ownership.                (21)
```

**QED.**
