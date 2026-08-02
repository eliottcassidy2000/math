---
id: THM-3046
title: "Quartic resolvent-root valuation realization of the binary-ternary clutch"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For four labelled distinct roots in a discretely valued splitting field,
  the valuations of the three cubic-resolvent root differences are exactly
  the three opposite-edge valuation sums.  Their parities are THM-3045's
  binary matching clutch, while their total modulo three is the ternary
  scalar and equals half the common quartic/resolvent discriminant valuation
  modulo three.  Thus the two exceptional primes have a literal quartic
  valuation realization.  The theorem does not supply an affine owner,
  a Keller restriction, a canonical root labelling, or an LRC consequence.
  The tropical Pluecker law forces the minimum matching sum to occur at
  least twice, yet all 24 clutch classes occur on honest 5-adic quartics.
source: codex-quartic-resolvent-valuation-clutch-2026-08-01
audit: >
  Independent hostile audit rederived the three indexed root-difference
  products, common discriminant valuation, exact THM-3045 quotient map,
  S4/V4 action, projector integrality equivalence, and prime-local split;
  checked all three rational 5-adic hostiles and the tropical Pluecker law;
  independently validated the CRT construction of all 24 honest quartic
  classes; required the final projected-restriction wording repair; and then
  reproduced normal, optimized, stored output, and both LF hashes exactly.
depends_on:
  - THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
  - THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary
script: 04-computation/quartic_resolvent_valuation_clutch_thm3046.py
output: 05-knowledge/results/quartic_resolvent_valuation_clutch_thm3046.out
script_sha256: 239faeeaa873f3ea892215274d33b81438423ca05033b4f9916ffdecf111f50f
output_sha256: 7daf9acd6fbd0e6c6a9b1dab6e28315cd304967e4c38ebd0f3b68d5da8b537d4
hash_basis: LF-normalized bytes
---

# THM-3046 -- quartic root valuations realize the binary-ternary clutch

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and exact statement

[THM-3045, the integral binary-ternary edge
clutch](THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md) computes
the exact quotient left when the three rational isotypic summands
`1+[22]+[31]` are intersected with the six-edge lattice of `K4`:

```text
F2[three perfect matchings] direct_sum F3_triv.             (1)
```

The theorem below identifies `(1)` on the root geometry of a quartic.

Let `(K,v)` be a discretely valued field, with `v(K*)=Z`, and let
`z_0,z_1,z_2,z_3` be four distinct elements of a valued splitting field on
which the chosen extension of `v` is still integer-valued.  Put

```text
f(T)=product_(i=0)^3 (T-z_i),
x_ij=v(z_i-z_j).                                           (2)
```

Use the three perfect matchings

```text
m_1=01|23,             m_2=02|13,             m_3=03|12, (3)
```

and the standard matching-resolvent roots

```text
u_1=z_0 z_1+z_2 z_3,
u_2=z_0 z_2+z_1 z_3,
u_3=z_0 z_3+z_1 z_2.                                    (4)
```

Then

```text
s_1:=v(u_2-u_3)=x_01+x_23,
s_2:=v(u_1-u_3)=x_02+x_13,
s_3:=v(u_1-u_2)=x_03+x_12.                              (5)
```

Consequently the THM-3045 clutch of the edge-valuation vector `x` is

```text
kappa=(s_1,s_2,s_3) mod 2,
tau  =(s_1+s_2+s_3) mod 3.                              (6)
```

If `R(U)=product_i(U-u_i)`, then

```text
Disc(R)=Disc(f),
v(Disc(f))=2(s_1+s_2+s_3),
tau=(v(Disc(f))/2) mod 3.                               (7)
```

The half in `(7)` is an integer in the fixed splitting-field valuation: the
discriminant is the square of the root Vandermonde.  Equations `(5)--(7)`
are relative to that chosen valuation normalization; no invariance under a
later ramification rescaling is claimed.

## 2. Root-difference proof

Direct expansion gives

```text
u_2-u_3=(z_0-z_1)(z_2-z_3),
u_1-u_3=(z_0-z_2)(z_1-z_3),
u_1-u_2=(z_0-z_3)(z_1-z_2).                             (8)
```

Every factor is nonzero because the quartic roots are distinct.  Applying
`v(ab)=v(a)+v(b)` proves `(5)`.  Multiplying the three identities in `(8)`
uses every edge of `K4` exactly once, up to sign.  Squaring therefore gives

```text
product_(i<j)(u_i-u_j)^2
  =product_(i<j)(z_i-z_j)^2,                            (9)
```

which is `(7)`.  This is the root-level mechanism behind the exact
quartic/resolvent discriminant identity in THM-2455 and THM-2598.

### Tropical restriction and full clutch realizability

The three products on the right side of `(8)` satisfy the Pluecker relation

```text
(z_0-z_1)(z_2-z_3)
 -(z_0-z_2)(z_1-z_3)
 +(z_0-z_3)(z_1-z_2)=0.                                (10)
```

The ultrametric inequality therefore forces the minimum of
`s_1,s_2,s_3` to occur at least twice.  This is the projected cluster-tree
restriction on the three matching sums.  The full six-edge valuation vector
obeys additional triangle ultrametric constraints.  The projected
restriction does **not** remove any class from `(1)`.

Indeed, prescribe `(kappa,tau) in F2^3 x F3`.  Choose a special coordinate
so the other two binary entries have a common value `c`; this is always
possible for three bits.  Put `t=c in {0,1}` and choose the unique
`a in {0,...,5}` satisfying

```text
a=kappa_special-c mod 2,                 a=tau mod 3.    (11)
```

Over the `5`-adic valuation, the following rational roots realize matching
sums `(t+a,t,t)`:

```text
t=a=0:       (0,1,2,3);
t=0<a:       (0,5^a,1,2);
0<t, a=0:    (0,5^t,2*5^t,1);
0<t, 0<a:    (0,5^(t+a),5^t,1).                        (12)
```

Relabelling the four roots places the special matching in any of the three
positions.  Since the sum of `(t+a,t,t)` is `3t+a`, `(11)` gives exactly the
prescribed binary parities and ternary total.  Thus all `24` clutch classes
occur on honest separable quartics even after imposing the tropical
Pluecker law.

## 3. Exact identification with the integral quotient

The six integers in `(2)` form an element of the edge lattice
`L=Z[E(K4)]`.  THM-3045's binary coordinate at a matching is the sum of its
two opposite-edge coordinates modulo two, and its ternary coordinate is the
sum of all six edge coordinates modulo three.  Equation `(5)` therefore
identifies the former with `s_i mod 2`, while

```text
sum_(i<j)x_ij=s_1+s_2+s_3                            (13)
```

identifies the latter with `(6)` and `(7)`.

Under relabelling by `S4`, the three entries of `kappa` are permuted through
the matching quotient `S4/V4=S3`, and `tau` is fixed.  Thus `(6)` is not only
an abstract isomorphism of groups: it is the exact equivariant quotient map
of THM-3045 evaluated on the quartic root-divisor vector.

Let `P_0,P_22,P_31` be THM-3045's rational isotypic projectors.  Their values
on `x` are all integral exactly when

```text
s_1,s_2,s_3 are even,             s_1+s_2+s_3=0 mod 3. (14)
```

Equivalently, all three resolvent-root difference valuations are even and,
in their presence, `12` divides `v(Disc(f))`.  At prime `2`, only their
three parity defects remain.  At prime `3`, only the half-discriminant
checksum remains.  At every prime at least `5`, the decomposition is already
integral.

This is the precise sense in which the same quartic object carries separate
binary and ternary tower grammars.  It is not a literal identification of
the Farey tree or a ternary triple tree with the two free factors of
`PSL_2(Z)`.

## 4. Sharp independence controls

Both pieces of `(6)` are necessary, even on honest quartic root packets.
Use the `5`-adic valuation on rational roots.

```text
roots                   (s_1,s_2,s_3)    kappa     tau   v(Disc)
(0,5,1,6)                 (2,0,0)          000       2        4
(0,25,1,26)               (4,0,0)          000       1        8
(0,25,5,1)                (2,1,1)          011       1        8. (15)
```

The last two packets have the same discriminant valuation and the same
ternary scalar but different binary matching data.  The first two have the
same binary matching data but different ternary scalars.  Hence:

- the discriminant valuation cannot recover the three matching parities;
- the matching parities cannot recover the ternary checksum;
- equality `Disc(R)=Disc(f)` alone does not split the integral edge lattice.

These are root-realizable stopping examples, not arbitrary lattice vectors.

## 5. Connection ledger and scope

The exact connection is:

```text
source:       labelled quartic root-difference divisor x in Z[E(K4)];
target:       F2[M] direct_sum F3 from THM-3045;
map:          opposite-pair valuation sums mod 2, total valuation mod 3;
preserved:    S4 relabelling, matching quotient, common discriminant;
destroyed:    root units, residue phases, affine presence, owner labels;
needed sidecar for geometry:
              a physical labelled root/owner realization and compatible
              divisor-to-source transport.                              (16)
```

In particular, this theorem does not constrain a hypothetical degree-four
Keller map merely from its graph quartic.  THM-2992 shows locally that even a
fixed signed-edge parity packet selects at most an unordered two-sheet block;
it does not choose an affine owner.  Nor does `(6)` turn a six-vertex
tournament, Feuerbach configuration, or LRC packet into quartic root data.

## 6. Exact verification

Run

```bash
python 04-computation/quartic_resolvent_valuation_clutch_thm3046.py
python -O 04-computation/quartic_resolvent_valuation_clutch_thm3046.py
```

Both executions byte-match the stored `19`-line transcript
`05-knowledge/results/quartic_resolvent_valuation_clutch_thm3046.out`.  The
companion uses explicit exceptions and no truth-bearing Python assertions.
It verifies `(8)--(9)` symbolically; exhausts all `6^6=46,656` denominator-six
residue vectors; checks the `24` uniform quotient fibres, projector
integrality, all `24` sheet permutations on an edge basis, the tropical
Pluecker identity, explicit `5`-adic representatives of all `24` clutch
classes, and the three root-realizable controls `(15)`.

## 7. Boundary ledger

```text
PROVED HERE:             root/resolvent difference factorization;
                         exact valuation clutch realization;
                         half-discriminant ternary checksum;
                         S4/V4 equivariance;
                         tropical Pluecker boundary and all-class
                         5-adic realizability;
                         projector integrality criterion;
                         root-realizable independence hostiles.

NOT PROVED:              affine owner or source realization;
                         Keller/Jelonek restriction;
                         a canonical labelling after quotienting roots;
                         tournament, Feuerbach, or LRC transfer;
                         literal identity of the modular free factors with
                         binary and ternary combinatorial trees.
```
