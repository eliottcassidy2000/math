---
id: THM-3539
title: "Fixed Keller newest-prime decomposition bound and quadratic LCA packet floor"
status: >
  PROVED + VERIFIED-EXACT for the wreath centralizer/normalizer, its complete
  predecessor-block and unordered-pair LCA orbit classification, and the
  resulting conditional packet formula for THM-3538.  The proved Keller data
  give only I<=D<=C_W(I), not D=C_W(I).  If the decomposition-group image has
  the full leaf-stabilizer point-and-pair orbits, the exponentially many
  carrier factors reduce to exactly n^2 valuation packets.  THM-3540 and
  THM-3542 prove this gate at n=2 and n=3; n>=4 and all-level prescribed-
  coordinate index zero remain OPEN.
source: codex/newest-prime-decomposition-sidecar/2026-08-16
depends_on:
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
  - THM-3538-fixed-keller-newest-prime-prescribed-coordinate-index-criterion
related:
  - THM-3531-fixed-keller-intrinsic-all-level-discriminant-square-class
  - THM-3537-fixed-keller-level-two-old-L-inertia-and-x-index
  - THM-3540-fixed-keller-depth-two-newest-prime-residue-saturation
  - THM-3542-fixed-keller-depth-three-newest-prime-residue-point-pair-saturation
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_newest_prime_decomposition_lca_packet_audit_20260816.py
output: 05-knowledge/results/keller_newest_prime_decomposition_lca_packet_audit_20260816.out
script_sha256: 8f9816389a701198a5809abf1f952c4b30393e908a9fed787c1d7a36beb356b3
output_sha256: 7155daf78e5f787c05f23bbb3429d0b3e330bdd54efeb4a9c3bc693fe2ef6294
semantic_sha256: 96b5ba9f5c0c66e29c5b40d20f77b11242de706a7c8b98c66647cddc283efbea
hash_basis: LF-normalized bytes
---

# THM-3539 -- the newest-prime inertia has a quadratic LCA shadow, conditionally

**PROVED GROUP THEORY + CONDITIONAL KELLER REDUCTION + VERIFIED-EXACT.**

Retain the fixed sporadic Keller inverse tower of THM-3535.  At depth `n`,
its Galois closure has full monodromy

```text
W_n=S_3 wr ... wr S_3                                  (1)
```

on the ternary words of length `n`.  Put

```text
r=n-1,                 B={0,1,2}^r,                  m=3^r.       (2)
```

The set `B` indexes the predecessor blocks in THM-3538.  Mark the unique
block `b_*` meeting `L=0` at the newest prime.  THM-3533 and THM-3535 give a
tame inertia generator `t` which is one transposition in the bottom `S_3`
over `b_*` and is the identity on every other leaf.

This theorem determines everything that follows from that group datum alone,
and no more.

## 1. Exact centralizer, normalizer, and fixed-sibling form

Write the last wreath layer as

```text
W_n=S_3^B semidirect W_r.                              (3)
```

If `tau` is the transposition supporting `t`, then

```text
C_(W_n)(t)
 = (C_(S_3)(tau) times product_(b!=b_*) S_(3,b))
     semidirect Stab_(W_r)(b_*).                       (4)
```

Since `C_(S_3)(tau)=<tau>`, this also gives

```text
N_(W_n)(<t>)=C_(W_n)(t).                               (5)
```

Indeed, a normalizer of the order-two subgroup must fix its unique nonidentity
element.  Let `ell_*` be the third leaf in the block `b_*`, namely the unique
sibling in that bottom block not moved by `t`.  The local stabilizer of
`ell_*` in `S_3` is exactly `<tau>`, so `(4)` gives the sharper identity

```text
C_(W_n)(t)=Stab_(W_n)(ell_*).                          (6)
```

This is the precise meaning of “centralizer/leaf stabilizer”: it is the
stabilizer of the third sibling in the ramified bottom block, not the
stabilizer of an arbitrary one of the many leaves fixed by `t`.

For completeness,

```text
|W_n|=6^((3^n-1)/2),
|C_(W_n)(t)|=|W_n|/3^n
             =2^((3^n-1)/2) 3^((3^n-1)/2-n).          (7)
```

Thus the conjugacy class has size `3^n`, exactly the number of bottom
transpositions: there are `3^(n-1)` bottom blocks and three transpositions in
each.

To prove `(4)`, write `g=(s_b)_b h` in `(3)`.  The conjugate `gtg^(-1)` is
supported over `h(b_*)`.  Equality with `t` therefore forces `h(b_*)=b_*`;
inside that block it forces `s_(b_*)` to centralize `tau`.  Every other
`s_b` is free.  These conditions are also sufficient.

## 2. What the Keller theorems prove about the decomposition group

Let `M_n/K_0` be the Galois closure and let `mathfrak P` be a prime over the
newest divisor

```text
p=(P_(n-1))                                             (8)
```

whose inertia is the supported transposition above.  Write `D=D_mathfrak P`
and `I=I_mathfrak P`.  The proved local anatomy gives exactly

```text
I=<t>.                                                  (9)
```

Inertia is normal in the decomposition group.  Conjugation by `D` therefore
acts on `I`; but `Aut(C_2)=1`.  Equations `(5)`--`(6)` now give the unconditional
bound

```text
<t> = I  <=  D  <=  N_(W_n)(I)
                 =  C_(W_n)(t)
                 =  Stab_(W_n)(ell_*).                (10)
```

This is the strongest decomposition-group conclusion presently implied by
THM-3530/3533/3535.  In particular, they do **not** prove

```text
D=C_(W_n)(t).                                          (11)
```

The missing coordinate is the residue extension in the exact sequence

```text
1 -> I -> D -> Gal(kappa(mathfrak P)/kappa(p)) -> 1.   (12)
```

THM-3538 passes to a strict henselization in order to split the predecessor
cover.  That operation separably closes the residue field and is deliberately
blind to `(12)`.  Full global monodromy supplies the ambient `W_n`, while the
supported inertia supplies the left end of `(10)`; neither determines the
residue group in the middle.

Both endpoints `D=I` and `D=C_(W_n)(t)` satisfy the currently proved abstract
group constraints.  This is a group-data hostile, not a claim that both occur
for the actual Keller divisor.  Deciding the actual `D` requires a residue-
field splitting theorem or an equivalent decomposition-cover computation.

## 3. Centralizer orbits on predecessor blocks

The bottom base group in `(4)` acts trivially on `B`.  Hence the induced
centralizer action on predecessor blocks is exactly

```text
H_r=Stab_(W_r)(b_*).                                   (13)
```

Identify `b_*` with `0^r`, and write `lambda(u,v)` for the depth of the least
common ancestor of two ternary words.  There are precisely `r+1=n` orbits on
`B`:

```text
{b_*},                                                  size 1;
S_d={b:lambda(b_*,b)=d},     0<=d<r,                   size s_d,
s_d=2*3^(r-d-1).                                      (14)
```

The stabilizer is transitive on each `S_d`: it may swap the two off-ray
children at depth `d` and then acts as the full ternary automorphism group in
the chosen off-ray subtree.  Different `d` cannot mix because rooted-tree
automorphisms preserve LCA depth.

## 4. Centralizer orbits on unordered block pairs

The `H_r`-orbits on `binom(B,2)` are the isomorphism types of the minimal
rooted subtree on the marked leaf `b_*` and an unordered pair `{u,v}`.  They
are completely classified as follows.

1. **Marked pairs.**  `{b_*,u}` with `u in S_d`, for `0<=d<r`.  There are
   `r` types, of size

   ```text
   s_d=2*3^(r-d-1).                                    (15)
   ```

2. **Unequal shells.**  Put
   `a=lambda(b_*,u)<b=lambda(b_*,v)`.  Then necessarily
   `lambda(u,v)=a`.  For every `0<=a<b<r` there is one orbit, of size

   ```text
   s_a s_b=4*3^(2r-a-b-2).                             (16)
   ```

3. **Equal shells.**  Put
   `lambda(b_*,u)=lambda(b_*,v)=d` and
   `c=lambda(u,v)`.  Every `0<=d<=c<r` gives one orbit.  Its size is

   ```text
   3^(2(r-d-1)),                    c=d;
   2*3^(2r-c-d-2),                  d<c.               (17)
   ```

For `c=d`, the two leaves enter the two different off-ray children at depth
`d`.  For `c>d`, there are two choices of their common off-ray child, then a
common prefix to depth `c`, an unordered choice of two of three children, and
two free descendant suffixes.  This proves `(17)`.  The same local
transitivity used in `(14)` proves that each displayed type is one orbit.

The pair-orbit census is therefore

```text
r + binom(r,2) + r(r+1)/2 = r(r+1)=n(n-1).             (18)
```

The sizes in `(15)`--`(17)` sum to `binom(3^r,2)`, providing an independent
partition check.

## 5. Conditional `n^2` reduction of the THM-3538 carrier

For one of THM-3538's integral observations `theta=y,z`, or `theta=u=1/x`
under its reciprocal chart gate, abbreviate

```text
A_b=h_theta(q_b),
R_(b,c)=Res(f_(theta,b),f_(theta,c)).                  (19)
```

Its exact local index is

```text
i_(theta,n)
 =sum_(b in B) v(A_b)+sum_({b,c} in binom(B,2)) v(R_(b,c)).       (20)
```

Let `bar D` be the image of `D` in `W_r`.  Equation `(10)` gives only

```text
bar D <= H_r.                                          (21)
```

The factors in `(19)` are equivariant under `bar D`, and `D` preserves the
chosen valuation.  Their valuations are therefore constant on `bar D`-
orbits.

Use the following explicit **decomposition saturation gate**:

```text
G_(2-orb): bar D has the same orbits as H_r on both
           B and binom(B,2).                           (22)
```

The stronger statement `bar D=H_r` is sufficient, and `(11)` is stronger
still.  Neither is currently proved.  Under `(22)`, choose one representative
factor for every type in `(14)`--`(17)`.  Then `(20)` becomes

```text
i_(theta,n)=v(A_*)+sum_d s_d v(A_d)
 +sum_d s_d v(R^* _d)
 +sum_(a<b) s_a s_b v(R^< _(a,b))
 +sum_d 3^(2(r-d-1)) v(R^= _(d,d))
 +sum_(d<c) 2*3^(2r-c-d-2) v(R^= _(d,c)).              (23)
```

There are exactly

```text
(r+1) + r(r+1) = (r+1)^2 = n^2                       (24)
```

representatives: `n` internal factors and `n(n-1)` cross-block resultants.
Since every valuation in `(20)` is nonnegative, the carrier is a unit iff all
`n^2` representatives in `(23)` are units.  Thus the exponential factor list
reduces to `O(n^2)` orbit types under `(22)`.

This is an orbit-count reduction, not a polynomial-time algorithm: producing
one representative resultant at depth `n` may itself remain expensive.

## 6. The unconditional quadratic floor and the lost information

Because `bar D` is only a subgroup of `H_r`, its orbits refine the LCA orbits.
Consequently the number `N_D` of valuation packets always satisfies

```text
n^2 <= N_D <= m+binom(m,2)
              =3^(n-1)(3^(n-1)+1)/2.                  (25)
```

The lower endpoint is attained by the full point-and-pair orbit gate `(22)`.
The upper endpoint is the trivial predecessor action, which is compatible
with the proved subgroup `D=I` at the group-data level because bottom inertia
acts trivially on `B`.  Hence the current theorems do not yield an
unconditional `O(n^2)` test.  Equation `(25)` is a **quadratic orbit floor**:
even the maximal allowed symmetry cannot collapse distinct rooted LCA types.

The packet map preserves exactly what `(23)` uses:

- the marked ancestry block;
- each block's LCA depth from it;
- for a pair, whether it contains the marked block and the three-leaf rooted
  LCA type; and
- under `(22)`, valuation and unit/nonunit status with the orbit multiplicity.

It destroys the individual ternary addresses, the choice of off-ray child,
the labelled strict-henselian section, the residue-field extension in `(12)`,
the actual unit residue and sign of a factor, and all distinctions inside one
LCA packet.  It also supplies no arithmetic reason that any representative in
`(23)` is a unit.

## 7. Exact boundary

THM-3539 proves no all-level coordinate-index equality.  In particular, it
does not prove:

1. `D=C_(W_n)(t)` or even `bar D=H_r` for the actual newest divisor;
2. the saturation gate `(22)`;
3. unitness of the `n^2` representatives at any new level `n>=5`;
4. a recurrence for their values or valuations;
5. transport to old primes, where THM-3537's ramified predecessor hostile
   invalidates the split newest-prime model; or
6. an arbitrary-Keller-map classification, `JC(2)`, `DC(2)`, or LRC(14).

THM-3540 subsequently proves the saturation gate `(22)` at `n=2`.  THM-3542
proves it at `n=3` by a good rational specialization whose point factors have
degrees `1,2,6` and whose injective pair-sum factors have degrees
`1,2,6,6,9,12`; the carrier has nine packets.  Neither theorem identifies the
full decomposition group, and `(22)` remains open for every `n>=4`.

The theorem instead isolates the exact missing sidecar: compute the residue
decomposition action `(12)`.  If it is point-and-pair saturated, THM-3538's
all-level problem becomes a quadratic family of arithmetic noncollision
tests; without that gate, exponential splitting can remain.

## 8. Exact companion

Reproduce with

```text
python -B 04-computation/keller_newest_prime_decomposition_lca_packet_audit_20260816.py
python -B -O 04-computation/keller_newest_prime_decomposition_lca_packet_audit_20260816.py
```

The companion checks `(7)` through depth seven, generates the full marked-
leaf stabilizer and recovers every point/pair orbit and every size in
`(14)`--`(17)` through predecessor depth five, checks the symbolic census
through `n=12`, and records the trivial-image hostile.  Ordinary and optimized
transcripts match the stored output exactly after LF normalization.

**QED.**
