---
id: THM-3056
title: "Tame quartic inertia clutch and binary-ternary index resonance"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a monic
  integral quartic over a strictly henselian DVR with tame cyclic splitting
  inertia, its six root-difference valuations split canonically into the
  first-order inertia-orbit skeleton plus a nonnegative projection-collision
  excess.  The total excess is exactly e times the monogenic order index,
  and its three perfect-matching sums give the exact correction to
  THM-3045/3046's F2[M]+F3 clutch.  Maximal orders have five explicit inertia
  fingerprints.  A sharp hostile shows that a transposition order of index
  one and a maximal three-cycle order have the same complete integer-normalized
  six-edge metric, matching sums, star sums, and clutch.  The missing scale is
  the base-normalized valuation (equivalently the tame ramification index),
  followed by a pointed fixed-sheet star and the full graph-order/affine-owner
  sidecars.  No quartic Keller branch is excluded.
source: codex-jc-quartic-valuation-2026-08-01
audit: >
  root/2026-08-01 independently checked the strict-henselian tame-cyclic
  reduction, all five faithful cycle types, the skeleton discriminant
  identity, the discriminant/index derivation of sum h=e*i, the fixed-sheet
  cross-resultant formula, both hostile discriminants and indices, and the
  owner-versus-fixed-sheet terminology boundary.  A follow-up audit derived
  every inertia-invariant matching cone and accepted the exact ambiguity
  classifier.  Ordinary, optimized, and stored transcripts byte-match after
  LF normalization.
depends_on:
  - THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch
  - THM-3046-quartic-resolvent-root-valuation-binary-ternary-clutch
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
  - THM-3042-subdirect-graph-order-common-quotient-and-singleton-owner-criterion
  - THM-3049-k4-matching-monomial-tropical-root-extraction-clutch
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/quartic_tame_inertia_clutch_index_resonance_thm3056.py
output: 05-knowledge/results/quartic_tame_inertia_clutch_index_resonance_thm3056.out
script_sha256: ba817e4a6bb9c5234a7774af1a5dc0c1cfbda0f8f04e8e64d7498be4e50f3d98
output_sha256: cb582d1e8d96a5fdb94857ef9bb2777cd02b9eb56eed3231713eedb8718c9a5b
hash_basis: LF-normalized bytes
---

# THM-3056 -- tame inertia has a scale-bearing matching-index clutch

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and statement

[THM-3045](THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md)
computes the exact `F2[M]+F3` quotient of the six-edge lattice of `K4`, and
[THM-3046](THM-3046-quartic-resolvent-root-valuation-binary-ternary-clutch.md)
realizes it through valuations of quartic and matching-resolvent root
differences.  THM-3046 deliberately fixes an integer normalization in the
splitting field and makes no claim that its residues survive a ramification
rescaling.  The theorem below identifies the missing scale, the exact order
tax hidden by that normalization, and the pointed coordinate required by
[THM-3038](THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary.md).

Let `R` be a strictly henselian DVR with fraction field `K`, normalized
valuation `v`, and residue characteristic different from `2` and `3`.  Let

```text
f(T) in R[T]                                               (1)
```

be monic, separable, and quartic.  Let `L/K` be its splitting field and
assume `L/K` is tame.  Because the residue field is separably closed, `L/K`
is totally ramified and cyclic.  Normalize its valuation `w` by

```text
w(L*)=Z,                   w restricted to K = e v,       (2)
```

where `e=[L:K]`.  Write `sigma` for a generator of `Gal(L/K)`, acting on the
labelled roots `z_0,...,z_3`.  Put

```text
x_ij=w(z_i-z_j),                                         (3)
```

and let `A=R[T]/(f)`, with normalization `Atilde` in its total fraction
algebra.  Define the monogenic order index

```text
i=length_R(Atilde/A).                                    (4)
```

For each edge `ij` of `K4`, define the first-order inertia skeleton

```text
c_ij=1  if i,j lie in the same nontrivial sigma-orbit,
c_ij=0  otherwise.                                      (5)
```

Then there are unique integers `h_ij>=0` such that

```text
x_ij=c_ij+h_ij,                                         (6)
```

and they obey the exact index identity

```text
sum_(i<j) h_ij=e i.                                     (7)
```

For the three perfect matchings

```text
m_1=01|23,            m_2=02|13,            m_3=03|12, (8)
```

put

```text
lambda_m=sum_(ij in m)x_ij,
lambda^sigma_m=sum_(ij in m)c_ij,
mu_m=sum_(ij in m)h_ij.                                 (9)
```

The vector `mu` is the **matching-index vector**.  It refines the scalar
order index exactly:

```text
lambda=lambda^sigma+mu,       mu_m>=0,
sum_m mu_m=e i.                                         (10)
```

Consequently THM-3045/3046's clutch satisfies

```text
kappa = kappa_sigma + mu                 in F2[M],
tau   = tau_sigma + e i                  in F3.          (11)
```

The binary correction needs the distribution `mu`, not only `i`.  The
ternary correction sees only its total.

If the monogenic order is maximal, `i=0`, then `h=mu=0`.  Up to permutation
of the three matching coordinates, the complete tame fingerprint is

| inertia type of `sigma` | `e` | normalized discriminant exponent `d` | `lambda^sigma` | `(kappa_sigma,tau_sigma)` |
|---|---:|---:|---:|---:|
| identity | 1 | 0 | `(0,0,0)` | `(000,0)` |
| transposition | 2 | 1 | `(1,0,0)` | `(100,1)` |
| double transposition | 2 | 2 | `(2,0,0)` | `(000,2)` |
| three-cycle | 3 | 2 | `(1,1,1)` | `(111,0)` |
| four-cycle | 4 | 3 | `(2,2,2)` | `(000,0)` |

Here `d=4-#orbits(sigma)`.  In arbitrary index, the intrinsic
base-normalized edge metric is

```text
xbar_ij=x_ij/e in Q,            sum_(i<j)xbar_ij=d/2+i. (12)
```

Thus the integer vector `x` without `e`, or equivalently without `xbar`, is
not a scale-bearing invariant of the base divisor.

There is also a pointed fixed-sheet gluing coordinate.  If `sigma` fixes `z_a`, then
`z_a in R`.  Writing

```text
f=(T-z_a)g_a,                    q_a=g_a(z_a)=f'(z_a),
t_a=sum_(j!=a)x_aj,                                     (13)
```

one has

```text
t_a=e v(q_a).                                           (14)
```

By THM-3038, the split monogenic singleton has gluing module `R/(q_a)`;
hence it splits precisely when `t_a=0`.  Formula `(14)` applies only after a
fixed sheet `a` has been identified.  It is not an affine-source owner test.

## 2. Proof of the skeleton and index formula

Every root of the monic integral polynomial is integral, so `x_ij>=0`.
Inertia acts trivially on the residue field.  Two distinct roots in the same
`sigma`-orbit therefore have the same reduction, so their difference has
positive integral valuation.  This proves `(6)` with `h_ij>=0`.

Let the orbit sizes of `sigma` be `ell_1,...,ell_r`.  The corresponding
factors of the quartic etale `K`-algebra are totally tamely ramified field
extensions of degrees `ell_j`.  Their maximal-order discriminant exponents
are `ell_j-1`, so

```text
d=sum_j(ell_j-1)=4-r.                                   (15)
```

The five possible cyclic cycle types on four letters give

```text
2 sum_(i<j)c_ij=e d.                                    (16)
```

This can be read directly from the table: the left half is respectively
`0,1,2,3,6`, while `(e,d)` is respectively
`(1,0),(2,1),(2,2),(3,2),(4,3)`.

The root discriminant and the order-index formula give

```text
2 sum_(i<j)x_ij
 =w(Disc(f))
 =e v(Disc(f))
 =e(d+2i).                                              (17)
```

Subtracting `(16)` from `(17)` proves `(7)`.  The three perfect matchings
partition the six edges, so summing `(9)` proves `(10)`.  Reducing `(10)`
modulo `2` and its coordinate sum modulo `3` proves `(11)`.  Nonnegativity
shows that `i=0` iff every `h_ij` vanishes, hence iff `lambda` equals its
fingerprint.  Division of `(17)` by `2e` proves `(12)`.

If `z_a` is fixed by the generator of the whole splitting group, then it is
in `K`, and integrality puts it in `R`.  Multiplicativity of `w` gives

```text
w(q_a)=sum_(j!=a)w(z_a-z_j)=t_a.                        (18)
```

Since `q_a in K`, equation `(2)` turns `(18)` into `(14)`.

## 3. Why the five maximal fingerprints are exact

For an inertia orbit of size greater than one, `c` is the complete graph on
that orbit.  Applying the three matching sums gives:

- one edge for a transposition, hence one matching value `1`;
- two opposite edges for a double transposition, hence one matching value
  `2`;
- a triangle for a three-cycle, and every perfect matching meets that
  triangle once, hence `(1,1,1)`;
- all six edges for a four-cycle, and every perfect matching contains two,
  hence `(2,2,2)`.

The identity has no skeleton edges.  Reduction modulo `2` and total modulo
`3` gives the last column of the table.  The table has a useful Keller-side
reading only in the maximal, correctly normalized local order.  Among the
fixed-point inertia types, a transposition has the one-hot fingerprint and a
three-cycle the all-one fingerprint.  Section 4 shows why dropping the
maximality/scale coordinate destroys this classifier completely.

There is nevertheless an exact all-index survivor.  Valuation invariance
forces `h` to be constant on the edge orbits of `sigma`.  Combining that
constraint with `(7)` gives the following necessary matching cones, up to
permuting matching coordinates; all displayed parameters are nonnegative
integers.

| inertia | complete matching vector `lambda` | order index |
|---|---|---:|
| identity | `(a,b,c)` | `a+b+c` |
| transposition | `(1+2a,b,b)` | `a+b` |
| double transposition | `(2+2a,2b,2c)` | `a+b+c` |
| three-cycle | `(1+i,1+i,1+i)` | `i` |
| four-cycle | `(2+2a,2+4b,2+2a)` | `a+b` |

For example, a transposition has two fixed edges and two two-edge orbits on
`E(K4)`.  The two repeated matching channels are therefore equal.  The
remaining matching excess is even because the total excess is `2i`.  A
three-cycle has one triangle orbit and one spoke orbit, so every matching
gets the same excess, necessarily `i`.  The other rows follow from the two
fixed/two paired edge orbits of a double transposition and the four-edge/two-
edge orbit split of a four-cycle.

These are necessary cones, not realizability theorems: residue units and the
ultrametric filtration can remove lattice points.  They still give a useful
nonmaximal classifier.  Among nontrivial inertia types with a fixed sheet:

```text
lambda non-diagonal                         => transposition;
lambda diagonal and even                    => three-cycle;
lambda diagonal and odd                     => ambiguous.            (18a)
```

At clutch resolution the same statement becomes

```text
transposition: kappa=(1,epsilon,epsilon), tau=1+2i mod3;
three-cycle:  kappa=(1+i mod2)(1,1,1),     tau=0.         (18b)
```

Thus a one-hot `kappa` forces the transposition lane, zero `kappa` forces the
three-cycle lane, and only `(kappa,tau)=(111,0)` remains ambiguous.  The
hostile in Section 4 realizes the smallest such ambiguity.

## 4. Sharp binary-ternary index resonance

Let `k` be algebraically closed of characteristic zero and put

```text
R=k[[t]],                       K=k((t)).                (19)
```

### 4.1 A transposition with index one

In `L_2=k((s))`, with `s^2=t`, take

```text
f_2(T)=(T^2-t)(T-t)(T-1),
(z_0,z_1,z_2,z_3)=(s,-s,t,1).                           (20)
```

The normalized splitting valuation has `w_2(s)=1`, so `e_2=2`.  Inertia
transposes `z_0,z_1` and fixes `z_2,z_3`.  In the edge order

```text
01,02,03,12,13,23,                                     (21)
```

the complete valuation vector is

```text
x^(2)=(1,1,0,1,0,0).                                   (22)
```

Moreover

```text
Disc(f_2)=4t^3(t-1)^6.                                  (23)
```

The maximal algebra is the product of one tame quadratic DVR and two copies
of `R`, with discriminant exponent `1`.  Thus `(23)` gives `i_2=1`.

### 4.2 A maximal three-cycle with the same entire metric

In `L_3=k((r))`, with `r^3=t` and a primitive cube root `zeta`, take

```text
f_3(T)=(T^3-t)(T-1),
(z_0,z_1,z_2,z_3)=(r,zeta r,zeta^2 r,1).                (24)
```

Now `e_3=3`, inertia cycles the first three roots, and the normalized edge
valuation vector is again

```text
x^(3)=(1,1,0,1,0,0)=x^(2).                             (25)
```

The factors in `(24)` are comaximal, the Eisenstein tame cubic order is
maximal, and

```text
Disc(f_3)=-27t^2(t-1)^2,              i_3=0.            (26)
```

Both packets therefore have

```text
lambda=(1,1,1),
star valuations=(2,2,2,0),
(kappa,tau)=(111,0).                                   (27)
```

This is stronger than a clutch collision: even the complete six-edge metric
and all four star sums coincide.  What differs is precisely the scale and
order data:

```text
xbar^(2)=x/2,             transposition,       index 1;
xbar^(3)=x/3,             three-cycle,         index 0. (28)
```

In `(10)--(11)`, the transposition skeleton `(1,0,0)` acquires matching-index
correction `(0,1,1)` of total `2=e_2 i_2`.  Its binary correction turns
`100` into `111`, while its ternary correction turns `1` into `0`.  This is
the exact binary-ternary **index resonance**.  The three-cycle packet has no
correction at all.

Therefore none of the following determines tame inertia or order index:

```text
integer-normalized matching valuations;
the F2[M]+F3 clutch;
all six integer-normalized edge valuations;
all four unpointed star valuations.                      (29)
```

The cheapest repair is to retain `e` or the base-normalized rational vector
`xbar`, then retain the inertia orbit labels.

## 5. One scalar index can move between fixed-sheet gluing and the complement

The transposition hostile `(20)` has two inertia-fixed sections, `a=t` and
`a=1`.  They distribute the same total order index in opposite ways.

For `a=t`, the complementary cubic and cross-resultant are

```text
g_t=(T^2-t)(T-1),
Disc(g_t)=4t(t-1)^2,               i_comp(t)=0,
q_t=g_t(t)=t(t-1)^2,               v(q_t)=1.             (30)
```

For `a=1`, they are

```text
g_1=(T^2-t)(T-t),
Disc(g_1)=4t^3(t-1)^2,             i_comp(1)=1,
q_1=g_1(1)=(t-1)^2,                v(q_1)=0.             (31)
```

Thus THM-3038's split index formula reads

```text
i_4=i_comp(a)+v(q_a)=1

a=t:                  (i_comp,v(q_a))=(0,1),
a=1:                  (i_comp,v(q_a))=(1,0).            (32)
```

The first fixed section is glued in the monogenic quartic order; the second
already splits.  Yet the total quartic index and the total standard-resolvent
index are both one in either pointed decomposition.  This gives a sharp
internal hostile to using a scalar discriminant/order tax as a fixed-sheet
splitting test.

The necessary pointed refinement is

```text
total index
  -> matching-index vector mu
  -> fixed-sheet star t_a/e=v(q_a)
  -> full graph-order discrepancy ideal
  -> affine-source regularity.                           (33)
```

The fourth step is THM-3042: extra graph generators may erase a nonunit
monogenic cross-resultant.  The fifth is THM-3038: even a finite-order
idempotent need not lie in the affine source.

## 6. Consequence for the planar quartic frontier

The closest proved Keller mechanism is THM-2633's fixed-branch constraint:
live `A4/S4` divisor inertia has fixed sheets, so the relevant nontrivial
cycle types are transpositions and three-cycles.  In a maximal tame
monogenic chart, Section 3 separates them exactly by `(1,0,0)` versus
`(1,1,1)`.  The resonance `(20)--(28)` proves that this separation is not
lawful for an arbitrary primitive inverse quartic.  A projection collision
of index one can make the transposition packet literally equal to the
maximal three-cycle packet.

Accordingly the next local computation for a hypothetical planar
degree-four Keller map must retain, at each relevant height-one divisor:

```text
1. the base valuation and ramification scale e;
2. the actual inertia orbit partition on the four normalization sheets;
3. the matching-index vector mu, not only the discriminant index i;
4. a labelled fixed sheet and its star valuation t_a/e;
5. the full graph-order common quotient/discrepancy ideal;
6. every original source-coordinate valuation on that trace.             (34)
```

THM-3049's tropical matching monomials supply item 3 only after items 1--2
are fixed; their root-extraction residues do not supply items 4--6.  The
counterfeit pair also explains why the primes `2` and `3` co-occur here
without yielding an exclusion: a binary ramification scale with one order
tax exactly imitates a ternary maximal orbit after integer normalization.

There is a direct warning for the open odd-Jelonek-exponent programme.  The
base discriminant exponent is

```text
v(Disc(f))=2 sum_m lambda_m/e=d+2i.                     (34a)
```

Thus its parity is the parity of the actual inertia discriminant exponent
`d`, not a function of the integer-normalized clutch.  In the resonance
hostile, the identical matching vector `(1,1,1)` and identical `tau=0` give
base exponents `3` in the transposition model and `2` in the three-cycle
model.  Consequently, once HYP-9027's proposed odd exponent is identified
with the actual integral local quartic order rather than a coefficient-chart
artifact, it would select the transposition lane in its stated 2-jet setting.
That law cannot be tested from THM-3046's clutch until the base ramification
scale and integral order have been restored.

This theorem decrements no planar-JC branch.  It does replace the vague
“missing valuation normalization” warning by the exact formulas `(7)`,
`(10)--(12)`, the pointed formula `(14)`, and the minimal hostile `(20)--(32)`.

## 7. Connection and boundary ledger

```text
source:       a monic integral tame quartic order at one base divisor;
target:       THM-3045/3046's three matching valuations and clutch;
map:          root-edge valuations -> inertia skeleton + excess h
              -> matching-index vector mu;
preserved:    exact discriminant/index tax and matching-channel location;
destroyed:    base scale if e is dropped, inertia labels, residue units,
              fixed-sheet choice, graph generators, affine incidence;
sidecars:     e or xbar, inertia partition, pointed star, full order, source
              coordinate valuations;
hostile:      (T^2-t)(T-t)(T-1) versus (T^3-t)(T-1).     (35)
```

```text
PROVED HERE:       tame skeleton/excess decomposition;
                   exact sum h=e*i and matching-index vector;
                   clutch correction and maximal inertia table;
                   pointed fixed-sheet cross-resultant valuation;
                   full-edge transposition/three-cycle resonance;
                   fixed-sheet-gluing/complement-index relocation.

NOT PROVED:        maximality of an unknown Keller graph order;
                   recovery of units or field-level square/cube roots;
                   full graph-order singleton membership;
                   affine-source ownership;
                   exclusion of A4, S4, G1, JC(2), or DC(2).              (36)
```

## 8. Exact companion

Run

```text
python3 04-computation/quartic_tame_inertia_clutch_index_resonance_thm3056.py
python3 -O 04-computation/quartic_tame_inertia_clutch_index_resonance_thm3056.py
```

Both modes LF-byte-match the stored transcript.  The companion enumerates
the five cyclic inertia types; checks `(16)` and every fingerprint; exhausts
all `3^6` nonnegative correction vectors for each type as an independent
lattice control of `(10)--(11)`; exhausts bounded assignments on every
inertia edge orbit to verify the necessary matching cones; and verifies both
hostile discriminants, edge metrics, matching and star vectors, order indices,
base scales, pointed cross-resultants, and the two complementary-cubic index
decompositions.
Every truth-bearing check uses an explicit runtime exception rather than a
Python assertion.

QED.
