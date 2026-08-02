---
id: THM-3095
title: "Marked S3 affine-lift half-face and oriented tetrahedral co-occurrence frame"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Every exact-order affine lift of a marked generating C2*C3 -> S3
  pair is uniquely a quartic owner, a matching direction, a cyclic direction
  orientation, and the split/full H1 bit.  The first mixed word obeys the
  exact formula w^2=eta d.  In the full branch the 24 S4 markings are the
  4*3*2 oriented tetrahedral co-occurrence frames, and the matching quotient
  forgets exactly the owner.  The carrier agrees with THM-3067 after
  orientation is forgotten, but the modular actions do not.  A labelled
  quartic frame selects a signed pair sum; the quotient retains only its
  squared resolvent channel.  No global affine owner, signed-edge descent,
  graph-quartic realization, Keller cofactor, or exclusion follows.
source: root-modular-resolvent-cooccurrence-2026-08-02
depends_on:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
  - THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss
  - THM-3083-exceptional-binary-point-ternary-direction-s4-tomography-clutch
  - THM-3092-modular-mixed-word-fingerprint-and-septimal-counterfeit-separation
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2775-modular-s4-to-weyl-d3-generator-frame-and-affine-parity-blindness
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
script: 04-computation/modular_marked_s3_affine_lift_half_face_thm3095.py
output: 05-knowledge/results/modular_marked_s3_affine_lift_half_face_thm3095.out
script_sha256: c5009bf45ae40a388d3e24aeb0bca6254afbda39c289c3dbee2cc4e0c594ff37
output_sha256: b151fcf86735600af483281613366f3c76aa07470c159f8a1f67184942e3fc1c
hash_basis: LF-normalized bytes
---

# THM-3095 -- the marked resolvent forgets the owner, not the half-face

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

## 1. The complete marked lift theorem

Let `X` be a four-element set and let `D` be its three perfect matchings.
The action on matchings is the standard quotient

```text
pi:S4 -> Sym(D)=S3,                    ker(pi)=V4.        (1)
```

Identify each `d in D` with the nonidentity element of `V4` whose two
transpositions are the two blocks of `d`.  Thus `d` is simultaneously a
matching and a translation direction on the affine four-point torsor.

Fix a marked generating pair

```text
(sbar,cbar) in S3 x S3,
ord(sbar)=2,          ord(cbar)=3,       <sbar,cbar>=S3. (2)
```

Put `wbar=sbar cbar`.  The product `wbar` is an involution of `D`, so it
fixes a unique matching.  Denote that matching/direction by

```text
d=d(sbar,cbar).                                          (3)
```

An **exact-order affine lift** of `(2)` is a pair `(s,c) in S4 x S4` with

```text
pi(s)=sbar,       pi(c)=cbar,       ord(s)=2,       ord(c)=3. (4)
```

The three-cycle `c` fixes a unique sheet `o in X`.  It also supplies a
cyclic orientation `epsilon` of the three directions.  Intrinsically,
`epsilon` belongs to the two-element torsor of cyclic orders on `D`; it is
not an absolute sign until one orientation has been chosen.  Equivalently,
the bijection

```text
X\{o} -> D,       x |-> the unique matching direction sending o to x (5)
```

transports the cyclic order of `c` on `X\{o}` to the cyclic order of `cbar`
on `D`.

There are exactly eight lifts `(4)`.  Four generate a point-stabilizer
`S3`, and four generate all of `S4`.  Put

```text
eta=0 in the split S3 branch,          eta=1 in the full S4 branch. (6)
```

Below, `eta d` denotes the identity of `V4` for `eta=0` and the nonzero
translation `d` for `eta=1`.

Then the complete invariant is

```text
(s,c) |-> (o,d,epsilon,eta),

{exact-order lifts of all marked pairs (2)}
    ~= X x D x Or(D) x F2,                               (7)
```

and the first mixed word `w=sc` satisfies the exact half-face law

```text
w^2 = 1                  if eta=0,
w^2 = d in V4\{1}        if eta=1.                       (8)
```

In particular, restricting `(7)` to `eta=1` gives an equivariant bijection

```text
Epi(C2*C3,S4)
   ~= {(o,d,epsilon):o in X, d in D, epsilon in Or(D)},
|Epi(C2*C3,S4)|=4*3*2=24.                                (9)
```

This refines THM-2595 and THM-3092.  The nonzero co-occurrence class does
not merely say that a half-face exists: the marked quotient already says
**which** of the three matching directions it is.

## 2. Explicit inverse and proof

Fix `(o,d,epsilon)`.  There is a unique three-cycle `c` fixing `o` and
inducing `epsilon` on `D`.  Put

```text
b=d(o),             a=c^(-1)(b),             e=c(b).     (10)
```

Then

```text
c=(a b e),                    d=(o b)(a e).               (11)
```

The two binary completions are

```text
s_0=(b e)=(b,c(b)),                    split branch,
s_1=(o a)=(o,c^(-1)(b)),               full branch.       (12)
```

Indeed,

```text
w_0=s_0c=(a e),                 w_0^2=1,
w_1=s_1c=(o a b e),             w_1^2=(o b)(a e)=d.      (13)
```

Both `w_0` and `w_1` induce the unique transposition of `D` fixing `d`.
Consequently `s_0` and `s_1` have the same image `sbar` in `(2)`.  The first
pair fixes `o` and generates the symmetric group on `X\{o}`.  The second
contains the cross transposition `(o a)` and its three `c`-conjugates, so it
generates `S4`.

Conversely, an order-three lift of `cbar` is a three-cycle and hence has one
fixed owner `o`.  For fixed `o`, the fibre of a quotient transposition has
exactly two order-two lifts: the inner and cross transpositions in `(12)`.
Thus `(12)` exhausts `(4)` and proves the bijection `(7)`.  It also identifies
`eta` with the `H^1(C2*C3,V4)=F2` class of THM-2595.

The count over all quotient markings is therefore

```text
6 marked S3 pairs * 4 owners * 2 branches = 48 lifts,
24 split + 24 full.                                      (14)
```

Every full direction occurs eight times, and every owner occurs six times.
Simultaneous conjugation by `S4` is free and transitive on the `24` full
triples.  Hence `(9)` is an `S4`-equivariant identification, not merely a
cardinality coincidence.

## 3. The exact forgetting diagram

For a fixed full marking, quotienting by `V4` acts on `(9)` by

```text
(o,d,epsilon) |-> (d,epsilon).                           (15)
```

The right side has `3*2=6` elements and is exactly the set of marked
epimorphisms `C2*C3 -> S3`.  Every fibre of `(15)` contains the four owners
once each.  Thus, **within the full branch**, the cubic matching quotient
forgets exactly the quartic owner and retains both the mixed-word direction
and its ternary orientation.

For all exact-order lifts, not only the full ones, the corresponding map is

```text
(o,d,epsilon,eta) |-> (d,epsilon),                       (16)
```

so the unpointed quotient also forgets the split/full cohomology bit.  The
same marked cubic pair therefore has four split and four full lifts.

Forgetting only `epsilon` in `(9)` leaves the twelve pairs

```text
(o,d) in X x D.                                          (17)
```

These are literally the directed tetrahedral flags `V x (V\{0})` of
THM-3067.  This is a carrier equivalence, **not an action equivalence**.
THM-3067 uses its binary modular move as a nonzero `V4` translation and
generates `A4`; here the binary generator is the transposition in `(12)`
and the full branch generates `S4`.  Importing the THM-3067 move law through
the common twelve-set would be a type error.

The projection ledger is

| retained data | size | exact meaning |
|---|---:|---|
| `(o,d,epsilon,eta)` | 48 | all exact-order affine lifts of full marked `S3` shadows |
| `(o,d,epsilon)` | 24 | full `S4` markings |
| `(o,d)` | 12 | oriented-edge/point-direction carrier, ternary orientation forgotten |
| `(d,epsilon)` | 6 | marked `S3` resolvent quotient, owner forgotten |
| `(o,epsilon)` | 8 | the eight three-cycles, binary completion forgotten |
| `d` | 3 | one matching/resolvent channel only |

## 4. The signed quartic survivor

Let a depressed separable quartic have roots labelled by `X`:

```text
f(T)=T^4+pT^2+qT+r,               sum_(x in X) alpha_x=0. (18)
```

A flag `(o,d)` selects the signed pair sum

```text
s_(o,d)=alpha_o+alpha_(d(o)),              u_d=s_(o,d)^2. (19)
```

The other block of the matching has sum `-s_(o,d)`.  Hence changing `o`
within one block leaves `s` fixed, while moving `o` to the other block
changes `s` to `-s`.  The square `u_d` is independent of the owner and is
the root of the matching resolvent

```text
S(U)=U^3+2pU^2+(p^2-4r)U-q^2                             (20)
```

labelled by `d`.  It is the translated standard resolvent `R_f(p+U)` of
THM-2598.  Thus the full co-occurrence frame chooses a **signed**
matching channel at a labelled fibre, whereas its marked `S3` quotient
chooses only the squared resolvent channel.

When `s=s_(o,d)` is nonzero, the exact quadratic blocks are

```text
a=(p+u_d+q/s)/2,                 b=(p+u_d-q/s)/2,

f(T)=(T^2-sT+a)(T^2+sT+b).                               (21)
```

Indeed the product differs from `f` in its constant coefficient by
`S(u_d)/(4u_d)`.  Formula `(21)` is the signed-edge factorization of
THM-2992, now reached from a marked modular frame.  It does not make `s`
an element of the base field or of a local integral order.

The smallest exact control labels the roots as

```text
(alpha_0,alpha_1,alpha_2,alpha_3)=(0,1,2,-3).             (22)
```

Take

```text
s=(01),             c=(123),             w=sc=(0123).
```

Then

```text
w^2=(02)(13),          s_(0,w^2)=2,          u_(w^2)=4.  (23)
```

Conjugating the marking by the kernel translation `(01)(23)` while keeping
the labelled quartic fibre fixed moves the ternary fixed owner from `0` to
`1`.  It preserves the quotient and `w^2`, but the selected sum becomes
`alpha_1+alpha_3=-2`.  The resolvent root `4` survives.  Simultaneously
relabeling the quartic roots would instead preserve the signed sum by
equivariance; that is not the forgetting map being tested here.  This is the
sharp owner/sign loss inside one full quotient fibre.

## 5. Cusp, grade-three, and Keller boundary

Suppose the marked `S3` pair is realized by genuine geometric meridians of
an `A2` cusp in a normalized cubic resolvent.  The quotient word `wbar` is
a transposition, and `d` is its unique fixed cubic sheet label.  In the full
quartic lift, `(8)` identifies the same label with the square of the
order-four quartic meridian.  In the split lift, the same cubic label remains
but the square is trivial.  This is exactly the two-branch cusp lift of
THM-3037, now resolved into owner, direction, orientation, and cohomology
coordinates.

The conclusion stops at that marked finite cover.  In particular:

1. The combinatorial owner `o=Fix(c)` is a sheet at a chosen labelled fibre.
   It is not a monodromy-invariant section, a retained point of an affine
   Zariski-main open, or the separated fixed factor required by THM-3064.
2. The signed value in `(19)` need not descend to the strict-henselian unit
   `s in R*` required by THM-2992.  Its Kummer phase, valuation, and block
   discriminants are additional data.
3. Neither `(7)` nor `(20)` supplies the primitive-element cofactor or its
   inverse-different unit.  THM-3064 shows that this cofactor is
   load-bearing for a Keller decoder.
4. The depressed cubic and cube-minus-square discriminant scaffold are
   universal quartic identities by THM-2455/2598, not grade-three Keller
   equations.  THM-2681 moreover proves that the specific THM-2473/1310
   cubic normalization on its principal chart cannot be the classical
   resolvent-root field of a dimension-three `S4` Keller map.

The typed connection is therefore

| field | value |
|---|---|
| source | an exact-order marked affine lift of a full `S3` matching action |
| map | `(s,c) -> (Fix(c), Fix(pi(sc)), cyclic direction order, H1 bit)` |
| target | an oriented tetrahedral flag and, in the full branch, a nonzero `V4` half-face |
| preserved | marked free-factor orders, ternary orientation, selected matching/resolvent label, split/full bit before quotient |
| destroyed by `S4 -> S3` | quartic owner and split/full bit; hence signed pair-sum phase |
| still absent physically | base-ring realization, signed-edge descent, inertia/index, affine owner, cofactor/inverse different |
| cheapest decisive test | the `48`-lift census plus one signed-root control `(22)--(23)` |

No degree-four Keller map, `G1`, `JC(2)`, `DC(2)`, SFC, GMC, or LRC
exclusion follows.

## 6. Sharp same-quotient hostile

The branch loss in `(16)` already occurs in the smallest labelled example.
With zero-based labels, take

```text
c=(123).
```

The two pairs

```text
split:  s_0=(23),       s_0c=(13),       (s_0c)^2=1,
full:   s_1=(01),       s_1c=(0123),     (s_1c)^2=(02)(13) (24)
```

have the same marked matching quotient.  The first fixes sheet `0` and
generates `S3`; the second generates `S4`.  Thus even the complete marked
cubic quotient does not determine whether its fixed direction is realized
as a quartic half-face.  The extra bit is exactly `eta`.

## 7. Exact companion and current status

Run

```bash
python 04-computation/modular_marked_s3_affine_lift_half_face_thm3095.py --check-stored
python -O 04-computation/modular_marked_s3_affine_lift_half_face_thm3095.py --check-stored
```

The self-contained companion uses explicit `require` gates and no Python
`assert` statements.  It checks:

1. the `S4 -> S3` matching quotient and exact `V4` kernel;
2. all six marked quotient pairs and all `48` exact-order affine lifts;
3. the `24+24` split/full census and `w^2=eta d` for every lift;
4. the explicit inverse `(10)--(12)` and singleton frame fibres;
5. all `24` full triples, twelve flags, six quotient fibres, and every
   simultaneous-`S4` equivariance identity;
6. four independent rational depressed-quartic root controls, all matching
   resolvent roots, every nonzero signed factorization `(21)`, and the
   `2 -> -2` origin hostile; and
7. the exact same-quotient split/full pair `(24)`.

Normal and optimized executions byte-match the stored sixteen-line
transcript after LF normalization.  Promotion remains blocked only on an
independent hostile audit of this candidate and its scope.
