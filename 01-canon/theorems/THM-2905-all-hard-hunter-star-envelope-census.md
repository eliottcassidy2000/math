---
id: THM-2905
title: "All-hard Hunter-star envelope census"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  The exact Hunter-star envelope closes
  2,964 of the 14,806 scalar-hard marked suffixes, including every prior
  direct pair-partition closure and 1,129 additional branches.  Six
  multi-hard bodies are new whole-root closures, taking the proved union to
  82 and the official seven-body residual to 3,350.
source: root-2026-07-29
depends_on:
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2899-all-root-ranked-suffix-scalar-census
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
  - THM-2903-one-hard-actual-h3-link-and-bad-triangle-closure
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2902-omission-six-ranked-h1-depth-one-closure
  - THM-2904-all-root-ranked-h1-hunter-closure
verification:
  - 04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.out
---

# THM-2905 -- all-hard Hunter-star envelope census

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

On the `14,806` scalar-hard marked suffixes of THM-2899, the exact
Hunter-star envelope from THM-2897 proves that `2,964` literal carriers
admit no five-cover.  This contains all `1,835` branches closed by
THM-2901's direct `q_5+2 beta_2<h` certificate and adds `1,129` branch
closures.

After recomposing all marked suffixes and subtracting the proved root union
through THM-2903, exactly six additional seven-body roots close.  Thus the
proved union has size `82`, and the official residual is

```text
3432-82=3350.                                             (1)
```

This is a strict finite-exact improvement, not a proof of LRC(14).

## 2. The Hunter-star envelope

Fix one scalar-hard marked suffix.  Let `C` be its literal carrier, `P` its
excluded gate prefix, and

```text
h=mu(C),
c(w)=mu(C intersect D_w),
q_1>=q_2>=...>=q_5
```

the five globally sealed singleton ranks over labels allowed outside `P`.
Let `beta_2` be THM-2901's attained global union cap for two distinct
allowed labels.

Suppose that an allowed five-set `Q` covers `C`.  Order its singleton
coverages as

```text
a_1>=a_2>=...>=a_5
```

and choose the `a_1` label as the centre of a star.  For each other label,
the centre--leaf intersection satisfies

```text
i_(1r)>=max(0,a_1+a_r-beta_2).                            (2)
```

Hunter's tree inequality on this star and `a_r<=min(a_1,q_r)` therefore
give

```text
U_C(Q)
 <=a_1+sum_(r=2)^5 min(a_r,beta_2-a_1)
 <=a_1+sum_(r=2)^5 min(a_1,q_r,beta_2-a_1).               (3)
```

The centre belongs to a distinct allowed pair, so `a_1<=beta_2`; also
`a_1<=q_1`.  Consequently

```text
G_5=max_(0<=a<=min(q_1,beta_2))
    [a+sum_(r=2)^5 min(a,q_r,beta_2-a)]                   (4)
```

is a valid upper bound on every five-label union.  The strict test

```text
G_5<h                                                    (5)
```

excludes a cover.

The objective in `(4)` is continuous and piecewise linear.  Its only
possible interior slope changes are

```text
a=q_r,          a=beta_2-q_r,          a=beta_2/2
                  for r=2,3,4,5.                         (6)
```

Thus the exact maximum occurs at `0`, at
`min(q_1,beta_2)`, or at a clipped value from `(6)`.  This is why the
postprocessing requires no new interval-union scan.

The same proof also gives the pointwise domination

```text
G_5<=min(q_1+q_2+q_3+q_4+q_5, q_5+2beta_2).              (7)
```

Hence every old direct closure must reappear in this census.

## 3. Exact all-hard census

The THM-2901 pair ledger records `q_1,q_2,q_3,q_5` but not `q_4`.  The
verifier therefore hash-pins and joins the complete THM-2899 top-five
ledger to the pair ledger by the exact tuple

```text
(body, rank, apex, excluded prefix).                      (8)
```

It also checks the stratum, carrier mass, component count, the four shared
singleton ranks, and the identities

```text
pair margin   =4h/7-beta_2,
direct margin =h-q_5-2beta_2.                            (9)
```

Every breakpoint in `(6)` is clipped to the full domain in `(4)`, including
both endpoints.  On every adjacent breakpoint cell, exact midpoint
linearity is checked independently of the maximization.  The locked result
is

```text
scalar-hard branches                                  14,806
Hunter-star strict closures                            2,964
old q_5+2beta_2 strict closures                        1,835
additional branch closures                             1,129
equality rows                                               0

stratum                    branches          G_5<h
low                           3,053             609
one                           7,853           1,571
both                          3,900             784.       (10)
```

None of THM-2901's `52` pair-cap exceptions satisfies `(5)`.  Therefore
the sharpened exact branch partition is

```text
Hunter-star closure                                   2,964
remaining finite H3-link obligations                 11,790
pair-cap exceptions                                      52
                                                     ------
                                                     14,806. (11)
```

Failure of `(5)` is only failure of a sufficient upper bound; it is not a
cover witness.

## 4. Whole-root recomposition

Exactly `16` bodies have every scalar-hard marked suffix closed by `(5)`.
Five are the old THM-2901 direct-terminal roots.  Five more belong to the
one-hard bank already closed by THM-2903; among those,
`(1,2,3,4,5,12,13)` is also a THM-2902 root.  The remaining six are new:

```text
(1,2,3,4,8,9,11),
(1,2,4,5,6,10,13),
(5,6,7,9,10,11,14),
(5,6,8,9,10,12,13),
(5,6,8,10,11,12,13),
(6,7,9,11,12,13,14).                                    (12)
```

The first five bodies in `(12)` have two scalar-hard suffixes; the last has
three.  Every other marked suffix closes by THM-2899, and THM-2892 supplies
the eight-body chamber.  Hence all six are whole-root closures.

The proved union before this theorem has size `76` by THM-2903.  Exact set
difference, rather than adding all `16` Hunter-terminal bodies, gives

```text
76+6=82,
```

which proves `(1)`.

## 5. Boundary and hostile controls

All theorem-bearing inequalities are strict.  There is no equality row.
The closest positive Hunter-star margin is

```text
1279/709320920
```

at

```text
E=(1,2,7,8,10,11,12), rank 7, apex 38.                  (13)
```

The closest failure is

```text
-233/73153080
```

at

```text
E=(1,4,5,8,9,11,13), rank 7, apex 28.                   (14)
```

The best of the `52` pair-cap exceptions still has negative margin

```text
-30403453/734894160
```

at

```text
E=(2,3,7,9,10,12,13), rank 2, apex 22.                  (15)
```

Thus the star correction substantially improves the H3-routed bank but
does not blur the separate exception family.

The faithful carrier here is a symmetric overlap-weighted graph on the five
candidate danger events.  The maximum singleton selects a proof-star centre,
but this is not an intrinsic orientation of the graph and does not define a
tournament.  The scalar envelope forgets which labels jointly attain the
rank and pair caps, all nonstar overlap, and every higher intersection; those
losses explain why `(5)` is sufficient rather than necessary.

## 6. Verification

The verifier uses exact rational arithmetic, hash-pins both inherited
ledgers, checks every join and identity in `(8)`--`(9)`, evaluates every
breakpoint in `(6)`, tests exact affine behavior between consecutive
breakpoints, and prints the branch and root consequence objects.  It uses
explicit guards rather than Python `assert`.  Ordinary and optimized
replays are byte-identical.

The complete semantic row digest is

```text
c1f60bdbc3669202dde60d8f44c9db887807b179cb42e02137d475c0f282d066.
```

Canonical artifacts:

```text
04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py
SHA-256 402498daeec59db0a9051c235ddcf0bdb4ee08250b0570e442bec6cf62782405

05-knowledge/results/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.out
SHA-256 c346cbce451b4d0104707b071c9874798e2cadc853102038b229be9ad8a6afe4
```

Reproduction:

```bash
python3 04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py
python3 -O 04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py
```

This theorem closes the six roots in `(12)`.  It does not close the
remaining `11,790` H3-link obligations, the `52` pair-cap exceptions, the
remaining `3,350` roots, or LRC(14). ∎
