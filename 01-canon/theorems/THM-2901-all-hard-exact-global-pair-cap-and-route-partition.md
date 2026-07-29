---
id: THM-2901
title: "All-hard exact global pair cap and route partition"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  On all 14,806 scalar-hard marked
  THM-2899 suffixes the global allowed two-comb cap is attained and exactly
  sealed.  Its direct partition certificate closes 1,835 branches and five
  whole roots; the remaining branch routes split as 12,919 finite H3
  pair-flag children and 52 exact pair-cap exceptions.
source: codex/lrc-j6-parity-seed-audit-2026-07-29
depends_on:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2898-unique-max-gate-five-parity-matching-closure
  - THM-2899-all-root-ranked-suffix-scalar-census
related:
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
  - THM-2900-flag-conditioned-rank-selective-partition-closure
verification:
  - 04-computation/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.ledger.out
---

# THM-2901 -- all-hard exact global pair cap and route partition

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Exact pair cap

Fix any of the `14,806` scalar-hard marked suffixes in THM-2899.  Retain
its literal carrier `C`, mass `h`, component count `r`, ordered excluded
prefix `P`, and allowed external labels.  Put

```text
c(w)=|C intersect D_w|,
q_j=the j-th largest allowed c(w),
gamma=(99/70)r/7,
beta_2=sup_{distinct allowed x,y}|C intersect (D_x union D_y)|.
```

THM-2899 gives `q_1<3h/7`, while THM-735(ii) gives

```text
c(w)<h/7+gamma/w.                                        (1)
```

Define

```text
W_2=floor(gamma/(3h/7-q_1))+1.                           (2)
```

If one endpoint of a pair is at least `W_2`, subadditivity and `(1)` put
its union strictly below `4h/7`.  Let `H_0` be the exact largest pair
union with both endpoints below `W_2`.  On every one of the `14,806`
branches the exact computation finds

```text
eta_2=H_0-q_1-h/7>0.                                     (3)
```

Set

```text
R_2=ceil(gamma/eta_2),       X_2=max(W_2,R_2).             (4)
```

For `w>=X_2`, every pair containing `w` has union

```text
<q_1+h/7+gamma/w <= q_1+h/7+eta_2=H_0.                   (5)
```

Thus the maximum over pairs with both endpoints below `X_2` is the
attained global value `beta_2`.  The strict inequality in `(1)` makes the
possible equality at the ceiling in `(4)` harmless.  This is an exact
self-seal, not a finite-window proposal.

Every finite head is evaluated by singleton-sum branch-and-bound followed
by exact literal two-comb unions.  The winning pair is independently
checked by subtracting both danger combs from `C`.  The exact ranges are

```text
55 <= W_2 <= 696,             139 <= R_2 <= 3182,
1 <= paid exact unions <= 155.
```

The computation pays `212,869` exact unions rather than enumerating all
`1,967,632,698` finite-head pairs.

## 2. Exact branch partition

THM-2897 gives three consequences of `beta_2`.

First,

```text
q_5+2beta_2<h                                             (6)
```

excludes a five-cover by partitioning its five labels into a least-covered
singleton and two disjoint pairs.  Inequality `(6)` holds on exactly
`1,835` branches.

Second,

```text
beta_2<4h/7                                               (7)
```

holds on exactly `14,754` branches.  THM-2893 then forces at least three
members of any hypothetical five-cover into the finite core

```text
H_3={w:c(w)>=(h-beta_2)/3}.                               (8)
```

After taking the direct closures first, the exact route partition is

```text
direct certificate (6)                                  1,835
finite H3 pair-flag / link-constrained three-cover      12,919
pair-cap exception                                          52
                                                         ------
                                                         14,806.  (9)
```

For the middle line one may enumerate pairs `L subset H_3`, since
`2<3`; no heaviness condition is imposed on `L`.  Put

```text
C_L=C minus D_L,       h_L=|C_L|,
J_L={w:c_(C_L)(w)>=h_L-beta_2}.                           (9a)
```

The heavy-link recursion in THM-2893 shows that all three labels of a
possible child cover lie in `J_L`: for each such label `w`, the parent
triple `L union {w}` is heavy and hence

```text
c_(C_L)(w)=U_C(L union {w})-U_C(L)>=h_L-beta_2.           (9b)
```

The child also retains `P union L` as forbidden labels.  Thus the `12,919`
rows are finite link-constrained child obligations, not closures.
An inherited parent-carrier cap may close some flags, but this is not
automatic; the same-cap deadlock applies to selected flags of full arity
`ell=s`, not to these `ell=2<s=3` pairs.  Exact child or inherited-cap
work must preserve the link and forbidden-prefix sidecars.

The `52` failures of `(7)` occur on `51` rank-one branches and one
rank-two branch, split by THM-2896 stratum as

```text
low 6,       one 32,       both 14.
```

Their closest strict failure is

```text
-417/3923920
```

at `E=(2,3,8,10,12,13,14)`, rank `1`, apex `22`; the most negative is

```text
-114553/5885880
```

at `E=(2,4,6,10,12,13,14)`, rank `1`, apex `22`.  Failure of the
sufficient cap is not evidence for a cover.

## 3. Five new whole-root closures

Combining `(6)` on the scalar-hard suffixes with THM-2899 on all earlier
suffixes closes exactly these five bodies:

```text
(1,2,3,4,6,11,12),
(1,2,3,5,6,10,13),
(1,2,4,5,6,12,13),
(1,3,4,5,6,7,14),
(5,7,8,10,11,13,14).                                    (10)
```

The last body has three hard suffixes; the other four have one each.
Every hard suffix of each body satisfies `(6)`.  The ordered hitting gate
therefore excludes every external six-cover.  THM-2892, already used in
THM-2899, supplies the eight-body chamber.

The five bodies in `(10)` are disjoint from the five scalar-terminal
bodies of THM-2899, the four THM-2895 bodies, and THM-2898's unique
maximum-gate body.  Hence the current proved union contains `15` roots,
and the official seven-body residual is

```text
3432-15=3417 roots.                                      (11)
```

Among the `65` roots with exactly one scalar-hard suffix, all `65` satisfy
`(7)` and four satisfy `(6)`.  This isolates a small, especially clean
test bank of `61` actual `H_3` pair-child problems.

## 4. Boundary and scope

The exact cap `beta_2` forgets compatibility and overlap between two
separately maximizing pairs.  The direct certificate `(6)` is sufficient,
not necessary.  The cutoff cardinalities recorded for `(8)` are upper
bounds on the actual core size; binomial counts from those cutoffs are
not actual flag counts.

This theorem proves the five root closures in `(10)`, the exact global
pair-cap atlas, and the route partition `(9)`.  It does not close the
`12,919` child obligations or the `52` pair-cap exceptions, close the
remaining `3417` roots, or prove LRC(14).

## 5. Verification

The verifier hash-pins the THM-2899 hard ledger and both exact interval
engines.  It reconstructs every literal carrier, verifies every winning
pair by independent residual subtraction, and uses explicit guards rather
than Python `assert`.  Locked ordinary and optimized replays are
byte-identical.  The semantic row digest is

```text
f3f63ac1953c0e2292488d91f59f831e0f7273b9c1eaeda32f74932d974e4ee4.
```

Canonical artifacts:

```text
04-computation/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.py
SHA-256 7ba8244d8fc78ebc0d9381e05d69ca53c849d6008ff9cfb43f0efcbb4b394f81

05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.ledger.out
SHA-256 5dea0eaa45dd52fbf1bef7cfcc328899a4789bc277b6e1e8ac2f4bdf192b85e4

05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.out
SHA-256 98b8ba171be1d38980e7271ef82e2bc1bde536afcf9864fa39138dbfbc93a3eb
```
