---
id: THM-795
title: Hamming-one full-residue rigidity at every lift height
status: PROVED (scale-free safe-interval and sheet-deck descent) + VERIFIED (twelve exact core rows and seven residual atoms)
source: codex-2026-07-14-S10
depends_on: []
related: [THM-724, THM-769, THM-770, THM-775, HYP-6775, HYP-6800, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_one_full_residue_rigidity_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_one_full_residue_rigidity_codex_S10.out
---

# THM-795 — Hamming-one full-residue rigidity at every lift height

Put

```text
delta=1/13,                    [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

### A. Every proper one-coordinate lift is loose

For `r in [12]` and every integer `k>=1`, define

```text
W_k(r)=([12]\{r}) union {r+13k}.                         (1)
```

Then

```text
M(W_k(r))>1/13.                                           (2)
```

The set in (1) is a complete nonzero residue transversal modulo 13.  It has
clearance exactly `1/13` at every nonzero thirteenth, so (2) is strict
looseness, not merely a lower bound.

### B. Uniform rigidity around every shallow AP dilation

Let `c>=1` with `13` not dividing `c`.  Fix `r in [12]`, and replace `cr` in
`c[12]` by a positive integer `w` satisfying

```text
w = cr  (mod 13).
```

Put

```text
A=(c[12]\{cr}) union {w}.                                (3)
```

If `M(A)=1/13`, then

```text
w=cr.                                                     (4)
```

Equivalently, the only tight point in Hamming radius one around any shallow
arithmetic-progression dilation is the dilation itself.

This is scale-free: neither `c` nor the replacement lift height is bounded.
It upgrades THM-770's height-twelve finite classification in one complete
unbounded direction.  It does **not** prove uniform shallow rigidity for
packets that differ in two or more coordinates.

## 1. A metric safe-interval lemma

For a finite positive set `P`, write

```text
phi_P(t)=min_(p in P)||pt||,          B=max(P).
```

If

```text
phi_P(t_0)=mu>delta,
rho=(mu-delta)/B,                                         (5)
```

then

```text
I=(t_0-rho,t_0+rho) subset {t:phi_P(t)>delta}.             (6)
```

Indeed, every `t->||pt||` is `p`-Lipschitz, so

```text
phi_P(t)>=mu-B|t-t_0|>delta
```

on (6).

Suppose `P union {w}` were tight.  Then every point of `I` would lie in the
closed danger set

```text
D_w={t:||wt||<=delta}.
```

Its connected components are closed teeth of length `2/(13w)`, separated by
nonempty gaps.  The connected interval `I`, of length `2rho`, must lie in one
tooth.  Therefore

```text
w <= T(P,t_0):=1/(13rho).                                 (7)
```

Thus every `w>T` is eliminated at once.

## 2. Eleven-runner cores and exact thresholds

For each `r`, put `P_r=[12]\{r}`.  The following table gives a strict core
witness, its exact clearance, the Lipschitz radius (5), and threshold (7).

| `r` | `t_r` | `mu_r=phi_(P_r)(t_r)` | `B_r` | `rho_r` | `T_r` |
|---:|---:|---:|---:|---:|---:|
| 1 | `1/14` | `1/7` | 12 | `1/182` | `14` |
| 2 | `7/15` | `2/15` | 12 | `11/2340` | `180/11` |
| 3 | `5/16` | `1/8` | 12 | `5/1248` | `96/5` |
| 4 | `4/17` | `2/17` | 12 | `3/884` | `68/3` |
| 5 | `7/18` | `1/9` | 12 | `1/351` | `27` |
| 6 | `3/19` | `2/19` | 12 | `7/2964` | `228/7` |
| 7 | `1/7` | `1/7` | 12 | `1/182` | `14` |
| 8 | `1/8` | `1/8` | 12 | `5/1248` | `96/5` |
| 9 | `1/9` | `1/9` | 12 | `1/351` | `27` |
| 10 | `1/10` | `1/10` | 12 | `1/520` | `40` |
| 11 | `1/11` | `1/11` | 12 | `1/858` | `66` |
| 12 | `1/12` | `1/12` | 11 | `1/1716` | `132` |

For `r<=6`, the clearance numerators, in increasing order of the members of
`P_r`, are respectively

```text
r=1:  2,3,4,5,6,7,6,5,4,3,2                  over 14;
r=2:  7,6,2,5,3,4,4,3,5,2,6                  over 15;
r=3:  5,6,4,7,2,3,8,3,2,7,4                  over 16;
r=4:  4,8,5,3,7,6,2,2,6,7,3                  over 17;
r=5:  7,4,3,8,6,5,2,9,2,5,6                  over 18;
r=6:  3,6,9,7,4,2,5,8,8,5,2                  over 19.
```

This directly verifies the first six values of `mu_r`.  For `r>=7`, no
member of `P_r` is divisible by `r`, because `2r>12`; hence every clearance
at `1/r` is at least `1/r`, and the runner 1 attains equality.  This verifies
the last six rows.  The remaining columns are exact substitutions in (5)--(7).

## 3. Proof of Part A

First take `r<=6`.  The replacement is `w=r+13k`.  By (7), every case with
`w>T_r` is loose.  The cases satisfying `w<=T_r` are exactly

```text
(r,k)=(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(6,2).         (8)
```

They have the following direct full-packet witnesses.

| `(r,k)` | `w` | witness `t` | `phi_(W_k(r))(t)` |
|---:|---:|---:|---:|
| `(1,1)` | 14 | `1/16` | `1/8` |
| `(2,1)` | 15 | `9/19` | `2/19` |
| `(3,1)` | 16 | `6/17` | `2/17` |
| `(4,1)` | 17 | `5/19` | `2/19` |
| `(5,1)` | 18 | `4/19` | `2/19` |
| `(6,1)` | 19 | `4/23` | `2/23` |
| `(6,2)` | 32 | `7/44` | `1/11` |

Every final-column value is strictly greater than `1/13`.  This closes all
`r<=6`.

Now let `r>=7`.  If `r` does not divide `k`, evaluate the full packet at
`t=1/r`.  Every core runner has clearance at least `1/r`.  Also

```text
||(r+13k)/r||=||13k/r||>=1/r,
```

because `13` is a unit modulo `r` and the residue is nonzero.  Thus the whole
packet has clearance at least `1/r>1/13`.

If `r|k`, then `k>=r`, so

```text
w=r+13k>=14r>T_r,
```

where the last inequality is immediate from the last six rows of the table.
The safe-interval lemma applies.  This exhausts every `r` and `k`, proving
Part A. ∎

## 4. The missing-owner splice deck

We prove Part B.  Assume (3) is tight.  Choose `a in [12]` with

```text
ar=1  (mod 13),
```

and consider the `c` preimages of `a/13` under the circle covering
`t->ct`:

```text
t_l=(a+13l)/(13c),             0<=l<c.                    (9)
```

At `t_l`, the phase of the deleted core runner `cr` would be `+1/13`, while
the surviving runner `c(13-r)` has phase `-1/13`.  Every other runner of
`cP_r` has clearance at least `2/13`.  Thus `c(13-r)` is the unique core
binder at this splice.

On a compatible real lift, move a sufficiently small distance `h<0` from
`t_l`.  The signed phase of `c(13-r)` moves below `-1/13`, so its circle norm
becomes strictly larger than `1/13`; the other ten clearances remain larger
than `1/13` by continuity.  Hence there is a one-sided interval

```text
J_l=(t_l-epsilon_l,t_l) subset G_(cP_r).                  (10)
```

If the full packet (3) is tight, its replacement runner must satisfy

```text
||wt||<=1/13
```

throughout `J_l`.  The danger set is closed and `t_l` lies in the closure of
`J_l`; therefore

```text
||w t_l||<=1/13                  for every l.             (11)
```

The phases in (11) form a translate of a uniform circular grid.  Indeed, the
step from `l` to `l+1` is `w/c`, and the grid order is

```text
D=c/gcd(c,w).                                               (12)
```

All `D` grid points lie in the one closed phase arc
`[-1/13,1/13]`, of length `2/13`.  If `D>=2`, a `D`-grid has largest
complementary gap `1/D`, so its shortest containing arc has length

```text
1-1/D >= 1/2 > 2/13,
```

a contradiction.  Hence `D=1` and

```text
c|w.                                                       (13)
```

Write `w=cu`.  Because `c` is a unit modulo 13 and `w=cr mod 13`,

```text
u=r mod 13.
```

Positivity gives `u=r+13k` for some `k>=0`.  If `k=0`, then `w=cr`, as
claimed.  If `k>=1`, then

```text
A=c W_k(r).
```

Multiplication by `c` is onto on the circle, so

```text
M(cW)=M(W)                                                 (14)
```

for every finite `W`.  Part A and (14) give `M(A)>1/13`, contradicting
tightness.  Therefore only `k=0` is possible, which proves Part B. ∎

## 5. What this removes from the uniform n=12 frontier

In the shallow full-residue chart, encode a packet by its twelve labelled lift
coordinates.  THM-770 proved AP rigidity through lift height twelve by a
finite exact CSP.  The present theorem removes, at **every** height and after
**every** unit dilation, the complete Hamming-one star around the AP locus.

Consequently any still-hypothetical non-AP shallow tight packet must:

1. differ from every AP dilation in at least two labelled coordinates; and
2. evade both a forced-splice sheet-deck descent and every one-coordinate
   core-safe interval.

This is genuine uniform progress on the `n=12` sporadic branch, but it is not
uniform emptiness.  The exact remaining shallow problem starts at Hamming
radius two, where two replacement danger sets can divide the splice germs
between them and the one-arc grid argument no longer applies separately.

That observation suggests the next recursive object: a bipartite incidence
graph between missing-owner splice sheets and replacement teeth, decorated by
the deck orders `c/gcd(c,w_i)` and by one-sided core-safe germ orientation.
The proof target is a Hall-type obstruction or a forced common divisor after
two or more replacement colours share the sheets.

## 6. Tournament Analysis and assumption challenge

The decisive vertices in Part B are not runners.  They are the splice sheets
`l in Z/cZ` together with the predicate “the replacement is dangerous on the
one-sided core-safe germ ending here.”  This preserves the deck order in
(12); quotienting directly to runner vertices destroys the `c` simultaneous
obligations.

For telemetry, assign each sheet the scalar danger margin

```text
eta_l=1/13-||w t_l||.
```

Orient a pair toward the larger margin, break equalities by a marked cyclic
cut on `l`, and take the resulting sorted order as the tie Hamiltonian path.
This always produces a transitive tournament: score histogram
`0,1,...,c-1`, zero directed cycles, singleton SCCs, and one Hamiltonian path.
Moving the cut is the switch/gauge on tied sheets.

The fingerprint is intentionally uninformative: it cannot determine whether
**all** margins are nonnegative, which is the predicate used in (11).
Fixed-circle-section vertices preserve the forced splice locations but lose
owner identity; gap vertices preserve width but lose the gcd deck; tooth
vertices preserve the replacement cover but lose the core-safe germ.  The
minimal faithful carrier is therefore the sheet--danger incidence deck with a
metric, oriented germ sidecar.  The bare tournament is only an ordering
telemetry layer.

## Exact replay

The `Fraction` verifier checks every core numerator vector, every rational
`rho_r` and `T_r`, the exact residual list (8), and all seven direct witnesses.
It labels the two infinite `r>=7` branches and the grid-length descent used in
Part B.  Its finite sheet-margin tournament census is telemetry only; no
finite enumeration is used to justify either uniform theorem.
