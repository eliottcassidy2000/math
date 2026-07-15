---
id: THM-816
title: Uniform looseness of the quartic order-three three-sheet interface
status: PROVED (uniform comb-discrepancy descent) + FINITE-EXACT (7,909-state recursive component closure)
source: codex-2026-07-15-S10/S11 quartic-s3 continuation
depends_on: [THM-810]
related: [THM-769, THM-772, THM-803, THM-807, HYP-6820]
verification:
  - 04-computation/lrc13_quartic_s3_recursive_comb_closure_codex_S11.py
  - 05-knowledge/results/lrc13_quartic_s3_recursive_comb_closure_codex_S11.out
  - 04-computation/lrc13_quartic_s3_global_endpoint_crosscheck_codex_S11.py
  - 05-knowledge/results/lrc13_quartic_s3_global_endpoint_crosscheck_codex_S11.out
---

# THM-816 — Uniform closure of the quartic order-three interface

Put

```text
F=F_13^*,              H=<5>={1,5,8,12}.
```

For a multiplicative coset `R=aH`, choose positive integers `u_r`, one for
each `r in R`, such that

```text
u_r=3r (mod 13),       3 does not divide u_r,
u_r=u_(-r) (mod 3).                                      (1)
```

Then the twelve-speed packet

```text
A_R=3([12] minus R) union {u_r:r in R}                   (2)
```

is uniformly loose:

```text
M(A_R)>1/13.                                             (3)
```

Consequently the exceptional order-three alternative of THM-810 is empty at
every lift height.  Any hypothetical tight four-coordinate perturbation of a
full-residue AP must take THM-810's common-scale alternative; the remaining
Hamming-four problem is the scale-one chart, not a quartic three-sheet
residual.

The congruence conditions make all twelve speeds in (2) distinct.  There are
three possible cosets and, by THM-810, four possible mod-three patterns.  No
height cutoff is assumed in the statement.

## 1. The exact owner-sheet carrier is a danger-comb cover

Write `delta=1/13` and, for a positive integer `w`, put

```text
D_w={t in R/Z: ||wt||<=delta}.                           (4)
```

For a coset `R`, let

```text
E_R={t:min_(q notin R) ||3qt||>delta}.                   (5)
```

This is a finite union of open rational intervals.  The full packet has no
strictly good time exactly when

```text
E_R subset union_(r in R) D_(u_r).                       (6)
```

Thus (6), not the scalar oriented-deck inequality, is the exact residual
predicate.  Under `tau=3t`, a point of `E_R` records one of the three lifted
owner-sheet obligations over the loose eight-speed quotient core.  Incidence
with a tooth of `D_(u_r)` records which exception blocks that obligation.
This is precisely the owner-sheet carrier isolated in THM-810, represented on
the original time circle so that all overlaps remain visible.

THM-810's lift-invariant `q=39` clock is consistent with this formulation,
not an obstruction to (3). At each of its eight points two opposite core
speeds have distance exactly `1/13`, so the point lies on the boundary of
`E_R`, not in its strict interior. The clock explains why perturbing a common
equality witness cannot prove looseness. The recursion below instead proves
that some open core-safe component survives all four lifted exception combs.

## 2. A sharp interval-comb discrepancy lemma

Let `E` be a union of `K` disjoint open intervals of total length `L`.  Put

```text
p=2/13.
```

For every positive integer `w`,

```text
|E intersect D_w| <= pL + p(1-p)K/w
                  = 2L/13 + 22K/(169w).                 (7)
```

### Proof

On one period of the normalized danger comb, the indicator has mean `p`.
A periodic primitive of `1_D-p` rises over the centered danger interval and
falls over its complement; the difference between its maximum and minimum is
exactly

```text
p(1-p).
```

After the substitution `x=wt`, the discrepancy on one arbitrary interval is
therefore at most `p(1-p)/w`.  Sum this inequality over the `K` components of
`E`.  Open/closed endpoints have measure zero, so the half-open convention is
irrelevant to (7).  This proves the lemma.  ∎

Suppose `m<=4` danger combs of speeds at least `v` cover `E`.  Sum (7) over
the `m` combs.  Since `mp<1`,

```text
(1-mp)L <= K p(1-p) m/v,

v <= B_m(E):=22mK/[13(13-2m)L].                         (8)
```

This converts every unbounded next speed into an exact bound determined only
by the current residual component union.

## 3. Recursive completeness of the finite closure

Order the four off-sheet speeds increasingly:

```text
v_1<v_2<v_3<v_4.
```

They are distinct because they have distinct nonzero residues modulo 13.
Starting with `E_0=E_R`, define

```text
E_j=E_(j-1) minus D_(v_j).                               (9)
```

If the full row were tight, then `E_j` would be covered by the `4-j`
remaining combs.  Applying (8) at every stage gives the necessary recursive
bound

```text
v_(j+1)<=floor(B_(4-j)(E_j)).                            (10)
```

Conversely, every positive solution of (1) lies in one unique progression

```text
u_r=b_r+39k,            k>=0,                            (11)
```

where `1<=b_r<=38` is fixed by the label and its negative-pair parity.  Hence
an increasing recursion which tries every unused label and every member of
(11) between the preceding speed and (10) contains the increasing ordering
of every hypothetical covering row.  A branch with no candidate can be
discarded by (8); a branch reaching four choices is tight exactly when its
residual interval union is empty.  This proves that the recursion is an
exhaustive reduction of the infinite problem, rather than a sampled height
search.

The exact root data are:

| missing coset `R` | quotient labels `[12] minus R` | `#components(E_R)` | `|E_R|` | first-speed bound |
|---|---|---:|---:|---:|
| `(1,5,8,12)` | `(2,3,4,6,7,9,10,11)` | 36 | `8893/45045` | 246 |
| `(2,3,10,11)` | `(1,4,5,6,7,8,9,12)` | 36 | `772/4095` | 258 |
| `(4,6,7,9)` | `(1,2,3,5,8,10,11,12)` | 42 | `443/2145` | 275 |

Every component is obtained by exact intersection with the safe bands

```text
((13k+1)/(13w), (13(k+1)-1)/(13w)),       0<=k<w.        (12)
```

The complete recursive audit over the three cosets and four parity patterns
is:

| chosen speeds | exact prefix states | minimum next bound | maximum next bound | no-candidate leaves |
|---:|---:|---:|---:|---:|
| 0 | 12 | 246 | 275 | 0 |
| 1 | 324 | 130 | 398 | 0 |
| 2 | 3,046 | 71 | 395 | 1,062 |
| 3 | 4,493 | 32 | 265 | 4,459 |
| 4 | 34 | — | — | — |

There are `7,909` states in total.  No prefix covers its residual safe set.
All 34 four-speed terminal states retain a nonempty open residual; even the
smallest terminal residual has exact measure

```text
41681/720720.                                            (13)
```

Thus (6) never holds, proving (3).  All comparisons and interval endpoints in
the replay are rational.  There is no floating point, time sampling, or
unproved numerical cutoff.  The search certificate digest is

```text
c603151d5d6606f0492bb6c484dc976fb57befa258d0ee823b9f24e16d7f0a15.
```

As an independent finite cross-check, the largest bound anywhere in the
recursive tree is 398.  A second implementation cuts the circle at every
threshold endpoint from every permitted progression value through 398, turns
the resulting 15,273--16,817 open cells into exact integer masks, and tests
the full Cartesian product rather than following the recursive tree.  It
checks 132,510 rows and again finds zero covers.  Its certificate digest is

```text
a02c30e784969c5865606455a38843f44ff49fb6a1b96faee2548a64204f8b83.
```

## 4. Tournament Analysis and assumption challenge

For telemetry on each of the 34 terminal rows, take the four exception combs
as vertices.  In the **raw gauge**, orient `u_i -> u_j` when

```text
|E_R intersect D_(u_i)| >= |E_R intersect D_(u_j)|,     (14)
```

using increasing speed as the fixed tie Hamiltonian path.  For the switched
**conditional gauge**, first remove the danger combs of the other two
vertices and compare the two marginal erosions of the remaining set.  This
switch tests whether pairwise priority survives the higher-order overlap
context.

Every raw and every conditional tournament is transitive: score histogram
`{0:1,1:1,2:1,3:1}`, zero directed triangles, four singleton SCCs, and one
Hamiltonian path.  Nevertheless the gauges have respectively 0, 1, 2, or 3
edge flips on

```text
15, 14, 4, 1
```

of the terminal rows.  The usual tournament fingerprints report no change
while the conditional incidence changes substantially.  They also omit the
absolute residual measure which decides (6).

This challenges the assumption that runners, residue labels, gaps, or fixed
safe components should be the final vertices.  Static components are not
closed under the operation: a new tooth can split one component into several.
The predicate-preserving object is the **dynamic residual interval union**
together with its bipartite incidence to individual labelled danger teeth.
Equivalently, its vertices are current owner-sheet component obligations and
comb teeth, and its arrows are actual threshold-cover incidences.  This
carrier preserves exactly whether a good time remains.  It deliberately
forgets the value of `M(A_R)` above `1/13`, which is not needed for uniform
looseness.  The tournament is useful overlap telemetry only after this richer
carrier has been retained.

## Exact replay

Run

```bash
python3 04-computation/lrc13_quartic_s3_recursive_comb_closure_codex_S11.py
python3 04-computation/lrc13_quartic_s3_global_endpoint_crosscheck_codex_S11.py
```

The canonical output is stored at the path in the frontmatter.
