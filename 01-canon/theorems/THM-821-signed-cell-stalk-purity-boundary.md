---
id: THM-821
title: Signed return-cell / deep-component stalk purity boundary
status: PROVED (pair-indexed uniform exact atomic factorization) + FINITE-EXACT at (13,5) (9,974 atomic stalks; explicit mixed fibres)
source: codex-2026-07-15-S10 deep signed-cell continuation
depends_on:
  - THM-803  # closed deep components and endpoint/cusp erosion selector
  - THM-817  # signed max-speed return cells and endpoint owners
related: [THM-774, THM-789, THM-807, HYP-6820]
verification:
  - 04-computation/lrc14_signed_cell_stalk_purity_codex_S10.py
  - 05-knowledge/results/lrc14_signed_cell_stalk_purity_codex_S10.out
---

# THM-821 — Signed return-cell / deep-component stalk purity boundary

Let `x>y` be positive odd exception speeds and put

```text
(a,b)=((x+y)/2,(x-y)/2),
Q_(x,y)(t)=||at||+||bt||,             theta=11/13.         (1)
```

For a ten-speed core `U`, retain THM-803's closed deep set and THM-817's
closed return cells:

```text
E_U={t:min_(u in U)||ut||>=1/11},
closure(R_U)=union_k R_k.                                  (2)
```

If `C` is one connected component of `E_U` and `R_k` is one signed return
cell, call

```text
S(C,k)=C+R_k                                               (3)
```

an **atomic erosion stalk**.  It succeeds when

```text
min_(t in S(C,k)) Q_(x,y)(t)>=theta.                       (4)
```

The sum in (3) is circular.  A proper arc is encoded exactly by its left
endpoint modulo one and its width.  Degenerate deep components are retained;
THM-803's deep set is closed.

## The result

### A. Exact pair-indexed numerical factorization

For each fixed pair in (1), the decorated sum state
`((x,y),S(C,k))` determines the atomic verdict (4).  Equivalently, the pair
index together with the exact input intervals `(C,R_k)` determines it.  On
each lifted representative of `S(C,k)`, the minimum is attained at an endpoint
or at a cusp

```text
j/(2a) or j/(2b),                    j in Z.               (5)
```

Thus, with the pair index retained, the exact selector argmin and the exact
signed margin

```text
mu_(x,y)(C,k)=min_(t in S(C,k))Q_(x,y)(t)-11/13           (6)
```

also determine the verdict.

This is a pair-indexed uniform atomic factorization, not a claim that the sum
arc or argmin determines a verdict after `(x,y)` is forgotten, or that a cell label,
component count, width, owner word, or unlabelled tournament determines the
LRC predicate, and not a proof that every admissible packet has a failing
atom.

The finite audit, explicit liar fibres, `U_n` family calculation, and
tournament telemetry below specialize to

```text
(x,y)=(13,5),                    (a,b)=(9,4).             (6a)
Q=Q_(13,5),                      mu=mu_(13,5).
```

### B. The finite-exact stalk audit

Reconstruct the `213` primitive rows in THM-817's deterministic cross-check
bank from seed `803807`.  Every row has one central return cell.  Their deep
components produce

```text
8,660 atomic stalks:       342 successes,       8,318 failures. (7)
```

Adjoin the complete disconnected row

```text
U_0=(1,2,3,4,7,500,503,504,505,506).                     (8)
```

It has `438` deep components and the three signed return cells `k=-1,0,1`,
so it contributes `1,314` atoms.  For each fixed label separately there are

```text
50 successes and 388 failures.                            (9)
```

The combined bank therefore has

```text
9,974 atoms:                  492 successes,
                            9,482 failures.               (10)
```

For a signature, a fibre is **mixed** when it contains both verdicts from
(4).  The exact table is:

| Forgetful signature | Fibres | Mixed | Largest mixed | Largest fibre |
|---|---:|---:|---:|---:|
| signed cell label `k` | 3 | 3 | 9,098 | 9,098 |
| exact return interval | 74 | 61 | 438 | 438 |
| return endpoint owners | 74 | 61 | 438 | 438 |
| full signed return cell | 74 | 61 | 438 | 438 |
| deep width | 2,718 | 65 | 668 | 668 |
| exact deep interval | 7,290 | 24 | 3 | 38 |
| deep endpoint owners | 3,519 | 234 | 144 | 144 |
| exact deep interval plus owners | 7,808 | 24 | 3 | 7 |
| cell label plus exact deep interval | 8,166 | 0 | 0 | 38 |
| exact return plus exact deep interval | 9,796 | 0 | 0 | 3 |
| full input stalk | 9,850 | 0 | 0 | 3 |
| sum width | 4,088 | 48 | 28 | 38 |
| exact circular sum arc | 9,796 | 0 | 0 | 3 |
| selector event shape | 3 | 2 | 4,804 | 4,804 |
| exact selector argmin | 9,086 | 0 | 0 | 60 |
| exact margin | 4,397 | 0 | 0 | 126 |

The zero-mixing result for `cell label + exact deep interval` is an empirical
fact about this bank.  It is **not** promoted to a general reconstruction
theorem: across cores, the label `k` alone need not determine the return-cell
endpoints.  In contrast, the zero mixing of the exact input intervals, exact
sum arc, exact argmin, and exact margin follows directly from Part A.

All `213` random rows fail the global containment obtained by requiring every
one of their atoms to succeed.  Hence row-level purity on this bank is
constant and vacuous; the nontrivial result is the atomic fibre audit (7)--
(10).

### C. Two exact liar fibres

The central-cell projection fails already within one random row.  Put

```text
U=(30,33,35,42,50,53,55,72,73,75),
R_0=[-2/10725,2/10725],             owners ((75),(75)).   (11)
```

The same complete signed cell has opposite verdicts on

```text
C_-=[1/330,2/165],
mu(C_-,0)=-17357/21450,

C_+=[309/803,106/275],
mu(C_+,0)=12049/156585.                                  (12)
```

Thus exact return position, width, sign, and return owners do not suffice.

Conversely, the exact deep interval alone fails as soon as returns disconnect.
Inside (8), take

```text
C=[49/132,2067/5566].                                    (13)
```

The negative and positive satellites give

```text
k=-1,  R_-=[-2/1001,-141/71500], owners ((7),(500)),
       mu(C,-1)=-557/12012;

k=+1,  R_+=[141/71500,2/1001],   owners ((500),(7)),
       mu(C,+1)=281/53625.                               (14)
```

The two sum arcs in (14) even have the same width
`2567/14610750`; their signed positions reverse the verdict.  This proves
that the apparent purity of the deep interval on the connected `213`-row
bank was accidental.

## Proof

For Part A, put `B=max(U)`.  Every connected component `C` is contained in one
`B`-safe band at threshold `1/11`, hence has width at most `9/(11B)`.
THM-817 puts every `R_k` inside one maximum-speed return tooth of width at most
`4/(143B)`.  Therefore

```text
width(C+R_k)<=9/(11B)+4/(143B)=121/(143B)<1.             (14a)
```

Thus the circular sum is a proper closed circular arc.  On every interval
avoiding the breakpoints (5), both triangular waves in (1) are affine.  Their
sum is affine, so its minimum on an arc occurs at an endpoint or at an internal
breakpoint.  This argument is uniform in the positive integers `a,b`, hence in
every odd pair `x>y`.  It is exactly THM-803's selector argument, now applied
before merging the raw `C+R_k` arcs.  Since

```text
E_U+closure(R_U)=union_(C,k)(C+R_k),                     (15)
```

global erosion containment holds precisely when every atom succeeds.

For Part B, the replay regenerates the `213` rows with Python's specified
Mersenne-Twister seed, then checks their bank digest before use.  It invokes
THM-817's max-speed cell constructor and independently compares its intervals
with THM-817's literal open-interval intersection.  It invokes THM-803's
closed deep-component constructor, but independently reconstructs every deep
endpoint owner by the equality `||ut||=1/11`.

For each pair `(C,R_k)`, the replay independently forms the lifted sum
interval, enumerates both endpoints and every point of (5) lying inside it,
computes the exact Fraction minimum, and records:

```text
(source row, C index, k,
 exact R endpoints, R owners,
 exact C endpoints, C owners,
 circular sum start and width,
 exact selector event and argmin,
 exact margin and verdict).                              (16)
```

Grouping all `9,974` records by each displayed signature gives the table.
Direct record lookup proves (11)--(14).  No floating point or sampled-circle
decision occurs.

## The `U_n` local family stalk

THM-817's family has

```text
B_n=506+360360n,
U_n=(1,2,3,4,7,B_n-6,B_n-3,B_n-2,B_n-1,B_n),
N_R(U_n)=3+1440n.                                        (17)
```

Since `360360=616*585`, every phase at

```text
s=479/616                                                (18)
```

is independent of `n`.  The exact local deep component containing `s` is

```text
C_n=[s-25/(616(B_n-3)),s],
left owner B_n-3,                  right owner B_n-2.     (19)
```

Indeed `B_n-2` has phase `10/11` at `s`, while the closest left constraint is
the `B_n-3` phase, at distance `25/(616(B_n-3))`.  In units of `1/616`, the
left phase slacks of the five high speeds with offsets `6,3,2,1,0` are

```text
436/(B_n-6), 25/(B_n-3), 504/(B_n-2),
367/(B_n-1), 230/B_n.                                    (19a)
```

The middle entry is strictly least for `B_n>=506`.  The least low-speed left
slack is `3/616`, from speed `4`, and is also larger.  This proves the left
endpoint uniformly; the zero right slack of `B_n-2` proves the right endpoint.

The central return cell and its owners are

```text
R_0=[-2/(143B_n),2/(143B_n)],       owners (B_n,B_n).    (20)
```

The nearest folded cusp to the right of `s` is `7/9`, at distance `1/5544`,
whereas the right extension in (20) is at most `1/36179`.  The left endpoint
is displaced by at most
`25/(616*503)+1/36179<17/616`, the distance to the nearest left cusp `3/4`.
Thus the whole sum stays in one affine cell.  There the two slopes are `-9`
and `+4`, so `Q` has slope `-5`.  Its minimum is therefore at the right
endpoint and equals

```text
Q(s+2/(143B_n))=69/616-10/(143B_n)<11/13.                (21)
```

Thus this explicit central stalk fails for every `n`.  The replay verifies
(19)--(21), and also the first and last positive satellite stalks, for the six
already emitted family rows `n=0,...,5`.  This is uniform local failure, not a
purity comparison: these selected family stalks all have the same verdict.

## Endpoint ownership, ancestry, and scope

For each chosen numerical predicate (4), exact endpoints of `C` and `R_k`
already determine the sum arc.  Endpoint owners can therefore be forgotten
*after* those intervals are given without changing this one verdict.  They
must nevertheless remain in the theorem-bearing carrier:

```text
(x,y) x [(signed cell, return endpoints and owners)
         x (deep component, deep endpoints and owners)]
    -> ((x,y), exact circular sum arc)
    -> selector event and mu_(x,y).                         (22)
```

Owners record which inequality moves an endpoint under speed replacement,
gcd descent, sheet transport, or recursive erosion.  The equality of a
single fixed arc after forgetting owners does not make that forgetting
functorial under those operations.

The finite audit proves none of the following:

- an all-size purity theorem for `cell label + deep interval`;
- a separator for globally successful versus globally failing rows;
- a mixed-fibre table or liar classification uniform in `(x,y)`;
- a bound on the number of deep or return components;
- a Čech, nerve, sheaf, or overlap-gluing theorem.

In particular, the replay emits raw pairwise Minkowski-sum stalks.  It does
not emit restriction maps, multiple-overlap data, or a cocycle equation, so a
Čech interpretation would add unproved structure.

## Tournament Analysis and challenged vertices

Tournament vertices are the sixteen forgetful signatures in the table, not
runners.  The first pair observable orders them lexicographically by

```text
(mixed-fibre count, largest mixed fibre, payload size),
```

and the switch puts payload size first.  Declaration order is the tie
Hamiltonian path.  Both completed tournaments are transitive, with score
histogram `(0,1,...,15)`, no directed triangle, sixteen singleton SCCs, one
Hamiltonian path, and `41` edge flips between gauges.

This telemetry expresses a real tradeoff—purity versus compression—but its
bare tournament is not the proof.  The exact mixed fibres (11)--(14) and the
factorization (15) are the theorem-bearing incidence data.  The challenged
assumption is that return cells, deep components, owners, widths, or event
types can be ranked separately and then used as a predicate.  Each such
projection has a mixed fibre already in the fixed `(13,5)` bank; their exact
pair-indexed signed incidence does not.

## Reproduction and certificate digests

```bash
python3 04-computation/lrc14_signed_cell_stalk_purity_codex_S10.py
```

The deterministic random-row bank has SHA-256

```text
303009db5bf61e2b5584f0664e740039aefda134e8ba80cf34de30cd897fcc71.
```

The full canonical certificate (all `9,974` decorated atoms, the fibre table,
both liar fibres, the six family stalk rows, and Tournament telemetry) has
SHA-256

```text
ebfe916e128a632c0e61bc5891df396f054b925747c56ca99ac7d3f252db4f8e.
```

The replay byte-matches its stored output, whose SHA-256 is

```text
bc7b53c51792cd790c2ebeb77b1cf8142a996a7f4aa26d261bc2d66b24275e16.
```
