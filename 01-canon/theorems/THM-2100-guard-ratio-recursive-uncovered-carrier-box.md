---
id: THM-2100
title: "The bounded guard-ratio terminal has a uniform recursive carrier box"
status: >
  PROVED. THM-2095 bounds the odd guard and one commensurate terminal speed
  by 14,396,844,825. Successively remove bounded danger combs from the guard
  complement. THM-2086's sharp overlap spectrum leaves a positive carrier
  floor at every one of the six stages, while one-dimensional interval
  discrepancy forces another remaining speed to be bounded by the carrier's
  component count. The exact recursion gives max(Q)<=534912686676376 and,
  through THM-2077, max(S)<=22822941298192042. This uniformly boxes the whole
  guard-ratio branch without coefficient templates; it does not discharge
  the box or prove LRC(14).
source: codex-2026-07-22-LRC-guard-ratio-recursive-carrier
depends_on:
  - THM-2077
  - THM-2086
  - THM-2095
related:
  - THM-1135
  - THM-2080
  - THM-2097
script: 04-computation/lrc14_guard_ratio_recursive_box_codex_20260722.py
output: 05-knowledge/results/lrc14_guard_ratio_recursive_box_codex_20260722.out
script_sha256: e03c8af9f8b0076dd4af403b15373ecf5bb5068e9827e0ed8a6ccd3a03593aee
output_sha256: be495283fb54c67d25c4464bab04e526568dd07e6c1bf49bd6b5739e9ac2385d
hash_basis: repository blobs with LF line endings
---

# THM-2100 -- a uniform recursive box for the guard-ratio branch

Retain a depth-four rank-seven terminal obstruction

```text
G_Q subset E_h,
Q={q_1,...,q_7} hereditarily primitive,
h odd,
```

with seven distinct terminal speeds, in THM-2081's notation. Suppose
THM-2087 Branch I holds. THM-2095 supplies a marked `q_1` and proves

```text
h,q_1<=H:=57(3^4*5^2*11*17*23*29)=14,396,844,825.     (1)
```

Then

```text
max Q<=534,912,686,676,376,                            (2)
max S<=22,822,941,298,192,042.                         (3)
```

The second bound is for the original thirteen-speed row.

## 1. Nested uncovered guard carriers

Write

```text
D_q={t:||qt||<1/14},
E_h={t:||ht||<1/7},
C_h=E_h^c.
```

For `K subset Q`, `|K|=k`, put

```text
A_K=C_h minus union_(q in K)D_q.                       (4)
```

The `7-k` unselected danger combs cover `A_K` pointwise. Otherwise a point
missed by them and by the selected combs would lie in `G_Q intersect E_h^c`,
contrary to guard containment.

Let `I(q,h)=measure(D_q intersect E_h)`. Since `measure(C_h)=5/7` and
`measure(D_q)=1/7`, the union bound gives

```text
measure(A_K)>=(5-k)/7+sum_(q in K)I(q,h).              (5)
```

THM-2086's exact odd-guard spectrum says `I(q,h)>=1/35`, except for

```text
q=6h,       I=1/42,
q=h/11,     I=2/77,                                   (6)
```

where the second exception exists only when integral. Distinctness permits
at most one occurrence of each exceptional ratio. Therefore

```text
delta_1=25/42,
delta_k=(5-k)/7+1/42+2/77+(k-2)/35,  2<=k<=6,          (7)
```

and `measure(A_K)>=delta_k`, where

```text
(delta_1,...,delta_6)
=(25/42,221/462,841/2310,577/2310,313/2310,7/330).    (8)
```

The final positive gap `delta_6=7/330` is load-bearing.

## 2. Exact component count

The guard boundary has `2h` points and the boundary of `D_q` has `2q`
points. Hence the boundary of `A_K` has at most

```text
2(h+sum_(q in K)q)                                    (9)
```

points. The set is closed, has positive measure by (8), and is not the whole
circle because it lies in `C_h`. After discarding isolated null points it is
therefore a union of at most

```text
R_K:=h+sum_(q in K)q                                  (10)
```

circular intervals with disjoint interiors. Coincident or retained closed
endpoints only merge components. Equivalently, cut inside the nonempty open
guard `E_h`; then no positive-length component wraps around the cut.

## 3. Harmonic selection lemma

For every circular interval `I` and positive integer `q`,

```text
measure(I intersect D_q)<=measure(I)/7+6/(49q).        (11)
```

Indeed, the primitive of the one-period mean-zero function
`1_(D_1)-1/7` has oscillation exactly `6/49`; scaling by `q` divides the
endpoint discrepancy by `q`. The bound is independent of the interval's
position and length.

Summing (11) over the components in (10) yields

```text
measure(A_K intersect D_q)
 <=measure(A_K)/7+6R_K/(49q).                          (12)
```

Let `m=7-k`. Coverage by the unselected combs and (12) imply

```text
measure(A_K)
 <=m measure(A_K)/7+(6R_K/49)sum_(q in Q minus K)1/q,
```

so

```text
sum_(q in Q minus K)1/q
 >=7k measure(A_K)/(6R_K)
 >=7k delta_k/(6R_K).                                 (13)
```

One of the `m` reciprocal terms is at least their average. Thus some
unselected speed satisfies

```text
q<=c_k R_K,       c_k=6(7-k)/(7k delta_k).             (14)
```

This is an unordered selection lemma: choose any witness supplied by (14),
adjoin it to `K`, and repeat. No size ordering or generic-position assumption
enters.

## 4. The exact recursion

Start with THM-2095's marked speed and put `U_1=H`. If the sum of the `k`
selected speeds is at most `U_k`, then `R_K<=H+U_k=:R_k`. Define

```text
b_(k+1)=floor(c_k R_k),       U_(k+1)=U_k+b_(k+1).     (15)
```

Exact arithmetic gives

| `k` | `delta_k` | `c_k` | `R_k` | `b_(k+1)` |
|---:|---:|---:|---:|---:|
| 1 | 25/42 | 216/25 | 28,793,689,650 | 248,777,478,576 |
| 2 | 221/462 | 990/221 | 277,571,168,226 | 1,243,418,355,401 |
| 3 | 841/2310 | 2640/841 | 1,520,989,523,627 | 4,774,568,778,091 |
| 4 | 577/2310 | 1485/577 | 6,295,558,301,718 | 16,202,606,721,059 |
| 5 | 313/2310 | 792/313 | 22,498,165,022,777 | 56,928,264,210,988 |
| 6 | 7/330 | 330/49 | 79,426,429,233,765 | 534,912,686,676,376 |

The selected-sum endpoint is `U_7=614,324,719,065,316`. Every earlier
selected speed is below the final entry, and the sixth step supplies the last
unselected speed. This proves (2).

## 5. Lift and exact scope

THM-2077 equation (13), for depth `r` and terminal maximum `B`, says

```text
max(S)<=2^(r+1)B max(1,(12-r)/(r+2)).                  (16)
```

Terminal size seven means `r=4`; hence

```text
max(S)<=32B*(8/6)=128B/3.                              (17)
```

Substituting the integer `B` from (2) and taking the floor gives (3).

The theorem turns THM-2097's templatewise finiteness into a uniform,
coefficient-template-free integer box for THM-2095's complete guard-ratio
branch. It does not enumerate or discharge those rows.

The useful objects are the nested carriers `A_K`, not a tournament on
runners. A selection-order tournament is transitive and forgets both
quantities doing the work: carrier mass and boundary complexity. The faithful
state is

```text
(measure(A_K), number of interval components of A_K),
```

with the exceptional ratio labels `6h` and `h/11` retained.

## Exact referee

The companion reconstructs all rational arrangement endpoints independently,
checks 1,280 overlap-spectrum rows, 316 carrier mass/component rows, 3,160
finite-union discrepancy rows, and every recurrence integer. Normal and
`python -O` transcripts byte-match and end in `PASS`. QED.
