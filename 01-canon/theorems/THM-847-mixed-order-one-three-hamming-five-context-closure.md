---
id: THM-847
title: Mixed order-one/order-three Hamming-five common-sheet contexts close by longest-component recursion
status: PROVED (uniform arbitrary-height closure of all 96 THM-823 mixed order-one/order-three contexts) + FINITE-EXACT (31,715-state Fraction endpoint certificate, standalone independent replay, and Tournament Analysis)
source: codex-2026-07-15-S10 frontier-object audit
depends_on: [LRC(<=13), THM-815, THM-816, THM-823]
related: [THM-837, THM-840, THM-844, THM-845, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_mixed_order_one_three_closure_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_five_mixed_order_one_three_closure_codex_S10.out
  - 04-computation/lrc13_hamming_five_mixed_order_one_three_context_closure_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_five_mixed_order_one_three_context_closure_codex_S10.out
---

# THM-847 — mixed order-one/order-three Hamming-five contexts close

Put

```text
delta=1/13,                 H={1,5,8,12} subset F_13^*.
```

Choose a multiplicative coset

```text
C=aH,
```

choose a label `b` outside `C`, and put

```text
R=C union {b},             P=[12] minus R.                 (1)
```

Give the two opposite pairs in `C` units in `{1,2}` which agree on each
pair:

```text
e_r=e_(-r) on C.                                            (2)
```

For `r in C`, let `u_r` be an arbitrary positive member of its labelled CRT
progression

```text
u_r=3r mod 13,             u_r=e_r mod 3.                  (3)
```

For the effective-order-one label, choose any `h_b>=1` and put

```text
u_b=3(b+13h_b).                                             (4)
```

## Theorem

Every packet

```text
B=3P union {u_r:r in C} union {u_b}                        (5)
```

is strictly loose:

```text
M(B)=max_t min_(w in B)||wt||>1/13.                        (6)
```

There are exactly

```text
3 quartet cosets * 8 labels outside the coset
                   * 4 opposite-pair unit words = 96      (7)
```

contexts in (1)--(4).  Thus (6) closes the complete mixed
order-one/order-three common-sheet branch isolated by THM-823, at arbitrary
lift height.

Together with THM-845 for the all-order-one branch and THM-844 for the
all-order-three branch, this closes every common-sheet survivor language in
THM-823's effective-order-at-most-twelve bank.  It does not classify
common-sheet rows with larger effective orders and does not close global
Hamming five or the `n=12` sporadic branch.

## 1. Why (5) is exactly the mixed normalized branch

THM-823 proves that if one of five replacement colours has effective order
one, then either all five orders are one or the other four orders are three
on a coset `C=aH`, with the order-one label `b` outside `C`.  Common-sheet
coverage of the quartet is possible exactly under (2).

In the mixed case the common sheet scale is three.  Divide the original AP
scale by its common gcd.  The seven retained labels become `3P`.  A quartet
replacement becomes (3), while the order-one speed is three times an integer
congruent to `b` modulo `13`, which is exactly (4).

The value `h_b=0` gives `u_b=3b`, so that coordinate has not been replaced.
It is the already closed order-three Hamming-four face of THM-816.  A proper
Hamming-five packet therefore starts at

```text
u_b=3b+39
```

and all five labelled progressions have common step `39`.  This proves that
(1)--(5) neither omit nor add a proper packet in the mixed branch.

The all-order-one alternative divides to an ordinary proper scale-one
Hamming-five packet and is closed by THM-845 (with zero-height coordinates
handled by the already closed smaller Hamming radii).

## 2. Exact recursive completeness

For a finite speed set `Q`, write

```text
E(Q)={t:||qt||>1/13 for every q in Q}.                     (8)
```

Numerically order the five replacements in a hypothetical covering row:

```text
v_1<v_2<v_3<v_4<v_5.                                      (9)
```

They are distinct because their labelled residues modulo `13` are distinct.
After a prefix, the exact state consists of

```text
(active open endpoint runs,
 remaining labelled progressions modulo 39,
 last chosen speed).                                      (10)
```

Before the last replacement a prefix has at most eleven speeds.  Settled
`LRC(<=13)` therefore gives maximin at least `1/12>1/13`, so its strict-safe
set has a positive-length component unless an earlier covering prefix has
already been recorded.

Suppose `m<=5` labels remain and a completion covers the current nonempty
residual.  It must cover a longest residual component of length `L_max`.
THM-815's one-interval discrepancy estimate gives the necessary next-speed
bound

```text
v_next <= floor(22m/[13(13-2m)L_max]).                    (11)
```

For comparison, if the residual has `K` components and total measure `L`,
the earlier global bound is

```text
v_next <= floor(22mK/[13(13-2m)L]).                       (12)
```

Since `L_max>=L/K`, (11) is no weaker than (12).

At every state the verifier enumerates, for every unused label, every member
of its progression larger than the last speed and at most the cap (11).
Thus it includes the next speed of every hypothetical covering row.  A state
with no such candidate cannot extend to a cover, while an empty residual at
any depth is exactly a covering prefix.  Numerical ordering gives every row a
unique recursive path.  This proves completeness independently of the
computed outcome.

Every residual is invariant under `t -> 1-t`.  Since only five of the twelve
labels are removed, `P` contains an even label, so `3P` contains an even speed
and `t=1/2` is not safe.  No residual component crosses the reflection fixed
point.  The replay may therefore work on `(0,1/2)`; this halves component
count and total measure, preserves `L_max`, and leaves (11)--(12) exact.

The strict-safe bands of one speed are represented literally as the open
rational intervals

```text
((13k+1)/(13u),(13(k+1)-1)/(13u)),       0<=k<u,           (13)
```

and every intersection is an exact two-pointer merge over
`fractions.Fraction` endpoints.  A common endpoint is never retained as
strict-safe interior.

## 3. Exact all-context outcome

The complete census is

| prefix depth | states | dead: no candidate |
|---:|---:|---:|
| 0 | 96 | 0 |
| 1 | 2,600 | 0 |
| 2 | 23,612 | 20,913 |
| 3 | 5,292 | 5,183 |
| 4 | 115 | 115 |

Thus

```text
total states                         31,715
covering prefixes                         0
depth-five terminal rows                  0
distinct replacement operators          378
minimum states in one context            122
maximum states in one context            733
largest longest-component cap            521
largest companion global cap             893.             (14)
```

The longest-component cap is strictly smaller than the global cap on every
one of the `31,715` visited states.  Since a hypothetical tight packet would
produce a covering path and the exhaustive tree has none, (6) follows
uniformly at every height.

The stored output emits all `96` context rows, including their CRT bases,
root geometry, node and dead-state counts, cap maxima, and a SHA-256 digest of
each recursive state word.  Their combined certificate is

```text
1cab41e8d32b09d93e9548d1baa486c0e33fb7979695d1b10725829c3a4aeb75.
```

The primary replay imports the independently frozen THM-844 interval engine
and pins its source hash.  A standalone second implementation separately
reconstructs the interval engine and the literal owner-sheet bank.  A third,
read-only audit enumerated all `63,360` marked mixed sheet presentations and
reproduced the same `96` contexts, complete state census (14), every candidate
list, and sampled closed-danger endpoint words with zero covers and zero
terminals.

## 4. Tournament Analysis and the faithful carrier

For telemetry, put the five labelled least allowed comb obligations at
tournament vertices.  Use

- root total-measure marginal erosion as the pairwise observable;
- root longest-component marginal erosion as the switch/gauge; and
- increasing `(least allowed speed,label)` as the tie Hamiltonian path.

Both gauges are scalar total orders.  In every context they therefore have
score histogram `(0,1,2,3,4)`, no directed triangle, five singleton SCCs,
and one Hamiltonian path.  Nevertheless they flip `395` edges in total, with
histogram

```text
flips per context:
{0:5,1:4,2:8,3:21,4:19,5:16,6:11,7:7,8:4,9:1}.           (15)
```

The tournament preserves a planning order only.  It destroys absolute
residual geometry, component incidence, the future progression action, and
cover truth.  The challenged assumption is that runners, residues, or even
the five combs alone should be the vertices of the theorem state.  The
faithful carrier is the evolving bipartite incidence

```text
(literal residual components) <-> (remaining labelled comb obligations), (16)
```

decorated by active endpoints and the last speed.  This is the same
operation-indexed lesson as THM-840: literal endpoints are Markov for monotone
addition, while replacement, deletion, and deck transport require the
labelled tooth ancestry.

## 5. Verification and scope

Run

```bash
python3 04-computation/lrc13_hamming_five_mixed_order_one_three_closure_codex_S10.py
python3 04-computation/lrc13_hamming_five_mixed_order_one_three_context_closure_codex_S10.py
```

and compare byte-for-byte with the stored output.  Use normal `python3`, not
`python3 -O`, because the frozen-census guards are assertions.  The frozen
hashes are

```text
parent source       6ee83b94b23f699939fea5ed53e5c50c7eb3f8a99aeac30375feabf3e0e8befb
primary source      8ab134ee85539f805378f9ea3d8e7b7829c0f6d965d172787fb7818bc1064b03
primary output      3d786b547f88673136600ca0b2a2a3ea46cbc6d0b9a94becda8a00cf88d65266
primary payload     d869bf92918ea378f45211b9564be1f9b3165371f02baa7acf4082875bdba4d2
independent source  81ed5a154b7deec91b4b9d3630e08dd9772cbeca6451b9fe37f2c942b81f1e27
independent output  48120979b9f8e81f4ab59ba533fa76b19490782c1217f3f65441dc7c1b8ea7e3
independent payload cbaa816724552a82f384ce9e56e8923120beb5b2fa40fbdb8f198eb2499b8b95
certificate         1cab41e8d32b09d93e9548d1baa486c0e33fb7979695d1b10725829c3a4aeb75
```

The replay uses no optimizer, floating point, sampled circle, or height
cutoff.  THM-847 closes the mixed order-one/order-three branch, but it does
not extend THM-823's common-sheet classification above effective order
twelve.  Those unbounded other-order languages, arbitrary-scale Hamming-six
and higher deck ramification, and the deep folded branches remain open.  In
particular, THM-847 does not prove global Hamming five, the `n=12` sporadic
branch, or LRC(14).
