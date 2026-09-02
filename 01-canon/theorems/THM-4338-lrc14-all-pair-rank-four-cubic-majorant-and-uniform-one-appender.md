---
id: THM-4338
title: "LRC(14) all-pair rank-four cubic majorant and uniform one-appender completion"
status: >
  PROVED RELATIVE TO THM-4150/4170/4231/4326/4333/4336 + FINITE-EXACT +
  FULL O2/O3 LEDGER INVARIANCE + SEVEN-PAIR CLEANROOM FLAT ENUMERATION;
  LRC(14) OPEN. Every eight-label body in the fixed thirty-label pool has
  rank-at-most-four retained mass strictly above 7/81 on every one of the
  181,194 THM-4231 residual pairs. A cubic OR majorant using only degree,
  codegree, and triple incidence closes the complete census. The resulting
  uniform one-appender cutoff is 13,737, or 12,274 after importing THM-4336's
  sharper two-pair controls. This is not arbitrary entry or LRC(14).
source: root + rank4_allpair / LRC14 continuation session, 2026-09-02
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4333-lrc14-rank-three-surplus-and-cofinal-third-tail-completion
  - THM-4336-lrc14-two-pair-rank-four-surplus-and-cofinal-one-appender
related:
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
artifact_root: 05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902
artifact_manifest: 05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/SHA256SUMS
artifact_manifest_sha256: 30e17681b798ddee10410b45489dab7f0c163e9ce948a93671a94270bd59f8de
audit: >
  PASS / ACCEPT WITH SCOPED INDEPENDENCE. Full O2 and O3 exact-all runs give
  byte-identical 181,194-row ledgers. The closed verifier pins the inherited
  pair and THM-4333 rank-three hashes, every ordered row, rank prefix, target
  identity, optimizer-body control, exact/economical extrema, and all
  appender cutoffs. A no-import cleanroom program independently reconstructs
  literal walls with an unfurled midpoint predicate and flat-enumerates all
  C(30,8)=5,852,925 bodies on seven stratified pairs; O2/O3/UBSan agree with
  the primary optimum and least mask. This is optimizer and wall independence
  on seven controls, not a second global pair census.
---

# THM-4338 -- all-pair rank-four cubic majorant and uniform one-appender completion

**PROVED RELATIVE TO THM-4150/4170/4231/4326/4333/4336 + FINITE-EXACT +
SCOPED INDEPENDENT AUDIT. LRC(14) REMAINS OPEN.**

## 1. Exact all-residual-pair statement

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                      (1)
```

Let `E_rem` be THM-4231's frozen ordered universe of `181,194` unordered
outsider pairs. Its raw-LF SHA-256 is

```text
9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1. (2)
```

For `Q={q,r} in E_rem`, put

```text
D=lcm{14v:v in P union Q}.                            (3)
```

Resolve every danger wall on the integer circle of circumference `D`. Retain
the open cells `A` on which both outsiders are safe, write `w_A` for the
integer width of `A`, and put

```text
F(A)={p in P: ||px||<1/14 on A}.
W_4(Q)=sum_(A:|F(A)|<=4) w_A,
L_4(K;Q)=sum_(A:|F(A)|<=4,F(A) intersect K=empty)w_A  (4)
```

for `K in binom(P,8)`. As in THM-4333/4336,

```text
0<=L_4(K;Q)/D<=mu(G_(K union Q)).                     (5)
```

> **All-residual-pair rank-four theorem.** For every `Q in E_rem` and every
> `K in binom(P,8)`,
>
> ```text
> L_4(K;Q)/D>7/81.                                    (6)
> ```

This closes the `181,192` rank-four pairs left open by THM-4336. It applies
to exactly

```text
181,194 * binom(30,8)=1,060,514,892,450               (7)
```

labelled fixed-pool ten-speed cores. It is not an entry theorem for arbitrary
speed rows.

## 2. The cubic coverage majorant

The new mechanism is the integer polynomial

```text
g(t)=14t-9 binom(t,2)+3 binom(t,3).                   (8)
```

On the only relevant multiplicities,

```text
t             0   1   2   3   4
g(t)          0  14  19  18  14
Delta g(t)       14   5  -1  -4.                     (9)
```

Therefore

```text
14 1_(t>0)<=g(t)       for 0<=t<=4.                  (10)
```

For a body `K`, define the cubic reward

```text
G_Q(K)=sum_(rank-at-most-four A) w_A g(|F(A) intersect K|). (11)
```

If `d_i,d_ij,d_ijk` are the degree, codegree, and triple-degree tensors of
the retained rank-at-most-four weighted failure hypergraph, then

```text
G_Q(K)=14 sum_(i in K)d_i
       -9 sum_({i,j} subset K)d_ij
       +3 sum_({i,j,k} subset K)d_ijk.                (12)
```

Thus the quartic incidence tensor is not needed. Equations `(10)--(12)` give
the uniform lower bound

```text
14L_4(K;Q)>=14W_4(Q)-G_Q(K).                          (13)
```

The decreasing marginals in `(9)` make `G_Q` submodular. At any partial body
`B`, every later marginal of a remaining label is at most its current
marginal. Hence the sum of the largest required current marginals is a sound
upper bound for every size-eight completion. An exact branch-and-bound may
therefore maximize `(12)` without enumerating all `binom(30,8)` leaves;
pruning only when this upper bound is strictly below the incumbent preserves
the least-mask tie.

For each pair, the program computes

```text
M_Q=max_(K in binom(P,8)) G_Q(K),
B_Q=14W_4(Q)-M_Q,
T_Q=81B_Q-98D.                                       (14)
```

Since `14*(7/81)D=98D/81`, positivity of `T_Q` proves `(6)` simultaneously
for every body on that pair.

This polynomial is tuned to the full rank-four interval: its endpoint values
at `t=1,4` are exactly `14`, while the interior surplus makes the objective
submodular. The plain degree bound replaces `g(t)` by `14t`; it discards the
overlap rebate and fails on most of the pair universe.

## 3. Finite-exact census and audit boundary

The exact all-pair census gives

```text
pair rows                                      181,194
plain-degree target positive                    29,595
cubic optimizer target positive                181,194
cubic optimizer target nonpositive                   0
exact-all search nodes                     344,646,646
exact-all prunes                           171,098,718
maximum nodes on one pair                       8,625  at Q=(78,218). (15)
```

Equivalently, an economical replay may accept the `29,595` degree-positive
pairs immediately and run the cubic optimizer only on the remaining
`151,599`; that replay uses `270,641,887` nodes and again has no failure.

The weakest normalized cubic lower certificate is unique at `Q=(50,70)`:

```text
D=91,205,797,082,400,
W_4=34,205,623,437,384,
M_Q=357,589,704,156,176,
B_Q=121,289,023,967,200,
least maximizing mask=0138c402,
labels={10,80,95,120,145,168,170,193},

B_Q/(14D)=5,227,975,171/55,037,980,998,
B_Q/(14D)-7/81=4,244,457,985/495,341,828,982>0.       (16)
```

The displayed body maximizes the cubic majorant. It need not minimize actual
`L_4`. Its directly integrated value is only a positive control; the
universal certificate is the maximum `M_Q` in `(14)`.

The global computation uses the independently audited event-sweep wall engine
from THM-4326 and checks the incremental, tensor, and direct-cell values at
every optimizer leaf. A separate cleanroom program imports no primary wall or
optimizer code: it rebuilds literal walls using an unfurled two-sided
midpoint predicate and flat-enumerates all

```text
binom(30,8)=5,852,925                                 (17)
```

bodies on seven stratified pair controls. These include the first and last
universe pairs, both THM-4336 controls, the economical cutoff hostile, and
interior samples. All seven optimum rewards and least masks agree with the
primary optimizer. This is scoped optimizer/wall independence on seven
controls, not a second global `181,194`-pair census.

The closed verifier pins the pair-universe and THM-4333 rank-three hashes,
checks every row in order, checks every rank-prefix and tick identity, and
recomputes all normalized extrema and appender cutoffs with exact integers.

## 4. Uniform one-appender completion

For `Q={q,r}`, put

```text
m_Q=B_Q/(14D),                C_Q=1891+q+r,            (18)
```

where `1891` is the sum of the eight largest labels of `P`. For
`U=G_(K union Q)`, `(5)`, `(13)`, and THM-4170 give

```text
mu(U intersect G_s)
 >=(6/7)mu(U)-6c(U)/(49s)
 >=(6/7)m_Q-6C_Q/(49s).                              (19)
```

The least integer cutoff certified by this pairwise cubic lower bound is

```text
s_0(q,r)=floor(54C_QD/(27B_Q-28D))+1.                (20)
```

The denominator is positive precisely because `(14)` has `T_Q>0`. Exact
enumeration of `(20)` over the complete pair ledger gives

```text
max_(Q in E_rem)s_0(Q)=13,737, uniquely at Q=(50,70). (21)
```

THM-4336 independently proves sharper actual `L_4` minima for `(50,70)` and
`(509,640)`, lowering their cutoffs to `6021` and `5295`. Replacing those
cubic certificates by the proved exact minima makes the new global maximum

```text
12,274, uniquely at Q=(50,140).                       (22)
```

> **Uniform one-appender corollary.** For every `Q in E_rem`, every
> `K in binom(P,8)`, and every integer `s>=12,274`,
>
> ```text
> mu(G_(K union Q union {s}))>4/63.                   (23)
> ```

All labels in `P union Q` are at most `769`, so `s` is fresh. The set in
`(23)` has eleven speeds. Consequently THM-4150 and dilation invariance give
a safe thirteen-speed row

```text
2d(K union Q union {s}) union {a,b}                   (24)
```

for every positive integer `d` and every two distinct positive odd integers
`a,b`.

This reduces THM-4333's uniform third-core-outside cutoff
`3,370,132,808` to `12,274` while using only one new core appender. It does
not make the fixed pool into a universal entry chart.

## 5. Connection contract and remaining obstruction

```text
source:       outsider-safe wall cells over the fixed thirty-label pool
target:       every eight-pool body on every THM-4231 residual pair
map:          discard failure rank >=5; dominate rank<=4 union coverage by
              the cubic degree/codegree/triple-degree polynomial; optimize
preserved:    exact grid, pair labels, rank<=4 mass, overlap incidence
              through triples, and a uniform lower certificate
destroyed:    quartic addresses, rank>=5 safe mass, arbitrary-row entry,
              component owners, and minority-anchor location
sidecar:      ordered pair universe, exact wall ledger, maximizing body,
              pairwise component bound, and strict cutoff arithmetic
decisive test: T_Q=81(14W_4-M_Q)-98D>0 for all 181,194 pairs
```

The rank-four residual-pair frontier is closed. The two major LRC(14)
obstructions remain logically separate:

1. **Entry:** no theorem sends an arbitrary ten-speed core into this fixed
   pool/pair chart while preserving the physical safety predicate.
2. **Minority anchor:** the `420|h`, `2+12` branch of THM-4330 still needs a
   located competition or occupation theorem; pair-only transition reuse is
   unbounded after tooth addresses are forgotten.

No independence assumption, stochastic extrapolation, or inference from the
seven cleanroom controls is used in the global finite proof.

## 6. Reproduction

From the repository root:

```bash
clang++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_rank4_cubic_majorant_allpair_probe_rank4allpair_20260902.cpp \
  -o /tmp/lrc14_rank4_cubic_majorant

/tmp/lrc14_rank4_cubic_majorant \
  05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/thm4231_remainder181194.csv \
  /tmp/rank4_cubic_majorant_exactall.csv \
  /tmp/rank4_cubic_majorant_exactall.out \
  --exact-all

diff -u \
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/rank4_cubic_majorant_exactall.csv \
  /tmp/rank4_cubic_majorant_exactall.csv

python3 -B \
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/verify_rank4_cubic_majorant_packet.py \
  | diff -u \
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/verify_rank4_cubic_majorant_packet.out -

clang++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_rank4_cubic_majorant_seven_pair_flat_audit_rank4allpair_20260902.cpp \
  -o /tmp/lrc14_rank4_seven_pair_flat

/tmp/lrc14_rank4_seven_pair_flat \
  05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/thm4231_remainder181194.csv \
  | diff -u \
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/seven_pair_flat_audit.out -

(cd \
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902 \
  && sha256sum -c SHA256SUMS)
```

The primary exact-all run takes about five to six minutes on the recording
machines. Its summary includes elapsed milliseconds and is therefore pinned
by the manifest rather than expected to be byte-identical across machines;
the mathematical CSV ledger is byte-invariant. Omitting `--exact-all` runs
the faster degree-then-cubic split.
