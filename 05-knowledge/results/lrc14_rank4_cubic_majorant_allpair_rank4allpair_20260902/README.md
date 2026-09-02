# LRC(14) all-residual-pair rank-four cubic-majorant packet

**Status: FINITE-EXACT computation supporting a theorem candidate.  It closes
the fixed-pool rank-four residual-pair census, not arbitrary entry and not
LRC(14).**

## Statement tested

Use the fixed thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.
```

Let `E_rem` be the exact ordered `181,194`-pair remainder of THM-4231.  Its
raw SHA256 is

```text
9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1.
```

For `Q={q,r} in E_rem`, resolve the walls of `P union Q` on circumference

```text
D=lcm{14v:v in P union Q}.
```

Retain cells on which `q` and `r` are safe.  If `F(A) subset P` is the pool
failure mask and `w_A` the integer width, put

```text
W_4(Q)=sum_(|F(A)|<=4) w_A,
L_4(K;Q)=sum_(|F(A)|<=4, F(A) intersect K=empty) w_A
```

for `K in binom(P,8)`.  The exact result is

```text
for every Q in E_rem and every K in binom(P,8),
L_4(K;Q)/D > 7/81.
```

Thus the packet extends THM-4336 from its two named pairs to all `181,194`
pairs in the frozen remainder.

## Cubic coverage majorant

If a retained cell has `t=|F(A) intersect K|` selected failures, then
`0<=t<=4`.  The integer cubic

```text
g(t)=14t-9 binom(t,2)+3 binom(t,3)
```

has the table

```text
t       0   1   2   3   4
g(t)    0  14  19  18  14
Delta      14   5  -1  -4.
```

Consequently

```text
14 1_(t>0) <= g(t),
```

and the decreasing increments make the induced weighted set function
submodular.  In incidence coordinates,

```text
G_Q(K)=14 sum_(i in K)d_i
       -9 sum_({i,j} subset K)d_ij
       +3 sum_({i,j,k} subset K)d_ijk.
```

Writing `M_Q=max_(|K|=8) G_Q(K)`, cellwise majorization gives, uniformly in
the body,

```text
14 L_4(K;Q) >= 14W_4(Q)-M_Q.                         (1)
```

This is the mechanism: rank-two and rank-three incidence force rebates for
multiple hits, so one does not need a quartic tensor to certify the rank-four
floor.  At a partial body, future marginals can only decrease.  The sum of
the largest required current marginals is therefore an admissible completion
upper bound for an exact cardinality-eight branch-and-bound.  Pruning is
strict below the incumbent, retaining the least-mask tie.

The cubic is a lower certificate for `L_4`; it is not asserted to equal the
coverage indicator.  The ledger's `direct_l4` is evaluated only at the body
maximizing `G_Q` and is a positive control, not a global `L_4` minimum.

## Complete results

The economical run first applies the raw union degree bound and invokes the
cubic optimizer only when that bound is inconclusive:

```text
pairs                                      181,194
degree-positive                             29,595
exact cubic fallbacks                      151,599
fallback positive                          151,599
fallback nodes                         270,641,887
fallback prunes                        134,308,799.
```

The frozen ledger is the stronger exact-all replay, which optimizes the cubic
on every pair:

```text
exact cubic positive                       181,194
nonpositive                                      0
nodes                                  344,646,646
prunes                                 171,098,718.
```

The weakest normalized exact-all cubic certificate is uniquely

```text
Q=(50,70),
D=91,205,797,082,400,
14W_4-M_Q=121,289,023,967,200,
(14W_4-M_Q)/(14D)=5,227,975,171/55,037,980,998,
surplus over 7/81=4,244,457,985/495,341,828,982.
```

The economical split also passes everywhere.  Its weakest normalized
certificate is `(499,636)`; this is a weaker degree certificate, not the
exact-all minimum.

## One-appender consequence

For each pair let

```text
m_Q=(14W_4-M_Q)/(14D),
C_Q=1891+q+r,
```

where `1891` is the largest eight-label sum in `P`.  THM-4170 gives

```text
mu(G_(K union Q union {s}))
 >= (6/7)m_Q-6C_Q/(49s).
```

The least integer cutoff certified by this row is

```text
s_Q=floor(54 C_Q D/(27(14W_4-M_Q)-28D))+1.           (2)
```

The maximum of `(2)` over the full ledger is `13,737`, uniquely at
`Q=(50,70)`.  THM-4336 already gives the sharper cutoff `6021` on that pair
(and `5295` on `(509,640)`).  Replacing those two rows by the stronger proved
values makes the uniform hybrid maximum

```text
s>=12,274,
```

uniquely at `Q=(50,140)`, whose strict pre-ceiling ratio is

```text
56,939,779,118,542,320 / 4,639,067,430,787.
```

Therefore, relative to THM-4150/4170/4336, every `Q in E_rem`, every
`K in binom(P,8)`, every integer `s>=12,274`, every positive integer `d`, and
every two distinct positive odd integers `a,b` give the safe thirteen-speed
row

```text
2d(K union Q union {s}) union {a,b}.
```

The cutoff exceeds every pool and outsider label (`max r=769`), so this set
has exactly eleven even body speeds and two distinct odd tails after transfer.

## Audit layers

1. The full exact-all ledger contains one ordered row per inherited pair,
   including the grid, all rank-zero-through-four masses, coarse bound, exact
   cubic certificate, maximizing body, optimizer-body direct `L_4`, and search
   totals.
2. The closed Python verifier pins the universe and ledger hashes, matches
   every rank-zero-through-three prefix against THM-4333, checks every integer
   identity and strict target, recomputes both census splits, orders normalized
   extrema with exact fractions, and derives every appender cutoff.
3. The independent C++ audit imports no primary source.  It reconstructs
   literal walls with an unfurled two-sided midpoint test and flat-enumerates
   all `C(30,8)=5,852,925` bodies, without pruning, at universe indices
   `0,20910,60000,120000,163336,181103,181193`.  These are pairs
   `(1,2)`, `(50,70)`, `(147,466)`, `(321,612)`, `(509,640)`, `(721,746)`,
   and `(766,768)`.  Every optimum reward and least mask matches the primary
   ledger.  O2, O3, and O1-UBSan transcripts are identical.
4. Full O2 and O3 exact-all runs produced byte-identical ledgers with SHA256
   `60dab8a471065dee132e61a5695e6a827616e0e45e7c6a67cd11b426fd86623a`.

## Reproduction

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
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/verify_rank4_cubic_majorant_packet.py

clang++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_rank4_cubic_majorant_seven_pair_flat_audit_rank4allpair_20260902.cpp \
  -o /tmp/lrc14_rank4_seven_pair_flat

/tmp/lrc14_rank4_seven_pair_flat \
  05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/thm4231_remainder181194.csv \
  | diff -u \
  05-knowledge/results/lrc14_rank4_cubic_majorant_allpair_rank4allpair_20260902/seven_pair_flat_audit.out -
```

Omitting `--exact-all` runs the economical degree-then-cubic split.  The
primary full replay takes about six minutes on the recording machine; the
seven-pair flat audit takes a few seconds.

## Scope and connection contract

```text
source:       pair-safe wall cells over the fixed thirty-label pool
target:       every eight-pool body on all 181,194 residual pairs
map:          replace the rank-four OR coverage cost by the cubic concave
              majorant g, then maximize the resulting submodular function
preserved:    exact grid, pair labels, rank<=3 incidence, uniform lower bound
destroyed:    exact rank-four coverage, higher-rank safe mass, cyclic address,
              owner and arbitrary-row entry
sidecar:      full pair ledger, maximizing body, direct optimizer-body control,
              component bound and exact cutoff arithmetic
decisive test: 81(14W_4-M_Q)-98D>0 for every residual pair
```

The computation proves a complete finite fixed-pool statement.  It does not
show that rank-four cells exhaust the safe set, identify the global full-mass
minimum, treat pairs outside `E_rem`, supply arbitrary-row entry, or prove
LRC(14).
