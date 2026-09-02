---
id: THM-4336
title: "LRC(14) two-pair rank-four surplus and cofinal one-appender completion"
status: >
  PROVED RELATIVE TO THM-4150/4170/4231/4333 + FINITE-EXACT +
  INDEPENDENT EXACT OPTIMIZER AND CLEANROOM WALL/CANDIDATE AUDIT; LRC(14)
  OPEN. For each of the THM-4333 residual pairs (50,70) and (509,640),
  every eight-label body in the fixed thirty-label pool has rank-at-most-four
  retained safe mass strictly above 7/81. Exact all-body enumeration and a
  direct-cell branch-and-bound agree on both minima. Consequently one fresh
  core appender s>=6021 or s>=5295, respectively, raises every such
  ten-label core above the 4/63 universal odd-tail threshold. Thus doubling
  the resulting eleven-speed body and adding any distinct positive odd pair
  gives a thirteen-speed safe row. This covers two named residual pairs, not
  the other 181,192 pairs, arbitrary entry, or LRC(14).
source: root + rank4_appender_probe / LRC14 continuation session, 2026-09-02
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4333-lrc14-rank-three-surplus-and-cofinal-third-tail-completion
related:
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4331-lrc14-safe-component-endpoint-denominator-odd-wall-escape
script: 04-computation/lrc14_rank4_two_appender_pair_probe_20260902.cpp
output: 05-knowledge/results/lrc14_rank4_two_appender_pair_probe_20260902.out
independent_audit_script: 04-computation/lrc14_rank4_two_appender_pair_probe_20260902_independent.py
independent_audit_output: 05-knowledge/results/lrc14_rank4_two_appender_pair_probe_20260902_independent.out
script_sha256: 91f885d72bf989ad79964cee92e0b7da0b1e1a9ee16d6df4d60c3c7ffc8e2f3f
output_sha256: bd48c1d4a3d66146f7d3ded13c911240e664ee5bef39bd2b83aa96da44627894
independent_audit_script_sha256: 4e69685eb7be35c4cdae81778a0c1d46217067893d4a425b4120caf03d085db1
independent_audit_output_sha256: fd342ff8ee0b4bb5c384209b2c779fae2fdd03d06d991dbe2bd591e73f7b3db7
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT WITH SCOPED INDEPENDENCE. A literal-midpoint C++20 wall
  program flat-enumerates all C(30,8)=5,852,925 bodies for each named pair.
  A direct weighted-cell maximum-coverage branch-and-bound on the same exact
  ledger independently returns the same optimum and least-mask tie in 1,841
  and 4,011 nodes. A no-import Python implementation reconstructs the walls
  by an unfurled two-sided midpoint test, checks the rank ledgers, remainder
  membership, both reported candidates, and the exact appender cutoffs. It
  does not independently enumerate the global optimum. C++ O2/O3/Werror and
  O0-UBSan outputs, and Python normal/optimized/hash-seeded outputs,
  reproduce their frozen transcripts.
---

# THM-4336 -- two-pair rank-four surplus and cofinal one-appender completion

**PROVED RELATIVE TO THM-4150/4170/4231/4333 + FINITE-EXACT + SCOPED
INDEPENDENT AUDIT. LRC(14) REMAINS OPEN.**

## 1. Exact rank-four statement

Retain the THM-4333 pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                       (1)
```

Fix one of the two outsider pairs

```text
Q in {{50,70},{509,640}}.                              (2)
```

Both pairs belong to THM-4231's frozen `181,194`-pair remainder. Put

```text
D=lcm{14v:v in P union Q}.                             (3)
```

Resolve all walls of `P union Q` on the integer circle of circumference
`D`. Retain each open wall cell `A` on which both labels in `Q` are safe,
write `w_A` for its integer width, and define its pool-failure set by

```text
F(A)={p in P:||px||<1/14 on A}.                        (4)
```

For `K in binom(P,8)`, define the rank-four retained mass

```text
L_4(K;Q)=sum_(A:|F(A)|<=4, F(A) intersect K=empty) w_A. (5)
```

Every retained cell lies in `G_(K union Q)`, while discarded higher-rank
cells may also be safe. Hence

```text
0<=L_4(K;Q)/D<=mu(G_(K union Q)).                      (6)
```

> **Two-pair rank-four theorem.** For every `K in binom(P,8)`,
>
> ```text
> L_4(K;Q)/D>7/81                                     (7)
> ```
>
> for both pairs in `(2)`. The exact minima are:

| `Q` | `D` | `min_K L_4(K;Q)` | `81L_4-7D` | least minimizing `K` |
|:--|--:|--:|--:|:--|
| `{50,70}` | `91,205,797,082,400` | `11,108,019,386,090` | `261,308,990,696,490` | `{10,80,95,120,143,145,168,170}` |
| `{509,640}` | `74,278,001,143,906,560` | `11,595,224,392,491,506` | `419,267,167,784,466,066` | `{8,80,88,95,145,168,170,286}` |

In reduced form, write the two minimum ratios as `m_Q`. Then

```text
m_(50,70)=100981994419/829143609840,
m_(509,640)=445970168941981/2856846197842560,          (8)

m_(50,70)-7/81
  =263948475451/7462292488560,
m_(509,640)-7/81
  =1791740033266949/25711615780583040.                 (9)
```

The root sum-of-largest-degrees screen is negative on both controls. Its
target ticks are `-431,061,626,070,252` and
`-40,171,801,599,658,044`, respectively. Thus rank-four overlap incidence,
not the marginal degree sequence, supplies `(7)`.

## 2. Exact finite proof

The primary computation constructs every literal wall in `(3)--(4)` and
classifies cells at exact integer midpoints. For failure rank at most four,
the inclusion--exclusion identity of THM-4333 becomes

```text
L_4(K)=W_4
 -sum_(i in K)d_i
 +sum_({i,j} subset K)d_ij
 -sum_({i,j,k} subset K)d_ijk
 +sum_({i,j,k,l} subset K)w_ijkl.                      (10)
```

The program evaluates `(10)` without pruning on all

```text
binom(30,8)=5,852,925                                  (11)
```

bodies for each pair. It then solves the equivalent direct weighted-cell
coverage problem by a separate exact branch-and-bound, using decreasing
marginals only as sound upper bounds. The two methods return the same
maximum covered weight and the same least-mask tie. In particular, the
strict positive ticks in the table prove `(7)` for the full finite universe,
not merely for the displayed minimizers.

As controls, the rank-zero-through-three prefixes are exactly

```text
22,084,503,304,648,
21,053,923,224,812,998,                                (12)
```

matching the audited THM-4333 packet. Direct integration of **all** safe
cells on the two rank-four minimizers gives

```text
14,349,106,369,428,
14,171,259,234,541,482.                                (13)
```

Those are full masses only for the two displayed minimizing bodies. They
are not asserted to be global full-mass minima, and no rank-at-least-five
cell is used in the proof of `(7)`.

The cleanroom Python audit imports no primary code. It independently
reconstructs the wall and rank ledgers with an unfurled two-sided midpoint
predicate, verifies both candidate masses, confirms both pairs occur in the
frozen remainder, and recomputes the transfer constants below. Global
optimizer independence comes from the direct-cell branch-and-bound in the
primary program; the Python audit does not repeat all `(11)` bodies.

## 3. One-appender LRC(14) family

The sum of the eight largest labels of `P` is

```text
290+286+264+252+240+193+190+176=1891.                 (14)
```

For `Q={q,r}`, put

```text
C_Q=1891+q+r,
gamma_Q=(6/7)m_Q-4/63.                                (15)
```

Thus

```text
Q={50,70}:   C_Q=2011,
 gamma_Q=118691847737/2902002634440,
 (6C_Q/49)/gamma_Q
   =714603342594960/118691847737<6021;

Q={509,640}: C_Q=3040,
 gamma_Q=703055796194263/9998961692448960,
 (6C_Q/49)/gamma_Q
   =3722062474903449600/703055796194263<5295.          (16)
```

> **Cofinal one-appender corollary.** For every `K in binom(P,8)`, put
>
> ```text
> s_0(Q)=6021  if Q={50,70},
>        5295  if Q={509,640}.                          (17)
> ```
>
> Then every integer `s>=s_0(Q)` satisfies
>
> ```text
> mu(G_(K union Q union {s}))>4/63.                    (18)
> ```
>
> Consequently, for every positive integer `d` and every two distinct
> positive odd integers `a,b`, the thirteen-speed row
>
> ```text
> 2d(K union Q union {s}) union {a,b}                  (19)
> ```
>
> is safe at threshold `1/14`.

**Proof.** Let `U=G_(K union Q)`. A safe set cut out by positive speeds has
at most the sum of those speeds many positive components, so `(14)` gives

```text
c(U)<=sum(K)+q+r<=C_Q.                                (20)
```

THM-4170's exact interval discrepancy inequality and `(6)--(8)` give

```text
mu(U intersect G_s)
 >=(6/7)mu(U)-6c(U)/(49s)
 >=(6/7)m_Q-6C_Q/(49s).                               (21)
```

The two ratios in `(16)` are nonintegral. Therefore `(17)` makes the final
quantity in `(21)` strictly larger than `4/63`, proving `(18)`. At the two
displayed cutoffs the exact lower bounds and surpluses are

```text
Q={50,70}:
 lower=41090163799831/647146587480120,
 lower-4/63=497192957/215715529160040;

Q={509,640}:
 lower=74714970194220413/1176544492478160960,
 lower-4/63=41197729678199/3529633477434482880.         (22)
```

Both cutoffs exceed every label in `P union Q`, so the set in `(18)` has
exactly eleven elements. Multiplication by `d` preserves Haar measure.
THM-4150 applied to `d(K union Q union {s})` therefore proves `(19)`.
Its body speeds are eleven distinct even integers and its two tails are
distinct odd integers, so no collision is hidden in the count. **QED.**

For comparison, THM-4333 gives a uniform cutoff `3,370,132,808` on every
pair in the full remainder from rank-three mass. The cutoffs `(17)` are
thousands rather than billions because `(7)` begins above `7/81`, whose
first `6/7` loss still leaves an absolute gap above `4/63`. This improvement
is proved only for the two named pairs.

## 4. Secondary two-appender wedge

The rank-four floor also yields a same-count two-appender family after one
pool label is dropped. Let `J in binom(P,7)` and extend it to an eight-set;
heredity of safe sets transfers `(7)` to `J`. Put

```text
C'_Q=1715+q+r,
delta_Q=m_Q-7/81,
epsilon_Q(s)=(6/7)delta_Q-6C'_Q/(49s),                (23)
```

where `1715` is the largest seven-label sum in `P`. If

```text
s>C'_Q/(7delta_Q),
t>max{s,(C'_Q+s)/(7epsilon_Q(s))},                    (24)
```

then two successive applications of `(21)` give

```text
mu(G_(J union Q union {s,t}))>4/63.                   (25)
```

The least integer passing the first inequality in `(24)` is `7412` for
`Q={50,70}` and `5872` for `Q={509,640}`. The second cutoff depends on
the chosen `s`; this is an ordered wedge, not a cofinal rectangle. The set
in `(25)` again has eleven elements, so doubling it and adding two distinct
odd tails gives thirteen relative speeds. Applying two appenders directly
to an eight-body would instead give fourteen relative speeds and is not an
LRC(14) corollary.

## 5. Connection contract and scope

```text
source:       pair-safe exact wall cells over the fixed thirty-label pool
target:       every eight-pool body on each of two named residual pairs
map:          discard pool-failure ranks at least five; optimize retained
              rank-at-most-four mass; apply one interval-discrepancy step
preserved:    exact retained mass and overlap incidence through rank four;
              pair labels; exact grid and normalized surplus
destroyed:    omitted full safe mass; component addresses and owners;
              all outsider pairs other than the two named controls
sidecar:      exact wall grid, failure masks, minimizing body, component
              bound C_Q, and nonintegral transfer threshold
decisive test: 81L_4(K;Q)-7D>0 for every one of 5,852,925 bodies
```

The theorem proves an infinite structured LRC(14) family, not arbitrary
entry and not LRC(14). It says nothing about the other `181,192` pairs in
THM-4333's remainder. The retained truncation is a lower bound, not the
complete safe mass. The bodies in `(13)` are not proved globally minimal.
The cutoffs are sufficient, not claimed minimal. No stochastic independence
or extrapolation from the two controls is used.

## 6. Reproduction

From the repository root:

```bash
clang++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_rank4_two_appender_pair_probe_20260902.cpp \
  -o /tmp/lrc14_rank4_two_appender_pair_probe
/tmp/lrc14_rank4_two_appender_pair_probe \
  | diff -u \
  05-knowledge/results/lrc14_rank4_two_appender_pair_probe_20260902.out -

python3 -B \
  04-computation/lrc14_rank4_two_appender_pair_probe_20260902_independent.py \
  | diff -u \
  05-knowledge/results/lrc14_rank4_two_appender_pair_probe_20260902_independent.out -
```

The audit paragraph in the frontmatter states the independence boundary.
