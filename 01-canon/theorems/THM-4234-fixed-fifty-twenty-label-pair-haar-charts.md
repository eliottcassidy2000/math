---
id: THM-4234
title: "Fixed-fifty pair charts and selected triple Haar chart"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191/4229 + FINITE-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. All sixty-six fixed-fifty
  twenty-label pair charts work for every positive outsider, proving
  chi_50>=20. All 220 omitted triples have a uniform cofinal tail for
  r>=589, and the maximin-surplus triple {132,176,264} has its whole finite
  remainder exhausted, proving chi_50>=21. The exact audits cover 2,100,384
  pair profiles and 924,232,608 pair literal checks, plus 4,084,080 triple
  profiles and 8,131,032 selected-triple literal checks, with no failure or
  equality. The other 219 triple finite remainders, arbitrary finite pair
  entry, physical entry, and LRC(14) remain OPEN.
source: codex-jc-lrc-niche-crossfeed-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4229-fixed-fifty-nineteen-label-petal-haar-charts
related:
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
primary_script: 04-computation/lrc14_fixed_fifty_pair_triple_haar_charts_thm4234.cpp
pair_output: 05-knowledge/results/lrc14_fixed_fifty_pair_haar_charts_thm4234.out
all_triple_output: 05-knowledge/results/lrc14_fixed_fifty_all_triple_limiting_charts_thm4234.out
selected_triple_output: 05-knowledge/results/lrc14_fixed_fifty_selected_triple_haar_chart_thm4234.out
primary_script_sha256: ffc88dea8a8696b8984a44c86b63b0e16dcc7bcd4a7f32626787b42c17c275f9
pair_output_sha256: adbb9806d18744cd3cd7c81bfad538cfa9bc5cbf567df3078803bcd70b07707d
all_triple_output_sha256: 43abb6969efb04564c36cf735c30836458f03f763a15092b7edac3d0b52e18af
selected_triple_output_sha256: 451cd54b9fff7f9ca5cb037adeaeaff8891d96ddfc59bab2bd7f4689daa1f464
pair_independent_script: 04-computation/lrc14_fixed_fifty_pair_haar_charts_thm4234_independent_audit.cpp
pair_independent_output: 05-knowledge/results/lrc14_fixed_fifty_pair_haar_charts_thm4234_independent_audit.out
pair_independent_script_sha256: f6cc2ef03ca6f294c21a444d531af34e057d2d4825ce6aa2468bdd1378ff1223
pair_independent_output_sha256: 61d03ef2eb78e7a44716688f778ad676a7aefb73fae725458c4b03bd3f45c1c1
triple_independent_script: 04-computation/lrc14_fixed_fifty_triple_haar_charts_thm4234_independent_audit.cpp
triple_independent_output: 05-knowledge/results/lrc14_fixed_fifty_triple_haar_charts_thm4234_independent_audit.out
triple_independent_script_sha256: 334fef6d7c0923b1286dd33ed6e2cd382731dfe0c19e5a6fee1a1b42700b8f4c
triple_independent_output_sha256: cc6b24d762bfcb9108b1d4df8a93b7689958ce9509a0d9d905b9c5e3dcd9b92b
selected_independent_script: 04-computation/lrc14_fixed_fifty_selected_triple_event_scatter_thm4234_independent_audit.cpp
selected_independent_output: 05-knowledge/results/lrc14_fixed_fifty_selected_triple_event_scatter_thm4234_independent_audit.out
selected_independent_script_sha256: 513ac8af58115211f0fa8aae991390dc348135264ffa1a2708a3400589279fdd
selected_independent_output_sha256: 259630cbe7464839a4e6e1dfa46c8d8028778c8a34453fce142131be56ca2c1c
compiler_invariance_output: 05-knowledge/results/lrc14_fixed_fifty_selected_triple_compiler_invariance_thm4234.out
compiler_invariance_output_sha256: d60da714b3be614a4ac2e9636cfca1428e06f8ec602b33ce0d5e2ec0c130c681
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary exact common-grid wall integrator and two
  clean-room simultaneous-event implementations agree on every essential
  pair and triple field. A third complement-scatter implementation audits
  the selected triple. Warning-clean compiler/thread variants agree as
  described below. Two pair and eight triple cutoff-body displays choose a
  different tied maximizer; every cutoff value and essential field agrees.
---

# THM-4234 -- fixed-fifty pair charts and selected triple Haar chart

**PROVED RELATIVE TO THM-4150/4156/4170/4191/4229 + FINITE-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

For a finite positive set `S`, write

```text
G_S={x in R/Z:min_(s in S)||sx||>=1/14},
alpha=4/63.                                             (1)
```

Retain the THM-4229 partition

```text
C={8,16,40,42,80,84,85,88,95,
   120,126,143,145,168,193,240,252,286},

O={10,15,20,30,60,63,132,170,176,190,264,290},
P=C union O.                                            (2)
```

> **Pair-chart theorem.** For every unordered pair `T in binom(O,2)`, every
> positive integer `r notin P union {50}`, and every
> `K in binom(C union T,9)`,
>
> ```text
> mu(G_(K union {50,r}))>=alpha.                        (3)
> ```
>
> **Triple-tail theorem.** For every `T in binom(O,3)`, every integer
> `r>=589`, and every `K in binom(C union T,9)`, equation `(3)` holds.
>
> **Selected-triple theorem.** Equation `(3)` holds for every positive
> `r notin P union {50}` when
>
> ```text
> T=T_*={132,176,264}.                                  (4)
> ```

Consequently the fixed-fifty universal chart number of THM-4211 satisfies

```text
chi_50>=21.                                             (5)
```

The closest proved mechanism is THM-4229's twelve separate nineteen-label
petals. The canonical hostile is THM-4207: separate marginal newcomer
certificates need not compose. The corrected near miss is MISTAKE-520:
THM-4156 supplies the hereditary zero-newcomer layer; THM-4203 is an
alternative E8 certificate, not that boundary. The least-used decisive
sidecar is the exact number of positive-length circular components of each
fixed-`50` base safe set.

The theorem is complementary to THM-4231. That theorem handles arbitrary
distinct `q,r>=1290` over the whole thirty-label pool. Here one newcomer is
fixed at `50`, the body chart is restricted, and every admissible finite
second newcomer is handled for all pairs and for `(4)`. Thus this theorem
occupies a genuine residual chart on the `q=50` ray.

## 2. Exact tail lemma and endpoint convention

Fix a required petal set `T` and core body `L`, and put

```text
V_(T,L)=G_(L union T union {50}),
mu(V_(T,L))=n_(T,L)/D,
D=91,205,797,082,400.                                  (6)
```

Let `c_(T,L)` be the number of positive-length circular components after
discarding finitely many endpoint-only components. This null modification
preserves every Haar measure used below. THM-4170's periodic discrepancy
estimate gives, for every positive integer `r`,

```text
mu(V_(T,L) intersect G_r)
 >=(6/7)(n_(T,L)/D)-6c_(T,L)/(49r).                   (7)
```

Define

```text
delta_(T,L)=54n_(T,L)-4D,
kappa(T,L)=ceil(54c_(T,L)D/(7delta_(T,L)))             (8)
```

whenever `delta_(T,L)>0`. Then `(7)--(8)` prove the target inequality for
every integer `r>=kappa(T,L)`. The cutoff is sufficient, not necessary or
minimal. Equality in the cutoff comparison is lawful because the target
inequality is non-strict.

Below each cutoff the programs integrate the joint wall arrangement on the
exact grid `lcm(D,14r)`. Thus the proof has two complementary parts:

```text
strict base surplus + component count   -> cofinal analytic tail;
literal grouped-wall integration        -> finite remainder.              (9)
```

Mass alone would not suffice: it loses the cyclic adjacency controlling the
error in `(7)`.

## 3. All sixty-six pair charts

For `T={p,q} in binom(O,2)`, bodies containing at most one of `p,q` are
already covered by THM-4229. The genuinely new bodies are exactly

```text
K=L union {p,q},             L in binom(C,7).          (10)
```

The complete limiting universe therefore has

```text
binom(12,2)binom(18,7)
=66*31,824
=2,100,384                                             (11)
```

labelled profiles. Exact integration gives

```text
strict pairs                         66 / 66
nonstrict profiles                   0
largest sufficient cutoff           563, at {170,290}
admissible (pair,r) literal rows     29,042
literal new-body checks              924,232,608
literal failures / equalities        0 / 0.             (12)
```

For each pair, every `L` in `(10)` has strict `delta_(T,L)>0`; `(8)` handles
the tail and the literal audit handles all permitted `r` below the pair's
largest cutoff. THM-4229 handles bodies containing at most one petal. This
proves `(3)`.

Equivalently, put a vertex at each label of `O` and join two vertices when
their twenty-label chart is universal in the sense of `(3)`. The exact
compatibility graph is the complete graph `K_12`. The observable is
symmetric, so this graph is explicitly **not a tournament** and no cosmetic
orientation is introduced.

## 4. The triple frontier

For `T in binom(O,3)`, the only new nine-bodies beyond the pair theorem are

```text
K=L union T,                  L in binom(C,6).          (13)
```

The exact limiting universe has

```text
binom(12,3)binom(18,6)
=220*18,564
=4,084,080                                             (14)
```

profiles. Every profile has strict limiting surplus. The largest sufficient
cutoff among all 220 triples and all their bodies is

```text
589, attained for T={170,264,290}.                     (15)
```

Together with the pair layers, equations `(7)--(8)` prove the triple-tail
theorem. This is a complete cofinal three-uniform hypergraph statement, not
a completion of all 220 finite remainders.

The triple maximizing its minimum limiting surplus is exactly `(4)`, with

```text
min_L delta_(T_*,L)=687,816,435,418,308,
max_L kappa(T_*,L)=470.                                (16)
```

It is called the **maximin-surplus triple**, not the unqualified “best”
triple: for example `{10,15,30}` has the smaller maximum cutoff `354`.

The literal selected-triple census exhausts every permitted `1<=r<470`:

```text
admissible r rows                  438
literal checks                    8,131,032
failures / equalities             0 / 0
closest r                         96
closest body                      {8,95,120,145,168,286}
closest scaled margin             1,168,598,590,832,358
closest denominator               182,411,594,164,800.  (17)
```

Bodies containing zero, one, or two labels of `T_*` inherit respectively
THM-4211, THM-4229, and the pair theorem. Equation `(17)` closes the
three-petal layer. This proves the selected-triple theorem and `(5)`.

Only the displayed triple has its finite remainder exhausted here. For the
other 219 triples, strict limiting surplus proves the cofinal tail `(15)`
but makes no claim below their individual cutoffs.

## 5. Pascal chart consequences

Fix a permitted outsider `r`. For a pair `T`, put `S=C union T`, so
`|S|=20`. Every eleven-subset of `S union {50,r}` is safe:

- THM-4156 handles subsets containing neither `50` nor `r`;
- THM-4191 handles subsets containing exactly one of them; and
- the pair theorem handles subsets containing both.

The lossless layer count is

```text
binom(20,11)+2binom(20,10)+binom(20,9)
=167,960+369,512+167,960
=705,432=binom(22,11).                                 (18)
```

Thus there are sixty-six separate twenty-two-label rank-eleven chart types.
For `S_*=C union T_*`, the same argument gives

```text
binom(21,11)+2binom(21,10)+binom(21,9)
=352,716+705,432+293,930
=1,352,078=binom(23,11).                               (19)
```

The internal nine-body split witnessing `chi_50>=21` is

```text
binom(18,9)+3binom(18,8)+3binom(18,7)+binom(18,6)
=48,620+131,274+95,472+18,564
=293,930=binom(21,9).                                  (20)
```

Safe-set monotonicity fills lower ranks in the same ambient chart. THM-4150
then supplies the universal distinct odd-tail transfer after common positive
scaling for every safe eleven-body. These are large structured LRC(14)
families, not a proof of arbitrary LRC(14).

## 6. Independent audit architecture

The primary implementation uses one exact common-denominator cell ledger.
It stores core failure masks, simultaneous petal-safe masks, and fixed-`50`
safety, then uses subset-zeta transforms to recover mass and
`starts-continuations` component counts. Its three modes produce the pair
ledger, all-triple limiting ledger, and selected-triple literal ledger.

The first clean-room implementation did not read the primary source. It
uses simultaneous grouped events and signed start/bridge deposits over the
eighteen-bit core. It independently reproduces all 66 pair records and
performs 198 direct body-local extremal replays. A fresh isolated O3 run
byte-matches the frozen output. For pair lanes `{20,132}` and `{132,170}` it
chooses a different displayed cutoff body from the primary among tied
maximizers; all cutoff values and consequence-bearing fields agree.

The second clean-room implementation independently reconstructs all 220
triple limiting ledgers and the selected finite remainder. It performs 440
direct triple base replays and one direct selected-literal replay. Eight
displayed cutoff-body witnesses differ from the primary only by tied
maximizers; every cutoff value and every essential field agrees. O0/O2/O3
and one-/four-thread streams byte-match.

A third implementation audits only `(4)` through an event/complement-scatter
representation. It reproduces all values in `(16)--(17)`, adds independent
body-local controls, and byte-matches under O0/O2/O3. No implementation uses
floating point, sampling, randomized order, or an inferred finite tail.

## 7. Reproduction and scope firewall

From the repository root:

```bash
g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_pair_triple_haar_charts_thm4234.cpp \
  -o /tmp/lrc14-thm4234-primary
OMP_NUM_THREADS=4 /tmp/lrc14-thm4234-primary
OMP_NUM_THREADS=4 /tmp/lrc14-thm4234-primary --triple-base
OMP_NUM_THREADS=4 /tmp/lrc14-thm4234-primary --selected-triple

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_pair_haar_charts_thm4234_independent_audit.cpp \
  -o /tmp/lrc14-thm4234-pair-independent
OMP_NUM_THREADS=4 /tmp/lrc14-thm4234-pair-independent

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_triple_haar_charts_thm4234_independent_audit.cpp \
  -o /tmp/lrc14-thm4234-triple-independent
OMP_NUM_THREADS=4 /tmp/lrc14-thm4234-triple-independent

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic \
  04-computation/lrc14_fixed_fifty_selected_triple_event_scatter_thm4234_independent_audit.cpp \
  -o /tmp/lrc14-thm4234-selected-independent
/tmp/lrc14-thm4234-selected-independent
```

Compare each stream with its matching frozen output. The primary
`--triple-base` and `--selected-triple` streams include the common 66-row
pair prefix by design.

This theorem does not unite all twelve petals into one chart, prove
`chi_50=21`, exhaust 219 triple finite remainders, replace the fixed center
`50`, enter an arbitrary LRC(14) counterexample, or prove LRC(14). The exact
next chart test for adding a fourth omitted label has three inherited
triple lanes plus one new quadruple lane:

```text
3binom(18,6)+binom(18,5)=55,692+8,568=64,260           (21)
```

profiles per candidate fourth label. This is a work universe, not a proved
four-petal chart. **QED.**
