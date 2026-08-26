---
id: THM-4229
title: "Fixed-fifty nineteen-label petal Haar charts"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191/4211 + FINITE-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. Each of the twelve displayed
  nineteen-label extensions of the fixed-fifty eighteen-label chart works
  separately for every positive outsider and every nine-body. All 525,096
  new limiting profiles have strict reserve, and 227,497,842 exact literal
  checks below the twelve sufficient cutoff maxima have zero failures and
  zero equalities. This proves chi_50>=19 and twelve separate
  twenty-one-label Pascal charts. The pair-union boundary formerly left open
  here is superseded positively by THM-4234.
source: codex-jc-lrc-niche-crossfeed-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
related:
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4234-fixed-fifty-twenty-label-pair-haar-charts
script: 04-computation/lrc14_fixed_fifty_nineteen_label_petal_haar_charts_thm4229.cpp
output: 05-knowledge/results/lrc14_fixed_fifty_nineteen_label_petal_haar_charts_thm4229.out
independent_audit_script: 04-computation/lrc14_fixed_fifty_nineteen_label_petal_haar_charts_thm4229_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_fixed_fifty_nineteen_label_petal_haar_charts_thm4229_independent_audit.out
script_sha256: 9d86f5ab9918ca8f4a7bd512711edb0779f5126c28d02a62f19e0cc908e00e1b
output_sha256: 98c7f7f3a2abbb5355493aadf4f559cd97a58d45c8529421ff6b88c48c571354
independent_audit_script_sha256: 42419e4751a64f3bddd490455b7245c8f0dc8b7392cee32b5439e734abcf2adb
independent_audit_output_sha256: bbbdfa0b89dc81eb14ddbe9a5e25f0f4e7c94b16490105bae659dc554ee1d9c5
hash_basis: raw LF bytes
primary_audit: >
  PASS. A nineteen-bit exact wall arrangement and subset-zeta transform
  reconstruct every new base mass and component count, derive every
  component-discrepancy cutoff, and exhaust all admissible literal outsiders
  below it. Warning-clean O2/O3 replays byte-match the frozen output.
independent_audit: >
  ACCEPT. A clean-room combined thirty-two-bit event-toggle machine treats
  the twelve petals as separate lanes over one eighteen-bit core, groups all
  coincident walls exactly, and independently reconstructs masses and
  components. It checks all 525,096 limiting profiles, all 227,497,842 new
  literal bodies, a deduplicated inherited-chart control, and 37 direct
  body-local extremal replays. O0/O2/O3 and one-/four-thread outputs
  byte-match.
---

# THM-4229 -- fixed-fifty nineteen-label petal Haar charts

**PROVED RELATIVE TO THM-4150/4156/4170/4191/4211 + FINITE-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and exact quantifiers

For a finite positive set `S`, write

```text
G_S={x in R/Z:min_(s in S)||sx||>=1/14},
alpha=4/63.                                             (1)
```

Retain

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},

C={8,16,40,42,80,84,85,88,95,
   120,126,143,145,168,193,240,252,286},

O=P\C={10,15,20,30,60,63,132,170,176,190,264,290}.    (2)
```

For `p in O`, put `C_p=C union {p}`.

> **Nineteen-label petal theorem.** For each `p in O` **separately**, every
> positive integer `r notin P union {50}`, and every
> `K in binom(C_p,9)`,
>
> ```text
> mu(G_(K union {50,r}))>=alpha.                        (3)
> ```

The word “separately” is load-bearing: `(3)` supplies twelve different
nineteen-label witnesses. It does not assert `(3)` for a union containing two
members of `O`.

In the notation of THM-4211,
`01-canon/theorems/THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart.md`,
the proved chart-number consequence is

```text
chi_50>=19.                                             (4)
```

This is a lower bound, not a maximality statement. By THM-4150,
`01-canon/theorems/THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer.md`,
for every positive `c` and every two distinct positive odd integers `a,b`,
there is an `x in R/Z` such that

```text
min_(v in 2c(K union {50,r}) union {a,b})||vx||>=1/14. (5)
```

Multiplication by `c` preserves Haar measure, so THM-4150 applies to `cH`
when `H=K union {50,r}`. Equation `(5)` gives infinite thirteen-speed
families; it is not full LRC(14).

## 2. Inheritance pass and proof mechanism

The closest proved mechanism is THM-4211's universal eighteen-label chart,
combined with the exact component-discrepancy estimate of THM-4170,
`01-canon/theorems/THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer.md`.
The canonical hostile is THM-4207,
`01-canon/theorems/THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number.md`:
for its fixed pair `(50,51)`, marginal repair certificates need not compose,
even when each newcomer works separately. The corrected near miss is
MISTAKE-520 in `01-canon/MISTAKES.md`: the complete pool was already Haar-safe
by THM-4156, so THM-4203's cardinality-at-most-seventeen result is an
alternative `E_8` certificate, not the zero-newcomer boundary. The least-used
sidecar retained here is the exact number of positive-length circular
components of each fixed-`50` base safe set, together with literal wall
addresses below its analytic cutoff.

The active connection contract is

```text
source:       a fixed-50 base V_K=G_(K union {50})
target:       V_K intersect G_r for every second outsider r
map:          (exact mass, exact component count) -> sufficient cutoff
preserved:    labelled body K, strict surplus, circular component count
destroyed:    r-specific phase alignment above the cutoff
sidecar:      literal joint-wall integration below the cutoff
decisive test: strict surplus plus exhaustive finite remainder.           (6)
```

## 3. Strict bases, exact cutoffs, and literal remainder

The bodies in `binom(C_p,9)` which do not contain `p` are the `48,620`
bodies of THM-4211. The new universe is exactly

```text
K={p} union L,       L in binom(C,8),
|binom(C,8)|=43,758.                                  (7)
```

Across the twelve values of `p`, `(7)` gives

```text
12*43,758=525,096                                     (8)
```

labelled new base profiles.

For one such `K`, put

```text
V_K=G_(K union {50}),
mu(V_K)=n_K/D,
c_K=# positive-length circular components of V_K,
D=91,205,797,082,400.                                  (9)
```

THM-4170 gives, for every positive integer `r`,

```text
mu(V_K intersect G_r)
 >=(6/7)(n_K/D)-6c_K/(49r).                           (10)
```

After scaling the target difference by `63`, define

```text
delta_K=54n_K-4D,
kappa(K)=ceil(54c_K D/(7delta_K)).                    (11)
```

Every one of the `525,096` values `delta_K` is strictly positive. Equations
`(10)--(11)` therefore prove `(3)` whenever `r>=kappa(K)`. The exact petal
ledger is:

| `p` | minimum `delta_K` | maximizing `max kappa` body | `max kappa` | literal rows | closest literal scaled margin; `r` |
|---:|---:|:---|---:|---:|:---|
|10|508193292182004|`85,95,145,168,193,240,252,286,10`|437|405|`377192499172398/91205797082400`; 70|
|15|480254240657412|`85,88,95,145,193,240,252,286,15`|456|424|`372745529447382/91205797082400`; 70|
|20|540263639214684|`8,88,95,145,193,240,252,286,20`|455|423|`444012607831338/91205797082400`; 70|
|30|494406802675044|`85,88,95,145,193,240,252,286,30`|452|420|`399927892144788/91205797082400`; 70|
|60|517087936862640|`85,88,95,145,193,240,252,286,60`|461|429|`814811803713810/182411594164800`; 96|
|63|486559331121600|`85,95,145,168,193,240,252,286,63`|439|407|`830320769823300/182411594164800`; 96|
|132|548283516621264|`8,88,95,145,193,240,252,286,132`|461|429|`464893690731876/91205797082400`; 36|
|170|482233002468372|`8,88,95,145,193,240,252,286,170`|502|470|`406817755161174/91205797082400`; 90|
|176|532859285925684|`40,85,95,145,193,240,252,286,176`|483|451|`443732236312050/91205797082400`; 48|
|190|546296777250552|`85,88,126,145,168,193,240,286,190`|452|420|`926718153725700/182411594164800`; 96|
|264|551704390034040|`85,95,145,168,193,240,252,286,264`|495|463|`465852503091024/91205797082400`; 48|
|290|530093788467192|`85,88,95,168,193,240,252,286,290`|490|458|`895853739199104/182411594164800`; 96|

For each `p`, the literal calculation exhausts every

```text
1<=r<max_K kappa(K),       r notin P union {50},       (12)
```

and every new body `(7)`. The total finite remainder contains

```text
5,199 admissible (p,r) rows,
227,497,842 literal new-body checks,
0 failures,
0 threshold equalities.                               (13)
```

Every closest margin in the table is strictly positive. Thus `(10)--(13)`
prove the new-body part of `(3)` for every outsider. THM-4211 supplies the
bodies not containing `p`; both present implementations also recompute them
as positive controls.

## 4. Independent exact audit

The primary program constructs each nineteen-label wall arrangement, stores
exact interval lengths by allowed failure mask, and uses integer subset-zeta
transforms to reconstruct masses, starts, and continuations. It repeats the
old eighteen-label control for each petal and performs

```text
252,775,380 inherited-chart control checks.            (14)
```

The clean-room audit did not read or import the primary source. It instead
uses one thirty-two-bit simultaneous event-toggle state for

```text
C | all twelve petals | 50 | r,
```

with thirteen independent output lanes over an eighteen-bit core. Coincident
walls toggle simultaneously on a checked common integer grid. Separate
transforms recover mass and the exact `starts-bridges` component count. A
body-local event sweep replays every limiting minimizer, every cutoff
maximizer including its component count, and every literal minimizer: `37`
direct extremal checks. The one-speed identity `mu(G_{1})=6/7` with one
component is a positive control.

The independent path deduplicates the inherited chart across petals, giving

```text
22,851,400 unique inherited literal controls,
0 failures,
0 equalities.                                          (15)
```

The difference between `(14)` and `(15)` is only repeated-versus-deduplicated
accounting. Both implementations agree on every theorem statistic, labelled
extremizer, and literal minimum. The independent source uses no randomized or
unordered structure; warning-clean O0/O2/O3 builds and one-/four-thread runs
all produce the frozen byte stream.

## 5. Twenty-one-label Pascal consequence

Fix one `p in O` and one permitted outsider `r`, and put

```text
A_(p,r)=C union {p,50,r},             |A_(p,r)|=21.    (16)
```

Every eleven-subset `H` of `A_(p,r)` is Haar-safe:

- if `H` contains neither `50` nor `r`, then `H subset P`, and the full
  hereditary consequence of THM-4156,
  `01-canon/theorems/THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer.md`,
  applies;
- if `H` contains exactly one of them, THM-4191,
  `01-canon/theorems/THM-4191-complete-full-pool-newcomer-haar-transfer.md`,
  applies to the other ten labels in `P`; and
- if `H` contains both, equation `(3)` applies to its other nine labels.

The three disjoint Pascal blocks have sizes

```text
binom(19,11)+2binom(19,10)+binom(19,9)
=75,582+184,756+92,378
=352,716=binom(21,11).                                 (17)
```

Safe-set monotonicity fills every nonempty subbody of rank at most eleven.
After arbitrary positive common scaling, THM-4150 transfers each rank-eleven
body to every pair of distinct positive odd tails. This gives twelve
separate twenty-one-label ambient Pascal charts, one for each `p`, not their
union.

MISTAKE-520 is essential to the first bullet: THM-4203,
`01-canon/theorems/THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion.md`,
remains a valid alternative repair-hypergraph certificate, but its former
size-seventeen boundary was superseded by THM-4156's already-proved complete
pool inequality. No novelty or proof obligation here relies on that obsolete
boundary.

## 6. Scale-one mod-six qualifier

The Pascal count `(17)` is a certificate count, not a claim that every row is
new relative to every inherited clock. At common scale `c=1`, the single
clock `x=1/12` handles any body whose labels are all nonzero modulo `6`, as

```text
||2h/12||=||h/6||>=1/6>1/14
```

for such a body label `h`, while every odd tail has distance at least `1/12`.
The chart `C` has eleven nonmultiples and seven multiples of `6`, and `50` is
a nonmultiple. Therefore the number `N` of nonmultiples in `A_(p,r)` is

```text
N=12+1_(6 does not divide p)+1_(6 does not divide r). (18)
```

Exactly `binom(N,11)` rank-eleven bodies are closed by this one clock. The
three possible values and the complementary counts are

| `N` | fixed-clock bodies | bodies defeating this clock |
|---:|---:|---:|
|12|12|352,704|
|13|78|352,638|
|14|364|352,352|

This is explicitly a **scale-one, single-clock** comparison. It neither
persists unchanged under arbitrary common scaling nor proves that the
complementary bodies are new relative to all earlier methods.

## 7. Separate-petal firewall and superseded frontier

Nothing above certifies a nine-body containing two distinct labels of `O`.
For an unordered pair `{p,p'} in binom(O,2)`, the next candidate chart is

```text
C_(p,p')=C union {p,p'},                 |C_(p,p')|=20. (19)
```

The present theorem already handles every nine-body in this chart containing
at most one of `p,p'`. Thus the genuinely new universe for that pair is

```text
K={p,p'} union L,       L in binom(C,7),
binom(18,7)=31,824.                                    (20)
```

There are

```text
binom(12,2)=66
```

such mixed-petal pairs and

```text
66*31,824=2,100,384.                                  (21)
```

new labelled base profiles before literal outsider scans. At the time of
this theorem, proving the analogue of `(3)` for these 66 separate
twenty-label charts was the next exact mixed-petal frontier; equation `(21)`
was only a work universe. THM-4234,
`01-canon/theorems/THM-4234-fixed-fifty-twenty-label-pair-haar-charts.md`,
subsequently completed every edge and advanced the frontier to triple finite
remainders. Nothing in the present theorem alone implied that later result.

## 8. Reproduction

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_nineteen_label_petal_haar_charts_thm4229.cpp \
  -o /tmp/lrc14-thm4229-primary
/tmp/lrc14-thm4229-primary

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_nineteen_label_petal_haar_charts_thm4229_independent_audit.cpp \
  -o /tmp/lrc14-thm4229-independent
OMP_NUM_THREADS=4 /tmp/lrc14-thm4229-independent
```

Compare each byte stream with its matching frozen output in
`05-knowledge/results/`. This completes the independently audited proof.
**QED.**
