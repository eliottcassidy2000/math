# Session 87: faster exact counts, with the old frontier retained

**Status: exact implementation improvement; not a new Burnside formula, new
state quotient, or tournament-class enumeration algorithm.** The target is
the complete integer prefix A000568(0),...,A000568(N), or a declared set of
vertex-color/attachment prefixes. None of these routines produces tournament
representatives, automorphism transporters, or metagraph edges.

## What changed

`degree_bucket_dp` is the recommended ordinary/colored counter. It retains
Session 22's degree-local odd-factorial arithmetic and precisely the same
future-response states and transitions. It instead stores

```
states[degree][response] = integer weight
```

rather than a dictionary keyed by `(degree,response)`, and interchanges the
loops from `degree,response,multiplicity` to `degree,multiplicity,response`.
This shares each transfer's metadata across the entire degree bucket,
removes degree from hashed keys, specializes the zero-multiplicity arithmetic,
and avoids a redundant zero-filter pass for **positive** ordinary/color weights.
At the final length-one stage, it evaluates the inherited unit-tail polynomial
before applying the common transfer factor. The optional `terminal_fold=False`
ablation retains the ordinary last stage.

The second interface, `palette_prefixes`, builds an integer three-tail kernel
once and reuses it for many palettes. This is an implementation of the old
THM-514 tail ladder, not a new compression theorem. Its one-time enumerative
build loses against the response DP for some smaller query batches; measured
break-even behavior is reported rather than hidden.

## Recovery of prior work

The benchmark does not present a weak Fraction program as the previous frontier.

| Earlier mechanism | Source retained or inspected | What this session inherits |
|---|---|---|
| March 7 odd partitions and integer common denominator | `04-computation/a000568_enum.py` | Neither odd-cycle filtering nor denominator clearing is new. |
| March 8 C/GMP v4 | `04-computation/a000568_gmp_enum_v4.c` | Precomputed multiplicity factors, exponent buckets, optional threads. We rerun one thread on the same machine. |
| March 17 pure-Python integer counter | `04-computation/a000568_turbo.py` | n!-weighted exact sum; rerun only function definitions, not unrelated top-level demos. |
| June 15 unit/prime-tail ladder | `01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md` | Integer K[m,t], soluble one-tail, three-tail divisibility coordinate. |
| June 15 divisor features | `01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md` | GCD/Euler-phi energy and incremental updates; no novelty claim. |
| October 21 packed future responses | `05-knowledge/results/tournament_enumeration_20261021a/check.py` | The real state quotient and minimal common odd denominator predate this session. |
| October 21 Session 22 degree-local DP | `05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py` | **Main performance baseline**, byte-identical SHA-256 `b76167ef82dda0b134681abce33e159d5dd80d15e92edb9e31ab9b490f6939ad`. |
| October 21 colored channels | `05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/channels.py` | F_n(a,b), hence C_n(2^b), already existed. |
| Session 86 reconstruction | `05-knowledge/results/session86_relations_attachment_carriers_20261028/recreate_a000568.py`, commit `6cd9f6468387e64f7ad30af77bae47f1ff8e27b4` | Readable independent controls and frozen reference terms, not the fastest historical engine. |

The benchmark includes a telemetry-light Session 22 control: transition counts
are added once per source state, rather than incremented inside every transfer.
Its source is constructed deterministically and hash-pinned in the receipt.
This separates part of the instrumentation cost from the execution-layout gain.

## Exact arithmetic and state proof

For an odd partition lambda with multiplicities m_l,

```
E(lambda) = sum_i (l_i-1)/2 + sum_(i<j) gcd(l_i,l_j)
z(lambda) = product_l l^m_l m_l!
a_n = sum_(odd lambda partitions n) 2^E(lambda)/z(lambda).
```

An even vertex cycle has an antipodal pair that a power reverses, so it fixes
no tournament. This does **not** justify dropping even cycles for ordinary
graphs or degree-even graphs.

At descending length l, used degree d, retain the future responses
`R_j=sum_(processed L) m_L gcd(L,j)`. Every digit is at most d<=N; width
`max(1,N.bit_length())` therefore prevents carries. Adding k l-cycles gives

```
d' = d+kl
R'_j = R_j+k gcd(l,j)
Delta E = k(l-1)/2 + l*k(k-1)/2 + k R_l.
```

Let `O(d)=odd(d!)` and `v(d)=v_2(d!)`. The normalized weight is
`W=O(d) sum 2^E/z` over histories merged at that state. The transfer is

```
F = O(d') / (O(d) l^k O(k))
W' += W F c^k 2^(Delta E-v(k)).
```

F is an integer: it is the odd part of
`d'!/(d! l^k k!)`, a labelled selection/cycle-count integer. The shift is
nonnegative: `binom(k,2)>=v_2(k!)`, and all other terms are nonnegative.
Each division is checked, not silently floored. A response coordinate j is
discarded only when `j>l-2` or `j>N-d'`; later degree only increases, so such
a cycle can never return. Regrouping states by degree changes none of this.

For k=0, F=1, shift=0, and no color factor changes. **Response expiration
still runs**, since different old responses can merge even on this branch.
Dropping the zero filter is safe only for the positive ordinary/color modes
implemented here; it is not exported to the inherited signed-core modes.

At l=1 the only response is the previous cycle count t. For each d,k, the
terminal calculation is

```
Q_d(2^k) = sum_t W(d,t) 2^(kt)
total[d+k] += F c^k 2^(binom(k,2)-v(k)) Q_d(2^k).
```

This is the old soluble one-tail, scheduled to avoid one dictionary write and
one common-factor multiplication per t. Telemetry still reports the same
mathematical state and transition counts for comparison; it is not a claim
that the folded implementation executes every old operation.

## Reusable three-tail kernel

For gamma with odd parts at least five, retain its mass m, length t, and
`u=#parts divisible by 3`. Define

```
K3[m,t,u] = sum_gamma (m!/z_gamma) 2^E_gamma.
```

Adding s threes changes energy by
`s(t+2u+1)+3 binom(s,2)` and centralizer by `3^s s!`. Therefore the ordinary
core kernel is recovered integrally by

```
K[m+3s,t+s] += (m+3s)!/(m! 3^s s!) * K3[m,t,u]
               * 2^(s(t+2u+1)+3 binom(s,2)).
```

Fold the three-tail **once**, not once per palette. For r=n-m,

```
n! C_n(c) = sum_(m,t) binom(n,m) K[m,t]
                    * 2^(binom(r,2)+rt) * c^(r+t).
```

At N=100 the K3 builder visits 98,937 residual cores and has 2,049 aggregate
cells; the resulting unit-tail kernel has 834 cells. The unfactored single-n
sum has 444,793 odd partitions. Those are exact work counts, **not** a claimed
speedup factor: kernel construction, transfer work, and integer sizes matter.
The response DP can outperform this enumerative kernel at larger bounds.

## Transfers to other counting problems

* **Vertex-colored tournaments:** `C_n(c)=sum 2^E c^length/z`. This is a
  genuine colored-object count for positive integer c. It does not enumerate
  the objects or recover their automorphism groups.
* **One new distinguished vertex:** its incidences with old vertices form a
  binary color on each old vertex. Hence rooted order-(n+1) tournaments are
  counted by `C_n(2)`. The equivalent fixed-point multiplier is `m_1`.
* **k ordered distinct named ports:** at n old vertices the count is
  `2^binom(k,2) C_n(2^k)`. The factor records the arcs among the named ports.
  At total order N the independent Burnside multiplier is `(m_1)_k`.
  This counts attachment orbits globally; it does not supply parent-specific
  transporter tables, class representatives, or extension adjacency.
* **Ordinary graphs:** the analogous sum runs over **all** partitions, with
  energy `sum floor(l_i/2)+sum gcd(l_i,l_j)`. Incremental energies and degree
  grouping transfer, but the tournament odd-only filter does not.
* **Royle-even versus degree-even:** the former is the automorphism orientation-
  parity condition; the latter requires every graph degree to be even. Literal
  graph enumeration through five gives respectively `1,1,1,2,4,12` and
  `1,1,1,2,3,7`. Their first count discrepancy is order four. The Royle-even
  equinumerosity supplies another consumer of the tournament total, not a
  new graph-class generation algorithm or a degree-even count.

These identities are inherited or independently derived/checked here, not
claimed as new literature. The actual improvement is execution and reuse.

## Exact modular alternative and rejected shortcuts

`verify.py` includes a deliberately small independent CRT control. It chooses
distinct verified primes p>n, evaluates the rational sum modulo p using
invertible z_lambda, and continues until their product exceeds the **proved**
bound `c^n 2^binom(n,2)`. This certifies unique reconstruction in the labelled
upper-bound interval. A guessed bit length or an asymptotic approximation is
not a reconstruction certificate. Repeating the partition/state work for
many primes can cost more than direct integers; no CRT speedup is claimed.

Full-state integer-key packing was also tested. It keeps the same states but
adds wide-key arithmetic and did not provide a reliable stronger result. The
degree-bucket implementation, not that exploratory variant, is recommended.
No truncation, rounding, approximate pruning, or imported sequence values enter
the production calculations. The known terms are compared only afterward.

## Validation and reproduction

```sh
python3 verify.py --receipt verification.json
python3 -O verify.py --receipt verification_optimized.json
python3 benchmark.py --receipt benchmarks.json
```

The scripts locate existing repository baselines by ancestry and verify their
SHA-256 hashes before import. No historical source snapshots need to be copied
into this packet. For an external staging directory, add
`--baseline-root /path/to/repository` (an exact-hash flat staging directory is
also accepted).

All Python code uses only the standard library. The optional historical C
baseline requires `cc` and GMP; the rest does not. Benchmark CPU times use
`process_time`, walls use `perf_counter`, and methods rotate/alternate order.
The shared host can change wall time. Memory probes use `tracemalloc` in
separate untimed runs and measure Python allocations, not whole-process RSS.

Controls include exact agreement with the strongest Session 22 prefix DP
through 200, identical stage state/transition rows, independent Fraction sums
through 50, the frozen OEIS terms through 50, literal all-permutation/all-bit
tournaments through five, 508 literal ordered-pair cycle-type checks through
14 including even-cycle rejection, literal color assignments, both ordered-root
formulas, literal ordinary/Royle-even/degree-even graph classes, and certified
CRT controls. An independent reviewer also checked the counting oracle through order 50 and the claimed implementation boundaries.

### Measured results

Running the benchmark writes `benchmarks.json` with all samples, output hashes, exact work counts,
source hashes, historical reruns, palette break-even cases, and memory probes.
The concise table below is generated from that receipt after completion.

Median process CPU seconds, **five runs per main method**. Every method
computed the same full prefix and had identical stage states/transitions.

| N | Session 22 original | Session 22 telemetry-light | Degree buckets, no terminal fold | Degree buckets + terminal fold | Gain vs telemetry-light |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.04835 | 0.04868 | 0.03681 | 0.03161 | 1.54x |
| 150 | 0.36081 | 0.35676 | 0.25988 | 0.23833 | 1.50x |
| 200 | 2.36547 | 2.32844 | 1.54850 | 1.45645 | 1.60x |

At N=200 the shared mathematical work is 2,994,814 transitions, 96,836
peak states, and 51,374 precomputed degree transfers. The first improvement
column is execution layout/specialization; the further terminal-fold saving
is polynomial evaluation before common multiplication. Neither is a smaller
state quotient. These are finite CPU observations, not asymptotic complexity
or universally-fastest claims.

Tracked Python peak memory was **1,629,200 to 928,456 bytes** at N=100
and **10,013,880 to 7,108,832 bytes** at N=150 (separate untimed probes).

| N | Palettes c=1,2,4,...,2^(q-1) | Shared K3, build included | Repeat new DP | Repeat prior channel DP |
|---:|---:|---:|---:|---:|
| 80 | 4 | 0.04412 | 0.05508 | 0.08165 |
| 80 | 16 | 0.12334 | 0.23513 | 0.33637 |
| 100 | 4 | 0.16873 | 0.13572 | 0.21048 |
| 100 | 16 | 0.35999 | 0.58106 | 0.85131 |

The four-palette N=100 case is an explicit **loss** versus repeating the
new DP. At sixteen palettes, shared K3 is about 1.91x faster at N=80 and
1.61x at N=100 than repeated new DP. This is a workload-specific reuse win,
not a replacement for the strongest scalar counter.

Historical same-machine reruns (three runs each): March Python turbo needed
median 2.90477 CPU seconds for a(100); the new DP computed the **whole**
prefix through 100 in 0.03336 seconds in that separate experiment. C/GMP v4
with one thread needed median 0.13142 CPU seconds for a(100). These historical
comparisons are context only: the meaningful current-baseline gain is the
**1.50–1.60x telemetry-light comparison above**, not the much larger ratio
against March code. C/GMP was also checked at 80 and 120.

The measurement receipt retains its original driver hash; a path-resolution-only
driver update is documented there. It also records removal of one accidental
trailing newline from three staging baseline copies, with canonical Git blob
IDs. Timed algorithms did not change; repository source hashes are now checked
against canonical bytes.

## Boundary and next step

For a scalar or modest number of palettes, use `degree_bucket_dp`. For many
palettes at moderate N, consult the shared-kernel benchmark before choosing
`palette_prefixes`. Neither removes response-state or partition growth.
Counts certify how many classes exist; speeding up canonical generation and
parent/port transport is a different task, addressed elsewhere in Session 87.
