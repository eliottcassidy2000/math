# Enumeration history recovered, with safe continuation targets

Source pin: `b8672af00bcc6d440a9a00a6915325faf39fcdff` (2026-10-28).
**INHERITED:** Burnside, unit/prime tails, future-response quotient, colored channels,
and root-switch theorems. **FINITE-EXACT:** the new marked-observer controls.
Historical timings below are source reports, not this session's measurements.

## 1. The actual algorithm lineage

| Date/source | Actual algorithm and output | Bottleneck / inheritance boundary |
|---|---|---|
| Mar 7–8, commits `8f81f8315104`, `d2a2deeb3745`, `aca54c0de4b8`; `04-computation/a000568_speedup.py`, `a000568_int_v2.py`, `a000568_enum.py` | Classical odd-partition Burnside sum; Fraction/divisor-profile DP; modular CRT; scaled-integer DP; finally direct grouped odd partitions with exact LCD arithmetic. | Counts, not tournament representatives. The early DP was superseded; repeated CRT traversals and huge batch weight arrays outweighed native arithmetic. |
| Mar 8, `a000568_gmp_enum_v3.c`, `a000568_gmp_enum_v4.c`, `a000568_results.txt` | GMP/C, precomputed `k^m m!`, atomic work queue, then buckets by exponent so large shifts occur only at final assembly. | Still one visit per odd partition. Source M1/8-thread timings: n100 0.2s, n150 5s, n200 200s; not portable benchmarks. |
| Apr 4, commit `4c088de5b50f` | Renewed vectorized/per-prime CRT and DP-CRT. | Source reports n100 50.9s versus GMP 0.15s; CRT is not automatically an acceleration. |
| Jun 15, THM-514 `the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`; THM-516 `burnside-core-kernel-phi-reframe.md` | Peel unit cycles; integer kernel `K[m,t]`; prime-tail divisibility statistics; totient/GCD quadratic form. | The O(n²) outer address does not make construction of its inner kernel polynomial. Full divisor profiles already encode the cycle multiplicities. |
| Oct 20, `tournament_future_boundary_20261020/README.md` | Merge histories with identical responses to every still-affordable future cycle. | Genuine avoidance of visiting every partition; keeps stage, degree, response, weight. No polynomial bound on all-order state growth. |
| Oct 21, `tournament_enumeration_20261021a`; `degree_local_tournament_enumeration_20261021/session22-enumeration` | Packed future responses; minimal uniform raw odd denominator `odd(N!)`; then degree-local `odd(d!)` normalization and precomputed integral transfers. | These are the strongest recovered scalar baselines. Same mathematical transitions, less arithmetic. The degree-local review reports 1.83–2.16× at ordinary bounds 80–200 relative to its telemetry-matched baseline. |
| Oct 21, `tournament_tail_kernels_20261021b` | Move cycle count out of the profile key into a polynomial payload, normalize `R_j-t`, fold the terminal functional before signed-excess expansion. | Fewer keys need not mean fewer coefficient cells or faster runtime. |
| Oct 27–28, Sessions79/86 | Aggregate prime cycle-index inversion/composition; root-switch operator and ordinary-child reconstruction. | Session86 counts through 50 / recomposition through 20 are independent checks, not numerical records or a new Burnside formula. |

The June negative control matters: `HYP-2523-a000568-delegates-to-odd-divisor-profiles-but-python-is-partition-bound.md`
measures n50 pairwise 0.0253s versus divisor-profile 0.0322s. At n100 an
average odd partition has only 5.666 distinct sizes. Rewriting local GCD work
was not the dominant improvement target.

## 2. Three output types must stay separate

1. **Scalar counts:** odd-partition Burnside / future-response DP. No class bank.
2. **Aggregate prime index:** `P plethysm (B-1)=1-B^{-1}-p1`.
   Retains permutation-cycle monomials; scalar specialization discards them.
   Recovering P from B and recomposing B is a structural consistency round trip,
   not an independent speed claim or prime-isotype census.
3. **Class generation:** aligned deletion parent plus attachment mask, then
   root-switch row support. `Q²=nQ`, `tr(Q)/n=rank_Q Q=rank_F2 Q=a_n`.
   It needs actual parent identifications and root transport. Counts or
   aggregate cycle indices cannot reconstruct its rows.

Source: `session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md`,
`session79_tournament_enumeration_20261027/{PRIME_CYCLE_INDEX,ROOT_SWITCH}.md`.
The n=0 count is1; Q formulas apply only to n>=1.

## 3. Useful synthesis: audit class-generation workload without generating classes

Let `C_m(c)=sum_(odd lambda|-m) 2^E c^ell(lambda)/z_lambda`.
For ordered k distinct marked vertices of an n-tournament,

```
R(n,k) = sum_(odd lambda|-n) (m1)_k 2^E/z_lambda
       = 2^binom(k,2) C_(n-k)(2^k).
```

Proof: remove the k distinguished fixed points from the permutation; the
remaining cycles choose one k-bit neighborhood pattern each, and the marked
vertices have `binom(k,2)` unconstrained mutual arcs. Equivalently,
`(m1)_k/m1!=1/(m1-k)!` reindexes the Burnside sum term by term.

This uses the existing `F_m(a,b)` channel engine at `(a,b)=(1,k)`:
`converse_tournament_count_channels_20261021/session22-converse-channels/channels.py`.
Rooted Burnside itself is already in
`session67_metagraph_enumeration_information_20261026/session67-recursion/verify_recursion.py`.
The contribution here is the connection and independent workload checks, not a
claim of a new colored-count algorithm.

In particular `p_n=R(n,1)=C_(n-1)(2)`. It gives pointed-bank size, whereas
`p_n/n` is not a_n (n3 gives 4/3, not 2). For even n the inherited binary
root-switch homology dimension is `p_n-2a_n`, obtainable without constructing Q.

New checks in `check_pointed_observers.py`: 266 direct marked/colored identities
through n=40, k<=6; 36 independent literal permutation-on-ordered-endpoint controls
through n=7; normal and `python -O` receipts byte-identical. Pointed counts n=1..10:
1,2,4,12,48,296,3040,54256,1716608,97213472. Even homology dimensions n=2..10:
0,4,184,40496,77747360. These evaluate the inherited theorem, not explicit Q matrices.

A next polynomial observer is
`sum 2^E/z product_l(1+t^l)^m_l`, counting pointed children by root outdegree.
For k ports, `sum_(c in {0,1}^k) x_c^l` retains joint neighborhood patterns.
Bounded searches found earlier fixed-group subset polynomials but no exact
same aggregate implementation; absence is not a priority proof. This is an
actionable payload extension, not a measured speedup.

## 4. Related sequences: reuse the operation, not an unjustified formula

March `burnside_enum_v2.c` already shared partition traversal, GCD interactions,
denominator factors and exponent buckets across tournaments, graphs, digraphs,
relations, and converse families. For simple graphs allow all cycle lengths
and replace the local term by floor(l/2); for digraphs double cross-GCD
interaction and use l-1 internally; loops add one cycle factor. Transfer the
future-response quotient when this pairwise interaction law is preserved,
but re-prove its integer normalization and orientation-monodromy filters.

Strong tournament counts already came from `S=1-1/A` in March
`a051337_from_a000568.py`: SCCs form an ordered sequence, not an unordered SET.
The current six-channel packet also supplies self-converse, strong-self-converse
and converse-merged counts. For anti-actions, endpoint-swap parity matters;
the Oct21 correction to THM-283 supersedes its older pair-orbit explanation.

The historical k-uniform hypergraph accelerator
`04-computation/k_uniform_fast_enum.c` uses divisor-signature Möbius inversion
instead of traversing every power of a long permutation. That mechanism is
transferable when only a bounded selection of cycle blocks is active. Its
current implementation is not safe over its advertised full k<=20 range:
see the explicit failure below.

## 5. Recovered hazards and repair boundary

* **Actual legacy parity bug:** `a000568_speedup.py`, blob
  `b6fb5dc726b70ee8a07c129ef7090255e750d800`, normalizes max_part before
  taking min(remaining,max_part). Its n5 generator emits
  `[5],[3,2],[1,1,1,1,1]`, and `a000568_davis(5)=10`, not12.
  Executed against exact source this session. Fix parity after taking min.
  This does not validate its unrelated legacy CRT path.
* **No exact rounding from relative asymptotics:** saved
  `a000568_asymptotic_exact_s26.out` fails n14/n15 and overflows float.
  At n30 even the single omitted 3-cycle contribution exceeds1 in absolute
  count. A relative-small tail is not an integer-rounding certificate.
* **Actual hypergraph truncation:** current C blob
  `3f890cd423631a3f8659c63d149e1790c08db2d5` accepts k<=20 but silently
  keeps only 10000 ordered compositions. k15 needs 16384. For cycle type
  `(3,2,1^12)` at n17 it drops `(1,2,12)` and the capped formula gives
  92 fifteen-subset orbits, versus93 by complementing to unordered pairs.
  Reproduced both by independent JavaScript evaluation and a compiled C harness.
  The independent audit additionally verifies that reversing the representation
  of these same cycle blocks restores 93: truncation causes representation-order
  dependence. Literal traversal of all 136 subsets at sizes 2 and 15 gives 93 each.
  Fail-fast capacity validation is the minimal safe repair. The proposed patch
  rejects the unsupported nontrivial capacity domain, and rejects malformed,
  overflowing integer arguments and k<1/n<0/threads<1 before allocations.
  Other bucket caps and fixed-width arithmetic remain unaudited; do not claim
  a general repair.

Immediate priorities: test/fix the parity generator; prohibit silent
composition truncation; benchmark against degree-local/terminal-fold baselines;
reuse the existing color channels for marker/workload audits; preserve full
partition payloads when the requested output is a prime cycle index.



## 6. Verified repair handoff

The two unified patches are `a000568_speedup.patch` and
`k_uniform_capacity_guard.patch`; they target the actual legacy repository
paths, not the modern counting engine. No Git writes were made by this lane.

`check_repository_legacy_enumeration.py` is the publication-ready,
repo-relative regression (standard library plus GCC/GMP). It needs no whole
original/fixed source snapshots. With the patches applied, run:

```sh
python3 PATH/TO/check_repository_legacy_enumeration.py --repo .
python3 -O PATH/TO/check_repository_legacy_enumeration.py --repo .
```

Prepublication runs use the documented source overrides. Both runs pass with
byte-identical JSON receipts: odd-partition bijection and direct counts 0..30
against the frozen degree-local kernel; Dirichlet counts 0..20; complete C
composition counts for every k=1..14; 17 small exact C controls; 10 unsupported
capacity/invalid-input rejections. The C probe also reproduces the old 92/93
failure using only a tiny explicitly historical capped generator. These tests
validate only the stated repairs, not the old CRT implementation or all
integer-overflow/bucket bounds in the hypergraph implementation.

Source SHA256 values are in `repository-repair-receipt.json`; the unchanged
Session22 baseline is `b76167ef82dda0b134681abce33e159d5dd80d15e92edb9e31ab9b490f6939ad`.
