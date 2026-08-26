---
id: THM-4163
title: "Order-eleven homogeneous-pair Johnson centrality"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every strong
  order-eleven tournament with a homogeneous pair has only central rational
  and exact-coset Johnson support-floor optimizers. The complete covering
  census consists of 9,355,949 strong order-ten quotient classes and all ten
  roots, hence 93,559,490 rooted presentations; all pass both gates and the
  minimum strict central-over-outer exact-coset margin is 380. Exactly 1,454
  presentations have a non-2 interior Johnson lattice, so parity rounding is
  not a valid proof shortcut. Prime order-eleven tournaments, actual response
  maximizers, target-class multiplicities, and order at least twelve remain
  outside scope.
source: codex-frontier-synthesis-creative-20260825ar
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula
  - THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4144-order-eleven-large-homogeneous-module-johnson-centrality
engine: 04-computation/tournament_order_eleven_pair_module_centrality_engine_thm4163.cpp
engine_shard0_output: 05-knowledge/results/tournament_order_eleven_pair_module_centrality_engine_shard0_thm4163.out
aggregate_script: 04-computation/tournament_order_eleven_pair_module_centrality_aggregate_thm4163.py
output: 05-knowledge/results/tournament_order_eleven_pair_module_centrality_thm4163.out
independent_audit_script: 04-computation/tournament_order_eleven_pair_module_centrality_independent_literal_audit_thm4163.cpp
independent_audit_output: 05-knowledge/results/tournament_order_eleven_pair_module_centrality_independent_literal_shard0_audit_thm4163.out
engine_sha256: 2cc8a6ef46a189db9b77c8cc929fd02659e571d0a63410c339f8949c17d8dec0
engine_shard0_output_sha256: 91240c194636daae05b7c0f81be43c862d554cadb0bb6d958a67b1affee99d85
aggregate_script_sha256: 6d5d1e0210573095d3bddaee2436f3dd5e9488365073364100468bd89a148e8b
output_sha256: 54f054d4ee6d2e5776da473aff0abb8f9c12e274789cb8bbca7be270f514abfb
independent_audit_script_sha256: f88d4b3414e2f5326b84ddea49c1f097f4409c7deb01048fac06e54d971b2d9b
independent_audit_output_sha256: 18e009fb99a10bbee7dab03d26244138aa07244c4ac1d8958b93a328af21c22e
quotient_schedule_sha256: 79c212eda5863fffb201a064398e9d425f6c21f0f14579f2d4483e160bfb2013
ordered_raw_summary_files_sha256: ec09112ccf6b7c527dcf194490628c184402f5674909a7b7f03e9530c1d68c39
boundary_manifest_sha256: dd1380ab4f6204f27c54a94a6402b8dd3c44a7f9f604da9556544e62e598fc3b
semantic_sum64: b93ec47d6452758c
semantic_xor64: 7ef3ce694e8d49fe
hash_basis: raw LF bytes
primary_audit: >
  PASS. A warning-clean quotient-once C++ engine reads nauty's complete
  strong order-ten representative stream in 1,024 residue shards, computes
  all ten rooted pair expansions by the THM-4162 tensor, and checks the exact
  rational cross-product and every Johnson coset. The aggregator verifies all
  shard files, rooted/quotient totals, 568 legitimate empty shards, the
  quotient schedule, zero failures, extrema, non-parity witnesses, two
  semantic accumulators, the ordered raw-summary digest, and a boundary-aware
  per-shard manifest. It rejects missing/extra entries, duplicate completion
  identifiers, and any departure from the exact shard grammar. Normal,
  optimized, and hash-seeded aggregation replays byte-match.
independent_audit: >
  ACCEPT. A separate literal order-eleven engine materializes each child,
  recomputes its endpoint DP, all 2^11 cut responses, moments, lattices, and
  coset floors. On the pinned 0/1024 shard it agrees row-semantically on all
  20,580 presentations with the tensor engine, including extrema, fourteen
  non-parity lattices, and zero rational/coset failures. Warning-clean builds,
  optimized/unoptimized controls, small-universe exhaustion, a fresh replay
  of all 1,024 nauty residue streams against the unsharded stream, and a
  clean-room aggregate audit leave no scope or arithmetic gap.
---

# THM-4163 -- order-eleven homogeneous-pair Johnson centrality

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and inheritance

Call `{a,b}` a homogeneous pair in a tournament `T` when every outside
vertex has the same orientation to `a` and `b`.

> **Theorem.** Let `T` be a strong tournament of order eleven containing a
> homogeneous pair. Every THM-4128 rational Johnson support-floor optimizer
> and every THM-4123 exact-coset support-floor optimizer of `T` lies on a
> central layer `t in {-1,1}`. The best central exact-coset floor exceeds the
> best outer floor by at least `380`.

The inheritance pass is:

- THM-4137 supplies every strong order-ten quotient isomorphism class;
- THM-4145 proves pair expansion preserves and reflects strongness and gives
  the exact rooted covering universe;
- THM-4162 is the closest mechanism: a parent-once endpoint tensor retaining
  the complete child cut field and exact layer lattice;
- THM-4133 is the canonical hostile: pair expansion can carry a central
  order-eleven parent to a noncentral order-twelve child, so no order-free
  inheritance statement is available;
- the least-used load-bearing sidecar is the full Johnson lattice, not its
  parity shadow. Section 6 exhibits `1,454` rows where that shadow is false.

This is support-floor centrality, not a statement about the layer containing
the actual largest ear response.

## 2. Exact covering universe

Contract a homogeneous pair in `T` to a distinguished vertex `r`. The
quotient `Q` has order ten and is strong by THM-4145. Conversely, label the
winner of the internal pair arc `a`, the loser `b`, and expand any rooted
strong `(Q,r)` to `P_r(Q)`. This recovers every target in the theorem.

THM-4137's complete order-ten generator stream has

```text
9,355,949 strong isomorphism classes.                  (1)
```

All ten roots therefore give the exact covering universe

```text
10*9,355,949=93,559,490 rooted presentations.          (2)
```

Root automorphisms and targets with several homogeneous pairs cause
duplication. That destroys target multiplicity but not a universal property;
`(2)` is not an order-eleven isomorphism-class count.

Homebrew nauty's `gentourng -q -c 10` stream was split by its exact
`res/1024` interface. The frozen schedule has `456` nonempty and `568` empty
residue shards, quotient-count digest

```text
79c212eda5863fffb201a064398e9d425f6c21f0f14579f2d4483e160bfb2013. (3)
```

The aggregator checks the count of every shard before accepting `(1)--(2)`.

## 3. Exact consequence objects and gates

For each rooted row, THM-4162 reconstructs the child Hamilton count `H`, the
`55` integer directed-exposure capacities, and the scaled packet

```text
(H,W2,D4x4,Chdx4).                                     (4)
```

Here `D4x4=4D_4` and `Chdx4=4C_hd` in THM-4128's notation. At order eleven,
the rational optimizer is central exactly when

```text
rho=2|Chdx4|/D4x4 < 1.                                 (5)
```

The engine tests `(5)` by integer cross-products, including equality as a
failure; no floating-point decision enters the proof.

For every cardinality `1<=m<=10`, THM-4162 computes from the capacity tensor

```text
A_m=sum_(|S|=m)F(S),       B_m=sum_(|S|=m)F(S)^2,
J_m=(B_m-H A_m)/(A_m-binomial(11,m)H),                 (6)
```

and the exact layer lattice `d_m` by Johnson-graph exchange gcds. With any
layer anchor `a_m`, the exact support floor is

```text
L_m=a_m+d_m ceil((J_m-a_m)/d_m),                       (7)
```

using the constant-layer convention when `d_m=0`. The exact-coset gate is

```text
max(L_5,L_6) > max_(m notin {5,6}) L_m.                (8)
```

Thus the tested objects are the theorem's actual rational and exact-coset
consequences, not proxy moments or parity-rounded approximations.

The integer widths are also certified. Partition the twelve positions of a
permutation into six adjacent pairs and independently swap within each pair.
Each orbit of size `64` contains at most one directed Hamilton path, so every
ear response satisfies

```text
F(S)<=12!/64=7,484,400.                                (8a)
```

Hence every endpoint count fits `uint32_t` (`11!<2^32`) and every layer
second moment is at most
`binom(11,5)*(12!/64)^2<2.6*10^16<2^63`. The signed
64-bit moment and packet arithmetic cannot overflow.

## 4. Complete finite-exact result

Across all `1,024` shards the exact totals are

```text
quotients                         9,355,949
rooted presentations            93,559,490
rational failures                        0
exact-coset failures                     0
presentations with inner lattice !=2 1,454.            (9)
```

The two row-order-independent semantic accumulators over the complete row
schema are

```text
sum64=b93ec47d6452758c,          xor64=7ef3ce694e8d49fe. (10)
```

Concatenating the raw shard-summary files in numerical shard order gives

```text
ec09112ccf6b7c527dcf194490628c184402f5674909a7b7f03e9530c1d68c39. (11)
```

The boundary-aware manifest

```text
dd1380ab4f6204f27c54a94a6402b8dd3c44a7f9f604da9556544e62e598fc3b (11a)
```

hashes `shard number || byte length || individual SHA-256` for every file.
The aggregator independently reconstructs `(9)--(11a)` and rejects a missing,
extra, malformed, duplicated-completion, or schedule-inconsistent shard.

The sharp rational row in this covering presentation set is

```text
child label:
1011111101111111111111111111111111110011111101110111110
quotient label:
101111110111111111111111111111110011111111111
root=4,
(H,W2,D4x4,Chdx4)=(2727,110646,3736874076,1052690016),
rho=2105380032/3736874076=175448336/311406173<1,
exact-coset margin=818.                                (12)
```

The minimum exact-coset margin occurs at

```text
child label:
1011011111111111111111111111111111111110110101110110000
quotient label:
101101111111111111111111111111111111101111110
root=3,
(H,W2,D4x4,Chdx4)=(243,23742,169212908,3515272),
rho=7030544/169212908,
central-minus-outer margin=380.                        (13)
```

Equations `(9)`, `(12)`, and `(13)` prove both strict gates and the advertised
uniform margin.

## 5. Independent literal route

The optimized engine never materializes a child `2^11` response vector: it
reuses the parent endpoint table for all roots, constructs only the
double-clone slice, and computes the moments/lattices from THM-4162. The
independent C++ audit takes the opposite route. It materializes each child,
runs a fresh order-eleven subset endpoint DP, enumerates every directed cut,
and computes every layer moment, gcd, and floor literally.

On the pinned nauty shard `0/1024`, the two implementations agree on all
`2,058` quotients and `20,580` roots. Both find zero failures, the same two
extrema, and fourteen non-parity lattice presentations. More strongly, their
full row schemas have identical commutative digests

```text
sum64=a3a1f225e6a370ee,          xor64=5ea64e1b98af2f79. (14)
```

The independent engine also rechecks that all `20,580` children are strong.
THM-4162 separately exhausts every rooted labelled quotient through order
five and compares the tensor and literal capacities entry by entry.

## 6. The non-parity lattice hostile

The tempting shortcut `d_m=2` on every interior layer survives all `202,012`
rooted labelled presentations through quotient order six and is false from
quotient order seven. The full order-eleven census finds `1,454` rooted rows
with a non-2 interior lattice. Its first frozen witness is

```text
child label:
1111111100111111100111111111111111110011110111111111111
quotient label:
111111110111111101111111111111110011101111111
root=9,                 H=21135,
bad_mask=1020,          max_inner_lattice=6.            (15)
```

Literal detail gives `d_m=6` for every `2<=m<=9` (and also for `m=1`), with
`d_10=30`. Thus ordinary odd/parity rounding would compute the wrong cosets;
the Johnson exchange lattice of THM-4162 is load-bearing even though the
centrality conclusion remains positive.

## 7. Replay and scope

Compile the warning-clean tensor and literal engines with

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/tournament_order_eleven_pair_module_centrality_engine_thm4163.cpp \
  -o /tmp/thm4163-pair-tensor

clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/tournament_order_eleven_pair_module_centrality_independent_literal_audit_thm4163.cpp \
  -o /tmp/thm4163-pair-literal
```

The pinned cross-check is

```bash
/opt/homebrew/bin/gentourng -q -c 10 0/1024 \
  | /tmp/thm4163-pair-tensor
/opt/homebrew/bin/gentourng -q -c 10 0/1024 \
  | /tmp/thm4163-pair-literal --quotient-stdin
```

For a full replay, make an isolated shard directory, run the tensor engine on
each `s/1024` for `0<=s<1024`, naming the outputs `shard_s.out`, then run

```bash
python3 -B \
  04-computation/tournament_order_eleven_pair_module_centrality_aggregate_thm4163.py \
  /path/to/shard-directory
```

Normal, optimized, and fixed-hash-seed aggregation streams byte-match the
canonical output. The front matter binds both engines, the pinned outputs,
the aggregate, schedule, full semantic accumulators, ordered raw-summary
digest, and boundary-aware manifest.

Together with THM-4144, this closes every order-eleven tournament having a
nontrivial homogeneous module. The remaining order-eleven structural stratum
is prime. This theorem does not deduplicate presentations, settle prime
tournaments, locate actual maximizing ears, prove an interval law, or extend
centrality to order twelve (which THM-4133 refutes). **QED.**
