# Aliquot-fiber square-sum quotients: shared ports and exact path lifting

Status: **PROVED** elementary obstruction and relation-composition repair; **FINITE-EXACT** minimality certificates and census. No novelty or literature-priority claim. No theorem about tournaments, aliquot-chain termination, or infinite arithmetic amplification.

## Inheritance and scope

This connects the user's inverse-aliquot fibers to the square-sum graph without confusing a square-deficit representation with a square-sum edge. Inherit the finite-fiber bound and `F13={27,35}` from [Session78](../session78_fourth_power_aliquot_20261027/ALIQUOT.md) and [ALIQUOT_FIBERS.md](ALIQUOT_FIBERS.md), and Session81 [aliquot observers](../session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md) and [observer grammar](../session81_counterexample_families_20261028/OBSERVER_GRAMMAR.md)'s descent test: an observable supports an operation precisely when the relevant output is constant on its fibers. The shared-port/clone-merger discipline supplies the hostile question: do separate existential witnesses use the **same named source** where two constraints meet? This is a shared logical mechanism, not a map from integers to tournaments.

Live board: complete inverse fibers / existential graph quotient / shared-source incidence / relational composition / observer descent. Closest mechanism: loss under a scalar quotient. Canonical hostile: 27 and 35 have the same aliquot sum. Corrected near miss: nonempty edges separately need not compose. Least-used coordinate: the actual source at the common fiber.

Let `s(n)` be the sum of positive proper divisors of positive n, and `F_t={n:s(n)=t}`. Only labels `t>=2` are considered. Such n are composite. If p is their least prime factor, `n/p>=sqrt(n)` is a proper divisor distinct from 1, so `s(n)>=1+sqrt(n)`. Therefore

`n <= (t-1)^2`.

For labels 2..N, a divisor sieve through `(N-1)^2` gives **complete fibers**, not truncated source samples. F1 is infinite and is excluded; F0 is also outside this claim.

Use the simple undirected source graph: distinct positive x,y are adjacent iff x+y is a perfect square. On a specified union of fibers, define the existential quotient on nonempty fiber labels: distinct a,b are adjacent iff some source in F_a is adjacent to some source in F_b. Internal fiber edges are irrelevant to a quotient path with distinct labels. A lift of the ordered distinct-label path `(t0,...,tk)` means sources `xi in F_ti` with every consecutive sum square. Distinct labels guarantee distinct sources.

## 1. Smallest three-label false path

The **complete** fibers are

`F4={9}, F10={14}, F13={27,35}`.

The induced source graph on their union has exactly two edges:

`9--27` (sum36), and `14--35` (sum49).

The other four unordered sums are 23, 44, 41, 62, none square. Its quotient is the path `4--13--10`, but there is no lift: the first edge requires middle source 27, the second 35. Thus the quotient connects two components of this induced graph. This is **not** a claim that 9 and 14 lie in different components of the graph on all positive integers, or of a larger union of fibers.

Even adjacency to the named source 9 fails to descend through s:

`s(27)=s(35)=13`, yet `9+27=36` is square and `9+35=44` is not.

The first failed implication is `(exists x,y E(x,y)) and (exists y,z E(y,z)) => exists x,y,z [E(x,y) and E(y,z)]`, where the two left-hand witnesses for y are independent.

**FINITE-EXACT minimum:** among three distinct labels >=2, the least possible maximum label of a quotient length-2 path with no lift is 13. Up to reversing the endpoints, the unique witness at that maximum is `(4,13,10)`. The complete fiber table needed to exclude all smaller maxima is:

| t | F_t |
|---|---|
|2,5|empty|
|3|4|
|4|9|
|6|6,25|
|7|8|
|8|10,49|
|9|15|
|10|14|
|11|21|
|12|121|
|13|27,35|

The checker exhausts all potential middle labels and unordered pairs of distinct endpoint labels using these complete fibers. This is a certified minimum for this exact predicate, not for every possible quotient failure.

## 2. Exact local repair, and why it must be propagated

For distinct labels a,b, retain the relation

`R_ab = {(x,y) in F_a x F_b : x+y is square}`.

Its shared-fiber incidence set is `P_b(a)={y in F_b : exists x in F_a, (x,y) in R_ab}`. Then

**`a--b--c` lifts iff `P_b(a) intersection P_b(c)` is nonempty.**

This follows by taking the same y for both edges. For the minimal witness the ports are `{27}` and `{35}`. Positive control: `4--13--14` lifts via 9,27,22, and also 9,27,169.

Merely testing every consecutive length-2 subpath is still insufficient for longer paths. An **actual arithmetic**, complete-fiber hostile is

`33--14--37--12`.

Its fibers and full consecutive relations are:

- `F33={45,87,247}`, `F14={22,169}`, `F37={155,203,299,323}`, `F12={121}`;
- `R_33,14={(87,169)}`;
- `R_14,37={(22,203),(169,155)}`;
- `R_37,12={(203,121)}`.

The first triple lifts as 87,169,155 (sums256,324); the second as 22,203,121 (sums225,324). But the whole word cannot lift: after its first edge the middle relation forces 169 to 155, while the final edge needs 203. These relations are all pairs, not selected edges with other possible lifts omitted.

**FINITE-EXACT minimum:** among four distinct labels >=2, with both consecutive three-label subpaths liftable but the full path not liftable, least maximum label is 37, with this unique witness up to reversal. The checker exhausts complete labels 2..100, hence in particular every smaller candidate. At N=100 there are 4,777 length-2 failures and 1,415 locally liftable length-3 failures, in both cases counting paths modulo reversal, not merely unordered sets of labels. These census counts are finite, not asymptotic.

## 3. Exact all-length repair

For an arbitrary fiber word with adjacent labels different, propagate named-source sets:

`S0=F_t0`,

`S_(i+1)={y in F_t(i+1): exists x in S_i with (x,y) in R_ti,t(i+1)}`.

Induction shows S_i consists exactly of possible last sources of lifts of the prefix. Consequently the whole word lifts iff S_k is nonempty. Equivalently, the Boolean matrix product of the consecutive R relations has a nonzero entry. Existentially replacing each matrix by one bit does **not** commute with composition: nonzero factors can have zero product.

If all labels are distinct this decides a simple source path. With repeated labels it decides a **walk**, and a simple-path requirement needs extra used-source memory; no claim of simple-path sufficiency in that case is made. Prescribed starting/ending sources are handled by singleton initial states and a final membership test. The full relation sidecar is sufficient; minimality of its storage is not asserted. General abstract disjoint-copy amplification is unnecessary here and would not by itself prove an arithmetic family of complete fibers.

The connection contract is explicit: source = named integers and square-sum adjacency; summary = aliquot labels and existential edges; preserved predicate = existence of each individual cross-fiber edge; destroyed information = compatibility of its endpoint with a neighboring edge; repair = source-resolved incidence and relation composition; cheapest decisive test = the four-source graph above. There is no intrinsic antisymmetric relation to justify a tournament construction.

## Reproduction and controls

Standalone Python 3 standard library, no external files or packages:

```
python3 check_fiber_ports.py
python3 -O check_fiber_ports.py --output /tmp/fiber_ports_optimized.json
cmp fiber_ports_results.json /tmp/fiber_ports_optimized.json
```

Default N=100 scans all sources 1..9801 twice, once by divisor sieve and once by independent divisor-pair enumeration. Adjacency is checked both by square generation and direct integer-square tests. Every length-2 and length-3 census candidate is checked against source-set propagation. An independent Cartesian-product lift oracle checks every four-distinct-label word through label 13, plus the long hostile: 11,881 words. Tests use explicit exceptions, not removable assertions. `--max-label` supports 37..500 (a resource guard, not a mathematical bound). Output is deterministic JSON beside the script unless overridden. Generated JSON is a reproducibility receipt, not a dependency. Normal and optimized runs agree byte-for-byte.

## Independent audit

A separate trial-divisor enumeration of every source1..9801 reproduced all complete fibers, direct square tests reproduced all edge relations, and independent relation composition reproduced both unique minimum witnesses and both census counts. Normal and optimized outputs agree. As an explicit ambient-scope control,9–27–22–14 is a valid source path once fiber14 is allowed: the earlier false-connectivity claim concerns only the selected induced union, not the entire square-sum graph.
