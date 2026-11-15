# Nonlocal chronological parent choices: a complete bounded target compiler

Here every newborn has exactly one backward edge to **any** old parent; no root-stability or last-vertex requirement is imposed. The reversed-edge graph is genuinely a tree. This is not the arbitrary abstract tree without its chronological ordering.

## Why an exact target never requires unbounded chronological padding

In the ordinal SCC order C_1→…→C_r of T, let q∈C_i be the selected parent. If q is the sink (terminal singleton), the new vertex lies just before it and no strong component changes. Otherwise the newborn and the entire suffix C_i,…,C_r form a single strong component: x→q reaches all those vertices, all of them except q beat x, and q can reach another suffix vertex. Earlier SCCs dominate this suffix and x, so are unchanged.

Starting from the singleton, there is therefore always a transitive prefix and at most one nontrivial, **terminal** strong component. Before any cyclic step, all births had to choose the original sink to remain transitive; their chronology has reversed-edge parents all equal to 0. At the first cyclic step, some of that original prefix is swallowed by the new final component. Later births can only enlarge or absorb more of this prefix. Any transitive prefix left at the end consists solely of old precyclic vertices, none ever used as a parent by a survivor (using it would swallow it). Remove them without losing the original vertex 0, which was the sink and is already in the selected strong suffix. The retained induced tournament still has a literal chronological one-parent history, and preserves H, since deleting the ordinal transitive prefix divides by 1. All count-1 cases have the singleton witness.

Consequently a proposed N>1 is realizable in this **entire** grammar iff it has a chronological strong witness on at most floor(1+3 log_5 N) vertices, by [Busch's bound](https://doi.org/10.37236/1141). This is a genuine all-order finite inverse, not a fixed-order convention. It leaves open whether this larger grammar covers the whole [Chan–Pak spectrum](https://arxiv.org/abs/2308.10214).

## Literal nonlocal parents at the three test values

A parent list is zero-indexed: its entries choose the parent of newborns 1,2,… respectively.

| Count | Minimum order in this grammar | One physical parent list |
|---:|---:|---|
| 23 | 6 | 0,0,0,1,1 |
| 33 | 6 | 0,0,0,1,2 |
| 35 | 7 | 0,0,1,0,4,0 |

The direct census visits all 4,037,914 chronological parent histories through order eleven, including those with transitive prefixes. It computes exact weighted starts by subset DP through order ten, applies the proved one-leaf formula to every order-eleven candidate, and also checks a deterministic sample of the last layer by independent full DP. At every odd target through **365**, it finds precisely the two exclusions 7 and 21 in the full physical leaf language. Exact all-order completeness of this finite decision range follows from the structural reduction and 5^(11/3)>365, not extrapolation. This explains why missing 23,33,35 under a local choice rule is a restriction on that rule, not on the full chronological tree of physical histories.

## Old transitive padding cannot hide a unit branch

If T has a nontrivial terminal strong core C and an old transitive prefix, every port inside C has its same defect as in C (all starts are confined to C and the unused prefix contributes Hamilton count 1). For a prefix port, let R be its nonempty later suffix containing C. A start at that port consists of the port followed by a chosen path of R, and the even starts require an odd choice of vertex set in R. Odd-cycle packing gives cancellation of the signed two-slot sum, so its doubled even weight is Q_2(R). Hence the defect is Q_2(R)−H(C)≥3H(C), because splitting a full path at every cut gives Q_2(R)≥(|R|+1)H(C), with |R|≥3. Such a vertex cannot be a unit port. This also shows arbitrary neutral early padding does not add a hidden branch to the chronological mixed-port census. In particular the [generic through-order-eight check](MOVABLE_PORTS.md) excludes a full depth-two unit binary subtree from any chronological source of H≤73: Busch puts its strong core at order at most eight, and padding cannot supply unit ports. For H=1 the scalar descendant 7 is already forbidden.

Replay: build `replay_parent_compiler.cpp` with C++17 and compare its printed receipt.
