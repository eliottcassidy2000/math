---
id: THM-902
title: THE RESIDUE-SIX TRIPLE-CERTIFICATE REDUCTION — an explicit rational 84-weight observable on unordered three-runner sector states pointwise dominates the negative K6 miss kernel; the remaining limiting inequality is the single three-speed bound q(a,b,c) <= 47/100
status: CLAIMED / CHECKPOINT STUB — exact pointwise reduction and primitive triple scan through 60 computed; independent verifier and relation-tail proof still being completed
source: codex-2026-07-16-S18
depends_on: [THM-891, THM-890, THM-898-fourrunner-boxhit-relation-stratified]
related: [THM-894, THM-896-level3-crossing, HYP-7024, HYP-7052]
verification: 04-computation/lrc14_residue6_triple_certificate_codex_S18.py -> 05-knowledge/results/lrc14_residue6_triple_certificate_codex_S18.out
---

# THM-902 — the residue-six triple-certificate reduction

This file reserves the exact theorem and filename before the long relation-tail audit.
The proved finite algebra and the still-open analytic step are deliberately separated.

Let `n=(n_0,...,n_6)` be the sector-count vector of the five moving runners and let
`M(n)={s in {1,...,6}: n_s=0}`.  Write

```text
L_6(n) = -K_6(M(n))
```

when `|M(n)|` is one or two, and `L_6(n)=0` otherwise.  There is an explicit symmetric
weight `beta:{0,...,6}^3 -> (1/100)Z`, constant on permutations, such that

```text
L_6(n) <= sum_{unordered triples of the five movers} beta(sector triple).       (1)
```

The 84 integer numerators of `beta`, in lexicographic
`combinations_with_replacement(range(7),3)` order, are stored literally in the verifier.
Equation (1) is a finite integer theorem: all 462 weak compositions of five into seven
sectors are checked, and the minimum cleared slack is zero.

For distinct positive speeds `a,b,c`, define

```text
q(a,b,c) = integral_0^1 beta(sec(a x),sec(b x),sec(c x)) dx.
```

Integrating (1) over a six-offset core gives the exact reduction

```text
-49 F_6(E) <= sum_{1<=i<j<k<=5} q(e_i,e_j,e_k).                                (2)
```

Consequently the single universal three-speed inequality

```text
q(a,b,c) <= 47/100                                                            (3)
```

would imply

```text
-F_6(E) <= 10(47/100)/49 = 47/490 < 0.097,
```

closing the sole limiting sign left by THM-891.

The current exact primitive scan covers `1<=a<b<c<=60` (`28,876` triples).  Its unique
maximum is

```text
q(1,4,7) = 81/175 = 0.462857... < 47/100.
```

This is finite-exact evidence, not yet a proof of (3).  The open analytic task is to
bound the zero-marginal cubic Fourier/relation remainder after the one- and two-coordinate
projections are discharged by THM-891's 21 pair rays.  The compact maximizer and every
floating LP support row found so far have relation height at most 21, agreeing with
THM-890 and the relation-stratified plateau in the explicit
`THM-898-fourrunner-boxhit-relation-stratified` file.

The proposed negation shortcut does not replace (3).  Reflection sends the exceptional
miss pairs `{1,5},{2,4}` to `{2,6},{3,5}`, while `K_1` is constant `-2` on every pair;
there is therefore no kernel-equivariant transfer from the negative residue-six face to
the already-closed positive residue-one face.

## Honest boundary

- [x] explicit rational triple observable;
- [x] pointwise domination on all 462 sector-count states;
- [x] exact primitive triple scan through largest speed 60;
- [x] negation-transfer proposal decided negatively at the kernel level;
- [ ] prove `q(a,b,c)<=47/100` for all primitive triples by relation stratification;
- [ ] compose the limiting closure with a core-uniform finite-`t` wall remainder.

