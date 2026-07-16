---
id: THM-904
renumber_note: was THM-902 at the first checkpoint; renumbered after the same pull exposed mac-mini's independently chosen THM-902 and boxeph's earlier HYP-7052
title: THE RESIDUE-SIX TRIPLE-CERTIFICATE REDUCTION — an explicit rational 84-weight observable on unordered three-runner sector states pointwise dominates the negative K6 miss kernel; the remaining limiting inequality is the single three-speed bound q(a,b,c) <= 47/100
status: PROVED exact pointwise reduction; FINITE-EXACT primitive triple scan through 60; universal 47/100 bound completed at VERIFIED-GRADE by THM-917, with the lambda_2/truncation remainder page still awaiting formal epsilon bookkeeping; propagation already proved exactly by THM-910
source: codex-2026-07-16-S18
depends_on: [THM-891, THM-899-lattice-law-boxhit-constants, THM-903-reflection-frame-residue6]
related: [THM-898-fourrunner-boxhit-relation-stratified, THM-905, THM-910, THM-912, THM-917, HYP-7024, HYP-7061]
verification: 04-computation/lrc14_residue6_triple_certificate_codex_S18.py -> 05-knowledge/results/lrc14_residue6_triple_certificate_codex_S18.out
---

# THM-904 — the residue-six triple-certificate reduction

This file reserves the exact theorem and filename before the long relation-tail audit.
The proved finite algebra and the still-open analytic step are deliberately separated.
`THM-910` independently closes the required limiting sign at `32/343`.  `THM-917`
subsequently completes the three-speed inequality below at verified-numerics grade;
only its explicit `lambda_2` remainder and truncation epsilon bookkeeping remain to be
written at fully formal analytic grade.

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

The proposed negation shortcut does not replace (3).  `THM-903` proves that time
reflection acts **within** each residue class, not from residue six to residue one.  In
the stationary-point frame the correct pinned reflection is `s -> 6-s`; its fixed
inner pairs are exactly `{1,5}` and `{2,4}`, the exceptional negative `K_6` support.
Thus the missing channel is reflection-even cubic data rather than a residue transfer.

## Honest boundary

- [x] explicit rational triple observable;
- [x] pointwise domination on all 462 sector-count states;
- [x] exact primitive triple scan through largest speed 60;
- [x] negation-transfer proposal decided negatively at the kernel level;
- [ ] prove `q(a,b,c)<=47/100` for all primitive triples by relation stratification;
- [ ] compose the limiting closure with a core-uniform finite-`t` wall remainder.
