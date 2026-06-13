# HYP-2444 - Pisano quotient controls the LRC14 single-stranger residual

**Status:** CONFIRMED for the one-stranger family; OPEN as a multi-stranger resource lemma.
**Source:** codex-2026-06-11-P6.
**Companions:** THM-491, THM-492, HYP-2436, HYP-2438, HYP-2443, THM-486.
**Computation:** `04-computation/lrc14_pisano_band_ladder_codex.py`;
stored output `05-knowledge/results/lrc14_pisano_band_ladder_codex.out`.

## Statement

In the exact family

```text
S(r) = 7*{1,...,12} union {r},  r <= 1092,
```

the obstruction to the plain band-1 shell horizon `q <= 27` is governed by
the antipodal unit quotient of `Z/27`.

The quotient `(Z/27)^*/+-` has 9 classes, matching

```text
pi(27) / pi(3) = 72 / 8 = 9.
```

The core `7*{1,...,12}` covers 8 of these 9 classes and misses exactly
`+-10`. Its shell-27 witness is the inverse multiplier class `+-8`.
Therefore a single stranger can block the core's last shell-27 unit twist
only by occupying the missing class `+-10`, or by collapsing to `0 mod 27`.

The exact finite computation sharpens this:

```text
plain q<=27 shell blockers in S(r): 8
their residues: r mod 13 = 0 and r mod 27 in {0, 10, 17}
old evaders after adding the B'(multiple-of-14) gate: 5
old evaders: 611, 702, 793, 962, 1053
```

Thus the residual is not just a list of exceptions. It is a two-coordinate
block:

```text
shell-27 Pisano quotient block  x  13-clock block.
```

## Stronger Family Closure

The same run found a sharper closure than the earlier band-2 phrasing.
For every valid primitive `S(r)` with `r <= 1092`:

```text
Q27 = {d*m : d | 14, m <= 27}
```

already contains a strict witness. In particular, the two rows with minimal
plain witness `q=41` are caught by the fibered shell `q=91`.

Exact row:

```text
valid rows: 936
minimal-q histogram: {13: 864, 27: 64, 40: 6, 41: 2}
Q27 coverage: 936/936
Q41 coverage: 936/936
B'(any runner) coverage: 936/936
B'(any) certificates: stranger target 792, core-runner target 144
```

So, inside the one-stranger family, the actual theorem is:

```text
Q27 alone closes the family; B'(any) is a parallel universal certificate.
```

This strengthens HYP-2438's evidence: the visible band-2 witnesses at
`q=40,41` are stress signals, but the fibered lattice already absorbs them
without enlarging the `m <= 27` horizon.

A tiny two-stranger stress probe also gives a useful negative lesson. Pairing
the 8 one-stranger plain-shell blockers over `7*{1,...,11}` produces 28
primitive rows, but all are immediately Q27-certified at `q=12`. Removing the
core runner `84` opens a low divisor clock. Thus a serious multi-stranger
search must preserve low-clock coverage while spending shell-27/Pisano
resources; simply pairing old residual residues is too easy.

HYP-2443 is the complementary marked-support view: it records which runners
pay the blocking cost at each denominator. HYP-2444 names the residue/Pisano
coordinate that the marked support is paying for in the single-stranger
residual.

## Proof Route

1. Prove the shell-27 quotient lemma by hand:
   `(Z/27)^*/+-` is one 9-cycle under multiplication by 2; the 7-core covers
   all classes except `+-10`; the inverse class `+-8` gives the core witness.
2. Convert the exact residue observation into a finite congruence lemma for
   `S(r)`: blocking all plain shells `q <= 27` forces
   `r mod 13 = 0` and `r mod 27 in {0,+-10}`.
3. Lift the closure from plain shells to the fibered lattice: show that the
   remaining two `q=41` rows admit the divisor-fiber witness `q=91`.
4. Generalize to HYP-2438 by replacing "one stranger" with a resource ledger:
   each additional stranger can block quotient classes, low clocks, or B'
   gaps, but those resources should be independent enough to bound the
   blockable band height.

## Tournament Analysis

For this session the tournament vertices are proof gates, not runners:

```text
Q41, band2, Q27, B'(any), band1, B'(mult14), shell27-residue, 13-clock.
```

The pairwise observable is

```text
(family coverage, old-evader coverage, shell-blocker coverage, -tie-rank).
```

The switch orients `A -> B` when `A` has the lexicographically larger tuple;
the tie Hamiltonian path is the declaration order. The resulting tournament
is transitive: score histogram `{0:1,1:1,...,7:1}`, no directed 3-cycles, and
one Hamiltonian path.

Assumption challenge: alternate vertex sets considered were runners, gaps,
fixed circle sections, section boundaries, residues, shell unit classes,
Pisano quotient classes, band-rung events, B' targets, and proof obligations.
The chosen proof-gate quotient preserves exact `S(r)` coverage and the residue
mechanism. It destroys arbitrary multi-stranger geometry, endpoint-owner data,
and nonlocal pressure dependencies; those are precisely the open part of
HYP-2438.
