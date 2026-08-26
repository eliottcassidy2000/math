---
id: THM-4154
title: "Mod-six fixed clock and Haar-pool inheritance correction"
status: >
  PROVED ELEMENTARY FIXED-CLOCK LEMMA + FINITE-EXACT LITERAL POOL AUDIT +
  INHERITANCE CORRECTION; NOT A NEW LRC FRONTIER. If no body label is
  divisible by 6, then x=1/12 gives every doubled body speed clearance at
  least 1/6 and every odd tail clearance at least 1/12. Hence every such
  dyadic odd-tail row has strict LRC(14) surplus at least 1/84. The P33,
  P40, and P43 pools of THM-4150/4152/4153 all satisfy this older divisor
  sieve, so the safety of every family counted there was already covered;
  the exact counts, abstract Haar transfer criteria, and exact geometry
  remain valid. Arbitrary bodies, parity-class entry, and LRC(14) are OPEN.
source: codex-lrc14-planar-jc-breakthrough-20260825
depends_on:
  - THM-366-lrc-small-denominator-divisibility-sieve
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
related:
  - THM-615-folding-identity-and-AP-even-part-confinement
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4152-second-tier-haar-finite-exception-pool40-odd-tail-transfer
  - THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
---

# THM-4154 -- mod-six fixed clock and Haar-pool inheritance correction

**PROVED ELEMENTARY FIXED-CLOCK LEMMA + FINITE-EXACT LITERAL POOL AUDIT +
INHERITANCE CORRECTION; NOT A NEW LRC FRONTIER.**

THM-366 already gives the relevant small-denominator witness: at LRC(14),
if no speed is divisible by `12`, then the unit phase `1/12` is an open
lonely witness. THM-2061 records the same obstruction after passing to the
dyadic two-odd-tail seam: a strict counterexample core must be
divisor-complete through `14`, in particular it must contain a multiple of
`6`. The lemma below specializes and stratifies the former and
constructively sharpens the latter by separating body and tail clearances.

## 1. General fixed-clock lemma

Let `H` be any nonempty finite set of positive integers satisfying

```text
6 does not divide h                    for every h in H.    (1)
```

Let `a,b` be positive odd integers. At the single physical phase

```text
x=1/12                                                       (2)
```

every doubled body speed obeys

```text
||2h x||=||h/6||.                                           (3)
```

The nonzero residue classes modulo `6` have exact distances

| `h mod 6` | `1,5` | `2,4` | `3` |
|:---:|:---:|:---:|:---:|
| `||h/6||` | `1/6` | `1/3` | `1/2` |

Therefore

```text
min_(h in H)||2h x|| >= 1/6.                              (4)
```

Every odd integer has residue `1,3,5,7,9`, or `11` modulo `12`, so

| odd residue mod `12` | `1,11` | `3,9` | `5,7` |
|:---:|:---:|:---:|:---:|
| distance at `x=1/12` | `1/12` | `1/4` | `5/12` |

Consequently

```text
||a x||>=1/12,                    ||b x||>=1/12.          (5)
```

Combining `(4)` and `(5)` proves

> **Fixed-clock lemma.** For every `H` satisfying `(1)` and every positive
> odd `a,b`,
>
> ```text
> min_(v in 2H union {a,b})||v/12|| >= 1/12
>                                    =1/14+1/84.           (6)
> ```

The body clearance is actually at least `1/6`; the full-row floor `1/12`
comes only from the odd-tail residue classes. Distinctness of `a,b` is not
needed for the inequality, but is retained in thirteen-speed applications.
This proves the lemma. **QED.**

## 2. Exact inheritance

Apply THM-366 with `n=14`, `m=12`, and unit `1 mod 12` to the speed set
`2H union {a,b}`. Condition `(1)` says that `12` divides no `2h`, and odd
tails are never divisible by `12`. Thus THM-366 already proves the common
clearance at least `1/12`; equations `(3)--(5)` identify the exact residue
mechanism and its two clearance strata.

Within THM-2061's normalized eleven-speed seam, the `N=6` clause of its
divisor-completeness conclusion says

```text
strict failure of 2H union {a,b}  =>  some h in H has 6|h. (7)
```

The contrapositive of `(7)` gives the same safety exclusion in that seam.
For arbitrary finite `H`, the direct proof and THM-366 give the stated
quantifiers. Hence `(6)` is a constructive specialization and quantitative
sharpening of existing canon, not a new reduction of the unresolved LRC(14)
core.

## 3. Literal audit of the three Haar pools

Write `P33`, `P40`, and `P43` for the explicit nested pools in
THM-4150, THM-4152, and THM-4153. The complete residue partition of `P33`
is

```text
residue 0: empty
residue 1: 1,19,25,31,43,73
residue 2: 2,8,20,32,38,50,62,80
residue 3: 51,69,75
residue 4: 4,10,16,34,40,58,64,76
residue 5: 5,17,23,29,41,47,53,71.                       (8)
```

The seven labels added to obtain `P40` have labelled residues

```text
67:1, 82:4, 86:2, 89:5, 93:3, 95:5, 141:3 mod 6,        (9)
```

and the three labels added to obtain `P43` satisfy

```text
111 mod 6=159 mod 6=285 mod 6=3.                        (10)
```

Thus all three pools have empty residue-zero class. Every eleven-subset of
each pool inherits `(1)`, and the fixed phase `(2)` proves the corresponding
odd-tail completion. The exact hereditary counts are

| pool | exact number of eleven-subsets |
|:---:|---:|
| `P33` | `binom(33,11)=193,536,720` |
| `P40` | `binom(40,11)=2,311,801,440` |
| `P43` | `binom(43,11)=5,752,004,349` |

These nested counts and the safe-set computations in THM-4150/4152/4153
remain factually correct. What is corrected is their proof significance:
the concrete pool-family safety corollaries were already covered directly by
THM-366; THM-2061's divisor pin records the same necessary condition in the
normalized seam, and `(2)` supplies the simpler common certificate. The
abstract Haar thresholds, cross-comb classifications, component-width
reductions, and compact-closed/proper-open equality mechanism are not
superseded.

## 4. Scope and information ledger

```text
source:       a dyadic body H with no label in 0 mod 6
target:       2H plus any positive odd tails
map:          h -> 2h, evaluated at the fixed phase x=1/12
preserved:    every individual residue and physical clearance
destroyed:    safe-set component geometry, Haar mass, and tail ratios
sidecar:      body labels mod 6 and tail labels mod 12
boundary:     adding one body label divisible by 6 kills this fixed clock
survivor:     THM-366/2061 may still close that row by another witness
decisive test: literal membership of H in the nonzero classes mod 6.       (11)
```

This theorem handles only the already isolated dyadic body plus odd-tail
parity seam. It does not prove entry into that seam, close bodies containing
multiples of `6`, handle mixed/even tails, or prove LRC(14).
In particular it corrects only the concrete `P33/P40/P43` novelty attribution,
not the abstract THM-4150/4152/4153 criteria, and it does not apply to
THM-4156's anchored bodies, which contain both `120` and `126`.
