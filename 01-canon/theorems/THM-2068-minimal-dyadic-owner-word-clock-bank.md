---
id: THM-2068
title: "Minimum THM-2066 owner-word sub-bank through quotient-core maximum 24"
status: >
  PROVED by exact bitset exhaustion. Among subsets of the THM-2066 clock
  window {15,...,34}, the minimum number of clocks whose complementary-owner-
  word obstruction closes all 59,880 primitive divisor-complete eleven-cores
  in {1,...,24} is seven. One minimum bank is
  {25,26,27,28,32,33,34}. Five clocks are dominated, every bank of at most six
  of the fifteen undominated clocks fails, and each clock in the displayed
  bank has an explicit private core. Minimality is only within the stated
  twenty-clock menu. This is exact finite atlas compression, not LRC(14).
source: codex-2026-07-21-LRC-minimal-owner-clock-bank
depends_on:
  - THM-2066
related:
  - THM-2061
  - THM-2060
  - THM-775
  - HYP-2230
script: 04-computation/lrc14_minimal_owner_clock_bank_codex_20260721.py
output: 05-knowledge/results/lrc14_minimal_owner_clock_bank_codex_20260721.out
script_sha256: 0ab3b5d55a70f32876908a6bd7ea6cb5ab6cb60863286eb8fc490d5809b9ebdb
output_sha256: 91c67726a5f8bfd9c7982e7f8f36f57224eba61007867815911797a3f62f6b8c
hash_basis: normalized repository blobs (LF)
---

# THM-2068 -- minimum owner-word sub-bank

Retain the THM-2066 dyadic seam

```text
S=2C union {x,y},    |C|=11,    x,y positive odd,       (1)
```

and its clock certificate

```text
R_N(C)=empty.                                           (2)
```

Here (2) means that the eligible odd-tail owner words on the labelled
weak-safe quotient packet `A_N(C)` contain no complementary pair. By
THM-2066, (2) excludes every strict counterexample over `C`.

## 1. The finite set-cover theorem

Let

```text
U={C subset {1,...,24}: |C|=11, gcd(C)=1,
     for every d=2,...,14 some c in C is divisible by d}. (3)
```

For `15<=N<=34`, define

```text
K_N={C in U:R_N(C)=empty}.                              (4)
```

Then

```text
|U|=59,880                                               (5)
```

and the minimum cardinality of a set `F subset {15,...,34}` satisfying

```text
U=union_(N in F) K_N                                    (6)
```

is exactly seven. In particular

```text
F_0={25,26,27,28,32,33,34}                              (7)
```

satisfies (6), while no set of six clocks in the stated window does.

The qualifier about the window is load-bearing. The theorem does not claim
that seven is minimum when arbitrary clocks outside `15,...,34` are allowed.

## 2. Exact cover data and dominance reduction

The individual cover sizes are

```text
15:0,     16:3502,  17:11424, 18:1786,  19:22741,
20:3452,  21:19218, 22:9950,  23:38446, 24:6009,
25:52894, 26:47849, 27:56680, 28:40078, 29:33111,
30:4194,  31:46674, 32:29128, 33:37831, 34:37645.       (8)
```

Exact bitset inclusion gives the strict or trivial dominations

```text
K_15 subset K_27,   K_16 subset K_32,
K_17 subset K_34,   K_18 subset K_26,
K_30 subset K_27.                                      (9)
```

Therefore any bank using a left-hand clock in (9) can replace it by the
corresponding right-hand clock without increasing its size or losing a core.
It is enough to search the fifteen undominated clocks

```text
19,20,21,22,23,24,25,26,27,28,29,31,32,33,34.          (10)
```

The exact exhaustive search checks respectively

```text
C(15,1)=15, C(15,2)=105, C(15,3)=455,
C(15,4)=1365, C(15,5)=3003, C(15,6)=5005              (11)
```

banks and finds none satisfying (6). At size seven the first covering bank in
lexicographic combination order is (7), encountered after `6,409` candidates.
This proves both assertions in Section 1.

The computation also finds `5,696` distinct nonempty clock-certificate
patterns among the cores. Thus the compression is not an artifact of a tiny
number of repeated rows.

## 3. Private witnesses for the displayed bank

For each `N in F_0`, the following core lies in `K_N` but in no `K_M` for
`M in F_0\{N}`:

```text
25: {1,2,3,4,5,9,11,13,14,20,24}
26: {1,2,3,4,8,9,10,12,13,14,22}
27: {1,2,3,4,8,9,10,11,12,13,14}
28: {1,6,8,9,10,12,13,14,15,22,23}
32: {1,2,9,10,13,14,17,19,22,23,24}
33: {1,2,3,9,11,13,14,17,19,20,24}
34: {1,3,8,9,10,12,13,14,21,22,23}.                  (12)
```

Hence the particular bank (7) is irredundant. Equation (11), rather than
private witnesses alone, proves global minimum cardinality inside the menu.

## 4. Meaning for the dyadic seam

THM-2066 used every clock from `15` through `34` as a convenient interval.
The present theorem shows that thirteen of those twenty clocks are dispensable
for its bounded normalized core atlas. The surviving bank is concentrated at
the upper end and includes the consecutive block `25,26,27,28` plus
`32,33,34`. This is evidence that the owner-word obstruction is detecting
phase resolution and complementary lift ownership rather than the small
divisor-completeness clocks themselves.

It is not evidence that only seven clocks will control unbounded cores. The
predicate (2) depends on the entire labelled safe packet, and new packet words
can appear as `max(C)` grows. A lawful extension must either prove a conductor
for these word languages or attach a descent that transports a failed core to
a smaller one.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that the correct finite vertices are runners,
clocks, or tail residues. For the optimization they are **core obligations**
`C in U`; a clock is a hyperedge `K_N` covering precisely those obligations it
discharges. This quotient preserves the exact question (6), but destroys the
owner words and residue labels needed to justify membership in each hyperedge.
Those data are therefore recomputed exactly, not inferred from cover counts.

Tournament Analysis is not a faithful carrier for set cover. Orienting clocks
by cover size, with clock value as the tie Hamiltonian path, gives a transitive
tournament and no directed-cycle or SCC information. Pairwise wins cannot
express whether the union of seven hyperedges covers every obligation. The
useful fingerprints are instead dominance inclusions (9), the 5,696
certificate-pattern histogram, and the private-core witnesses (12).

## 6. Exact referee

The companion script uses Python integers as exact bitsets. It independently
rebuilds every safe packet and eligible owner word, enumerates all
`C(24,11)=2,496,144` cores, retains exactly the `59,880` members of (3), checks
that every certificate pattern is nonempty, proves each dominance inclusion,
exhausts (11), checks (7), and emits the private witnesses (12). The stored
output ends in `PASS`; the same output is obtained under `python -O`. QED.
