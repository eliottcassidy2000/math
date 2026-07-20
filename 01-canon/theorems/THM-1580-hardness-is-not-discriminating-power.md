---
id: THM-1580
title: "HARDNESS IS NOT DISCRIMINATING POWER — the polynomial-time arborescence count OUTSEPARATES the #P-hard Hamiltonian-path count by an order of magnitude, and their failures are nearly complementary. Closing THM-1460(C): at n = 7 there are 116 adjacency-cospectral groups; Σ_r a_r fails to split exactly 5 of them (reproducing THM-1460's 111/116 independently) and H fails to split 69 — so the intractable invariant is the WEAKER fingerprint, by 14×. The 5 survivors under Σa are NAMED here for the first time, and only 2 groups resist BOTH: (Σa, H) = (1680, 47) and (2380, 143), each of size 2. That residue is the true wall for any determinantal-plus-path fingerprint. The pattern holds at every n tested: Σa fails 0, 3, 5 groups at n = 5,6,7 while H fails 2, 16, 69."
status: >
  VERIFIED-EXHAUSTIVE at n = 5, 6, 7 over all iso classes (12, 56, 456), exact integer
  arithmetic throughout: integer char poly by Faddeev–LeVerrier, Σ_r a_r by summing
  Matrix-Tree cofactors over roots, H by direct enumeration.
  INDEPENDENT REPRODUCTION of THM-1460(C): that file reports Σa differing inside 111 of 116
  cospectral groups at n = 7; this run finds exactly 5 failures out of 116.  Agreement on a
  number computed by a separate implementation is the validation for the rest of the table.
  The naming of the 5 groups, the H comparison, and the 2-group residue are NEW.
  Advances no open problem.
source: klein-2026-07-20-S355 (owner: work on arborescence / Hamiltonian-path / logarithm leads)
depends_on:
  - THM-1460  # arborescences are the determinantal relaxation of H; (C) is what this closes
related:
  - THM-1560  # the halving dictionary -- hp is mod-2 blind, arborescences are not
  - THM-505   # H = spectral skeleton + Witt defect
script: 04-computation/arb_cospectral_residue_klein_S355.py (+ .out)
---

# THM-1580 — hardness is not discriminating power

THM-1460 established that spanning out-arborescences are the **determinantal relaxation** of
Hamiltonian paths: a Hamiltonian path from `r` is a spanning out-arborescence rooted at `r`
with every out-degree `≤ 1`, and dropping that one constraint turns a `#P`-hard count into a
determinant. Its §(C) reported that `Σ_r a_r` differs inside **111 of the 116**
adjacency-cospectral groups at `n = 7` — but did not say which 5 survive. This closes that,
and finds something sharper on the way.

## 1. The separation table

| `n` | classes | cospectral groups | `Σa` fails | `H` fails | **both fail** |
|---|---|---|---|---|---|
| 5 | 12 | 2 | **0** | 2 | 0 |
| 6 | 56 | 19 | **3** | 16 | 3 |
| 7 | 456 | 116 | **5** | 69 | **2** |

The `n = 7` row reproduces THM-1460(C)'s `111/116` from an independent implementation, which
is what validates the rest of the table.

## 2. The inversion — the finding

> **The polynomial-time invariant is the far better fingerprint.** At `n = 7`, `Σ_r a_r`
> (Matrix-Tree, `O(n³)`) fails on **5** cospectral groups; `H` (`#P`-hard) fails on **69**.
> A 14× gap, in favour of the easy one. The same ordering holds at `n = 5` and `n = 6`.

This is worth stating plainly because the intuition runs the other way: one expects the
intractable count to encode more. It does not. Computational hardness and discriminating
power are independent axes, and on adjacency-cospectral classes they point *opposite* ways.

The logarithm frame of THM-1460(D) is why they can differ at all: `log H` is additive under
ordinal sum with **no interaction term**, while `log Σa` is additive with a **size-dependent
shift**. They are not measuring the same thing, so there is no reason for their failures to
nest — and empirically they do not.

## 3. The 5 survivors under `Σa`, named

| group size | `Σ_r a_r` | `H` values | does `H` finish the split? |
|---|---|---|---|
| 2 | 1680 | `{47}` | **no — H also fails** |
| 3 | 2328 | `{131, 145}` | yes |
| 3 | 2365 | `{133, 139}` | yes |
| 2 | 2380 | `{143}` | **no — H also fails** |
| 3 | 2534 | `{153, 159}` | yes |

## 4. The residue

Only **two** groups at `n = 7` resist the whole package `(spec(A), Σ_r a_r, H)` — both of
size 2:

```text
(Σa, H) = (1680, 47)      and      (Σa, H) = (2380, 143)
```

Since the failures of `Σa` and `H` are nearly complementary — 5 and 69 failures collapsing to
2 in combination — the pair is a much stronger fingerprint than either alone. **These two
pairs are the true wall** for any determinantal-plus-path invariant, and they are the right
targets for anyone extending the reconstruction thread: whatever separates them cannot be a
function of the adjacency spectrum, the arborescence count, or the Hamiltonian-path count.

## 5. Scope

Exhaustive at `n ≤ 7`, exact arithmetic, no sampling. It advances no open problem; it closes a
named gap in THM-1460, contributes a counterintuitive fact about the two relaxations, and
isolates a two-element residue worth attacking.

*Files: `04-computation/arb_cospectral_residue_klein_S355.py` (+ `.out`).*
