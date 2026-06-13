---
source: opus-2026-06-02-S568 (remote-control)
status: PROOF PROGRAM + unified target (verified across even n) + per-n routes; leverages HYP-2091 parity. Not full proofs.
tags: [LRC, proof-program, even-n, floor-tight, transversal, parity, HYP-2091, pinch, shield, spectral-gap, n14]
---

# A proof program for LRC at even n up to 14: the floor-tight-transversal target

**Prompt (user):** leverage HYP-2091; spend a long session; plan creatively;
read the repo; try to get proofs for each even n up to 14.

This is the careful plan. It (1) states the **single unified target** that every
project lens reduces to, (2) **verifies** it across even n=4..14, (3) gives the
**cleanest route per even n** using HYP-2091's parity structure, and (4) lays out
the **creative attack on n=14**. Honest: 4–12 are theorems in the literature; the
genuine deliverable is the unified localisation and the n=14 route, not new full
proofs of 4–12.

## 0. Status of the even cases (n = total runners)

| n | status | literature method |
|---|---|---|
| 4 | PROVED | Betke–Wills/Cusick (geometry of numbers) |
| 6 | PROVED | Bohman–Holzman–Kleitman (averaging) |
| 8 | PROVED | Rosenfeld 2025 (finite-check sieve) |
| 10 | PROVED | Trakulthongchai 2025 (sieve) |
| 12 | PROVED | Sungkawichai–Trakulthongchai 2026 (sieve+poly) |
| **14** | **OPEN** | — |

## 1. The unified target (every lens collapses here) — VERIFIED

`M(S)=max_t min_i‖v_i t‖`; LRC(n) ⟺ `M(S) ≥ 1/n`. Split by measure:

- **positive measure ⟹ `M>1/n` ⟹ lonely.** The IGNORE regime: all incommensurate
  (Weyl), all low-resonance, the round/orbit/transitive side. Trivial.
- **measure zero = the WORRY set.** Here `M=1/n` (tight) **or** `M<1/n`
  (counterexample). So:

> **LRC(n) ⟺ every measure-zero speed set is FLOOR-TIGHT (`M=1/n`), never below.**

**Verified (`lrc_unified_target_even_n_s568.py`, exhaustive small boxes):** for
n=4,6,8,10,12,14, **0 counterexamples** — every measure-zero config has `M=1/n`,
and (except the n=8 sporadics) is a **perfect antipodal transversal mod `2n-1`**:

| n | primitive sets | measure-0 (worry) | counterexamples | perfect transversal |
|---|---|---|---|---|
| 4 | 79 | 1 | **0** | 1/1 |
| 6 | 461 | 2 | **0** | 2/2 |
| 8 | 1716 | 3 | **0** | 1/3 (+2 sporadic) |
| 10 | 5005 | 1 | **0** | 1/1 |
| 12 | 12376 | 1 | **0** | 1/1 |
| 14 | 8568 | 1 | **0** | 1/1 |

So the worry-set is **tiny and structured** (1–3 configs per box), and the target
is: *prove the worry-set never dips below `1/n`.* Every project thread is a
different description of this same set:

| lens | the worry-set is… |
|---|---|
| measure (S564) | the measure-zero configs |
| resonance/orbit (S563) | the resonance-maximal closed orbits |
| Burnside (S565) | the dual-Burnside **fix side** (self-converse) |
| strong (S566) | the **regular rotational encirclement** |
| polygon (S567/**HYP-2091**) | the even-n **polygon/dihedral** face (m odd) |
| spectral gap (oracle S552, HYP-2084) | the **perfect antipodal transversals** mod 2n-1 |

## 2. Leverage of HYP-2091 (the parity carrier)

HYP-2091/S567: even n ⟹ `m=n-1` odd ⟹ the worry-set lives on the **clean polygon
ladder** (rotational `R_m`, dihedral `D_{2m}`); `n→n+2` preserves it. Consequence
for the program:

- The even-n worry-set is *uniformly* the rotational encirclement — one shape, all
  even n, climbing the `n+2` ladder. A proof that works at one rung should be
  `n+2`-recursive.
- The transversal description (`mod 2n-1`) and the polygon description (`R_{n-1}`)
  are the **outside/necklace** (A000016) reading of the same object — the worry-set
  is "all outside, no mesh." So the target only needs the *boundary* data.

## 3. The cleanest route per even n

- **n=4, 6 (elementary range).** The worry-set is 1–2 transversals (AP, and the
  `{2}`-flip sporadic at n=6). Route: the **pinch lemma** (S557) gives `M=r/s`; for
  these the binding pair sums to `n` with `gcd=1`, so `s=n`, `r=1`, `M=1/n` exactly
  — *provided no other runner is closer*, i.e. **no multiple of `n`** (THM-369 /
  C′). Closing n=4,6 = proving the 1–2 explicit transversals are floor-tight, a
  finite check the literature already implies. **Reproducible, not new.**
- **n=8 (composite `2n-1=15`).** Worry-set = 1 transversal **+ 2 non-transversal
  sporadics** (`{1,2,3,4,5,7,12}`…). These escape the transversal description (the
  non-unit-pair hole, S559). Route needs the **non-transversal closure** — exactly
  the gap that recurs at n=14.
- **n=10, 12 (composite/prime `2n-1`).** Worry-set = 1 transversal (the AP) in box;
  proven in the literature by the sieve. Route: the AP is floor-tight via the
  `(1,n-1)` pinch at `t=1/n` (no multiple of `n`).
- **n=14 (the frontier).** Worry-set in box = 1 transversal (AP); V* (the
  distance-1 sporadic, speed 24) and the non-transversal sporadics live just
  outside. The full worry-set = perfect transversals mod 27 **∪** non-unit-pair
  non-transversals. Route below.

## 4. The creative attack on n=14

Three converging sub-attacks, all aimed at "worry-set is floor-tight":

**(A) The transversal branch (spectral-gap chain, oracle S552/S553, HYP-2084).**
Reduce to the `2^{n-1}` perfect antipodal transversals mod `2n-1`; show each has
`M=1/n` (not below). HYP-2084 verified this in bounded boxes (every below-`2/(2n-1)`
row is floor-tight). *Open:* prove for all transversals. Parity lever: each
transversal is a flip-set `F` off the AP (S553/S572); prove `M(F)=1/n` by the
unit-shell witness clock (a missed unit antipodal shell gives a `2/(2n-1)`
witness, so only the *all-lower* transversal = AP and the rare sporadic stay at the
floor).

**(B) The non-transversal branch (the composite-`2n-1` hole, S559/S563-apex).**
`2n-1=27=3³` composite ⟹ non-unit antipodal pairs ⟹ sporadic tight non-transversals
(as at n=8). Close with a **second witness clock**: lift `mod 2n-1 → mod p(2n-1)`
to restore invertibility of the non-unit shell (S559 handoff), or the **r/p
higher-lift** of the pinch sieve (S562). This is the `2q`-apex residual (HYP-2063).

**(C) The strong/encirclement branch (S566).** Show the observer always reaches
sole-source; the only doubtful configs are forced into the **regular rotational
encirclement** (the worry-set), and a regular rotational encirclement of `m=n-1`
(odd) runners always leaves a `≥2/n` gap. Moon's theorem frames it: the runners
Hamiltonian-cycle around the observer; uniform (regular) spacing is the extremal,
and it touches the floor but never closes (the `1/n` wall).

**The one lemma that would close it:** *a perfect antipodal transversal mod `2n-1`
(and its non-unit-pair cousins) has `M(S)=1/n` exactly.* (A) gives the unit part,
(B) the non-unit part, (C) the geometric meaning. The pinch lemma supplies the
`r/s` form; the shields (THM-396/397) supply the obstruction calculus.

## 5. The honest gap and why it is small

Every lens has **localised** LRC(14) to: *the perfect-transversal-and-cousins
worry-set is floor-tight.* That set is `2^{13}`-ish transversals plus a thin
non-transversal layer — a **finite, highly-symmetric** target (regular rotational
polygons), not the continuum. The remaining work is the two clock-lemmas (A unit-
shell, B non-unit lift). This is dramatically smaller than "check all speed sets,"
and it is `n+2`-uniform (HYP-2091), so a proof at the transversal level would cover
**all even n at once**, not just 14.

## 6. Honest status

No new full proof of any even n here. The contributions: (i) the **single unified
target** (measure-0 ⟹ floor-tight) and its **verification across even n=4..14**
(0 counterexamples); (ii) the **collapse of all six lenses onto one worry-object**
(perfect transversal / regular rotational encirclement); (iii) the **HYP-2091
parity carrier** making the target `n+2`-uniform; (iv) the **three-branch attack**
on n=14 with the precise closing lemma. The closable-now cases (4,6,10,12) are
finite checks of 1–2 explicit transversals (literature-implied); 8 and 14 need the
non-unit-pair clock lemma.

**Artifacts:** `04-computation/lrc_unified_target_even_n_s568.py` (+`.out`).
Builds on HYP-2091/S567 (parity), S552/S553/HYP-2084 (spectral gap/transversals),
S559/HYP-2063 (apex/non-unit hole), THM-396/397 (shields), S563-566 (orbit/
Burnside/strong), S557 (pinch), THM-369 (sieve). New: **HYP-2093**.
