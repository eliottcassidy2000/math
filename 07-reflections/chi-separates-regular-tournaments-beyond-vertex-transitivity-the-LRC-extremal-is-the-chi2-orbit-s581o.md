---
source: oracle-2026-06-03-S581o
status: computation + answer (chi separates regular tournaments beyond VT and beyond cyclicity; the LRC tight orbit is the chi=2 rotational R_m, not Paley)
tags:
  - lonely-runner
  - regular-tournament
  - dichromatic-number
  - chi
  - paley
  - vertex-transitive
  - extremal-family
---

# χ Separates Regular Tournaments Beyond Vertex-Transitivity — and the LRC Extremal Is the χ=2 Orbit, Not Paley

User: among the maximally-cyclic (regular) tournaments, does χ add anything beyond
vertex-transitivity — i.e. is there a tight LRC config that is regular but not the Paley/AP
orbit, and does its χ differ?  Answer (χ = the **dichromatic number** = min #colors so each
class is acyclic/transitive; a tournament subset is acyclic iff it has no 3-cycle):

## 1. "Paley/AP" is *two* orbits, and only χ tells them apart

A regular tournament has every score `(m−1)/2` — "maximally cyclic." Two facts force a finer
invariant:
- **All regular tournaments on `m` vertices have the SAME number of 3-cycles** (it is a
  function of the score sequence, which is constant for regular ones): m=7 → 14, m=9 → 30,
  m=11 → 55. So *cyclicity cannot distinguish regular tournaments.*
- The symmetric ones are all **vertex-transitive** (VT), so VT cannot distinguish the
  symmetric ones either.

Yet the **rotational `R_m`** (connection set `{1,…,(m−1)/2}` = the AP at its n-clock tight
time) and the **Paley `QR_m`** (`m ≡ 3 mod 4`) are **non-isomorphic** — and **χ separates
them** (`lrc_regular_tournaments_chi_vs_vt_s581.py`, `..._s581b.py`):

| m | χ(R_m, AP-rotational) | #3-cyc | aut | χ(Paley QR_m) | aut | iso? |
|---|---|---|---|---|---|---|
| 5 | 2 | 5  | 5  | (n/a, 5≡1) | – | – |
| 7 | **2** | 14 | 7  | **3** | 21 | **no** |
| 9 | 2 | 30 | 9  | 3 (QR-circ) | 81 | no |
| 11| **2** | 55 | 11 | **4** | – | no |

> **χ(R_m) = 2 for every m, while Paley has χ = 3 (m=7), 4 (m=11).** `R_7` and `QR_7` are both
> VT, both regular, both 14 three-cycles, both self-converse — *indistinguishable by
> vertex-transitivity or cyclicity* — but `χ` is `2` vs `3`. So **χ adds strictly beyond
> vertex-transitivity (and beyond the 3-cycle count): it is the right finer invariant on the
> maximally-cyclic family.**

And there is even more structure χ must carry: m=7 has **3** regular tournaments (only 2 are
circulant), m=9 has 15 — so regular orbits beyond Paley/AP exist from m=7 on.

## 2. The LRC extremal is the χ=2 (least-cyclic) orbit, *not* Paley

`R_m` is the **dichromatic-number-2** regular tournament — the regular tournament *closest to
transitive* (it splits into just two transitive parts). Paley (χ=3) is genuinely *more*
robustly cyclic. The LRC tight family (S576o) realizes **only the χ=2 orbit**
(`..._s581b.py`, runner tournament at the tight time, ties resolved by ±ε):

- **n=6 AP, n=8 AP:** the `−ε` resolution is **exactly `R_m`** (regular, `==R_m`, χ=2).
- **n=8 sporadics `{1,2,3,4,5,7,12}`, `{1,4,5,6,7,11,13}`:** χ=2, near-`R_m` (tie-degenerate,
  not exactly regular) — **but never Paley.**
- **No tight config through n=8 realizes the χ=3 (Paley) orbit or any other regular orbit.**

> **The LRC tight regular orbit is the MINIMALLY-cyclic one (`R_m`, χ=2 = the AP) — not the
> maximally-symmetric Paley (χ=3).** Being χ=2 ("barely cyclic," nearest to transitive) is
> exactly why the AP sits at the loneliness boundary `M=1/n`: it is the regular tournament
> least able to spread the runners. Paley, more cyclic, corresponds to *looser* speed sets.
> **χ=2 is a candidate characterization of the LRC tight regular orbit.**

## Direct answer to the question
- **Does χ add beyond vertex-transitivity?** **Yes, decisively** — `R_m` and Paley are both VT
  with equal 3-cycle counts, and only χ separates them (2 vs 3,4).
- **Is there a tight config that is regular but not the AP/Paley orbit?** Among the tight
  family checked (n≤8): **no** — every tight config is the χ=2 `R_m` orbit (the AP); the
  χ=3 Paley orbit and the other regular orbits are *not* realized by tight configs.
- **Does its χ differ?** The tight orbit is uniformly **χ=2**; the regular orbits the LRC
  extremal does NOT touch are exactly the **χ≥3** ones. So χ both (i) refines VT on the
  regular family and (ii) appears to *characterize* the tight orbit (χ=2).

## Honest limits / next
- The third (non-circulant) regular tournament on 7, and the 12 non-circulant ones on 9, were
  not χ-computed (full regular enumeration is slow) — so "`R_m` is the *unique* χ=2 regular
  tournament" is **conjectured, not proven.** Verifying it would make "χ=2 ⟺ AP-rotational
  orbit" a clean characterization.
- The tight family is checked only through n=8 (bounded speeds). The conjecture to test next:
  **every LRC tight config (all n) realizes a χ=2 (near-`R_m`) tournament; no tight config is
  Paley-type.** Equivalently: regular tournaments with χ≥3 give *strictly lonely* speed sets.
- Connects to: HYP-2091 (the rotational `R_m` is the clean-polygon extremal), S576o/S577o (the
  worry-set / tie-wall), and the Barajas-Serra circular-chromatic LRC method (the χ here is
  the tournament-native dichromatic analogue).

## Artifacts
```
04-computation/lrc_regular_tournaments_chi_vs_vt_s581.py          (+.out)
04-computation/lrc_tight_config_regular_orbit_chi_s581b.py        (+.out)
```
Related: HYP-2091, S576o (even-ladder worry), S577o (tie-wall), vt-trienerment-polygon-rigidity-s589.
