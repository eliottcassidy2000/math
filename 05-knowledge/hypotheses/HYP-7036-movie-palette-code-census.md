# HYP-7036 — The movie palette code [n,k,d] census: conjecture refuted-as-stated, replaced by the rate law

**Status:** RESOLVED (death-star-2026-07-16-S22; owner directive: compute the census).
Script/results: `movie_palette_code_census_deathstar_S22.py` / `.out` (17 cores, 6 species).

## Structure theorem (proved, one paragraph)

Each wall is crossed exactly once per period ⟹ each wall lies on exactly one edge of the
movie multigraph ⟹ the palette map (cycle space → F₂^walls) is INJECTIVE and the palette
code **is the cycle code of the movie multigraph** — a graphic code with
**[n, k, d] = [#events, #events − #states + 1, girth]** (no self-loops; parallel edges
give d = 2).

## The census (17 cores)

| species | example | [n, k, d] | rate k/n | κ |
|---|---|---|---|---|
| exact dilate | 5·[2,3,5,7,11,13] | [1260, 1009, 2] | **0.801** | 3 |
| consecutive @ scale | [50..55] | [2114, 1002, 2] | **0.474** | 4 |
| consecutive @ scale | [30..35] | [1274, 246, 2] | 0.193 | 4 |
| far bank | [0..5,50] | [378, 43, 2] | 0.114 | 3 |
| random | #0 | [5390, 295, 2] | 0.055 | 3 |
| planted gcd-sub | 24+48=72 | [3724, 127, 10] | 0.034 | 3 |
| consecutive small | [20..25] | [840, 22, 37] | 0.026 | 4 |
| AP | [0..6] | [84, 1, 84] | 0.012 | 3 |
| planted sporadic | 23+48=71 | [4032, 15, 337] | 0.004 | 3 |
| generic incoherent | c = 10/20/30 | [n, 1, n] | ~0.001 | 3–5 |

## Findings

1. **The S21 conjecture (d = coherence meter) is REFUTED as stated:** κ is FLAT (3–5)
   across the whole taxonomy — thirteen-element-scale integers always carry tiny relations
   (Probe B's lesson, recurring at 6 elements) — while d spans 2..3458. No monotone d ↔ κ.
2. **The rate k/n is the coding-theory face of coherence** (dilate 0.80 > consecutive@scale
   0.47/0.19 > far 0.11 > random 0.06 > incoherent 0.00x), with the sharp refinement:
   **only SUBGROUP-TYPE relations (gcd/dilate substructure) drive recurrence** — the
   gcd-24 planted core (rate 0.034, d = 10) vs the isolated-relation core 23+48=71 (rate
   0.004, d = 337). Sporadic additive relations do NOT force returns: the S21 transference
   duality is one-directional at finite scale (Fourier-side deviations ≠ state-side
   recurrence unless the relation generates a rational sublattice).
3. **d is a recurrence-ONSET indicator**, with three regimes: d = n (k = 1, "Hamiltonian
   movies", zero proper recurrence — the [n,1,n] full-cycle code: AP, GW-ish, small/
   incoherent cores); a middle band (d = 10..434: few long fundamental loops); d = 2
   (parallel edges) once recurrence is dense. Fragile as a distance, informative as a
   phase marker.
4. **The Moore-bound connection lands on the dichotomy**: incoherent movies sit at the
   cycle-space critical point (avg degree barely > 2 — outside the Alon–Hoory–Linial/
   2607.14068 density regime); coherent movies enter the dense regime where girth-forcing
   applies (and indeed collapse to d = 2). The paper's machinery applies exactly to the
   coherent species — one more face of the seam dichotomy. (Caveat: the .out's "moore~"
   display column uses a clamped log and is invalid near avg-degree 2; excluded from
   conclusions.)

## Named next steps

(i) The rate law deserves a formula: k/n for dilates is computable exactly (period
collapse: k = n − n/c·... the c-fold cover structure) — derive and test; (ii) the
subgroup-vs-sporadic relation dichotomy should be pushed back through THM-890: which
relation TYPES carry the Fourier deviation mass (conjecture: sporadic relations carry
Fourier mass but no recurrence — testable by comparing THM-890 weights with return counts);
(iii) the middle-band d (few-long-loops movies) as a certificate: d large ⟹ no short
expressive cycles ⟹ palette vectors well-separated ⟹ a decidable incoherence certificate
cheaper than Fourier scans.

-> HYP-7027, THM-889/890, T1540, 2607.14068 reflection; death-star-S22.
