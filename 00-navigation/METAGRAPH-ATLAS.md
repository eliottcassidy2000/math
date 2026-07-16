# THE TOURNAMENT METAGRAPH ATLAS
*(klein-2026-07-15-S313 cont.3 — consolidated map of every charted aspect of G_n, each with its
canonical source; update in place when an aspect moves. "Level k" = c3 = # cyclic triangles =
tie-splits from transitive; x = (n³−n)/3 − 8k.)*

## 1. Vertices and quotients
- **|V(G_n)| = A000568**: 1, 1, 2, 4, 12, 56, 456, 6880, … Burnside over all-odd cycle types only
  (even cycles fix nothing — the Feit–Thompson face). Merged |V| = (A000568 + SC)/2; SC = 2, 2, 8,
  12, 88, 176 (n = 3..8). Sources: CLAUDE.md canon; THM-852 (freeness), ORBITS = 2, 3, 22, 101.
- **Even-graph dual E_n**: V = 2, 3, 7, 16, 54 (A002854); χ grows faster; bridge matrix full rank.
  Source: `07-reflections/even-graphs-as-first-class.md`.

## 2. The axis (levels) — COMPLETELY MAPPED
- **Level lattice**: full step-8 progression floor → (n³−n)/3 (THM-866, tie-split walk); k(T) = c3(T),
  each tie-split kills exactly one 3-cycle (HYP-6948); climb increments T_⌊(n−1)/2⌋ (THM-867/868).
- **Level widths (classes per level k = 0..k_max)**, exact n ≤ 7 (atlas census, this session):
  - n=4: 1,2,1 | n=5: 1,3,2,2,3,1 (palindromic!) | n=6: 1,4,4,4,10,8,12,8,5
  - n=7: 1,5,7,8,17,23,39,40,59,61,79,52,47,15,3 (peak 79 at k = 10, NOT at the regular end)
  - Palindromy holds n ≤ 5 only. n = 8 level census: opus-S316/THM-866 .outs.
- **Landau corona** (which lattice shells the tournament world fills): trivial zone quadratic
  ((n²−1)/4 + 1 shells), then a contiguous corona of width (n−6)(n²−1)/24 (odd) / (n−6)(n²−4)/24
  (even) up to the ceiling; onset x* = 2n²+6 / 2n²+n, witness = two-champion duel. **⌊n/2⌋+1
  REFUTED** (n=8 coincidence). Source: THM-869.
- **Lattice geometry**: n=8 sits in E8 (Σ=0 half-coset slice; floor = 70 roots) — THM-868; n=16
  sits in D₁₆⁺ (Milnor's partner, NOT split-frame E8⊕E8); residue laws = discriminant-form
  arithmetic, signature n−1 mod 8 — THM-872. n=24: OPEN (Niemeier selection).

## 3. Symmetry — the odd/solvable world
- **|Aut| always odd ⟹ solvable (Feit–Thompson)**: THM-868 §5. Histograms (exact, this session):
  n=5: {1:7, 3:4, 5:1}; n=6: {1:41, 3:12, 5:2, 9:1}; n=7: {1:399, 3:47, 5:4, 7:1, 9:4, 21:1}
  (21 = QR₇ Paley, the multiplier law Z₇⋊Z₃).
- **Fermat-rung rigidity**: rotational tournaments at n = 3, 5, 17, 257 have Aut = Z_n exactly;
  n=17 has exactly 16 rotational classes (free Z₁₆ torsor, no canonical unit = the sedenion
  degeneration); n=9 anti-rigid: C₉{1,3,4,7} = C₃[C₃], Aut = Z₃≀Z₃ (81 = n²). Source: THM-871.
- **A₅ cannot enter by symmetry**; it enters (i) by lattice (E8 = McKay(2I), THM-868), (ii) as
  monodromy of the 5 tetrahedral blocks of the icosian frame (2I/2T, kernel ±1 — verified), (iii) as
  base monodromy on labelled boxes: 24 orbits on the n=5 box (2×12 + 8×20 + 14×60, no fixed
  tournament), 560 on the 6-axis box, 4,495,872 on the n=8 (5+3) box. Deck groups of all class
  fibrations are solvable — any unsolvable monodromy must come from base loops, never decks.
- **Kakeya**: needle directions of A₅ = the 31 icosahedral axes; K(A₅) = 15 = K_odd(A₅) (the
  even axes are free); K(2I) = 30 (±-descent). Source: THM-870.

## 4. Edges and layers (the flip structure)
- Waggly = all Hamming layers d = 1..m; wiggly = d = 1; blue/black = d = m (CLAUDE.md canon).
- H-gradient: mostly increasing, level edges 0,0,1,15,136 (n=3..7), H-decreasing exist n ≥ 7.
- Mod-16 selection on the d = m flip layer: THM-790 (a FLIP-layer rule, not a level rule).
- Spine/ribs/sea (SC-SC / SC-NS / NS-NS): THM-A/B/C; sea dominates (96% at n=8).

## 5. H (Hamiltonian paths) and the 2-adic tower
- H = 1 + 2^d hypotenuse formula; tilings·|Aut| = H (LEM-003); H-spectrum per class: kps atlas.
- Digit tower: digit 0 constant (Rédei); digit 1 = c3-governed to n = 6, dies n = 7; the locker
  tournament refutation H(D_11) ≡ 3 (THM-865). "The OCF is the Toda tower" (opus-S316).

## 6. Truncation grammar (the shadow calculus this all sits in)
- Vandermonde tail: polygonal = rank-1 truncation of Pascal; Moser/A000127; D-sequence; the (r,g)
  master cancellation + missing-region law (deficit exactly 1 at d₀ = (g+2)(r+2)−1); Brown
  completeness everywhere. Sources: HYP-6911/6912/6946/6947, drafts S313/S313c2.
- Corona = the grammar's Landau instance (THM-868 §3, THM-869); Farey-14/golden and the OCF tower
  are the other instances (S315-S317).

## 7. Named gaps (the unmapped aspects)
1. Corona onset minimality for ALL n (convexity argument to formalize — THM-869).
2. n = 24 Niemeier selection (THM-872 §4).
3. Level-width sequences: OEIS status; palindromy break mechanism at n = 6.
4. Exact metagraph antichain width at n ≥ 7 (old formula fails; distinct from level width).
5. A₅ base-monodromy on a NATURAL class fibration (the THM-848 split (H,K) fibre is the candidate).
6. K(G) for the other odd-relevant groups (Z_q⋊Z₃ multiplier groups; the 16 sedenion classes as a
   "direction set"?).
7. The 16 sedenion classes: full invariant table (c4 partially separates: 1105..1428).

## 8. Spectral columns (klein-S315/S316: census, reciprocity, coning)

**What the spectra see (exact censuses):**

| n | classes | cpA-tie groups (classes in ties) | distinct cpA | distinct cpK |
|---|---------|----------------------------------|--------------|--------------|
| 6 | 56      | (census S315c2)                  | —            | 6            |
| 7 | 456     | 116 (404 = 89%)                  | 152          | 11           |
| 8 | 6880    | 1460 (6817 = 99.1%)              | 1523         | 50           |

- **cpK is a FUNCTION of cpA** (THM-924: cpK(y) = 2^{n−1}[cpA((y−1)/2) + (−1)^n cpA((−y−1)/2)]);
  walk moments and skew d-moments are cpA-determined too (walk reciprocity det(zI−A+J) = (−1)^n cpA(−1−z)).
  Distinct-cpK sequence **1, 2, 2, 6, 11, 50** (n=3..8) = the symmetrization collapse about Re z = −1/2. New sequence.
- **Splitter ranking at n=7** (S315c2): τ_in 116/116 > cpL 111 > scores 85 > H 47 > |Aut| 1 > cpK 0 (now a theorem).
- **Invisibility = the cone stratum** (THM-918 construction + THM-925 completeness): at n=8, the
  one-eyed panel (cpA,cpL,H,τ_in) has EXACTLY 27 invisible pairs = the manufactured cones; the
  two-eyed panel (+τ_out) has EXACTLY 4 = the double-blind four (H=23,29,31,43; all with sink AND
  source; τ_in = τ_out on them). Zero wild (non-cone) invisibles at n=8. Tower: invisible pairs at
  every n ≥ 8 (cone transform laws; sink = black hole crushes τ_in to det(L+I), source = mirror scales by n).
- **Open:** n=9 completeness (first wild pair?); OEIS submission for 1,2,2,6,11,50; per-level
  spectral-tie distribution vs the H-gradient.
