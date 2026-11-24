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

- [Hamilton path viewpoints, two-edit squares and exact line partners](../05-knowledge/results/chanpak_tiling_viewpoint_20261123/README.md). On a strong tournament the generated viewpoint group is symmetric except for the triangle and the five-vertex cyclic triangle substitution. Every composable pair of unit count flips at heights at least five shares a fixed-path board. Physical placements in one rigid class can reach exponentially many heights on blue lines and on black lines.

## 5. H (Hamiltonian paths) and the 2-adic tower
- H = 1 + 2^d hypotenuse formula; tilings·|Aut| = H (LEM-003); H-spectrum per class: kps atlas.
- **H by level, exact to n = 11** (THM-4889, [packet](../05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md)): level floors φ_n(k) nondecreasing in k; the global floor φ(k) = 1, 3, 5, 9, 11, 15, 23, 33, 41, 51 (k ≤ 9 exact);
  zone law φ_n(Λ(n−1)+1) > max H(n−1) (values inherited from n−1 live at levels ≤ Λ(n−1), new values above); top-level separation φ_n(Λ(n)) > ψ_n(Λ(n)−1) at n ≤ 10
  (CONFIRMED at order 11 (φ_11(55) = 91,395 > ψ_11(54) = 89,357: every regular 11-tournament beats every non-regular one) at n = 11); the unattained top band = inter-level gap + top-level sparsity; flips raising c3 raise H in 98 % of cases (n ≤ 9); H(T)/max_v H(T−v) ∈ [1, 5.35]. OPEN laws: HYP-9481.
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
3. Level-width sequences: OEIS status (the n = 7 row is not in OEIS, checked 2026-11-06); palindromy break mechanism at n = 6. Exact widths for n = 9, 10, 11 (classes, strong classes and distinct H per level) are in THM-4889's packet (`data/least_order_1x.txt` tails).
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
- **n=9 (THM-931, exhaustive, 191536 classes exact):** one-eyed panel 264 groups/277 pairs;
  extended panel 35 pairs = 27 SUSPENSION-manufactured (cone/cocone commute: the two predicted
  27-families coincide as suspension images) + **8 WILD pairs = CHIRALITY pairs (T, T^op)**:
  non-SC, cpL self-dual, tau_in = tau_out, near-regular scores; cpA-deck and every
  reversal-invariant invariant (c3, c4, c5, |Aut|) tie AUTOMATICALLY; the full-panel deck
  (per-card tau) splits all 8. Blindness bifurcates: LAUNDERED (sink/source) + CHIRAL.
- **Sequences:** distinct cpK 1,2,2,6,11,50 (n=3..8) — NOT in OEIS (checked 2026-07-16),
  submission-ready. Wild-pair counts by n: 0 (n<=8), 8 (n=9).
- **Open:** n=10 wild hunt (is every wild pair chiral? 98M extensions, overnight-feasible or
  chirality-targeted); per-level spectral-tie distribution vs the H-gradient.

## Chan–Pak constructions and exact inverses

- **PROVED + FINITE-EXACT:** [Lonely Runner filler circles](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/FILLER-CIRCLES.md): a literal 61-tick model has exact distance 5/61, and a primitive seven-axis subbank gives a certified 4.338 ppm relative window at fixed Jacobian. [Joint CRT and prime grid](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/FILLER-CRT.md): small mod-77 wheel and six exact prime starts on a 37-vertex layout. The full 14-axis bank forces common profile content and cannot itself give prime starts.

- **PROVED + FINITE-EXACT:** [Third-slot Jacobian fibers](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/THIRD-SLOT-JACOBIAN.md): first `Q3` height read at eight fixed sites; two literal prime sevenfold squares share the same Jacobian, and a verified 14-pair filler bank gives `(t+1)^14` distinct exact starts at fixed curvature. No third-slot contribution to mixed gain survives the census through eleven fixed sites; the primitive frames have a structural same-boundary explanation.

- **PROVED + FINITE-EXACT:** [Minimal scalable frame census](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/TEMPLATE-CENSUS.md): at six fixed sites only A/B are variable-width prime-capable; the complete seven-site scan adds reads `11Q+9D`, `29Q+15D`, `48Q-D`, with pure interaction `4D`, `3D`, `64Q` and examples outside both earlier frames. [Jacobian leading-term vetoes](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/LEADING-TERM-VETO.md): bundled charts impose a rank-one partial degree matrix and common binomial leading edge in the exceptional direction.

- **PROVED + FINITE-EXACT:** [Glued filler capacity](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/MODULE-CAPACITY.md): arbitrary ordinal joins of one-apex fillers plus any fixed bank cannot cover primes in the scalable frames; without extra bank their first scalar and prime escapes are 689 and 1747. [Four bundled Jacobian ports](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/BUNDLED-JACOBIAN.md): constant determinant suffices whenever the controlled edge graphs have matching number at most four, via partial degree eight and a shear for the residual 6/8 pair.

- **PROVED + FINITE-EXACT:** [Unrestricted Jacobian squares](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/UNRESTRICTED-JACOBIAN.md) and [scalable two-apex modules](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/TWO-APEX-SCALABLE.md): every incident pair has mixed height at least twice the outside height in the coherent gauge, the disjoint 31→217 square gives an integral Keller map, and 103→721 has a monotone strong disjoint square. Two arbitrary-module constructions have lower heights `6Q+5D`, `15Q+D` and upper exactly sevenfold; the first jointly controls height and curvature residues and reaches every fixed-modulus scale window of width `O(X^theta)`.

- **PROVED + FINITE-EXACT:** [Complete exact one-apex inverse](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/ONE-APEX-INVERSE.md): all two-contact sevenfold patterns decided with `O(log p)` Boolean candidates; this chart reaches only `O(sqrt(X) log X)` counts up to X. Opposite-side 97→679 and mixed 103→721 first occur at orders 11 and 12; at 103, 227 and 1399 no monotone marked route exists within this chart, at any order.

- **PROVED + FINITE-EXACT:** [Two-edit follow-up](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/README.md#two-edit-follow-up): prime-to-own-sevenfold compression now has cap `B(7p)+2` (at most `+1` on the split multicore branch); a literal 43→13→301 square realizes the two-pad support surplus. The physical wheel accepts any unit cofactor phase and, after one common gcd, its eight distinct form readings are pairwise coprime at sufficient cutoff. Transitive-spine one-apex contact flips have no one-edit sevenfold move and first have a two-edit eligible move 37→259 at host order ten; the same-side exterior chart has an exact but sparse inverse.

- [Lonely Runner wheel and exact two-edit prime ports](../05-knowledge/results/chanpak_lrc_two_edit_ports_20261123/README.md): distance exactly two to own (7p) for every eligible prime through 71, a labelled unit square with nonisomorphic equal-height middles, and simultaneous small-prime avoidance on the literal fork.
- [Planar prime ports and factor-seven seams](../05-knowledge/results/chanpak_planar_prime_ports_20261123/README.md): an exact inverse chart can miss every multiple of seven; smallest literal prime-to-own-sevenfold arc has order nine, and the first eleven eligible primes admit no such edge at any order.
- [Square-sum lifts, signed Collatz and tournament arcs](../05-knowledge/results/square36_signed_collatz_20261123/README.md). Of 118 square-sum pairs among 72 heights, 101 admit one flip. Five labels of an exponent-61 six-cycle are attained.
- [Signed chart and 19th-vertex decoder](../05-knowledge/results/chanpak_signed_square36_20261123/README.md) · [Kaprekar host at base 765764](../05-knowledge/results/chanpak_kaprekar_decimal_20261123/EVEN-BASES-THROUGH-16.md) · [Edit-radius bounds](../05-knowledge/results/chanpak_edit_radius_20261123/EDIT-RADIUS.md).
- [Five Heegner midpoint flips](../05-knowledge/results/chanpak_heegner_selector_20261122/README.md) · [Balanced windows](../05-knowledge/results/chanpak_heegner_absolute_targets_20261123/README.md).

[Order-118 E61 bracket, optimal 7:10 direction and prior constructions](../05-knowledge/results/chanpak_odd_square_spines_20261122/README.md). Exact E61 crossing and Chan–Pak OPEN.

- [Arc-flip geometry at 7 and 21](../05-knowledge/results/chanpak_gap_cubical_geometry_20261122/README.md) proves the core collapse for cutoffs 15–23.
- [2–3–11 brackets and Collatz clocks](../05-knowledge/results/prime_brackets_collatz_2311_20261122/README.md) excludes E61 traces with at most two odd steps at all depths.
- [Two-clock Farey inverse](../05-knowledge/results/chanpak_farey_two_three_clocks_20261122/README.md) covers arbitrary-flank bands and incompatible binary and ternary lifts.
- [Complementary placement](../05-knowledge/results/chanpak_farey_complementary_exact_20261122/README.md) excludes the full adaptive cube-root window on spines 4 and 8.
- [Clocks](../05-knowledge/results/odd_tables_signed_carries_20261111/SYNTHESIS.md); [more](S111-plus.md).
- [Exact U5-free capacity](../05-knowledge/results/ramsey_neighbourhood_lifts_20261030/README.md).
- [THM-5313](../05-knowledge/results/chanpak_switching_divided_powers_20261122/README.md) · [Exponent-31 flip](../05-knowledge/results/chanpak_e31_physical_20261122/README.md) · [E61 inverse and exponents 3–23](../05-knowledge/results/chanpak_e61_fibonacci_target_20261122/README.md).
- THM-1805: triangle Vandermonde. THM-1862/1936: join products for `H`, signed Rédei; additive `c3`.
- THM-2183/2195, 3121–3248: ordinal metrics and resolvents via THM-5313 above. Q4 radius/non-P OPEN.
