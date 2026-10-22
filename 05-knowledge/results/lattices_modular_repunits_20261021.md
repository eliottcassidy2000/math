# Lattices, modular forms and repunits: tournament scores in II_{25,1}, Leech deep-hole theta series, three PDF referee reports (opus-2026-10-21-S52)

**Owner prompt (2026-10-21).** The prompt listed these topics:
- 1729, the 24 transverse dimensions of the bosonic string, q^(1/24) and algebraic geometry;
- Smale's h-cobordism theorem;
- the Heegner numbers 1, 2, 3, 7, 11, 19, 43, 67, 163;
- why 31, 63, 127 are rare, and all prior repunit work;
- 25 = 24 + 1, Lorentzian versus Euclidean, the even unimodular lattice II_{25,1};
- the 24 Niemeier lattices, deep holes and sphere packing, and "the same in 8 dimensions";
- {7,21} against {8,24};
- E8: 248 = 8 + 3·80, and its Langlands self-duality;
- Δ = (E4³ − E6²)/1728, E8 = E4², E10 = E4E6, E14 = E4²E6.

Three PDFs by the "SRFP311T1 Collaboration" came with it, treated as data: `nebe_voronoi_rigidity.pdf`,
`mclaughlin_slice_arithmetic.pdf` and `e8.pdf`.

**Concurrent work.**
- Session27 received the same prompt with four other PDFs and wrote the main synthesis:
  [`session27_repunits_lattices_20261022/SYNTHESIS.md`](session27_repunits_lattices_20261022/SYNTHESIS.md). It covers:
  - the repunit catalogue and the measure-aware rarity law;
  - θ/η characters and 1729 mechanisms;
  - II_{9,1} null quotients and the E8 deep hole;
  - the {7,21} → {8,24} type obstruction;
  - the Smale/Kervaire–Milnor source check;
  - an independent check of the norm identity in §1.1.
- mac-mini owns the class-number/Heegner, circulant, QR-code → E8/Golay/Leech and {7,21} ↔ {8,24} lanes
  (THM-4751 / HYP-9362; MSG-1913, MSG-1914).

This note records what S52 adds. Status tags used throughout: PROVED, KNOWN (cited), FINITE-EXACT, NUMERICAL,
NUMEROLOGY, REFUTED.

**Scripts.** All are in `04-computation/experiments/lattices_repunits_20261021/` and each prints ALL CHECKS PASSED:
- `pdf_referee.py`, `lorentzian_scores.py`, `lattice_cusps.py`, `numerology_tests.py`;
- `bitdepth_trace_test.py` (cypari, about 25 minutes);
- `hyp9363_holes.py` and `hyp9363_deep_hole_theta.py` (cypari, numba);
- `hyp9363_weil_attempt.py` (a failed sufficient condition, recorded).

The independent audit is in `audit_S52/` (REPORT.md plus five scripts); MISTAKE-661 records the errors it found.

## 1. Tournament scores in II_{25,1} (THM-4752)

For a 25-vertex tournament T with scores s, put v_T = (s | 70) ∈ II_{25,1}. Write ρ = (0, 1, …, 24 | 70) for Conway's
Weyl vector and ν = (½²⁵ | 5/2).

1. **Norm and the Leech cusps (PROVED).** v_T² = −2c3(T). So v_T is null exactly for the transitive tournaments, the
   25! permutations of ρ, i.e. the Leech cusps.
   - 0² + … + 24² = 70² is Watson's cannonball, and n = 25 is the only order with a nonzero null transitive vector.
   - Announced in MSG-1913 and checked independently by Session27, which adds a hostile example: equal v_T with
     H = 657 and 643.
   - Earlier repo sightings: `ogg_triangular_chowla_20260927.md` had named (0, …, 24; 70) "the Leech vector in
     II_(25,1)" as NUMEROLOGY. The variance form Σ(s − μ)² = n(n² − 1)/12 − 2c3 is THM-785 / THM-4739.
2. **Backward length is the Weyl-vector pairing (PROVED).** For every ranking π, BL_T(π) = −v_T·ρ_π.
   - The score rankings minimise it (rearrangement; Kano–Sakamoto, Networks 13 (1983) 143–151).
   - The minimum is U(T) = Σ Landau slacks = c3 + ½|s↑ − ι|²; this closed form is elementary and possibly folklore
     (Monsuur, EJOR 2005, unread).
   - U is the horospherical height of v_T over the Leech cusps. The minimising cusp is also the nearest, but at
     squared distance 2(U − c3), not U.
3. **The A24 cusp (lattice side KNOWN, tournament reading new).** The 24 adjacent transpositions e_i − e_{i+1} of score
   ranks and r* = (−½, ½²³, 3/2 | 5/2) are Leech roots. They form an Ã24 cycle and sum to ν. So:
   - ν·ρ = −25 = −h(A24), and ν^⊥/ν = N(A24);
   - these 25 are exactly the Leech roots orthogonal to ν (a two-line lemma from r_λ·r_μ ≤ 0).

   This is the entry 2ν (height 25) of Borcherds' table of norm-0 vectors
   (https://math.berkeley.edu/~reb/lattices/norm0; thesis, arXiv math/9911195, Lemmas 2.2.4, 2.4.2, 2.5.1), read in
   tournament coordinates. Only ν, of the two S₂₅-fixed null rays, lies in the closed Leech chamber; the other is
   (−1²⁵ | 5).
4. **Every tournament has an A24 class (PROVED).** v_T·ν = −25 for all T.
   - The class of T in N(A24) is its centred score vector s − 12 ∈ A24, with |s − 12|² = 1300 − 2c3.
   - The transitive class (12, …, −12) is the Weyl vector ρ_{A24}, and 1300 = 2h(h + 1).
   - A 25-vertex digraph with m arcs and 5 | m lands in the glue class [−m mod 25].
   - At n = 25, centred scores in A_{n−1} are THM-872's odd-n case. The new step is the canonical ambient Niemeier
     lattice ν^⊥/ν. **The n = 24 question of THM-869/872 stays open.**
5. **The light cone (PROVED).** By Landau's theorem, as in THM-1825-q-redei, the score vectors are the lattice points of
   conv(S₂₅ρ).
   - The 25! vertices are Leech cusps, and every other score vector is timelike.
   - S₂₅ = W(A24) acts by relabelling.
   - Score vectors are the weights of V_δ, since Σ_T x^score = s_δ (camion_busch §6.1; Gessel 1979).

**What the lift loses.** H, all cycles beyond triangles, and the arcs. The dictionary says nothing about the gaps 7 and
21 or about repunits.

## 2. The Leech lattice seen from a deep hole (HYP-9363, FINITE-EXACT; new conjecture)

The audit counted Leech roots by their level r·ν = −k and found 25, then 5·p(k − 1) for k ≤ 22 (p is the partition
function). These are the counts N_k of Leech points at squared distance 2 + 2k/25 from the A24 hole. I conjectured an
η-quotient before computing further:

    θ_{Λ+c}(τ) = 5 η(τ)²⁵/η(τ/25) + 25 Δ(τ)  =  θ_{N(A24)+ρ/25}(τ) + 25 Δ(τ).

The second form follows from two facts: η(τ)²⁵/η(τ/25) is the 25-core generating function (Garvan–Kim–Stanton), and
the five glue cosets of N(A24) shifted by ρ/25 have equal θ.

**A24 check.** Exact enumeration of 20,219,385 Leech points (Leech lattice built from scratch as {ρ, ρ′}^⊥) confirms
it through q^{1+60/25}, i.e. 61 coefficients.
- This includes the predicted first departure from 5p(k − 1): N₂₅ = 7275, not 7875, coming from Δ's −24q².
- The run 5p(k − 1) for k ≤ 24 holds because every partition of fewer than 25 is a 25-core.

**Generalisation.** For the pure A-types A_{h−1}^k (k(h − 1) = 24):

    θ_{Λ+c} = h^{k/2} (η(τ)^h/η(τ/h))^k + (24 + k) Δ.

The coefficient is consistent because 24h/(h − 1) = 24 + k. Tested on two more holes, both exact:
- **A1²⁴ (c = v/2, |v|² = 8):** 48, 4096, 97152, 1130496, 8400704, 45785088 through q^{3.5}. The η-term counts sums of
  24 triangular numbers.
- **A2¹² (c the centroid of a √6-triangle with 36 nearest points):** through q^{1+8/3}.

**General conjecture.** θ_{Λ+c} − θ_{N+ρ_N/h} = (24 + k)Δ for every deep hole. **No proof yet.**
- Weil-invariance of the coset indicator on M = {ρ, ν}^⊥, with D = (Z/25)², FAILS as a sufficient condition
  (`hyp9363_weil_attempt.py`).
- The identity needs relations θ_{M+(−1,v)} = θ_{M+(v,−1)} among different cosets.
- Literature searches (deep-hole theta series, η-quotients, t-cores) found nothing.

## 3. Cusps in ranks 8, 16, 24 (KNOWN, reproduced FINITE-EXACT; `lattice_cusps.py`)

Vinberg's algorithm, run from scratch with the finite-volume criterion, gives:
- **II_{9,1}:** 10 roots (E10) and one cusp Ẽ8, z = (1, 0⁸ | 1); ρ² = −1240 and ρ·z = −30 = −h(E8).
- **II_{17,1}:** 19 roots and two cusps, Ẽ8 + Ẽ8 (z = (1⁹, 0⁸ | 3)) and D̃16 (z = (1, 0¹⁶ | 1)). These are the
  E8×E8 and Spin(32)/Z₂ heterotic lattices. THM-480 and THM-872 already met D16⁺ rather than E8².

Siegel's mass formula checks exactly in ranks 8, 16 and 24.
- Rank 16: 691/277667181515243520000; the 691 is B₁₂'s.
- Rank 24: the 24 Niemeier lattices with |Aut N| = |W(R)|·|G₁|·|G₂| from SPLAG Table 16.1.
  - The audit showed every |G_i| entry is pinned by the identity.
  - The check caught my recall error: G₂(A₂¹²) is M₁₂, not M₁₁.
- θ_N = E₁₂ + (24h − 65520/691)Δ.

**E8 holes and Milnor's drums.**
- Shallow hole (5/6, 1/6⁷): squared distance 8/9 to exactly 9 points (a simplex).
- Deep hole e₁: distance 1 to 16 points (a cross-polytope; Session27).
- E8² and D16⁺ have the same θ (480 roots, 61920 norm-4 vectors) but different root systems.

**"The same in 8 dimensions".** In II_{9,1}, E10's unique affine subdiagram Ẽ8 is the only cusp, i.e. E8. In
II_{17,1} there are two cusps; in II_{25,1}, 24 (23 deep holes plus Leech). The E8 holes correspond to the
maximal-rank subsystems D8 and A8, with E8 = D8⁺ = A8⁺.

## 4. Referee reports on the three supplied PDFs

| PDF | Claim | Verdict |
|---|---|---|
| nebe_voronoi_rigidity | Extremal weight-24 theta E4⁶ − 1440E4³Δ + 125280Δ²; shells 52,416,000 / 39,007,332,000 / 6,609,020,221,440 (the first is the HYP-2433 shell); density (3π/2)²⁴/24! = 2.31783e−8 | **Correct** (exact q-series; norm-12 shell 437,824,977,408,000 in the audit) |
| | Vectors of norm < 2μ are strictly Voronoi-relevant, so **at least** 6,648,079,969,440 facets (Theorem 4.1, a lower bound; norm-12 relevance left open) | **Correct as a lower bound.** The lemma is folklore: \|v+2u\|² − \|v\|² ≥ 2(2μ − \|v\|²) > 0 |
| | x + y (x ⊥ y minimal) of norm 12 is not strictly relevant | **Correct**: x + y ≡ x − y mod 2Λ |
| | Cohn–Elkies bound Δ_CE(48) ≈ 3.01e−8 and a "23 % gap" | **Inconsistent with published LP numerics.** Interpolating Table 1 of Afkhami-Jeddi–Cohn–Hartman–de Laat–Tajdini (arXiv:2006.02560) gives ≈ 4.5e−7, about 20× the extremal density (NUMERICAL; the same interpolation reproduces Leech at d = 24 to 2e−4) |
| mclaughlin_slice_arithmetic | Even rank-2 forms of det 15 with min ≥ 4 are [4, ±1, 4]; this is the norm form of p₂ = (2, ω) in Q(√−15), h = 2 | **Correct, elementary** |
| | Stab(L15) = (McL:2)×2; (L15⊕U)⊕K_McL ⊂ II_{25,1} of index 15 (Nikulin) | **Correct given Conway's 2-2-3 theorem**; the audit's constructed Leech gives 552 = 2·276 |
| | Prop. 6.1: D = −15 is the only fundamental discriminant with a reduced non-principal form a = c = 2 | **True.** The proof skips b ∈ {0, ±2}, but those give the non-fundamental D = −16, −12 |
| | (3,1) spacetime and vacuum selection | The paper itself flags these as open |
| e8 | ζ_E8(10) = 240·2⁻⁵ζ(5)ζ(2) = (5/4)π²ζ(5) | **Correct, classical** |
| | "Bit Depth": y² = x³ + Ax + B mod 2⁵²¹ − 1 with A = int(−ζ(5)·10¹⁶⁰), B = int(qπ²·10¹⁶⁰) has a_p = 5 at q = 5/4 (the paper's own trace line is commented out) | **REFUTED.** PARI SEA: every a_p for q ∈ {5/4, 1, 3/2, 2, 240/96} at 15 and 200 digits is a generic 259–262-bit trace. The audit also proves a_p ≠ ±5 without SEA: [p + 1 ∓ 5]P ≠ O, and one order is certified |
| | Proton/electron ratio 1836.15 "emerges" | Unsupported; NUMEROLOGY |
| | Promised "odd zeta values for Λ24" | See the Leech identities below |

**The Leech analogue** (NUMERICAL, from the classical splitting
ζ_Λ(s) = (65520/691)·2^{−s/2}·[ζ(s/2)ζ(s/2 − 11) − L(Δ, s/2)]):
- **s = 10 (analytic continuation).** ζ(−6) = 0 removes the Eisenstein part, so
  ζ_Λ(10) = −(65520/691)·2⁻⁵·L(Δ,5) = −1.97552397007037, with L(Δ,5) = 0.666709188434003643.
  The value is purely cuspidal: no odd ζ appears.
- **s = 26 (convergent).** ζ_Λ(26) = (65520/691)·2⁻¹³·[ζ(13)ζ(2) − L(Δ,13)] = 0.00749948328893. The non-critical
  value L(Δ,13) = 0.99720986806 is present.
- Both values were checked by independent routes in the audit (θ-continuation and E4³ − 720Δ, agreeing to 1e−40).

## 5. Modular constants, im J and exotic spheres (KNOWN; `numerology_tests.py`)

- **The 24 of q^{1/24}** is |Mp₂(Z)^ab| = 24, since η's multiplier e^{2πi/24} generates the character group. Compare
  |SL₂(Z)^ab| = 12 = |Pic(M₁,₁)| (Mumford).
- **1/η²⁴** counts the 24 transverse oscillators (26 = 24 + 2).
- **744 = 720 + 24 = 3·248 = dim E8³**, the roots plus the rank in θ_{E8³}/η²⁴.
- **Adams.** |im J| in stem 4s − 1 is denom(B₂ₛ/4s) = 24, 240, 504, 480, 264, 65520. The q-coefficients of E₂ₛ are
  their reciprocals up to sign (the same Bernoulli number in two normalisations).
  - So E4 = θ_E8 carries 240 as both the E8 roots (the kissing configuration) and π₇ˢ = Z/240.
  - The 691 in θ_Leech = E₁₂ − (65520/691)Δ is num(B₁₂/24), which divides |bP₂₄| = 2¹⁰·2047·691 (exotic 23-spheres).
- **Kervaire–Milnor.** |bP₄ₘ| = 2^{2m−2}(2^{2m−1} − 1)·num(4|B₂ₘ|/m) gives 28, 992, 8128, 261632, 1448424448,
  67100672. This matches MISTAKE-563's corrections, Milnor 2011 Theorem 5 and Session27's source check.
  - The Mersenne factors 2^{2m−1} − 1 are the odd-exponent repunits; 63 is never one.
  - The Euclid factor is perfect at m = 2, 3, 4, 7, 9, 10.
  - 3 ∤ 2^{2m−1} − 1, the Φ₂(2) = 3 arithmetic of THM-4749.
- **Smale.** Kervaire–Milnor compute Θₙ, the group of h-cobordism classes of homotopy spheres. Smale's h-cobordism
  theorem identifies Θₙ with the diffeomorphism classes Sₙ for n ≥ 5 (Milnor 2011, footnote 2). The E8 plumbing (σ = 8)
  bounds the generator of bP₈ = Z/28. Reading "smale h curve order" as h-cobordism is an interpretation, also
  Session27's.
- **Link to repunit H-values: ANALOGY only** (binary doubling). No map exists. Atiyah's 1987 paper on log η and framings
  is relevant context but was not re-read here.

## 6. Numerology tests (`numerology_tests.py`)

- **The Paley cube pattern (NUMEROLOGY).** H(P_q)/|Aut| = 1, 9, 1729 at q = 3, 7, 11 (tabled in THM-4747), so value −
  1 = 8m³. It fails at q = 19, 23, 27, 31, the same place as the prior split-prime reading
  (`07-reflections/the-1729-resonance-is-isolated…`, Session27 §6).
- **HYP-9273** (S33/S34; pre-registered q = 31 pass) re-checks at all seven q, including q = 27 with |Aut| = 3|Aff|.
- **Heegner (KNOWN).** The j-values are −15³, …, −640320³, and e^{π√163} − 640320³ − 744 = −7.4993e−13. The lane is
  mac-mini's; see also THM-4666 and S17.
- **Perfect codes (KNOWN).** Binary 1-perfect codes exist exactly at lengths 2^r − 1; Golay satisfies
  1 + 23 + 253 + 1771 = 2¹¹.
- **Projective-plane sizes q² + q + 1 = 111_q against the strong H-spectrum below 2100** (723,544 classes, reproducing
  S51).
  - Only 7 = |PG(2,2)| and 21 = |PG(2,4)| are missing; these are the gaps (THM-1370).
  - All other sizes up to 1893 are attained.
  - The factor-2 band test has little power: it contains 98.6 % of all odd d. The audit's regression on log d and
    [3|d], [9|d], [5|d] finds no projective-plane outlier; 133 has z = +2.74, an excess. The 0.67 at 31 is the
    3-divisibility effect.
- **A lawful design-side map {7,21} → {8,24} (KNOWN design theory; recorded by mac-mini as HYP-9362 (d), NUMEROLOGY
  for H; my pointer was MSG-1914; the chains themselves go back to S80, S88, S100).**
  - An extension of a projective plane of order q is a 3-(q² + q + 2, q + 2, 1) design. Its block count
    (q² + q + 2)(q² + q + 1)/(q + 2) forces (q + 2) | 12, i.e. q ∈ {2, 4, 10}, and q = 10 does not exist (Lam 1989).
  - So {7, 21} = {111₂, 111₄} are the point sets of the only extendable projective planes.
  - They extend to S(3,4,8) (extended Hamming → E8) and to S(5,8,24) via S(3,6,22) and S(4,7,23) (Golay → Leech,
    using the Z/4 Hensel lift; Bonnecaze–Solé–Calderbank 1995).
  - As a statement about H-values this is NUMEROLOGY.

## 7. Why 31, 63, 127 are rare: the merged answer

Sources: S51 (THM-4748, THM-4749, HYP-9359), Session24, Session27. S52 adds only §6.
1. **Construction complexity.** 2^k + 1 needs one reversed arc. 2^k − 1 needs three for k ≥ 5; 15 = 2⁴ − 1 needs two.
   - The Fibonacci fan realises 2^k − 1 for every k ≥ 4.
   - The fan and its converse are the only realisations with at most three back arcs among strong classes with
     H < 2100, for k = 6, 9, 11 (FINITE-EXACT).
2. **Degree law (Session24, THM-4749).** N(n,h) is eventually, for n ≥ ⌊3 log₅ h⌋, a polynomial of degree w(h), the
   longest factorisation into strong values.
   - 31 and 127 are prime, and 63 factors only through the gaps 7 and 21, so w = 1 for all three; 33, 65, 129 have
     w ≥ 2.
   - N₆₃(n) = 2(n − 7).
   - The comparison can go the other way: 255 is commoner than 257.
   - Session27's measure caveat: this concerns N(n,h) as n → ∞, not strong-class counts such as S(31)/S(33).
3. **Arithmetic support (HYP-9359, NUMERICAL).**
   - Strong-class counts favour 3 | d and 9 | d, almost all through one-vertex gluings and order-3 automorphisms.
   - For odd k, 3 divides 2^k + 1 and never 2^k − 1.
4. **Open.** S₇(63) = 0, and the H ≡ 7 (mod 8) depletion in the second octave.

## 8. Prior repo work merged

Repunits, Mersenne and 1729 are covered by Session27's catalogue (113,808 paths) and S51. For lattices, modular forms
and topology, an S52 catalogue sweep (3102 files) found the following.

**Score lattices.**
- THM-868 (E8 at n = 8, A8 at n = 9), THM-869 / THM-872 (D16⁺, D24⁺, the open n = 24 question).
- THM-785 / THM-4739 (score variance), THM-1825-q-redei (Landau = permutohedron), THM-207, camion_busch §6.1.
- Session27's LORENTZIAN_SCORE_BRIDGE and SCORE_ROOT_COLLISION.

**Codes.**
- THM-481: Paley gauge codes ê8, Golay g24, extremal [48,24,12]; q = 71 fails.
- THM-480: the skew tower gives d16⁺, not e8².
- THM-489: P24 ↦ 16Δ; plain Construction A of Golay gives N(A1²⁴).
- The S43 P_23 → Golay → M24 → Leech chain; HYP-2433 (extremal first shells).

**Class numbers.** THM-4666 (class numbers as Paley ranking excesses), THM-4558, and S17 Heegner-7.

**1729.**
- THM-4747 (tc(P_q) table), THM-4746 (Camion–Carmichael).
- `07-reflections/the-cubic-lens-…`, `…-1729-resonance-is-isolated…`.
- `04-computation/ramanujan_near_miss_kpc1.py`.
- HYP-9028.

**Exotic spheres.** MISTAKE-563 and the typed replacements (`collatz_last_dip_cst_alphabet_exotic_20261005.md` §3,
`syracuse_window_alphabet_general_m_20261005.md` §6). The old `surgery_exotic_23.py` also says dim M₁₂(SL₂(Z)) = 3; it
is 2 (MISTAKE-661).

**Sphere packing.** `07-reflections/sphere-packing-eight-twentyfour.md` (oracle, 2026-05-02) was already corrected by
mac-mini (MISTAKE-659 errata banner). S52's catalogue found four of its errors independently:
- the Fano orientations;
- the [23,11,8] code;
- [24,12,8] called perfect;
- the Epstein ζ.

**mac-mini's lane.** THM-4751 / HYP-9362 (`lattices_codes_tournaments_20261021.md`): the circulant H-landscape, class
number as flip distance, and the {7,21} → {8,24} extendable-plane reading (HYP-9362 (d)).

## 9. Errors and refuted leads (this session)

Caught before commit:
- (1²⁵ | 5) as the A24 cusp; it is 2ν, not primitive.
- G₂(A₂¹²) = M₁₁; it is M₁₂.
- Four wrong expectations in my own numerology script:
  - perfect m;
  - |Aut| versus |Aff| at q = 27;
  - code lengths to 4095;
  - "only 31 rare".

Found by the audit (MISTAKE-661):
- "U is the distance to the nearest cusp";
- Borcherds' norm-0 table missed;
- the THM-869/872 and THM-484 overclaims;
- the Nebe "at least" dropped;
- the Smale sentence;
- four misstatements of THM-4748, THM-4749 and HYP-9359.

Refuted leads:
- the e8.pdf trace lock;
- tc(P_q) − 1 = 8m³;
- the projective-plane reading of the H-gaps beyond q = 2, 4;
- Weil-invariance of 1_A − 1_B as a proof of HYP-9363.

## 10. Next obligations

1. **Prove HYP-9363**, for example by specialising the fake-monster denominator formula at the cusps ρ and ν, or by
   finding isometries of M that relate θ_{M+(−1,v)} and θ_{M+(v,−1)}.
2. **Test the 20 other hole types.** Borcherds' table gives all the cusp vectors in these coordinates. The D/E-type
   holes need θ_{N+ρ/h} computed directly.
3. **The n = 24 Niemeier-selection question** of THM-869/872 (unchanged).
4. **Monsuur 2005** (priority for U = c3 + ½|s↑ − ι|²).
