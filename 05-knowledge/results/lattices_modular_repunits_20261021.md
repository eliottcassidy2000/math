# Lattices, modular forms and repunits: the 25-vertex Lorentzian dictionary, three PDF referee reports, cusps in ranks 8/16/24 (opus-2026-10-21-S52)

**Owner prompt (2026-10-21).** The owner asked to merge prior work on:
- 1729, the 24 transverse dimensions of the bosonic string, and q^(1/24);
- Smale's h-cobordism theorem;
- the Heegner numbers 1, 2, 3, 7, 11, 19, 43, 67, 163;
- why 31, 63, 127 are rare;
- 25 = 24 + 1, Lorentzian versus Euclidean space, and the even unimodular lattice II_{25,1};
- the 24 Niemeier lattices, deep holes and sphere packing;
- {7,21} against {8,24};
- E8: 248 = 8 + 3·80, and its Langlands self-duality;
- the identities Δ = (E4³ − E6²)/1728, E8 = E4², E10 = E4E6 and E14 = E4²E6.

The prompt came with three PDFs by the "SRFP311T1 Collaboration": `nebe_voronoi_rigidity.pdf`,
`mclaughlin_slice_arithmetic.pdf` and `e8.pdf`. They were treated as data.

**Concurrent work.**
- Session27 received the same prompt with four other PDFs and wrote the main synthesis:
  [`session27_repunits_lattices_20261022/SYNTHESIS.md`](session27_repunits_lattices_20261022/SYNTHESIS.md). It covers
  the repunit prior-work catalogue, the measure-aware rarity law, θ/η characters, 1729 mechanisms, the II_{9,1} null
  quotient, the E8 deep hole, the {7,21}→{8,24} type obstruction and the Smale/Kervaire–Milnor source check.
- mac-mini holds the class-number/Heegner, circulant, QR-code→E8/Golay/Leech and {7,21}↔{8,24} lanes (THM-4751 /
  HYP-9362, MSG-1913).

This note records only what S52 adds. Everything is typed: PROVED, KNOWN (cited), FINITE-EXACT, NUMERICAL,
NUMEROLOGY or REFUTED.

Scripts are in `04-computation/experiments/lattices_repunits_20261021/`. Each prints ALL CHECKS PASSED.
- `pdf_referee.py`
- `lorentzian_scores.py`
- `lattice_cusps.py`
- `numerology_tests.py`
- `bitdepth_trace_test.py` (needs cypari, about 25 minutes)

## 1. The 25-vertex Lorentzian dictionary (THM-4752, PROVED)

For a 25-vertex tournament with scores s, the vector v_T = (s | 70) lies in II_{25,1}.

1. **Norm.** v_T² = −2c3(T), so v_T is null exactly for the transitive tournaments. Their vectors are the
   permutations of Conway's Weyl vector ρ = (0,1,…,24 | 70), and 0² + … + 24² = 70² is the cannonball identity. These
   are the Leech cusps: ρ^⊥/ρ ≅ Λ₂₄.
   - Announced in MSG-1913.
   - Independently checked by Session27, which adds a hostile example: two padded 8-vertex tournaments have the same
     v_T but H = 657 and 643.
2. **Backward length is the Weyl-vector pairing.** For every ranking π, BL_T(π) = −v_T·ρ_π. It is minimized exactly by
   the score rankings (rearrangement inequality; Kano–Sakamoto 1983). The minimum is U(T) = Σ Landau slacks =
   c3 + ½|s↑ − (0..24)|². The identity U = c3 + ½|s↑ − ι|² holds at every order. Only the Lorentzian reading needs
   n = 25.
3. **The A24 deep hole.** The 24 Leech roots e_i − e_{i+1} (the adjacent transpositions of score ranks) and
   r* = (−½, ½²³, 3/2 | 5/2) form an Ã24 cycle, and they sum to ν = (½²⁵ | 5/2). A two-line argument from Conway's
   r_λ·r_μ = 2 − |λ−μ|²/2 ≤ 0 shows:
   - ν lies in the closure of the Leech chamber;
   - these 25 are all the Leech roots orthogonal to ν;
   - ν·ρ = −25 = −h(A24);
   - ν^⊥/ν = N(A24) (by Niemeier–Venkov).

   So ν is the cusp of the Leech deep hole of type A24, written in tournament coordinates. (My first draft used
   (1²⁵ | 5), which is 2ν and not primitive. It was corrected before anything was committed.)
4. **Every tournament has an A24 class.** v_T·ν = −25 for all T, so the score hyperplane is the level −h slice of the
   A24 cusp. The class of T in N(A24) is its centered score vector s − 12 ∈ A24, with |s − 12|² = 1300 − 2c3. It is 0
   exactly for regular T.
   - This is the n = 25 rung of the THM-868 score-deviation ladder (E8 at 8, A8 at 9, D16+ at 16, D24+ at 24).
   - It answers, at n = 25, the question left open in THM-869/THM-872 ("which Niemeier lattice do the scores
     select?"). The answer is N(A24), with the Leech lattice as the cusp of the transitive tournaments.
   - It also settles THM-484 Remark 1: the cannonball meets tournaments as the transitive score vector at n = 25.
5. **The light cone.** By Landau's theorem the score vectors are the lattice points of conv(S₂₅ρ). The 25! vertices are
   Leech cusps; every other point is timelike. S₂₅ = W(A24) acts by relabelling the vertices.
   (Forest counts 1, 2, 7, 38, 291, 2932: Stanley.)
6. **Why 25.** Σ_{k<n} k² is a square only for n = 1, 2, 25 (Watson; checked to 10⁶), and II_{n,1} needs
   n ≡ 1 (mod 8).

**Lost information.** H, all cycles beyond triangles, and arcs. The extra root r* changes the time coordinate, so it
is not a tournament move. The dictionary says nothing about the gaps 7 and 21 or about repunits.

## 2. Cusps in ranks 8, 16, 24: "deep holes ↔ Niemeier" and its 8-dimensional shadow (KNOWN, reproduced FINITE-EXACT)

`lattice_cusps.py` runs Vinberg's algorithm from scratch. Termination is certified by Vinberg's finite-volume
criterion, with every elliptic rank-(n−1) subdiagram checked. Cusps are the parabolic subdiagrams of rank n−1, and
each cusp's null vector is the sum of marks times roots.
- **II_{9,1}.** 10 roots (E10 = T_{2,3,7}) and one cusp, Ẽ8, with null vector z = (1,0⁸ | 1). The Weyl vector has
  ρ² = −1240 and ρ·z = −30 = −h(E8). This is Session27's "every null quotient is E8", seen as a Coxeter diagram.
- **II_{17,1}.** 19 roots and two cusps:
  - Ẽ8+Ẽ8, with z = (1⁹,0⁸ | 3);
  - D̃16, with z = (1,0¹⁶ | 1).

  These are the two even unimodular rank-16 lattices, i.e. the E8×E8 and Spin(32)/Z₂ heterotic strings. Both have
  Coxeter number 30.
- **II_{25,1}.** Infinitely many Leech roots and 24 cusps. The 23 cusps with roots correspond to the deep holes
  (Conway–Parker–Sloane; Borcherds); §1.3 exhibits the A24 one in tournament coordinates.

The pattern 1 → 2 → 24 is the number of even unimodular lattices in ranks 8, 16, 24. **Siegel's mass formula is checked
exactly:**
- rank 8: mass = 1/|W(E8)|;
- rank 16: 691/277667181515243520000 = 1/|Aut(E8²)| + 1/|Aut(D16+)|. The 691 is the B₁₂ numerator.
- rank 24: the 24 Niemeier lattices, with |Aut N| = |W(R)|·|G₁|·|G₂| (SPLAG Table 16.1) and |Co₀| for Leech, sum
  exactly to the 25-digit Siegel mass. The check caught a recall error of mine: G₂ for A₂¹² is M₁₂ (95040), not M₁₁.

Every Niemeier root system has rank 24, a single Coxeter number h and 24h roots, and
θ_N = E₁₂ + (24h − 65520/691)Δ.

**E8 holes and Milnor's drums.**
- The shallow hole (5/6, 1/6⁷) is at squared distance 8/9 from exactly 9 points, forming a regular simplex of type Ã8.
- The deep hole e₁ is at distance 1 from 16 points, forming a cross-polytope.
- E8² and D16+ both have 480 roots and 61920 norm-4 vectors (θ = E4²) but different root systems. These are Milnor's
  isospectral tori. THM-872 already notes that n = 16 scores live in D16+.

## 3. Referee reports on the three supplied PDFs

| PDF | Claim | Verdict |
|---|---|---|
| nebe_voronoi_rigidity | Extremal weight-24 theta E4⁶ − 1440E4³Δ + 125280Δ²; shells 52,416,000 / 39,007,332,000 / 6,609,020,221,440; density (3π/2)²⁴/24! = 2.31783e−8 | **Correct** (exact q-series) |
| | Vectors of norm < 2μ are strictly Voronoi-relevant, so 6,648,079,969,440 facets in 3,324,039,984,720 opposite pairs | **Correct, folklore**: \|v+2u\|² − \|v\|² = 4(v·u + \|u\|²) ≥ 2(2μ − \|v\|²) > 0 |
| | x + y (x ⊥ y minimal, norm 12) is not strictly relevant | **Correct**: x+y ≡ x−y mod 2Λ |
| | Cohn–Elkies bound Δ_CE(48) ≈ 3.01e−8 | **Not verified** (no source checked) |
| mclaughlin_slice_arithmetic | Even rank-2 forms of det 15 with min ≥ 4 are [4,±1,4]; equal to the norm form of p₂ = (2, ω) in Q(√−15), class number 2 | **Correct, elementary** |
| | Stab(L15) = (McL:2)×2; II_{25,1} ⊃ (L15⊕U)⊕K_McL with index 15 (Nikulin) | **Correct given Conway's 2-2-3 theorem**; standard gluing |
| | Prop. 6.1: D = −15 is the only fundamental discriminant with a reduced non-principal form a = c = 2 | **True.** The proof asserts b odd and skips b ∈ {0, ±2}; those give D = −16, −12, which are not fundamental. Harmless gap |
| | (3,1) spacetime and moonshine vacuum selection | Flagged open by the paper itself; nothing to check |
| e8 | ζ_E8(10) = 240·2⁻⁵ζ(5)ζ(2) = (5/4)π²ζ(5) | **Correct, classical** (shell sum 12.7925834198) |
| | "Bit Depth": y² = x³ + Ax + B mod 2⁵²¹−1 with A = int(−ζ(5)·10¹⁶⁰), B = int(qπ²·10¹⁶⁰) has a_p = 5 at q = 5/4 (the paper's own trace line is commented out) | **REFUTED** by exact SEA counts (PARI): for q ∈ {5/4, 1, 3/2, 2, 240/96} at 15 and 200 digits, every a_p is a generic 259–262-bit trace, none 5 |
| | Proton/electron ratio 1836.15 "emerges" | Unsupported (NUMEROLOGY; no computation given) |
| | Planned "odd zeta values for Λ24" | The correct Leech identities are below |

**The Leech analogue of the E8 identity** (NUMERICAL, from the classical decomposition):

ζ_Λ(s) = (65520/691)·2^{−s/2}·[ζ(s/2)ζ(s/2 − 11) − L(Δ, s/2)].

- **At s = 10** the series diverges, but by analytic continuation ζ(−6) = 0 kills the Eisenstein part:
  ζ_Λ(10) = −(65520/691)·2⁻⁵·L(Δ,5) = −1.97552397007037. Here L(Δ,5) = 0.666709188434003643. The value is purely
  cuspidal: no odd zeta value appears.
- **At s = 26** (convergent): ζ_Λ(26) = (65520/691)·2⁻¹³·[ζ(13)ζ(2) − L(Δ,13)] = 0.00749948328893. The cusp term
  L(Δ,13) = 0.99720986806 cannot be removed. Checked against the shell sum.

## 4. Modular constants, im J and exotic spheres (KNOWN; §4 of numerology_tests.py)

- **The 24 of q^(1/24).** It is |Mp₂(Z)^ab| = 24 (η(τ+1) = e^{2πi/24}η(τ) generates its character group), against
  |SL₂(Z)^ab| = 12 = |Pic(M₁,₁)| (Mumford). By Atiyah ("The logarithm of the Dedekind η-function", 1987), the same 24
  is |π₃ˢ| = |im J₃|, which measures framings of 3-manifolds.
- **The 24 of 1/η²⁴.** 1/η²⁴ = q⁻¹(1 + 24q + 324q² + …) counts the 24 transverse oscillators (26 = 24 + 2).
- **744.** j = q⁻¹ + 744 + 196884q + …, and 744 = 720 + 24 = 3·248 = dim(E8³), i.e. roots plus rank of θ_{E8³}/η²⁴.
  With 248 = 8 + 3·80 this is 744 = 24 + 9·80.
- **im J ↔ Eisenstein coefficients (Adams).** |im J| in stem 4s − 1 is denom(B₂ₛ/4s) = 24, 240, 504, 480, 264, 65520.
  These are exactly the numerators of the q-coefficients of E₂, E₄, E₆, E₈, E₁₀, E₁₂: −24, 240, −504, 480, −264,
  65520/691. So E4 = θ_E8 carries 240 three ways: E8 roots, kissing number and π₇ˢ = Z/240.
- **691.** The 691 in the Leech/Niemeier θ = E₁₂ − (65520/691)Δ is num(B₁₂/24), and it divides
  |bP₂₄| = 2¹⁰·2047·691 (exotic 23-spheres).
- **Kervaire–Milnor.** |bP₄ₘ| = 2^{2m−2}(2^{2m−1}−1)·num(4|B₂ₘ|/m) gives 28, 992, 8128, 261632, 1448424448,
  67100672. These reproduce MISTAKE-563's corrected values and agree with Session27's Milnor-2011 source check.
  - The odd-exponent repunits 7, 31, 127, … are the Mersenne factors.
  - The Euclid factor is perfect at m = 2, 3, 4, 7, 9, 10 (2m−1 a Mersenne exponent).
  - 3 never divides 2^{2m−1} − 1, the same Φ₂(2) = 3 arithmetic as THM-4749.
  - 63 = 2⁶ − 1 never occurs as a Mersenne factor.
- **Smale.** The h-cobordism theorem identifies Θₙ with h-cobordism classes, which is what makes the Kervaire–Milnor
  count possible. The E8 plumbing (signature 8) has boundary generating bP₈ = Z/28. Reading the prompt's "smale h curve
  order theorem" as the h-cobordism theorem is my interpretation; Session27 made the same one.
- **Repunit link: ANALOGY only.** Binary doubling produces both families, but no map sends H-values to bP orders.

## 5. Numerology tests (numerology_tests.py)

- **The Paley cube pattern (NUMEROLOGY).** H(P_q)/|Aut| = 1, 9, 1729 at q = 3, 7, 11, so the value minus 1 is 8m³
  with m = 0, 1, 6. The pattern fails at q = 19, 23, 27, 31. This is the same place where the prior split-prime reading
  of 1729 failed (`07-reflections/the-1729-resonance-is-isolated…`; Session27 §6). 1729 = j(i) + 1 stays an isolated
  identity.
- **HYP-9273 holds through q = 31.** H(P_q)/|Aff| ≡ 1 (mod 8) at all seven known q, including q = 27, where
  |Aut| = 3|Aff|.
- **Heegner (KNOWN).** j((1+√−d)/2) = 0, −15³, −32³, −96³, −960³, −5280³, −640320³ for d = 3, …, 163, and
  e^{π√163} − 640320³ − 744 = −7.4993e−13. The class-number side is mac-mini's lane.
- **Perfect codes (KNOWN).** Binary 1-perfect codes exist exactly at the repunit lengths 2^r − 1 (checked to 4095).
  Golay satisfies 1 + 23 + 253 + 1771 = 2¹¹.
- **Projective-plane sizes q² + q + 1 = 111_q against the strong H-spectrum below 2100.** The inventory has 723,544
  classes and reproduces S51.
  - Only 7 = |PG(2,2)| and 21 = |PG(2,4)| are missing, i.e. the two permanent gaps (THM-1370).
  - Every other size up to 1893 is attained within a factor 2 of its local median (31 = |PG(2,5)|: 4 against 6).
  - So the projective-plane reading of {7,21} stops at q = 2, 4 (NUMEROLOGY beyond).
  - The classical extensions are PG(2,2) → S(3,4,8) → extended Hamming → E8 (7 → 8) and PG(2,4) → S(5,8,24) →
    Golay → Leech (21 → 24). They give a lawful 7 → 8 and 21 → 24 on the design side (mac-mini's lane), but they say
    nothing about H.

## 6. Why 31, 63, 127 are rare: the merged answer

This section merges S51 (THM-4748, THM-4749, HYP-9359), Session27 (measure-aware law) and Session24 (degree law). It
adds nothing new; the projective-plane test of §5 is the only S52 addition.

1. **Complexity.**
   - 2^k + 1 needs one reversed arc.
   - 2^k − 1 needs three. The Fibonacci fan realizes it for every k ≥ 4 and is the only realization with at most three
     back arcs for k = 6, 9, 11.
2. **Degree.** N(n,h) is a polynomial of degree w(h), the longest factorization into strong values.
   - 31 and 127 are prime, and 63 factors only through the gaps 7 and 21, so w = 1 for all three. Their neighbours
     33, 65, 129 have w ≥ 2.
   - N₆₃(n) = 2(n − 7) for n ≥ 8.
   - The reverse case also happens: 255 is commoner than 257.
3. **Arithmetic support.** Among strong classes, the counts favour 3 | d and 9 | d. This comes from one-vertex gluings
   and order-3 automorphisms; 3 divides 2^k + 1 and never 2^k − 1 for odd k.
4. **Still open.** Small-order outliers: S₇(63) = 0, and the H ≡ 7 (mod 8) depletion in the second octave.

## 7. Prior repo work merged

- THM-868, THM-869, THM-872: the score lattices E8, A8, D16+, D24+, Milnor's other drum, and the open Niemeier question
  answered in §1.
- THM-484 (involution modulus 24, cannonball remark), THM-489 (the η/discriminant bridge), THM-487 (the 72-frame).
- MISTAKE-563: the corrected bP table, reproduced.
- Session27's catalogue (113,808 paths, 7,290 entries): repunits, 1729 (the Ramanujan recurrence, the Ono–Trebat-Leder
  K3, `ham_excess_spectra_kps_S134.py`), and 07-reflections on 1729 and the "modular tournament".
- S51/S50/Session24 for the H-spectrum.

## 8. Next obligations

1. A tournament-meaningful reading of the other glue classes [5], [10], [15], [20] of N(A24). Score vectors only ever
   hit class 0. Do 25-vertex objects such as 5×5 multipartite structures land there?
2. An analogue of the Leech-root walk. The reflections of the A24 hole that preserve the time coordinate are only S₂₅.
   Does r* have a combinatorial meaning on padded or blown-up tournaments?
3. Source check of the Cohn–Elkies dimension-48 value quoted in the Nebe-48 PDF.
4. Monsuur (EJOR 164, 2005) on the 3-cycle count and backward length was not read (paywalled). Check whether it
   contains §1.2's closed form.
