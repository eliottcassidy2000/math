# Tournaments, Collatz and 3x−1, and the arithmetic of 2 and 3

Session mac-mini-2026-11-15-synthesis (wave 38 of this packet). **No open problem is solved here, and Chan–Pak
remains OPEN.**

This is a companion to the concurrent
[mod-8 dictionary](../chanpak_mod8_synthesis_span3_20261115/SYNTHESIS.md) of klein-2026-11-15-synthesis, which
answers the same owner prompt. Read that file for:
- (ℤ/8)* = ⟨−1, 3⟩ and its three characters;
- the Zolotarev signs of the Collatz multipliers;
- the Paley sign rules;
- Artin's constant for one-cycle doubling tables;
- a first open-problem table.

This file adds:
- the composition algebra, including a new three-block law;
- the 2-adic ladder beyond mod 8, including the metagraph reading mod 16;
- the merge with the earlier multiplication-table work (THM-4792);
- the archimedean obstruction;
- a local–global plan for Chan–Pak;
- further bridges to open problems.

Every statement is labelled PROVED, KNOWN (cited), FINITE-EXACT, NUMERICAL, CONJECTURE or ANALOGY.

## 0. The picture in one paragraph

Hamiltonian-path counts are values of multilinear forms in port vectors. Under the basic composition laws these forms
carry the arithmetic of ℤ[√−2], the ring in which 2 ramifies (2 = −(√−2)²) and 3 splits (3 = (1 + √−2)(1 − √−2)).
The 3-cycle is the prime 1 + √−2: H(C₃) = 3 = N(β(point)). Collatz is dynamics built from the same two primes 2 and 3.

The two subjects meet at every place of ℚ:
- **at 2:** Rédei's parity, the character χ₋₂ and the 2-adic norm group of ℚ₂(√−2);
- **at the odd primes:** the multiplication tables mod p, whose 3x+1 / 3x−1 alternation is the sign of halving on
  the Paley tournament;
- **at ∞:** angles in ℤ[√−2], where norms are thin and bilinear forms are not.

Chan–Pak then reads as a local–global problem for a thin multilinear family, and Collatz as a dynamics problem in the
same local rings. That is a shared shape, not a reduction of one problem to the other.

## 1. One object, three composition laws

For a tournament X:
- P_k(X) counts ordered k-tuples of disjoint nonempty paths covering X; the port vector is v_X = (P_1, P_2, …);
- the slot polynomial f_X(x) = Σ_k P_k(X)·C(x, k) has the Grinberg–Stanley odd-cycle expansion
  Σ_C 2^|C| x^(n−d(C));
- H = f(1), and Q = Π₂ = f(2) counts ordered pairs of possibly empty disjoint paths
  ([pair-line census](../chanpak_pair_line_census_20261115/PAIR_LINE_STRUCTURE.md) §4).

| composition | count | algebra | status |
|---|---|---|---|
| ordinal sum A ⇒ B | f_(A⇒B) = f_A·f_B pointwise | split torus: (H, Q) multiply coordinatewise | PROVED |
| point-cyclic join A → B → • → A | Re(β_A β̄_B), β = (H + P₂) + H√−2 | definite form of discriminant −8, Gram [[3,1],[1,1]] | PROVED: THM-5172 (C), THM-5222 |
| diagonal join (B ≅ A) | N(β_A) = (H + P₂)² + 2H² | norm form of ℤ[√−2] | PROVED |
| three-block cyclic join | T(v_A, v_B, v_C), symmetric trilinear; tridiagonal Jacobi Gram J(v_C) with J_aa = p_(a−1) + 3p_a + p_(a+1) and J_(a,a+1) = p_a + p_(a+1) | a cubic form generalising the binary one | **PROVED here, THM-5242 (A)**; 400/400 brute-force checks |
| double cone a → B → c → a | P_D(q) = q[P(q+1) + qP(q) − P(q−1)] | a linear operator on path-cover polynomials | PROVED, [LLT.md](../chanpak_generic_lines_llt_20261115/LLT.md) |

Two readings of the table:
- **Positivity selects the definite form.** Because Q ≥ 2H (wave 35), every realised β lies in a cone avoiding the
  isotropic directions of the split (Pell) forms. That is why cyclic joins carry the anisotropic ℤ[√−2] arithmetic
  and never a Pell supply.
- **The three-block law is symmetric.** H(cyc(A, B, C)) = H(cyc(B, A, C)), so reversing the triangle never changes
  the count. With two points it is the triangle reader 3H + P₂; for blocks of order at most 2 its cube is
  3(P₁ + P₂)(P₁² + P₂²).

## 2. The 2-adic ladder: where "mod 8" lives

| level | tournament side | Collatz side | status |
|---|---|---|---|
| mod 2 | H is odd (Rédei 1934) | 3x ± 1 is even for odd x (a forced halving) | KNOWN |
| mod 4 | H ≡ 1 + 2·(number of directed odd cycles) (odd-cycle expansion) | v₂(3x + 1) = 1 iff x ≡ 3 (mod 4) | PROVED / elementary |
| mod 8 | N(β) ≡ 3 for every diagonal join; the primes x² + 2y² are those in ker χ₋₂ = ⟨3⟩ (klein) | 3^k ≡ 1, 3 (mod 8); v₂(3x + 1) is fixed by x mod 8 (klein) | PROVED |
| mod 16 | at odd order N(β_A) ≡ 3H² + 8·HC(A), so **the norm reads the parity of the closure number** of A's node in the tiling metagraph (THM-4729): its count of closed tilings, i.e. Hamiltonian cycles | v₂(3x + 1) ≥ 4 iff x ≡ 5 (mod 16) | **PROVED here, THM-5242 (B)**; 300/300 checks |
| mod 32 | at even order N(β_A) ≡ 3H² + 16·c_(n−1)(A) | the next 2-adic digit of x | PROVED here |
| 2-adic limit | every diagonal join is a 2-adic norm from ℚ₂(√−2) | **the closure of the Collatz group ⟨2, 3⟩ is exactly N(ℚ₂(√−2)^×)**; negation, which exchanges 3x+1 and 3x−1, is the non-norm class: (−1, −2)₂ = −1 | **PROVED here, THM-5242 (D)** |
| cycles | — | every cycle denominator 2^K − 3^k with K ≥ 3 is a non-norm. Only the trivial cycle 1 → 4 → 2 → 1 (and its 3x−1 mirror) has a norm-class denominator | PROVED here |

**Reading.** "Mod 8 in Collatz" and "mod 8 in the norm form" are one structure: the norm residue symbol of
ℚ₂(√−2).
- Its kernel contains both Collatz multipliers, 2 = N(√−2) and 3 = N(1 + √−2).
- On the tournament side, the generators are the point (Q = 2) and the 3-cycle (H = 3).
- The mirror −1 that turns Collatz into 3x−1 is exactly what the norm group misses.

Going one digit further:
- the tournament norm starts reading cycle structure, through the metagraph's closed tilings;
- the Collatz argument starts reading deeper valuations.

ANALOGY: László Rédei proved both the parity of H (1934) and, with Reichardt, the 4-rank law for class groups of
quadratic fields. Higher 2-power ranks are governed by fields such as ℚ(√−2), ℚ(i) and ℚ(√2) (Rédei symbols; Smith's
work on 2^∞-Selmer groups). The ladder above is the tournament-side analogue of such a 2-power tower. No theorem is
claimed.

## 3. The odd places: multiplication tables, Paley tournaments, circulant cycles

1. **The centre of the table (THM-4792, PROVED).** For odd m = 2k + 1 the centre labels of the multiplication table
   are k = −1/2 and k + 1 = 1/2, and their products are ±1/4. The map 3x−1 sends 1/4 to −1/4, and 3x+1 sends it back.
   The halving-shift alternation of the edge tables of 3x+1 and 3x−1 exists iff −1 ∈ ⟨2⟩ mod m (kind-pasteur
   Theorem A; Moree 1997). It holds for all primes ≡ 3, 5 (mod 8), never for primes ≡ 7 (mod 8), and depends on
   higher power residues for primes ≡ 1 (mod 8).
2. **Alternation is the Paley sign (THM-5242 (C), PROVED).** For a Paley prime p ≡ 3 (mod 4), these are equivalent:
   - the tables of 3x+1 and 3x−1 alternate mod p;
   - halving reverses the Paley tournament P_p;
   - p = x² + 2y²;
   - p ≡ 3 (mod 8).

   klein proved the Paley and x² + 2y² half. The alternation half is the merge with the owner's earlier
   multiplication-table thread. Checked on all 332,397 Paley primes below 10^7. **Every prime diagonal-join count is a
   prime of this kind.**
3. **Densities.**
   - By Hasse (1966) the alternating primes, those with ord_p(2) even, have density 17/24. We count 0.70817 below
     10^7 (KNOWN, FINITE).
   - The primes whose doubling table is a single cycle have Artin's density A = 0.3739… under GRH (Hooley). klein
     measured 2A within primes ≡ 3, 5 (mod 8).
4. **Circulants are rational Collatz cycles (THM-4849, PROVED).** A circulant tournament on ℤ_q is a rational cycle of
   the Terras map with denominator 2^q − 3^((q−1)/2). By THM-5242 (D) this denominator is ≡ 5 (mod 8) at every Paley
   prime q: it lies in the non-norm class, opposite to the counts of diagonal joins.
5. **Class numbers (THM-4666, PROVED).** For Paley primes, the ranking excess of the natural order is a class number.
   The pair-block order has excess (1 − 2(2/p))·h(−p). That is 3h(−p) exactly at the x² + 2y² Paley primes and
   −h(−p) at the others: the same mod-8 split, appearing in the Euler factor 2 − χ(2) of Dirichlet's formula.
6. **Paley data (FINITE-EXACT, computed independently by klein and here).** H(P_p)/(p(p−1)/2) = 1, 9, 1729,
   6,857,869,865 for p = 3, 7, 11, 19. Q(P_p) is also divisible by p(p−1)/2, because automorphisms act freely on path
   pairs too. The Paley tournament P₇ attains the minimum of Q/H at order 7.

## 4. The archimedean place: norms are thin, bilinear forms are not (HYP-9792)

- **Diagonal joins are a density-zero source of primes.**
  - In the complete range below 4,190,211 (blocks of order ≥ 11 have P₂ ≥ 2^11 − 2), only 26 of the 73,940 primes
    ≡ 3 (mod 8) are diagonal joins: 3, 11, 227, 35963, …
  - The reason is an angle obstruction. The ratio x/y = Q/H − 1 is at least about 2n − 5 at order n (exact minima for
    n ≤ 10), and the size of p forces n ≍ log p / log log p.
  - So the reachable primes lie in a cone of angle O(log log p / log p) in ℤ[√−2], a vanishing share by Hecke.
  - This archimedean obstruction comes on top of the mod-8 reciprocity obstruction klein compared with Apollonian
    sub-families.
- **The bilinear family is locally free.**
  - Over the explicit two-apex sources
    ([TWO_APEX_COARSE_PLACEMENT.md](../chanpak_pair_line_census_20261115/TWO_APEX_COARSE_PLACEMENT.md)), (H, Q) mod ℓ
    is within total variation 0.005–0.013 of uniform at order 12 for ℓ = 5, 7, 11, 13, halving per spine vertex.
  - At ℓ = 2, 3 the proved laws are the only constraints.
  - The point-cyclic join Re(β_A β̄_B) is a bilinear (Type II) form over ℤ[√−2], the shape that Friedlander–Iwaniec
    (X² + Y⁴) and Heath-Brown (x³ + 2y³) needed to capture primes in thin sets.
- **Small blocks are not enough.** With both blocks of order ≤ 8, point-cyclic joins miss 351 of the 500 odd numbers
  below 1,000 and 90% of those in [10^5, 2·10^6]. Coverage has to come from large or combined families, as in the
  placement packets.

## 5. Where this sits among open problems (additions to klein's table)

| open problem | meeting point | kind |
|---|---|---|
| Hasse density of even order / Artin's conjecture (base 2) | alternation of the 3x±1 tables ⟺ ord_p(2) even, density 17/24 (Hasse, unconditional); one-cycle doubling tables (Artin, GRH) | proved bridge; the Artin part conditional |
| Ramanujan's form x² + y² + 10z² (Ono–Soundararajan, under GRH) | another "every odd number except a finite list" theorem. Its 18 exceptions are 3, 7, 21, 31, 33, 43, 67, 79, 87, 133, 217, 219, 223, 253, 307, 391, 679, 2719 (verified to 20,000). **That 7 and 21 appear in both lists is a coincidence; no claim.** The shared shape is local data plus analytic equidistribution plus finite computation, with GRH entering through Siegel zeros | analogy |
| Duke–Schulze-Pillot (ternary forms, large N, ineffective) | the likely shape of a Chan–Pak proof: every sufficiently large odd N by an analytic argument, then a finite computation | analogy |
| Bhargava–Hanke 290 theorem | a finite universality certificate. The 168 carrier forms covering [333, 10^9] (HYP-9722) play this role empirically | analogy |
| Tunnell's theorem / BSD / congruent numbers | Tunnell counts x² + 2y² + 8z² and x² + 2y² + 32z². The binary core x² + 2y² is the diagonal-join norm. Conductor 32 also appears in the Paley K6 − P6 sum (THM-4653, HYP-9276) | meeting point of forms; analogy |
| Friedlander–Iwaniec, Heath-Brown (primes in thin sets) | point-cyclic joins are bilinear over ℤ[√−2]. Type I input: HYP-9792 (ii). Type II input: open | method transfer (proposed) |
| Hadamard / skew-Hadamard conjecture | Paley tournaments give skew Hadamard matrices of order p + 1. Doubly regular tournaments are equivalent to skew Hadamard matrices (THM-5136/5137). The same P_p carry the alternation sign of §3 | meeting point |
| Gauss class numbers (solved: Heegner–Baker–Stark) | Paley ranking excess = h(−p); the mod-8 factor 1 − 2(2/p) (THM-4666) | proved bridge |
| Furstenberg ×2×3 / effective ×a×b (Bourgain–Lindenstrauss–Michel–Venkatesh) | the multiplication tables by 2 and 3 mod m are the finite ×2×3 systems. ⟨2, 3⟩ ∩ (squares) acts on P_p by automorphisms | analogy / method |
| Collatz (3x+1) and 3x−1 | the 2-adic norm-group reading: ⟨2, 3⟩ never contains −1 2-adically, but contains it mod p for every p ≡ 3, 5 (mod 8). A local–global contrast; no reduction | proved identities; the rest analogy |

## 6. What the synthesis suggests for Chan–Pak

Read Chan–Pak as a local–global problem for multilinear tournament families.

| place | condition on counts | source |
|---|---|---|
| 2 | H odd; Q ≡ 0 (mod 4), and Q ≡ 0 (mod 8) at even order; fibre steps 12 or 24 | Rédei; pair-line census §4; wave 37 |
| 3 | Q ≡ (−1)^n H (mod 3) | THM-4756 |
| ℓ ≥ 5 | no obstruction (the two-apex family equidistributes) | HYP-9792 (ii) |
| ∞ | placement within O(N^(3/8)) of any large N | [target-address packet](../chanpak_target_address_20261115/README.md) |

What is missing is the **global step**: exact absorption of the archimedean residual. Wave 38 tested the
fibre-swap repair of wave 37 against the two-apex construction (README, wave 38).
- **Exact but coupled.** Modules inside a two-apex source are read multilinearly, never additively. The mixed term
  was nonzero in 240/240 tests, about 9% of a single swap.
- **Coarse.** A module swap moves the count by a constant fraction of itself: 180 module pairs give 25 values spread
  over about 3·10⁹ steps.

**Design principle (thin readers).** An exact repair needs a reader whose covers that split the repair block have
essentially unique completions elsewhere. The splitting multiplier must be O(1), as in the triangle reader (M = 1).
Hosts in which the rest of the tournament multiplies the split covers turn every swap into a coarse move.

**Two consequences.**
- The norm route (diagonal joins) cannot reach most primes (HYP-9792 (i)).
- The bilinear route needs (a) a thin reader for the residual and (b) a Type II estimate for
  Re(β_A β̄_B) over two independent explicit families. Its Type I input, equidistribution modulo small primes, is
  already visible.

## 7. Receipts

| claim | script (in `scripts/`) | output |
|---|---|---|
| triangle law, symmetry, Gram matrices | `w38_triangle_form.py` (with `w38_hq.c`) | 400/400; Jacobi matrices printed |
| mod-16/32 law, Paley equivalence, Hasse, 2-adic norm group, Paley data, diagonal-join census | `w38_norm_bridges.py` | 300/300; 0 exceptions in 332,397 primes; 0.70817; True; 26 of 73,940 |
| two-apex equidistribution | `w38_two_apex_equidistribution.py` | TV table in HYP-9792 |
| bilinear coverage by small blocks; Ramanujan control | `w38_bilinear_cover.py` | misses listed; the 18 exceptions |
| fibre widths (exact and sampled); two-apex additivity | `w38_fibre_widths.py`, `w38_sample_pairs.c`, `w38_sampled_fibres.py`, `w38_two_apex_additivity.py`, `w38_two_apex_coupled_fill.py` | README wave 38 |
