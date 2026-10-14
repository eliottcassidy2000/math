# Numbers given tournament structure

**Contents.**
- Fermat primes are the primes whose circulant tournaments are all rigid. Wieferich primes are where doubling-invariant tournaments on `Z_(p^2)` explode (THM-4644); 3511 is a Paley prime.
- Aliquot parity flips only at squares, so mixed-parity aliquot cycles need squares (HYP-9273 A, which contains the opposite-parity amicable problem).
- Odd aliquot cycles are frustrated: Dickson's 1913 theorem read as signed-cycle balance, the structural reason odd lengths are rare.
- `{deficient, perfect, abundant}` reads as `{hyperbolic, Euclidean, spherical}` at four divisors, along the `(2,3,k)` chain through 6 and 42.
- Paley graphs on Fermat primes and Ramsey.
- Paley Hamiltonian orbit counts `≡ 1 (mod 8)` (HYP-9273 B).

2026-10-14, session opus-2026-10-14-S33 (worktree `codex/session-number-tournament-bridges-20261014`).

Owner's prompt (verbatim): "keep working the most precise and meaningful hypotheses you can conjecture in your continual pursuit toward proving a major open result, think about merging in how {3,5,27,257,65537} might relate with the platonic solids or other groups of 5 we have studied and also how 1093,3511 wieferich primes relate with 7,21 as tournament hamiltonian path count forbidden values. consider also aliquot sums, draw creative analogies between {deficient, abundant, perfect} and other triples we have studied previously. think about how sociable cycles of length 1 (perfect) 2 (amicable) 5 and 28 have been found, but not 3. keep trying to draw creative analogy between number theory and our tournament graph theory novel results and techniques developed for discrete math, try to give the numbers appropriate structure to unlock problems. maybe prove that an amicable pair is impossible with one number even and the other odd. in the reproduced Parts graph, the radius-two cross-edges alone allow four colours. Adding just six cross-edges of one specific geometric type makes the SAT solver report that four colours are impossible. This isolates a small interaction inside a known five-chromatic graph"

**Status.**

- **PROVED:**
  - THM-4644 ([canon file](../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md)):
    - Fermat primes are exactly the primes all of whose circulant tournaments are rigid;
    - Wieferich primes are exactly the primes at which the circulant tournaments on `Z_(p^2)` invariant under the odd part of doubling jump from `2^((p−1)/d')` to `2^((p+1)(p−1)/(2d'))`.
  - The aliquot parity lemma and its corollaries (section 3; the amicable case is classical).
  - The cycle identities `Σ σ(n_i) = 2 Σ n_i` and `Π σ(n_i)/n_i >= 2^k`.
  - The frustration lemma for odd cycles (section 3b; Dickson's mechanism) and its reading through the unoriented incidence matrix of `C_k`.
- **FINITE-EXACT:**
  - all aliquot cycles with every member `<= 10^8`;
  - the double flips for odd `a <= 10^9`;
  - the frustration census on those cycles and on seven published cycles of lengths 6, 8 and 9;
  - THM-4644's checks (brute force for `p <= 11`, stabilizers for `p <= 31`, `E` for `p < 4000`);
  - `H(P_q)` for `q = 3, 7, 11, 19, 23` and `H(P_27) mod 8·1053`;
  - all circulant tournaments for `p <= 19`;
  - the doubly regular tournament of order 15;
  - the Paley clique numbers;
  - primary pseudoperfect numbers below 60000.
- **CITED:**
  - Burnside (1901); Alspach (via THM-871);
  - Wieferich (1909), Meissner (1913), Beeger (1922), PrimeGrid (OEIS A001220);
  - Greenwood–Gleason (1955), Evans–Pulham–Sheehan (1981), Angeltveit–McKay (2024);
  - Gioia–Vaidya (1967; rediscovered by D'Abramo 2005);
  - Hagis (1970);
  - Iannucci–Luca (2009); Pollack (2011);
  - Dickson (1913);
  - Kobayashi–Pollack–Pomerance (2009); Pollack, "A remark on sociable numbers of odd order" (Lemma 5, the sharp gcd restriction);
  - Bratley–Lunnon–McKay (1970); Moews and Pedersen;
  - Chernykh (OEIS A066873);
  - Poulet (1918);
  - Rédei (1934), Moon (1966);
  - Hurwitz (1893), Klein (1879).
  - From this repository: THM-871, THM-307, THM-200/THM-029, THM-1370, THM-4523, THM-4631, THM-4641, the S10 Fermat/Platonic verdict, the camion/Busch gap note, the arithmetic-seams note (2026-09-21), and mac-mini's Hadwiger–Nelson type-M note.
  - **The same-day codex packets on the same owner themes** (section 1b): the [dyadic aliquot three-cycle frontier](aliquot_dyadic_frontier_20261014e.md), the [guarded core compiler](aliquot_guarded_cycles_20261014e.md), [tournament path defects](tournament_path_defects_20261014e.md), and [Fermat/Wieferich order lifts](fermat_wieferich_order_20261014e.md).
- **OPEN, with evidence:** HYP-9273 ([file](../hypotheses/HYP-9273-aliquot-square-avoidance-paley-hamiltonian-orbits-number-tournament-trichotomies.md)): (A) aliquot square avoidance; (B) Paley orbit counts `≡ 1 (mod 8)`.
- **ANALOGY / NUMEROLOGY:** typed where they occur.
- **Not yet audited** (section 9).

## 0. Answers in brief

1. **"Maybe prove that an amicable pair is impossible with one number even and the other odd."** This is a classical **open problem**. What can be proved, and what we add:
   - **Structure (Gioia–Vaidya 1967; rediscovered by D'Abramo 2005).** If `σ(m) = σ(n) = m + n` with `m` even and `n` odd, then `σ(n)` and `σ(m)` are odd. Since `σ(k)` is odd iff the odd part of `k` is a square, `n` is an odd square and `m` is a square or twice a square. Conversely, any amicable pair with a member in `□ ∪ 2□` has opposite parity. The problem is therefore exactly: **is there an amicable pair of squares or twice-squares?**
   - **Literature.** Iannucci and Luca (2009) proved that `s(n^2)` is a square for at most `x/(log x)^(3/2+o(1))` values `n <= x`; Pollack (2011) notes that the same holds for twice a square. Pollack (2011) also proved that amicable `l`-full numbers are `o(x^(1/l))`. Hagis (1970) proved that coprime opposite-parity pairs have both members above `10^60`. Chernykh's distributed search has finished all pairs with smaller member below `10^20`, and every known pair has equal parity.
   - **Our reduction.** The odd member is `a^2`, so we search over `a`. Every odd `a <= 10^9` was checked, so there is no opposite-parity pair with odd member `<= 10^18`. This range is inside Chernykh's; it is an independent check, not a new bound.
     - Only 4476 values have `s(a^2) ∈ □ ∪ 2□` ("double flips"), far below the Iannucci–Luca bound (about `10^7` at `10^9`).
     - 4458 of them are primes. For a prime, `s(a^2) = a + 1` is far too small to close a pair.
     - The 18 composite ones grow slowly: 1, 2, 3, 6, 13, 15, 17, 18 up to `10^2, …, 10^9`. None closes a pair.
   - **The generalisation (HYP-9273 A).** Every aliquot cycle avoids `□ ∪ 2□`. Our census of all cycles with members `<= 10^8` finds no exception.
2. **Sociable cycles of length 1, 2, 5, 28 but not 3.**
   - The known lengths are 1, 2, **4**, 5, **6**, **8**, **9**, 28 (Wikipedia's table, October 2024). They comprise 5398 four-cycles, five 6-cycles, four 8-cycles, and one each of lengths 5, 9 and 28. No cycle of length 3 or 7 is known.
   - **The sharpest repository statement on 3-cycles** is codex's same-day theorem. A primitive 3-cycle made entirely of numbers `2^a p` would need three distinct exponents, the second-smallest at least 129 ([note](aliquot_dyadic_frontier_20261014e.md)).
   - **The structural answer for odd lengths (section 3b).** Dickson (1913) proved that a sociable cycle of odd length `k > 1` cannot be regular, i.e. the gcd of its members cannot be a unitary divisor of each member.
     - Kobayashi–Pollack–Pomerance (2009) reprove this and bound odd-order cycles by `x/(log_2 x)^(1/2+o(1))`. Pollack later improved the bound to `x/(log x)^(1+o(1))` using the sharp gcd form (his Lemma 5).
     - The proof is a sign-frustration argument:
     - if `q | σ(n_i)` for every member, then `n_(i+1) = σ(n_i) − n_i ≡ −n_i (mod q)`, so the aliquot map is a sign flip mod `q` along the cycle;
     - an odd closed walk then forces `n ≡ −n`.
   - **The tournament dictionary (ANALOGY with an exact core).**
     - In tournament and graph language this is Harary balance: an all-negative cycle is balanced iff it is even. Equivalently, an odd cycle is not 2-colourable.
     - Over `Q`, the cycle equations `n_i + n_(i+1) = σ(n_i)` are the unoriented incidence system of the cycle graph `C_k`. For odd `k` (non-bipartite) the solution is unique, with `n_i = ½ Σ_j (−1)^j σ(n_(i+j))`. For even `k` it has a free alternating direction.
   - **Data.**
     - 229 of the 231 amicable pairs, 9 of the 10 four-cycles, and all 6 published 6- and 8-cycles have an odd prime dividing every `σ(n_i)` but not every member, which no odd cycle can have.
     - Both known odd cycles (lengths 5 and 9) have `gcd σ(n_i)` equal to `8` and `2`.
   - **What this does and does not explain.** It explains why odd lengths are much rarer than even ones. It does **not** single out 3 (5 and 9 exist).
   - **Parity and abundancy give no obstruction to 3-cycles.** The parity lemma allows all-even, all-odd, or exactly two members in `□ ∪ 2□`. Abundancy needs only `AAD` or `ADD`.
3. **`{3, 5, 17, 257, 65537}` and groups of five.** (We read the prompt's 27 as 17 = `F_2`; 27 = `3^3` is not a Fermat number.)
   - **Already in the repository ([arithmetic-seams note](arithmetic_seams_20260921_tournaments.md) §5, 2026-09-21; codex's [path-defect packet](tournament_path_defects_20261014e.md)).** The five numbers themselves are Hamiltonian path counts. Reversing the extreme arc of a transitive tournament gives `H(E_n) = 1 + 2^(n−2)`, which is `3, 5, 17, 257, 65537` at `n = 3, 4, 6, 10, 18`.
   - **New, exact (THM-4644 part 1).** The Fermat primes are exactly the primes at which every circulant tournament on `Z_p` is rigid (`Aut = Z_p`). So "infinitely many Fermat primes?" is a tournament question. At a non-Fermat prime, an odd-order multiplier group always gives a symmetric circulant tournament. The Paley primes `p ≡ 3 (mod 4)` are the maximally symmetric case. This is the repository's old "Paley maximally symmetric versus Collatz rigid" contrast, now as a dichotomy of primes.
   - **Exact, classical.** The Paley graphs on the Fermat primes 5 and 17 are the unique extremal graphs for the only two non-trivial diagonal Ramsey numbers known: `P_5 = C_5` for `R(3,3) = 6`, and `P_17` for `R(4,4) = 18`. `P_257` has clique number 7 (computed) and is not extremal.
   - **The Platonic five (S10 verdict, unchanged).** The two fives share only `{3, 5}`. What is exact is the triangle `{3, 4, 5}`: Schläfli numbers = constructible polygons up to 5 = field sizes `q` with `PSL(2, q)` (or `PGL(2, 3)`) the rotation groups.
   - **Collatz (S32).** At a Fermat prime `F_k`, the `F_k x + 1` collector weights have denominator `2^(2^(k+1)) − 1 = F_0 F_1 ⋯ F_k`. The period-32 clock `2^32 − 1` is the product of the five known Fermat primes.
4. **1093, 3511 and `{7, 21}`.**
   - **New, exact (THM-4644 part 2).** Wieferich primes are exactly the primes at which the circulant tournaments on `Z_(p^2)` invariant under the odd part `g` of doubling explode:
     - 1093: `2^6564` instead of `2^12`;
     - 3511: `2^3512` instead of `2^2`;
     - for every other prime below 4000 there is no jump.
     - This is the tournament reading of the doubling-orbit jump already in THM-4631 (`c_2 = p c_1`) and codex's split/grow tower. Codex's same-day path-defect packet gives a different, path-level criterion: `p` is Wieferich iff `p` divides the number `(2^p − 2)/p` of free marked-path orbits of `E_(p+2)`. See section 1b.
   - **Exact.** 3511 is a Paley prime (`ord_3511(2) = 1755 = (p−1)/2`, `p ≡ 3 (mod 4)`), so its doubling-invariant circulant tournaments are the Paley tournament and its reverse, exactly as `<2> = {1, 2, 4}` gives `P_7` at 7.
   - **Exact.** `1093 = (3^7 − 1)/2` is a base-3 repunit of length 7, with `ord_1093(3) = 7`.
   - **Exact, repository (THM-307).** The Burnside count of tournaments mod `p` involves the Euler quotient of 2, which vanishes exactly at Wieferich primes.
   - **The link to the forbidden Hamiltonian-path values `{7, 21}` is NUMEROLOGY** (section 6).
   - **The honest kinship is a second-digit lock (ANALOGY).** Rédei (`H` odd) and Fermat (`2^(p−1) ≡ 1 (mod p)`) are first-digit facts. The gap 21 is a mod-4 lock, and Wieferich primes are where the mod-`p^2` digit vanishes.
5. **`{deficient, perfect, abundant}` and other triples (section 4).**
   - **The exact anchor.** `σ(n)/n = Σ_(d | n) 1/d` for every `n`, so `σ(n)/n − 2 = Σ_(d | n, d > 1) 1/d − 1`. When `τ(n) = 4` there are exactly three divisors `d_1 < d_2 < d_3` above 1. The right side is then the angle excess (over `π`) of the triangle with angles `π/d_i`: deficient, perfect, abundant are hyperbolic, Euclidean, spherical.
   - Only `6 = 2 · 3` is perfect, which is the Euclidean `(2, 3, 6)`. No `n` with `τ(n) = 4` is abundant, so the spherical side is empty.
   - **The `(2, 3, k)` chain.** `k = 3, 4, 5` are spherical (the Platonic rotation groups of orders 12, 24, 60). `k = 6` is Euclidean, `1/2 + 1/3 + 1/6 = 1`, the perfect 6. `k = 7` is the first hyperbolic case, `χ = −1/42`, with `42 = 2 · 3 · 7` the next primary pseudoperfect number; this gives Hurwitz's `84(g − 1)` and Klein's 168.
   - **The other triples are typed ANALOGY** (section 4).
6. **The Parts-graph remark.** Owner's summary of mac-mini's Hadwiger–Nelson type-M note; SOLVER verdicts, no DRAT certificate. The rigid radius-2 hinge is a small certificate of non-4-colourability, and the shared theme is the minimal obstruction (ANALOGY):
   - an odd cycle is the minimal certificate of non-2-colourability; in aliquot form it is Dickson's frustration (section 3b);
   - `H = 7` is impossible because its only conflict graph `K_3` (three pairwise-conflicting odd cycles) never occurs (THM-200/THM-029; all `n` for `{7, 21}` in THM-1370);
   - an opposite-parity amicable pair needs squares on both sides.
7. **The precise hypotheses** (HYP-9273):
   - **(A)** every aliquot cycle avoids `□ ∪ 2□`;
   - **(B)** every Paley tournament `P_q` has `H(P_q)/(q(q−1)/2) ≡ 1 (mod 8)`. Six cases hold; broader versions are refuted (all circulants at `p = 19`; the doubly regular order-15 tournament). The first test is `q = 31`.

## 1. THM-4644: Fermat and Wieferich primes through circulant tournaments (PROVED + FINITE-EXACT)

Statement and proof are in the [canon file](../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md).

**The mechanism.**
- By Burnside's theorem a circulant tournament of prime order has `Aut = Z_p ⋊ M`, where `M` is the multiplier stabilizer, of odd order (`−1` reverses every arc).
- When `(Z/p)^*` is a 2-group, its only involution `−1` lies in every non-trivial subgroup, so `M = 1`. When `p − 1` has an odd factor, a union of cosets of the odd subgroup gives a symmetric tournament.
- On `Z_(p^2)` the odd part `g` of doubling has `ord_(p^2)(g) = d'` exactly at Wieferich primes, instead of `p d'`. So the number of `g`-orbits on the units, and with it the number of `g`-invariant tournaments, grows by a factor `p` in the exponent.
- The same orbit count gives THM-4631's Wieferich jump in fast-mode multiplicity: 1093 has 3 modes at level 1 and 3279 at level 2; 3511 has 2 and 7022.

**Data** (`.out` §1–2):
- **Brute force.** `|Aut| = p` for every circulant tournament at `p = 3, 5`; `{7, 21}` at 7; `{11, 55}` at 11.
- **Stabilizers.** Trivial for all circulants exactly at `p = 3, 5, 17` among `p <= 31`.
- **The exponent `E`**, against the value `2^((p−1)/d')` it would take for a non-Wieferich prime:

| `p` | `ord_p(2)` | `d'` | `ord_(p^2)(g)` | `E` | non-Wieferich value |
|---|---|---|---|---|---|
| 3 | 2 | 1 | 3 | 2 | 2 |
| 7 | 3 | 3 | 21 | 2 | 2 |
| 31 | 5 | 5 | 155 | 6 | 6 |
| 127 | 7 | 7 | 889 | 18 | 18 |
| **1093** | 364 | 91 | **91** | **6564** | 12 |
| **3511** | 1755 | 1755 | **1755** | **3512** | 2 |

## 1b. Relation to the same-day codex packets (CITED, repository)

Codex worked on the same owner themes on the same day. To avoid double counting:

| topic | codex packet | this note |
|---|---|---|
| Fermat numbers in tournaments | `H(E_n) = 1 + 2^(n−2)` (from the arithmetic-seams note): the five Fermat numbers as Hamiltonian path counts | THM-4644 part 1: Fermat **primes** as the primes where every circulant tournament is rigid (the converse direction is new; the forward direction is THM-871) |
| Wieferich primes | marked-path orbit count `(2^p − 2)/p` divisible by `p`; split/grow tower `C_k = ((p−1)/d) p^(min(k−1, s−1))`; Fermat primes are never Wieferich (classical, Crandall–Dilcher–Pomerance) | THM-4644 part 2: the same orbit jump counted as `g`-invariant circulant tournaments on `Z_(p^2)`, including the non-unit orbits; part 3: 3511 is a Paley prime |
| aliquot 3-cycles | PROVED: dyadic-times-prime 3-cycles need three distinct exponents, second-smallest `>= 129`; affine core compiler; the `I + P` distinction; Pollack's gcd restriction replayed on the 5-cycle | the parity lemma and HYP-9273 (A); the census statistics of section 3b over all 255 cycles; the signed-graph (balance) reading |
| abundancy | `Σ σ = 2 Σ n`, `Π σ/n > 2^L`, the minimum node abundant and the maximum deficient | the same identities (found independently), plus the `τ(n) = 4` triangle reading |

## 2. Fermat primes and the groups of five (PROVED + CITED + NUMEROLOGY)

* **The tournament meaning of "Fermat prime" (THM-4644 part 1).**
  * Rotational tournaments at `F_k` up to isomorphism: the multiplier group `Z_(p−1)` acts freely, giving `2^((p−1)/2)/(p − 1)` classes, namely `1, 1, 16, 2^120, 2^32752`. THM-871 computed 16 at `p = 17`.
  * The Fermat and Paley primes are the two ends of a symmetry spectrum; the Collatz graph (`Aut(N, T) = 1`, THM-4523) is rigid like the Fermat rungs.
* **Ramsey.** Clique numbers by exhaustive branch and bound; the graphs are self-complementary, so `α = ω`.

  | graph | `ω` | bound obtained |
  |---|---|---|
  | `P_5` | 2 | `R(3,3) > 5` |
  | `P_13` | 3 | `R(4,4) > 13` |
  | `P_17` | 3 | `R(4,4) > 17` |
  | `P_29` | 4 | `R(5,5) > 29` |
  | `P_37` | 4 | `R(5,5) > 37` |
  | `P_41` | 5 | `R(6,6) > 41` |
  | `P_257` | 7 | `R(8,8) > 257` |

  * `P_5` and `P_17` are the unique extremal graphs for `R(3,3) = 6` and `R(4,4) = 18` (CITED).
  * That these are Fermat primes is NUMEROLOGY: `43 <= R(5,5) <= 46` (upper bound Angeltveit–McKay 2024), nowhere near 257.
* **The Platonic five.** S10's verdict stands. The exact common structure is `{3, 4, 5}` and `PSL(2,3) = A_4`, `PGL(2,3) = S_4`, `PSL(2,5) = PSL(2,4) = A_5`.
* **The `F_k x + 1` collector** (S32 Proposition 2): `d = ord_(F_k)(2) = 2^(k+1)`, so the weights `2^(d−m)/(2^d − 1)` have denominator `F_0 F_1 ⋯ F_k`. EXACT and trivial.
* **About 27.** If 27 was intended, it is `3^3`, the Collatz record seed with 111 steps. It is not a Fermat number; NUMEROLOGY.

## 3. Aliquot parity, amicable pairs and sociable cycles (PROVED + FINITE-EXACT + CITED)

**Lemma (classical).** `σ(n)` is odd iff the odd part of `n` is a perfect square, i.e. iff `n ∈ □ ∪ 2□`.

*Proof.* `σ(2^e) = 2^(e+1) − 1` is odd, and for odd `m`, `σ(m) ≡ τ(m) (mod 2)`. Divisors of `m` pair as `d ↔ m/d` except `√m`. ∎

**Corollaries (PROVED).**
1. **Parity of the aliquot sum.** `s(n) = σ(n) − n ≡ n + [n ∈ □ ∪ 2□] (mod 2)`: the aliquot map preserves parity except at squares and twice-squares.
2. **Even count.** Along any aliquot cycle the parity changes an even number of times, so the number of members in `□ ∪ 2□` is even.
3. **Amicable pairs.** For an amicable pair `(m, n)`, `m + n = σ(m)` is odd iff `σ(m)` is odd. So the pair has opposite parity iff both members are in `□ ∪ 2□`, iff one is (Gioia–Vaidya 1967).
4. **Arithmetic mean.** For a `k`-cycle, `Σ σ(n_i) = Σ (n_i + n_(i+1)) = 2 Σ n_i`.
5. **Geometric mean.** With `r_i = n_(i+1)/n_i` (so `Π r_i = 1`), `Π σ(n_i)/n_i = Π (1 + r_i) >= Π 2√r_i = 2^k`, with equality iff all `r_i = 1` (a perfect number). So every non-trivial cycle has an abundant member and a deficient one.

**Census** (`aliquot_parity_squares_20261014.out`; exact σ-sieve to `10^8`; cycles found exactly when all members are `<= 10^8`):
- **248 cycles:** 5 perfect, 231 amicable pairs (54 of all cycles are all-odd), ten 4-cycles, Poulet's 5-cycle (12496) and 28-cycle (14316). There is no 3-cycle.
- **Parity and squares.** Every cycle is parity-homogeneous; one 4-cycle is all odd: `28158165 → 29902635 → 30853845 → 29971755`. No cycle contains a square or a twice-square. Corollary 4 holds.
- **Abundancy patterns:** `AADD` (8 up to rotation), `ADDD` (1), `ADAD` (2) for the 4-cycles; `AADDD` for the 5-cycle.

**Double flips** (odd `a <= 10^9` with `s(a^2) ∈ □ ∪ 2□`). There is no int64 overflow: for odd `a <= 10^9`, `σ(a^2)/a^2 < Π_(3 <= p <= 23) p/(p−1) < 3.06`, so `σ(a^2) < 3.1·10^18 < 2^63`.

| `a <=` | `10^1` | `10^2` | `10^3` | `10^4` | `10^5` | `10^6` | `10^7` | `10^8` | `10^9` |
|---|---|---|---|---|---|---|---|---|---|
| double flips | 2 | 7 | 17 | 39 | 95 | 242 | 607 | 1634 | 4476 |
| composite | 0 | 1 | 2 | 3 | 6 | 13 | 15 | 17 | 18 |
| `x/(log x)^1.5` | 2.9 | 10 | 55 | 358 | 2560 | 19500 | 155000 | 1.26·10^6 | 1.06·10^7 |

- **The prime ones** are `a = 2b^2 − 1` (and `a = 3`), with `s(a^2) = a + 1`.
- **The composite ones** are `49` (`s = 20^2`), `287 = 7·41`, `4081 = 7·11·53`, and mostly semiprimes `pq`. 11 of the 18 have a prime factor `p` with `p + 1 = 2b^2` (`p = 7, 71, 127, 199, 449, 647`), i.e. a prime double flip; `71` occurs three times. The largest ratio `s(a^2)/a^2 = 0.8299` occurs at `a = 949818597 = 3·7^2·13·601·827`.
- **For fixed `p`,** the condition on `q` in `s(p^2 q^2) = (p+1)q^2 + (p^2+p+1)(q+1) ∈ □ ∪ 2□` is a conic in `q`. When `p + 1` lies in the same class as the target (square, or twice a square), the conic degenerates and has finitely many solutions. For `p = 71` (`72 = 2·6^2`), the twice-square case `y^2 − (6q)^2 = 5113(q+1)/2` with `t = y − 6q` gives `q = (2t^2 − 5113)/(5113 − 24t)` for `51 <= t <= 213`. The integer values are `q = 409, 3391` (prime) and `85625` (not prime), checked. Otherwise it is a Pell conic, with at most `O(log x)` solutions up to `x`.
- **Result.** No double flip satisfies `σ(s(a^2)) = σ(a^2)`.

**Status.**
- The opposite-parity amicable problem is OPEN; every known amicable pair has equal parity.
- Sociable 3-cycles are unknown. Bratley–Lunnon–McKay searched to `5·10^7` (1970). Moews and Moews searched exhaustively for cycles of length greater than 2 whose member preceding the largest one is below `1.54·10^11`; Pedersen extended this to `4.27·10^11` (1997).
- HYP-9273 (A) is the parity-homogeneity statement.

## 3b. Odd aliquot cycles are frustrated (PROVED + CITED + FINITE-EXACT)

**Lemma (frustration; the mechanism of Dickson 1913, KPP 2009 Theorem E; Pollack's Lemma 5 is the sharp form).** Let `n_0 → ⋯ → n_(k−1) → n_0` be an aliquot cycle, `G = gcd(n_i)` and `Q = gcd(σ(n_i))`.
- If `k` is odd, then `odd(Q) | G` and `Q | 2G`. Pollack shows `Q | G` unless `2 ∥ Q`.
- If `k` is even, then `n_i ≡ (−1)^i n_0 (mod Q)`.

*Proof.* `n_(i+1) = σ(n_i) − n_i ≡ −n_i (mod Q)`. Around the cycle, `n_0 ≡ (−1)^k n_0`. For odd `k` this gives `Q | 2n_0`, and the same for every `i`. ∎

**Reading.**
- **Linear algebra.** The cycle equations `n_i + n_(i+1) = σ(n_i)` are `B^T n = σ` for the unoriented incidence matrix `B` of the cycle graph `C_k`. `det(I + P) = 1 − (−1)^k`, where `P` is the cyclic shift.
  - **Odd `k`** (`C_k` not bipartite): the system is uniquely solvable over `Q`, with `n_i = ½ Σ_(j=0)^(k−1) (−1)^j σ(n_(i+j))`, and invertible mod every odd `q`.
  - **Even `k`:** it is solvable iff the alternating sum of the `σ(n_i)` vanishes, which it always does on a cycle, with the free direction `(1, −1, 1, …)`.
  - Codex's [guarded-core packet](aliquot_guarded_cycles_20261014e.md) records the same `I + P` distinction, and notes that it treats lengths 3 and 5 alike.
- **Signed graphs.** With every cycle edge negative (`n ↦ −n` mod `Q`), Harary's balance theorem says the cycle is balanced iff `k` is even. So the aliquot map, wherever `σ` is divisible by `q`, is a 2-colouring constraint, and odd cycles are its obstructions.
- **Dickson's theorem.** An odd cycle (`k > 1`) whose gcd `a > 1` is unitary in every member would force `σ(a)/a ∈ {1, 2}`. KPP use this together with "`σ(n)` is divisible by every small prime power for almost all `n`" to show that irregular cycles, hence all odd-order cycles, are rare.

**Census** (`aliquot_dickson_frustration_20261014.out`; all 248 cycles with members `<= 10^8`, plus seven published cycles recomputed and closed: 6-cycles from `21548919483`, `90632826380`, `1771417411016`, `4773123705616`; 8-cycles from `1095447416`, `1276254780`; the 9-cycle from `805984760`). Every cycle satisfies the lemma.

| length | cycles | regular (Cohen) | an odd prime divides every `σ(n_i)` but not every member |
|---|---|---|---|
| 2 | 231 | 189 | 229 |
| 4 | 10 | 8 | 9 |
| 6 | 4 | 3 | 4 |
| 8 | 2 | 2 | 2 |
| 28 | 1 | 0 | 0 |
| **5** | 1 | **0** (Dickson) | **0** (forced): `G = 8`, `Q = 8` |
| **9** | 1 | **0** (Dickson) | **0** (forced): `G = 2`, `Q = 2` |

**What the table shows.**
- Even cycles almost always use the sign flip modulo some odd `q`. For amicable pairs the median odd excess `odd(Q)/gcd(odd(Q), G)` is 13041, and 94% have `3` in it.
- An odd cycle cannot do this, so it must avoid every common odd divisor of the `σ(n_i)`.
- Members of 2- and 4-cycles are "σ-typical": `3 | σ(n)` for 98.5% of them against 85.1% of all `n <= 10^8`, and `9 | σ(n)` for 93.6% against 57.4%.
- This is why odd lengths are rare. It does not single out 3: lengths 5 and 9 exist, and the long even 28-cycle also has no odd excess.

## 4. Triples (EXACT anchors, ANALOGY elsewhere)

| `{deficient, perfect, abundant}` | the other triple | exact anchor |
|---|---|---|
| `σ(n)/n` below, at, above 2 | triangle `(π/d_1, π/d_2, π/d_3)`: hyperbolic (defect), Euclidean, spherical (excess) | `τ(n) = 4`: `σ(n)/n − 2 = Σ_(d > 1) 1/d − 1` (a tautology; the content is the reading); only 6 is perfect (↔ `(2,3,6)`); none is abundant below `10^5` or at all (`1/p + 1/q + 1/pq <= 1`, `1/p + 1/p^2 + 1/p^3 < 1`) |
| — | `(2,3,k)`: spherical (Platonic) `k <= 5`, Euclidean 6, hyperbolic from 7 | `χ(2,3,6) = 0` with 6 perfect; `χ(2,3,7) = −1/42` with 42 primary pseudoperfect; primary pseudoperfect numbers `2, 6, 42, 1806, 47058` |
| aliquot cycles | perfect on average | `Σ σ(n_i) = 2 Σ n_i`; geometric mean of abundancy `>= 2` |
| Collatz words | `2^L > 3^l` (positive cycles), `=` (impossible), `<` (negative cycles) | S32 carry sign; the four cycles are near-perfect words |
| `px + 1` maps | contracting `p < 4`, critical `p = 4` (absent: `p` odd), expanding `p > 4` | the drift `log p − 2 log 2` |
| tournaments | vertex score below, at, above `(n−1)/2` | regular = every vertex "perfect" |
| signed cycles | balanced (even) / frustrated (odd) | section 3b (a dichotomy, not a triple) |

**NUMEROLOGY (law of small numbers).** Sylvester's sequence `2, 3, 7, 43, 1807` (`e_(n+1) = e_n^2 − e_n + 1`) gives the Egyptian fractions of 1 behind the pseudoperfect numbers. It contains the Heegner numbers 7 and 43 (S17). Also `1807 = 13 · 139`, where 139 is the Collatz denominator `3^7 − 2^11` of the `−17` cycle.

## 5. Paley Hamiltonian paths and HYP-9273 (B) (FINITE-EXACT + OPEN)

| `q` | `H(P_q)` | `q(q−1)/2` | orbits `H/|Aff|` | mod 8 |
|---|---|---|---|---|
| 3 | 3 | 3 | 1 | 1 |
| 7 | 189 | 21 | 9 | 1 |
| 11 | 95095 | 55 | 1729 | 1 |
| 19 | 1172695746915 | 171 | 6857869865 | 1 |
| 23 | 15760206976379349 | 253 | 62293308207033 | 1 |
| 27 | (mod `8 · 1053`) | 351 (affine) | — | 1 |

**`q = 27`.** Under the full `Aut` (order 1053, with Frobenius) the orbit count is `3 (mod 8)`; under the affine group it is 1.

**Controls** (`.out` §3):
- **All circulant tournaments on `Z_p`**, residues of `H/|Aut|` mod 8:

  | `p` | `p mod 4` | residues |
  |---|---|---|
  | 3, 7, 11 | 3 | all 1 |
  | 19 | 3 | `{1: 128, 3: 126, 5: 96, 7: 162}` |
  | 5 | 1 | all 3 |
  | 13 | 1 | mixed |
  | 17 | 1 | mixed |

- **The doubly regular `T_4`** (order 15, `|Aut| = 21`): orbit count `9444525 ≡ 5 (mod 8)`.

**Conclusion.** The congruence is not a circulant phenomenon (refuted at 19) and not a doubly regular one (refuted at 15). It is Paley-specific among the cases computed. The orbit count's oddness is Rédei's theorem plus the odd order of the affine group; the hypothesis concerns the next two binary digits. The orbit count `1729` at `q = 11` is the taxicab number (NUMEROLOGY).

## 6. Wieferich, {7, 21}, and second-digit locks (ANALOGY + NUMEROLOGY)

**Exact facts:**
- THM-4644 part 2.
- THM-4631's multiplicity jump.
- THM-307's Euler quotient.
- `1093 = R_7(3)` with `ord_1093(3) = 7`.
- 3511 is a Paley prime.

**The second-digit-lock analogy (ANALOGY).**
- On the tournament side, the OCF `H = 1 + 2α_1 + 4α_2 + ⋯` reads `H` in base 2. Rédei is the first digit; the gaps lock at the second (21 at `n = 6` is a mod-4 lock).
- On the number side, `2^(p−1) = 1 + p·q_p(2) + ⋯` reads `2^(p−1)` in base `p`. Fermat is the first digit; Wieferich primes are where the second digit vanishes.
- We know no map between the two.

**NUMEROLOGY:**
- `{7, 21} = {|V(P_7)|, |Aut(P_7)|}`.
- Both Wieferich primes are `≡ 1 (mod 3)` and quadratic residues mod 7 (`1093 ≡ 1`, `3511 ≡ 4 (mod 21)`). Each condition holds for half the primes, so the joint coincidence has prior probability `1/16`.
- `13 | p − 1` for both, and 13 is a Wilson prime.
- `364 = C(14, 3)`.
- The `H = 7` obstruction is a triangle `K_3` of conflicting 3-cycles, while aliquot 3-cycles are unknown: a rhyme on the number 3, with no mechanism.

## 7. Hostile controls

1. **Amicable parity.** Not proved: the problem is open. The reduction is classical (Gioia–Vaidya). Our search range is inside Chernykh's complete enumeration, so it is a check, not a new bound.
2. **Dickson's frustration** explains the rarity of odd lengths, not the absence of 3. Since 5 and 9 occur, no argument of this kind can forbid 3 by parity of length alone. Codex makes the same point for the linear `I + P` level. The frustration lemma itself is classical (Dickson; Pollack's Lemma 5); what is new here is the census statistics and the balance reading.
3. **HYP-9273 (A)** has no proof and no mechanism beyond the thin-set heuristic. **(B)** rests on six cases; its two broader versions were refuted in this session.
4. **THM-4644 is elementary.** Its value is the reformulation. Neither part solves the underlying problem (infinitude of Fermat or of Wieferich primes).
5. **The Parts remark** is the owner's summary of mac-mini's SOLVER result, not re-run here.

## 8. Open threads

1. **HYP-9273 (B) at `q = 31`.** The direct dynamic programme needs about 33 GB. Alternatives: count orbit representatives using the free affine action, or compute the OCF digits `α_1 mod 4` and `α_2 mod 2` for Paley tournaments by character sums.
2. **Composite double flips.**
   - The data (18 up to `10^9`) and the conic structure suggest polylogarithmic growth, which would make an enumeration of all odd-member candidates far beyond `10^20` feasible via Pell recursions.
   - A completeness argument is needed: every composite `a` has some prime `p | a`, and the conic in the cofactor must be handled uniformly.
3. **A finer frustration for `k = 3`.** Dickson's argument uses one modulus at a time. A 3-cycle needs, for each small odd prime `q`, a member whose `σ` is prime to `q` (or all members divisible by `q`). Joint constraints (mod 4 and mod 9 together with the abundancy pattern `AAD`/`ADD`) might give a quantitative heuristic for the 3-cycle count versus the 5-cycle count.
4. **Tournament reformulations of other prime classes** (Mersenne, Wilson, Wall–Sun–Sun) in the style of THM-4644.

## 9. Audit

Pending: an independent adversarial audit will be recorded here.

## 10. Reproduction and references

**Reproduction:**
- `python3 04-computation/experiments/number_tournament_bridges_20261014.py` (about 3 minutes, 4 GB);
- `python3 04-computation/experiments/aliquot_parity_squares_20261014.py 1e8 1e9` (about 2 minutes, 5 GB);
- `python3 04-computation/experiments/aliquot_dickson_frustration_20261014.py` (about 1 minute, 2 GB; reads the JSON of the previous script).

**References.**
- **Groups and tournaments:**
  - W. Burnside, *Proc. LMS* (1901), transitive groups of prime degree;
  - B. Alspach (1970), via THM-871;
  - L. Rédei (1934); J. W. Moon, *Topics on Tournaments* (1968).
- **Wieferich primes:**
  - A. Wieferich (1909); W. Meissner (1913); N. G. W. H. Beeger (1922);
  - OEIS A001220 (PrimeGrid completed the search to `2^64` in December 2022).
- **Ramsey:**
  - R. E. Greenwood and A. M. Gleason, *Canad. J. Math.* 7 (1955);
  - R. J. Evans, J. R. Pulham and J. Sheehan, *J. Combin. Theory B* 30 (1981);
  - V. Angeltveit and B. D. McKay, `R(5,5) <= 46` (2024).
- **Amicable pairs:**
  - A. A. Gioia and A. M. Vaidya, Amicable numbers with opposite parity, *Amer. Math. Monthly* 74 (1967) 969–973, correction 75 (1968) 386;
  - G. D'Abramo, arXiv:math/0501402 (2005);
  - P. Hagis, Lower bounds for relatively prime amicable numbers of opposite parity, *Math. Comp.* 24 (1970);
  - D. E. Iannucci and F. Luca, Squares whose aliquot sum is a square, *JP J. Algebra Number Theory Appl.* 13 (2009) 17–26;
  - P. Pollack, Powerful amicable numbers, *Colloq. Math.* 122 (2011) 103–123;
  - S. Chernykh's list (OEIS A066873: 4197267 pairs with smaller member below `10^20`).
- **Sociable numbers:**
  - L. E. Dickson, Theorems and tables on the sum of the divisors of a number, *Q. J. Math.* 44 (1913) 264–296;
  - M. Kobayashi, P. Pollack and C. Pomerance, On the distribution of sociable numbers, *J. Number Theory* 129 (2009) 1990–2009;
  - P. Pollack, A remark on sociable numbers of odd order (pollack-math.net/socodd.pdf): Lemma 5, and the bound `x/(log x)^(1+o(1))` for fixed odd order;
  - P. Poulet (1918);
  - P. Bratley, F. Lunnon and J. McKay, *Math. Comp.* 24 (1970);
  - D. and P. C. Moews; J. O. M. Pedersen (search to `4.27·10^11`, 1997);
  - Wikipedia, "Sociable number" (table as of October 2024).
- **Hurwitz and Klein:** A. Hurwitz (1893); F. Klein (1879).
- **This repository:**
  - THM-871; THM-307; THM-200, THM-029, THM-1370; THM-4523; THM-4631;
  - the arithmetic-seams note ([tournaments](arithmetic_seams_20260921_tournaments.md) §5, the Fermat family `H(E_n) = 1 + 2^(n−2)`);
  - codex 2026-10-14: [dyadic three-cycle frontier](aliquot_dyadic_frontier_20261014e.md), [guarded core compiler](aliquot_guarded_cycles_20261014e.md), [tournament path defects](tournament_path_defects_20261014e.md), [Fermat/Wieferich order lifts](fermat_wieferich_order_20261014e.md);
  - THM-4641 ([note](repair_algebra_minus_one_and_wilson_dictionary_20261014.md));
  - the S10 Fermat/Platonic verdict ([note](gilbreath_fermat_platonic_20260926.md));
  - the camion/Busch gap note ([note](camion_busch_gaps_polyhedra_collatz_20261002.md));
  - mac-mini's Hadwiger–Nelson type-M note ([note](hn_type_m_hinge_rotations_and_forced_radii_20261014.md));
  - the S17 and S32 notes.
