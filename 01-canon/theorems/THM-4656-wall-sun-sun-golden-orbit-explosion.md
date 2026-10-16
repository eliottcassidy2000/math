---
id: THM-4656
title: "Wall-Sun-Sun primes in the circulant-tournament dictionary of THM-4644: the golden (metallic) unit is a Wieferich base exactly at (metallic) Wall-Sun-Sun primes (KNOWN), so the unit-invariant Cayley tournaments on Z_(p^2) (split primes) and on Z[mu]/(p^2) (inert primes) explode exactly there (verified at the Pell primes 13, 31); golden rigidity holds exactly when the Pisano period is a power of 2, i.e. p | L_(2^j) (OEIS A074714), including every Mersenne prime 2^q - 1 with q = 2 or q = 3 (mod 4) (the necessity half of Lucas-Lehmer with seed 3); golden Paley primes = Fibonacci-primitive-root primes = 3 (mod 4)"
status: >
  The tournament orbit-explosion reading (parts 1-2, assembled from THM-4644's orbit argument with the metallic unit mu_a in
  place of 2) is PROVED and was not found in the literature (audit E). KNOWN inputs, cited: part 1's unit criterion (p is an
  a-Wall-Sun-Sun prime iff mu_a is a Wieferich unit at p: Bragman-Rowland 2022, arXiv:2202.00704, Thm 10; McIntosh-Roettger
  2007); Wall 1960 (pi(p^2) in {pi(p), p pi(p)}); pi(p^2) = pi(p) iff p^2 | F_(p-(5/p)) iff L_p = 1 mod p^2 (Andrejic 2006;
  McIntosh-Roettger 2007; Sun-Sun 1992 gave the Fermat-Last-Theorem connection); no Wall-Sun-Sun prime below 2^64
  (PrimeGrid 2020-2022); the invariant-subset count (Akin 2018, arXiv:1806.11241, as in THM-4644); the part 3 list (OEIS
  A074714); Fibonacci primitive roots (Shanks 1972, OEIS A003147); Lucas-Lehmer; pi(n) <= 6n (Freyd-Brown 1992, AMM E3410).
  Part 3's equivalences and the Mersenne statement are PROVED (elementary). FINITE-EXACT: brute-force orbit counts for 27
  cases (a = 1, 2, 3), including the metallic Wall-Sun-Sun primes p = 13, 31 (a = 2), where the jump is observed; audit E's
  independent Burnside counts (76 cases plus p = 191, 241, 643, 1546463); the golden-rigid list for p < 10^6. Audited
  (audit E, MISTAKE-615).
  Extension 2026-10-15: the formerly FINITE-EXACT full split-ring remark is PROVED at every p^k, together with the inert
  all-depth count, by Theorem U1 in unit_lift_clock_20261015a.md; 61 full-ring/Burnside checks cover 2,442,949 elements.
source: mac-mini-2026-10-14-tourn
depends_on: [THM-4644 (Fermat and Wieferich primes as circulant-tournament phenomena, opus S33), Wall 1960, Bragman-Rowland 2022, McIntosh-Roettger 2007]
verification:
  - 04-computation/experiments/wall_sun_sun_20261014/wss_orbits.py (+ .out)
  - 04-computation/experiments/wall_sun_sun_20261014/golden_fermat_check.py (+ .out)
  - 04-computation/experiments/audit_E_20261014/wss_audit.py, golden_rigidity_audit.gp, golden_paley_fpr.py, wss_fullring_split.py (+ .out)
  - 04-computation/experiments/unit_lift_clock_20261015a.py (+ 05-knowledge/results/unit_lift_clock_20261015a.json)
note: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
related: [THM-4644, THM-4631, THM-4641 (Wilson / {2,4,8} / Fermat dictionary), HYP-9278 (S35: rigidity circuits), THM-4657 (opus S36, ceded to this theorem; independent confirmation at 40 split pairs)]
---

# THM-4656: Wall–Sun–Sun primes in the circulant-tournament dictionary

## Setting

- For an integer `a >= 1`, let `μ = μ_a` be a root of `x² − a x − 1`. So `μ₁ = φ` is the golden ratio and `μ₂ = 1 + √2` is the silver ratio. Put `D = a² + 4` and `R_k = Z[μ]/(p^k)`.
- For a prime `p ∤ 2D`, let `π_a(p)` be the order of `μ` in `R_1^*`. This is the period of `U_n(a, −1)` mod `p`, and `π_1` is the Pisano period.
  - If `p` splits, i.e. `(D/p) = +1`, then `π_a(p)` divides `p − 1`.
  - If `p` is inert, then `π_a(p)` divides `2(p + 1)`.
- `π_a(p²)` is either `π_a(p)` or `p·π_a(p)` (Wall). The prime `p` is an **`a`-Wall–Sun–Sun prime** when `π_a(p²) = π_a(p)`. For `a = 1` these are the Wall–Sun–Sun primes, and there is none below `2^64`.
- A Cayley tournament on an abelian group `G` of odd order is a set `S ⊂ G ∖ {0}` with `S ∩ (−S) = ∅` and `S ∪ (−S) = G ∖ {0}`.

## Statement

1. **(Split primes: Wieferich units; KNOWN criterion, new reading.)** Let `p` split and let `μ₁ ∈ Z/p²` be a root of `x² − a x − 1` (Hensel lift).
   - **KNOWN** (Bragman–Rowland 2022, Thm 10; McIntosh–Roettger 2007): `p` is an `a`-Wall–Sun–Sun prime iff `μ₁^(ord_p μ₁) ≡ 1 (mod p²)`, i.e. iff `μ₁` is a Wieferich base at `p`.
   - **Assembly with THM-4644 part 2.** THM-4644 part 2 therefore holds verbatim with the base 2 replaced by `μ₁`. Put `d = ord_p μ₁`, `d'` its odd part and `g = μ₁^(d/d')`. The circulant tournaments on `Z_(p²)` with `gS = S` number `2^E`, where

     `E = (p − 1)/d'` normally, and `E = (p + 1)(p − 1)/(2d')` at `a`-Wall–Sun–Sun primes.
2. **(Inert primes: orbit explosion on `Z[μ]/(p²)`.)** Let `p` be inert, `π = π_a(p)`, `π'` its odd part and `g = μ^(π/π')`. The Cayley tournaments on `(R_2, +) ≅ (Z/p²)²` with `gS = S` number `2^E`, where

   `E = (p + 1)(p² − 1)/(2π')` normally, and `E = (p² + 1)(p² − 1)/(2π')` at `a`-Wall–Sun–Sun primes.

   So the exponent jumps by the factor `(p² + 1)/(p + 1) ≈ p` exactly at Wall–Sun–Sun-type primes. This mirrors THM-4644's factor `(p + 1)/2` at Wieferich primes.

   *Full-ring extension (PROVED, 2026-10-15).* For split primes the same orbit count `(p² − 1)(p^(2−ε) + 1)/π'` holds on the full ring `Z[μ]/(p²) ≅ (Z/p²)²`, with no orbit fixed by negation. This was previously only FINITE-EXACT here (audit E: 81 split cases with `p < 200`, and split exceptions 31, 191, 643, 1546463). [The all-depth proof, Theorem U1](../../05-knowledge/results/unit_lift_clock_20261015a.md), supplies the missing zero-divisor argument: the odd-order generator has norm one, so its two split components are reciprocal and every `g^d-1` has equal component valuations. The same proof gives all `p^k` orbit counts in both split and inert rings, with lift-delay coordinate `h=v_p(mu^pi-1)`. Its independent full-ring/Burnside census covers 61 cases and 2,442,949 elements.
3. **(Golden rigidity.)** Let `p ∤ 10`. The following are equivalent:
   - every subgroup `1 ≠ H <= ⟨φ⟩ ⊂ R_1^*` contains `−1`, so no Cayley tournament on `(R_1, +)` is invariant under a nontrivial power of `φ`;
   - `π(p)` is a power of 2;
   - `p` divides a Lucas number `L_(2^j)`, and then `π(p) = 2^(j+2)`.

   These primes start `3, 7, 47, 127, 1087, 2207, 4481, 21503, 34303, 119809, 524287` (`p < 10^6`). The list is KNOWN (OEIS A074714).

   **Mersenne primes.** Every Mersenne prime `M_q = 2^q − 1` with `q = 2` or `q ≡ 3 (mod 4)` is golden-rigid, with the maximal period `π(M_q) = 2^(q+1)`.
   - Such `M_q` is inert, since `M_q ≡ 2, 3 (mod 5)`.
   - For inert `p`, `φ^(p+1) = N(φ) = −1`. With `p + 1 = 2^q` this gives `φ^(2^q) = −1`.
   - This is the necessity half of the Lucas–Lehmer test with seed 3.
   - For `q ≡ 1 (mod 4)`, `M_q` splits and is not golden-rigid.

   **The Fermat mirror.**
   - Golden rigidity (`ord φ` a power of 2) is the golden mirror of the prime divisors of Fermat numbers (`p | 2^(2^j) + 1`, i.e. `ord_p 2` a power of 2).
   - The mirror of the Fermat primes themselves (`p − 1 = 2^k`, so the split torus `F_p^*` is a 2-group) is the inert primes with `p + 1 = 2^q` (the nonsplit torus containing `φ` is a 2-group). Those are the Mersenne primes above.
4. **(Golden Paley primes.)** Let `p ≡ 3 (mod 4)` split in `Q(√5)`. Exactly one root `φ₁` of `x² − x − 1` is a quadratic residue, since the product of the roots is `−1`.
   - If `φ₁` generates `QR_p`, the circulant tournaments on `Z_p` invariant under `φ₁` are exactly the Paley tournament `P_p` and its reverse (THM-4644 part 3 with base `φ₁`).
   - `φ₁` generates `QR_p` iff the other root `φ₂ = −φ₁^(−1)` is a primitive root. So these are exactly the primes `≡ 3 (mod 4)` with a Fibonacci primitive root (KNOWN: Shanks 1972; OEIS A003147).
   - Below 2000 there are 52 of them, starting `11, 19, 31, 59, 71, 79, 131, …`.
   - At `p = 11`, `φ₁ = 4` has order 5 and `⟨4⟩ = {1, 3, 4, 5, 9} = QR_11`. So `P_11` is the golden tournament on `Z_11`.
   - Its 55 arcs equal `|Aut P_11| = 11·5`. For every Paley tournament of prime order `p`, `Aut(P_p) = {x ↦ ax + b : a ∈ QR_p}` acts regularly on the `C(p, 2)` arcs.
   - That this 55 is also `F_10 = F_(π(11))` is a COINCIDENCE: by Luo's theorem, 55 is the largest triangular Fibonacci number.

## Proof

**1.**
- The equivalence is the cited criterion. Directly: with `μ₂ = −μ₁^(−1)`, `π_a(p^k) = lcm(ord_(p^k) μ₁, ord_(p^k) μ₂)`. Since `−1` has order prime to `p`, `μ₂` lifts with the same `p`-power as `μ₁`. Hence `π_a(p²) = π_a(p)` iff `ord_(p²) μ₁ = ord_p μ₁`.
- THM-4644's proof of part 2 uses only `ord_(p²)(base) ∈ {d, pd}` and the odd order of `⟨g⟩`. So it applies verbatim.

**2.**
- `R_2` is local, with maximal ideal `pR_2` and `|R_2| = p⁴`.
- `⟨g⟩` has odd order `π' p^ε`, where `ε = 0` iff the prime is `a`-Wall–Sun–Sun.
- It acts freely on the `p²(p² − 1)` units, giving `(p² − 1)p^(2−ε)/π'` orbits.
- On `pR_2 ∖ {0} ≅ F_(p²)^*` it acts through `g mod p`, of order `π'`, giving `(p² − 1)/π'` orbits.
- Negation fixes no orbit, since `−1 ∉ ⟨g⟩`. So the count is `2^(N/2)`, with `N` the total number of orbits.

**3.**
- `⟨φ⟩` is cyclic of order `π(p)`.
- **If `π(p) = 2^k`:** every nontrivial subgroup contains the unique involution `φ^(2^(k−1))`, and that involution is `−1`.
  - Split case: `φ^(π/2)` is a pair `(x, y)` with `x, y ∈ {±1}`. Its components multiply to `N(φ)^(π/2) = (−1)^(π/2) = 1`, because `4 | π` (`π = 2^k` with `k >= 2`). So `x = y = −1`.
  - Inert case: `F_(p²)` is a field, whose only involution is `−1`.
- **If `π(p)` has an odd prime factor `r`:** the subgroup of order `r` avoids `−1`, and unions of its cosets give invariant tournaments, as in THM-4644 part 1.
- **Lucas form:** `L_n = φ^n + φ^(−n)` for even `n`. So `p | L_(2^j)` iff `φ^(2^(j+1)) = −1` iff `ord(φ) = 2^(j+2)`.
- The Mersenne case is shown in the statement.

**4.**
- The product of the two roots is `−1`, a nonresidue since `p ≡ 3 (mod 4)`. So exactly one root is a residue.
- If that root generates `QR_p`, of odd order `(p−1)/2`, its cosets are `QR_p` and `−QR_p`.
- **Fibonacci primitive root.** Suppose `φ₁` has odd order `(p−1)/2`. Then `φ₂ = −φ₁^(−1)` has order `lcm(2, (p−1)/2) = p − 1`. Conversely, if `φ₂` is a primitive root, `φ₁ = −φ₂^(−1) = φ₂^((p−1)/2 − 1)` has order `(p−1)/gcd(p−1, (p−1)/2 − 1) = (p−1)/2`.

**Checks** (`wss_orbits.out`, `golden_fermat_check.out`; audit E's independent counts in `audit_E_20261014/`):
- Brute-force orbit counts equal the formulas in all 27 cases (`a = 1, 2, 3`). These include the silver Wall–Sun–Sun primes:
  - `p = 13` (inert): 4080 orbits, against 336 otherwise;
  - `p = 31` (split): 64 orbits, against 4.
- `a`-Wall–Sun–Sun primes found:

  | `a` | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  |---|---|---|---|---|---|---|---|
  | primes | none below `3·10^5` | 13, 31 (known: also 1546463) | 241 | 3 | 3, 11 | 191, 643 | 5 |

- Audit E recounted by Burnside, with periods from direct iteration: 76 small cases (`a = 1..7`) and `p = 191, 241, 643, 1546463`, with no mismatch.
- `π(p) = 2^k ⟺ p | L_(2^j)`: no mismatch for `p < 10^6`. This was checked by `golden_fermat_check.py`, the direct modular check in `wss_orbits.py` (3), and audit E with PARI.
- `π(M_q) = 2^(q+1)` for `q = 2, 3, 7, 19, 31, 107, 127`. ∎

## Reading

- **The dictionary.** The two classical "Wieferich" conditions on the units `2 ∈ Z` and `φ ∈ Z[φ]` correspond to two orbit explosions of invariant circulant tournaments. The unit criterion itself is KNOWN; the tournament reading is the new part.
  - 2 is the generator behind the Collatz and doubling maps. `φ` is the generator behind Fibonacci recursion and the Pisano period.
  - The jump factor is `(p+1)/2` for 2, and `(p² + 1)/(p + 1)` for inert `φ`: a "quadratic" Wieferich phenomenon.
- **Rigidity duality.**
  - Fermat primes (`p − 1 = 2^k`) make the split torus a 2-group. Mersenne primes (`p + 1 = 2^q`, inert) make the golden nonsplit torus a 2-group.
  - The divisors of Fermat numbers (`ord_p 2` a power of 2) are mirrored by the golden-rigid primes (`ord φ` a power of 2).
- **The Pisano bound `π(n) <= 6n`** (KNOWN, Freyd–Brown 1992; equality iff `n = 2·5^k`) is the ramified prime 5 at work. `φ ≡ 3 (mod 5)` is a double root, and `5 | F_5` makes the 5-part of `π` grow by a full factor 5 at each level. See the note.
