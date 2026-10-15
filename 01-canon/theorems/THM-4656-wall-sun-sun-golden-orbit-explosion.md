---
id: THM-4656
title: "Wall-Sun-Sun primes are the golden Wieferich primes of THM-4644: the unit-invariant Cayley tournaments on Z[mu]/(p^2) explode exactly at metallic Wall-Sun-Sun primes (verified at the Pell primes 13, 31); golden rigidity holds exactly when the Pisano period is a power of 2, i.e. p | L_(2^j), which includes every Mersenne prime 2^q - 1 with q = 3 (mod 4) (Lucas-Lehmer with seed 3); golden Paley primes"
status: >
  PROVED (elementary; the same orbit argument as THM-4644 part 2, with the metallic unit mu_a in place of 2). KNOWN inputs:
  Wall 1960 (pi(p^2) in {pi(p), p pi(p)}); the equivalences of Wall-Sun-Sun primes (Sun-Sun 1992; pi(p^2) = pi(p) iff
  p^2 | F_(p-(5/p)) iff L_p = 1 mod p^2); no Wall-Sun-Sun prime below 2^64 (PrimeGrid 2020-2022); Lucas-Lehmer with seed 3.
  FINITE-EXACT: brute-force orbit counts for 27 cases (a = 1, 2, 3), including the metallic Wall-Sun-Sun primes p = 13, 31
  (a = 2), where the jump is observed; the golden-Fermat list for p < 10^6.
source: mac-mini-2026-10-14-tourn
depends_on: [THM-4644 (Fermat and Wieferich primes as circulant-tournament phenomena, opus S33), Wall 1960, Sun-Sun 1992]
verification:
  - 04-computation/experiments/wall_sun_sun_20261014/wss_orbits.py (+ .out)
  - 04-computation/experiments/wall_sun_sun_20261014/golden_fermat_check.py (+ .out)
note: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
related: [THM-4644, THM-4631, THM-4641 (Wilson / {2,4,8} / Fermat dictionary), HYP-9278 (S35: rigidity circuits)]
---

# THM-4656: Wall–Sun–Sun primes as golden Wieferich primes

## Setting

- For an integer `a >= 1`, let `μ = μ_a` be a root of `x² − a x − 1`, so `μ₁ = φ` is the golden ratio and `μ₂ = 1 + √2` is the silver ratio. Put `D = a² + 4`, and let `R_k = Z[μ]/(p^k)`.
- For a prime `p ∤ 2D`, let `π_a(p)` be the order of `μ` in `R_1^*`. This is the period of `U_n(a, −1)` mod `p`; `π_1` is the Pisano period.
  - If `p` splits, that is `(D/p) = +1`, then `π_a(p)` divides `p − 1`.
  - If `p` is inert, then `π_a(p)` divides `2(p + 1)`.
- `π_a(p²)` is either `π_a(p)` or `p·π_a(p)` (Wall). The prime `p` is an **`a`-Wall–Sun–Sun prime** when `π_a(p²) = π_a(p)`. For `a = 1` these are the Wall–Sun–Sun primes; none exists below `2^64`.
- A Cayley tournament on an abelian group `G` of odd order is a set `S ⊂ G \ {0}` with `S ∩ (−S) = ∅` and `S ∪ (−S) = G \ {0}`.

## Statement

1. **(Golden Wieferich bases; split primes.)** Let `p` split and let `μ₁ ∈ Z/p²` be a root of `x² − a x − 1` (Hensel lift). Then `p` is an `a`-Wall–Sun–Sun prime iff `μ₁^(ord_p μ₁) ≡ 1 (mod p²)`, i.e. iff `μ₁` is a Wieferich base at `p`.

   Consequently THM-4644 part 2 holds verbatim with the base 2 replaced by `μ₁`. Put `d = ord_p μ₁`, `d'` its odd part, and `g = μ₁^(d/d')`. The circulant tournaments on `Z_(p²)` with `gS = S` number `2^E`, where

   `E = (p − 1)/d'` normally, and `E = (p + 1)(p − 1)/(2d')` at `a`-Wall–Sun–Sun primes.

2. **(Inert primes: orbit explosion on `Z[μ]/(p²)`.)** Let `p` be inert, `π = π_a(p)`, `π'` its odd part, and `g = μ^(π/π')`. The Cayley tournaments on `(R_2, +) ≅ (Z/p²)²` with `gS = S` number `2^E`, where

   `E = (p + 1)(p² − 1)/(2π')` normally, and `E = (p² + 1)(p² − 1)/(2π')` at `a`-Wall–Sun–Sun primes.

   So the exponent jumps by the factor `(p² + 1)/(p + 1) ≈ p` exactly at Wall–Sun–Sun-type primes, mirroring THM-4644's factor `(p + 1)/2` at Wieferich primes.

3. **(Golden rigidity, the mirror of Fermat primes.)** Let `p ∤ 10`. The following are equivalent:
   - every subgroup `1 ≠ H <= ⟨φ⟩ ⊂ R_1^*` contains `−1`, so no Cayley tournament on `(R_1, +)` is invariant under a nontrivial power of `φ`;
   - `π(p)` is a power of 2;
   - `p` divides a Lucas number `L_(2^j)`, and then `π(p) = 2^(j+2)`.

   These primes start `3, 7, 47, 127, 1087, 2207, 4481, 21503, 34303, 119809, 524287` (`p < 10^6`). Every Mersenne prime `M_q = 2^q − 1` with `q ≡ 3 (mod 4)` is among them: `M_q` is inert, with `2(M_q + 1) = 2^(q+1)`, so `π(M_q)` is a power of 2. In fact `π(M_q) = 2^(q+1)`, the maximum: `M_q | L_(2^(q−1))`, which is the Lucas–Lehmer test with seed 3.
   - Fermat primes (`p − 1 = 2^k`) are the rigid primes for the split torus `F_p^*`.
   - Mersenne primes (`p + 1 = 2^q`) are the rigid primes for the nonsplit torus containing `φ`.

4. **(Golden Paley primes.)** Let `p ≡ 3 (mod 4)` split in `Q(√5)`. Exactly one root `φ₁` of `x² − x − 1` is a quadratic residue, since the product of the roots is `−1`. If it generates `QR_p`, then the circulant tournaments on `Z_p` invariant under `φ₁` are exactly the Paley tournament `P_p` and its reverse (THM-4644 part 3 with base `φ₁`). Below 2000 this happens for 52 primes, starting `11, 19, 31, 59, 71, 79, 131, …`.
   - At `p = 11`, `φ₁ = 4` has order 5, and `⟨4⟩ = {1, 3, 4, 5, 9} = QR_11`. So `P_11` is the golden tournament on `Z_11`. Its 55 arcs equal `|Aut P_11| = 11·5`; for every Paley tournament, `Aut(P_p)` acts regularly on the `C(p,2)` arcs.
   - That this 55 is also `F_10 = F_(π(11))` is a COINCIDENCE: by Luo's theorem, 55 is the largest triangular Fibonacci number.

## Proof

**1.**
- With `μ₂ = −μ₁^(−1)`, we have `π_a(p^k) = lcm(ord_(p^k) μ₁, ord_(p^k) μ₂)`.
- `−1` has order prime to `p`, so `μ₂` lifts with the same `p`-power as `μ₁`. Hence `π_a(p²) = π_a(p)` iff `ord_(p²) μ₁ = ord_p μ₁`.
- THM-4644's proof of part 2 uses only `ord_(p²)(base) ∈ {d, pd}` and the odd order of `⟨g⟩`. So it applies verbatim.

**2.**
- `R_2` is local, with maximal ideal `pR_2` and `|R_2| = p⁴`.
- `⟨g⟩` has odd order `π' p^ε`, where `ε = 0` iff the prime is `a`-Wall–Sun–Sun.
- It acts freely on the `p²(p² − 1)` units, giving `(p² − 1)p^(2−ε)/π'` orbits.
- On `pR_2 \ {0} ≅ F_(p²)^*` it acts through `g mod p`, which has order `π'`, giving `(p² − 1)/π'` orbits.
- Negation fixes no orbit, because `−1 ∉ ⟨g⟩`.
- So the count is `2^(N/2)`, where `N` is the total number of orbits.

**3.**
- `⟨φ⟩` is cyclic of order `π(p)`.
- If `π(p) = 2^k`, every nontrivial subgroup contains the unique involution `φ^(2^(k−1))`. That involution is `−1`, the only element of order 2 that lies in `⟨φ⟩` in both the split and the inert case (`φ₂ = −φ₁^(−1)`).
- If `π(p)` has an odd prime factor `r`, the subgroup of order `r` avoids `−1`, and unions of its cosets give invariant tournaments, as in THM-4644 part 1.
- For the Lucas form, `L_n = φ^n + φ^(−n)` for even `n`. So `p | L_(2^j)` iff `φ^(2^(j+1)) = −1` iff `ord(φ) = 2^(j+2)`.
- The Mersenne case follows from `2^q ≡ 3 (mod 5)` for `q ≡ 3 (mod 4)`, so `M_q ≡ 2 (mod 5)` is inert.

**4.**
- The product of the two roots is `−1`, which is a nonresidue since `p ≡ 3 (mod 4)`. So exactly one root is a residue.
- If that root generates `QR_p`, of odd order `(p−1)/2`, then its cosets are `QR_p` and `−QR_p`.

**Checks** (`wss_orbits.out`, `golden_fermat_check.out`):
- Brute-force orbit counts equal the formulas in all 27 cases (`a = 1, 2, 3`). These include the silver Wall–Sun–Sun primes:
  - `p = 13` (inert): 4080 orbits, against 336 otherwise;
  - `p = 31` (split): 64 orbits, against 4.
- `a`-Wall–Sun–Sun primes found:
  - `a = 1`: none below `3·10^5`;
  - `a = 2`: 13, 31 (known: 13, 31, 1546463);
  - `a = 3`: 241;
  - `a = 6`: 191, 643.
- `π(p) = 2^k ⟺ p | L_(2^j)`: no mismatch for `p < 10^6`.
- `π(M_q) = 2^(q+1)` for `q = 3, 7, 19, 31, 107, 127`. ∎

## Reading

- **The dictionary is exact.** The two classical "Wieferich" conditions on the two units `2 ∈ Z` and `φ ∈ Z[φ]` correspond to two orbit explosions.
  - 2 is the generator behind the Collatz and doubling maps; `φ` is the generator behind Fibonacci recursion and the Pisano period.
  - The jump factor is `(p+1)/2` for 2 and `(p² + 1)/(p + 1)` for inert `φ`, a "quadratic" Wieferich phenomenon.
- **Rigidity duality.** Fermat primes (`p − 1 = 2^k`) make the split torus a 2-group. Mersenne primes (`p + 1 = 2^q`, with `q ≡ 3 (mod 4)` so that `p` is inert) make the golden nonsplit torus a 2-group, and the Lucas–Lehmer test is exactly the statement that `φ` then has maximal order.
- **The Pisano bound `π(n) <= 6n`** (equality iff `n = 2·5^k`) is the ramified prime 5 at work. `φ ≡ 3 (mod 5)` is a double root, and `5 | F_5` makes the 5-part of `π` grow by a full factor 5 at each level. See the note.
