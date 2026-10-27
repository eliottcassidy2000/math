# The Eisenstein grid as a half-shifted copy of the first grid; 13 and 37; sums of four fourth powers

2026-10-27, mac-mini-2026-10-27-slab. **Status:** THM-4824 (PROVED parts plus FINITE-EXACT checks) and HYP-9428 (OPEN),
with a typed table of the numerical connections. Independent audit pending at the pre-audit checkpoint.

## Owner's directive (verbatim)

"consider the idea that the eisenstein integers are just another coordinate grid on top of the first grid, but centered
at (1/2,1/2), or something along those lines for inspiration. consider prior work on 13 and 37 and consider that of the
nonnegative integers up through 50, 13 natural numbers can be written as the sum of 4 4th powers. 0 to4,16 to 19, 32 to
34,48 to 49 =5+4+3+2, while 37 cannot. the lone 50 + ranges of 11,12,13. consider previous work connecting
5,10,25,50,100 and 11,12,13 in the aliquot related work."

## 0. What the repository already had

- **13 and 37.** In the session69 power/aliquot audit
  (`session69_prime_ports_power_aliquot_20261026/session69-arithmetic-audit/REPORT.md`, §6), the "square-later" integers
  give `(X, qualifying, even failures) = (100, 37, 13)`.
  - **Definition.** Write n = 2^a·u with a ≥ 1 and u odd. Then n is square-later iff u > 2^a.
  - **The split.** Of the 50 evens up to 100, 37 qualify and 13 fail. This is the owner's 13 + 37 = 50, with 50 here
    counting evens.
- **5, 10, 25, 50, 100.** The mac-mini-2026-10-25-pyth note `fresh_remarking_quotient_and_divisor_arrangements_20261026.md`
  §6 and `power_doubling.out` record that "each odd non-power b₀ gives a grid {2^i b₀^j}".
  - **The five numbers.** For b₀ = 5 the grid contains 5, 10, 25, 50, 100 at (i, j) = (0,1), (1,1), (0,2), (1,2), (2,2).
  - **How they connect.** They are joined by doubling and squaring: 10 = 2·5, 25 = 5², 50 = 2·25, 100 = 10² = 4·25.
- **11, 12, 13.**
  - **Aliquot fibres** (`aliquot_fibres_to_45_20261027.out`): s⁻¹(11) = {21}, s⁻¹(12) = {121 = 11²} (a unique
    prime-square source), s⁻¹(13) = {27, 35}.
  - **Lonely runners.** The session66 aliquot packet's copy of CORE-PAPERS cites "Eleven, twelve, and thirteen lonely
    runners" (Sungkawichai–Trakulthongchai, arXiv 2604.23906, PREPRINT).
- **Not found.** No prior statement tying {5, 10, 25, 50, 100} to {11, 12, 13} in the aliquot work was found, in any
  file type, agent letters included. The nearest items are the ones listed above.

## 1. The two grids, made exact (THM-4824)

**The decomposition.** In the coordinates (u, v) of (u + v√−3)/2, Z[ω] = {u ≡ v (mod 2)}.
- **The first grid.** Z[√−3] is the set with u, v even.
- **The second grid.** The set with u, v odd is a copy of the first, centred at (½, ½) in its units. Its norms are always
  odd.
- **The units.** They split 2 + 4: ±1 in the first grid, and ±ω, ±ω² in the shifted grid.

This is the owner's picture, and it is the mechanism of the merge debt:

| lattice | unit charges | catalytic merge debt |
|---|---|---|
| first grid Z[√−3] alone | ±1 (collinear) | g(n) = n(n−1) on the real line; catalysts are useless; nothing else is reachable (THM-4824 (B); this is THM-4820's f(n,0,0)) |
| glued with the (½,½) copy, i.e. Z[ω] | 6 | ⌊Q/3⌋ plus the exceptional set (THM-4820) |
| square grid Z[i] | 4 | **zero** at every point (THM-4824 (A1)) |

So the half-shift glue is what makes the catalysis possible, and the 60°/120° geometry of the glue is what keeps it from
being total.

**A critical point.** In the glued family {u ≡ v (mod 2)}, Q = u² + c v², the Eisenstein grid c = 3 is exactly the shape
at which the centred points become as short as the first grid's basis vector. It is also the last shape at which they
are minimal.
- **Below c = 3:** four glue units.
- **At c = 3:** six units.
- **Above c = 3:** two collinear units.

**Debt by aspect.**
- c = 1 is the square grid: zero.
- c = 2 has the floor 2⌊Q/8⌋, with its own exceptions.
- c = 3 has ⌊Q/3⌋.
- c = 3/2 and 5/2 have measured ratios of about 0.13 and 0.28.

A tempting angle law fits three of these points and is refuted by the other two (HYP-9428 (2)).

**Beyond the plane.**
- **Glue in other dimensions.** The 8-dimensional half-shift glue E₈ = D₈ ∪ (D₈ + (½)⁸) is free, and so are the
  checkerboards D_n.
- **The conjecture.** A₂ is the only root lattice with a positive floor (HYP-9428 (1)). Its distinguishing property is
  that no two of its roots are perpendicular.

**Parity cosets.** The four cosets of 2Z[ω] refine the two grids. Every tight split of x is a point of x + 2Z[ω] on two
or four norm shells (THM-4824 (E)). This is the all-type form of Session75's parity-shell criterion, and it gives a
finite, purely arithmetic seed test. For example, 51 is a seed.

## 2. Sums of four fourth powers (zeros allowed)

**Below 81.** The only fourth powers available are 0, 1 and 16. So the representable set is the triangle
{16a + b : a + b ≤ 4}, with 15 points.
- **The composition simplex.** These are the compositions of 4 into the three types 0, 1, 16. This is the same simplex
  as THM-4820's charge vectors (p, q, r), here weighted by 0, 1, 16 instead of 1, ω, ω².
- **Up to 50.** The representable positives number 4 + 4 + 3 + 2 = 13: namely 1–4, 16–19, 32–34 and 48–49. The other
  37 are not representable.
- **Gaps.** The gaps have lengths 11, 12, 13, 14 = 16 − 5 + a, for a = 0, 1, 2, 3. Within 1..50 the last gap is the
  lone 50.
- **Mod 16.** Odd bases give 1 (mod 16) and even bases 0 (mod 16), so a sum of four fourth powers is ≡ 0, …, 4 (mod 16).
  The eleven residues 5, …, 15 never occur. This is the two-grid (even/odd) obstruction, and it is why 79 = 4·16 + 15
  needs 19 fourth powers.

| N | 5 | 10 | 25 | 50 | 100 | 200 | 500 | 1000 | 10⁴ |
|---|---|---|---|---|---|---|---|---|---|
| representable in 1..N | 4 | 4 | 8 | **13** | 21 | 30 | 53 | 96 | 540 |
| not representable | 1 | 6 | 17 | **37** | **79** | 170 | 447 | 904 | 9460 |

**The zero-sum quartic identity** (THM-4824 (F), classical). If a + b + c = 0, then

    a⁴ + b⁴ + c⁴ = 2(a² + ab + b²)² = (2/9) Q(a + bω + cω²)².

So the fourth-power sums of zero-sum triples are twice squared Eisenstein norms. Of the 13 representable numbers up to
50, the ones of this form are 2, 18 and 32. The zero-sum condition is the catalyst condition 1 + ω + ω² = 0 read on the
coefficients.

## 3. 13 and 37, typed

| statement | type |
|---|---|
| 13 and 37 are the two smallest primes ≡ 1 (mod 12), the primes that split in both grids, Z[i] and Z[ω]: 13 = 2² + 3² = 1 + 3·2², 37 = 1² + 6² = 5² + 3·2² | IDENTITY |
| the primitive Pythagorean triangles with a leg 12 have hypotenuses exactly {13, 37}: 5-12-13 and 35-12-37 | IDENTITY (via (c−b)(c+b) = 144) |
| the primitive 120° triangles (c² = a² + ab + b²) with a side 7 have longest sides exactly {13, 37}: 7-8-13 and 7-33-37 | IDENTITY (via (2c−2b−7)(2c+2b+7) = 147) |
| the same pair arising in both grids | FINITE-EXACT: among the primes ≡ 1 (mod 12) below 200 (13, 37, 61, 73, 97, 109, 157, 181, 193), {13, 37} is the only pair sharing both a Pythagorean side and a 120° side (`counts.py` §6). For example, 61 and 109 share the leg 60 but no 120° side. No general law is claimed (NUMEROLOGY) |
| 13² = 11² + 3·4² and 37² = 13² + 3·20²: the squares of the first-grid primes above 13 and 37 give the chain 11 → 13 → 37 | IDENTITY each; the chain is NUMEROLOGY |
| 13 = #Z[ω] with Q ≤ 3, the star number S₂; 37 = #Z[ω] with Q ≤ 9 = S₃ = centred hexagonal H₄; also 13 = #Z² with x² + y² ≤ 4 and 37 = #Z² with x² + y² ≤ 10 | IDENTITY each; common ball counts of the two grids are 1, 13, 37, 61, 97, 109, 121, 169, … (all ≡ 1 mod 12) |
| 13 of 1..50 are sums of four fourth powers, and 37 are not | IDENTITY (§2) |
| session69: 37 square-later evens and 13 failures up to 100 | IDENTITY |
| the two 13/37 splits are the same phenomenon | NUMEROLOGY. The sets differ (failures: 2, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64, 80, 96; representable: 1–4, 16–19, 32–34, 48, 49), and at doubled cutoffs the splits are 81/19 and 21/79 |
| 79 non-representable numbers in 1..100, which is also the number needing g(4) = 19 fourth powers; 19 square-later failures up to 200; 8, 13, 21 representable at 25, 50, 100 (Fibonacci), breaking at 200 (30) | NUMEROLOGY |

## 4. 5, 10, 25, 50, 100 and 11, 12, 13

| n | Gaussian norm | Eisenstein norm | square-later | s(n) |
|---|---|---|---|---|
| 5 | yes | no | odd | 1 |
| 10 | yes | no | qualifies (u = 5 > 2) | 8 |
| 25 | yes (= 3² + 4²) | yes (= Q(5)) | odd | 6, perfect |
| 50 | yes (1² + 7² = 5² + 5², the first two-way sum of two positive squares) | no | qualifies | 43 |
| 100 | yes | yes (= Q(10)) | qualifies | 117 |
| 11 | no | no | odd | 1 |
| 12 | no | yes (= Q(2√−3)) | fails (u = 3 ≤ 4) | 16 |
| 13 | yes | yes | odd | 1 |

- **Which grid sees which numbers.** The two grids separate the doubling chain of 5: the Gaussian grid sees all five
  numbers, the Eisenstein grid only the squares 25 and 100. Of 11, 12, 13, the Eisenstein grid sees 12 and 13 and the
  Gaussian grid only 13.
- **Square-later.** 10, 50 and 100 all qualify, and 12 is one of the 13 failures up to 100.
- **Aliquot.** The aliquot fibres of 11, 12, 13 have the shapes pq (21), p² (121 = 11²), and p³, pq (27, 35). These are
  the shapes of THM-4805 (A'), which already appear together at 6, 7, 8.

All of this is IDENTITY at the level of individual facts. No structural link between the two triples is claimed.

## 5. What is new, and what is open

- **New, PROVED (THM-4824).**
  - The square grid, every D_n and E₈ are free.
  - The first grid alone is collinear.
  - The Eisenstein grid is the critical (tie) member of the half-shift glued family.
  - The aspect-2 floor 2⌊Q/8⌋, and the window lemma.
  - The parity-coset shell lemma.
- **Classical.** The zero-sum quartic identity, and the structure of four fourth powers below 81.
- **Open (HYP-9428).**
  - A₂ is the only paying root lattice.
  - The glued-family ratio κ(c).
  - The aspect-2 exceptional set (on-axis primes ≡ 5 mod 8).
  - The 13/37 coincidences are typed numerology.

## Reproduction

All scripts are in `04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/`:
- `counts.py`: fourth powers, square-later, lattice balls, the quartic identity, the triangle facts;
- `grid_debt.py 200`: the square grid, the first grid, and Z[ω];
- `root_lattice_debt.py 40` and `root_lattice_debt2.py`: the glued family and the root lattices;
- `e8_check.py`;
- `aspect_ratio_debt.py`, `aspect2_exceptions.py`, `aspect2_axis.py`.

Each takes seconds to a few minutes.
