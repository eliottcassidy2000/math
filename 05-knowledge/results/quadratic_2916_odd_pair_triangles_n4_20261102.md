# x -> x^2 - 29/16 against the n = 4 metagraph, primes modulo even numbers until +-1, and the jagged triangles (odd pairs; PRIMEGAME's trial-division region)

mac-mini-2026-11-02-triangles, 2026-11-02. All work committed directly to `main` and pushed after each step.

## 0. The owner's prompt (pasted, verbatim)

"Consider prior work on the hourglass atom on 4 vertices. also on {10,20,39} and the concept of primes modulo even
numbers until they achieve the value of 1 or -1. up until that point, they have gone through the ways for primes to
sum to eachother, with 2 connecting each even number with the one prior, and then patterns in the odds connected to
eachother. consider fractran and turing machines and cellular automata and selecting pairs of odd numbers from this
jagged but triangular like region, consider what primegame and other representations could look like in these formats
and get creative and free in your pattern seeking. connect back with previous work and think how the below connects
with the 8 tilings in the n=4 tournament metagraph condensing down to the 4 isomorphism class nodes in the metagraph
(or 3 in the merged metagraph): for the map x ↦ x² − 29/16, we can classify every orbit starting at an odd integer
divided by 4. Exactly eight such starting values are bounded; all enter the cycle −7/4 → 5/4 → −1/4 → −7/4"

## 1. Prior work (credited)

- **x -> x^2 - 29/16.** THM-4139 (codex, 2025-08-25): the complete rational preperiodic graph is the eight
  quarter-integers plus the exact 3-cycle, so there is no rational 6-cycle. This is classical (Poonen 1998); c = -29/16
  is the standard example with a rational 3-cycle. The owner's statement is this theorem restricted to odd/4 starts.
- **"Hourglass."** In the tiling literature here, hourglass cycles are HYP-3053's bridge-to-bridge cycles of the
  diagonal-layer graph: C(n-2, 2) of them, so exactly one at n = 4. The Pythagorean hourglass of THM-4773 / THM-4780
  is a different object.
- **{10, 20, 39}.** These are the clique-mean carrier dimensions of the concurrent packet `chart_clique_norm_20261103`
  (three colours, forbidden clique sizes (3,3,5)), at n = 3, 4, 5:
  1 + 2·C(n,2) + sum_{k=3}^{4} C(n,k) + #{t_i <= n}. Reproduced here as 10, 20, 39 (then 69 at n = 6).
  - Coincidences only (NUMEROLOGY): 10 and 20 are even numbers next to a prime (11; 19), and 39 is the square-sum hole
    of THM-4851.
- **Goldbach machinery.** GOLDBACHGAME (THM-4856); the 25-state machine and its CA and FRACTRAN forms (THM-4858); the
  convolution picture of prime pairs (session65).

**Typing.**
- STRUCTURAL: a predicate-preserving map.
- NUMEROLOGY: equal counts or labels.
- GENERIC: forced, or true of every object of the class.

## 2. Results

1. **The 8 -> 4 -> 3 chain is real for the map.** On the eight bounded points B = {±1/4, ±3/4, ±5/4, ±7/4} (pre-registered
   Q1, Q2, PASSED):
   - f maps B onto 4 points {-7/4, -5/4, -1/4, 5/4}, with fibres the pairs {x, -x};
   - f maps those onto the 3-cycle, with fibre shape (2,1,1);
   - depths to the cycle: 0 for {-7/4, 5/4, -1/4}, 1 for {-5/4, 1/4, 7/4}, 2 for {±3/4}.

   On numerators, the map is the integer map g(a) = (a^2 - 29)/4 on odd a. Its orbit is bounded iff |a| <= 7
   (checked for |a| <= 100,001; for |a| >= 9, g(a) > |a|). The eight bounded numerators {±1, ±3, ±5, ±7} are the units
   mod 16.
2. **Against the n = 4 metagraph** (M1-M3 PASSED).
   - **The metagraph chain.** 8 tilings -> 4 classes, with fibres {5,1,1,1} (TT4 1; source+C3 1; C3+sink 1; the
     strong 4-tournament 5) -> 3 merged nodes, with fibres {5,2,1}. There are 4 lines.
   - **Verdict: NUMEROLOGY.**
     - The counts 8, 4, 3 agree.
     - The first-stage fibre shapes differ: {5,1,1,1} against {2,2,2,2}.
     - The second-stage shape (2,1,1) is forced for any map of 4 onto 3 (GENERIC).
     - The complement lines and the ±x pairs are both fixed-point-free involutions on 8 points, hence conjugate
       (GENERIC).
     - No bijection carries the class line multigraph (a star at the strong class, with a loop there) onto f's
       functional graph on f(B) (a triangle with a pendant edge) (M3).
   - **Further failed structural reading.** (Z/16)^× ≅ Z/2 × Z/4 is not GF(2)^3, the tiling cube.
   - **What survives.** The shared shape: an 8-set, folded by a fixed-point-free involution and one further merge. In
     the metagraph that merge is the converse; in the map, the f-fold of a 4-point set onto the cycle.
3. **Primes modulo even numbers until +-1** (P1 PASSED; post-hoc proposition proved and checked). For even m, scan the
   odd primes in increasing order, recording residues mod m, until the first prime P(m) = ±1 (mod m).
   - **Bound.** P(m) < m^2 for every even m <= 20,000. P(m) = m ± 1 for 41.8% of them, and the largest P(m)/m is 17,
     at m = 16012.
   - **Proposition (PROVED, elementary).** Let R(m) be the set of residue pairs {r, m - r} that the scan hits and G(m)
     the Goldbach partitions of m into odd primes.
     - G(m) ⊆ R(m) always.
     - If m - 1 or m + 1 is prime, R(m) = G(m). The scan then stops at m ± 1, having seen exactly the odd primes below m.

     This is the owner's "up until that point, they have gone through the ways for primes to sum to each other". It
     is exact when m has a prime neighbour.
   - **Open observation (FINITE-EXACT, not proved).** For all 5,820 even m <= 20,000 with no prime neighbour,
     R(m) ⊋ G(m). The extra pairs are wrap-around residues of primes in (m, P(m)). A proof for all m would need a
     Goldbach-type existence statement.
   - **The prime 2.** It has residue 2, and the pair (2, m - 2) joins m to the previous even number m - 2: "2
     connecting each even number with the one prior". Since m - 2 is even, 2 never enters G(m) for m > 4.
4. **The jagged triangles and their machine forms.**
   - **The odd-pair triangle.** Rows m = 2k; entries are odd a <= k, marked when a and m - a are both prime (printed
     in `quadratic_triangle_n4.out`).
     - As a cellular automaton, row m + 2 is row m with its second factor shifted one cell: χ(m+2-a) = χ(m - (a-2)).
       So the triangle is the space-time diagram of a shift register run against the static prime indicator, i.e. the
       convolution r(m) = sum_a χ(a) χ(m - a) (session65's carrier). GENERIC.
     - GOLDBACHGAME selects, in each row, the first marked cell from the left (its least prime split). The Leng
       automaton (THM-4858) sweeps the same rows in unary.
   - **PRIMEGAME's region.** PRIMEGAME sweeps the jagged triangle {(n, d) : D(n) <= d <= n - 1}: row n has n - D(n)
     trial divisors, and exactly the prime rows run down to d = 1 (THM-4853's cost is a sum over this region).
     - These are "numbers modulo descending moduli until the remainder is 0". The owner's scan is "primes modulo
       increasing even moduli until the residue is ±1". Both select pairs from a jagged triangular region by a
       first-hit rule.
     - In machine terms each is a nested counter loop: FRACTRAN 2-cycles, or Turing scanners.

## 3. Reproduction

Run these from `04-computation/experiments/circulant_fractran_goldbach_20261101/`:

```text
python3 quadratic_triangle_n4.py     # < 1 s; Q1-Q2, M1-M3, P1, P3 (pre-registered) and the triangles
python3 residue_pairs_posthoc.py     # ~15 s; the residue-pair proposition and the open observation (post-hoc)
```

## 4. Audit

(pending)
