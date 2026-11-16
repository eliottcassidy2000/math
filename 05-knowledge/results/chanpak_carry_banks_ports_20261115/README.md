# Carries between the dyadic and the triadic step: bank rigidity, port spectra, and 27, 63, 1473, 1729 in the descent tree

Session `mac-mini-2026-11-15-carry`, 2026-11-15. **The Chan–Pak conjecture is OPEN; nothing here proves it.**
Not independently audited (owner directive of 2026-11-13). Every statement carries its status.

`H(T)` = number of Hamiltonian paths; a *count* is an odd number that is `H` of a tournament; the *half-count* is
`n = (H - 1)/2`. `B + x` is a base tournament `B` with one new vertex `x` (the *head*), whose out-set is `X`.

## 0. Results at a glance

| # | statement | status |
|---|---|---|
| 1 | **Cut form.** `H(B + x) = H(B) + 2 Σ_{a in X, b not in X} W(a,b)`, where `W(a,b)` sums `H(B - V(P))` over the paths `P` from `a` to `b` with an even number of vertices. Every pair of head contacts interacts: the second difference is `-2(W(e,f) + W(f,e)) <= -2 H(B - e - f)`. | PROVED (one line from Claim A, THM-4876 (1)); brute-force check |
| 2 | **Bank rigidity.** If `2^k` head patterns `X0 ∪ Y`, `Y ⊆ S`, give distinct counts inside `L` consecutive odd numbers, then `H(B - e - f) <= 2(L - 1)` for all `e, f` in `S`. For a complete bank (`L = 2^k`): `H(B - e - f) <= 2^(k+1) - 2`. Independent low digits of one head exist only where deleting any two of them collapses the base. | PROVED |
| 3 | Complete 2-banks (4 consecutive counts from 2 contacts) exist on 1, 2, 5, 9, 15, 21 bases of orders 4–9; **no complete 3-bank exists on any base of order at most 9**, although 336 strictly submodular quadratic bijections of the 3-cube onto 8 consecutive integers exist abstractly. | FINITE-EXACT |
| 4 | Without the cube structure a single base carries long runs: 24 consecutive odd counts at order 8 (2763..2809, base count 657), 29 at order 9 (107..163, base count 5), 54 on a 12-vertex tape with one 4-block (337..443). The pure tape (one apex) gives 4, 5, 6, 6, 7, 9, 9 at orders 8–14. | FINITE-EXACT |
| 5 | **Port spectra** `Π(d)` (defects `δ` with `2d + δ` a one-vertex child of a realiser of `d`), complete for the 50 counts `d <= 104`: both digits `±1` for 26 counts; no `-1` at 11, 29, 35, 39, 43, 47, 53, 55, 63, 67, 71, 77, 85, 87, 101; no `+1` at 3, 13, 25, 31, 43, 73, 75, 83, 85, 87, 93. | FINITE-EXACT |
| 6 | **The commutator carry.** The triadic step `T: h -> 3h` (ordinal sum with a triangle) and the dyadic step `D: h -> 2h - 1` (tail at a port of defect -1) satisfy `D(T(h)) - T(D(h)) = 2`: the two orders give consecutive odd counts `6h - 3`, `6h - 1`. Both are physical for `h = 3, 5, 9, 15, 17, 19, 23, 25, 27, 31, 33`; a port is missing at 11, 13, 29. | identity PROVED; realisations FINITE-EXACT (`h <= 33`) |
| 7 | **8 and 9.** The dyadic ray `D^i(3) = 2^(i+1) + 1` and the triadic ray `T^j(3) = 3^(j+1)` from the triangle meet exactly once, at `9 = D²(3) = T(3)` (`2^a + 1 = 3^b` only for `(1,1)`, `(3,2)`); in half-counts `4 = 100₂ = 11₃`. This is the only relation among the numbers `2^a + 1`, it is why the floor ladder of THM-5262 climbs the powers of three to `243` and stops, and 21 is the first count whose two semigroup parents (`T(7)`, `D(11)`) are both blocked. | PROVED (arithmetic); consequences as tagged in section 4 |
| 8 | **27 and 63.** Writing `S` for the signed descent, `S³(63) = 27` (63 → 47 → 35 → 27 → 21: three rising steps on the chain of 21). The complete parent set `P(63) = {1, 11, 25, 33, 37}` and `27 ∉ P(41)` show that no count-27 tournament lies inside a count-63 one at any depth (agreeing with `coalescent_27_63_20261116`); `63 = 2·33 - 3`. | PROVED (arithmetic) / FINITE-EXACT |
| 9 | **1729 and 1472.** `1729 = 12³ + 1` has half-count `864 = 2^5·3^3`: it is the fifth point of the dead ray above 55 (half-count 27), and `1729 = H` of the one-apex word `11101001000`. `1473` has half-count `2^5·23` (ray above 47), `2945` has half-count 1472. They descend to 21 through the two branches: 1729 via 41, 1473 via 27; `1729 - 1473 = 2^8 = 2^5(55 - 47)`. | PROVED (arithmetic, literal check) |

Not claimed: a complete low bank at large counts, a selector for arbitrary targets, any density statement.

## 1. Sources and credit

* `coalescent_27_63_20261116` (another session): no 27-path tournament embeds in a 63-path witness; exact rewrites;
  square-sum ports for 8, 9, 1472, 1729; the golden square-sum carry (`s ≡ 1 mod b`).
* `physical_radix_carry_20261116` (another session): a physical radix-four rollover `(P,1) → (P,3) → (P+1,1)` on a
  transitive tape with one head, and `BUDGET.md`: the same masks cannot fill a larger low bank; "find a new source
  of independent low digits".
* `six_tutte_collatz_20261116` (the signed descent, written `S` here), THM-5262 (my previous note: floor, ladder, descent tree,
  parent sets), THM-4876 (1) (Claim A) and Session214 (floor).

This note answers the bottleneck question in one direction (sections 2–3: what one head can and cannot do, on any
base) and gives the arithmetic frame in which "carry" has an exact meaning for tournaments (sections 4–5).

## 2. The cut form and the rigidity of banks

**Cut form (1).** Claim A: `H(T) - H(T - v) = 2 Σ_C H(T - V(C))` over the odd cycles through `v`. A cycle through
the head is `x → a ⇝ b → x` with `a` in `X`, `b` outside `X` and a path `a ⇝ b` in `B` with an even number of
vertices. Hence the formula of row 1: the count is `H(B)` plus twice a directed cut of the weighted digraph `W`
on `V(B)`. `W` depends only on `B`: *a small graph that encodes the arithmetic of every head over this base.*
For the transitive tape, `W(a,b) = 2^(b-a-2)` (`b >= a + 2`), `1` (`b = a + 1`), and the half-count of a one-apex
word is `a - J(a)`: the binary value of the marks minus their pair energy.

Since `a → b` or `b → a`, the two-vertex path gives `W(e,f) + W(f,e) >= H(B - e - f) >= 1`. For `e, f` outside `X`:

    H(B + x_{X+e+f}) - H(B + x_{X+e}) - H(B + x_{X+f}) + H(B + x_X) = -2 (W(e,f) + W(f,e)) <= -2 H(B - e - f).

`scripts/cut_form_check.py`: 3,108 (base, out-set) pairs, no failure.

**Bank rigidity (2), PROVED.** Four counts inside `L` consecutive odd numbers differ by at most `2(L - 1)`, so a
second difference is at most `4(L - 1)` in absolute value. Therefore `H(B - e - f) <= W(e,f) + W(f,e) <= 2(L - 1)`
for every pair of contacts in the bank. So a bank of `2^k` consecutive counts at height `N` needs a base in which
*every two bank contacts form a feedback pair*: removing them leaves fewer than `2^(k+1)` Hamiltonian paths,
however large `N` is. On a transitive tape this forces the contacts into a window of `k + 3` consecutive tape
vertices (compare the budget bound of `physical_radix_carry_20261116`, proved there for its masks); the statement
here holds for every base.

**Census (3).** `scripts/bank.c` searches every base for the largest `k` with a context `X0` and a set `S` of `k`
contacts whose `2^k` patterns give `2^k` consecutive odd counts.

| base order | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|
| bases (isomorphism classes) | 4 | 12 | 56 | 456 | 6,880 | 191,536 |
| with a complete 2-bank | 1 | 2 | 5 | 9 | 15 | 21 |
| with a complete 3-bank | 0 | 0 | 0 | 0 | 0 | 0 |

Examples of 2-banks: the strong 4-tournament (counts 9, 11, 13, 15), a base of count 9 (89..95), a base of count 39
(209..215). The absence of 3-banks is not forced by the algebra: `scripts/quad_bank3.py` lists all 4,608
quadratic bijections of the 3-cube onto 8 consecutive integers; 336 have all three interactions negative, e.g.

    f(y) = (-y1 + 2 y2 + 4 y3) - (y1 + y2 + y3)²      with values 0, -2, 1, 3, -3, -1, 2, -4.

So tournaments could in principle carry three-digit banks; through order 9 none does.

**Runs (4).** `scripts/ext_runs.c`: for each base the longest run of consecutive odd numbers among its `2^n`
extension counts. Order 8: maximum 24 (2763..2809 on a base of count 657). Order 9: maximum 29 (107..163, on the
4-block with five tape vertices); best run by base count: 5 for the transitive base, 10 for count 3, 29 for
count 5. `scripts/tape_block.c` continues the tape families: the pure tape has runs 4, 5, 6, 6, 7, 9, 9 at orders
8–14 (about one more per two tape vertices); one 4-block on a tape of order 12 carries 54 consecutive counts
(337..443) on a single base, and over all placements the family covers 173 consecutive odd numbers
(1903..2247) at order 12 and 175 (5057..5405) at order 13. Long runs exist, but they are not cubes of contacts:
the low digit is spread over many contacts, as in the three-row table of `RADIX4.md`.

## 3. Port spectra

`Π(d)` = the odd `δ`, `|δ| <= 15`, such that some realiser of `d` has a one-vertex extension of count `2d + δ`.
`out/port_spectra_104.txt` lists `Π(d)` for the 50 counts `d <= 104` (complete: realisers of counts up to 223 have
strong cores of order at most 10). Sizes range from 5 to 16; `Π(15)` is every odd number from -5 to 15,
`Π(23)` every one from -9 to 15, `Π(11) = {-3, 1, 3, 7, 11, 13, 15}`, `Π(3) = {-1, 3, 7, 9, 11}`.
In half-counts a child is `2n + j` with `j = (δ + 1)/2`: a binary expansion whose digit set depends on the count.
A missing count is a missing digit at every smaller count; the nearest ones are `+1` at 3 (for 7) and `-1` at
11 (for 21).

## 4. The dyadic and the triadic step

`T: h -> 3h` is always physical (ordinal sum with a triangle). `D: h -> 2h - 1` is physical at `h` iff
`-1 ∈ Π(h)`; for a tail vertex at a port `u` of `A` the defect is `D_A(u) - E_A(u)` (THM-5262). In half-counts
they are `n -> 3n + 1` and `n -> 2n`.

**Commutator (6).** `T(D(h)) = 6h - 3`, `D(T(h)) = 6h - 1`. The unit of the low digit is the commutator of the two
scales. Literal pairs: (15, 17), (27, 29), (51, 53), (87, 89), (99, 101), (111, 113), (135, 137), (147, 149),
(159, 161), (183, 185), (195, 197). Longer words give the Collatz constants: a word with `i` letters `D` and `j`
letters `T` sends `h` to `2^i 3^j h - Σ_l 2^(i-l) 3^(t_l)`, `t_l` = number of `T` after the `l`-th `D`.

**The semigroup.** On half-counts the two steps generate the affine semigroup `<2n, 3n + 1>`. Its orbit of 0 has
at most `c X^σ` elements below `X`, `2^(-σ) + 3^(-σ) = 1`, `σ = 0.787885` (word count; measured exponents
0.755, 0.760, 0.763, 0.768 at `X = 10^4..10^7`): a self-similar set of density zero. Of its 33 counts up to 209,
26 are reached physically by `T` and (-1)-tails from the one-point tournament; 57, 105, 173, 201 are blocked by
missing ports of defect -1 at 29, 53, 87, 101, and 113, 171, 209 depend on the blocked 57 and 105.

**8 and 9 (7).** `2^a + 1 = 3^b` forces `(a, b) = (1, 1)` or `(3, 2)`: for `a >= 3` the left side is `1 mod 8`,
so `b = 2c` is even and `2^a = (3^c - 1)(3^c + 1)` is a product of two powers of two differing by 2, whence
`c = 1`; `a = 2` gives 5. So:

* the fans `3, 5, 9, 17, ...` and the triangle powers `3, 9, 27, ...` share only 9: the count 9 is `D²(3)` and
  `T(3)`, with witnesses the fan `1000` and two triangles;
* among the numbers `φ_a = 2^a + 1` the only multiplicative relation is `φ_3 = φ_1²`; hence
  `3^j = φ_a φ_b` exactly for `j <= 4` (3, 9, 27 = 3·9, 81 = 9·9), and the ladder of THM-5262 (C)
  (`2m - 1` inside `3m`) climbs 1 → 3 → 9 → 27 → 81 → 243 with half-counts 0, 1, 4, 13, 40, 121 and has no
  instance at `m = 243`;
* the half-count 10 of 21 is `2·5 = 3·3 + 1`: both semigroup parents exist arithmetically and both are blocked
  (7 is not a count; 11 has no port of defect -1). The half-count 3 of 7 has no semigroup parent at all.

The square-sum carry of `GOLDEN.md` is a unit `s ≡ 1 (mod b)` in a cyclic group; the tournament carry found here
is a commutator of two affine maps. Each is exact because of one small identity (`a² = 1 + bq` there,
`(3·2n + 2) - (2·3n + 1) = 1` here), and in each the irrational ratio (`φ` there, `log 3 / log 2` here) only
explains why the coincidence does not repeat: after 9 the two rays never meet again.

## 5. 27, 63, 1473 and 1729 in the descent tree

`out/descent_numbers.txt`. The signed descent (THM-5262, section 5) has the rising chain
`21 ← 27 ← 35 ← 47 ← 63 ← 83 ← ...` and the branch `21 ← 41 ← 55` whose rising child 55 is a dead ray.

* **63 → 47 → 35 → 27.** Three rising steps; none can be induced (THM-5262 (A)), and the whole relation is not
  induced either: the one-vertex deletions of a count-63 tournament have counts 1, 11, 25, 33, 37, and a chain
  27 → x → 63 of two strict extensions would need `x = 41` with `27 ∈ P(41)`, which fails. Its nearest parent is
  the fan count 33. One-apex witnesses: `27 = H(100010)`, `63 = H(1001100)`.
* **The Collatz meeting.** `62 = T²(27)` and 63 have the same numbers of odd and even shortcut steps to 182, with
  constants differing by `3^6`; so `27 + 2^10 j` and `63 + 2^8·9 j` both reach `182 + 3^8 j` (checked `j < 200`).
  For a count `h` the first shortcut step `(3h + 1)/2` is the increment floor; it is attained only at `3 → 5`,
  so these orbits leave the realisable extensions at the first step (`27 ∉ P(41)`, `63 ∉ P(95)`).
* **1729.** Half-count `864 = 2^5·3^3`: `1729 = D^5(55)` (dyadic step) on the dead ray `55, 109, 217, 433, 865, 1729`, where 55
  is the count of half-count 27, the classical long-orbit start. Descent: `... → 55 → 41 → 21`. Witnesses: the
  one-apex word `11101001000` (12 vertices, literal check), and the ordinal products `13·133`, `19·91` (the two
  cube decompositions of `PORTS.md`); `7·247` is not available.
* **1472.** As a half-count it is the count `2945 = 5·19·31`; the neighbour `1473 = 2·736 + 1` has half-count
  `2^5·23` and lies on the ray `47, 93, 185, 369, 737, 1473, 2945`. Descent: `... → 47 → 35 → 27 → 21`.
  1473 is not a one-apex count. So `1729 - 1473 = 2^8 = 2^5 (55 - 47)`: the two numbers are the fifth dyadic
  ancestors of the rising children of 41 and 35, one through each branch above 21; the `257` of `PORTS.md`
  is this `2^8` plus the parity offset.

## 6. Open

1. A complete 3-bank on any base (HYP-9822 (a)); the growth of single-base runs with the order.
2. A class of pointed tournaments closed under tail extension in which both digits `±1` survive (HYP-9822 (b)).
3. Chan–Pak.

## Files

`scripts/`: `cut_form_check.py`, `bank.c`, `ext_runs.c`, `tape_block.c` (inputs: out-mask lines; `tourn_masks.py`
converts `gentourng` output), `quad_bank.py`, `quad_bank3.py`, `port_spectra.py` (reads `lift_near_order9.txt` of
`chanpak_signed_descent_lifts_20261115/out`), `descent_numbers.py`. `out/`: the outputs quoted above.
