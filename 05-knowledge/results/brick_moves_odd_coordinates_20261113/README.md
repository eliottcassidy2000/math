# Brick moves, the two integer coordinates of the odd numbers, {3,10} and {4,11}, and FRACTRAN

mac-mini-2026-11-13-cuts (fifth prompt) · **Chan–Pak and HN OPEN.**

Status: **FINITE-EXACT** (§1). **ELEMENTARY** (§2). **PREDICTIONS** (§4), each with its test. Script:
[`scripts/brick.py`](scripts/brick.py), which needs the `hl*.txt` data of the busch_staircase note.

## 1. The 12-vertex test: the brick law, corrected and confirmed

Last session predicted that two strong-component compositions of n vertices share a path count exactly when a Busch
brick, 3 or 5, moves between components. **As stated, this is FALSE.**
- At 11 vertices, 10 of the 71 shared counts across (4,7), (5,6) and (3,8) are not {3,5}-moves. For example,
  165 = 5·33 = 11·15 trades a 5 for an 11.
- At 12 vertices, 61 of the 468 shared counts across (4,8), (5,7), (6,6) and (3,9) are not {3,5}-moves.

**Corrected law (FINITE-EXACT, n = 11 and 12, all two-component compositions).** Every shared count v = xy = x′y′ is
an S-unit re-association: x/x′ has numerator and denominator built only from primes in S = {3, 5, 11, 13}. These are
exactly the primes of the order-5 strong spectrum {9, 11, 13, 15}. There are 0 exceptions among 71 + 468 shared
pairs. {3, 5} alone suffices for the comparisons (4,7)/(3,8), (4,8)/(6,6), (4,8)/(3,9) and (5,7)/(6,6). The others
need an 11 or a 13 to move: for example 165 = 5·33 = 11·15, and 561 = 33·17 = 3·187, where 33 = 3·11 is a 6-core
count. This is consistent with P1 below.

**FRACTRAN reading.** A two-component tournament is a two-register machine. A collision is one FRACTRAN step
(x, y) → (x·u, y/u), where u is a fraction over the "brick primes" of the small strong spectra. This is the same
move set as Conway's PRIMEGAME, whose state lives in the exponents of a few primes. It links to the prior
[PRIMEGAME/metagraph](../metagraph_primegame_20261102/TWO_BIRTH_THEOREM.md) work and to the FRACTRAN word compiler
of [coverage_machine_induction](../coverage_machine_induction_20261103/README.md). The Chan–Pak question becomes:
which register totals x·y are reachable when x and y are strong counts? The bricks 3, 5, 11, 13 are the fraction
alphabet.

## 2. Two integer coordinates of the odd numbers (ELEMENTARY)

The odd numbers are the coset 1 + 2Z. There are two natural integer coordinates on it:
- m = (N − 1)/2, base point 1, so 0 ↦ 1. This is the OCF cycle content: H = 1 + 2m, where m = c1 + 2c2 + 4c3 + …
  sums packings of odd cycles.
- r = (N + 1)/2, base point −1, so 1 ↦ 1. This is the balanced-ternary coordinate of the clone-clock notes, where the
  three branches are r ↦ 3r − 1, 3r, 3r + 1.

The involution N ↦ −N swaps them: m(−N) = −r(N), r − m = 1 and m + r = N. In the ring Z_2, the odd numbers are the
unit group Z_2^× = {±1} × (1 + 4Z_2). This group is the Galois group of the 2-power cyclotomic tower Q(ζ_{2^∞}), in
which the Fermat blocks live (2^j + 1 is a product of cyclotomic values Φ_d(2)). The two coordinates differ by the
sign factor, which is complex conjugation.

**The exceptions in each coordinate.**
- m-gaps {3, 10} = {T_2, T_4} (triangular numbers). Here m = 3 is the universal gap: it is a gap of I(Ω(T), b) in
  every base checked (THM-5070), because no tournament has exactly three pairwise-meeting odd cycles. And m = 10 is a
  gap in every base except 4 and 8 (base 3 misses N = 31); see the chanpak_waves note for the proof at b ≥ 9.
- r-gaps {4, 11}. In balanced ternary these are "11" and "11T": 4 = 3·1 + 1 is the +1 branch that a single vertex
  cannot take, and 11 = 3·4 − 1 is the −1 branch from the dead node 4. In both coordinates the two gaps are 7 apart,
  since 21 − 7 = 2·7.
- Their sums are 3 + 4 = 7 and 10 + 11 = 21.

**A copy of N starting at 11, with eight hidden values in a 3-group and a 5-group.** In m-coordinates the attained
set is N minus {3, 10}. From 11 on it is a full copy of the natural numbers; 11 is r(21) and H = 23. Behind that copy
sit exactly eight attained counts greater than 1, namely 3, 5, 9, 11, 13, 15, 17, 19. Split by divisibility by 3:
- **{3, 9, 15}**, a group of 3: the multiples of the brick 3, with m ≡ 1 mod 3.
- **{5, 11, 13, 17, 19}**, a group of 5: the brick 5 and the primes 11, 13, 17, 19.

The gap 21 lies in the 3-class (m = 10 ≡ 1 mod 3), and the gap 7 in the other class. The eight hidden values, cut
3 + 5, carry the two Busch bricks as their class sizes.

## 3. How these pieces connect

- **Galois coordinate to Fermat blocks.** r and m differ by complex conjugation in Z_2^×. The Fermat blocks
  2^j + 1 are cyclotomic, and 7 = Φ3(2) is the first cyclotomic value lying outside the 2-power tower; that is the
  Fermat-blind prime of the staircase note.
- **Brick primes to FRACTRAN.** The law of §1 says composition collisions are finite-alphabet fraction moves. The
  hidden eight values of §2 are exactly the counts generated before the copy begins.

## 4. Predictions toward proofs (each with its test)

- **P1, brick law for all n.** For two-component compositions (a, b) and (c, d) of n with a ≤ b and c ≤ d, every
  shared count is an S-unit re-association with S = primes(Sp(a)) ∪ primes(Sp(c)). *Test:* n = 13 with
  (4,9), (5,8), (6,7).
- **P2, m = 3 is universal.** No tournament has Ω equal to K3 plus isolated vertices. So 1 + 3b is never
  I(Ω(T), b), for every b ≥ 2. *Test:* prove it via pancyclicity of the strong core (Moon), then check b ≤ 12.
- **P3, the m = 10 pattern.** For b ≤ 12, m = 10 is a gap exactly when b ∉ {3, 4, 8, …}. *Test:* rerun THM-5070's
  chains for b = 9–12 and predict the pattern from b mod 3 and b mod 4.
- **P4, staircase runs.** Strong order 10 contains every odd number in roughly [140, 14000]. *Test:* sample order-10
  tournaments, together with the strong-ear augmentation of THM-4744.
- **P5, shallow cuts.** Sources of order ≤ 10 cover every odd N in [69, 10^7]. *Test:* run `cover_cuts.py` with
  order-10 port pairs computed in C.
- **P6, ten interiors.** The paired-port HN architecture still has no forcing selector at 10 interiors. The first
  survivors, if any, appear at 11 or 12 interiors and contain a Moser-spindle core. *Test:* `interiors_sat.py 10`
  after C-side filtering.
- **P7, the Fermat join sequence.** For every r ≥ 2, Φ3(2^(2^r)) = 1 + F_r·2^(2^r) is a single Fermat-block cut
  whose source is the Fermat block R_(2^r) itself. *Test:* the port spectrum of R_4 and its padding (the case 273).
- **P8, the balanced-ternary dead branch.** In r-coordinates the dead branch 1 → 4 → 11 stops at depth 2.
  All of 11's children (32, 33, 34, i.e. H = 63, 65, 67) are attained through other parents, and no r ≥ 12 is a
  gap. *Test:* this is equivalent to Chan–Pak; check it against the existing 2^32 verification and give a structural
  proof for the subtree of 11.
- **P9, FRACTRAN closure.** The set of products x·y over two-component compositions of order ≤ n, closed under
  brick moves, equals the full order-n non-strong spectrum. *Test:* n ≤ 12 from these data.
- **P10, square sums.** The number of odd square-sum 4-cycles below X grows like c·X·log X (two-square
  representations). *Test:* count to 10^4.
