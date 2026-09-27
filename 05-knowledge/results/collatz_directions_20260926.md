# Directions toward a Collatz proof: the no-descent fractal, why no prefix rank exists, the lowest-set-bit map, the two-place series, the Mahler sibling, and irrational versus dyadic hierarchies

**Session:** opus, `collatz-directions-20260926`, 2026-09-26.
**Owner's directive:** "consider the snippets for directions to head towards
a Collatz proof and prioritize creative mathematical reasoning, not
re-verification": (a) the reset lane's obstruction, "a rank across successive
excursions that survives regenerated precision and pays the original source
threshold; the strongest obstruction involves two complete excursions, both
above the original source, with arbitrarily large combined growth; the
construction regenerates arbitrarily deep precision around the same negative
cycle, so tracking the cycle label cannot resolve it"; (b) the forest tiling
results (`17` polygon multisets, `21` cyclic stars with counts `1, 3, 7, 10`
by valency `6, 5, 4, 3`; six stars that cannot extend; periodic tilings with
two stars and unbounded orbit count `k_N`; the Sierpinski selection
`Q = 2Q + {00, 10, 01}` with `3^r` centres per `2^r x 2^r` block and unbounded
hierarchical depth; `F = U + S` with nested versus disjoint classes; the
`14`-edge refinement graph).
**Inherits:** [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)
(thin divergence, summable reciprocals, the two-place corollary on the
Bernstein series), [THM-4495](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md)
(exact order of the no-descent counts), [THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md),
the S11/S12 notes (residual, shadows), the repository's Mahler `3/2`
frontier ([`MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md`](../../00-navigation/MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md):
THM-3848, THM-4072, THM-4074, THM-4077, THM-4082), the forest tiling note
(`forest_20260926_tilings.md`), and the S10 Fermat-tower note.

**Status: PROVED small propositions (the no-descent fractal `E_inf` is closed
with box dimension `h*` and is the closure of the fixed points of the
no-descent words, the height-minimal points of the negative rational cycles;
no rank that is a function of the parity prefix can decrease at every
excursion; the lowest-set-bit conjugacy) + CITED (the reciprocal-sum form of
no-divergence is THM-4476's Corollary 8) +
FINITE-EXACT (`E_inf` meets the negative integers exactly in `{-1, -5, -17}`
up to `10^6`; the lowest-set-bit map reaches a power of two for every
`n <= 10^5`; dimension estimates to `k = 2000`; the lower approximations of
`log_2 3` and their gaps) + DIRECTIONS marked as such + INDEPENDENTLY AUDITED
(HAS GAPS -> repaired, section 7: wording overclaims, no numerical error).
Collatz OPEN; nothing here is a proof.** Script: `04-computation/experiments/collatz_directions_20260926.py`
-> `.out`.

## 0. The answer in one paragraph

No rank that is a function of finite precision can decrease at every
excursion, and this is a one-line theorem (section 1, Proposition 4); the
lane's controller carries the full integer in its registers, so it is not in
that class, and the proposition says what it must use, not that it fails: the set `E_inf` of 2-adic integers whose parity
word never has a coefficient descent is a closed set of box dimension
`h* = 0.94996` in `Z_2`, it is the closure of the negative rational cycle
points at their height-minimal rotation (`-1, -5, -17` and the fixed points
`x_w` of every other no-descent word; `-7` and `-25` are cycle points outside
`E_inf`), every one of its cylinders contains positive integers, and the
conjecture implies `Z^+ ∩ E_inf = empty`, with the converse under Terras's
equality `sigma = sigma_inf` (Proposition 3). "Regenerated precision around the
same negative cycle" is the density of cycle points in `E_inf`; "two
excursions with unbounded combined growth" is the positive dimension. So a
rank must use the whole integer, never a prefix alone, and the only such
rank we know is the stopping time itself. What is left are global
inputs of three kinds, none of which is a rank: (i) an Archimedean input on
the real shadow `xi = n C_inf` of a divergent orbit, whose adaptive Mahler
sequence `xi 3^l / 2^(d_l)` sits within `eta_l >= 1/3` of the odd iterates
(section 4; the repository's Mahler frontier is stuck at exactly the same
wall, "couple the exact residue classes to an Archimedean or genuinely
all-depth `2`-adic obstruction"); (ii) a transversality input distinguishing the finite binary
expansions of integers inside `E_inf`, which the lowest-set-bit form
`r -> 3r + lsb(r)` makes explicit (section 2: the conjecture says this map
always reaches a power of two); (iii) the two-place series (section 3: the
conjecture's no-divergence half is "the reciprocal sum of every orbit
diverges", THM-4476's Corollary 8). The tiling snippets describe a dyadic,
self-similar hierarchy (the Sierpinski selection is the Lucas set of odd
binomials, the Gilbreath tower of S10); the Collatz hierarchy is the
continued fraction of `log_2 3`, irrational, with integer cycles only at the
first three lower approximations `1/1, 3/2, 11/7` (gaps `1, 1, 139`): the
heuristic reading (section 5, marked) is that no finite family of certificate
shapes tiles the residual because the self-similarity that would close the
tiling does not exist.

## 1. The obstruction, exactly (PROVED small propositions)

Write the parity word of `x in Z_2` in `T`-coding (`T(x) = x/2` or
`(3x+1)/2`), and its height after `j` steps `h_j = o_j log_2 3 - j` with `o_j`
the number of odd steps (coefficient no-descent at every step iff `h_j >= 0`
for all `j`). Let

```text
E_inf = { x in Z_2 : h_j(x) >= 0 for every j >= 1 }.
```

**Proposition 1 (closed, dimension `h*`).** `E_inf` is closed; the number of
length-`k` cylinders meeting it is exactly THM-4495's `W_k = Theta(2^(h* k)
k^(-3/2))`, so its box-counting dimension in the 2-adic metric is
`h* = h(log_3 2) = 0.9499555`. *Proof.* Closed: each condition `h_j >= 0` is a
condition on the first `j` letters, i.e. on a finite union of cylinders. A
length-`k` word with all prefix heights `>= 0` extends to an infinite such
word (repeat it: `h_(mk + i) = m h_k + h_i >= 0`), and every infinite parity
word is the word of exactly one point of `Z_2` (the parity-vector bijection
of Terras, Lagarias and Bernstein; the finite form is THM-4476 step 1), so
the cylinders meeting `E_inf` are exactly the no-descent words, counted by
`W_k`; the dimension is
`lim log_2 W_k / k = h*`. ∎ Numerically `log_2 W_k / k = 0.8496, 0.9084, 0.9384,
0.9435` at `k = 60, 200, 1000, 2000`, and `log_2 W_k / k + 1.5 log_2 k / k =
0.9973, 0.9657, 0.9534, 0.9517`, both converging to `0.94996` at the rate the
polynomial factor dictates.

**Proposition 2 (the negative copy is dense).** The rational points of
`E_inf` include the fixed points `x_w = S_w/(2^A - 3^p)` of every no-descent
word `w` (all negative), and these are dense in `E_inf`: the cylinder of any
finite no-descent word `u` contains `x_u` (the point with word `u^infinity`
is fixed by `T^k`, hence equals `x_u`; `u^infinity` is no-descent by the
extension argument). `E_inf` is therefore the closure of these fixed points,
which are the height-minimal points of the negative rational cycles; the
other points of a cycle (`-7`, `-25`, ...) are not in `E_inf`. Among the
negative integers the members of `E_inf` with `|n| <= 10^6` are exactly
`-1, -5, -17` (script part 4, by the criterion "the orbit never decreases in
absolute value"; the audit confirmed the list by the height criterion
directly). For negative integers a coefficient descent is an actual descent
(`S_j > 0`), and the converse holds step by step (`|U(m)| < |m|` iff
`v >= 2`); the multi-step converse was verified for `|n| <= 10^6`, not proved.

**Proposition 3 (the conjecture as transversality).** `Z^+ ∩ E_inf = empty`
iff every positive integer has finite coefficient stopping time; with
Terras's equality `sigma = sigma_inf` (verified to `10^7`, THM-4512) this is
the conjecture. Positive integers are dense in `Z_2`, so every cylinder of
`E_inf` contains positive integers: at every finite precision the integers
are indistinguishable from the fractal.

**Proposition 4 (no bounded-precision rank).** There is no function `Phi` of
finite parity prefixes, with values in a well-ordered set, that decreases at
every complete excursion (a rise followed by a local descent) of every
no-descent word: the words `(1, 1, 2)^N` are no-descent (prefix heights `0.585, 1.170,
0.755` within a period, period height `3 log_2 3 - 4 = +0.755`) with `N`
complete excursions (two rises, then `v = 2`, a local descent that stays
above the start), so
`Phi(empty) >= N` for every `N`. Precisely, the excluded class is: functions of the finite parity prefix,
with well-ordered values, strictly decreasing across every rise-then-local-
descent of every no-descent word. Not excluded (audit): ranks that depend on
the integer `n` (the lane's controller carries `k, r, b` with `b = (n+5)/2^L`,
the whole integer, so its rank is of this kind); ranks required to decrease
only at excursions that return below the *original source*, which are
vacuous on `E_inf` and are the closest reading of "pays the original source
threshold"; ranks allowed to increase between excursion ends are excluded.
If "complete" is read as "returns below its own base", the witness
`(1, 1, 1, 1, 2, 3)^N` (heights `0.585, 1.170, 1.755, 2.340, 1.925, 0.510`,
period `+0.510`, fixed point `-697/217`) replaces `(1, 1, 2)^N`. So a prefix
rank cannot pay the source; a rank on the integer must, and for every `K`
there are `n = n' mod 2^K` with different rank: this is one reading of the
lane's "regenerated precision" and "two excursions with unbounded combined
growth", and it is why prefix certificates can only ever cover residue
classes (the density sieve, whose residual has THM-4495's order, S11).

## 2. The lowest-set-bit map (PROVED conjugacy, FINITE-EXACT)

Let `lsb(r)` be the lowest set bit of `r` as a power of two and
`rho(r) = 3r + lsb(r)`. Starting from `r_0 = n` odd,

```text
r_l = 2^(d_l) m_l = 3^l n + S_(l-1),     m_l = the Syracuse orbit,  d_l = v_1 + ... + v_l,
```

because `3 (2^d m) + 2^d = 2^d (3m + 1) = 2^(d + v) U(m)` (checked for all odd
`n <= 2000`, `l <= 200`). So the Collatz map is "triple and add your lowest
bit", without ever dividing; the orbit reaches `1` iff `r_l` is eventually a
power of two (then `rho(r) = 4r` for ever; every `n <= 10^5` gets there, at
most `129` steps; `27` reaches `2^70` after `41` steps, its `41` odd steps and
`70` halvings), and diverges iff the odd part `r_l / lsb(r_l)` tends to
infinity, in which case `r_l ~ xi 3^l` with `0 < xi < infinity` (full rate,
THM-4476) while `lsb(r_l) = 2^(d_l)` with `d_l >= l`. The carry is the
`3`-weighted sum of the orbit's own past lowest bits, `S_(l-1) = sum_(t<l)
3^(l-1-t) lsb(r_t)`, and the divergence question becomes: **can `3^l n` plus
the `3`-weighted sum of its own past lowest set bits stay divisible by a power
of two within a subexponential factor of itself, for ever?** The two copies
are both visible: `r_l = S_(l-1) mod 3^l` (the `3`-adic side reads the carry)
and `r_l = 0 mod lsb(r_l)` (the `2`-adic side reads the halvings). This is the
natural home of the lane's carry-defect potential; it says nothing new by
itself, but it puts the "regenerated precision" on the table as the size of
`lsb(r_l)` relative to `r_l`.

## 3. The two-place series (CITED + corollary)

THM-4476's "two places" corollary: for the word `d` of a divergent positive
orbit the real series `R(d) = sum_l 2^(d_l)/3^(l+1)` converges, while its
2-adic value is `R_2(d) = -n`; for the `3n-1` sheet (the negative copy) with
`R_2(d) = n > 0` one has `R(d) < n`. In the notation of this thread,
`S_(l-1)/3^l = n (C_l - 1)` with `C_l = prod_(j<l) (1 + 1/(3 m_j))` (checked
exactly on `27`, `-1`, `-5`): for an orbit reaching `1` the product diverges
and so does the real series (`27`: partial sums `0.33, 0.56, 0.85, ..., 2.14`
and growing); for the cycle points the real and 2-adic values coincide
(`-1`: `0.9999...`; `-5`: `3.63, 3.72, 3.78, ... -> 5`); for a divergent
orbit the real value would be `n (C_inf - 1)` with `C_inf < infinity`.

**Corollary (reciprocal-sum form).** The no-divergence half of the conjecture
is: `sum_j 1/m_j(n) = infinity` for every positive integer `n`. (Divergent
implies summable by THM-4476; summable implies `m_j -> infinity`.) So the
conjecture asks that the Bernstein series of every positive integer's word
diverge in `R` while converging 2-adically to `-n`; for the rational cycle
points the two completions agree, and a divergent integer orbit would be a
word whose two limits are a positive real and a negative integer. This is
the exact form of the local-global gap named in S11.

## 4. The Mahler sibling and the Archimedean object (DIRECTION)

Mahler's `Z`-numbers (`{xi (3/2)^n} < 1/2` for all `n`) reduce to the integer
map `g -> ceil(3g/2)` with a residue condition: the halving is rigid (one per
step) and the real number `xi` is recovered from the integer orbit. A
divergent Collatz orbit is the adaptive-halving sibling: with
`xi = n C_inf = lim 2^(d_l) m_l / 3^l`,

```text
xi 3^l / 2^(d_l) = m_l + eta_l,     eta_l = m_l (C_inf / C_l - 1) >= 1/3,
```

so the real sequence `xi 3^l / 2^(d_l)` sits within `eta_l` of the odd
integers `m_l`, with the exponents `d_l` chosen by the orbit itself. The
repository's Mahler frontier is stuck at the same wall as the reset lane:
THM-4072 (finite-state obstruction), THM-4074 (arbitrarily long reset
runways after which every finite carry word beginning in `1` can be
programmed, in the denominator-19 family: the exact analogue of "regenerated
precision"), THM-4077/THM-4082 (a 2-adic tangent
isometry that produces no `Z`-number), and its live task reads "couple the
exact residue classes to an Archimedean or genuinely all-depth 2-adic
obstruction". The safe-tail shift there has binary-ultrametric Hausdorff
dimension `log_2(3/2) = 0.585` (THM-3848); the Collatz no-descent fractal has
`h* = 0.950`; the two are different objects (THM-3848's shift is not the
unknown `Z`-language), so the comparison is a loose analogy, and Mahler's
countability theorem (at most one `Z`-number per unit interval) has no
Collatz content because the Collatz starting points are already integers.
The direction: the only all-`n` structural theorems about `xi (3/2)^n` that
exist (Flatto-Lagarias-Pollington: the fractional parts spread over an
interval of length at least `1/3`; Dubickas) are Archimedean digit-spreading
arguments that survive arbitrary depth. A Collatz version would bound the
spread of `{eta_l}` or of the shadow error from below along any divergent
orbit and contradict the convergence `C_l -> C_inf`. I do not know how to do
it; it is where the lane's rank must be replaced by a global invariant.

## 5. Irrational versus dyadic hierarchies (the tiling snippets; DIRECTION, heuristic)

* The Sierpinski selection `Q = 2Q + {(0,0), (1,0), (0,1)}` is the set of
  `(i, j)` with `i & j = 0`, i.e. the odd binomial coefficients (Lucas), the
  Pascal-mod-2 support: `3^r` centres per `2^r x 2^r` block is the
  `3^(K-1-m)` count of the Gilbreath tower (S10). The owner's phrase "any finite collection of bounded-radius observations
  can remain fixed while the distance to the nearest hexagon grows
  unboundedly" (the forest note says it as "any finite family of
  bounded-radius local features can be fixed while the needed global address
  information remains unbounded") is the tiling form of Proposition 4:
  local data cannot see the hierarchy. The two-star alphabet `3^6, 3^4.6`
  with unbounded orbit count `k_N` is the tiling form of "two letters
  (`v = 1` rise, `v >= 2` fall) with unbounded global structure".
* The dyadic hierarchy closes because `2 Q + {...}` is exact
  self-similarity. The Collatz hierarchy is the continued fraction of
  `log_2 3`: best lower approximations `A/p = 1/1, 3/2, 11/7, 19/12, 84/53,
  569/359`, gaps `3^p - 2^A = 1, 1, 139, 7153, 4.0e22, 2.1e168`. Among words with `p <= 12`, integer negative cycles exist only at the
  first three (gap `1` twice, by `3 - 2 = 1` and `9 - 8 = 1`, the only
  solutions of `3^p - 2^A = 1`, and the sporadic `139 | S`), none at `19/12`
  (S12 audit); whether any exist at `84/53` or beyond is not settled by
  anything cited here. The shadows of the other lower approximations are
  rational points of `E_inf` with denominators growing like the gap.
  Heuristic, not a result: because `log_2 3` is irrational, no depth of the
  hierarchy repeats the previous one exactly, so no `Q = 2Q + F` describes
  the Collatz residual, only an infinite sequence of ever-worse rational
  approximations; the rigorous content behind "no finite bank covers `Z^+`"
  is only that `W_k >= 1` for every `k` (the word `1^k`), which does not use
  irrationality.
* The counts `1, 3, 7, 10` and `F = U + S` (`10, 7, 3, 1`, nested versus
  disjoint) have no Collatz reading I can defend beyond the Mersenne
  coincidence `1, 3, 7`; the `14`-edge refinement graph is a tiling fact.
  I record them as the owner's, unlinked.

## 6. What a proof would have to contain (ranked, honest)

1. **An Archimedean input** on `xi = n C_inf` (section 4): the one family of
   all-`n` arguments about `3^l` and powers of two that exists in the
   literature. Hardest, most promising in shape.
2. **A transversality input**: a property of finite binary expansions that
   the recursion `r -> 3r + lsb(r)` cannot preserve for ever inside `E_inf`
   (section 2). Nothing known separates integers from the other points of
   `E_inf` except finiteness of the expansion, and that interacts with the
   carry only through the recursion itself.
3. **Diophantine inputs** handle cycles (`2^A - 3^p | S`, Steiner,
   Simons-de Weger) and have no divergence analogue: the gap `2^A - 3^p` is
   replaced by a moving target.
4. **Density and dimension** are done up to constant factors (THM-4495,
   THM-4498, THM-4499, THM-4476) and cannot be pushed to the pointwise statement
   (Propositions 3-4).

Proposition 4 says what the lane's "rank across excursions" must be built
from (the integer, which its controller already carries), not that it fails;
the owner's two-excursion construction (an external snippet, not in the
repository) is one instance of the same phenomenon, and the note's own proof
uses `(1, 1, 2)^N`. The shell-revisit statement of S9/S10 (HYP-9161 is
REFUTED as stated for arbitrary finite segments; its whole-orbit form is
OPEN) remains the one statement in this thread that is pointwise and not
obviously in that class, because it constrains the orbit of one integer
through its own address, not through a prefix.

## 7. Independent audit (2026-09-26)

Auditor subagent: `04-computation/experiments/collatz_directions_20260926_audit.py`
(sha256 `020264cee1f5b76a2cde3fe0ec1089ac6782e859b374626509c6a808dd41523b`) ->
`05-knowledge/results/collatz_directions_20260926_audit.out`
(sha256 `a67fc04393c71134e9d703c0cae2048fcc37d4e7be3fa10c4b8c51cdaff3fc66`; 53 checks, all pass).
CONFIRMED: every number (the `W_k` values and dimension estimates, the
Spitzer identity to `k = 300`, the fixed points of all no-descent words to
`k = 12`, `E_inf` meets the negative integers in `{-1, -5, -17}` by the height
criterion, the lowest-set-bit identities, `129` steps at `n = 77031`, the
lower approximations and gaps, `139 | S` for `-17`, the two-place identities
to `l = 60` and the exact cycle values, Terras's equality `T`-coded to
`10^6`); the four propositions under the note's definitions; the citations.
CORRECTED (applied above): the Terras caveat dropped from section 0 and the
wiring; the claim that the lane's rank is in the excluded class (its
controller carries the whole integer; Proposition 4 excludes prefix ranks
only, and the precise excluded class is now stated); "closure of the
negative rational cycle points" (only their height-minimal rotations, i.e.
the fixed points of no-descent words); the two-excursion snippet presented
as a proof; integer cycles beyond `p <= 12` unsupported; HYP-9161's REFUTED
status; the reciprocal-sum form relabelled CITED (THM-4476 Corollary 8); the
dimension comparison with THM-3848's shift labelled a loose analogy; the
irrationality argument labelled heuristic; `THM-4074` "beginning in `1`";
two quotations; the "optimal in order" gloss; the parity-bijection citation;
the multi-step converse for negative integers marked verified, not proved.
Verdict: HAS GAPS -> repaired; no numerical error.
