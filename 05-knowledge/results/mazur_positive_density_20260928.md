# Mazur's positive-density log-time convergence, digested: the harmonic mass of a Syracuse tree layer is `3^n mu_n(seed)`, the negative cycles are the resonances of the 3-adic Syracuse law (exact spike profile), and the seed-1 mass is a computable test of the theorem (`limsup n^(1/6) H_n(1) > 0` is forced; `H_18(1) = 0.41` and falling)

**Session:** opus, `collatz-poset-dag-20260927` (S19 continuation), 2026-09-28.
**Owner's directive:** "thoroughly digest the attached pdf for any possible
connections to our work we can leverage or extend towards proofs, and spend
a long session doing as much math as you can"; "keep heading in these
directions and new ones you synthesize"; Kaprekar's `6174`, `495`, the
2-digit cycle `9, 81, 63, 27, 45`, the perfect cuboid; the pasted "inverter"
text (valuation sequences to 2-adic values; periodic implies rational; open
whether an unrepeating sequence can have a rational value; a natural-number
value would refute Collatz).
**The paper:** M. Mazur, *Explicit Positive-Density Collatz Convergence in
Logarithmic Time*, v2, September 2026 (AI-assisted; ProofAtlas platform;
accompanying Lean 4 development). Companion (same platform, same supporting
source revision): *Positive lower density of Collatz predecessors*.
**Inherits (cited):** Tao 2022 (Syracuse random variables, fine-scale mixing,
Proposition 1.14; from memory), Terras 1976 (the bijection between residue
classes mod `2^k` and parity words), Krasikov–Lagarias 2003 (`x^0.84`
predecessors of `1`; from memory), Lagarias 1985 / Bernstein–Lagarias 1996
(the 2-adic conjugacy map and the Periodicity Conjecture; from memory),
THM-4476 (thin divergence; its corollaries (2), (3), (8) on the Bernstein
value), THM-4514 (spine blocks, two-place descent tree), the barrier atlas
(`collatz_procgen_20260922_barrier_atlas.md`, typing rule of §0), the STICKY
note (S16), the Artin-corrections note (S18: residue mod `3` of a visited
value is the parity of the arriving valuation).

**Status: CITED (Mazur's Theorem 1.1 and the companion theorem: Lean-checked
on the platform, review pending, no independent replay by us, automated
reviews only) + PROVED (Theorem A, the harmonic-mass identity and its
recursion; Theorem B, cycle resonances, lower-bound form; Theorem C, the
seed-1 test) + FINITE-EXACT (the 3-adic Syracuse law to level 18, rational
to level 6; every elementary component of the paper checked at small size,
P1–P8) + OBSERVED (the exact spike profile `D` to four decimals at levels
17 and 18; the decline of the seed-1 mass; linear growth of the second
moment; convergence of the entropy deficit) + CONDITIONAL (absolute
continuity of the 3-adic Syracuse law, from the fine-scale estimate as
stated) + OPEN (`lim H_n(1) > 0`?) + DIRECTION. Audit: pending (section 11
will hold the record).**

Scripts and outputs (all in this session's worktree):
`04-computation/experiments/mazur_positive_density_20260928.py` (P1–P8),
`mazur_harmonic_mass_20260928.py` (exact rational law, tree layers),
`mazur_harmonic_mass_deep_20260928.py` (FFT law to level 18; run with the
level as argument), `mazur_seed1_test_20260928.py` (Theorem C's numbers),
`kaprekar_cuboid_20260928.py`; outputs `05-knowledge/results/mazur_*.out`,
`kaprekar_cuboid_20260928.out`.

---

## 0. What is new, in one screen

1. **The paper's object is the repo's object.** Mazur's reference density
   `rho_q = (2/3) 3^q mu_q` is the level-`q` law of Tao's Syracuse random
   variable, and his weighted inverse histories from a seed `M` with weights
   `omega(w) = 3^d 2^(-A(w))` are exactly the *harmonic mass* of the
   Syracuse tree of `M`. **Theorem A (PROVED):** for every odd integer `y`
   not divisible by `3` (either sign) and every `n >= 1`,
   `3^n mu_n(y mod 3^n) = sum_{x odd : S^n(x) = y} 3^n 2^(-A_n(x))`,
   the sum over the depth-`n` predecessors of `y` (which have the sign of
   `y`). So the abstract transfer-operator quantity at the residue of the
   seed *is* the weighted count of the seed's actual inverse histories; the
   paper's "seed by residue averaging" is the choice of the tree with the
   largest such mass among the trees of `R_j = (4^j - 1)/3`.
2. **The negative cycles are the resonances of the 3-adic Syracuse law
   (Theorem B, PROVED).** A Syracuse cycle of length `k`, total valuation
   `A`, through `y_0` gives `mu_{km}(y_0 mod 3^{km}) >= 2^(-Am)`, i.e.
   `rho >= (2/3)(3^k/2^A)^m`. Negative cycles have `3^k > 2^A`: the class
   of `-1` carries `rho_n(-1) = 0.9748 (3/2)^n` (the word `1^n` lands there
   exactly: `(3/2)^n - 1 = -1 mod 3^n`), the `{-5, -7}` cycle grows like
   `(9/8)^(n/2)`, the seven-cycle of `-17` like `(2187/2048)^(n/7)`. The
   limit density has singularities of order `|y - y_0|_3^(-s)`, `s = 1 -
   (A/k) log_3 2 = 0.3691, 0.0536, 0.0086`. Positive cycles (`2^A > 3^k`)
   are anti-resonant: the trivial cycle contributes `(3/4)^n`. The relative
   profile of the `-1` spike on its forward rational closure is exact:
   `D(-1/2) = 1/2`, `D(-1/4) = D(1/8) = D(11/16) = D(49/32) = 3/4`,
   `D(-1/8) = D(1/16) = D(11/32) = 3/8`, `D(-1/16) = 3/16`, all reproduced
   to four decimals at levels 17 and 18 (path sums `2^(l - A)`).
3. **The seed-1 mass tests the theorem (Theorem C, PROVED).** With `H_n(1)
   := 3^n mu_n(1 mod 3^n)`: positive lower density of `{x : tau(x) <= C ln
   x}` for *any* `C` forces `limsup_n n^(1/6) H_n(1) > 0`. The exact
   sequence is `1, 8/7, 1376/1387, 0.928, 0.964, 0.955, 0.860, 0.775,
   0.697, 0.637, 0.591, 0.537, 0.500, 0.473, 0.459, 0.435, 0.421, 0.410`
   (`n = 1..18`), falling about `3%` per level since level 8, with
   `n^(1/6) H_n = 0.665` at `n = 18`. A positive limit near `0.35` is
   consistent with the data; so is decay to `0`. The paper does not claim
   anything about the seed `1`; its residue averaging is precisely what
   avoids this quantity. This is the sharpest computable stake we have on
   the theorem's *conclusion*, independent of its astronomical constants.
4. **The 3-adic Syracuse law looks absolutely continuous but not
   square-integrable.** Its entropy deficit `n ln 3 - H(mu_n)` converges
   (`1.06, 1.08, ..., 1.173` at `n = 10..18`, increments shrinking
   geometrically; information dimension `1`, `E[rho log rho] -> 0.84` on
   units), while `E[rho_n^2]` over units grows linearly (`0.31 n + 0.8`).
   The fine-scale distances `||mu_18 - lift(mu_m)||_1` decrease in `m`
   (`0.85` at `m = 1` to `0.35` at `m = 9`), the qualitative content of
   Mazur's (2.3) and Tao's Proposition 1.14, whose explicit constant here
   is `2^(2^6536)`-sized and says nothing at these levels. As stated, (2.3)
   *is* the `L^1`-Cauchy property of the densities, hence absolute
   continuity of the 3-adic law (CONDITIONAL on (2.3)).
5. **Typing.** Mazur's theorem is a residue-averaging mechanism: it
   overcomes DRIFT and STICKY (size-free Terras weights, exact inverse-orbit
   sums at every scale), and is blind to SHEET, DEFECT, DIM (positive
   density says nothing about the exceptional set; the companion theorem
   gives every target `a`, `3 ∤ a`, a positive-density basin, so extra
   cycles would be invisible to it). It replaces the Krasikov–Lagarias
   `x^0.84` row by positive density, *if* it stands.
6. The owner's seeds: the "inverter" is the Bernstein 2-adic map and the
   open question is the Periodicity Conjecture, on which THM-4476 already
   has the partial result (aperiodic words of slow growth have no
   odd-denominator rational value); Kaprekar's routine has the exact
   invariant of the nine-iteration typology (the digit multiset, a finite
   state space), which is why `495`, `6174` and the 2-digit cycle `9 ->
   81 -> 63 -> 27 -> 45` are theorems by enumeration; the perfect cuboid
   remains open (Euler bricks `(44,117,240)`, `(240,252,275)` below 300).

---

## 1. The paper: statement, platform status, what was and was not verified

**Theorem 1.1 (Mazur).** The rational `c` and the integer `X_0` defined in
Section 8 are positive and, for `X >= X_0`,
`#{n < X : tau(n) <= (523/50) ln n} >= c X`,
where `tau(n)` is the number of steps of `n -> n/2` (even), `3n+1` (odd) to
reach `1`. Corollary 1.2: the same with `262/25 = 10.48`. The constant
`523/50 = 10.46` is the central-word time constant `3/log(4/3) = 10.428`
(one Syracuse step of mean valuation `2` costs three Collatz steps and
shrinks by `3/4`) with rounding losses.

**The formal statement, as read on the platform** (main file, commit
`830b9d3f`): `step`, `fastConvergent C n := 0 < n ∧ ∃ k, k <= C log n ∧
step^[k] n = 1`, `countBelow C X := Nat.count (fastConvergent C) X`, and
`positive_density (X) (hX : cutoff <= X) : densityConstant * X <=
countBelow (523/50) X`, with `densityConstant_pos`, `cutoff_pos`, and the
`262/25` variant. This matches Theorem 1.1 exactly (count of positive
integers strictly below `X`; natural logarithm; the map (1.1)).

**Platform status (read 2026-09-28):** "Lean checked, recorded build
passed, unfinished proof steps: none"; publication "Review pending"; 386
files, 61,169 lines; toolchain `leanprover/lean4:v4.30.0-rc2`, mathlib
`5450b53e`; axioms `propext`, `Classical.choice`, `Quot.sound`; the
platform states that its public hosting does not accept the Lean result as
such; the reviews shown are by AI systems (the platform names a Claude model)
with independence not asserted. The companion page (*Positive lower density
of Collatz predecessors*) is marked "Accepted formalization", 388 files /
61,800 lines, pinned on the same supporting revision, with theorems
`predecessors_positive_lower_density`,
`nonconvergence_positive_lower_density_of_counterexample`,
`universal_reaches_one_of_arbitrarily_dense_convergence`.

**Not done by us:** an independent replay of the Lean build (the source
package is a download, which this session does not initiate without the
owner), a check that the Lean development's internal statement of the
mixing estimate (2.3) is what the paper says, and any human review. The
paper itself says (Section 9) that the formal development "includes a
proof of the finite-distribution mixing estimate (2.3)" and that the
correspondence between the exposition and the declarations is "provided in
the accompanying theorem-specific source package". Our typing is therefore
**CITED, review pending**: the elementary layer is verified below, the
analytic core (2.3) with its numerical coefficient (Section 8.1) is taken on
the platform's word.

---

## 2. The argument's skeleton

1. **Finite residue spaces and transfer operators.** `G_t = Z/3^t`. For a
   word `w = (a_1..a_d)` the affine inverse map induces an injection `F_w :
   G_t -> G_{t+d}`, `F_w(z) = 3^d 2^(-A) z + sum_j 3^(j-1) 2^(-(a_1+...+a_j))`,
   and the weighted transfer `(T_w g)(y) = 3^d 2^(-A) sum_{F_w(z) = y} g(z)`.
   **Lemma 2.1:** `<|T_w g|>_{t+d} = 2^(-A(w)) <|g|>_t`, with or without
   absolute values (the fibre is a single point; the `3^d` cancels).
2. **Reference density.** `mu_q` = law of Tao's Syracuse random variable on
   `G_q` (`Y_0 = 0`, `Y_{k+1} = 2^(-A)(3 Y_k + 1)`, `A` geometric with
   `P(A = a) = 2^(-a)`); `rho_q = (2/3) 3^q mu_q`, `<rho_q> = 2/3`, zero on
   nonunits.
3. **Fine-scale mixing (2.3):** `||rho_q - rho_m ∘ pi_{q,m}||_q <= (2/3)
   C_A m^(-A)` for `1 <= m <= q`, "the finite-distribution form of the
   fine-scale estimate in Tao's analysis [Proposition 1.14]"; used with the
   fixed exponents `6` and `2`; Lemma 8.1 gives the closed coefficient `C`
   of Section 8.1 for exponent `6`.
4. **Deterministic residue spread (3.9):** distinct endpoints of bounded
   histories in one class mod `3^q` are spaced at least `3^q` apart; equal
   endpoints share their word; so the number of central histories landing
   in a class is bounded by a polynomial in the generation.
5. **Seed by residue averaging (Proposition 4.2):** the seeds `R_j = (4^j -
   1)/3` (`S(R_j) = 1` with valuation `2j`) permute the residues mod `3^q`
   as `j` runs over `0..3^q - 1`; averaging the weighted terminal sums over
   the residue classes and picking a class at least as good as the mean
   gives a fixed odd seed `M >= 16 b_0`, `3 ∤ M`.
6. **Terminal selection (Proposition 5.2) and the comparison of two residue
   moduli (Proposition 5.3),** then **coverage (Proposition 6.3):** the
   source charge `x omega(w) <= M` (each history's source times its weight
   is at most the seed), so a terminal weighted sum bounded below at scale
   `X` is a count of distinct sources below `X`.
7. **Constants (Section 8):** `A* = 6409`, `E* = 2170`, `L* = 280`, then
   towers: `D_exp = 2^(8192 A* 2^(3 E*))` (about `2^(2^6536)`), `C* = (32
   A* D*)^A*` (`log_2 log_2 C* ≈ 6548`), `N = 20000(ceil(log_2 F) + 64)`
   with `F = 2467 b 16^b (C+1)` (`log_2 N ≈ 6563`), and `c^(-1)` beyond
   `2^(2^(2^6535))`. The paper says these are "an explicit certificate, not
   a claim of numerically usable constants".

---

## 3. Elementary components verified (`mazur_positive_density_20260928.py`)

| part | claim of the paper | check | result |
|---|---|---|---|
| P1 | seeds `R_j = (4^j-1)/3` are odd, `S(R_j) = 1` with valuation `2j`; `v_3(4^j - 1) = 1 + v_3(j)`; `R_0..R_{3^q-1}` permute `Z/3^q` | exact, `q <= 5` | holds |
| P2 | inverse history `R_{i+1} = (2^{a_{i+1}} R_i - 1)/3`; source `x = (2^A R - C)/3^d`, `C = sum 3^(j-1) 2^(A - (a_1+..+a_j))`; `3^d x < 2^A R` | 252 histories | holds |
| P3 | Lemma 2.1 with and without absolute values; existence of an integral source iff `R mod 3^d ∈ Im F_w` | 200 cases; `R < 3000`, words `<= 3`, valuations `<= 4` | holds |
| P4 | `rho_q`: mean `2/3`, zero on nonunits | exact rationals, `q <= 5` | holds; `ell^1` distances between consecutive lifts `0.317 .. 0.198` (slow at tiny levels) |
| P5 | mean identity `<sum_w T_w rho_k> = (2/3) p(V)` (`p(V) = 225/256`); residue pullback | seed `M = R_41`; the composed value `4.7143` reproduced exactly | holds |
| P6 | endpoint spread (3.9): same class mod `3^q` implies spacing `>= 3^q`; equal endpoints share the word | enumeration | holds |
| P7 | Section 7: `tau(x) = d + A + tau(M)`; `A log 2 - d log 3 <= log(x/M) <= A log 2 - d log 3 + sum_j log(1 + 1/(3 x_j))`; the constant `3/log(4/3)` | 300 + 300 sampled histories from unit seeds | identities exact; `tau/log x` means `5.35 / 4.81` (seed-diluted); seed-free `(d+A)/(A log 2 - d log 3)` means `6.19 / 5.40` at mean valuations `2.39 / 2.56`; `3/log(4/3) = 10.428` at `A = 2d` |
| P8 | the constants | evaluated | `D_exp ≈ 2^(2^6536)`, `log_2 log_2 C* ≈ 6548`, `log_2 N ≈ 6563`, `c^(-1) > 2^(2^(2^6535))` |

Empirically `54.22%` of `n < 10^6` satisfy `tau(n) <= 10.46 ln n`; the
theorem's lower bound `c` is the reciprocal of a triple exponential.

---

## 4. Theorem A: the harmonic mass of a Syracuse tree layer

**Notation.** `S(x) = (3x+1)/2^{v_2(3x+1)}` on odd integers; `A_n(x)` the
sum of the first `n` valuations; `T_n(y) = {x odd : S^n(x) = y}` (all
depth-`n` predecessors, loops included); `T_n^first(y)` those with first
arrival at depth `n`; `mu_n` the level-`n` law of the Syracuse random
variable; `rho_n = (2/3) 3^n mu_n`.

**Theorem A (PROVED).** For every odd integer `y` with `3 ∤ y` and every
`n >= 1`:
`H_n(y) := sum_{x ∈ T_n(y)} 3^n 2^(-A_n(x)) = 3^n mu_n(y mod 3^n) = (3/2) rho_n(y)`,
and every `x ∈ T_n(y)` has the sign of `y`. If `3 | y` both sides are `0`.

*Proof.* Fix a word `w` of length `n`, `A = A(w)`. (i) The odd integers
whose first `n` Syracuse valuations are `w` form one residue class mod
`2^(A+1)` (Terras; induction on `n`: `v_2(3x+1) = a` is one class mod
`2^(a+1)`, and `x mod 2^(a_1+...+a_n+1)` determines `x_1 mod
2^(a_2+...+a_n+1)` bijectively). (ii) For such `x`, `3^n x + C_w = 2^A
S^n(x)` with `C_w = sum_{j=1}^n 3^(n-j) 2^(a_1+...+a_{j-1})` (odd), so
`S^n(x) ≡ C_w 2^(-A) = Y_n(w) mod 3^n`. (iii) Conversely, if `y ≡ Y_n(w)
mod 3^n` then `x := (2^A y - C_w)/3^n` is an odd integer (`2^A y` even,
`C_w` odd) lying in the class of (i) (its residue mod `2^(A+1)` is `(2^A -
C_w) 3^(-n)`, the same as for any element of that class, because `y` is
odd), hence `x` has word `w` and `S^n(x) = y`. So the words `w` with
`Y_n(w) ≡ y` are in bijection with `T_n(y)`, `w ↔ x`, and `2^(-A(w)) =
2^(-A_n(x))`. (iv) `mu_n(y mod 3^n) = sum_{w : Y_n(w) ≡ y} 2^(-A(w))`.
(v) Sign: `S` maps negatives to negatives (`3x+1 <= -2` for `x <= -1`), so
a predecessor of `y > 0` is positive; likewise for `y < 0`. ∎

**Corollary A1 (exact class means).** `Y_n mod 3 = 2^(-a_n) mod 3`, so
`mu_n(2 mod 3) = 2/3`: the mean of `rho_n` over the units is `1`, over the
class `1 mod 3` it is `2/3`, over the class `2 mod 3` it is `4/3`
(asserted in the deep script at every level).

**Corollary A2 (the layer recursion, PROVED).** Writing `H_n^{(r)}(y)` for
the part of `H_n(y)` carried by nodes `x ≡ r mod 3`:
`H_{n+1}(y) = H_n^{(1)}(y) + 2 H_n^{(2)}(y)`.
(A node `x ≡ 1` has the children `(2^a x - 1)/3`, `a` even, total weight
`sum 3 · 2^(-a) = 1`; a node `x ≡ 2` has `a` odd, total weight `2`; a
node `x ≡ 0` is a leaf.) For the seed `1`: the depth-1 layer is `{(4^k -
1)/3}` with weights `3 · 4^(-k)`, classes `1, 2, 0` for `k ≡ 1, 2, 0 mod
3`, so `H_1 = 1` splits as `16/21, 4/21, 1/21` and `H_2 = 16/21 + 8/21 =
8/7`, which is the exact `9 mu_2(1)`. Exact values: `H_1 = 1`, `H_2 =
8/7`, `H_3 = 1376/1387`; then `0.927627, 0.964301, 0.955420` (`n = 4, 5,
6`), all matched by direct enumeration of the tree with truncated
valuations. The first-arrival masses are much smaller at small depth
(`1/4, 0.393, 0.135, 0.184, 0.269, 0.232` for `n = 1..6`) because the loop
`1 -> 1` carries `(3/4)^n` and the all-arrival mass is the geometric
smoothing `H_n = sum_{j <= n} (3/4)^(n-j) H_j^first` (`H_0^first = 1`).
The first-arrival layers are far from equidistributed mod `3` (class-0
shares `0.19, 0.74, 0.26, 0.23, 0.29, 0.50`), and the recursion reproduces
each next first-arrival mass exactly.

**Reading.** Mazur's weight `omega(w) = 3^d 2^(-A(w))` on an inverse
history of the seed `M` is the summand of `H_d(M)`; his source charge `x
omega(w) <= M` is the identity `1/x = 3^d 2^(-A) · prod_{j<d}(1 +
1/(3 x_j))` (section 6) read as `x · 3^d 2^(-A) = M / prod(...) <= M`.
Lemma 2.1 is the statement that the harmonic mass is conserved in the mean
over residues; Theorem A says the *actual* mass of the tree of a seed is
the reference density at the seed's residue. The whole positive-density
argument is a lower bound on the harmonic mass of some tree, transported
to a count of small sources by the size identity.

---

## 5. Theorem B: cycle resonances and the exact spike profile

**Theorem B (PROVED, lower-bound form).** Let `y_0` lie on a Syracuse cycle
(on `Z`) of length `k` with word `w` and total valuation `A`. Then for
every `m >= 1`, `mu_{km}(y_0 mod 3^{km}) >= 2^(-Am)` and `rho_{km}(y_0) >=
(2/3)(3^k/2^A)^m`. Since `y_0 = C_w/(2^A - 3^k)` with `C_w > 0`, negative
cycles have `3^k > 2^A` (growing spikes), positive cycles `2^A > 3^k`
(decaying). *Proof.* `y_0 ∈ T_{km}(y_0)` with word `w^m`; apply Theorem A.
For the loop `-1` directly: `Y_n(1^n) = sum_{j<n} 3^j 2^(-j-1) = (3/2)^n -
1 = (3^n - 2^n)/2^n ≡ -1 mod 3^n`. ∎

**Consequences.** `sup rho_n >= (2/3)(3/2)^n`, so no level-uniform bound
on the reference density exists; in the 3-adic limit the density (when it
exists, section 7) has the singularities `|y - y_0|_3^(-s)`, `s = 1 -
(A/k) log_3 2`: `s = 0.3691` at `-1`, `0.0536` at `-5` and `-7` (`k = 2,
A = 3`), `0.0086` at the seven points `-17, -25, -37, -55, -41, -61, -91`
(`k = 7, A = 11`). A hypothetical positive cycle would have `s < 0`,
extremely close to `0` (its mean valuation exceeds `log_2 3` by the cycle
constraints), an anti-resonance nearly flat over the accessible levels: no
lever on the cycle problem here, but a clean statement of what a cycle
does to the 3-adic law.

**Numerics (FFT law, levels `<= 18`).** The maximal atom of `rho_n` on the
units is the class of `-1` at every level (`1.33, 2.10, 3.20, 4.85, 7.32,
11.0, 16.6, 24.9, ...`, `1441` at `n = 18`), with `rho_n(-1)/(3/2)^n ->
0.9748` (the loop word contributes exactly `(2/3)(3/2)^n`; the rest is the
negative tree of `-1`). The atoms of the `{-5,-7}` cycle grow (`rho_18(-5)
= 6.3`, `rho_18(-7) = 7.6`, against `1` for a typical unit); `-17` sits at
`2.0–2.5`.

**The spike profile (exact, OBSERVED to four decimals, PROVED as a lower
bound).** Define, for a 3-adic integer `x`, `D(x) := lim_n rho_n(x mod
3^n) / rho_n(-1 mod 3^n)`. The recursion `rho_{n+1}(z) = 3 sum_{a ≡ eps(z)
mod 2} 2^(-a) rho_n((2^a z - 1)/3 mod 3^n)` (the parents' form of the law;
asserted at `z = 1` at every level) and the growth `3/2` per level of the
spike give, on the forward rational closure of `-1` under all `S_a(y) =
(3y+1)/2^a`,
`D(z) = 2 sum_{a} 2^(-a) D((2^a z - 1)/3)`, `D(-1) = 1`, `D = 0` off the closure,
i.e. `D(z) = sum over backward paths from z to -1 of 2^(l - A)` (each step of
valuation `a` costs `2^(1-a)`; `a = 1` is free). Hence `D(-1/2) = 1/2`
(`a = 2` from `-1`), `D(-1/4) = 2(2^(-1) · 1/2 + 2^(-3) · 1) = 3/4`, and
the all-ones forward orbit of `-1/4`, `x_j = 3^(j+1)/2^(j+2) - 1 = -1/4,
1/8, 11/16, 49/32, ...`, all have `D = 3/4`; `D(-1/8) = D(1/16) = D(11/32)
= 3/8`, `D(-1/16) = 3/16`. Observed at levels 17 and 18: `0.5000, 0.7500,
0.7500, 0.7500, 0.7500, 0.3750, 0.3750, 0.3750, 0.1875`. The second tier of
atoms at level `n` is therefore the set of classes `3^(j+1)/2^(j+2) - 1 mod
3^n`, `j < n`, the "shadows" of `-1` (they agree with `-1` to level `j+1`),
each at `3/4` of the top atom, which is what the atom lists show at every
level (`n = 8`: `-1: 24.9`, then five classes at `18.7`, all of the form
`-1 + 2^(-(j+2)) 3^(j+1)`). Equality in the profile is a statement about
the limit and is left OBSERVED; the lower bound `rho_n(x) >= (2/3)(3/2)^n
· 2^(depth - A_0)` for a point at depth `depth` and valuation `A_0` in the
closure is PROVED by Theorem A applied to the word `1^(n-depth)` followed by
the path.

**Global statistics of the law (levels 1–18).**

| `n` | `rho_n(1)` | median of `rho_n` on units | `E[rho_n^2]` | entropy deficit `n ln 3 - H(mu_n)` | share of units with `rho < 0.1` | `||f_n - f_{n-1}||_1` |
|---|---|---|---|---|---|---|
| 6 | 0.6370 | 0.587 | 2.667 | 0.935 | 0.043 | 0.263 |
| 10 | 0.4250 | 0.533 | 3.911 | 1.062 | 0.058 | 0.178 |
| 14 | 0.3155 | 0.511 | 5.162 | 1.131 | 0.066 | 0.130 |
| 18 | 0.2736 | 0.499 | 6.419 | 1.173 | 0.070 | 0.097 |

`E[rho_n^2]` (over units; equals `(2/3) 3^n P(Y_n = Y_n')`) grows by `0.31`
per level: the law is not square-integrable in the limit, its density has
a tail `P(rho > t) ~ t^(-2)` fed by the spikes' forward closures. The
entropy deficit converges (increments `0.021, 0.018, ..., 0.0088` at `n =
11..18`, ratio `0.9`), so the information dimension is `1` and `E_units[rho
log rho] -> 0.84`. The `L^1` distance between consecutive Haar densities
decays like `n^(-0.93)`. The median density of a unit class settles at
`0.50`, with `7%` of the classes below `0.1`: the tree of a typical seed is
thin at many depths.

**Fine-scale distances.** `||mu_N - lift(mu_m)||_1` at `N = 18`: `0.854,
0.739, 0.650, 0.580, 0.523, 0.474, 0.430, 0.390, 0.354` for `m = 1..9`
(at `N = 17`: `0.851, ..., 0.345`; the dependence on `N` is weak). This is
the quantity bounded in (2.3) by `(2/3) C_A m^(-A)`, decreasing in `m` as
it should; the explicit `C` of Lemma 8.1 says nothing below `m` of the
order of `2^(2^6536)`.

---

## 6. Theorem C: the seed-1 test of positive-density log-time convergence

**Theorem C (PROVED).** Let `F_C = {x >= 1 : tau(x) <= C ln x}` and `H_n
:= 3^n mu_n(1 mod 3^n)`. If `F_C` has positive lower natural density for
some `C > 0`, then `limsup_n n^(1/6) H_n > 0`; indeed the Cesàro means of
`n^(1/6) H_n` are bounded below. In particular Mazur's Theorem 1.1 implies
`limsup_n n^(1/6) H_n > 0`, and `H_n = o(n^(-1/6))` would refute it (for
every constant `C`, not only `523/50`).

*Proof.* (a) Lower natural density `c` beyond `X_0` gives, by partial
summation, `sum_{x ∈ F_C, x < X} 1/x >= c ln(X/X_0)`. (b) Write `x = 2^k
m`, `m` odd. Then `tau(m) <= tau(x) <= C ln X`, and the Syracuse depth
`n(m)` (number of odd terms before `1`) is at most `tau(m)`. Hence
`sum_{x ∈ F_C, x<X} 1/x <= 2 + 2 sum_{1 <= n <= C ln X} G_n`, `G_n := sum_{m ∈ T_n^first(1)} 1/m`.
(c) For `m ∈ T_n^first(1)` with orbit `m = x_0, ..., x_n = 1` and word of
total valuation `A`, telescoping `2^{a_j} x_j = 3 x_{j-1} (1 + 1/(3
x_{j-1}))` gives the exact identity
`1/m = 3^n 2^(-A) Pi(m)`, `Pi(m) = prod_{j<n} (1 + 1/(3 x_j))`.
The `x_j` (`j < n`) are `n` distinct odd positive integers, so `sum_{j<n}
1/x_j <= sum_{k<=n} 1/(2k-1) <= 1 + (1/2) ln n`, and `Pi(m) <= e^(1/3)
n^(1/6)`. Therefore `G_n <= e^(1/3) n^(1/6) H_n^first <= e^(1/3) n^(1/6)
H_n` (first-arrival nodes are among all depth-`n` predecessors). (d)
Combining, with `N = C ln X`: `sum_{n <= N} n^(1/6) H_n >= (c /(2 e^(1/3)
C)) N - O(1)`, so `limsup n^(1/6) H_n >= c/(2 e^(1/3) C) > 0`. ∎

**The discount is mild in practice.** Over all odd `x <= 2·10^6`: `max
Pi(x) = 1.2531` (at `x = 993`), mean `1.167`, and `Pi(x)/(e^(1/3)
n^(1/6)) <= 0.764`; `Pi(9) = 1.2486` (`sum 1/x_j = 0.68`), `Pi(27) =
1.1989`. So `G_n` and `H_n^first` agree within a factor `1.25` on the
tested range, and the `n^(1/6)` is a proof convenience, not the truth.

**The data.**

| `n` | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 |
|---|---|---|---|---|---|---|---|---|---|
| `H_n(1) = 3^n mu_n(1)` | 1.143 | 0.928 | 0.955 | 0.775 | 0.637 | 0.537 | 0.473 | 0.435 | 0.410 |
| `n^(1/6) H_n(1)` | 1.28 | 1.17 | 1.29 | 1.10 | 0.94 | 0.81 | 0.73 | 0.69 | 0.66 |

Every ratio `H_{n+1}/H_n` from `n = 7` to `n = 17` is below `1` (`0.90` to
`0.975`; the last three `0.948, 0.967, 0.975`). By Corollary A2 the ratio
is `f_n^{(1)} + 2 f_n^{(2)} = 1 + f_n^{(2)} - f_n^{(0)}`: the layers of
the tree of `1` have for eleven consecutive levels more harmonic mass on
the leaves (`0 mod 3`) than on the doubly fertile class `2 mod 3`. The
"parents" form `rho_{n+1}(1) = 3 sum_j 4^(-j) rho_n(R_j)` with the other
seeds' values (`rho_18`: `R_2 = 5: 0.283`, `R_4 = 85: 0.227`, `R_5 = 341:
1.160`, `R_7: 0.047`, `R_8: 0.999`) gives the fixed-point extrapolation
`rho_∞(1) ≈ 0.237`, i.e. `H_∞(1) ≈ 0.35`, if those have converged (they
drift too). **Verdict: OPEN.** A positive limit near `0.35` and a slow
decay to `0` are both consistent with eighteen levels. The next level costs
`9×` (level 18 took two minutes and `20 GB`); level 19 or 20 is feasible on
this machine with a chunked FFT and would be the cheapest new evidence on
the theorem's conclusion available anywhere. This is the first point in the
thread where a Lean-checked positive-density claim meets an exact sequence
it must dominate.

**What Theorem C does not say.** Nothing about Mazur's proof; nothing about
the full tree of `1` without the log-time restriction (the companion
theorem's `a = 1` case has no depth cutoff, so no test of this kind
follows); nothing if `H_n` converges to a positive number.

---

## 7. Absolute continuity, typing, and what changes for the repo

**Absolute continuity (CONDITIONAL, derivation ours).** Let `f_n(y) = 3^n
mu_n(y mod 3^n)` be the Haar density of the uniform lift of `mu_n` to
`Z_3`; the laws are consistent (`Y_n mod 3^m` has the law of `Y_m`), so
`f_m = E[f_n | level m]` and `(f_n)` is a nonnegative martingale of mean
`1`. The estimate (2.3), `||f_q - f_m||_{L^1(Haar)} = ||mu_q - lift
mu_m||_1 <= C_A m^(-A)` uniformly in `q >= m`, is exactly the `L^1`-Cauchy
property; so under (2.3) `f_n -> f_∞` in `L^1`, and the 3-adic Syracuse law
is `f_∞ dHaar`, absolutely continuous, with `||f_∞ - f_m||_1 <= C_A
m^(-A)`. Its density is unbounded (Theorem B), infinite on the forward
rational closure of every negative cycle, not in `L^2` (the second moment
grows linearly), with finite entropy `int f log f`. We did not verify that
Tao's Proposition 1.14 is literally (2.3); the paper says so, and the
numerics of section 5 are consistent with it. If it is, absolute
continuity of the 3-adic Syracuse law is a corollary of Tao 2022 and
presumably known to its author; we have not seen it stated.

**Atlas typing (conclusion rule of atlas §0).** Conclusion: a positive
lower density of integers reaching `1` in logarithmic time.
* **DRIFT: O.** The conclusion is false for `5x+1` (its cycles' basins are
  conjecturally null); the proof uses the negative drift through the
  central restriction `A ≈ 2d` and the time constant `3/log(4/3)`.
* **STICKY: O.** The transfer weights `3^d 2^(-A)` are exact inverse-orbit
  sums at every scale (the paper: "these transfers agree with the actual
  inverse-orbit sums"; Theorem A), i.e. size-free memorylessness, which the
  aliquot control lacks. The pattern of atlas §7 holds: a residue-averaging
  mechanism.
* **SHEET: B.** The conclusion holds for `3x+1` on the negatives (each
  negative cycle conjecturally attracts a positive density in log time);
  the companion theorem makes this explicit: every target `a` with `3 ∤ a`
  has a positive-density basin, so an extra cycle would have one too.
* **DEFECT: B.** A null planted set does not change a positive-density
  conclusion.
* **INTEGRAL: B** (the transfer operators act on `Z/3^t` and the Terras
  classes; a rational cycle in `Z_(2)` is invisible).
* **UNIFORM: B** (presumably runs unchanged for `3x+k`, `3 ∤ k`, with
  seeds `(4^j - k)/3`; not checked).
* **DIM: B.** Positive density leaves the whole non-descending set
  unresolved; no descent certificate is produced.

**For the repo.** (i) The Krasikov–Lagarias `x^0.84` row (predecessors of
`1`) would be replaced by positive density, and the "entry lane" coverage
of S15/S16 would get a certificate that is not class-based (the seed is a
single integer, the histories are actual orbits); the divergence half
(THM-4476/4499, DRIFT) is untouched, as the typing says. (ii) The
companion's two elementary corollaries are worth keeping in view: a single
non-convergent `a` implies a positive-density set of non-convergent
integers (its predecessors), and convergence on a set of upper density `1`
implies convergence everywhere. With THM-4476 (every divergent orbit is
thin), a counterexample would have to come with a positive-density *union*
of thin divergent orbits or cycle basins; no contradiction, a shape. (iii)
The 3-adic Syracuse density is a new object for the thread: its
resonances are the negative cycles, its typical value is `1/2`, its dead
zones (`7%` of classes below `0.1`) are the residues whose trees are thin,
and the S18 observation that `80%` of visited values above the start range
are `2 mod 3` is the class mean `4/3` of Corollary A1 seen from the orbit
side. (iv) Theorem C is the leverage the owner asked for: a Lean-checked
claim converted into a computable inequality on an exact sequence.

---

## 8. The inverter (the owner's pasted text) and the Periodicity Conjecture

The map from a valuation sequence `a = (a_1, a_2, ...)` to the unique
2-adic integer with that Syracuse word is Bernstein's conjugacy
`Phi(a) = - sum_{j>=1} 2^(a_1+...+a_{j-1}) / 3^j ∈ Z_2`
(the 2-adic limit of `x = (2^A x_n - C)/3^n`); THM-4476 writes it as
`R_2(d)` in the cumulative variables `d_l`. Facts: (1) periodic (or
eventually periodic) `a` gives a rational `Phi(a)` with odd denominator
(geometric series); (2) conversely a rational `p/q` (`q` odd) has the
word of the integer `p` under `x -> (3x+q)/2^v`; so "an unrepeating sequence
with rational value" is the same thing as "a divergent integer orbit of some
`3x+k` map, `k` odd", and a natural-number value is a divergent `3x+1`
orbit: the pasted claim is right, and the pasted open question is the
**Periodicity Conjecture** (Lagarias 1985; Bernstein–Lagarias 1996; from
memory): `Phi^(-1)(Q ∩ Z_2)` is exactly the set of eventually periodic
words, equivalently every `3x+k` orbit on `Z` is eventually periodic. (3)
The repo's partial result is THM-4476, corollaries (2), (3), (8): if `Phi(a)
= n` is a positive integer for an aperiodic `a`, the real series `R(a)`
converges strictly below `n` and the orbit diverges at full rate `c = n
prod(1 - 1/(3 m_l)) > 0`; no such orbit has `m_j <= K j^alpha` for all `j`
with `alpha < 1/h* = 1.05268`, hence **every aperiodic word whose
minus-sheet discrepancy stays within `C log_2 j` with `C < 1.05268` (or
plus-sheet with `C < 0.02634`) has a Bernstein value that is not an
odd-denominator rational**; and an orbit is eventually periodic iff its
reciprocal sum diverges, so the Periodicity Conjecture is equivalent to
`sum_i 1/|T_k^i(n)| = ∞` for every odd `k` and every `n`. In the language of
this note: a rational value of an aperiodic word is a point of `Z_(2)`
whose *2-adic* word never repeats while its *3-adic* Syracuse law (section
5) is that of an integer orbit; Theorem B shows what a *periodic* word does
to the 3-adic side (a resonance or anti-resonance at the cycle point); an
aperiodic rational would be a point of `Z_3` carrying, level by level, the
mass `3^n 2^(-A_n)` of a single divergent orbit, which THM-4476 (5) makes
`o(1)`: its discrepancy tends to `-∞`. So the two sides agree: an aperiodic
rational value is a vanishing atom, never a resonance.

---

## 9. Kaprekar's routine and the perfect cuboid (the owner's numbers)

Kaprekar's map on `d`-digit strings (descending digits minus ascending
digits, leading zeros kept): `d = 2`: the single cycle `9 -> 81 -> 63 ->
27 -> 45 -> 9` (`= 9 · (1, 9, 7, 3, 5)`), reached by all 90 non-repdigit
starts; `d = 3`: the fixed point `495`; `d = 4`: `6174`; `d = 5`: two
4-cycles and the 2-cycle `{53955, 59994}`; `d = 6`: a 7-cycle and the fixed
points `631764`, `549945`. Every image is divisible by `9`, and the map
factors through the multiset of digits: a finite state space per `d`. In
the nine-iteration typology of the STICKY note ("provability tracks the
presence of an exact invariant"), Kaprekar is the extreme case: the
invariant is a finite quotient, so every statement about it is a theorem by
enumeration, and nothing transfers to Collatz beyond the typology's row.
`6174 = 2 · 3^2 · 7^3`, `495 = 3^2 · 5 · 11`; no structural link to the
thread was found and none is claimed. Euler bricks with edges below 300:
`(44, 117, 240)`, `(240, 252, 275)`; neither has an integral space
diagonal; the perfect cuboid is open (searches far beyond this range are in
the literature; from memory). Its relevance to the thread is the same as
the Ogg/triangular seed of S17: a source of three-way integer constraints,
not a Collatz technique.

---

## 10. Reproduction and status table

```
python 04-computation/experiments/mazur_positive_density_20260928.py   > 05-knowledge/results/mazur_positive_density_20260928.out
python 04-computation/experiments/mazur_harmonic_mass_20260928.py      > 05-knowledge/results/mazur_harmonic_mass_20260928.out
python 04-computation/experiments/mazur_harmonic_mass_deep_20260928.py 17 > 05-knowledge/results/mazur_harmonic_mass_deep_20260928.out
python 04-computation/experiments/mazur_harmonic_mass_deep_20260928.py 18 > 05-knowledge/results/mazur_harmonic_mass_deep18_20260928.out   # 20 GB, 2 min
python 04-computation/experiments/mazur_seed1_test_20260928.py         > 05-knowledge/results/mazur_seed1_test_20260928.out
python 04-computation/experiments/kaprekar_cuboid_20260928.py          > 05-knowledge/results/kaprekar_cuboid_20260928.out
```

| item | status |
|---|---|
| Mazur Theorem 1.1, Corollary 1.2; companion theorem and corollaries | CITED (platform Lean-checked; review pending / accepted formalization; no independent replay) |
| P1–P8 elementary components | FINITE-EXACT, all consistent with the paper |
| Theorem A (harmonic mass = `3^n mu_n(seed)`), Corollaries A1, A2 | PROVED; FINITE-EXACT to level 6 (rational) and 18 (float) |
| Theorem B (cycle resonances; `-1` exact; singularity exponents) | PROVED (lower bounds); OBSERVED growth constant `0.9748` |
| spike profile `D` on the closure of `-1` | PROVED as lower bound; equality OBSERVED to four decimals (levels 17, 18) |
| Theorem C (positive-density log-time convergence forces `limsup n^(1/6) H_n(1) > 0`) | PROVED |
| `lim H_n(1) > 0`? | OPEN; `H_18 = 0.410`, ratios `< 1` for eleven levels; extrapolation `≈ 0.35` |
| absolute continuity of the 3-adic Syracuse law | CONDITIONAL on (2.3) as stated (Tao Prop. 1.14 per the paper); numerics consistent |
| second moment linear, entropy deficit convergent, median `0.50` | OBSERVED (levels `<= 18`) |
| inverter = Bernstein map; open question = Periodicity Conjecture; THM-4476 partial result | CITED + PROVED (repo) |
| Kaprekar cycles, `9 | image`, digit-multiset invariant; Euler bricks | FINITE-EXACT; perfect cuboid OPEN (CITED, from memory) |

**Next probes.** Level 19–20 of the law (chunked FFT), for `H_n(1)`; the
class split `f_n^{(0,1,2)}` of the tree of `1` at depth `n` from the joint
law of `(Y_n mod 3^(n+1), A mod 2)`; the exact profile theorem (`D` as the
Green's function of the transfer operator at the `-1` eigen-singularity);
a proof or refutation of `liminf H_n(1) > 0`, which, by Theorem C, is the
cheapest possible confrontation with a positive-density claim.
