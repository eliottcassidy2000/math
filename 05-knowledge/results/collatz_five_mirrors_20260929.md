# Five papers as structural mirrors of the 3-adic Syracuse law: the Fourier profile of Tao's Syracuse random variable (primitive maxima at the powers of two, decay ratio rising to `0.87` per level, scale-invariant Fourier mass `0.47`), the same-length collision threshold, and carry-polynomial reciprocity under word reversal

**Session:** opus, `collatz-poset-dag-20260927` (S20), 2026-09-29.
**Owner's directive:** "merge ideas related to the linked papers below into
the current mathematical exploration process and extend the ideas of past
work in the repo toward Collatz related proofs", with five arXiv papers of
2026-09-26/27: Viaclovsky, *A two-parameter family of complex structures on
`S^6`* (2609.33785); Narode, *Necessary and sufficient condition for
reciprocal polynomials to be monogenic* (2609.33766); Merca, *The two-block
odd partition function* (2609.33770); Dalfó–Fiol–Reyes, *A note on
three-quarters circulant digraphs* (2609.33718); Lyu, *Polynomial
compressibility and forbidden oriented forests* (2609.33600). Per the
standing pattern for fresh papers (memory: structural analogy, the shape of
the hard direction, not subject overlap), each paper is read for the shape
of its argument and matched to the thread's frontier (S19: Theorems A–C of
[`mazur_positive_density_20260928.md`](mazur_positive_density_20260928.md);
THM-4476, THM-4514; the barrier atlas §§0, 7, 8).
**Inherits (cited):** Tao 2022 (Syracuse random variables; fine-scale
mixing, Proposition 1.14; from memory, as restated by Mazur's (2.3)), Mazur
2026 (Lemma 8.1: "primitive Fourier coefficient bound `C* h^(-6409)` at level
`h`"), Terras 1976, S19 Theorem A (harmonic mass), the S19 Fourier
consistency (`Y_n mod 3^m` has the law of `Y_m`).

**Status: PROVED (Propositions 1–5: level-independence of the primitive
Fourier profile; the implication (2.3) ⟹ super-polynomial decay of the
primitive maxima; no uniform one-step gap of the geometric Gauss sums; the
same-length spread lemma; carry reciprocity under reversal) + FINITE-EXACT
(the Fourier profile of the law to level 17 in float64; the collision table;
the reversed cycles; the Gauss-sum table) + OBSERVED (the maxima sit at the
powers of two `2^(h+3)`, `2^(h+4)`; the decay ratio rises from `0.65` to
`0.87` and levels; the Fourier mass per conductor level is `0.47`) +
DIRECTION. No Collatz proof step. Audit: pending (section 8 will hold the
record).**

Scripts and outputs:
`04-computation/experiments/collatz_five_mirrors_20260929.py` →
`05-knowledge/results/collatz_five_mirrors_20260929.out` (profile to level
14, Gauss sums, collisions, reciprocity);
`collatz_five_mirrors_fourier_deep_20260929.py` → `..._fourier_deep_...out`
(profile to level 17; `8 GB`, one minute);
`collatz_five_mirrors_costsplit_20260929.py` → `..._costsplit_...out` and
`..._costsplit14_...out` (the resonant coefficient decomposed by total
cost and by last valuation at levels 10 and 14).

---

## 0. What is new, in one screen

1. **The Fourier profile of the 3-adic Syracuse law is a single sequence
   `M(h)`, and Mazur's mixing estimate is a statement about it.** With
   `mu_hat_n(t) = E e(t Y_n/3^n)`, consistency gives `mu_hat_n(3^(n-h) u) =
   mu_hat_h(u)` for every unit `u mod 3^h` (Proposition 1, PROVED), so the
   maximal coefficient of conductor exactly `3^h`, `M(h) := max_u
   |mu_hat_h(u)|`, does not depend on the level it is read at. The
   fine-scale estimate (2.3) of the Mazur digest (Tao's Proposition 1.14 in
   Mazur's restatement) implies `M(h) <= C_A (h-1)^(-A)` for every `A`
   (Proposition 2, PROVED): the primitive Fourier maxima must decay faster
   than any polynomial. Computed to `h = 17` (FINITE-EXACT, float64):
   `0.5774, 0.3779, 0.2522, 0.1770, 0.1293, 0.0961, 0.0759, 0.0609, 0.0480,
   0.0383, 0.0319, 0.0265, 0.0221, 0.0191, 0.0163, 0.0144, 0.0125`; the
   maximum sits at `u = 2^s` with `s = h + 3` or `h + 4` (`2^17` at `h =
   14`, `2^21` at `h = 17`); the ratio `M(h)/M(h-1)` rises from `0.65` to
   about `0.87` and levels off over `h = 14..17` (`0.867, 0.851, 0.885,
   0.868`). A fixed power law is excluded by the steepening (local exponents
   `1.76, 1.86, 1.99, 2.08` on the doublings `5→10, 6→12, 7→14, 8→16`); a
   geometric decay at rate about `0.87` per level fits levels `11–17` and
   would satisfy (2.3) with room. No tension with Tao's proposition, whose
   explicit constants are vacuous here; but this is the first measurement of
   what the mixing estimate actually bounds.
2. **Scale invariance of the Fourier energy.** The Fourier mass at conductor
   level `h`, `sum_(cond t = 3^h) |mu_hat(t)|^2`, is `0.667, 0.476, 0.462`
   and then `0.464 .. 0.472` for `h = 4..17` (slowly increasing): every
   3-adic scale carries the same energy, which is exactly the linear growth
   of the second moment of the density found in S19 (`(3/2)·0.31 = 0.466`).
   The typical coefficient at level `h` has `|mu_hat|^2 3^h = 0.70`:
   square-root cancellation, as for random phases; the slow sup-norm decay
   comes from the powers of two alone.
3. **Where the resonance comes from.** At `h = 10` the resonant coefficient
   (`t = 2^12`, `0.0383`, against `0.0026` for a generic unit) is carried by
   the words of total cost `A = 12..17` (`60%` of `|.|` coherent), far below
   the typical cost `2h = 20`, and by the last valuations `1, 2, 3`
   (`0.0225, 0.0120, 0.0051`): the low-cost words, exponentially rare but
   phase-coherent, are the obstruction to fast sup-norm mixing. At `h = 14`
   the same decomposition gives the same picture (`A = 19..23`, coherent
   fraction `0.56`, last valuation `1` carrying `57%`), and the
   contribution-weighted cost ratio is stable, `A/h = 1.49` at `h = 10`,
   `1.48` at `h = 14`: the resonance is a fixed large-deviation family,
   whose mass rate at `A/h = 1.48` is `e^(-0.093 h)` (ratio `0.91` per
   level); the observed `0.87` is that rate times the coherence loss.
4. **One-step geometric Gauss sums have no uniform gap (Proposition 3,
   PROVED).** `G_j(t) = c_j sum_r 2^(-r) e(t 2^(-r)/3^j)` satisfies `sup_t
   |G_j(t)| -> 1` (`0.577, 0.582, 0.789, 0.887, 0.944, 0.971, 0.986, 0.994,
   0.997, 0.9985, 0.9993, 0.9997` for `j = 1..12`, at `t = 2^(j+1)`), while
   the mean over units settles at `0.5430`; the 2-adic reading `e(t m_r/2^r
   + t/(2^r 3^j))`, `m_r = -3^(-j) mod 2^r`, is exact. Mixing is a
   multi-step (renewal) phenomenon, as in Tao's proof, never a spectral gap
   of one step.
5. **The same-length spread lemma (Proposition 4, PROVED)** — the Moore
   bound of the circulant paper: two distinct words of the same length `d`
   with costs `A, A'` and `A + A' <= (n - d) log_2 3 + 1` land on distinct
   classes mod `3^n`. The actual first collisions in the tree of `1` need
   cost sums `29..63` at depths `2..4` for `n = 4..12`, three to five times
   the bound: the "tessellation obstruction" of the paper (optimal tiles that
   fail to tile) has no counterpart; the Syracuse residues stay injective far
   beyond the Moore regime.
6. **Reciprocity (Proposition 5, PROVED):** with the exclusive-prefix carry
   `C_w(u,v) = sum_j u^(d-j) v^(a_1+...+a_(j-1))` (the affine identity `3^d x
   + C_w = 2^A y`) and the inclusive one `C'_w`, word reversal is polynomial
   reciprocity, `C_(rev w)(u,v) = u^(d-1) v^A C'_w(1/u, 1/v)`; the 2-adic
   source `x ≡ -C_w 3^(-d)` and the 3-adic target `y ≡ C_w 2^(-A)` are the
   two reciprocal evaluations of one integer. Reversal is an involution on
   rational cycles that preserves `(k, A)`; it fixes the cycles of `-1` and
   `{-5,-7}` (a rotation) and sends the seven-cycle of `-17` to the rational
   point `-13801/139`, not an integer (FINITE-EXACT).
7. **The other two mirrors are shape only** (section 1): Viaclovsky's
   contracting affine lift and monodromy relation are the affine IFS
   `S_a(y) = (3y+1)/2^a` and the cycle equation, with the difference that a
   Collatz cycle is a fixed point of a non-identity affine map; Lyu's
   "bounded local cliques do not bound compressibility; only forests do" is
   the atlas's lesson that local word constraints do not tame the global
   stopping time.

---

## 1. The five papers and their mirrors

**Viaclovsky (2609.33785).** Complex structures on `S^6` from the rational
elliptic surface `y^2 = x^3 + t x + 1` (fibers `III*, I_1, I_1, I_1`;
Mordell–Weil generated by `P = (0, 1)`): the line bundle `M = O(P - O)`
minus its zero section, quotiented by a fiberwise linear lift `Φ = λ Φ_0`
of translation by `-2P`, contracting for small `λ`, gives a torus
fibration; Mumford fillings at the nodal values and a logarithmic
transform at infinity (an order-four free affine action) produce a multiple
fiber; the integral monodromies satisfy `T_* T_1 T_2 T_3 = I`, `T_*^4 = I`,
`(T_i - I)^2 = 0`; a gauge-invariant integer `d` in the winding data of the
elliptic logarithm is forced to `d^2 = 1` by primitivity of the local
monodromy image. *Shape of the hard direction:* build a global object as a
quotient by a contracting affine lift of a translation, then control the
degenerate fibers integrally. *Mirror:* the Syracuse maps `S_a(y) = (3y +
1)/2^a` are contracting affine maps of `Z_3` (ratio `|3|_3 = 1/3`), the
3-adic Syracuse law is the stationary law of this IFS (S19 §4), and the
inverse histories are the "quotient by the lift". The cycle equation
`3^k y_0 + C_w = 2^A y_0` is the sphere relation for the affine monodromy
around a cycle, with linear part `3^k/2^A ≠ 1`: a Collatz cycle is a fixed
point of a non-identity affine map, whereas in the torus fibration the
linear parts multiply to the identity and the invariant lives in the
translation part (`d`). The primitivity used to pin `d` is the Terras
primitivity of the word class (one class mod `2^(A+1)`), the step (i) of
S19's Theorem A. Verdict: shape only; nothing to compute.

**Narode (2609.33766).** A reciprocal polynomial `f(x) = x^n g(x + 1/x)`
is monogenic iff `g` is monogenic and `f(1)`, `f(-1)` are individually
squarefree (the earlier condition "`f(1) f(-1)` squarefree" is not
necessary: `x^4 + 3x^2 + 1`). *Shape:* halve the degree by the symmetry `x
↔ 1/x`, and locate the whole obstruction at the two fixed points `x = ±1`.
*Mirror:* word reversal is polynomial reciprocity of the carry
(Proposition 5); the fixed points of reversal are the palindromic words,
and the cycles fixed by reversal up to rotation are `-1` and `{-5,-7}` but
not the seven-cycle; the two evaluations `C_w(3,2)/3^d` (2-adic source
class) and `C_w(3,2)/2^A` (3-adic target class) are the "`x` and `1/x`" of
one integer. No degree-halving analogue was found: the carry has no
symmetric variable. Verdict: one exact identity, one finite fact.

**Merca (2609.33770).** `a(n)`, the signed count of partitions into
exactly two part sizes each with odd multiplicity (a double Lambert series
`sum q^(n_1)/(1+q^(2n_1)) · q^(n_2)/(1+q^(2n_2))`), equals `(σ(n) -
λ(n))/4 - σ(n/2)/2 + σ(n/4)` (Jacobi's two- and four-square counts), hence
`a(n) >= 0`; the generating function factors as the `pod` product times a
theta series over triangular numbers with signs `(-1)^(T_(r-1))`; `a(12n +
11) ≡ 0 mod 3`. *Shape:* a signed sum whose positivity is invisible in its
definition becomes obvious in a cancellation-free arithmetic description.
*Mirror:* the 3-adic law `mu_n(y) = sum_w 2^(-A(w)) [Y_n(w) ≡ y]` is an
`n`-fold geometric Lambert sum (weights `2^(-a)`, 3-adic phases); its
Fourier expansion `mu_n(y) = 3^(-n) sum_t e(-ty/3^n) mu_hat_n(t)` is the
signed description, and S19's Theorem A (`3^n mu_n(y)` = the harmonic mass
of the predecessors of `y`) is the cancellation-free one; the Lambert
factors are the one-step Gauss sums of section 4, evaluated through the
other prime by `e(t m_r/2^r + t/(2^r 3^j))`. The thread's positivity
question, `liminf H_n(1) > 0`, is the seed-`1` value of that description;
its layer recursion `H_(n+1) = H^(1) + 2H^(2)` (S19, Corollary A2) is the
analogue of Merca's divisor recursion. Verdict: the frame is right; a
closed formula for `H_n(1)` is not in reach (it would be a formula for the
tree of `1`).

**Dalfó–Fiol–Reyes (2609.33718).** On the circulant digraph `CD(N, a, b)`
allow the step pairs `(+a,+b), (+a,-b), (-a,+b)` but not `(-a,-b)`: a
sector-constrained distance, planar layers of size `3ℓ + 1`, Moore-type
bound `N(k) = (3k^2 + 5k + 2)/2` on the order at diameter `k`, not
attained for `k > 1` because the optimal tiles do not tessellate the plane;
admissible generator pairs from the lattice of translations (`u_i a + v_i
b ≡ 0 mod N`, `N = det`). *Shape:* a restricted step set on a cyclic group,
a Moore bound, and an obstruction (tiling) separating the counting optimum
from realizability. *Mirror:* the frequency-side recursion of the 3-adic
law, `mu_hat_n(t) = sum_(a>=1) 2^(-a) e(t 2^(-a)/3^n) mu_hat_(n-1)(t 2^(-a)
mod 3^(n-1))` (verified exactly), is a twisted geometric walk on the cyclic
unit group `<2> ≅ Z/(2·3^(n-1))` — a weighted circulant with generator
`2^(-1)` — with the sector restriction that the parity of `a` is fixed by
the class of the target mod `3`; the ball growth is exponential
(`2^(A-1)` words of cost `A`), not planar; the Moore bound is Proposition 4
and it is far from attained (section 5): there is no tiling obstruction,
the residues of same-length words stay injective long past the bound.
Verdict: the sharpest of the five; it produced the Fourier profile of
section 2.

**Lyu (2609.33600).** For acyclic oriented graphs, the compressibility
`τ(H)` (least `n` such that `H` maps into every tournament of order `n`)
satisfies `p(H) <= τ(H) <= r_tr(p(H))` (longest path, transitive Ramsey
number), and the upper bound is attained by connected graphs of arbitrary
girth with oriented clique number `3`: bounding local cliques does not
give polynomial control; a forbidden graph gives a polynomially
`τ`-bounded class only if its underlying graph is a forest; the alternating
orientation of `P_4` gives `τ = p` when triangle-free. *Shape:* local
constraints do not tame a global Ramsey-type parameter; only tree-like
forbidden structure does. *Mirror:* the barrier atlas: bounded valuations,
forbidden runs and residue constraints on the valuation word never bound
the stopping time (words with all `a_j <= 2` and mean below `log_2 3` exist;
the atlas's UNIFORM/DIM rows), and the one constraint that does organise
the orbit is tree-like, the descent tree `D(m)` of THM-4514 (every value
has a unique parent below it). `τ >= p` is "time to `1` at least the
depth". The repo's tournament thread (Paley `T_p`, flip-rank; THM-640) is
a side link only. Verdict: remark.

---

## 2. The Fourier profile of the 3-adic Syracuse law

**Definitions.** `mu_n` is the law of `Y_n` on `Z/3^n` (`Y_0 = 0`, `Y_(k+1)
= 2^(-A)(3 Y_k + 1)`, `P(A = a) = 2^(-a)`); `mu_hat_n(t) = sum_y mu_n(y)
e(ty/3^n)` (any sign convention; only moduli are used). For `t ≠ 0` write
`t = 3^(n-h) u` with `3 ∤ u`; `h` is the conductor level.

**Proposition 1 (level independence, PROVED).** `mu_hat_n(3^(n-h) u) =
mu_hat_h(u)` for every unit `u mod 3^h`, `1 <= h <= n`. *Proof.* `e(3^(n-h)
u Y_n/3^n) = e(u (Y_n mod 3^h)/3^h)` and `Y_n mod 3^h` has the law `mu_h`
(consistency: `Y_n mod 3^h` is the function `F_h` of the last `h`
valuations, S19 §7 and the audit's item 13). ∎ Hence `M(h) := max_(u unit
mod 3^h) |mu_hat_h(u)|` is intrinsic; the script confirms the same profile
at `n = 4, 6, ..., 14, 17`.

**Proposition 2 ((2.3) forces super-polynomial decay, PROVED).** If
`||mu_q - lift(mu_m)||_1 <= C_A m^(-A)` for `1 <= m <= q` (Mazur's (2.3),
the finite-distribution form of Tao's Proposition 1.14), then `M(h) <= C_A
(h-1)^(-A)` for all `h >= 2`. *Proof.* Take `q = h`, `m = h-1`. The lift
of `mu_(h-1)` is constant on the cosets of `3^(h-1) Z/3^h Z`, so its
Fourier coefficient at a primitive `u` vanishes (`sum_(k=0)^2 e(uk/3) =
0`), and `|mu_hat_h(u)| = |sum_y (mu_h - lift)(y) e(uy/3^h)| <= ||mu_h -
lift mu_(h-1)||_1`. ∎ Mazur's Lemma 8.1 is the case `A = 6` with his `C`;
his Section 8.1 describes the underlying "primitive Fourier coefficient
bound `C* h^(-6409)` at level `h`", i.e. a bound on `M(h)` with `C*` a
three-fold exponential tower.

**The profile (FINITE-EXACT, float64 FFT of the level-17 law; identical at
every level where compared).**

| `h` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `M(h)` | .5774 | .3779 | .2522 | .1770 | .1293 | .0961 | .0759 | .0609 | .0480 | .0383 | .0319 | .0265 | .0221 | .0191 | .0163 | .0144 | .0125 |
| ratio | – | .655 | .667 | .702 | .730 | .743 | .789 | .803 | .789 | .797 | .835 | .828 | .834 | .867 | .851 | .885 | .868 |
| argmax `u` | 1 | 2^2 | 2^3 | 2^4 | 2^5 | 2^6 | 2^8 | 2^9 | 2^10 | 2^12 | 2^13 | 2^14 | 2^16 | 2^17 | 2^18 | 2^20 | 2^21 |
| Fourier mass at level `h` | .667 | .476 | .462 | .464 | .466 | .466 | .465 | .466 | .467 | .467 | .468 | .469 | .470 | .470 | .471 | .471 | .472 |
| typical `|mu_hat|^2 3^h` | 1.00 | .714 | .692 | .696 | .698 | .699 | .698 | .699 | .700 | .701 | .702 | .703 | .704 | .705 | .706 | .707 | .708 |

`M(1) = 1/sqrt 3` exactly (`|(1/3) e(1/3) + (2/3) e(2/3)|`). The maximum
is always at a power of two, `u = 2^s` with `s = h + 3` or `h + 4` for `h
>= 7`, i.e. `t/3^h = 2^s/3^h ≈ 0.03`: the character `e(2^s y/3^h)` is not
a low frequency in the real sense; it is the frequency that turns the
3-adic inverses `2^(-a)` of the small valuations into small integers
`2^(s-a)` (the one-step Gauss sums of section 4 show the same alignment).
The values `|mu_hat_17(2^s)|`, `s = 0..60`, form a single bump centred at
`s = 21` (`0.0125`), with `0.0111` at `s = 20`, `0.0106` at `23`, and the
generic size `3·10^(-4)` (`= 0.84 · 3^(-17/2)`) outside `s = 13..27`.

**Decay.** Local exponents from the doublings `M(h)/M(2h)`: `1.76` (5→10),
`1.86` (6→12), `1.99` (7→14), `2.08` (8→16): steepening, so no fixed
power law; `3.75 h^(-2)` fits `7 <= h <= 14` to `4%` and undershoots by
`2–4%` at `15–17`. The ratios rise from `0.65` to `0.87` and are flat over
`14..17` at `0.87 ± 0.02` (they fluctuate with the `s = h+3`/`h+4`
alternation of the argmax). Reading: a geometric decay at rate about
`0.87` per level, which is compatible with Proposition 2 (super-polynomial)
with room, and with Tao's proposition, whose constants say nothing below
`m` of the order of `2^(2^(2^8697))`. OBSERVED, seventeen levels; the
asymptotic regime is not established.

**Scale invariance.** By Parseval, `sum_t |mu_hat_n(t)|^2 = 3^n sum_y
mu_n(y)^2 = (3/2) E_units[rho_n^2]`; the second moment grows by `0.31` per
level (S19), so the Fourier mass per conductor level is `(3/2)·0.31 =
0.466`, the table's constant. Every 3-adic scale carries the same Fourier
energy: the law is "`1/f`" in the 3-adic scale, its density is in no `L^p`
for `p >= 2` if this continues (S19 §7, OBSERVED), and the fine-scale
`ell^1` distances `||mu_N - lift mu_m||_1` (`0.85 .. 0.35` for `m = 1..9`
at `N = 18`) are dominated by the bulk of square-root-size coefficients,
not by the maxima: `d(m, N) >= M(m+1)` (Proposition 2's inequality) is far
from tight (`0.74` against `0.25` at `m = 2`).

**The resonance, decomposed (`h = 10`, `t = 2^12`).** The joint law of
`(Y_10 mod 3^10, A)` (cost coordinate to `50`) splits the coefficient
`0.0383` into contributions by total cost: `A = 14: 0.0168` (mass
`0.044`), `A = 16: 0.0144` (mass `0.076`), `A = 13: 0.0095`, `A = 15:
0.0072`, `A = 17: 0.0044`, then `A = 12, 11, 19, 18, 10` below `0.004`;
`sum_A |.| = 0.0658`, coherent fraction `0.58`. By last valuation: `a = 1:
0.0225`, `a = 2: 0.0120`, `a = 3: 0.0051`, `a = 4: 0.0019`. At `h = 14`
(`t = 2^17`, `0.0191`; `collatz_five_mirrors_costsplit14_20260929.out`):
`A = 21: 0.0083`, `A = 19: 0.0062`, `A = 20: 0.0044`, `A = 22: 0.0043`, `A
= 23: 0.0034`, coherent fraction `0.56`, `a = 1: 0.0110`, `a = 2: 0.0054`,
`a = 3: 0.0023`. The contribution-weighted mean of `A/h` is `1.49` at `h
= 10` and `1.48` at `h = 14`, against the mass-weighted mean `2.00`. So the
resonance is carried by a fixed family, the words of total cost `A ≈ 1.48
h` with small last valuations: exponentially rare (the rate function of the
cost, `I(α) = α H(1/α) - α ln 2` in nats for `A = α h`, `H` the entropy,
gives `I(1.48) = -0.093`, i.e. mass `e^(-0.093 h)`, ratio `0.91` per
level) but phase-coherent; the observed sup-norm ratio `0.87` is this rate
times a coherence loss of about `0.96` per level. This is the mechanism behind the slow
sup-norm decay and the reason it must be renewal-theoretic (Tao's
Fourier–renewal method, Mazur's §8.1 "renewal bounds"): no single step has
a gap (section 4), and the obstruction is a large-deviation family of
words, whose rate (`≈ 0.87` per level) is what a sharp version of (2.3)
would have to compute.

**Direction.** A proof that `|mu_hat_h(2^s)| <= C r^h` for some `r < 1`
uniformly in `s`, together with the square-root behaviour of the generic
coefficients, would give (2.3) with a usable constant; the cost
decomposition suggests the route: bound the coherent low-cost family by its
large-deviation mass and show the remaining words cancel at the
square-root rate. This is the analytic core of Tao's proposition restated
as a concrete inequality on one explicit sequence.

---

## 3. The same-length spread lemma (the Moore mirror)

**Proposition 4 (PROVED).** Let `w ≠ w'` be words of the same length `d`
with total costs `A, A'`. If `Y(w) ≡ Y(w') mod 3^n` then `A + A' > (n - d)
log_2 3 + 1`. Equivalently, same-length words with `A + A' <= (n - d) log_2
3 + 1` land on distinct classes mod `3^n`, and (by S19's Theorem A) the
corresponding depth-`d` predecessors of any two integers in one class mod
`3^n` are distinct nodes.

*Proof.* `Y(w) = C_w/2^A` with `C_w = sum_(j=1)^d 3^(d-j) 2^(a_1+...+a_(j-1))`
odd. The congruence gives `3^n | C_w 2^(A') - C_(w') 2^A`. The difference is
nonzero: equality would force `A = A'` (both carries odd) and `C_w =
C_(w')`, and `(d, C_w)` determines `w` (`a_1 = v_2(C_w - 3^(d-1))`, then
recurse on `(C_w - 3^(d-1))/2^(a_1)`). And `0 < |C_w 2^(A') - C_(w') 2^A|
< 2^(A+A') 3^d / 2` since `C_w <= 2^A (3^d - 1)/2`. Hence `3^n < 2^(A+A')
3^d/2`. ∎

**The actual first collisions in the tree of `1`** (nodes of depth `d`
with valuations `<= 24`, reduced mod `3^n`; minimal `A + A'` among the
colliding pairs at the first colliding depth):

| `n` | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| first depth `d` | 2 | 2 | 3 | 3 | 3 | 3 | 4 | 4 | 4 |
| minimal `A + A'` | 29 | 31 | 32 | 32 | 56 | 56 | 58 | 61 | 63 |
| bound `(n-d) log_2 3 + 1` | 4.2 | 5.8 | 5.8 | 7.3 | 8.9 | 10.5 | 10.5 | 12.1 | 13.7 |

The bound is loose by a factor `3–5`: unlike the three-quarters digraphs,
where the counting optimum is never realised for diameter `> 1` because
the optimal tiles do not tessellate, the Syracuse residue walk has no
tiling obstruction at all; its same-length words stay injective far beyond
the Moore regime, and the mixing of section 2 happens only once the
exponentially many words of typical cost exhaust the `2·3^(n-1)` units
(depth `≈ n log 3/log(4/3) = 3.8 n`).

---

## 4. One-step geometric Gauss sums (the Lambert factors)

`G_j(t) := c_j sum_(r=1)^(L_j) 2^(-r) e(t 2^(-r)/3^j)`, `L_j = 2·3^(j-1)`
the order of `2` mod `3^j`, `c_j = 1/(1 - 2^(-L_j))`; this is the
one-step factor of the frequency recursion at frequency `t` when the
previous coefficient is `1`.

**Proposition 3 (PROVED).** (i) `G_j(t) = c_j sum_r 2^(-r) e(t m_r/2^r +
t/(2^r 3^j))` with `m_r = -3^(-j) mod 2^r` (the 3-adic inverse of `2^r`
read through the 2-adic inverse of `3^j`: `2^(-r) mod 3^j = (1 + m_r
3^j)/2^r`). (ii) `sup_(t unit) |G_j(t)| -> 1` as `j -> ∞`: for `t = 2^s`
with `2^s <= 3^j/2^K`, all `r <= s` give phases `e(2^(s-r)/3^j)` within
`2π 2^(-K)` of `1`, so `|G_j(2^s)| >= 1 - 2π 2^(-K) - 2^(-s)`.
FINITE-EXACT: `sup |G_j| = 0.577, 0.582, 0.789, 0.887, 0.944, 0.971,
0.986, 0.994, 0.997, 0.9985, 0.9993, 0.9997` (`j = 1..12`), at `t = 1, 8,
8, 16, 32, 64, 128, 256, 18659, 2048, 4096, 8192`; the mean of `|G_j|` over
units converges to `0.5430`; `2.9%` of the units have `|G_j| > 0.9`. The
2-adic reading agrees to `10^(-9)`.

So the frequency walk has no one-step spectral gap, exactly as the
three-quarters digraph has no one-step metric: the contraction of the
primitive coefficients is a property of the composed walk with the
renewal structure of the costs (section 2), and the powers of two are the
frequencies along which the walk is slowest.

---

## 5. Carry reciprocity and reversed cycles (the reciprocal mirror)

**Proposition 5 (PROVED).** For a word `w = (a_1..a_d)` with `A = sum
a_j`, put `C_w(u,v) = sum_(j=1)^d u^(d-j) v^(a_1+...+a_(j-1))` and
`C'_w(u,v) = sum_(j=1)^d u^(d-j) v^(a_1+...+a_j)`. Then `C_(rev w)(u,v) =
u^(d-1) v^A C'_w(1/u, 1/v)`. *Proof.* `u^(d-1) v^A C'_w(1/u,1/v) = sum_j
u^(j-1) v^(a_(j+1)+...+a_d)`; with `i = d+1-j` this is `sum_i u^(d-i)
v^(a_d + ... + a_(d+2-i)) = C_(rev w)(u,v)`. ∎ (The exclusive carry itself
is not reciprocal to the reversed one; the shift by one valuation is the
whole content. Checked on 200 random words; `w = (1,2)`: `C_w = 5`, `C'_w
= 14`, `C_(rev) = 7 = 3·8·C'_w(1/3,1/2)`.)

**Reading.** The affine identity `3^d x + C_w(3,2) = 2^A y` gives the
2-adic class of the source, `x ≡ -C_w 3^(-d) mod 2^(A+1)`, and the 3-adic
class of the target, `y ≡ C_w 2^(-A) mod 3^d`: one integer `C_w(3,2)`,
evaluated with the two prime powers in the denominator — the reciprocal
pair. Reversal exchanges the roles of the two carries and, for a cycle
(`x = y = y_0 = C_w/(2^A - 3^k)`), maps the rational cycle point `y_0` to
`y_0' = C_(rev w)/(2^A - 3^k)`, preserving `k`, `A` and the linear part
`3^k/2^A`: an involution on the rational cycles of the INTEGRAL barrier
(atlas §0). FINITE-EXACT: `-1 -> -1`; `{-5,-7}`: `(1,2) -> (2,1)`, a
rotation, `y_0' = -7`; the seven-cycle of `-17`: `(1,1,1,2,1,1,4) ->
(4,1,1,2,1,1,1)`, not a rotation, `C_(rev) = 13801`, `y_0' = -13801/139`,
a rational cycle of `3x+139`... i.e. of the map `x -> (3x + 139)/2^v` on
integers, and not an integer cycle of `3x+1`. Reversal therefore does not
preserve integrality; the integer cycles are not closed under it, which is
one more way the negative cycles are special objects of the 3-adic law
(S19, Theorem B) rather than of the word combinatorics alone.

---

## 6. Typing and what changes for the repo

* No barrier is typed by this note: it produces no termination or density
  statement. Its content is a measurement of the object the fine-scale
  mixing estimate bounds (section 2), two elementary lemmas (Propositions 4,
  5), and one negative structural fact (Proposition 3).
* **For the mixing estimate.** (2.3) ⟹ `M(h) = o(h^(-A))` for every `A`
  (Proposition 2, PROVED). The measured `M(h)` to `h = 17` decays with a
  ratio rising to `0.87`; if the ratio stays below `1`, the decay is
  geometric and (2.3) holds with a usable constant, which no published
  proof provides (Mazur's `C` is a three-fold tower, S19 audit). The
  concrete inequality to prove is `|mu_hat_h(2^s)| <= C r^h`; the cost
  decomposition names the family that carries the coefficient.
* **For the S19 seed-1 test.** Unchanged; the Fourier side explains why the
  `ell^1` distances of S19 decay slowly (bulk square-root coefficients at
  every level, constant Fourier energy per level) while remaining Cauchy.
* **For the atlas.** Proposition 4 is a residue-spread lemma of the kind
  Mazur's (3.9) uses, in the weaker same-length form; Lyu's theorem is a
  reading of the UNIFORM/DIM rows (local word constraints never bound the
  stopping time); the reversal involution is a symmetry of the INTEGRAL
  row's objects.

---

## 7. Reproduction and status table

```
python 04-computation/experiments/collatz_five_mirrors_20260929.py              > 05-knowledge/results/collatz_five_mirrors_20260929.out
python 04-computation/experiments/collatz_five_mirrors_fourier_deep_20260929.py 17 > 05-knowledge/results/collatz_five_mirrors_fourier_deep_20260929.out   # 8 GB, 1 min
python 04-computation/experiments/collatz_five_mirrors_costsplit_20260929.py 10 > 05-knowledge/results/collatz_five_mirrors_costsplit_20260929.out
python 04-computation/experiments/collatz_five_mirrors_costsplit_20260929.py 14 > 05-knowledge/results/collatz_five_mirrors_costsplit14_20260929.out   # 5 GB
```

| item | status |
|---|---|
| Proposition 1 (level independence of the primitive profile) | PROVED; checked at levels 4–17 |
| Proposition 2 ((2.3) ⟹ `M(h) <= C_A (h-1)^(-A)`) | PROVED |
| `M(h)`, `h <= 17`; argmax at `2^(h+3)`, `2^(h+4)`; ratios `0.65 -> 0.87` | FINITE-EXACT (float64 FFT); the geometric reading OBSERVED |
| Fourier mass per level `0.466`, typical `|mu_hat|^2 3^h = 0.70` | FINITE-EXACT; the identity with the second-moment slope PROVED (Parseval) |
| cost/last-valuation decomposition of the resonance at `h = 10` | FINITE-EXACT |
| Proposition 3 (Gauss sums: 2-adic reading; no uniform gap) | PROVED; table FINITE-EXACT |
| Proposition 4 (same-length spread lemma) and the collision table | PROVED; FINITE-EXACT |
| Proposition 5 (carry reciprocity); reversed cycles | PROVED; FINITE-EXACT |
| the five mirrors (section 1) | DIRECTION / remark |

**Next probes.** The joint `(Y_h, A)` decomposition at `h = 14` (memory
`3^14 × 60` doubles) to see whether the resonant cost ratio `A/h` drifts;
a large-deviation lower bound `M(h) >= c e^(-I(α) h)` from the coherent
family; the profile to `h = 18` (`20 GB`); the reversal involution on the
rational cycles of `3x + k` for small `k` (which rational cycles are
reversal-symmetric).
