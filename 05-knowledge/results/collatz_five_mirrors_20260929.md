# Five papers as structural mirrors of the 3-adic Syracuse law: the Fourier profile of Tao's Syracuse random variable (primitive maxima at the powers of two, decaying like the no-descent probability, `0.46 P_h`, on the powers of two to level 120, scale-invariant Fourier mass `0.47`), the same-length collision threshold, and carry-polynomial reciprocity under word reversal

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
2026 (Lemma 8.1, and the paragraph following it in §8.1: "primitive Fourier
coefficient bound `C* h^(-6409)` at level `h` for characters not factoring
through level `h - 1`"), Terras 1976, S19 Theorem A (harmonic mass), the
S19 Fourier consistency (`Y_n mod 3^m` has the law of `Y_m`), THM-4476 and
THM-4487 (the thin-divergence exponent `h* = h(log_3 2) = 0.94996` and the
no-descent tilt).

**Status: PROVED (Propositions 1–5: level-independence of the primitive
Fourier profile; the implication (2.3) ⟹ super-polynomial decay of the
primitive maxima; no uniform one-step gap of the geometric Gauss sums; the
same-length spread lemma; carry reciprocity under reversal; Proposition 6,
the reversal-invariant trace; the closure of the powers of two under the
frequency recursion) + VERIFIED (the Fourier profile of the law to level
18 in float64, agreeing with the rational law to level 6, with an
independent forward recursion to level 14 to `4·10^(-7)`, and with the
closed recursion on the powers of two to six digits) + FINITE-EXACT (the
collision table; the reversed cycles; the Gauss-sum table; the traces) +
OBSERVED (the maxima sit at the powers of two `±2^s`, `s - h` growing
slowly; the decay ratio rises from `0.65` to `0.89` over eighteen levels;
the Fourier mass per conductor level is `0.46–0.47`; on the powers of two
the coefficient is `0.44–0.48` times the no-descent probability over the
levels `20..120`) + CONJECTURAL (`M(h) ≍ P_h`: the maximal primitive
coefficient decays at the no-descent rate `3^(h*-1)`, a necessary condition
for (2.3), not the mixing estimate itself) + DIRECTION. No Collatz proof step. Audited SOUND WITH CORRECTIONS (eighteen applied; section 8); section 2b, the `3x+k` reversal survey and Proposition 6 were added after the audit's snapshot, and section 2b was then audited separately (fifteen corrections applied; section 8).**

Scripts and outputs:
`04-computation/experiments/collatz_five_mirrors_20260929.py` →
`05-knowledge/results/collatz_five_mirrors_20260929.out` (profile to level
14, Gauss sums, collisions, reciprocity);
`collatz_five_mirrors_fourier_deep_20260929.py` → `..._fourier_deep_...out`
(profile to level 17; `8 GB`, one minute);
`collatz_five_mirrors_costsplit_20260929.py` → `..._costsplit_...out` and
`..._costsplit14_...out` (the resonant coefficient decomposed by total
cost and by last valuation at levels 10 and 14);
`collatz_five_mirrors_fourier_deep_20260929.py 18` → `..._fourier_deep18_...out`
(the profile to level 18); `collatz_five_mirrors_reversal_20260929.py` →
`..._reversal_...out` (reversal on the integer cycles of `3x+k`);
`collatz_five_mirrors_tracesum_20260929.py` → `..._tracesum_...out` (the
reversal-invariant trace, Proposition 6).

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
      than any polynomial. Computed to `h = 18` (VERIFIED, float64):
   `0.5774, 0.3779, 0.2522, 0.1770, 0.1293, 0.0961, 0.0759, 0.0609, 0.0480,
   0.0383, 0.0319, 0.0265, 0.0221, 0.0191, 0.0163, 0.0144, 0.0125, 0.0112`;
   the maximum sits at `u = ±2^s` with `s - h = 0` for `h <= 6` and `s - h
   = 1, 2, 3, 4, 5` on `h = 7–9, 10–12, 13–15, 16–17, 18` (OBSERVED; `2^17`
   at `h = 14`, `2^23` at `h = 18`); the ratio `M(h)/M(h-1)` rises from
   `0.65` to `0.89` with no plateau (window means `0.78, 0.82, 0.87` over
   `h = 6–9, 10–13, 14–18`). Eighteen levels do not separate a geometric
   decay from a shifted power law (`C (h + 2.5)^(-2.5)` fits `7–18` to
   `1.8%`, the audit's finding; the readings separate near `h = 31`). The
   separation comes from the powers of two themselves: they are a closed
   family under the exact frequency recursion (`2^j -> 2^(j-a)`), so
   `mu_hat_h(2^j)` is computable without the law (section 2b). That
   recursion reproduces every FFT value to six digits and continues to `h =
      120`, where `max_j |mu_hat_h(2^j)|` decays geometrically up to a
   polynomial prefactor (`C h^(-1.1) r^h`, `r ≈ 0.944`; level ratio `0.906
   -> 0.930` over `h = 30..120`; local exponents `2.1, ..., 6.1` on
   successive doublings, growing linearly with slope `0.083 ≈ |log_2
   0.9465|`, the signature of a geometric    law) and is `0.44–0.48` times the
   no-descent probability `P_h = P(a_1 + ... + a_j < j log_2 3 for all j <=
   h)` over the levels `20..180` (`0.44–0.58` over `20..300`, section 2c). `P_h` decays at the rate `3^(h* - 1) =
   0.9465` (PROVED, standard large deviations), `h* = h(log_3 2)` the
   thin-divergence exponent of THM-4476: on the powers of two, which carry
   the maximum wherever that was checked, the    maximal Fourier coefficient
   decays like the no-descent probability (OBSERVED to `h = 300`, with the
      exponential rate `0.947–0.950` by the admissible fits, the no-descent
   rate at its lower edge); that the
   decay rate of `M(h)` is the no-descent rate `3^(h*-1)` is CONJECTURAL
   (it needs the maximum to stay on the powers of two), and it concerns the
   maximal primitive coefficient, not the `ℓ^1` distances of (2.3). This is
   consistent with (2.3) (geometric beats every polynomial) and is the
   first measurement of the necessary condition that Proposition 2 extracts
   from the mixing estimate.
2. **Scale invariance of the Fourier energy.** The Fourier mass at conductor
   level `h`, `sum_(cond t = 3^h) |mu_hat(t)|^2`, is `0.667, 0.476, 0.462`
   and then `0.464 .. 0.472` for `h = 4..18` (slowly increasing); with
   Proposition 1 the mass at level `h` is exactly `(3/2)(E_units[rho_h^2] -
   E_units[rho_(h-1)^2])`, and the S19 second-moment increments `0.308 ->
   0.314` give `0.462 -> 0.472`: every 3-adic scale carries nearly the same
   energy, which is the linear growth of the second moment found in S19.
   The typical coefficient at level `h` has `|mu_hat|^2 3^h = 0.70`:
   square-root cancellation, as for random phases; the slow sup-norm decay
   comes from the powers of two alone.
3. **Where the resonance comes from.** At `h = 10` the resonant coefficient
   (`t = 2^12`, `0.0383`, against `0.0026` for a generic unit) is carried by
   the words of total cost `A = 12..17` (`60%` of `|.|` coherent), far below
   the typical cost `2h = 20`, and by the last valuations `1, 2, 3`
      (`0.0225, 0.0120, 0.0051`) (FINITE-EXACT at `h = 10` and `14`; at `h =
   14`: `A = 19..23`, coherent fraction `0.56`, last valuation `1` carrying
   `57%`). The contribution-weighted cost ratio is stable at the two levels,
   `A/h = 1.49` and `1.48`, and the cost family `A ≈ 1.48 h` has
   large-deviation mass `e^(-0.093 h)` (ratio `0.91` per level). The reading
   that the low-cost coherent words carry the coefficient is DIRECTION at
      these levels and is consistent with section 2b: to `h = 120` the
   coefficient at the powers of two is proportional to the no-descent
   probability, whose critical family has `A/h -> log_2 3 = 1.585` (the
   resonant exponent is `s = h log_2 3 - 6 ± 1` there).
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
   classes mod `3^n`. In the tree of `1` with valuations `<= 24` the first
   same-length collisions mod `3^n` appear at depths `2..4` with minimal
   cost sums `29..63` for `n = 4..12`, `4.4–7` times the bound; with
   valuations `<= 70` collisions appear at depth `2` already for `n <= 9`
   and the minimal colliding sum at depth `d` decreases with `d`, but at
   every depth `d <= 5` and `n <= 12` it exceeds the bound by a factor of at
   least `4.4`: the Syracuse residues stay injective well beyond the Moore
   regime, and no counterpart of the paper's tessellation obstruction was
   found (an analogy, not a theorem).
6. **Reciprocity (Proposition 5, PROVED):** with the exclusive-prefix carry
   `C_w(u,v) = sum_j u^(d-j) v^(a_1+...+a_(j-1))` (the affine identity `3^d x
   + C_w = 2^A y`) and the inclusive one `C'_w`, word reversal is polynomial
      reciprocity, `C_(rev w)(u,v) = u^(d-1) v^A C'_w(1/u, 1/v)`; the 2-adic
   source `x ≡ -C_w 3^(-d) mod 2^A` and the 3-adic target `y ≡ C_w 2^(-A)
   mod 3^d` are the two reciprocal evaluations of one integer. Reversal is an involution on
   rational cycles that preserves `(k, A)`; it fixes the cycles of `-1` and
   `{-5,-7}` (a rotation) and sends the seven-cycle of `-17` to the rational
   point `-13801/139`, not an integer (FINITE-EXACT). On the maps `3x + k`,
   `k <= 41` odd, reversal permutes the integer cycles of fixed `(d, A)`:
   its fixed points are exactly the rotation-symmetric words, and it has
   genuine pairs — the 5-cycles of `3x+13` with minima `227 ↔ 259` and `251
   ↔ 287`, the 3-cycles of `3x+37` with minima `23 ↔ 29` (and their triples
   for `3x+39`); the seven-cycle family and the long cycles of `3x+5` (`d
   = 17`) and `3x+13` (`d = 15`) have rational partners (FINITE-EXACT,
   orbits verified; section 5). The trace of a rational cycle is
   reversal-invariant (Proposition 6, PROVED): the pairs have equal element
   sums (`2499`, `125`), and the seven-cycle and its rational partner both
   sum to `-327`.
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
3-adic Syracuse law is the stationary law of this IFS (S19 §§2, 7), and the
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
its layer recursion `H_(n+1) = H^(1) + 2H^(2)` (S19, Corollary A2) is the analogue of Merca's linear recurrence (his Corollary 1.5). Verdict:
remark; `mu_n` is positive by definition (a sum of the weights `2^(-A(w))`),
so Merca's positivity phenomenon has no counterpart, and the open `liminf
H_n(1) > 0` is an asymptotic question about a positive sequence; a closed
formula for `H_n(1)` is not in reach (it would be a formula for the tree of
`1`).

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
`2^(-1)`, all steps `a >= 1` allowed (the parity restriction `a ≡ ε(z) mod
2` lives on the space side, in the parents' recursion of S19 §5); the ball
growth is exponential
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

**The profile (VERIFIED, float64 FFT of the level-18 law, `20 GB`;
identical at every level where compared, and reproduced to six digits by
the closed recursion of section 2b).**

| `h` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `M(h)` | .5774 | .3779 | .2522 | .1770 | .1293 | .0961 | .0759 | .0609 | .0480 | .0383 | .0319 | .0265 | .0221 | .0191 | .0163 | .0144 | .0125 | .0112 |
| ratio | – | .655 | .667 | .702 | .730 | .743 | .789 | .803 | .789 | .797 | .835 | .828 | .834 | .867 | .851 | .885 | .868 | .894 |
| argmax `u` | 1 | 2^2 | 2^3 | 2^4 | 2^5 | 2^6 | 2^8 | 2^9 | 2^10 | 2^12 | 2^13 | 2^14 | 2^16 | 2^17 | 2^18 | 2^20 | 2^21 | 2^23 |
| Fourier mass at level `h` | .667 | .476 | .462 | .464 | .466 | .466 | .465 | .466 | .467 | .467 | .468 | .469 | .470 | .470 | .471 | .471 | .472 | .472 |
| typical `|mu_hat|^2 3^h` | 1.00 | .714 | .692 | .696 | .698 | .699 | .698 | .699 | .700 | .701 | .702 | .703 | .704 | .705 | .706 | .707 | .708 | .708 |

`M(1) = 1/sqrt 3` exactly (`|(1/3) e(1/3) + (2/3) e(2/3)|`). The maximum
is always at a power of two up to sign (`|mu_hat(-u)| = |mu_hat(u)|`),
`u = ±2^s` with `s - h` stepping up by one every two or three levels (`s -
h = 1` at `h = 7`, `5` at `h = 18`), so `2^s/3^h` decreases (`0.12` at `h = 7`, `0.027` at `h = 14`, `0.022`
at `h = 18`) until it settles near `2^(-6)` (section 2b: `s = h log_2 3 -
6 ± 1` for `20 <= h <= 120`): the
character `e(2^s y/3^h)` is not a low frequency in the real sense; it is
the frequency that turns the 3-adic inverses `2^(-a)` of the small
valuations into small integers `2^(s-a)` (the one-step Gauss sums of
section 4 show the same alignment).
The values `|mu_hat_17(2^s)|`, `s = 0..60`, form a single bump centred at
`s = 21` (`0.0125`), with `0.0111` at `s = 20`, `0.0106` at `23`, and the
generic size `7·10^(-5)` (`= 0.84 · 3^(-17/2)`) outside `s = 13..27`; at
level 18 the bump is centred at `s = 23` (`0.0112`).

**Decay (levels `<= 18`).** Local exponents from the doublings
`M(h)/M(2h)`: `1.76` (5→10), `1.86` (6→12), `1.99` (7→14), `2.08` (8→16),
`2.10` (9→18): steepening, so no pure power law `C h^(-α)`; but a shifted
power law `C (h + 2.5)^(-2.52)` fits `7 <= h <= 18` to `1.8%` (the audit's
fit; local exponents `1.86, 1.94, 2.01, 2.06, 2.10`, ratios `0.854–0.881`
at `14..18`), better than a geometric `C r^h` (`r = 0.86`, `4.3%` on
`11–18`), and it predicted the level-18 value from the levels below
(`0.0110` against the geometric `0.0105`; measured `0.0112`). `3.75
h^(-2)` fits `7 <= h <= 14` to `4%` and `M(h)` lies `2–4%` below it at
`15–18`. The ratios rise from `0.65` to `0.89` and are still rising over
`14..18` (mean `0.873` against `0.823` over `10..13`; the `±0.02`
fluctuation follows the steps of `s - h`). So the eighteen FFT levels do
not decide between a geometric decay (compatible with Proposition 2) and a
polynomial one (which would contradict it); the two readings separate by a
factor `2` near `h = 31`, and Tao's constants say nothing below `m` of the
order of `2^(2^(2^8697))`. Section 2b decides it for the powers of two.

**Scale invariance.** By Parseval, `sum_t |mu_hat_n(t)|^2 = 3^n sum_y
mu_n(y)^2 = (3/2) E_units[rho_n^2]`, and with Proposition 1 the mass at
level `h` is exactly `(3/2)(E_units[rho_h^2] - E_units[rho_(h-1)^2])`; the
second-moment increment is `0.308 -> 0.314` over `h = 3..17` (S19: `0.31`),
so the mass per level is `0.462 -> 0.472`, slowly increasing. Every 3-adic
scale carries nearly the same Fourier energy: the law is "`1/f`" in the 3-adic scale, its density is in no `L^p`
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
level) but phase-coherent; at `h <= 18` the sup-norm ratio `0.87` lies below this
mass rate `0.91` (a coherence loss, if the decay is geometric; DIRECTION,
the quotient `0.87/0.91` is not a measurement). Section 2b shows, to `h =
120`, that the coefficient at the powers of two is `0.46` times the
no-descent probability, whose critical family has `A/h -> log_2 3 = 1.585`
and rate `0.9465`: the finite-level family `A ≈ 1.48 h` is the beginning of
that critical family. This is the candidate mechanism behind the slow
sup-norm decay and the reason it must be renewal-theoretic (Tao's
Fourier–renewal method, Mazur's §8.1 "renewal bounds"): no single step has
a gap (section 4), and the obstruction is a large-deviation family of
words, whose contribution is what a sharp version of (2.3) would have to
bound.

**Direction.** A proof that `|mu_hat_h(2^s)| <= C r^h` for some `r < 1`
uniformly in `s`, together with the square-root behaviour of the generic
coefficients, would give (2.3) with a usable constant; the cost
decomposition suggests the route: bound the coherent low-cost family by its
large-deviation mass and show the remaining words cancel at the
square-root rate. This is the analytic core of Tao's proposition restated
as a concrete inequality on one explicit sequence.

**2b. The powers of two to level 120 (added after the audit's snapshot;
`collatz_five_mirrors_powers_of_two_20260929.py`).** The frequency
recursion `mu_hat_n(t) = sum_(a>=1) 2^(-a) e((t 2^(-a) mod 3^n)/3^n)
mu_hat_(n-1)(t 2^(-a) mod 3^(n-1))` sends `t = 2^j` to the frequencies
`2^(j-a)`: the family `{2^j mod 3^n : j ∈ Z}` (the inverse powers included)
is closed under it (PROVED, trivial). With the valuation truncated at `a <=
40` (error `<= 2^(-40)` per level) and `mu_hat_0 ≡ 1`, the values `m_n(j)
:= mu_hat_n(2^j mod 3^n)` follow for every `j` and `n` from a recursion
over about `5000` exponents per level, in Python big integers for the
residues and float64 for the phases, without the law itself. Validation:
`max_j |m_h(j)|` reproduces every FFT value `M(h)`, `h <= 18`, to six
digits (`0.577350, 0.377924, ..., 0.011187`) with the same argmax `s - h =
1..5` for `9 <= h <= 18` (for `h <= 8` the family covers all units and
the maximum is `M(h)` by definition). Beyond, `max_j |m_h(j)|` is a lower
bound for `M(h)` (an equality if the maximum stays on the powers of two,
true wherever checked):

| `h` | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 | 110 | 120 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `max_j |m_h(j)|` | 8.88e-3 | 3.29e-3 | 1.38e-3 | 6.16e-4 | 2.83e-4 | 1.33e-4 | 6.45e-5 | 3.19e-5 | 1.59e-5 | 8.06e-6 | 4.12e-6 |
| argmax `s - h` | 6 | 12 | 17 | 23 | 29 | 35 | 40 | 46 | 52 | 58 | 64 |
| ratio to the level below | .905 | .906 | .919 | .931 | .926 | .921 | .929 | .938 | .938 | .934 | .930 |
| `P_h` (no-descent probability) | 1.93e-2 | 7.07e-3 | 2.99e-3 | 1.33e-3 | 6.15e-4 | 2.83e-4 | 1.40e-4 | 6.94e-5 | 3.45e-5 | 1.79e-5 | 9.31e-6 |
| `max_j |m_h(j)| / P_h` | .461 | .465 | .463 | .463 | .460 | .470 | .462 | .459 | .462 | .449 | .443 |

(`P_h = P(a_1 +
... + a_j < j log_2 3` for all `j <= h)` for i.i.d. geometric(1/2)
valuations, by dynamic programming.) The decay is geometric with a
polynomial prefactor: `C h^(-1.1) r^h`, `r = 0.944`, fits `h = 40..120` to
`0.9%` (a pure geometric fit has ratio `0.930` and residual `12%`; a
shifted power law is not a well-posed competitor over a bounded range,
since with a large shift it tends to a geometric law), and the local
exponents on successive doublings, `2.1` (10→20), `2.3` (15→30), `2.7`
(20→40), `3.1` (25→50), `3.5` (30→60), `4.4` (40→80), `5.3` (50→100),
`6.1` (60→120), grow linearly in `h` with slope `0.083` per level and
intercept `1.1`, as `h |log_2 r| + β` for `C h^(-β) r^h` (`|log_2 0.9465|
= 0.079`); a shifted power law would saturate them at its exponent. The
audit's shifted power law, which fit the eighteen FFT levels better than a
geometric law, fails beyond `h ≈ 30`. The resonant exponent is `s = h log_2
3 - 6 ± 1` for every `20 <= h <= 120`: the resonant frequency is a fixed
fraction `2^s/3^h ≈ 2^(-6)` (`0.009–0.019`) of the modulus, and `s/h`
approaches `log_2 3 = 1.585` (`1.53` at `h = 120`).

*Is the maximum on the pure powers of two?* Every family `u 2^j` (`u` a
fixed unit) is closed under the recursion in the same way, so the
comparison can be run for the small odd multipliers `u <= 49`, `3 ∤ u`
(seventeen families; `collatz_five_mirrors_multiplier_families_20260929.py`):
at every level `h <= 60` the pure family `u = 1` gives the largest
coefficient, with `pure/max = 1.0000` throughout (VERIFIED; the other
sixteen families never exceed it). Units with larger multipliers are not
covered, so "the maximum stays on the powers of two" remains OBSERVED (to
`h = 18` over all units, to `h = 60` over these families).

**The no-descent proportionality (OBSERVED).** `max_j |m_h(j)| / P_h`
lies in `0.44–0.48` at every level `20 <= h <= 120` (window means `0.459,
0.465, 0.456` over `20–40, 41–80, 81–120`; `0.481` at `h = 41`, `0.443` at
`120`, the coefficient's prefactor decaying slightly faster than `P_h`'s).
The rate of `P_h` is the large-deviation rate of the valuation sum at the
critical slope `log_2 3`: `P_h^(1/h) -> e^(-I(log_2 3)) = 3^(h* - 1) =
0.94650` with `h* = h(log_3 2) = 0.94996` the binary entropy of `log_3 2`,
the entropy that governs the thin divergent orbits of THM-4476 (`N(X) <= K
X^(h*+ε)`) and the tilt of THM-4487 (PROVED, standard: `P_h <= P(S_h < h
log_2 3) <= e^(-I h)` by Chernoff, and the matching lower bound by tilting
the valuations to mean `log_2 3 - ε`; the cheapest way to stay under the
critical line is to walk along it). The convergence is slow because of a
polynomial prefactor, `P_h e^(I h) ≈ h^(-1.3)` over `40..120` (expected
`h^(-3/2)`): the single-level ratios `P_h/P_(h-1)` oscillate with the
integer part of `j log_2 3` (`0.903` at `70`, `0.949` at `120`), and the
geometric-mean ratio is `0.9325` over `61..120` and `0.9366` over
`101..120`, approaching `0.9465` from below. Reading (DIRECTION): the
maximal Fourier coefficient of the 3-adic Syracuse law at the powers of
two is, up to a constant near `0.46`, the probability that a random
valuation word of length `h` never descends — the coherent family is the
no-descent set, whose members' phases `e(2^(s-T_j)/3^j)` (`T_j` the suffix
sums, `2^(s-T_j)` small against `3^j` on the no-descent set) stay near `1`,
the last few valuations fixing the constant, while the descending words
cancel. Consequences: (i) by Proposition 2's inequality the `ℓ^1` distances
of (2.3) satisfy `||mu_q - lift mu_m||_1 >= max_j |m_(m+1)(j)| ≈ 0.46
P_(m+1)` (PROVED / OBSERVED), so (2.3) cannot hold with a rate faster than
the no-descent rate; whether it holds with that rate is a statement about
the bulk coefficients, which dominate the distances (section 2), and is
not decided here; (ii) the decay rate of the maximal primitive coefficient
and the thin-divergence exponent would be the same number (CONJECTURAL for
`M(h)`: it needs the maximum to stay on the powers of two, and a proof of
the square-root cancellation over the descending words); (iii) Mazur's
`C* h^(-6409)` would then be a polynomial statement about a geometric
quantity. What is PROVED here: the closure and exactness of the recursion,
`M(h) >= max_j |m_h(j)|`, the FFT agreement to `h = 18`, the rate of
`P_h`; what is OBSERVED: the decay law and the proportionality to `h =
120`.

**2c. The obligations worked (S21, 2026-09-29; added after both audits;
`collatz_five_mirrors_powers_of_two_20260929.py 300`,
`collatz_five_mirrors_rate300_20260929.py`,
`collatz_five_mirrors_coherent_20260929.py`,
`collatz_five_mirrors_coherent120_20260929.py`).**

*The rate to level 300 (VERIFIED numerics, OBSERVED law).* The closed
recursion continued to `h = 300` (`M(300) = 7.42·10^(-11)`; exponent window `12000` per level; the truncation `a <= 40` is not self-certifying here, its a priori bound `300·2^(-40) = 2.7·10^(-10)` exceeding `M(300)`; the audit's `a <= 60` run certifies the values to `3·10^(-4)` relative and agrees with the `a <= 40` run to `2·10^(-12)`): the local exponents on doublings keep growing
linearly, `3.10, 5.27, 7.27, 9.19, 11.07, 12.97` for `25→50, 50→100,
75→150, 100→200, 125→250, 150→300`, least-squares slope `0.0785` per level (two-point `0.077` between the doublings `50→100` and `150→300`) against `|log_2 3^(h*-1)| = 0.079`, i.e. `r = 0.947`; the geometric-mean level ratio of `M` is
`0.9383, 0.9418, 0.9427, 0.9433` over `100–200, 150–300, 200–300,
250–300`, that of `P_h` `0.9375, 0.9404, 0.9413`, both rising toward
`0.9465`. A free fit `C h^(-β) r^h` over `100..300` gives `β = 1.62, r =
0.9489` for `M` (residual `1.0%`) and `β = 1.29, r = 0.9461` for `P_h`
(`2.3%`); with `r` pinned to `3^(h*-1) = 0.9465` the fits have `β = 1.14`
(`M`, residual `5.9%`) and `β = 1.38` (`P_h`, `2.5%`): over a finite range
`β` and `r` trade against each other, so the exponential rate of `M` is `0.947–0.950` by the admissible estimators (fits with residual at most `2.5%` have `β = 1.5–1.75` and `r = 0.948–0.950`; the least-squares doubling slope `0.0785` gives `0.947`; the pure geometric `0.941` has residual `20%`, and the fit pinned to `3^(h*-1) = 0.9465` needs `β = 1.14` at six times the free fit's residual): the no-descent rate sits at the lower edge of the bracket, `0.15–0.35%` below the fitted rates, a gap that 300 levels cannot attribute to the prefactor or to the rate; the identification of the rate with the no-descent rate stays OBSERVED only in this weaker sense. The ratio `M/P_h` stays in `0.44–0.48` up to `h = 180` and then drifts up (`0.500,
0.519, 0.560, 0.579` at `h = 200, 240, 280, 300`): over `20 <= h <= 300` the ratio moves within `0.44–0.58`; the rise by a factor `1.26` over `180..300` is what a rate `0.2%` per level above the no-descent rate produces, and whether `M` and `P_h` have the same exponential order (the conjecture `M ≍ P_h`) or `M` decays slightly slower is not decided by these levels, which lean slightly against the conjecture; "`0.46 P_h`" is a description of the first two hundred levels, not an identity. The resonant exponent is `s = h log_2 3 - 5.3 .. 7.2` at every level `20..300` (`-6.2 .. -7.1` at the multiples of 25; `-6.49` at `h = 300`).

*The mechanism, tested directly (VERIFIED numerics; the reading
DIRECTION).* A dynamic programme over `(total cost A, level j, prefix sum
P)` with the exact phase factors `e((2^(s - A + P_(j-1)) mod 3^j)/3^j)`
splits the resonant coefficient by the excursion of the prefix sums above
the critical line, `E(w) = max_j (P_j - j log_2 3)`. The strict no-descent words (`E < 0`) carry `54%, 41%, 35%` of `|mu_hat_h(2^s)|` at `h = 40, 80, 120`, a share falling by about `0.1` per forty levels (their own coherent fraction is `0.25, 0.19, 0.15`: their phases are
spread), so the descending words do not cancel; they add roughly in phase.
The words with `E < c` reproduce the coefficient as a vector to `10%` from `c = 5, 9, 15` on (to `5%` from `c = 6, 10, 20`) at `h = 40, 80, 120`: the band width grows roughly like `h/8`. The magnitude alone is within `4%` for every `c >= 5` (`|coh_6|/|full| = 1.008, 1.007, 0.997`; at `c = 4`: `0.92, 0.93, 0.89`), but at `h = 120` it wanders up to `1.12` (`c = 12`) while the vector remainder is `0.19–0.29` for `c = 5..9` (`c = 16`: `0.065`; `c = 24`: `0.002`): the deep-descending words' net contribution is a rotation whose size decays with `c` more slowly as `h` grows. Within the band the phases are far from aligned (coherent
fraction `0.02, 0.015, 0.013` at `c = 6`): the coefficient is a small
residue of the band's mass (`P_h^(6)/P_h = 24, 30, 34`), not a sum of
aligned terms. The typical coefficient is at the square-root scale: random units at `h = 30` have rms `|mu_hat| ≈ 5·10^(-8)` (twelve units: `0.8·10^(-8) .. 8·10^(-8)`, rms `3·10^(-8)`; the audit's hundred: rms `5.4·10^(-8)`) against the level average `0.84·3^(-15) = 5.9·10^(-8)`; the median is only `0.4` of the rms (`2.3·10^(-8)`) because the coefficient distribution over the units is skewed (exact at `h = 14`: median `0.62` rms, `38%` of the units below half the rms), `10^5` below the resonant window (`3.3·10^(-3)`).

*The excursion profile (added after the third audit's snapshot;
`collatz_five_mirrors_excursion_profile_20260929.py`).* The same dynamic
programme with the maximal excursion bin `m(w) = max_j floor(P_j - j log_2
3)` as a state gives the contribution of each bin. At `h = 40, 80, 120`
the cumulative sum over `m <= 4` reaches `0.98, 0.99, 0.98` of the
coefficient's magnitude, with vector remainders `0.056, 0.053, 0.195`; the
per-bin contributions decay with `m` while the bin masses grow (the bin
`m = 20` holds `1.4, 14, 35` times the strict no-descent mass and
contributes `3·10^(-4), 3·10^(-3), 2·10^(-2)` of the coefficient), and
the bins beyond `m = 4` carry a rotating remainder of `5–20%` whose
cumulative effect decays with `m` non-monotonically (`h = 120`: `0.19,
0.29, 0.19, 0.13, 0.09, 0.024` at `m = 4, 6, 8, 10, 14, 20`). So the
carrier of the coefficient is the near-critical band `m <= 4` at every
level, and the width needed to bring the remainder below `10%` grows
slowly with `h` (`4, 3, 14`, the criterion being sensitive to the phase
rotation of the deep bins).

*What this leaves as the proof-shaped targets.* With `s = h log_2 3 -
6`, the phase of a word at level `j` is `e((2^(s - T_j) mod 3^j)/3^j)`,
`T_j` the suffix cost, which is near `1` iff `A - h log_2 3 + 6 <= P_(j-1) <= A - h log_2 3 + j log_2 3 + 6 - K` (on the band `A - h log_2 3` lies between `-∞` and `c`; below the lower end the exponent is negative and the phase is a scrambled inverse power): (T1) *cancellation of deep descents* — the sum of `2^(-A) e(2^s Y_h/3^h)` over the words whose excursion exceeds `c` is at most `ε(c)` times the band sum, with `ε(c) -> 0` uniformly in `h` (for fixed `h` the statement is empty; the data show the uniform version needs `c` of order `h/8` on `40..120`, so the target may have to be stated with `c = c(h) = o(h)`); (T2) *the band residue* — the sum over the words with excursion below the same `c(h)` has modulus of the exponential order of `P_h` (its mass is `≍ P_h` by large
deviations; the content is that the spread phases leave a residue of the
same order). (T1) is a statement about the scrambling of `2^m mod 3^j`
for `m` beyond `j log_2 3`, the same equidistribution that Tao's
Fourier–renewal method organises; (T2) is a non-cancellation statement
inside the near-critical band. Both are OPEN; neither is a Collatz step by
itself, but together they would turn the necessary condition of
Proposition 2 into a two-sided description of the sup-norm side of the
mixing estimate, with the thin-divergence exponent `h*` as its rate.

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
C_(w')`, and `(d, C_w, A)` determines `w` (`C_w` does not involve `a_d`:
`a_1 = v_2(C_w - 3^(d-1))`, recurse on `(C_w - 3^(d-1))/2^(a_1)` for `a_2,
..., a_(d-1)`, then `a_d = A - a_1 - ... - a_(d-1)`). And `0 < |C_w 2^(A')
- C_(w') 2^A|
< 2^(A+A') 3^d / 2` since `C_w <= 2^A (3^d - 1)/2`. Hence `3^n < 2^(A+A')
3^d/2`. ∎

**The actual collisions in the tree of `1`** (nodes of depth `d`, reduced
mod `3^n`; the session's run used valuations `<= 24` and reported the
minimal `A + A'` at the first colliding depth, `29, 31, 32, 32, 56, 56, 58,
61, 63` at `d = 2, 2, 3, 3, 3, 3, 4, 4, 4` for `n = 4..12`; the audit
re-ran the search with valuations `<= 70` and found that the first
colliding depth and the minimal sum both depend on the cap). Minimal
colliding `A + A'` at depth `d` with valuations `<= 70` (exact for every
entry `<= 70 + 3d`; `*` = upper bound; `–` = none found):

| `n` | `d = 2` | `d = 3` | `d = 4` | `d = 5` | bound at `d = 2` | bound at `d = 3` |
|---|---|---|---|---|---|---|
| 4 | 29 | 22 | 26 | 29 | 4.2 | 2.6 |
| 5 | 31 | 32 | 29 | 29 | 5.8 | 4.2 |
| 6 | 49 | 32 | 29 | 33 | 7.3 | 5.8 |
| 7 | 73 | 32 | 36 | 40 | 8.9 | 7.3 |
| 8 | 99* | 56 | 43 | 44 | 10.5 | 8.9 |
| 9 | 99* | 56 | 55 | 44 | 12.1 | 10.5 |
| 10 | – | 70 | 58 | 50 | 13.7 | 12.1 |
| 11 | – | 80* | 61 | 52 | 15.3 | 13.7 |
| 12 | – | 127* | 63 | 59 | 16.9 | 15.3 |

At every depth `d <= 5` and `n <= 12` the minimal colliding `A + A'`
exceeds the bound `(n - d) log_2 3 + 1` by a factor of at least `4.4`
(`4.4–7` on the first-depth entries). Unlike the three-quarters digraphs,
where the counting optimum is never realised for diameter `> 1` because
the optimal tiles do not tessellate, no tiling obstruction of the Syracuse
residue walk was found (an analogy, not a theorem); its same-length words
stay injective well beyond the Moore regime, and the depth-`d` layers with
valuations `<= 70` already cover every unit class mod `3^n` from `d ≈ n/2`
on for `n <= 8` (the entropy count `4^d ≈ 3^n` gives `d ≈ 0.79 n`).

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
`2π 2^(-K)` of `1`, so `|G_j(2^s)| >= (1 - 2^(-s)) cos(2π 2^(-K)) - 2^(-s)
>= 1 - 2π 2^(-K) - 2^(1-s)` (head of modulus at least `(1 - 2^(-s))
cos(2π 2^(-K))`, tail at most `2^(-s)`; with `2^(-s)` in place of
`2^(1-s)` the bound fails, e.g. `|G_9(2)| = 0.207`).
FINITE-EXACT: `sup |G_j| = 0.577, 0.582, 0.789, 0.887, 0.944, 0.971,
0.986, 0.994, 0.997, 0.9985, 0.9993, 0.9997` (`j = 1..12`), at `t = ±2^(j+1)
mod 3^j` (`1, 8, 8, 16, 32, 64, 128, 256, 18659 = -2^10, 2048, 4096,
8192`); the mean of `|G_j|` over
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
2-adic class of the source, `x ≡ -C_w 3^(-d) mod 2^A` (the class mod
`2^(A+1)` is `(2^A - C_w) 3^(-d)`, `y` being odd), and the 3-adic
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

**Reversal on the integer cycles of `3x + k` (FINITE-EXACT; added after
the audit's snapshot, `collatz_five_mirrors_reversal_20260929.py`).** For
the map `x -> (3x + k)/2^v` on odd integers the affine identity is `3^d x
+ k C_w = 2^A y`, the cycle point of a word is `y_0 = k C_w/(2^A - 3^d)`,
and the reversed word's point is `y_0' = k C_(rev w)/(2^A - 3^d)`; the
survey takes every cycle reached from the odd starts `|x| <= 5·10^4`
(orbits capped at `10^15`), reverses its word from its minimal element,
and verifies the reversed orbit step by step (integrality of the affine
fixed point alone is not enough: the valuations must be exactly the
reversed ones). Result for odd `k <= 41`: reversal permutes the integer
cycles of fixed `(d, A)`; its fixed points are exactly the cycles whose
word is a rotation of its reverse (`RIV` in the output), and it has genuine
transpositions: for `3x+13` the 5-cycles of cost `A = 8` pair as `{227,
347, 527, 797, 601}` (word `(1,1,1,2,3)`) `↔ {259, 395, 599, 905, 341}`
(word `(3,2,1,1,1)`, read from `905`) and `251 ↔ 287`, while `211, 283,
319` are self-dual; for `3x+37` the 3-cycles of cost `6` pair as `{23, 53,
49}` (`(1,2,3)`) `↔ {29, 31, 65}` (`(2,1,3)`); `3x+39` repeats the `3x+13`
pairs tripled (`681 ↔ 777`, `753 ↔ 861`). Every cycle of `3x+1` is
self-dual except the seven-cycle, whose partner is rational; the long
cycles of `3x+5` (`187, 347`, `d = 17`) and `3x+13` (`131`, `d = 15`)
also have rational partners, and the `k`-multiples of the `3x+1` cycles
appear for every `k` (`3(kx) + k = k(3x + 1)`). Reversal is thus a
symmetry of the rational cycle set (the INTEGRAL row's objects) that
sometimes descends to the integers: an involution on the integer cycles
with fixed points and pairs, and the seven-cycle of `-17` is the smallest
`3x+1` orbit it moves off the integers.

**Proposition 6 (the trace of a cycle is reversal-invariant; PROVED, added
after the audit's snapshot).** For a word `w` let `S_w(u,v) = sum_j
C_(rot_j w)(u,v)` over the `d` cyclic rotations. Then `S_(rev w) = S_w` as
polynomials; consequently the elements `y_j = k C_(rot_j w)(3,2)/(2^A -
3^d)` of the rational cycle of `w` for `x -> (3x+k)/2^v` and those of the
reversed cycle have the same sum. *Proof.* `S_w(u,v) = sum_(i=1)^d
u^(d-i) B_(i-1)(w; v)`, where `B_m(w; v) = sum_j v^(b_(j,m))` and `b_(j,m)`
is the sum of the `m` cyclically consecutive entries starting at position
`j`; the multiset `{b_(j,m)}_j` is the same for `w` and `rev w` (a
cyclically consecutive block of `rev w` is a reversed block of `w`), so
every `B_m` is reversal-invariant. ∎ Checked on 300 random words with
rational `(u, v)` (`collatz_five_mirrors_tracesum_20260929.py`): the pairs
above have equal sums, `227 + 347 + 527 + 797 + 601 = 259 + 395 + 599 +
905 + 341 = 2499` and `23 + 53 + 49 = 29 + 31 + 65 = 125`, and the
seven-cycle's rational partner (`-13801/139, -2579/139, -3799/139,
-5629/139, -4187/139, -6211/139, -9247/139`) sums to `-327 = -17 - 25 -
37 - 55 - 41 - 61 - 91`. The sum of squares is not invariant. So the
Narode mirror gives one exact invariant of the reversal involution — the
trace — and the reversal pairs are pairs of integer cycles of equal trace,
which is a strong constraint on any candidate family.

---

## 6. Typing and what changes for the repo

* No barrier is typed by this note: it produces no termination or density
  statement. Its content is a measurement of the object the fine-scale
  mixing estimate bounds (section 2), two elementary lemmas (Propositions 4,
  5), and one negative structural fact (Proposition 3).
* **For the mixing estimate.** (2.3) ⟹ `M(h) = o(h^(-A))` for every `A`
  (Proposition 2, PROVED). The eighteen FFT levels do not decide between a
  geometric and a polynomial decay of `M(h)`; the closed recursion on the
    powers of two (section 2b) shows `max_j |mu_hat_h(2^j)|` decaying to
  `h = 120` like `C h^(-1.1) r^h` with `r ≈ 0.944` (fit to `0.9%`; a pure
  geometric fit gives ratio `0.930` with `12%` residual), the doubling
  exponents growing linearly with slope `0.083 ≈ |log_2 0.9465|`:
  geometric with a polynomial prefactor, at the no-descent rate `3^(h*-1)
  = 0.9465` asymptotically (OBSERVED). If the maximum stays on the powers
  of two (true wherever checked), `M(h)` decays like `P_h`, which is the
  necessary condition for (2.3) that Proposition 2 extracts (a geometric
  bound on `M(h)` does not by itself give (2.3), whose `ℓ^1` distances are
  dominated by the bulk coefficients), and the conjectured sharp form of
  the coefficient bound is `M(h) ≍ P_h`, the no-descent probability of
  THM-4476/THM-4487: the decay rate of the maximal primitive coefficient
  and the thin-divergence exponent would be one number. The concrete
  inequality to prove is `|mu_hat_h(2^s)| <= C P_h`, i.e. square-root
  cancellation over the descending words.
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
python 04-computation/experiments/collatz_five_mirrors_fourier_deep_20260929.py 18 > 05-knowledge/results/collatz_five_mirrors_fourier_deep18_20260929.out   # 20 GB, 2 min
python 04-computation/experiments/collatz_five_mirrors_reversal_20260929.py     > 05-knowledge/results/collatz_five_mirrors_reversal_20260929.out
python 04-computation/experiments/collatz_five_mirrors_tracesum_20260929.py     > 05-knowledge/results/collatz_five_mirrors_tracesum_20260929.out
python 04-computation/experiments/collatz_five_mirrors_powers_of_two_20260929.py 120 > 05-knowledge/results/collatz_five_mirrors_powers_of_two_20260929.out   # 3 min
python 04-computation/experiments/collatz_five_mirrors_multiplier_families_20260929.py 60 49 > 05-knowledge/results/collatz_five_mirrors_multiplier_families_20260929.out   # 10 min
python 04-computation/experiments/collatz_five_mirrors_powers_of_two_20260929.py 300 > 05-knowledge/results/collatz_five_mirrors_powers_of_two300_20260929.out   # 15 min
python 04-computation/experiments/collatz_five_mirrors_rate300_20260929.py        > 05-knowledge/results/collatz_five_mirrors_rate300_20260929.out
python 04-computation/experiments/collatz_five_mirrors_coherent_20260929.py       > 05-knowledge/results/collatz_five_mirrors_coherent_20260929.out   # 10 min
python 04-computation/experiments/collatz_five_mirrors_coherent120_20260929.py    > 05-knowledge/results/collatz_five_mirrors_coherent120_20260929.out   # 10 min
python 04-computation/experiments/collatz_five_mirrors_excursion_profile_20260929.py 40 80 120 > 05-knowledge/results/collatz_five_mirrors_excursion_profile_20260929.out   # 15 min
```

| item | status |
|---|---|
| Proposition 1 (level independence of the primitive profile) | PROVED; checked at levels 4–18 |
| Proposition 2 ((2.3) ⟹ `M(h) <= C_A (h-1)^(-A)`) | PROVED |
| `M(h)`, `h <= 18`; argmax at `±2^s`, `s - h = 0..5`; ratios `0.65 -> 0.89` | VERIFIED (float64 FFT, matching the closed recursion to six digits); the decay law is not decided by these levels |
| closure of the powers of two under the frequency recursion; `M(h) >= max_j |mu_hat_h(2^j)|` | PROVED |
| `max_j |mu_hat_h(2^j)|` to `h = 120`: `C h^(-1.1) r^h`, `r ≈ 0.944` (pure geometric ratio `0.906 -> 0.930`), doubling exponents `2.1 -> 6.1` with slope `0.083 ≈ |log_2 0.9465|`; `= 0.44–0.48 P_h` over `20..120` | VERIFIED (exact recursion, truncation `2^(-40)`, float64 error `<= 1.1·10^(-10)`, independently recomputed; 30-digit check at `h = 30, 60`); the proportionality OBSERVED |
| the rate of `P_h`: `P_h^(1/h) -> 3^(h*-1) = 0.9465` | PROVED (Chernoff and tilting) |
| the maximum over the seventeen families `u 2^j`, `u <= 49` odd, `3 ∤ u`, is on the pure powers of two at every level `h <= 60` | VERIFIED |
| `M(h) ≍ P_h`, rate `3^(h*-1) = 0.9465` (a necessary condition for (2.3), measured; not the `ℓ^1` estimate itself) | CONJECTURAL; the levels `180..300` lean slightly against it (rate `0.947–0.950`) |
| the recursion to `h = 300`: doubling exponents `3.1 -> 13.0` growing linearly (slope `0.0785`), rate `0.947–0.950` (`3^(h*-1)` at the lower edge), `M/P_h ∈ [0.44, 0.58]` rising after `180` | VERIFIED (certified by the audit's `a <= 60` run); the rate identification OBSERVED in the weak sense, the data leaning slightly against `M ≍ P_h` |
| the resonant coefficient split by excursion above the critical line (`h = 40, 80, 120`): strict no-descent words carry `54%, 41%, 35%`; the band `E < c` reproduces the coefficient to `10%` from `c = 5, 9, 15` at `h = 40, 80, 120` (magnitude within `4%` from `c = 5`); deep descenders cancel slowly, the band width growing like `h/8`; random units at `h = 30` sit at the square-root scale, `10^5` below the resonance | VERIFIED (exact DP and recursion); the mechanism reading DIRECTION; targets (T1), (T2) OPEN |
| Fourier mass per level `0.462 -> 0.472`, typical `|mu_hat|^2 3^h = 0.70` | VERIFIED; the level-by-level identity with the second-moment increments PROVED (Parseval + Proposition 1) |
| cost/last-valuation decomposition of the resonance at `h = 10` | FINITE-EXACT |
| Proposition 3 (Gauss sums: 2-adic reading; no uniform gap) | PROVED; table FINITE-EXACT |
| Proposition 4 (same-length spread lemma) and the collision table | PROVED; FINITE-EXACT |
| Proposition 5 (carry reciprocity); reversed cycles | PROVED; FINITE-EXACT |
| reversal on the integer cycles of `3x+k`, `k <= 41`: an involution with fixed points (rotation-symmetric words) and pairs (`3x+13`: `227 ↔ 259`, `251 ↔ 287`; `3x+37`: `23 ↔ 29`) | FINITE-EXACT (orbits verified; cycles reached from `|x| <= 5·10^4`) |
| Proposition 6 (the cycle-sum polynomial, hence the trace of a rational cycle, is reversal-invariant; the pairs have equal traces `2499`, `125`; the seven-cycle and its partner both sum to `-327`) | PROVED; FINITE-EXACT |
| the five mirrors (section 1) | DIRECTION / remark |

## 8. Audit record (2026-09-29)

**Auditor:** an independent session (own script
`04-computation/experiments/collatz_five_mirrors_20260929_audit.py`, output
`collatz_five_mirrors_20260929_audit.out`, report
`collatz_five_mirrors_20260929_audit.md`, 69 numbered claims): the law
recomputed exactly to level 6 and in float64 to level 14 by an independent
forward recursion; `M(h)` to `10^(-7)` for `h <= 14`; the masses and a
sharpened Parseval identity (level by level); the Gauss-sum table; the
collision table re-run with valuations `<= 70`; the cost decomposition; the
cycles, the `3x+139` orbit and the traces; the (2.3) and §8.1 quotations
against the S19 note; the five paper summaries against the texts (the
Viaclovsky paper against its arXiv page). The note was extended during the
audit (the `h = 14` decomposition, the level-18 profile, the `3x+k`
reversal survey, Proposition 6); those parts received quick checks only,
and section 2b (the powers of two to level 120) was written after the
audit and is not covered by it.

**Verdict: SOUND WITH CORRECTIONS.** Propositions 1, 2, 4, 5, 6 hold;
every recomputed number reproduces; the quotations and mirrors are
accurate. Eighteen corrections, all applied:

1. **The decay reading was overstated.** "Levels off", "flat over
   `14..18`", "would satisfy (2.3) with room" and the conditional "if the
   ratio stays below `1`" were not supported: a shifted power law `C (h +
   2.5)^(-2.5)` fitted the eighteen levels better than a geometric law and
   would have contradicted Proposition 2; the readings separate only near
   `h = 31`. Sections 0, 2 and 6 now say so. (The question was then
   settled for the powers of two by the closed recursion of section 2b,
   which reaches `h = 120`: geometric, at the no-descent rate.)
2. **The argmax law was misdescribed:** `s - h = 0, 1, 2, 3, 4, 5` in steps
   of two or three levels, not "`h + 3` to `h + 5`"; `2^s/3^h` decays
   geometrically; the argmax is defined up to sign.
3. **The collision table measured its valuation cap:** with valuations `<=
   70` the first colliding depth drops and the minimal colliding sum
   decreases with depth; the looseness factor is at least `4.4`, not
   `3–5`; the "depth `3.8 n`" sentence had no derivation and was wrong
   (coverage of the unit classes from `d ≈ n/2`). The table was replaced by
   the audit's.
4. **Two proof slips repaired:** Proposition 3(ii)'s intermediate
   inequality had tail `2^(-s)` instead of `2^(1-s)` (178 violations,
   `|G_9(2)| = 0.207` against `0.499`); Proposition 4's recovery step needs
   `(d, C_w, A)`, not `(d, C_w)` (`C_w` does not involve `a_d`).
5. **The 2-adic source class** is mod `2^A`, not `2^(A+1)`.
6. **The resonance mechanism** was a heuristic stated as fact, its
   "coherence loss `0.96`" a quotient; now DIRECTION.
7. **Mirrors:** the three-quarters "sector restriction" sits on the space
   side, not the frequency side; Merca's positivity has no counterpart
   (`mu_n` is positive by definition) and his "divisor recursion" is a
   linear recurrence.
8. **Slips and labels:** the generic coefficient size is `7·10^(-5)`, not
   `3·10^(-4)`; the Fourier mass per level is `0.462 -> 0.472`, not a
   constant; "undershoots" was inverted; the `C* h^(-6409)` bound is in the
   paragraph after Lemma 8.1; float64 profiles are VERIFIED, not
   FINITE-EXACT; the header said "level 17".

MISTAKES entry: MISTAKE-549. Not checked by the auditor: the literal text
of Tao's Proposition 1.14; levels 15–18 and the `h = 14` decomposition
(the session's FFTs); the completeness of the `3x+k` cycle survey.

**Follow-up audit of section 2b (2026-09-29, same auditor; own numpy
implementation of the closed recursion, brute-force check against the word
sum for `n <= 4`, a rigorous float64 error analysis (`ℓ^∞` contraction:
error `<= n(2^(-40) + 41ε) = 1.1·10^(-10)` at `h = 120`), a 30-digit
recomputation at `h = 30, 60`, an own dynamic programme for `P_h`, three
decay models, the multiplier families to `h = 60`; files
`collatz_five_mirrors_20260929_audit2.py/.out/.md`, 37 claims): SOUND WITH
CORRECTIONS, fifteen applied.** Holds: closure and exactness, every value,
offset, ratio and local exponent, the FFT agreement, `P_h` and the theorem
`P_h^(1/h) -> 3^(h*-1)`, `3^(h*-1) = e^(-I(log_2 3))`. Corrections: (1)
the implication in consequence (i) and in section 6 ran the wrong way —
the family gives a lower bound on the `ℓ^1` distances of (2.3) (necessary,
not sufficient); (2) the decay has a polynomial prefactor (`C h^(-1.1)
r^h`, `r ≈ 0.944`; the linear growth of the doubling exponents with slope
`0.083 ≈ |log_2 0.9465|` is the evidence for the no-descent rate; the
shifted-power residual was a grid artefact; the `P_h` ratios oscillate
and their geometric mean approaches `0.9465` from below); (3) "`0.46 ±
0.02` at every level" was a band the data leave once (`0.481` at `h =
41`): `0.44–0.48`; (4) wording: the `P_h` rate is PROVED, the mechanism is
DIRECTION with phases near `1` involving the suffix sums, "the sup-norm
mixing rate is the no-descent rate" was an overreach; (5) new fact: `s = h
log_2 3 - 6 ± 1` for `20 <= h <= 120`, amending the first audit's
"`2^s/3^h` decreases geometrically". Not checked: the prefactor exponent `3/2` of `P_h`; `M(h) ≍ P_h` beyond the sampled windows.

**Third audit, section 2c (2026-09-29, same auditor; own recursion to `h = 300` with `a <= 40` and a certified `a <= 60` run, reversed-order and noise-injection runs, own `P_h`, own excursion-split DP validated at `c = ∞`, a hundred random units at `h = 30` and the exact coefficient distribution at `h = 12, 14`; files `collatz_five_mirrors_20260929_audit3.py/.out/.md`, 25 claims): SOUND WITH CORRECTIONS, eleven applied.** Holds: every level-300 value (to `4·10^(-7)`), all fits and ratios, the phase factorisation with its inverse-power reading, the split's completeness, the random units. Corrections: (1) the rate bracket had the pure-geometric misfit as its lower end; the admissible fits give `0.947–0.950`, the no-descent rate at the lower edge, and the `M/P_h` rise over `180..300` leans against `M ≍ P_h`; (2) "the band within 6 bits reproduces the magnitude" was a coincidence of a rotated vector at `h = 120`: as a vector the band needs `c = 5, 9, 15` for `10%` (`h/8`), so (T1) must be uniform in `h` with `c = c(h)`; (3) the level-300 certification (`a <= 40` not self-certifying); (4) `0.44–0.48` not `0.46 ± 0.02`, `s - h log_2 3 = -5.3 .. -7.2`, the strict share's trend, rms against median for the random units, the two-sided alignment condition. The excursion-profile paragraph was added after this audit's snapshot and is not covered by it.

**Next probes.** The targets (T1), (T2) of section 2c (cancellation of
deep descents; the band residue), replacing the earlier "coherent
no-descent family with fixed initial phases", which section 2c refutes as
a mechanism (the strict no-descent words carry a third to a half of the coefficient over `h = 40..120`, a decreasing share, and their phases are spread); whether the argmax of `|mu_hat_h|` over all units stays on the
powers of two beyond `h = 18` (a chunked FFT at `h = 19`, or a search over
`±2^s u` for small units `u`); the constant `0.46` as a computable
expectation over the first valuations; the reversal involution on the
rational cycles of `3x + k` beyond the integer ones, and whether the
reversal pairs of `3x+13` and `3x+37` begin an infinite family.
