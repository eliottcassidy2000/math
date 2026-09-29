# Five papers as structural mirrors of the 3-adic Syracuse law: the Fourier profile of Tao's Syracuse random variable (primitive maxima at the powers of two, decaying at the no-descent rate `3^(h*-1)` on the powers of two — proved to be its upper rate conditional on square-root cancellation at the negative powers of two, S22 — scale-invariant Fourier mass `0.47`), the same-length collision threshold, and carry-polynomial reciprocity under word reversal

**Session:** opus, `collatz-poset-dag-20260927` (S20, S21, S22), 2026-09-29.
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
levels `20..120`) + CONJECTURAL (`M(h) ≍ P_h`: the maximal primitive coefficient decays at the no-descent rate `3^(h*-1)`, a necessary condition for (2.3), not the mixing estimate itself) + DIRECTION. S22 (section 2d): PROVED — Lemma G (a quantitative one-step gap on half the units, with the resonant set named), the exponent-walk identity, Lemma R' (the renewal bound `|c_J| <= mass_h(J) Ñ_(J-1)`), the Chernoff mass law `mass_h(J) <= e^(-hI) e^(θ*δ) (log_2 3 - 1)^(-J)`; CONDITIONAL — Theorem C: if the weighted negative-power norm `Ñ_n` decays at any rate below `log_2 3 - 1 = 0.585`, the resonant window decays at least at the rate `3^(h*-1)` (an upper bound) with a constant prefactor; VERIFIED — `Ñ_n` below the critical rate to level 320 (rate `0.568–0.574`, constant `1.58`), at the Parseval scale for `n <= 122` and again from `136` on but lifted to `8.7` times it at `n = 124..134` by a travelling ridge (the ridges are the multiplier-family coincidences `u · 2^Q ≡ ∓1 mod 3^(n_0)`, OBSERVED, seeds exact, the growth phase of that one OPEN), the bound term by term; OBSERVED — `J_eff ≈ 1.6 √h`; the hypothesis H is CONJECTURAL, and (T2) reduces to one constant `Σ_J γ_J` provided the limits `γ_J` exist (DIRECTION). No Collatz proof step. Audited SOUND WITH CORRECTIONS (eighteen applied; section 8); section 2b, the `3x+k` reversal survey and Proposition 6 were added after the audit's snapshot, and section 2b was then audited separately (fifteen corrections applied; section 8).**

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
`collatz_five_mirrors_tracesum_20260929.py` → `..._tracesum_...out` (the reversal-invariant trace, Proposition 6); S22: `collatz_five_mirrors_gap_20260929.py` → `..._gap_...out` (Lemma G over all units to `n = 9`, the Ramanujan average, the scrambled-level split), `collatz_five_mirrors_renewal_20260929.py` → `..._renewal_...out` (the exponent profile with negative exponents, the negative family, the `J`/`F` decompositions, the Gauss-sum regimes), `collatz_five_mirrors_renewal_bound_20260929.py` → `..._renewal_bound_...out` (Lemma R' term by term, `Ñ_n`, the bound to `h = 300`), `collatz_five_mirrors_renewal_cert_20260929.py` → `..._renewal_cert_...out` (the `a <= 60` certification, `J_eff`, the mass law, the bound to `h = 2000`); `collatz_five_mirrors_negfamily160_20260929.py` → `..._negfamily160_...out` (`Ñ_n` to `160`, the bump diagnostic); `collatz_five_mirrors_ridges_20260929.py` → `..._ridges_...out` (the ridge surface to `200` in a `600`-window); `collatz_five_mirrors_ridge_track_20260929.py` → `..._ridge_track_...out` (the dominant ridge to `320`, Gauss sums and digit runs, `Ñ_n` to `320`); `collatz_five_mirrors_gamma_20260929.py` → `..._gamma_...out` (the `J`-terms at the scale `h^(-3/2) e^(-hI)` to `h = 200`); `collatz_five_mirrors_multiplier_law_20260929.py` → `..._multiplier_law_...out` (the amplitudes and exponents of the families `u 2^j`, `u <= 127`, at `h = 20..50`); `collatz_five_mirrors_floor_angle_20260929.py` → `..._floor_angle_...out` and `collatz_five_mirrors_floor_angle_J_20260929.py` → `..._floor_angle_J_...out` (the top part's real-angle coherence by floor flag, crossing exponent and `J`).

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
0.519, 0.560, 0.579` at `h = 200, 240, 280, 300`): over `20 <= h <= 300` the ratio moves within `0.44–0.58`; the rise by a factor `1.26` over `180..300` is what a rate `0.2%` per level above the no-descent rate produces, and whether `M` and `P_h` have the same exponential order (the conjecture `M ≍ P_h`) or `M` decays slightly slower is not decided by these levels, which lean slightly against the conjecture (S22, section 2d: under the hypothesis H the exponential rate is `3^(h*-1)`, so the lean concerns the prefactor only; without H it stands); "`0.46 P_h`" is a description of the first two hundred levels, not an identity. The resonant exponent is `s = h log_2 3 - 5.3 .. 7.2` at every level `20..300` (`-6.2 .. -7.1` at the multiples of 25; `-6.49` at `h = 300`).

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

### 2d. S22: the renewal structure of the resonant coefficient — a one-step gap lemma, the exponent walk, a renewal bound (Lemma R'), and the sharp rate `3^(h*-1)` conditional on square-root cancellation at the negative powers of two

**Directive.** "keep going with T1 and T2 toward a proof." Worked in the
same worktree; scripts
`04-computation/experiments/collatz_five_mirrors_gap_20260929.py`,
`..._renewal_20260929.py`, `..._renewal_bound_20260929.py`,
`..._renewal_cert_20260929.py`, outputs `05-knowledge/results/collatz_five_mirrors_gap_20260929.out`,
`..._renewal_20260929.out`, `..._renewal_bound_20260929.out`,
`..._renewal_cert_20260929.out`. Constants used throughout:
`m := log_2 3 = 1.58496`, `p := log_3 2 = 0.63093`, `q := 1 - p = 0.36907`,
the tilt `θ* := ln(2q) = -0.30362` (the exponential tilt of the geometric
law `2^(-a)` to mean `m`), `I := θ* m - ln(m - 1) = 0.054979`, and the
identity `e^(-I) = 3^(h*-1) = 0.946505` (PROVED: `ln 3^(h*-1) = -ln p -
(q/p) ln q - ln 3` and `-I = -(1/p) ln(2q) + ln q - ln p` agree term by term,
using `1/p = m` and `m - 1 = q/p`). The number `1/(m - 1) = 1.70951` is
the growth factor of the mass law below.

**Result in one paragraph.** The prefix-excursion coordinates of section 2c
mixed two different scramblings. Reading the phase product *from the top
level down*, the exponent `κ_j = s - T_j` (`T_j = a_j + ... + a_h` the
suffix cost) performs a downward random walk with i.i.d. geometric steps,
and the level-`j` phase `ω_j(κ_j) = e((2^(κ_j) mod 3^j)/3^j)` is near `1`
in the *corridor* `0 <= κ_j < j m - K`, a wrapped power of two below the
floor (`κ_j > j m - K`), and a modular inverse power above the ceiling
(`κ_j < 0`). Three things are proved: (i) a quantitative one-step gap
lemma (Lemma G) that names the frequencies where the geometric Gauss sum
can approach `1` — the corridor is exactly where it does; (ii) an exact
renewal identity `mu_hat_h(2^s) = Σ_J c_J` over the number `J` of bottom
levels above the ceiling, with the bound `|c_J| <= mass_h(J) Ñ_(J-1)`
(Lemma R'), where `mass_h(J)` is the explicit probability of exactly `J`
ceiling levels and `Ñ_n = Σ_(m>=1) 2^(-m) |mu_hat_n(2^(-m) mod 3^n)|` is a
weighted norm of the law's Fourier coefficients on the *negative* powers of
two; (iii) the Chernoff mass law `mass_h(J) <= e^(-hI) e^(θ*δ) (m-1)^(-J)`
at `s = hm - δ`. Consequently (Theorem C, CONDITIONAL): if
`Ñ_n <= C ρ^n` with `ρ < m - 1 = 0.58496`, then `|mu_hat_h(2^s)| <= C' e^(θ*δ) e^(-hI)` for every `h` and every `0 <= s <= hm` — the resonant window
of the powers of two decays at least at the no-descent rate `3^(h*-1)`
with a *constant* prefactor — geometric with an explicit constant on this family and under H, where the cited unconditional bound on all of `M(h)` is polynomial with a tower constant (`C* h^(-6409)`); the comparison is between a conditional bound on one family and an unconditional bound on the full level. The hypothesis is VERIFIED to level `320`: the negative-power coefficients sit at the Parseval scale up to structured ridges (below), `Ñ_n 3^(n/2) ∈ [0.055, 1.56]` for `n <= 80` (`<= 1.29` for `20 <= n <= 80`), rate `0.574` per level over `20..80` and `0.568–0.574` over the windows to `320`, against the critical `0.585`; the margin is `1.3%`, since `3^(-1/2)/(m-1) = 0.98699`, and the ridges are the reason it cannot be called safe. Term by term the inequality `|c_J| <= mass_h(J) Ñ_(J-1)` holds with worst ratio `0.84` (`h = 40`), `0.90` (`h = 80`), `0.91` (`h = 120`), and `0.92` over all `h <= 81` and all `0 <= s <= h log_2 3` (fourth audit), and the bound `Σ_J mass_h(J) Ñ_(J-1)` equals `0.18, 0.14 e^(-hI)` at `h = 40, 100` (computed) and `0.11, 0.08 e^(-hI)` at `h = 300, 2000` under the hypothesis with an extrapolated constant (for `h >= 200` the terms with `J - 1 > 80` are not negligible without H), so the S21 rate bracket
`0.947–0.950` is read as the rate `0.9465` with a slowly varying
prefactor. What (T1) and (T2) become: the deep-ceiling part of (T1) is
Lemma R' plus the hypothesis; the floor part of (T1) is OBSERVED negligible
(the words with any floor level carry `< 1%` of the coefficient at `h = 40,
80`) without a proof; (T2) is reduced to the non-cancellation of a
convergent series of `≈ 1.6 √h` comparable rotating terms `c_J` (`J_eff =
10, 14, 19` at `h = 40, 80, 120`), whose per-mass weights scale roughly
like `1/h` (the ballot factor of the top part staying in the corridor),
which is how `M(h) ≍ h^(-3/2) e^(-hI) ≍ P_h` would arise. No Collatz step;
the reductions are the deliverable.

**Lemma G (one-step gap; PROVED).** Let `t` be a unit mod `3^n`, `0 < t <
3^n`, `θ := t/3^n ∈ (0,1)`, and `ξ := -t 3^(-n) ∈ Z_2` with binary digits
`ξ_1, ξ_2, ...` (`ξ_i = 0` for `i <= v_2(t)`). Then, exactly,
`G_n(t) = Σ_(r>=1) 2^(-r) e(θ/2^r + (ξ mod 2^r)/2^r)` — the residue `t
2^(-r) mod 3^n` equals `(t + m_r 3^n)/2^r` with `m_r = ξ mod 2^r` the unique
representative in `[0, 2^r)` (from `0 <= t + m 3^n < 2^r 3^n`) — and
`|G_n(t)| <= 1/4 + (5/16 + cos(2πφ)/4)^(1/2)`, `φ = (2ξ_2 - ξ_1 - θ)/4`
(the two leading terms in closed form, the tail bounded by `1/4`). Hence:
`|G_n(t)| <= (1 + √5)/4 = 0.8090` whenever the two leading digits of `ξ`
differ, i.e. `v_2(t) = 1` or (`t` odd and `t ≡ (-1)^(n+1) mod 4`); for
`4 | t` the bound is `1/4 + (5/16 + cos(πθ/2)/4)^(1/2)`, below `1` unless
`θ -> 0`; for `t` odd with `t ≡ (-1)^n mod 4` it is `1/4 + (5/16 +
sin(πθ/2)/4)^(1/2)`, below `1` unless `θ -> 1`. The gap classes are half the units up to one (`6562` against `6560` at `n = 9`). Checked against the modular sum to `5·10^(-16)` (`n <= 5`) and over all units to `n = 9`: the two gap
classes never exceed `0.672`, the class `4 | t` with `θ > 1/2` (equivalently `t` odd, `ξ ≡ 3 mod 4`, `θ <= 1/2`) never exceeds `0.846`, and only the classes `4 | t, θ -> 0` and `t odd, θ -> 1`
approach `1` (`0.9967` at `n = 9`), which is Proposition 3(ii) with the resonant set named: `t = ±2^v u` with `v >= 2` and `2^v u` small (an odd resonant `t` is `t = 3^n - 4u'` with `4u'` small; the involution `t -> 3^n - t` pairs `{v_2(t) = 1}` with `{t odd, ξ ≡ 1 mod 4}` and `{4 | t}` with `{t odd, ξ ≡ 3 mod 4}`); at `n = 9` the units with `|G_9| > 0.99` are `±2^8, ±2^9, ±2^10, ±2^11, ±5·2^8`. On the family
`2^k` the phases are `θ = 2^(k - nm)`, so the gap fails precisely in the
corridor `k < nm - K`, and for the resonant `s = hm - δ` the level-`h`
Gauss sum is `1 - O(2^(-δ))`. A corollary used below: for the frequency
`2^(-m) mod 3^n` the gap holds iff the binary digits of the 2-adic number
`-3^(-n)` at positions `m+1` and `m+2` differ.

**The exponent walk (PROVED, an identity).** For a word `w = (a_1..a_h)`,
`Y_h = Σ_j 3^(h-j) 2^(-T_j)` (in `Z[1/2]`), so `e(2^s Y_h/3^h) = Π_(j=1)^h
ω_j(κ_j)` with `κ_j := s - T_j` and `ω_j(k) := e((2^k mod 3^j)/3^j)`
(inverse powers for `k < 0`); this is the phase factorisation of section
2c read from the top: `κ_h = s - a_h`, `κ_(j-1) = κ_j - a_(j-1)`, a
downward walk with i.i.d. geometric steps started at `s`, and no bridge
condition. For `A <= s` the real-number identity `e(2^s Y_h/3^h) =
e(2^(s-A) C_w(3,2)/3^h)` holds with the carry polynomial of section 5
(checked at `h = 6`, all words with letters `<= 4`, to `2·10^(-11)`): the
phase is the real number `2^(s-A) Σ_j 2^(P_(j-1) - jm)` and no 2-adic
scrambling occurs; for `A > s` the factor `2^(s-A)` is a modular inverse.
Ramanujan: `(1/L_n) Σ_(k mod L_n) ω_n(k) = μ(3^n)/L_n = 0` for `n >= 2`
(`-1/2` at `n = 1`; checked to `n = 7`): a walk with uniformly distributed
exponents would give a vanishing coefficient, so the whole coefficient is
a correlation effect between the walk and the phases.

**The corridor and the three regimes (VERIFIED).** At `s = hm - δ` the
exponent walk starts `δ` bits below the critical line (`κ_h - hm = -δ -
a_h`) and moves away from the critical line toward the ceiling `κ = 0` at `2 - m = 0.415` bits per level, so a typical word stays in the corridor from the top until the
exponent turns negative at `T_j > s`, i.e. for the bottom `J ≈ 0.21 h`
levels. Three regimes for the one-step Gauss-sum modulus `|G_n(2^k)|`:
corridor `0 <= k <= nm - 6`: mean `0.81, 0.93, 0.96, 0.97` at `n = 10, 20,
30, 40`; ceiling `k ∈ [-40, -1]`: mean `0.525, 0.545, 0.562, 0.551`,
geometric mean `0.48–0.53`; floor `k ∈ [nm + 1, nm + 40]`: mean
`0.553, 0.527, 0.552, 0.515` — the scrambled exponents behave like the
random units of section 4 (mean `0.543`). The profile `|mu_hat_n(2^k)|`
over `k ∈ [-60, nm + 60]` at `n = 20, 40, 60, 80` (renewal output, C1):
the maximum sits at `k = nm - 5.7, -6.4, -6.1, -6.8`; on the *ceiling side*
`k ∈ [-60, -1]` the coefficients sit at the Parseval scale `3^(-n/2)` (rms
`1.3·10^(-10)`–`2.9·10^(-10)` against `3^(-20) = 2.9·10^(-10)` at `n = 40`;
`1.7`–`2.7·10^(-20)` against `8.2·10^(-20)` at `n = 80`) — the negative family is closed under the recursion (exponents only decrease) and every remaining level is scrambled, so away from the *ridges* described below the coefficients are generic; on the *floor side* `k = floor(nm) + d` they are much larger, `≈ 0.1 · 2^(-d) M(n)` (the ratio to `2^(-d) M(n)` is `0.06–0.14` for `-2 <= d <= 12` at `n = 40, 80`; `7.7·10^(-6), 2.5·10^(-6), ..., 4.6·10^(-9)` for `d = 0..8` at `n = 80`; max `4·10^(-12)`, rms `7·10^(-13)` at `d ∈ [21, 60]` against `3^(-40) = 8·10^(-20)`), because a word can leave
the floor zone in one large first step of cost `2^(-d)` and then run the
corridor. This asymmetry is the reason the renewal bound is organised
around the ceiling and not the floor.

**Lemma R' (renewal bound; PROVED).** For every `h >= 1` and every integer `s >= 0`, with `J(w) := #{j : κ_j < 0} = #{j : T_j > s}` (the set is `{1..J}`
since `T_j` decreases in `j`) and `c_J := E[1_(J(w)=J) Π_j ω_j(κ_j)]`, so
that `mu_hat_h(2^s) = Σ_(J=0)^h c_J`:
`|c_J| <= mass_h(J) · Ñ_(J-1)`, `mass_h(J) := P(T_(J+1) <= s < T_J)`,
`Ñ_n := Σ_(m>=1) 2^(-m) |mu_hat_n(2^(-m) mod 3^n)|` (`Ñ_(-1) = Ñ_0 = 1`).
Proof. `J = 0` is `A <= s`, and `|c_0| <= P(A <= s) = mass_h(0)`. For `J >=
1` write `u = (a_(J+1)..a_h)`, `a = a_J`, `v = (a_1..a_(J-1))`; `J(w) = J`
iff `T(u) <= s < T(u) + a`; put `κ := s - T(u) >= 0`. For `j < J`, `κ_j =
(κ - a) - (a_j + ... + a_(J-1))`, the exponent walk of the word `v` at level
`J - 1` started at `κ - a < 0`, so `E_v Π_(j<J) ω_j(κ_j) = mu_hat_(J-1)(2^(κ-a)
mod 3^(J-1))`. Hence `c_J = Σ_u 2^(-T(u)) 1_(T(u)<=s) Π_(j>J) ω_j(κ_j(u)) ·
Σ_(a>κ) 2^(-a) ω_J(κ - a) mu_hat_(J-1)(2^(κ-a))`, and `|c_J| <= Σ_(κ>=0)
P(T_(J+1) = s - κ) Σ_(m>=1) 2^(-κ-m) |mu_hat_(J-1)(2^(-m))| = [Σ_κ P(T_(J+1) =
s - κ) P(a_J > κ)] · Ñ_(J-1) = mass_h(J) Ñ_(J-1)`. ∎ (For `J = h` the top part is empty, `mass_h(h) = 2^(-s)`; for `s < 0` the identity `c_h = mu_hat_h(2^s)` still holds but `mass_h(h)` as defined is `0`, which is why the lemma is stated for `s >= 0`.) Exact evaluation: `mass_h(0) = P(S_h
<= s)`, `mass_h(J) = Σ_(t<=s) P(S_(h-J) = t) 2^(-(s-t))` with `S_N` a sum of
`N` geometric variables, computable by the recursion `P(S_(N+1) = t) =
(P(S_N = t-1) + P(S_(N+1) = t-1))/2`. Checked term by term against the
exact DP (renewal-bound output, b): at `h = 40`, `s = 57`, the ratios
`|c_J|/(mass_h(J) Ñ_(J-1))` are `0.016, 0.029, 0.053, 0.136, 0.273, 0.439,
0.171, 0.512, 0.732, 0.484, 0.468, 0.206, 0.695, ...`, worst `0.840`; at `h = 80`, `s = 120`, worst `0.810` over `J <= 30` (the DP's range; over all `J` the worst is `0.904` at `J = 64`, and `0.913` at `J = 64` for `h = 120`); the inequality was also checked at every `h <= 81` and every `0 <= s <= floor(h log_2 3)`, worst ratio `0.918` (fourth audit); the DP masses agree with the exact mass law to the DP's cost window (`J <= 20` complete to `< 1%`).

**The mass law (Chernoff form PROVED; local form VERIFIED).** For `0 <= J
<= h` and `s = hm - δ`: `mass_h(J) <= P(T_(J+1) <= s) <= e^(-θ* s) Z(θ*)^(h-J)
= e^(-hI) e^(θ* δ) (m - 1)^(-J)`, with `Z(θ) = E e^(θ a) = (e^θ/2)/(1 -
e^θ/2)`, `Z(θ*) = q/(1-q) = m - 1`, the tilted mean `1/(1 - e^(θ*)/2) = m`
(Chernoff at the fixed tilt `θ* < 0`; `e^(-θ* hm) Z^h = e^(-hI)` by the
definition of `I`). The exact masses normalised by `e^(-hI) e^(θ* δ)` grow
by the factor `1.35, 2.05, 1.90, 1.76, ...` per unit of `J` at `h = 40`,
`1.72, 1.71, 1.71, 1.71, 1.70, ...` at `h = 1000` (renewal-cert output,
iii): the growth factor `1/(m-1) = 1.7095` is exact in the limit, and at
finite `h` it is damped by the Gaussian window of the tilted local limit
theorem, centred at `J ≈ δ/m ≈ 4` with width `≈ 0.61 √h` in `J`, which is
where the `√h` of the next paragraph comes from. Since `1/(m-1) ·
3^(-1/2) = 0.98699`, the series `Σ_J mass_h(J) Ñ_(J-1)` is *marginally*
convergent when the negative family sits exactly at the Parseval scale:
this is the sharp reason the rate question was delicate at `h <= 300`.

**Theorem C (CONDITIONAL; PROVED given the hypothesis).** Hypothesis H:
`Ñ_n <= C ρ^n` for all `n >= 0` with some `ρ < m - 1 = 0.58496`. Then for all `h >= 1` and all integers `0 <= s = hm - δ`, `0 <= δ <= hm`:
`|mu_hat_h(2^s)| <= Σ_J mass_h(J) Ñ_(J-1) <= e^(-hI) e^(θ* δ) [1 + C (m-1)^(-1)
/ (1 - ρ/(m-1))] =: C' e^(θ* δ) (3^(h*-1))^h`.
With `ρ = 3^(-1/2)` the bracket is `1 + 131.4 C`. Two consequences: the
exponential rate of the resonant window is *exactly* `3^(h*-1)` (an upper
bound at this rate here; the lower bound is (T2)), and the low side of the
resonance decays at least like `e^(θ* δ) = 0.738^δ` (observed `≈ 0.85^δ`
near the peak at `h = 80`, steepening further out: `2.5·10^(-13)` rms on
`k ∈ [0, nm/2]` against the peak `6.5·10^(-5)`). The hypothesis is what the
data say (VERIFIED, not proved): `Ñ_n 3^(n/2) = 1.00, 0.93, 0.92, 0.95,
0.90, 0.82, 0.56, 0.56, 0.45, 0.28, 0.48, 0.67` for `n = 1..12`, then in `[0.055, 1.56]` to `n = 80` (max `1.56` at `n = 16`, min `0.055` at `n = 78`; `<= 1.29` for `20 <= n <= 80`);
least-squares rate `0.5736` per level over `20..80` (`0.557` over `40..80`), below `3^(-1/2) = 0.5774` and below the critical `0.585`; the sup over `m <= 60` of the single coefficients has `N_n 3^(n/2) ∈ [0.67, 3.77]` for `n <= 80`. Extended to `n = 320` (`_negfamily160_`, `_ridge_track_` outputs): the rate of `Ñ_n` is `0.5678` over `20..320`, `0.5710` over `160..320`, `0.5744` over `200..320`, and `(log_2 3 - 1)^(-n) Ñ_n <= 1.58` over `20..320` (the maximum at `n = 130`), so H holds numerically with `C = 1.6`, `ρ = 0.585`, to level `320`; but `Ñ_n 3^(n/2)` is *not* uniformly of order one: it drops to `0.002–0.1` on most levels past `140` and spikes to `7.97` at `n = 128` (`N_n 3^(n/2) = 13.4` there, `3.6–13` over `n = 88..130`). The spikes are structured, not noise — the ridges of the next paragraph.
The negative family was certified against the valuation truncation by an `a
<= 60` rerun (`Ñ_n` to `1.3·10^(-11)` relative, single coefficients to
`1.2·10^(-10)`, `n <= 80`), and the tail `m > 60` of `Ñ_n` is at most
`2^(-60)`, which changes the bound at `h <= 300` by less than `10^(-18)`.
For `h <= 150` every `Ñ_(J-1)` that matters in the bound `B_h(s) := Σ_J mass_h(J) Ñ_(J-1)` is a computed number (the terms with `J - 1 > 80` carry at most `Σ_(J>=82) mass_h(J) <= 6.7·10^(-7) e^(-hI)` at `h <= 150` even with `Ñ <= 1`; the tail `m > 60` of each `Ñ_n` is bounded by `2^(-60)`), so the computed inequalities `|mu_hat_h(2^(s*))| <= B_h(s*) <= 0.213 e^(-hI)` (`40 <= h <= 81`; `<= 0.164 e^(-hI)` for `82 <= h <= 150` with `Ñ_n` to `n = 120`; `0.25` at `h = 20`, above `e^(-hI)` for `h <= 19`) and `max_(0<=s<=hm) |mu_hat_h(2^s)| <= max_s B_h(s) <= 0.872 e^(-hI)` (`40 <= h <= 81`; the maximum of `B_h(s)` sits at `s = floor(hm)`, where `e^(θ*δ) ≈ 1`) rest on Lemma R', the exact masses and float64 arithmetic only (fourth audit, output H). The bound against the measured maxima
(renewal-bound output, d): `M(h)/bound = 0.070, 0.043, 0.039, 0.028, 0.021,
0.019, 0.016, 0.010, 0.010` and `bound/e^(-hI) = 0.180, 0.178, 0.136, 0.141,
0.148, 0.120, 0.105, 0.124, 0.110` at `h = 40, 60, 80, 100, 120, 150, 200,
250, 300`; `M(h)/e^(-hI) = 0.0125 -> 0.0011` over the same levels, a
prefactor `≈ h^(-1.2)` that the sup-based bound cannot see.

**The ridges of the negative family (OBSERVED; mechanism identified, VERIFIED at the seeds).** In a window of `450–600` negative exponents followed to level `320` (`collatz_five_mirrors_ridges_20260929.py`, `_ridge_track_`), the surface `v_n(m) = 3^(n/2) |mu_hat_n(2^(-m))|` is a generic background of order one crossed by *ridges*: lines `m = m_0 - 1.3 (n - n_0)` along which `v` is `2–30`. The dominant one is born at `n = 11..15` at `m = 414..409` with `v = 3.1, 4.4, 6.1, 8.4, 12.2`, peaks at `29` (`n = 29`), is `8–25` to `n = 96`, `4–6` at `n = 115..125`, `1–2` at `n = 130..160`, `0.3–0.9` at `n = 165..300`, and reaches the small exponents (`m* = 43, 27, 18` at `n = 300, 310, 320`) with `v = 0.05–0.25`; its slope is `-1.33, -1.26, -1.16, -1.42, -1.24, -1.30` per level over the successive fifty-level stretches. Mechanism: along its birth the one-step Gauss sum is `|G_n(2^(-m*))| = 0.981, 0.992, 0.996, 0.994, 0.997` in the class (`t` odd, `ξ ≡ 3 mod 4`, `θ = 0.84..0.94`) and the 2-adic expansion of `3^(-n)` has a run of `11, 13, 14, 16, 17` equal digits starting at positions `413, 411, 410, 408, 407` (upper end fixed at `423`; the control position `m* + 40` has runs of `2–10` and `|G| = 0.37–0.95`). By the duality behind Lemma G this is the statement that the residue `2^(-423) mod 3^n` is the *small negative integer* `-55` for every `n <= 15`: `v_3(55 · 2^423 + 1) = 15` (`2^(-423) ≡ 3^n - 55 mod 3^n` for `n = 9..15`, `≡ 3^15 - 55 mod 3^16` at `n = 16`). So at the exponents `-423 + j` the negative family *is* the multiplier family `-55 · 2^j` of section 2b at the levels `n <= 15`, and the ridge is that family's resonance: the corridor for a multiplier `u` is narrower by `log_2 u` bits, so the shift `δ -> δ + log_2 u` in the mass law predicts an amplitude `≈ M(n) u^(θ*/ln 2) = M(n) u^(-0.438)` and a resonant exponent shifted by `-log_2 u`. Tested on the closed families `u 2^j` for `u = 5..127` at `h = 20, 30, 40, 50` (`collatz_five_mirrors_multiplier_law_20260929.py`): the exponent shift holds within `±1.5`; the amplitude is only a trend, least-squares exponent `-0.52, -0.55, -0.51, -0.47` with prefactor `0.9, 0.9, 0.7, 0.6`, and a scatter of a factor four that depends on the multiplier's fine structure (`u = 11, 13, 37, 47, 49` give `0.20–0.33` of the pure maximum, `u = 41, 43, 55, 65` give `0.05–0.08`; the ridge's `0.20 M(15)` at `u = 55` is `0.06–0.14` at the higher levels). Once the coincidence ends (the 3-adic digit of `2^(-423) + 55` at position `15` is non-zero) the ridge is a remnant carried along the walk: it decays at `0.5675` per level over `n = 34..300` (`v` from `25` to `0.3`), slightly *faster* than the Parseval rate. The other ridges are the same phenomenon with `u = 1` (the wrapped copies of the pure resonance: `2^(-480) ≡ 2^6 mod 3^6`, `2^(-154) ≡ 2^8 mod 3^5`, from the periods `L_6 = 486`, `L_5 = 162`; the latter continues as `2^(-154) ≡ 13 mod 3^n` for `n <= 7`, `v_3(13 · 2^154 - 1) = 7`, the family `13 · 2^j`) and with other small `u`. The spike of `Ñ_n` at `n = 124..134` (`Ñ_n 3^(n/2) = 0.90, 3.37, 7.97, 8.67, 2.38, 1.80`, maximum `8.67` at `n = 130`; `[0.029, 0.31]` for `81 <= n <= 122` and `0.03–0.35` at the even levels `136..160`) is the level-`5` ridge arriving at `m <= 3` (its path `m ≈ 155 - 1.3 (n - 5)`: `m* = 63, 48, 24, 11, 3, 1` at `n = 84, 96, 112, 124, 128, 130`); at `n = 128` the leading nine equal digits of `-3^(-128)` (`3^128 ≡ 1 mod 2^9`) make the one-step sum coherent (`|G| = 0.967`), but along the rest of its path the one-step sums are generic (`0.17–0.74`, fourth audit), so its growth from `1.0` to `13.4` times the Parseval scale over `n = 84..128` (a local rate `0.61` against `0.577`) is a multi-level coherence whose mechanism is OPEN — the dominant remnant, by contrast, decays at `0.5675`. The extrapolation constant `1.3` used for the bound at `h = 2000` is exceeded by this wave (up to `8.7`); the bound at `h <= 300` is unaffected (the mass at `J - 1 ∈ [124, 134]` is below `10^(-7)` at `h = 300`), and at `h ≈ 600`, where the mass peak reaches `J ≈ 126`, those terms would enter with a factor up to `7`, still a constant. Reading: a ridge is born whenever `u · 2^Q ≡ ∓1 mod 3^(n_0)` for a small `u`, lives at the family's resonant amplitude, `0.05–0.35` of `M(n)` for `u <= 127` (growing against the Parseval scale by `1.64` per level), for `n <= n_0`, then decays as a remnant; since `v_3(u 2^Q ∓ 1) >= n_0` has frequency about `3^(-n_0)` over the pairs `(u, Q)`, the largest ridge that a window of `W` exponents and multipliers `u <= U` can show has `n_0 ≈ log_3(W U)` and amplitude `≈ (W U)^0.45` times the Parseval scale — polynomial, not exponential, which is why H can survive them: the prediction is `Ñ_n <= C n^(1/2) 3^(-n/2)` or so, and the observed maxima (`1.56` at `n = 16`, `3.5` at `n = 32`, `13.4` at `n = 128`, on the single coefficients) are consistent with a slow polynomial growth. Nothing here is proved; the mechanism is exact at the seeds, and the amplitude of a family is a trend with large scatter.

**What the hypothesis is.** `Ñ_n` is the weighted `ℓ^1` norm of the vector
`(mu_hat_n(2^(-m)))_(m>=1)`, which evolves by the closed linear maps
`(L_n v)(m) = Σ_a 2^(-a) ω_n(-m-a) v(m+a)` with `ω_n(-m-a) = e((-3^(-n) mod
2^(m+a))/2^(m+a) + 2^(-m-a) 3^(-n))` — the phases are the reversed binary
digits of the 2-adic number `-3^(-n)` — so H says that the products `L_n
··· L_1 1` contract at a rate below `m - 1` in that norm. At `m = 1` the
coefficient is `E[(-1)^(Y_n) e(Y_n/(2·3^n))]`, the parity of the residue
`Y_n ∈ [0, 3^n)` twisted by a slow phase: H is a square-root-cancellation
statement for one explicit character sum over the 2-adic digits of `3^(-n)`,
of the `×2 ×3` kind, and is CONJECTURAL. Lemma G gives a gap at the step
`(n, m)` iff two consecutive digits of `-3^(-n)` differ, which holds at about
half the steps but does not by itself control the product (the vector `v`
is not constant in `a`). (From memory, not verified against the source:) the Fourier–renewal method would give at best a geometric bound with an unusably small `c` at these frequencies; H needs the sharp constant.

**The `J`-series and (T2) (VERIFIED numerics; the reading DIRECTION).**
The remainders `|Σ_(J'<=J) c_J' - full|/|full|` drop below `10%` at `J_eff =
10, 14, 19` for `h = 40, 80, 120` (`≈ 1.6 √h`; `h/8` and `2.7 ln h .. 4 ln h`
fit worse), which replaces the excursion band width `c = 5, 9, 15` of section
2c: the carrier is the Gaussian `J`-window of the mass law, and its width
grows like `√h` at least until `(1/(1 - 0.987))^2 ≈ 6000` levels. The terms
are comparable and rotating: `|c_J|/|full| = 0.02, 0.05, 0.10, 0.26, 0.52,
0.81, 0.26, 0.59, 0.44, 0.20, 0.10` with arguments `+2.2, +1.8, +3.1, +0.2,
+1.9, -0.9, +2.1, +0.4, +2.3, -1.3, +0.6` for `J = 0..10` at `h = 40`, and
the same pattern at `h = 80, 120` (`0.65, 0.64` at `J = 5, 7`; `0.59, 0.64,
0.62` at `J = 5, 7, 8`). The per-mass weights `|c_J|/mass_h(J)` times `h` are
`1.94, 1.89, 1.48` (`J = 4`), `1.86, 1.97, 1.60` (`J = 5`), `0.62, 0.76, 0.66`
(`J = 7`), `0.35, 0.46, 0.41` (`J = 8`) at `h = 40, 80, 120`: consistent with
a `1/h` law within the spread (exponents between `-1.5` at `J = 2` and
`-0.8` at `J = 9`), not decisive. The reading: the top part must keep its
real angle `Σ_(j>J) 2^(κ_j - jm)` from spreading, which forces the bridge
from the top to stay in the corridor — a ballot event of probability `≍
1/h` — and the words that leave it (`F >= 1` floor levels, `K = 3`) carry `2.4%, 2.0%` of the mass at `h = 40, 80`; their largest single term is `2.3%, 1.3%` of the coefficient and their vector total `0.4%` (the `F = 0` words give `100.4%, 99.65%` of it; the sum of their moduli is `6.5%, 3.4%`; renewal output, C3). So `M(h) ≈ e^(-hI) · h^(-1/2)` (local CLT) `· h^(-1)` (ballot) `· |Σ_J
γ_J|` with `γ_J` the `J`-terms per unit of that scale, i.e. `M ≍ h^(-3/2)
e^(-hI) ≍ P_h` iff the series `Σ_J γ_J` converges to a non-zero limit: that
limit is (T2) in its final form, one complex number. Its convergence is H;
its non-vanishing is open, and nothing here bounds it below. The S21 observation `M/P_h ∈ [0.44, 0.58]` is this constant read at finite `h`, still drifting because the `J`-window is still widening. Direct test (`collatz_five_mirrors_gamma_20260929.py`, `h = 40, 80, 120, 160, 200`): `P_h / (h^(-3/2) e^(-hI)) = 6.8, 8.1, 9.0, 9.2, 9.2` (the no-descent probability has exactly the `h^(-3/2)` prefactor, constant `≈ 9.2`); `|mu_hat_h(2^s*)|/P_h = 0.463, 0.462, 0.443, 0.462, 0.500`; `|mu_hat_h(2^s*)|/(h^(-3/2) e^(-hI)) = 3.16, 3.76, 3.97, 4.26, 4.60`, still growing slowly; the terms `γ_J(h) = c_J/(h^(-3/2) e^(-hI))` have arguments that settle from `h = 120` on (`J = 5`: `-0.86, -0.67, -0.47, -0.50, -0.54`; `J = 12`: `0.82, 1.04, 1.28, 1.28, 1.27`) and moduli that converge only for `J <= 5` (`J = 5`: `2.56, 2.43, 2.34, 2.28, 2.22`; `J = 4`: `1.63 -> 1.19`), while `J >= 7` keep growing (`J = 8`: `1.38, 2.19, 2.47, 2.61, 2.67`; `J = 17`: `0.015, 0.57, 1.80, 3.24, 4.53`) because the Gaussian window needs `h >> J^2`: the constant `Σ_J γ_J` is not computable from `h <= 200`, only its first terms are.

**The floor part of (T1) and the ballot factor, measured (VERIFIED at `h = 40, 80, 120`; added after the audit's snapshot).** Conditional on the crossing exponent `κ_c = κ_(J+1)` the top part (levels `j > J`, real phase `e(Σ_(j>J) 2^(κ_j)/3^j)`) and the bottom part are independent, so the words with a floor level cancel iff their top-part sums `W(κ_c, F = 1) = Σ 2^(-cost) e(angle_top)` are small against their mass. Exact DP with the phase applied at the top levels only (`collatz_five_mirrors_floor_angle_20260929.py`, `_floor_angle_J_`): over all words, the top-part coherence `|W|/mass` is `0.88, 0.90, 0.83` for `F = 0` and `0.077, 0.062, 0.055` for `F >= 1` at `h = 40, 80, 120`, the same at every `κ_c = 0..8` (`0.88 ± 0.01` against `0.07–0.08` at `h = 40`), and `|W(F>=1)|/|W(F=0)| = 0.002, 0.0015, 0.003`: a floor level spreads the real angle by a factor `12–16` in coherence, which is the mechanism, measured. Conditional on `J` as well: the words with `J = 0` (`A <= s`) all carry a floor level, because the corridor `0 <= κ_j < j log_2 3 - 3` is empty at `j = 1` and a single point at `j = 2` (a definitional artefact of `K = 3` at the lowest levels: `2^(κ_j)/3^j` there is a large real angle, not a wrapped power), and their coherence is `0.016, 0.004, 0.002`; for `J = 1..4` the `F = 0` fraction of the `J`-class times `h` is `8.1, 8.6, 7.7` / `19.6, 23.3, 21.8` / `27.7, 36.3, 35.5` / `32.5, 46.0, 47.1` at `h = 40, 80, 120` — the ballot factor `≍ 1/h` of the top bridge staying above the floor, stable in `h` for the small `J` that carry the coefficient (for `J >= 8` the fraction saturates and the factor disappears, which is why the typical words show none) — while the `F = 0` coherence at fixed small `J` decreases slowly (`0.50, 0.34, 0.26` at `J = 4`) and the `F >= 1` coherence is `<= 0.06`, `<= 0.013`, `<= 0.005`. So `c_J ≈ mass_h(J) · (const/h) · (slowly decreasing coherence) · (bottom factor ≈ 3^(-J/2))`, and the measured `c_J/(h^(-3/2) e^(-hI))` at `J = 4` (`1.63, 1.38, 1.28` for `h = 40, 80, 120`) matches the product `mass(F = 0) · coherence` (`0.054, 0.016, 0.008` in units of `e^(-hI)`, i.e. `≈ h^(-1.7)` over this range) to the accuracy of the window. What remains unproved for the floor part is a smoothness statement for the law of the real angle of a bridge that dips below the floor; the numbers above are its content.

**What changed for (T1) and (T2).** (T1) *ceiling part*: for every `J_0`,
the words with at least `J_0` ceiling levels contribute at most `Σ_(J>=J_0)
mass_h(J) Ñ_(J-1)`, PROVED, and under H at most `C' e^(θ* δ) e^(-hI) ·
(ρ/(m-1))^(J_0)` uniformly in `h` — an *absolute* cancellation at the scale
`e^(-hI)`; relative to the coefficient (`≍ h^(-3/2) e^(-hI)` if (T2)) the
uniform version needs `J_0 ≍ (3/2) ln h / |ln 0.987| ≈ 115 ln h`, so the
band width in section 2c's sense is `O(log h)` asymptotically, after the
`√h` regime. (T1) *floor part*: OBSERVED negligible and the mechanism measured (a floor level divides the top-part coherence by `12–16`; the ballot factor `≍ 1/h` of the small-`J` classes verified), no proof; the mechanism is the spread of a real angle, not 2-adic scrambling. (T2): reduced
to `Σ_J γ_J ≠ 0`. The S21 rate bracket `0.947–0.950` is read as the rate `0.9465` with a non-power prefactor, and the "levels `180..300` lean against `M ≍ P_h`" of section 2c is, under H, no longer evidence against the *exponential rate* `3^(h*-1)` (Theorem C forbids any rate above `0.9465` under H, so the fitted `0.947–0.950` must be a prefactor effect); it remains evidence about the *prefactor* (the rise of `M/P_h` by `1.26` over `180..300`), which Theorem C does not address, and without H the S21 reading stands.
Obligations: (1) prove H in any form with `ρ < 0.585` — the first genuinely
new target, a single closed family, generic frequencies, the 2-adic digits
of `3^(-n)`; (2) compute the limits `γ_J` (the top-part bridge sums converge
as `h -> ∞` at fixed `J`) and the constant `Σ_J γ_J` to settle (T2)
numerically; (3) done to `n = 320` (rate `0.568–0.574`, constant `1.58`); next, the fine structure behind the multiplier families' amplitudes (a factor four of scatter around the `u^(-0.5)` trend), the ridge inventory as a function of `v_3(u 2^Q ∓ 1)`, and whether the polynomial bound on the ridges can be proved from the mass law of section 2b; (4) the
floor part of (T1) as a smoothness statement for the law of the real angle.

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
python 04-computation/experiments/collatz_five_mirrors_gap_20260929.py            > 05-knowledge/results/collatz_five_mirrors_gap_20260929.out   # 3 min
python 04-computation/experiments/collatz_five_mirrors_renewal_20260929.py        > 05-knowledge/results/collatz_five_mirrors_renewal_20260929.out   # 1 min
python 04-computation/experiments/collatz_five_mirrors_renewal_bound_20260929.py  > 05-knowledge/results/collatz_five_mirrors_renewal_bound_20260929.out   # 1 min
python 04-computation/experiments/collatz_five_mirrors_renewal_cert_20260929.py   > 05-knowledge/results/collatz_five_mirrors_renewal_cert_20260929.out   # 1 min
python 04-computation/experiments/collatz_five_mirrors_negfamily160_20260929.py 160 > 05-knowledge/results/collatz_five_mirrors_negfamily160_20260929.out   # 1 min
python 04-computation/experiments/collatz_five_mirrors_ridges_20260929.py 200 600  > 05-knowledge/results/collatz_five_mirrors_ridges_20260929.out   # 1 min
python 04-computation/experiments/collatz_five_mirrors_ridge_track_20260929.py 320 450 > 05-knowledge/results/collatz_five_mirrors_ridge_track_20260929.out   # 1 min
python 04-computation/experiments/collatz_five_mirrors_gamma_20260929.py 40 80 120 160 200 > 05-knowledge/results/collatz_five_mirrors_gamma_20260929.out   # 4 min
python 04-computation/experiments/collatz_five_mirrors_multiplier_law_20260929.py > 05-knowledge/results/collatz_five_mirrors_multiplier_law_20260929.out   # 2 min
python 04-computation/experiments/collatz_five_mirrors_floor_angle_20260929.py 40 80 120 > 05-knowledge/results/collatz_five_mirrors_floor_angle_20260929.out   # 2 min
python 04-computation/experiments/collatz_five_mirrors_floor_angle_J_20260929.py 40 80 120 > 05-knowledge/results/collatz_five_mirrors_floor_angle_J_20260929.out   # 2 min
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
| `M(h) ≍ P_h`, rate `3^(h*-1) = 0.9465` (a necessary condition for (2.3), measured; not the `ℓ^1` estimate itself) | CONJECTURAL; the levels `180..300` lean slightly against it (rate `0.947–0.950`); under H (section 2d) the exponential rate is `3^(h*-1)` and the lean is a prefactor statement |
| the recursion to `h = 300`: doubling exponents `3.1 -> 13.0` growing linearly (slope `0.0785`), rate `0.947–0.950` (`3^(h*-1)` at the lower edge), `M/P_h ∈ [0.44, 0.58]` rising after `180` | VERIFIED (certified by the audit's `a <= 60` run); the rate identification OBSERVED in the weak sense, the data leaning slightly against `M ≍ P_h` (under H the lean is about the prefactor only; section 2d) |
| the resonant coefficient split by excursion above the critical line (`h = 40, 80, 120`): strict no-descent words carry `54%, 41%, 35%`; the band `E < c` reproduces the coefficient to `10%` from `c = 5, 9, 15` at `h = 40, 80, 120` (magnitude within `4%` from `c = 5`); deep descenders cancel slowly, the band width growing like `h/8`; random units at `h = 30` sit at the square-root scale, `10^5` below the resonance | VERIFIED (exact DP and recursion); the mechanism reading DIRECTION; targets (T1), (T2) OPEN (superseded in part by section 2d) |
| Lemma G (one-step gap `<= 0.809` on half the units: `v_2(t) = 1` or `t` odd with `t ≡ (-1)^(n+1) mod 4`; the resonant set `t = ±2^v u` small named) | PROVED; checked over all units to `n = 9` |
| the exponent-walk identity; the real-phase identity for `A <= s`; the Ramanujan average `0` | PROVED (checked at `h = 6`, `n <= 7`) |
| Lemma R' (`|c_J| <= mass_h(J) Ñ_(J-1)`) and the Chernoff mass law `mass_h(J) <= e^(-hI) e^(θ*δ) (log_2 3 - 1)^(-J)`, `e^(-I) = 3^(h*-1)` | PROVED; term by term VERIFIED at `h = 40, 80` (worst ratio `0.84`) |
| Theorem C: H (`Ñ_n <= C ρ^n`, `ρ < 0.585`) ⟹ `|mu_hat_h(2^s)| <= C' e^(θ*δ) (3^(h*-1))^h` for all `h`, `s <= h log_2 3` | PROVED given H; H CONJECTURAL |
| the negative family: `Ñ_n 3^(n/2) ∈ [0.12, 1.56]` for `n <= 80`, rate `0.574` over `20..80`, certified by `a <= 60`; to `n = 320`: rate `0.568–0.574`, `(log_2 3 - 1)^(-n) Ñ_n <= 1.58`; at the Parseval scale to `n = 122` and from `136` on, not between (`8.67` times it at `n = 130`, audit 4) | VERIFIED |
| the ridges of the negative family: lines `m = m_0 - 1.3 (n - n_0)` of amplitude `2–30` times the Parseval scale; the dominant one is the multiplier family `-55 · 2^j` (`v_3(55 · 2^423 + 1) = 15`), the others wrapped copies of the pure resonance (`u = 1`) and `13 · 2^j` (`v_3(13 · 2^154 - 1) = 7`); the families' amplitudes `0.05–0.35 M(h)` with a `u^(-0.5)` trend and a scatter of four, their exponents shifted by `-log_2 u`; remnants decay at `0.5675` per level; polynomial size heuristic | OBSERVED; the seeds VERIFIED (exact residues); the size heuristic CONJECTURAL |
| the floor mechanism: top-part coherence `0.88–0.90` without a floor level, `0.055–0.077` with one (`h = 40..120`); the `F = 0` fraction of the classes `J = 1..4` is `(8, 22, 35, 47)/h` (the ballot factor) | VERIFIED; no proof |
| `P_h = (9.2 ± 0.1) h^(-3/2) e^(-hI)` at `h = 120..200`; `|mu_hat_h(2^s*)|/P_h = 0.44–0.50` to `h = 200`; the `J`-terms' phases settle from `h = 120`, their moduli converge only for `J <= 5` | VERIFIED; the constant of (T2) not computable at these levels |
| the bound `Σ_J mass_h(J) Ñ_(J-1) = 0.14–0.18 e^(-hI)` at `h = 40..150` (computed; `0.11–0.12 e^(-hI)` at `h = 200..300` under H), `M/bound = 0.07 -> 0.01`; `F = 0` words carry `99.6–100.4%`; per-mass weights `≈ 1/h` | VERIFIED; `J_eff = 10, 14, 19 ≈ 1.6 √h` OBSERVED (four points with `h = 20`); the readings DIRECTION |
| the profile asymmetry: ceiling side `3^(-n/2)`, floor side `≈ 0.1·2^(-d) M(n)` | OBSERVED |
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

**Fourth audit, section 2d (2026-09-29, a new auditor; own code: the law itself to `h = 10` by forward DP against the closed recursion, the closed recursion on the negative family to `n = 160` with `a <= 40/50/60` and a 30-digit `mpmath` check, a complete walk DP for the `c_J` to `h = 120`, the exact masses, the bound `B_h(s)` at every `h <= 150` and every `0 <= s <= floor(h log_2 3)`; files `collatz_five_mirrors_20260929_audit4.py/.out/.md`, 20 claims): SOUND WITH CORRECTIONS, twenty applied.** Holds: Lemma G (the 2-adic reading as an integer identity, the two-term bound, the digit-change criterion), the exponent-walk and real-phase identities, the Ramanujan average, Lemma R' (for `s >= 0`) and its proof, the Chernoff mass law with `e^(-I) = 3^(h*-1)`, Theorem C's derivation, every recomputed number (`Σ_J c_J` against the closed recursion to `10^(-12)`, against the law to `10^(-16)`; `|m_h(s)| <= B_h(s)` at every level and exponent tested, worst `0.918`). Corrections: (1) three errors in Lemma G's prose — the inequality for the class `4 | t` was reversed (the class with `θ <= 1/2` *contains* the resonance), the odd resonant frequencies are `3^n - 4u'` not `-(odd)`, and the gap classes are half the units up to one; (2) `s >= 0` in Lemma R' and Theorem C; (3) the sentence "for `h <= 81` ... `M_res(h) <= 0.18 e^(-hI)`" was false as stated (undefined quantity, a constant read off three levels, and true only at the argmax exponent) — replaced by the computed statement `B_h(s*) <= 0.213 e^(-hI)` for `40 <= h <= 81`, `max_s B_h(s) <= 0.872 e^(-hI)`; (4) the bound values at `h >= 200` depend on H (the terms with `J - 1 > 80`); (5) the floor-side prefactor is `0.1`, not `1`, and the `d`-indexing was off by one; (6) the floor mass at `h = 80` is `2.0%` not `7.8%` (the difference was the DP's cost window) and "`< 2.3%, < 1.3%`" were the largest single terms, not totals; (7) the ranges of `Ñ_n 3^(n/2)` (`[0.055, 1.56]`) and `N_n 3^(n/2)` (`[0.67, 3.77]`); (8) the status paragraph ("at least at the rate"; `J_eff` OBSERVED; (T2)'s reduction presupposes the limits `γ_J`); (9) the S21 withdrawal holds only under H and only for the exponential rate, the prefactor evidence stands; (10) the comparisons with the Fourier–renewal method were rhetorical or from memory; (11) the auditor's own extension of `Ñ_n` to `n = 160` found the wave at `n = 124..134` (`8.7` times the Parseval scale), independently of this session's ridge work, and correctly notes that its growth phase runs along generic one-step sums (mechanism OPEN). The ridge, gamma, multiplier-law and floor-angle paragraphs were added after this audit's snapshot and are not covered by it.

**Next probes (S21; see section 2d for the S22 state).** The targets (T1), (T2) of section 2c (cancellation of
deep descents; the band residue), replacing the earlier "coherent
no-descent family with fixed initial phases", which section 2c refutes as
a mechanism (the strict no-descent words carry a third to a half of the coefficient over `h = 40..120`, a decreasing share, and their phases are spread); whether the argmax of `|mu_hat_h|` over all units stays on the
powers of two beyond `h = 18` (a chunked FFT at `h = 19`, or a search over
`±2^s u` for small units `u`); the constant `0.46` as a computable
expectation over the first valuations; the reversal involution on the
rational cycles of `3x + k` beyond the integer ones, and whether the
reversal pairs of `3x+13` and `3x+37` begin an infinite family.
