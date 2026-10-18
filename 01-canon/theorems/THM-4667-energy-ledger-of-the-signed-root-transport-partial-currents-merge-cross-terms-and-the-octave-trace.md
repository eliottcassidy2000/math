---
id: THM-4667
title: "Signed ROOT transport energy ledger: sibling-pair fibres, ROOT projection with its unrooted complement, unconditional current-budget sufficiency, and the primitive affine Fourier recurrence"
status: >
  PROVED (elementary bookkeeping; scope repaired 2026-10-17 by clock_duality_oct15): statements 1-5 and 8-10 (statement 8 checked by brute force on all seeds below
  2 x 10^6; statement 10 checked to machine precision against clock_duality's exact recursion for d <= 4).
  FINITE-EXACT: statement 6 (the ideal energy ladder to d = 6 in exact integers; d = 7..12 in float64 with a rigorous
  truncation bound) and the complete energy profiles of every N < 2^23 in
  statement 7 (exact incremental algorithm, invariant checked against the direct current at every N). NUMERICAL: the
  cross-term anatomies and the operator norms in statement 7. The uniform bound that would close the tail-debt gap is
  HYP-9288 (OPEN; calibrated as beyond Krasikov-Lagarias). The sibling pair and its opposite clock parities are
  clock_duality_oct15's R5. The Syracuse bridge of statement 10 was found concurrently by clock_duality_oct15 in measure
  form (root_positive_transport_20261017g.md, equation (5)). The complete fibre decomposition, the energy (l2) reading,
  the ROOT/shell identity and the additive-Fourier form are new in the repository.
  Repair (clock_duality_oct15): the ROOT projection retains an unrooted energy complement, and the current-budget
  implication survives without universal convergence. Primitive affine phases require the source/terminal boundary
  kernels before comparison with ordinary Syracuse offsets. No novelty claim for the bookkeeping itself.
session: opus-2026-10-17-S41
source: 05-knowledge/results/root_transport_energy_20261017.md
scripts:
  - 04-computation/experiments/root_transport_energy_ledger_20261017.py, root_transport_energy_profile_20261017.py (+ .out, .json)
  - 04-computation/experiments/root_transport_incremental_trace_20261017.py (+ _17/_20/_22 .out and .json; per-N arrays .npz committed for _12/_17 only, the 2^21/2^23 arrays are regenerated in 38 s / 161 s)
  - 04-computation/experiments/root_transport_octave_trace_20261017.py (+ .out, .json)
  - 04-computation/experiments/root_transport_depth_ladder_20261017.py, root_transport_crossterm_anatomy_20261017.py (+ .out, .json)
  - 04-computation/experiments/root_transport_energy_ladder_banded_20261017.py (+ .out, .json; d <= 12 in clock_duality's positive gauge)
  - 04-computation/experiments/root_transport_sibling_anatomy_20261017.py, root_transport_symbol_norm_20261017.py, root_transport_fourier_check_20261017.py (+ .out, .json)
related:
  - clock_duality_oct15 root_braid_logistic_20261017f.md (R1 odd-source fold, R5 sibling fold, sections 3-3c: transport, R7, R8 depth ladder gamma_d), marked_contacts_root_current_20261017e.md (energy kernel G_X, formula (4)), root_positive_transport_20261017g.md (concurrent: P1 positive gauge, equation (5) Syracuse factorisation = measure form of statement 10, P3 uniform l1 remainder for d <= m/16, P4 slow packet)
  - THM-4658 (two-octave energy transport), HYP-9283 (budgets; calibration (B), (C)), HYP-9284 (odd-step suppression), THM-4659/4660, THM-4565 (S20 coupled-orbit covariance)
  - HYP-9288 (energy cone)
---

# THM-4667 — Energy ledger of the signed ROOT transport

## Setting

`T` is the shortcut map. Let `f(n) = (-1)^clock(n)` if n reaches1 (clock = first-hit number of `T`-steps), and `f(n)=0` otherwise. Thus `f(Tn)=-f(n)` globally, including the even ROOT cycle. No universal convergence is assumed. Put
`w(n) = 3[n ≡ 2 (mod 3)] - 1`, `H(M) = Σ_(n<=M) w(n) f(n)`, `R(N) = H(4N) - H(N)`. Following clock_duality_oct15:
- **Sources.** A source is a primitive odd *seed* `s` (`s ≢ 5 (mod 8)`). The map `B(x) = 4x + 1` sends the other odd
  numbers to seeds by `seed(x)`, which preserves `f`.
- **Initial coefficients.** `c_0(s) = Σ_(j: B^j(s) <= 4N) K_N(B^j(s))`, with
  `K_N(m) = w(m)[m <= 4N] + (-3χ(m) - w(m))[m < 2N]` (the two threshold events of R1).
- **Odd transport.** `s -> P(s) = seed((3s + 1)/2^(a(s)))`, `a(s) = ν_2(3s + 1) ∈ {1, 2}`, with sign
  `σ(s) = (-1)^(a(s))`. ROOT `1` is fixed. Write `(Tc)(u) = Σ_(P(s) = u) σ(s) c(s)` and `c_d = T^d c_0`.
- **Energy.** `E_d(N) = Σ_u c_d(u)^2`.

## Statements

1. **(Invariant.)** `Σ_u c_d(u) f(u) = R(N)` for every `d`.
2. **(ROOT-projected partial currents.)** For `R_u^(d)(N) = Σ_(s: P^d(s) = u) c_0(s) f(s)`,
   `R_u^(d)(N)=f(u)c_d(u)`. Consequently
   `E_d(N)=Σ_u R_u^(d)(N)^2+Σ_(u not reaching ROOT) c_d(u)^2`.
   The partial current runs through the depth-`d` class of `u`. Equality with the first sum alone requires that the
   unrooted complement vanish; it holds in the certified finite census, but is not assumed for arbitrary sources.
3. **(Merges.)** `E_(d+1) - E_d = 2 Σ_u Σ_(s < s' in P^(-1)(u)) σ(s)σ(s') c_d(s) c_d(s')`.
   - The fibre is `P^(-1)(u) = {(2^(a_j) B^j(u) - 1)/3 : j >= 0, B^j(u) ≢ 0 (mod 3)}` (clock_duality's formula (10)).
   - Its elements are seeds at scales about `(2^(a_j)/3) 4^j u`.
4. **(End and unconditional sufficiency.)** When every source in the finite support reaches ROOT, eventually
   `E = R(N)^2`, and `R(N)^2 = E_0(N) + 2 Σ_d Cross_d(N)` with the finite merge sum through stabilization.
   Without that assumption, the ROOT coefficient still eventually equals `R(N)`; hence
   `R(N)^2 <= liminf_d E_d(N) <= sup_d E_d(N)`. Therefore `sup_d E_d(N) <= C N` unconditionally gives
   `|R(N)| <= sqrt(C N)`. Convergence of full energy to `R(N)^2` is not supplied without control of the other components.
5. **(Octave trace.)** `T_d(X) := Σ_(X <= N < 2X) E_d(N) = Tr(T^d G_X T^(*d))`.
   - `G_X = Σ_N c_0(N) c_0(N)^T` is clock_duality's energy kernel aggregated on seeds.
   - Always `D_X = Σ_N R(N)^2 <= liminf_d T_d(X)`. If all initially supported sources in the finite octave
     reach ROOT, then eventually `T_d(X)=D_X`. This equality holds for the certified finite census.
   - A bound `sup_d T_d(X) <= C X^2` gives `D_X <= C X^2`. Through THM-4658 this gives a budget of the same type for the
     census energy `E_X`. By HYP-9283(B) that certifies the Mercator pole (constant budgets up to 1806.6, logarithmic
     budgets with `b` of order `10^3`), which forces `limsup log B(X)/log X >= 0.99821`.
6. **(Ideal energy ladder, FINITE-EXACT.)** In clock_duality's ideal recursion R8,
   `K_d(r, b) = Σ_j σ_j K_(d-1)(s_j(r), b - o_j)` with `o_j = a_j - 1 + 2j`.
   - The band-weighted energy `E_d = Σ_b 2^(-b) Σ_r K_d(r, b)^2 / (4·2^d (1 - 2^(-P_d)))` equals the actual `E_d(N)/N`
     to leading order. Its diagonal part is exactly `E_(d-1)`: each old residue is hit once per `(j, branch)`, and
     `(1/2)·(1 + 1/2)·Σ_j 4^(-j) = 1`. So the energy changes only by cross terms.
   - **Values.** For `d = 0, ..., 6`: `E_d = 4.857143, 3.388403, 2.575382, 2.137430, 1.915079, 1.778419, 1.693330` and
     `γ_d = 2.285714, 1.784095, 1.448878, 1.228070, 1.079843, 0.967086, 0.877335`. The first four `γ_d` are as in R8; the
     rest are new.
   - **Cross terms.** `-1.4687, -0.8130, -0.4380, -0.2224, -0.1367, -0.0851`, all negative: the ideal transport
     dissipates energy.
   - **Depths 7-12.** These come from the banded positive-gauge computation (float64, rigorous band tail
     `< 3·10^(-25)`), whose `γ_7..γ_10` agree with clock_duality's certified values to 12 digits.
     - `E_7..E_12 = 1.634126, 1.590908, 1.556295, 1.527640, 1.501754, 1.479358`.
     - Cross terms `-0.0592, -0.0432, -0.0346, -0.0287, -0.0259, -0.0224`.
     - `γ_11 = 0.588440`, `γ_12 = 0.548051`.
   - **Rate.** The consecutive cross-term ratios are 0.55, 0.54, 0.51, 0.62, 0.62, 0.70, 0.73, 0.80, 0.83, 0.90,
     0.87. The decay is slower than geometric, and whether `E_d` has a positive limit is OPEN.
   - **Anatomy (`d <= 5`).** The cross terms concentrate at lag 1, which by statement 8 is the sibling pair.
7. **(Checks, FINITE-EXACT / NUMERICAL.)**
   - *Every `N < 2^23`* (incremental algorithm, 120,409,136 tree nodes; `R(N) = H(4N) - H(N)` verified at every `N`):
     - **Octave traces.** For every octave `X = 2^2, ..., 2^22` the trace `T_d(X)/X^2` attains its maximum at `d = 0`
       (`T_0/X^2 -> 51/7`) and ends at `D_X/X^2`. At `X = 2`, `T_5 = 4 > T_0 = 2`.
     - **Pointwise.** `max_N sup_d E_d(N)/N = 10.06`, at `N = 7910`, where `R^2/N = 9.08`. In the octaves
       `2^18..2^21` every profile peaks at `d = 0`.
     - **Ideal regime.** The actual `E_6/N` averages 1.693 in every octave `2^18..2^22`, matching the ideal value.
     - **Shells.** `Σ_(X <= N < 2X) max_d c_d(1)^2 <= 5.20 X^2` for every octave `2^9..2^22`. This is statement 9's shell
       partial sums with the maximum inside the sum.
   - *Sibling anatomy* (octaves `2^10..2^12` in full, samples at `2^16`):
     - INTRA (statement 8) carries 94-100% of the energy change at `d <= 2`, and the bulk afterwards.
     - The termination climb at `2^12` (trace 1.98 -> 3.98) is INTRA with positive sign.
   - *Operator norm* of the ideal transport on last-digit-balanced vectors, computed as the sup over `|z| = 1` of its
     band symbol: 1.283, 1.324, 1.335, 1.346 for `d = 1..4`. So the ideal transport is not a contraction.
8. **(Sibling factorisation.)** The seeds split into three classes:
   - **heads** `s ≡ 1 (mod 8)`, with `a = 2`;
   - **tails** `2s + 1 ≡ 3 (mod 16)` of the heads, with `a = 1` and `P(2s + 1) = P(s)`;
   - **singletons** `x ≡ 7 (mod 8)` or `x ≡ 11 (mod 16)`, with `a = 1` and `x = (2P(x) - 1)/3`.

   Every fibre `P^(-1)(u)` is a disjoint union of sibling pairs `{s_k, 2s_k + 1}` at scales spaced by the factor 64,
   together with the singleton `(2u - 1)/3` when `u ≡ 2 (mod 3)`. Hence
   - `(Tc)(u) = -[u ≡ 2 (3)] c((2u - 1)/3) + Σ_k (c(s_k) - c(2s_k + 1))`;
   - `E_(d+1) - E_d = INTRA_d + INTER_d`, with `INTRA_d = -2 Σ_(s ≡ 1 (8)) c_d(s) c_d(2s + 1)`;
   - in partial-current coordinates this is `2 Σ_(s ≡ 1 (8)) R_s^(d) R_(2s+1)^(d)
     - 2 Σ_(s ≡ 1 (8), s not reaching ROOT) c_d(s)c_d(2s+1)`; siblings belong to the same basin;
   - INTER_d collects the cross terms between different blocks of one fibre, at scale ratios at least 32.
9. **(ROOT and shells.)** Let `depth(x)` be the number of canonical primitive transports from `x` to ROOT, or infinity if none. This is not by itself a first-hit shortcut clock. Then
   - `c_d(1) = Σ_(x: depth(x) <= d) c_0(x) f(x)`;
   - `E_d(N) = c_d(1)^2 + Σ_(u ≠ 1) c_d(u)^2`.

   The ROOT coefficient is the partial sum of the current over the stopping-time shells `{depth = d'}`, `d' <= d`.
10. **(Fourier form of the ideal transport.)** Put `o = a - 1 + 2j` (so `A=o+1=a+2j` is the homogeneous exponent of the primitive branch, with terminal braid retained), and let `K̂(ξ, b) = Σ_r K(r, b) e(-rξ/3^(d+1))`. Then
    `K̂_d(ξ, b) = Σ_(o=0)^(b) (-1)^(o+1) e(-μ_o ξ/3^(d+1)) K̂_(d-1)(2^(-(o+1)) ξ mod 3^d, b - o)`, with
    `μ_o = 2^(-(o+1)) - (1 - 4^(-j))/3`.
    - Unfolding composes the **primitive** affine maps `P_A(r) = (3r + c_A)/2^A` with `c_A = (3 + 2^a - 2^A)/3`. Their
      translations are not the ordinary Syracuse offsets for the same letter word.
    - After the sign gauge, the branch law is geometric(1/2), and the sign is constant on each slice `Σ o_i = n`.
    - With the band retained, the exact positive-kernel relation is `P = C A`, `S = A C`, `P^d = C S^(d-1) A`:
      - `A` chooses `(3r + 1)/2^a`, `a = 1, 2`, with probabilities 2/3 and 1/3 and band increment `a - 1`;
      - `C` chooses `B^(-j)` with probability `(3/4)4^(-j)` and band increment `2j`;
      - `S` is ordinary geometric Syracuse.

      The boundary kernels cannot be omitted when transferring a Fourier estimate.
    - This bridge was found concurrently by clock_duality_oct15 in measure form (root_positive_transport_20261017g.md,
      equation (5)). Tao's Proposition 1.14 concerns the offset alone, while the band readout keeps the exponent.

## Proof

1. **Invariant.** `f(s) = σ(s) f(P(s))`. One odd shortcut step plus `a(s) - 1` halvings is `a(s)` clock steps, and `seed`
   strips a braid whose shortcut-clock difference is even. More precisely, for actual states `x_i=B^(j_i)s_i`,
   the physical clock through d odd transitions is `Σ_i(a_i+2j_i)+2j_0-2j_d`; retain both boundary braid indices.
   A first-hit controller must separately stop its literal trajectory. The sign relation is exact even off ROOT. So `Σ_u (Tc)(u) f(u) = Σ_s σ(s) c(s) f(P(s)) = Σ_s c(s) f(s)`; the initial
   value is R1. ROOT: `P(1) = 1` and `σ(1) = +1` (the cycle `1 -> 2 -> 1`).
2. **Partial currents.** By induction `c_d(u) = Σ_(P^d(s) = u) σ_d(s) c_0(s)`.
   The invariant says `f(s)=σ_d(s)f(u)`, so `R_u=f(u)c_d(u)`. Now `f(u)^2` is the indicator of the ROOT basin.
   Splitting the coefficient energy into that basin and its complement proves statement2. The inverse formula
   `σ_d(s)=f(s)f(u)` is valid only when `u` reaches ROOT.
3. **Merges.** Expand `(Σ_(s ∈ P^(-1)(u)) σ c)^2`. For the fibre: `(3s + 1)/2^a = B^j(u)` for some `j`, with `a` forced by
   `B^j(u) mod 3`.
4. **End.** On finite initial support, the rooted sources have a finite maximum entry depth. Unrooted sources
   never enter1, so after that maximum `c_d(1)=R(N)` regardless of the other components. This proves the unconditional
   liminf inequality. If all sources are rooted, the entire vector eventually equals `R(N) δ_1`, as in every computed range.
5. **Octave trace.** By linearity `E_d(N) = c_0(N)^T (T^d)^T T^d c_0(N)`. Sum over `N` and use the cyclicity of the
   trace. Applying statement4 to the finitely many rows in an octave gives the unconditional liminf comparison.
   The calibration is THM-4658's descent `sqrt(e_(4X)) <= (sqrt(e_X) + sqrt(δ_X))/2 + sqrt(3/X)` iterated along
   `X_j = 4^j X_0`, followed by HYP-9283(B).
6. **Ideal ladder.** For fixed `j` and branch class, `r -> s_j(r)` is a bijection from the `3^d` residues `r` with
   `B^j(r)` in that class onto `Z/3^d`. Shifting the band by `o_j` multiplies the band weight by `2^(-o_j)`. Summing
   `2^(-o_j) = 4^(-j) 2^(-(a_j - 1))` over `j` and the two classes (factor `1 + 1/2`), together with the density change
   `2^(-1)` between depths, gives the diagonal factor 1. The values come from the exact integer recursion, validated
   against clock_duality's rational `γ_0, ..., γ_3`, and from the anatomy script.
7. **Checks.** Computation (scripts above). The incremental algorithm rests on two facts:
   - `c_0(N + 1) - c_0(N)` is supported on `seed(4N + 1)`, `seed(4N + 3)` and `seed(2N + 1)`;
   - `T` maps a point mass at a seed to a signed point mass at its image.

   So `c_d(N)` for all `d` is updated along three tree paths. Sources at ROOT stay there, and ROOT is kept as a separate
   column per depth.
8. **Sibling factorisation.** Heads and tails:
   - For a head `s ≡ 1 (mod 8)`: `3s + 1 ≡ 4 (mod 8)`, so `a = 2`. Put `t = (3s + 1)/4`. Then `3(2s + 1) + 1 = 2(4t + 1)`,
     so `a(2s + 1) = 1` and `(3(2s+1) + 1)/2 = B(t)`, and `seed(B(t)) = seed(t)`. This is clock_duality's R5.

   The rest of the fibre:
   - **Classes.** A seed `x ≡ 3 (mod 8)` is a tail iff `x ≡ 3 (mod 16)`: for `x = 8k + 3`, `(3x + 1)/2 ≡ 4k + 5 (mod 8)`.
     The remaining seeds have `(3x + 1)/2 ≢ 5 (mod 8)`: it is `≡ 3` or `7` for `x ≡ 7 (mod 8)`, and `≡ 1` for
     `x ≡ 11 (mod 16)`. That is already a seed, so `x = (2P(x) - 1)/3` and `P(x) ≡ 2 (mod 3)`.
   - **Indexing.** In the fibre of `u`, index the elements by `o = a - 1 + 2j`. The element with index `o` exists iff
     `B^j(u) ≡ 2^(-a) (mod 3)`, i.e. iff `u + j ≡ 2^(-a) (mod 3)`, since `B(x) ≡ x + 1 (mod 3)`.
   - **Pairs.** For each `u` the admissible `o` form two residue classes mod 6 (adjacent, or `{0, 5}`). Indices `o`
     (odd, `a = 2`) and `o + 1` give the pair `(s, 2s + 1)`. The step `o -> o + 6` multiplies the scale by `4^3 = 64`.
   - **Singleton.** The only unpaired index is `o = 0`, present iff `u ≡ 2 (mod 3)`, with element `(2u - 1)/3`.
   - **Energy.** Expanding `(Tc)(u)^2` gives INTRA and INTER. On rooted siblings `f(2s+1)=-f(s)` gives the partial-current expression; off ROOT both f-values vanish, so the coefficient complement must be retained.
   - **Check.** Brute force on all seeds below `2·10^6` (pairs for heads below `10^6`; all 375,000 singletons).
9. **ROOT and shells.** `c_d(1) = Σ_(P^d(x) = 1) σ_d(x) c_0(x)`. Here `P^d(x) = 1` iff `depth(x) <= d` (ROOT is fixed with
   sign `+1`), and `σ_d(x) = f(x) f(1) = f(x)` by statement 2.
10. **Fourier form.** Write `v_o(r) = K_(d-1)(s_o(r), b - o)` on its domain. It is the push-forward of `K_(d-1)(·, b - o)`
    under the injective affine map
    - `F_o(s) = B^(-j)((3s + 1) 2^(-a)) = 3·2^(-(o+1)) s + μ_o (mod 3^(d+1))`.

    Then `v̂_o(ξ) = e(-μ_o ξ/3^(d+1)) K̂_(d-1)(2^(-(o+1)) ξ mod 3^d)`. The sign is `(-1)^a = (-1)^(o+1)`. Iterating, the
    phases compose into the primitive affine offset of that word, with band shift `Σ o_i`.
    At `A=3`, `P_3(r)=(3r-1)/8` differs from ordinary `S_3(r)=(3r+1)/8` (values1 and8 at r=0 mod9).
    The probability and band-preserving factorization follows from `U_a B^(-j)=U_(a+2j)` and is proved in
    [positive ROOT transport, section2](../../05-knowledge/results/root_positive_transport_20261017g.md).
    Its terminal clock correction is essential:11 ->17 ->13=B(3) uses three shortcut steps but homogeneous
    letters1,4 of totalfive. Under the sign gauge the branch law has geometric(1/2) letters. Check: statement 10 reproduces R8 to
    `5·10^(-10)` relative `3·10^3` at `d = 4` (`root_transport_fourier_check_20261017.out`).

## Reading

- **The square-root law as energy conservation.** If the relative signs at every merge were fair independent coins, the
  energy would be a nonnegative martingale. Doob gives `Pr(sup_d E_d >= lambda E_0) <= 1/lambda`; this
  is a probability estimate in the model, not a deterministic Collatz bound. A normalized Rudin-Shapiro butterfly
  conserves energy exactly (the unnormalized two-output transform doubles it). By statement 8 the Collatz merge
  keeps only its difference channel, so the energy change of a pair is `((c - c')^2 - (c + c')^2)/2`.
- **Finite ideal dissipation.** Negative aggregate sibling cross terms are observed in the computed ideal regime.
  About two thirds of the initial energy is dissipated by depth6, consistent with HYP-9284 and S37's cancelled octave
  memory. Neither an all-depth sign nor control of INTER follows from those finite values. The operator itself expands some balanced vectors (statement 7), so the dissipation belongs to the actual
  vectors.
- **Open uniform targets.** The sampled large-current rows are ROOT-dominated at late depth. This suggests a
  shell maximal inequality, but it does not identify the complete all-scale gap: uniform ideal cross terms,
  discrete-to-density matching and the nonterminal energy must also be controlled. HYP-9288 implies the shell
  budget; the converse is not proved.

## Correction lineage

The 2026-10-17 clock_duality audit repairs the ROOT-projection equality and the ordinary-Syracuse phase identification
from a430c818f9. All displayed finite census values and the primitive Fourier recurrence are retained.
[Outer packets, section7](../../05-knowledge/results/outer_packet_golden_clocks_20261017h.md) proves the unconditional
survivor and adds the all-depth source-rounding energy bound. The MISTAKES entry records the omitted component and
boundary register; the abstract hostile is not a claimed positive-integer Collatz counterexample.
