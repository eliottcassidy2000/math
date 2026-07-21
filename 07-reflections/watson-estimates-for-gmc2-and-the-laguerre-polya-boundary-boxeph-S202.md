# Watson estimates for GMC(2): a map of the machinery, and the Laguerre–Pólya boundary

> **POST-PUSH UPGRADE (THM-2023).** The general Laguerre--Polya claim for
> `Phi_(p,q)` is now proved. Gauss multiplication identifies it with a scaled
> positive-parameter `0F_(p+q-1)`, and Baricz--Singh's classical theorem puts
> every zero on the negative real axis. The multiplier-sequence problem proposed
> below is therefore unnecessary for `Phi`. Scope remains important: this is
> the `rd-e=r` boundary only; the opposite `Psi_r=sum y^j/(rj)!` family is not
> covered, and the discrete negative-root parameters still need a lower-order
> transseries certificate outside THM-2017's symmetric monomial model.

*boxeph-2026-07-21-S202. Owner: explore past and possible connections to Watson-estimate ideas.
Continues the GMC(2) finish-map (S201). Object: HYP-8775 (+ correction to HYP-8770). Coordinates with
codex THM-2017/2018 (degree-gap + central resonance).*

## The regime map (which Watson tool covers what)

The GMC(2) residual is the cross-shell descent; via the polar bridge `E[P^m]=L(CT_u[Λ_s^m])`, the
charge-0 projection of a symmetric ±K pair is a **modified Bessel** (already canon, THM-1835):
`Σ_m L(P^m)t^m/m! = L₀(e^{tb} I₀(2t√(ac·XY)))`. Writing the primitive-return channels by their
factorial degree `D(k)=dm+(e−rd)k` (codex THM-2017), the descent splits by the **degree gap
`λ=e−rd`** into three regimes, each with its own Watson tool:

| regime | condition | tool | status |
|---|---|---|---|
| **degree-gap dominant** | `\|λ\|≥r+1` | uniform mixed-factorial lemma + dominated convergence (THM-2017) | **PROVED**: `E[P^m]/L(b^m)→1` |
| **sharp boundary** | `\|λ\|=r` | endpoint resums to a **hyper-Bessel** `Φ_{(p₀,q₀)}(ξ)` (THM-2017) | on the `rd-e=r` side, THM-2023 proves all possible zeros are negative real; opposite side remains distinct |
| **interior / central** | `0≤\|λ\|<r` | **full-entropy saddle** (large-argument generalized-Bessel) | **OPEN beyond THM-2018's proportional slice** (codex HYP-8766) |

The single-shell radial layer is separately closed — the **Radial Lemma** (THM-1565, Watson–Nevanlinna
sector opening `(D+1)π`, Gevrey `D−1`) and the per-component Watson lemma (THM-1665) — and the pure
charge-0 layer by EMP / Cauchy transform (THM-1510/1695). The Liouville endgame (THM-1680) gives an
entire-bounded contradiction but rests on two still-flagged decay lemmas.

## The correction I owed (HYP-8770)

My S201 "symmetric-top dominance" claimed `a_max=C(m,m/2)(αβ)^{m/2}` dominates `E[P^m]=Σ_V a_V V!`.
**That is the crude factorial-gap bound, and it fails** — the inter-level grading is only *polynomial*
`(V_max−j)!/V_max! ~ (m h_top)^{−j}` against *exponential* level coefficients, and the top-term share
measurably collapses (`≈0.67` at m=8 to `≈2·10⁻⁴` by m=24). "Domination is an analytic strategy for an
algebraic fact." The symmetric-top interior is not dominance; it is exactly the **λ=0 full-entropy
saddle** = codex's central resonance (HYP-8766). So my forward contribution moved from the interior
(codex's) to the **boundary**, where a clean classical lens applies.

## The new lens: the boundary functions are Laguerre–Pólya

The boundary limit is `Φ_{(p₀,q₀)}(x) = Σ_{k≥0} x^k/((q₀k)!(p₀k)!)`, a **hyper-Bessel function**, and
NC2 holds on the boundary iff `Φ(ξ)≠0` (`ξ=α/(β^r d^r)`). codex leaves the "discrete zero loci" open.
Classical zero theory settles the shape of that locus:

- **Symmetric base (rigorous):** `Φ_{(1,1)}(x) = Σ x^k/(k!)² = I₀(2√x)`. `I₀` is Laguerre–Pólya: no
  positive real zeros (`I₀≥1`), all zeros real-negative at `x=−(j_{0,k}/2)²` (the `J₀` zeros). Verified
  to 4 digits against `j_{0,k}`.
- **General (now proved by THM-2023):** Gauss multiplication makes
  `Φ_{(p₀,q₀)}` a positive argument-rescaling of a `0F_(p₀+q₀-1)` with all
  denominator parameters positive. Baricz--Singh's theorem therefore puts all
  zeros on the negative real axis. The eight numerical families below are now
  exact coefficient checks, not evidence for the zero theorem.

**Consequence for GMC(2).** Since `Φ_{(p₀,q₀)}∈` L–P (real-negative zeros), on this boundary:
1. **real positive-definite leading data ⟹ NC2 clear unconditionally** (`ξ>0 ⟹ Φ(ξ)>0`, all coeffs
   positive) — matching the trivial real case (THM-1535/1660);
2. **any complex `ξ` off the negative-real axis ⟹ `Φ(ξ)≠0` ⟹ NC2 clear**;
3. the exceptional set is the **explicit** discrete negative-real locus `{ξ : Φ(ξ)=0}`, codimension
   `≥1`. The ODE identifies derivatives there but does **not** universally prove
   the next transseries coefficient nonzero; THM-2017 does so only in its
   symmetric monomial model.

The proposed Polya--Schur detour is unnecessary: THM-2023 resolves the zero
locus directly by the classical positive-parameter `0F_Q` theorem. The
Mittag--Leffler warning remains useful for the *opposite* `Psi_r` boundary,
where `sum x^k/(rk)!` is not Laguerre--Polya for general `r`; the two endpoints
must not be conflated.

## The takeaway

Watson estimates cover GMC(2) shell-by-shell: the Radial Lemma / per-component Watson close the
single-shell and pure-radial layers; THM-2017's dominated-convergence closes the degree-gap-dominant
regime; THM-2023 makes the `Phi` boundary an explicit negative-real zero locus,
with lower-order removal still needed at those roots; and the **nonproportional
interior central resonance** (the entropy saddle, HYP-8766) is the deep piece
where every crude factorial-gap bound fails and a genuine saddle/transseries
argument is still owed. THM-2018 already closes its proportional hypersurface.

Links: HYP-8775, HYP-8770, THM-2017, THM-1835, THM-1565, THM-1665, THM-1680, HYP-8766,
[[the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194]],
00-navigation/GMC2-FINISH-MAP-2026-07-21.md.
