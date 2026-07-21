# Watson estimates for GMC(2): a map of the machinery, and the Laguerre–Pólya boundary

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
| **sharp boundary** | `\|λ\|=r` | endpoint resums to a **hyper-Bessel** `Φ_{(p₀,q₀)}(ξ)` (THM-2017) | clear iff `Φ(ξ)≠0` — **the zero-loci (this note)** |
| **interior / central** | `0≤\|λ\|<r` | **full-entropy saddle** (large-argument generalized-Bessel) | **OPEN** (codex HYP-8766/HYP-8771) |

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
saddle** = codex's central resonance (HYP-8771). So my forward contribution moved from the interior
(codex's) to the **boundary**, where a clean classical lens applies.

## The new lens: the boundary functions are Laguerre–Pólya

The boundary limit is `Φ_{(p₀,q₀)}(x) = Σ_{k≥0} x^k/((q₀k)!(p₀k)!)`, a **hyper-Bessel function**, and
NC2 holds on the boundary iff `Φ(ξ)≠0` (`ξ=α/(β^r d^r)`). codex leaves the "discrete zero loci" open.
Classical zero theory settles the shape of that locus:

- **Symmetric base (rigorous):** `Φ_{(1,1)}(x) = Σ x^k/(k!)² = I₀(2√x)`. `I₀` is Laguerre–Pólya: no
  positive real zeros (`I₀≥1`), all zeros real-negative at `x=−(j_{0,k}/2)²` (the `J₀` zeros). Verified
  to 4 digits against `j_{0,k}`.
- **General (strong evidence):** `Φ_{(p₀,q₀)}` has **all-real, all-negative** smallest zeros and
  **log-concave** (Turán) coefficients for every tested `(p₀,q₀)` — (1,1),(1,2),(1,3),(2,2),(2,3),
  (1,4),(3,4),(3,5). This is the Laguerre–Pólya signature.

**Consequence for GMC(2).** If `Φ_{(p₀,q₀)}∈` L–P (real-negative zeros), then on the boundary:
1. **real positive-definite leading data ⟹ NC2 clear unconditionally** (`ξ>0 ⟹ Φ(ξ)>0`, all coeffs
   positive) — matching the trivial real case (THM-1535/1660);
2. **any complex `ξ` off the negative-real axis ⟹ `Φ(ξ)≠0` ⟹ NC2 clear**;
3. the exceptional set is the **explicit** discrete negative-real locus `{ξ : Φ(ξ)=0}`, codimension
   `≥1`, removed one order down by the hyper-Bessel ODE `θ²Φ=ξΦ` (codex).

So codex's "boundary zero-loci open" becomes a **named classical problem**: prove `Φ_{(p₀,q₀)}∈` L–P.
The right machinery is **Pólya–Schur / multiplier sequences** (is `{1/(ak)!}` a multiplier sequence?
— equivalently `Σ x^k/(k!(ak)!)∈` L–P), anchored by the **Mittag–Leffler reality theorem**
(`Σ x^k/Γ(ak+1)` has only real-negative zeros iff `a≤2`) and **Schur's Hadamard-composition theorem**
(the Hadamard product of two L–P⁺ functions is L–P⁺). Caveat: single reciprocal-factorial functions
`Σ x^k/(ak)!` are *not* L–P for `a≥3` (Mittag–Leffler `a>2`), so the reality of the *product*
`Φ_{(p₀,q₀)}` is genuinely a multiplier-sequence statement, not a Hadamard corollary — and where the
higher zeros of `Φ` might turn complex is precisely codex's **inner resonance band**.

## The takeaway

Watson estimates cover GMC(2) shell-by-shell: the Radial Lemma / per-component Watson close the
single-shell and pure-radial layers; THM-2017's dominated-convergence closes the degree-gap-dominant
regime; the **boundary is a Laguerre–Pólya problem** (this note — evidence + the `I₀` base + the
Pólya–Schur/Mittag–Leffler anchor), reducing codex's zero-loci to an explicit removable locus; and the
**interior central resonance** (the entropy saddle, HYP-8771) is the one deep piece where every crude
factorial-gap bound provably fails and a genuine saddle/transseries argument is still owed.

Links: HYP-8775, HYP-8770, THM-2017, THM-1835, THM-1565, THM-1665, THM-1680, HYP-8766,
[[the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194]],
00-navigation/GMC2-FINISH-MAP-2026-07-21.md.
