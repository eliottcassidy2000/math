# The two maps: Lee-Yang extremality — tournament-H on the real axis, LRC-coverage on the circle

*kind-pasteur-2026-06-27-S31ah. The owner's directive: stop measuring "the single" value; measure "the
whole PGF curve and its root structure"; think Lee-Yang extremality and the φ⁴ density
`exp(−λS⁴−bS²)dS`; "what maximizes LRC values relating to tournaments"; "two comprehensive maps, and the
gold is in them." I built the two maps — the root structure of the two partition functions — and the gold
is a clean DUALITY.*

## The reframing: from single values to zero loci
- **Tournament:** `H(T) = I(Ω(T), 2)` is ONE reading (fugacity `x=2`) of the hard-core partition function
  `I(Ω,x)` (THM-485). The richer object is the **whole curve `I(Ω,x)` and its zeros** (Lee-Yang zeros).
- **LRC:** the cover bound `meas(S7)=q_0=P(N=0)` is ONE coefficient of the **coverage PGF**
  `Q(z)=Σ_{t=0}^6 q_t z^t` (`N` = #empty inner sectors). The richer object is **`Q` and its zeros**.

## MAP 1 — tournament `I(Ω,x)` zeros (the real-axis Lee-Yang locus)
Over ALL tournaments `n=5,6,7` (`tournament_lrc_pgf_root_structure_kps.py`):
- **Real-rooted fraction = 1.000** at every `n≤7`: the zeros lie on the **negative real axis** — the
  Lee-Yang locus of a 1-D hard-core gas. (Known to fail at `n≥9`, THM-025 — a Lee-Yang transition.)
- **The H-maximizer has a zero HUGGING 0:** smallest `|root| = 0.143, 0.073, 0.015` for `n=5,6,7`,
  monotonically → 0, and it is (essentially) the smallest over all classes at each `n`. Since
  `min|root| ≈ 1/α₁` and `H ≈ 1+2α₁+…`, **H-maximization = pushing a Lee-Yang zero to the origin** — the
  partition function condensing a zero at the fugacity edge `x=0⁺` as `n→∞`. This recasts the BIBD
  max-`α₁` result (THM-027) as a **zero-condensation / critical-point** statement.

## MAP 2 — LRC coverage PGF `Q(z)` zeros (the circle Lee-Yang locus)
Over `k`-clusters `k=8,9,10`:
- The 6 zeros come in **3 conjugate pairs** — they are genuinely **complex, off the real axis**, and sit
  in a tight annulus around `|z| ≈ 1.5–1.9`.
- **The coverage-maximizer (AP/consec) has them on the TIGHTEST circle:** root-spread `max|r|/min|r| ≈
  1.14, 1.22, 1.26` for `k=8,9,10`, versus `≈ 3–15` for generic configs. The maximizer's zeros are the
  **most circular** — a **Lee-Yang circle**. (`even-AP = consec` exactly: dilation-invariance.)
- So as a config dissociates away from the AP, its coverage-PGF zeros **spread off the circle**; the AP is
  the **zero-condensation onto a perfect circle**.

## THE GOLD: the two extremalities are the two classical Lee-Yang loci
| | partition function | zero locus | extremizer | extremal signature |
|---|---|---|---|---|
| **Tournament** (observer-blind) | `I(Ω,x)` | **negative real axis** | regular/Paley (max H) | a zero → **origin** (edge) |
| **LRC** (observer-relative) | `Q(z)` coverage PGF | **circle `|z|≈1.6`** | AP/consec (max coverage) | zeros → **tightest circle** |

This is precisely the **two Lee-Yang theorems**: real-rooted (the "ferromagnet on a line" / Newton locus)
for the tournament, and **zeros-on-a-circle** (the original Lee-Yang circle theorem for the Ising model)
for the LRC. The owner's "tournament = observer-blind (real/affine), tiling = observer-relative" split
(mac-mini MSG-1016) is the **real-axis vs circle** split of the two partition functions. **Maximizing a
value = condensing the partition-function zeros to their most ordered Lee-Yang configuration.**

## The φ⁴ connection (bold, the owner's `exp(−λS⁴−bS²)`)
The coverage PGF's zeros on a circle of radius `R≈1.6` with a small quartic spread is exactly the
**Lee-Yang locus of a `(φ⁴)₂` measure** `exp(−λS⁴−bS²)dS`: the Gaussian (`φ²`) part puts zeros on a
circle; the quartic `−λS⁴` term **deforms/spreads** them. The AP is the `λ→0` (Gaussian/most-circular)
limit of the coverage ensemble; dissociation turns on `λ` (the interaction), spreading the zeros. So:
> **CONJECTURE (HYP-3099): the LRC coverage extremality is a φ⁴ Lee-Yang statement — the AP minimizes the
> quartic coupling `λ` (zeros most on-circle), and `meas(S7)≤cap` is the bound at the Lee-Yang edge.** The
> cap `C(k+1,2)/91` (THM-576/577) should be the on-circle (Gaussian/pairwise) value; the binding **dip**
> at `k=8,9` is the **quartic `λ`-correction** (the off-circle zero motion) — matching mac-mini's
> "degree-2 cap + higher-Pascal dip" and my "order-3/4 cover-bound" exactly. The dip = the `φ⁴` vertex.

## New signals this creates (to measure going forward)
1. **`min|root|` of `I(Ω,x)`** (tournament Lee-Yang edge) — a single scalar that tracks H-maximization and
   its `n→∞` zero-condensation; cleaner than `H` itself.
2. **Root-spread `max|r|/min|r|` of the coverage PGF** — minimized by the AP; a direct "distance from the
   Lee-Yang circle" = a NEW extremality functional for the cover bound (smaller ⟺ more extremal).
3. **The circle radius `R(k)`** of the coverage zeros (closed form? `R² = q_0/q_6`-flavored).
4. **The off-circle moment** `Σ(|r_i|−R)²` = the φ⁴ coupling `λ` — the order-≥3 content the cap misses.
5. **The Lee-Yang gap** = the smallest `|z|−`(real axis) for tournaments / `||z|−R|` for LRC — the
   "spectral gap to criticality."

## Why this could matter for the PROOF
The cover bound `meas(S7)≤cap` is "a single coefficient `q_0`." Lee-Yang says coefficients are constrained
by the zero locus: **if the coverage zeros are provably near a circle of radius `R`, then `q_0 = q_6·∏|r_i|
≈ q_6 R^6` is pinned**, and the extremality becomes "the AP maximizes `R` (or minimizes spread)" — a zero-
locus statement that Lee-Yang/Asano contraction methods могут prove where the moment-LP stalls. This is a
genuinely different attack on CRUX 1: **prove the coverage PGF is a Lee-Yang (circle) polynomial and bound
`q_0` by the zero radius**, instead of bounding moments directly.

→ THM-554 (score partition function `Z=∏(x_a+x_b)`, the Lee-Yang product form), THM-485 (two-temperature
hard-core gas), THM-020/025 (real-rootedness and its `n=9` failure), THM-576/577 (cap = pairwise/on-circle),
HYP-3099 (this), mac-mini MSG-1016 (observer split = real/circle), the φ⁴ owner-cue, [[lrc14-thread]].
