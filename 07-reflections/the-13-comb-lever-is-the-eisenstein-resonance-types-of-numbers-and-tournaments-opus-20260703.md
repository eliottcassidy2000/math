# The 13-comb lever is the Eisenstein resonance — and rational/irrational time, tournament type, and covering certificate are one trichotomy

*opus-2026-07-03-S52. The owner asked me to work the 13-spaced comb lever (OPEN-Q-108, the wide-far/slow-base
residual) creatively — Vitali, the rational/irrational and odd/even dualities, and how types of numbers relate
to tournaments. The lever turns out to be a single arithmetic fact — `14` is a primitive 6th root of unity
mod `183 = Φ₆(14)` — and that fact organizes the whole covering ↔ tournament ↔ number-type picture.*

## The lever is `14 = e^{iπ/3}` mod 183 (verified)

The covering-min binds at `t* = 14/183`, and `183 = Φ₆(14) = 14² − 14 + 1 = 3·61`. So **14 is a root of the
6th cyclotomic polynomial mod 183** — a primitive 6th root of unity. Its powers mod 183 are the six
Eisenstein units:
```
    14¹=14,  14²=13,  14³=−1,  14⁴=−14,  14⁵=−13,  14⁶=1        (mod 183)
```
Three consequences, all exact:
- **`13 = 14²` mod 183**, and **`14·13 = 14³ = −1` mod 183** — so `‖13·t*‖ = 1/183`. A **13-spaced comb**
  `{w, w+13, …, w+13(r−1)}` therefore has phases at `t*` forming an AP with step `1/183`, spanning
  `(r−1)/183 ≤ 12/183 = 0.0656 < 1/14`. **The whole comb fits inside one danger-arc width** — it is "as safe
  as a single runner." This is the lever, and it is nothing but the resonance `14·13 ≡ −1`.
- **The deep well `{1,…,12, 182}` (182 = 13·14) is lonely ONLY at `t* = 14/183`**: `M(S) = 14/183`, attained
  at denominator `183`, and the best witness with `q ≤ 45` gives min-distance `1/15 < 1/14` — **not lonely**.
  The small-`q` census provably fails; the certificate is the single Eisenstein rational `14/Φ₆(14)`.
- **The binding pair is the antipode.** At `t*`, runner 1 sits at `14/183` and runner `182 ≡ −1` sits at
  `−14/183 = 169/183` — antipodal. The `Z/2` of the order-6 group (`14³ = −1`) is *literally* the binding pair.

## Types of numbers ↔ tournaments ↔ covering certificates (one trichotomy)

The **Cayley transform** makes the bridge exact: a tournament `T` gives a skew `S = A − Aᵀ`, then
`U = (I+S)⁻¹(I−S) ∈ SO(n)` with eigenvalues **on the unit circle**. The *number-type of the eigenvalue* is
the tournament type — and it is the *same* number-type as the covering witness `t`:

| number type of `e^{2πi t}` | tournament (Cayley eigenvalue) | covering time `t` | certificate |
|---|---|---|---|
| **root of unity** (`t ∈ ℚ`) | **transitive / reducible** (eigenvalues = roots of unity) | **resonant, pinned** | rational — the **census** |
| ↳ **Eisenstein** (`t = a/Φ₆(m)`) | the **6th-root/antipode** sublattice | the deep-well **binding** `14/183` | the **`183` = Φ₆(14) certificate** (this lever) |
| **quadratic irrational** (`i√p`) | **Paley / regular** (skew spectrum `±i√p`, Gauss sum) | resonant at a prime scale | the heptagon `√−7`, `√21` |
| **generic irrational** | **mixing / regular** (equidistributed) | free, equidistributed | the **sweep** (THM-608) |

So the two-sided architecture (bounded-magnitude **rational census** vs large-magnitude **irrational sweep**)
is the *coarse* version of a **finer trichotomy on the type of `t`**: small-`q` rational (census), **deep
cyclotomic rational** (the Eisenstein Φ₆ lever — invisible to bounded-`q`), and irrational (sweep). **The
13-comb lever is the middle mode the two-sided split skipped over.**

## The rational/irrational duality, refined (with kps-S31)

kps-S31 (`base_floor_of_cite`) proved the clean half: *the good set is open, so a real (possibly irrational)
lonely point forces a rational one* — the census direction. This lever supplies the other half's structure:
**which** rational. The deep well shows it is **not** always small-`q`; it can be the **cyclotomic** `Φ₆(14)`.
So "rational" bifurcates:
- **generic small-`q`** (Farey, the census, `q ≤ ~45`) — the *equidistribution-adjacent* rationals;
- **special cyclotomic `Φ_k(n)`** (here `Φ₆(14) = 183`, Eisenstein) — the *resonant* rationals where the
  covering-min actually binds, and where the danger arcs tile the circle into `Φ₆(14)` equal pieces.

The extremal (hardest) configurations live at the **cyclotomic** rationals, exactly the ones a bounded-`q`
search steps over — which is why mac-mini's HYP-4040 proves *no uniform arithmetic band* closes it: the band
would have to reach `q = Φ₆(max-speed) → ∞`.

## Odd/even = the antipode inside the Eisenstein

`order(14 mod 183) = 6 = 2 × 3`, and this factorization *is* the odd/even ⊕ Eisenstein split:
- the **`2`** (`14³ ≡ −1`) is the **antipode / complement** — the `Z/2` of the merged metagraph
  (`R = ` complement), realized as the binding pair `(1, 182)` at `±14/183`;
- the **`3`** (`14² ≡ 13`, `3 | 183`) is the **Eisenstein cube-root** `ω` — the `ι-even` side (`√−3`).

And `14 = 2·7`: the **`2`** is that same antipode, the **`7`** is the **heptagon / Paley** prime (`√−7`, the
`ι-odd` certificate, THM-503). So the covering binds where the **even/Eisenstein `√−3`** (in `183 = 3·61`)
meets the **odd/heptagon `√−7`** (in `14 = 2·7`) — and their cross is `√21 = √(3·7)`, the OPEN-Q-108 residual
of my S27–S30 biquadratic work. **The covering-min lives at the compositum `ℚ(√−3, √−7)`**: Eisenstein-3 in
the denominator `Φ₆(14)`, heptagon-7 in the modulus `14`.

## Vitali / three-gap (the measure geometry)

At the **rational** `t*`, the phases `{v·t* mod 1}` land on the finite lattice `(1/183)ℤ` — the circle is
partitioned into `183 = Φ₆(14)` equal arcs (a *rational Vitali cover* — a genuine partition, no choice
needed). The 13-spaced comb occupies one length-`(r−1)/183` sub-arc; the three-gap theorem for `{k·13t*} =
{k/183}` degenerates to a **single** gap `1/183` (uniform). At an **irrational** `t`, the same phases
equidistribute (Weyl) and the three-gap theorem gives the genuine `≤ 3` distinct gaps — the
Vitali/`ℝ-mod-ℚ` (non-measurable-flavored) regime the free sweep exploits. **The census counts on the finite
Eisenstein lattice; the sweep integrates over the equidistributed circle.** The lever is the statement that
on the finite lattice `(1/Φ₆(14))ℤ`, a covering family's danger arcs still leave the deep-well gap open.

## The mechanism for `r ≥ 2` (the open renormalization piece)

`t* = 14/183` sits at distance `1/(13·183)` from `1/13` — the *obvious* 13-resonance, where the comb
**collapses** to a single point (`(w+13k)/13 ≡ w/13`) that a covering comb (`13 | w`) pins **at the origin**
(dangerous). `t*` shifts off `1/13` by exactly one lattice step to **lift the collapsed comb off the origin**
into the safe band, while the base `{1,…,12}` — `1`-Lipschitz with slack `13/2562` — stays safe. So the
general lever is a **resonant sweep**: sweep `t` in the tiny window around `1/13` where the comb stays tight
(`|13t| < 1/(14(r−1))`), placing the fast center in the safe band. This is a THM-608 variant keyed to the
**phase** spread `(r−1)/183` (small by resonance) rather than the **speed** spread `13(r−1)` (large) — which
is exactly why THM-608's condition (ii) (speed-spread `D`) does not cover it, and a phase-spread version does.

## Status

- **Verified (exact):** `14` is a primitive 6th root mod `183 = Φ₆(14)`; `14·13 ≡ −1` (comb resonance);
  comb span `(r−1)/183 < 1/14`; deep well lonely **only** at `14/183` (census `q ≤ 45` fails at `1/15`);
  binding pair `(1,182)` antipodal; base slack `13/2562`. (`comb13_eisenstein_resonance…S52.py`.)
- **Already formalized (r = 1):** the deep-well `14/183` certificate — kps `LadderPackData.deepWell_lonely`.
- **The insight (this session):** the lever = the Eisenstein resonance; **covering time-type = tournament
  type = number type** (a trichotomy: small-`q` census / cyclotomic Eisenstein / irrational sweep); odd/even =
  antipode `Z/2` inside Eisenstein `Z/6`; the compositum `ℚ(√−3,√−7)` (183 = 3·61 meets 14 = 2·7) with cross
  `√21`.
- **Open (route):** the `r ≥ 2` **phase-spread** THM-608 variant (a resonant sweep in the `1/13` window) — the
  next Lean target for the renormalization's resonant branch.

Related: OPEN-Q-108; THM-608 (HYP-4045) + the tower (HYP-4046); kps `deepWell_lonely` / `base_floor_of_cite`
(HYP-4044 chain, the rational/irrational duality); mac-mini HYP-4040 (no uniform band = `Φ₆(max-speed) → ∞`);
klein HYP-3791 (non-13 combs equidistribute = the irrational side); the biquadratic `ℚ(√−3,√−7)` / `√21`
(S27–S30); the Cayley transform (S30, two-geometries-one-algebra); THM-503 (seven-vanishing). Script:
`04-computation/comb13_eisenstein_resonance_opus_20260703_S52.py`.
