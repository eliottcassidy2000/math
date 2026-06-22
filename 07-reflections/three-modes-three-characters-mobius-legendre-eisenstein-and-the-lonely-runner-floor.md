# Three recursion modes = three characters (Möbius, Legendre, Eisenstein), and the lonely-runner floor is their principal term

*kind-pasteur-2026-06-22-S31q. The owner asked how coprime density / Euler totient / multiplicative
functions relate the three tournament recursion modes `A+B+C−D−E−F+G`, `A+B−C`, `A+B−C+D−E−F+G`. They are
three arithmetic characters, and they organize the LRC(14) floor at the apex prime 7.*

## The three modes ARE three characters
The three recursion modes (Mode A: `n→n−1` hypotenuse; Mode B: `n→n−2` Cayley–Dickson binary; Mode C:
`n→n−3` Eisenstein ternary, THM-291/HYP-2689) are three **generators**, and a recursion over them is an
inclusion–exclusion over the `2³−1 = 7` nonempty subsets of `{A,B,C}`, graded `3+3+1` (THM-549, the
half-tiling recurrence). Decoding the owner's sign words over `A..G = 1..7`
(`lrc_three_mode_character_decomposition_kps.py`):

| word | signs | `+` on | `−` on | character |
|---|---|---|---|---|
| `A+B+C−D−E−F+G` | `+++---+` | singletons `{1,2,3}`, triple `{7}` | pairs `{4,5,6}` | **Möbius μ** (incl–excl) |
| `A+B−C+D−E−F+G` | `++-+--+` | `{1,2,4,7}` = QR(7) | `{3,5,6}` = NQR(7) | **Legendre χ₇** |
| `A+B−C` (+ `S_ω=A+ωB+ω²C`) | `++-` | doubling orbit `{1,2,4}` | — | **Eisenstein χ₃ / ω** |

So: **Möbius (additive incl–excl) ↔ Legendre (quadratic, apex prime 7) ↔ Eisenstein (cubic, `2n−1=27=3³`)**.
These are the multiplicative-function avatars of the three modes.

## The totient is the Möbius mode (the copy rule, HYP-2882)
The bridge to coprime density is the **copy rule** `Σ_{d|n} c(d) = n ⟹ c = μ * id = φ` — *the Euler
totient is the Möbius transform of the identity* (verified `n=1..8`). On a `q`-grid `φ(d)` is exactly the
number of residues of exact denominator `d`. So the Möbius mode is not a metaphor: it is literally the
totient/coprime-density structure (`1/ζ(2) = Σ μ(n)/n²`) that the LRC floor is built from
(`zeta2-governs-the-lonely-runner-floor.md`). And `H(T) = Π H(strong atoms)` is the same Möbius/Euler
product on the tournament side (codex-S100): **`φ`-multiplicativity ↔ `H`-multiplicativity, both primitive
packets before the scalar quotient.**

## The floor decomposed: principal (μ) + apex (χ₇) + sporadic (χ₃)
The apex-7 witness floor is the totient sum `floor = Σ_{b<7} φ(b)·2(7−b)/(7b) = 146/35 ≈ 4.171` (per V),
a sum of **positive** terms. Splitting by the apex character `χ₇(b)`:

> `floor = 4.171` = QR-part `2.857` + NQR-part `1.314`; the **χ₇ oscillation `QR−NQR = 1.543`** is the
> apex-7 *bias* (the doubling orbit `{1,2,4}` is weighted more), and `χ₃` splits the `27=3³` Eisenstein
> part. The principal totient term (`4.171`) **dominates** the χ₇ oscillation (`1.543`): `osc/floor < 1`.

Across the whole composite family LRC(2q) (`lrc_floor_character_qtrend_kps.py`, q = 3,5,7,11,13,17,19,23):
`osc/floor` stays `< 1` for every q (`0.60, 0.23, 0.37, 0.37, 0.07, 0.03, 0.26, 0.20`), **largest at the
smallest `q=3` (LRC(6) is the tightest, matching the repo's q-uniform margin)**. The principal never
loses to the character oscillation.

## What this says for proof and disproof
- **Proof side:** the floor is manifestly positive (a totient/Möbius sum of positive widths), and
  q-uniformly the principal dominates the apex χ₇ oscillation — so the floor never vanishes and never
  flips sign. This re-derives the `3/π²` floor's positivity *with its apex structure exposed*: the hard
  prime-flavored cases are exactly where the totient principal is thin (few small resonance denominators),
  not where the sign flips.
- **Disproof side:** a counterexample would need the χ₇ (apex) or χ₃ (sporadic-27) **character
  oscillation to overwhelm the Möbius principal** and drive `floor < cap`. The computation shows the
  oscillation is bounded by the principal (`osc/floor < 1` always), so the characters can *bias* the
  floor but cannot *flip* it. The disproof has no room in the character oscillation — it is pushed back,
  as before, onto the CAP comparison (the prime-thin margin), i.e. OPEN-Q-108.
- **The unification (`14 = 2·7`, `2·14−1 = 27 = 3³`):** the LRC(14) is governed by exactly the three
  characters of its arithmetic — `μ` (the AP / consec-extremality / coprime principal), `χ₇` (the apex
  prime, `H=7` forbidden, the parity-dual scar), and `χ₃` (the Eisenstein `3³` sporadic = the GW tiler,
  HYP-2138). The three tournament recursion modes and the three LRC characters are the same three objects.

## Net
Not a new proof — it identifies the owner's three modes as **Möbius / Legendre / Eisenstein characters**,
shows (via `φ = μ*id`) that the *coprime-density floor is the Möbius mode*, decomposes the apex-7 floor
into principal + χ₇ + χ₃, and verifies q-uniformly that the principal dominates the character oscillation
(so no disproof lives there). The residual is unchanged — the cap comparison over covering clusters — but
it now wears the right clothes: the difficulty of `n` is the thinness of its Möbius principal relative to
its prime characters, exactly the `ζ(2)`-floor reflection's "arithmetic of n is the difficulty of n,"
made into a character decomposition.

→ HYP-2882 (copy rule / packet duality), THM-291 (Mode B), HYP-2689 (Mode C Eisenstein), THM-549
(3+3+1 half-tiling), `zeta2-governs-the-lonely-runner-floor.md`, HYP-2856 (3/π² floor), HYP-2138
(sporadics iff 2n−1 composite), `the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md`,
[[lrc14-thread]].
