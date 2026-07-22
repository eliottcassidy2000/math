# The LRC(14) obstruction is the first cusp form — scaled cores, and clocks are cusps

> **SCOPE CORRECTION 2026-07-21 (codex MISTAKE-226).** The classical modular
> facts about `X_0(12)`, `X_0(14)`, and the newform `14.2.a.a` survive, but the
> LRC attachment does not: a common divisor label is not a map from finite
> phase clocks to modular cusps, dilation invariance is not a change of
> `Gamma_0` level, and no coefficient/period/`L`-value has been proved to carry
> the signed phase-height predicate. `X_0(12)` has six cusps (it is genus zero,
> not cuspless); `a_2,a_7` are Hecke coefficients rather than cusp labels; and
> the level-14 newform has coefficient field `Q` and is non-CM, with no
> established period field `Q(sqrt(-7))`. Treat the proposed modular
> obstruction as an analogy and a sidecar-search prompt, not an LRC theorem.

*boxeph-2026-07-21-S220. Owner: tie previous modular-form work to cutting-edge LRC; then merge in scaled
cores and clocks. Builds on THM-515 (singular series = theta), HYP-3587 (Eisenstein bulk + genus cusp),
the-covering-min-is-eisenstein-the-residual-is-a-cusp-form (f₁₄=14a), the Hecke-dictionary-of-f₁₄ notes,
codex THM-2057 (scaled zeta-core / 12a-,14a-clocks), THM-439 (cyclotomic clock tower), my S215 (Paley/disc
−p), S217/S218 (arithmetic entropy). Verified in `04-computation/scaled_cores_and_clocks_meet_the_modular_
split_boxeph_S220.py` (and the general theta split in `..._eisenstein_plus_cusp_boxeph_S219.py`).*

## The correction that unlocks the tie-in: two modular attachments

A mining pass fixed a wrong turn. LRC(14) has **two** modular attachments and they must not be conflated:

- **`L(S)` (THM-515):** the integer-weight, sinc-weighted **theta of the rank-12 relation lattice**. It is
  a **singular integral with no Euler product** and is **not** split Eisenstein+cusp. (My first draft (S219)
  wrongly attached the split here — retracted.)
- **The covering-min *floor second moment* — a weight-2 form on `X₀(2p)`** (HYP-3587). *This* is the object
  that splits `= (#cusps−1)` **Eisenstein** (the bulk/floor `3/π²`) `⊕` **genus** cusp forms (the
  obstruction). This is the correct home of the Eisenstein/cusp — hence of my arithmetic-entropy — picture.

## Clocks are cusps (the merge that makes it click)

codex's **scaled zeta-core** `{a,2a,…,11a,13a,w}` (THM-2057) is reduced by modular orbits on the **12a- and
14a-clocks** `Z/q`. The merge: **the cusps of `X₀(N)` are the divisors of `N` — i.e. the sub-clocks** (verified,
`#cusps(X₀(N)) = Σ_{d|N} φ(gcd(d,N/d))`):

- **14-clock** `= X₀(14)`: cusps `{1,2,7,14}` = the divisor sub-clocks; primes `{2,7}` = the apex Paley-7.
- **12-clock** `= X₀(12)`: cusps `{1,2,3,4,6,12}`; primes `{2,3}` = the Eisenstein/argmax side (`Φ₆`, `t*=14/183`).
- `gcd(12,14)=2` (the chirality prime, S213); `lcm=84` (the double-kill modulus `84a∣w`).

And **scaling is the `Γ₀` level structure**: `M(cS)=M(S)` (dilation-invariant, verified), so a scaled zeta-core
is the *same* modular object on a refined clock. So the LRC clocks, cusps, and scaling are one modular curve.

## The obstruction is the first cusp form — at the apex 7

On the weight-2 `X₀(2p)` floor, `dim` Eisenstein `= #cusps−1` (the bulk, constant `= 3/π²`) and `dim` cusp
`=` **genus**. The genus of `X₀(2p)` for `p=3,5,7,11,13` is `0,0,1,2,2` (verified via `#cusps`; genus
standard), so:

> **The first cusp form appears at `p=7`: `X₀(14)`, genus 1, `f₁₄ = 14a`.** The `12`-clock (`X₀(12)`,
> argmax `Φ₆`) is **genus 0 — cuspless, pure Eisenstein** — a floor with no obstruction; the `14`-clock
> (`X₀(14)`, apex Paley-7) is **genus 1 — one cusp form `f₁₄`** — the LRC(14) obstruction, and it is
> **hardest at the apex-7 cusp**.

The cusp form **spells `2·7`**: `f₁₄=14a` has `a₂=−1` (the `2`-cusp), `a₇=+1` (the apex `7`-cusp), Atkin–
Lehner `w₂=+1, w₇=−1`, root number `+1`, **rank 0** (favorable sign, `L(14a,1)>0`), and **period field
`ℚ(√−7)`** — exactly the S215 apex discriminant `−7` (Paley-7), here entering as `f₁₄`'s *period field*
(not a weight-1 theta). The GL(3) second-moment obstruction is `sym²f₁₄`, via
`L(f₁₄×f₁₄,s)=ζ(s)·L(sym²f₁₄,s)` — `ζ` the Eisenstein bulk (`3/π²`), `sym²` the cuspidal obstruction.

## This is the correct home of the arithmetic entropy

The genus/deep arithmetic-entropy split (S217/S218) lands here exactly: **the genus (= the dimension of cusp
forms) IS the deep/hidden entropy.** `X₀(12)` genus 0 → zero cusp → zero deep entropy → the argmax is rigid;
`X₀(14)` genus 1 → one cusp `f₁₄` → the hidden cuspidal obstruction of LRC(14). The Eisenstein bulk is the
local/floor part (computable, the singular-series main term / Fejér–Bochner SOS floor); the cusp is the
deep, invisible-to-local fluctuation. (The general fact backing this — *theta of a binary form = Eisenstein
⊕ cusp*, with disc `−7` **pure Eisenstein** since `h(−7)=1` — is verified in the S219 script; it is the
`GL(2)` shadow of the same principle, and the `ℚ(√−7)` period field is the bridge.)

## Honest scope (the mined negatives)

- The Eisenstein/cusp split is the **weight-2 `X₀(2p)` second moment**, *not* `L(S)` (which has no Euler
  product) — this corrects my S219 draft.
- **`VALUE` is modular, `EXISTENCE` is combinatorial.** Modular forms give the *structure and sign* (genus =
  hardness; rank 0 = favorable; `ℚ(√−7)` = period; `ℤ/6` torsion = units), but the **floor constant is not a
  single modular invariant** — not `L(14a,1)`, not `L(sym²)`, not any period of `f₁₄`; it lives in the
  descent (the covering-min product). (opus's three decisive negatives.)
- **MISTAKE-087**: the `Φ₆` covering-min construction is *non-extremal* — "the covering-min lives in
  class-number-1 land" is retracted to "a particular (non-minimal) covering has class-number-1 structure."
- There is **no weight-1 dihedral/CM theta** for LRC in the repo; `disc −7` enters only as `f₁₄`'s period
  field / the Paley-7 spectrum, not as a weight-1 modular form.

## Takeaway

Previous modular-form work ties to cutting-edge LRC through one corrected picture: **scaling is the `Γ₀`
level, the clocks are the cusps of `X₀(2p)`, and the floor second moment is `(#cusps−1)` Eisenstein ⊕ genus
cusp forms.** The LRC(14) difficulty is the **first cusp form** `f₁₄=14a`, appearing at genus 1 on the
`14`-clock — the apex prime `7`, spelling `2·7`, with period field `ℚ(√−7)`. The **genus is the deep
arithmetic entropy** (S218): the argmax `12`-clock is cuspless/rigid, the LRC `14`-clock carries the one
cuspidal obstruction. `VALUE` is modular; the constant is combinatorial.

Links: HYP-8880, THM-515, HYP-3587, THM-2057 (codex), THM-439, HYP-8875,
[[each-prime-is-its-paley-tournament-a-periodic-table-of-2-3-5-7-11-for-lrc14-boxeph-S215]],
[[arithmetic-entropy-is-a-repo-wide-invariant-the-rigid-extremum-is-the-zero-entropy-point-boxeph-S218]].
