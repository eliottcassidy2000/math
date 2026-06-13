---
source: monad-explorer-2026-06-06-S700
status: clean elementary theorem (AP tight times = units mod n, verified all n in 5..23)
        + structural reframe (even n = two negation fixed points; the second apex n/2)
        extending S679 from prime/odd n to even n.  LRC(14) itself NOT proved.
tags: [LRC, n14, n19, negation-involution, fixed-point-residual, euler-phi, units,
       parity, apex, second-apex, Vstar, sieve, reframe, exploration]
---

# The second apex: even-n LRC and the two fixed points of negation

S679 (`lrc19-sieve-apex-involution`) gave a beautiful reframe of the Lonely Runner
Conjecture at the **prime** frontier n=19: the hard core is the *fixed-point
residual of the negation involution* `σ: a ↦ −a` on `Z/n`, with the apex (a
sieve-forced multiple of n, residue 0) sitting at the unique fixed point. The
table it drew — LRC@prime-n, Rédei reversal, V* sign gauge, κ-even — unified four
threads around "the hard part is the fixed-point set of an involution."

That reframe is **parity-special**, and this is the gap I chased. Negation on
`Z/n` has fixed-point set `{a : 2a ≡ 0 (mod n)}`:

- **n odd:**  `Fix(σ) = {0}` — one apex.
- **n even:** `Fix(σ) = {0, n/2}` — **two** apices.

n=14 is the smallest even unproved case, and `Fix(σ₁₄) = {0, 7}`. The second
fixed point `n/2 = 7` is exactly the structural feature S679's prime picture
cannot see, and (I claim) it is the same object the codex sessions kept circling
as "parity" / the "V* seam" (HYP-2218 / HYP-2230 / HYP-2231) without naming it as
a negation fixed point.

## The clean theorem underneath (verified, elementary, all n)

For the canonical near-counterexample `AP_n = {1,…,n−1}` (gap `M = 1/n`):

> **The simultaneous tight times on the division grid are exactly the units:**
> `j/n` is tight (every runner `≥ 1/n`) **iff `gcd(j,n)=1`.**
> So the tight-time set is `{ j/n : j ∈ (Z/n)^× }`, of size **`φ(n)`**.

*Proof.* `‖a·j/n‖ = dist(aj mod n)/n ≥ 1/n` iff `aj ≢ 0 (mod n)`. If `d=gcd(j,n)>1`
then `a = n/d ∈ {1,…,n−1}` gives `aj ≡ 0` (that runner sits on the origin); if
`gcd(j,n)=1` then `aj ≢ 0` for every `a ∈ {1,…,n−1}`. ∎

Verified exactly (Python `Fraction`) for n = 5,7,11,13,14,19,21,22,23: `tight ==
units` in every case. For n=14 this gives `{1,3,5,9,11,13}/14` — **exactly** the
"shared best times" list of HYP-2231, now with a one-line reason.

Two finer facts, also verified:

- **Binders come in negation pairs.** At a tight time `j/n`, the runners at the
  floor `1/n` are exactly `a ≡ ±j⁻¹ (mod n)` — the negation orbit `{b, −b}`.
  Since `b = j⁻¹` is a unit and `n/2` is a *non*-unit for `n>2`, a binder is
  **never** the self-conjugate speed.
- **The self-conjugate runner n/2 is inert-then-lethal.** For even n the speed
  `n/2` lies *inside* `AP_n`. At every unit time it sits at distance `1/2`
  (maximal — it never binds); at every even `j` it sits **on the origin**
  (`‖(n/2)(j/n)‖ = ‖j/2‖ = 0`). So `n/2` is the runner that deletes the even-j
  division times, and `t=1/2 = (n/2)/n` is killed by the whole even-speed
  sublattice.

## Why this sharpens (and partly corrects) S679

S679 wrote "apex = the sieve multiple of n = the unique negation fixed point 0."
For **prime** n those two notions *coincide* — the only fixed point is `0`, which
is also the `q=n` sieve-coverage residue, and it lives *outside* the AP. For
**even** n they **separate** into two genuinely different objects:

| fixed point | role at even n | lives where |
|---|---|---|
| `0`   | the `q=n` sieve apex (a forced multiple of n); kills *all* `j/n` | outside AP (external) |
| `n/2` | the AP's own self-conjugate runner; kills the *even-j* `j/n`, inert at unit times | **inside** AP |

So the even-n hard core is a **two-point** fixed set, and only one of the two
points (`0`) is a sieve obstruction; the other (`n/2`) is an internal structural
feature of the tight config. S679's identification "apex = fixed point" was a
prime-only accident of the two coinciding.

## The shell-side duality (why the codex machinery can't see it)

The pair-sum / worry-set modulus is `C = 2n−1`, which is **always odd**. Hence
negation on the shell torus `Z/C` always has the single fixed point `0`,
*regardless of the parity of n*. **The parity asymmetry lives entirely on the
speed torus `Z/n`, never on the shell torus `Z/C`.** This is exactly why the
`Res₂₇` carry/owner machinery (HYP-2164 → HYP-2253), which works mod `C=27`, has
to re-import parity as a *separate* coordinate (HYP-2230's `v mod 2 = r+k`): the
second apex `n/2` is invisible in the shell and must be carried by hand.

## Why n=19 is a cleaner lab than n=14 (matches S678/HYP-2254)

S678 ranked n=19 as the cleanest local-carrier lab (`n` prime, `C=37` prime) and
n=14 as the awkward one. The reframe explains it structurally:

| n | parity | `|Fix(σ_n)|` | φ(n) tight times | C = 2n−1 |
|---|---|---|---|---|
| 19 | odd  | **1** | 18 (all `j`)            | 37 (prime) |
| 14 | even | **2** | 6  `{1,3,5,9,11,13}`    | 27 = 3³ |

n=14 carries *two* obstructions the prime case lacks: the second negation fixed
point `n/2=7` **and** the composite shell `C=27=3³` (the `gcd` strata `1,3,9` of
the codex work). Both are absent at n=19. (Note: composite *odd* n=21 also drops
φ to 12, but keeps `|Fix|=1`; the extra fixed point is a pure **parity** effect,
orthogonal to the φ-drop.)

## An honest, optimistic handoff

The inert-then-lethal behaviour of `n/2` suggests a reduction, not just a
description: at every tight time the self-conjugate runner is pinned at the
*maximal* distance `1/2`, so it is **slack, never binding**. If one can show the
`n/2` runner stays slack throughout the off-grid recovery (not just on the grid),
then the even-n core collapses onto the *same* single-apex residual as the odd
case — i.e. the only genuinely hard fixed point is `0` (the sieve apex), and the
S606 twisted-involution cure that S679 wants for odd n would suffice for even n
too. That is the well-posed next object:

> **Conjecture (S700 handoff).** For even n, the self-conjugate runner `n/2` can be
> deleted from any tight/near-tight `Z/n`-config without lowering `M` below `1/n`
> off the division grid — reducing even-n LRC to the odd-n (single-apex) residual.

Concretely for n=14: is every near-floor config's `M` controlled by the
`{0}`-apex residual once the slack `7`-runner is set aside? This connects the
second apex directly to HYP-2253's "primitive apex debt" target and to the
`Vstar` seam (V* keeps speed 7 and the `3+24=27` pair-sum, and has the **same**
six unit tight times — verified).

## Status

- **Proved/verified (elementary, all n in 5..23):** AP tight times = units mod n
  (`φ(n)` of them); binders are negation pairs `±j⁻¹`; `Fix(σ_n) = {0}`/`{0,n/2}`
  by parity; the self-conjugate `n/2` runner is inert at unit times, lethal at
  even times; the shell modulus `C=2n−1` is always odd.
- **Reframe delivered:** even-n LRC has a *two-point* negation fixed set; the
  second apex `n/2` is the named parity obstruction, internal to the AP, invisible
  to the shell torus — sharpening (and correcting the parity-blindness of) S679.
- **Not proved:** LRC(14). The off-grid recovery is untouched; the handoff
  conjecture (slack-`n/2` reduction even→odd) is open.

Artifacts: `04-computation/lrc_negation_fixed_points_even_n_s700.py`,
`05-knowledge/results/lrc_negation_fixed_points_even_n_s700.out`, HYP-2259, T753.
Builds on S679 (HYP-2255), HYP-2231 (the {1,3,5,9,11,13} walls), HYP-2230 (parity
carry), HYP-2254 (n=19 cleaner), HYP-2164 (V*), THM-401 (pair-sum modulus C=2n−1),
and the S606 twisted-involution program.
