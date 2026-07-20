---
id: THM-1500
title: "THE GMC MASTER THEOREM: BOTH COUNTEREXAMPLES ARE ONE CONSTRUCTION, THE m! IS FORCED, AND n=3 IS MINIMAL FOR IT. Set P = (1+Z)(W − g(Z)U), Q = Z with U ≥ 0 independent of a standard complex Gaussian pair (Z,W). Then E[Pᵐ] = m!·[sᵐ](1+s)ᵐ·Φ(−s·g(s)) in ONE line, where Φ is the MGF of U — so everything is controlled by the single series f(s) := Φ(−s g(s)). (B) f IS UNIQUELY FORCED: by Lagrange inversion Σₘtᵐ[sᵐ]f(s)(1+s)ᵐ = f(t/(1−t))/(1−t), so demanding vanishing for every m ≥ 1 gives f(s) = f(0)/(1+s) and NOTHING ELSE; the mixed moment is then E[QPᵐ] = m!·[sᵐ]s(1+s)^{m−1} = m! AUTOMATICALLY. The m! is a theorem, not a coincidence of either example. (C) THE MASTER EQUATION Φ(−s g(s)) = 1/(1+s) with U = χ²_d/2 forces 1 + s·g(s) = (1+s)^{2/d}, i.e. g(s) = ((1+s)^{2/d} − 1)/s — recovering g = 2+s at d=1 (the owner's (2+Z), DERIVED not guessed) and g = 1 at d=2 (the S133 four-term example). The two known counterexamples are the SAME construction at d = 1 and d = 2. (D) MINIMALITY: g is a polynomial iff 2/d ∈ ℤ₊, i.e. d ∈ {1,2}, so n = 2+d ∈ {3,4} and n=3 is minimal; n=2 needs d=0, giving c·s·g(s) = log(1+s), and log(1+s)/s is not a polynomial — THIS FAMILY CANNOT REACH GMC(2). (E) This REFUTES my own HYP-8330 speculation that the alternating-binomial collapse might be the only mechanism and might force even n. The dimension table is corrected: GMC(N) is false for every N ≥ 3."
status: >
  (A) The GMC(3) counterexample is VERIFIED-EXACT two independent ways: the complex
  formalism (E[Z^aW^b] = a!δ_ab, E[U^k] = (d/2)_k) for m = 1..10, and a full real-coordinate
  expansion in a,b ~ N(0,1/2), T ~ N(0,1) for m = 1..6, giving E[Pᵐ] = 0 and E[QPᵐ] = m!
  with zero imaginary part. It genuinely lives in THREE real Gaussians.
  (B) PROVED (Lagrange–Bürmann with φ(s) = 1+s) and verified, including a perturbation
  spot-check: altering any single Taylor coefficient of f breaks the vanishing immediately.
  (C) PROVED and verified — d=2 re-derived and re-confirmed against the S133 example.
  (D) PROVED within the stated family. It is a statement about THIS construction, NOT a
  proof that GMC(2) holds — see Honest scope.
source: mac-mini-2026-07-20-S134 (owner supplied the GMC(3) counterexample and asked for the
  better next question on "binomial collapse", explicitly aiming at GMC(2)/(3)/(4))
supersedes_in_part: THM-1480   # dimension table N>=4 -> N>=3; HYP-8330 speculation refuted
depends_on:
  - THM-1480  # the GMC(4) counterexample, verified and reduced to four terms
related:
  - THM-1435  # "THE WITNESS IS NOT PRODUCED"
  - THM-1300  # the JC(3) counterexample the chain rests on
script: 04-computation/gmc3_master_theorem_macmini_S134.py (+ .out)
---

# THM-1500 — the GMC master theorem

**One line.** Both known Gaussian-moment counterexamples are the *same* construction with a
different `U`; the generating function that makes them work is **uniquely forced**; and that
uniqueness is exactly what caps the construction at three real Gaussians.

## (A) The GMC(3) counterexample is correct

With `Z = (X+iY)/√2`, `W = Z̄`, `U = T²/2`:

> `P = (1+Z)(W − (2+Z)U)`, `Q = Z` — a 5-term quartic in three real Gaussians.

Verified exactly `m = 1..10` in the complex formalism, and independently re-expanded into
real coordinates (`a, b ~ N(0,½)`, `T ~ N(0,1)`, 10 monomials) giving `E[Pᵐ] = 0` and
`E[QPᵐ] = 1,2,6,24,120,720` with **zero imaginary part**. It genuinely lives in three
independent real standard Gaussians.

## (B) The controlling series, and why it is forced

Take the construction in general: `P = (1+Z)(W − g(Z)U)`, `Q = Z`, with `U ≥ 0` independent
of `(Z,W)`. Using `E[W^r F(Z)] = r!·[s^r]F(s)` and writing `Φ` for the **moment generating
function** of `U`, one line handles all `m` at once:

> **`E[Pᵐ] = m!·[sᵐ] (1+s)ᵐ · Φ(−s·g(s))`**

So everything depends on the single series `f(s) := Φ(−s g(s))`. Now apply Lagrange–Bürmann
with `φ(s) = 1+s` (so `σ = t/(1−t)`, `φ' = 1`):

`Σ_m t^m [s^m] f(s)(1+s)^m = f(t/(1−t))/(1−t)`.

Demanding `[sᵐ](1+s)ᵐf(s) = 0` for **every** `m ≥ 1` says that series is the constant `f(0)`,
i.e. `f(t/(1−t)) = f(0)(1−t)`. Substituting `u = t/(1−t)`:

> **`f(s) = f(0)/(1+s)` — and nothing else works.**

Verified, including a perturbation check: bumping any single Taylor coefficient of `f`
destroys the vanishing at once (`f + s⁰, s¹, s², s³` all break it immediately).

**The `m!` is then automatic**, not a feature of either example:

`E[QPᵐ] = m!·[sᵐ] s(1+s)^{m−1} = m!·[s^{m−1}](1+s)^{m−1} = m!·1`.

## (C) The master equation, and both examples from it

> **`Φ(−s·g(s)) = 1/(1+s)`**

For `U = χ²_d/2` — i.e. `U` built from `d` real Gaussians — `Φ(x) = (1−x)^{−d/2}`, so

`(1 + s g(s))^{−d/2} = (1+s)^{−1}` ⟹ `1 + s g(s) = (1+s)^{2/d}` ⟹ **`g(s) = ((1+s)^{2/d} − 1)/s`**

| `d` | `n = 2+d` | `2/d` | `g(s)` | polynomial? |
|---|---|---|---|---|
| 1 | **3** | 2 | **`2 + s`** | yes |
| 2 | **4** | 1 | **`1`** | yes |
| 3 | 5 | 2/3 | `2/3 − s/9 + 4s²/81 − …` | **no** |
| 4 | 6 | 1/2 | `1/2 − s/8 + s²/16 − …` | **no** |

> `d = 1` gives `g = 2+s` — **exactly the `(2+Z)` in the owner's construction, derived rather
> than guessed.** `d = 2` gives `g = 1` — **exactly the four-term S133 example**
> `P = (1+W)(W̄ − |Z|²)`. The two counterexamples are one construction at `d = 1, 2`.

## (D) Minimality: why three, and why not two

`g` is a **polynomial** iff `2/d` is a positive integer, i.e. **`d ∈ {1,2}`**. Since the total
real dimension is `2` (for `Z,W`) plus `d`:

> **`n ∈ {3,4}`, and `n = 3` is minimal for this family, attained uniquely at `d = 1`.**

`n = 2` would require `d = 0`, i.e. `U` constant `= c`. Then `Φ(x) = e^{cx}` and the master
equation becomes `c·s·g(s) = log(1+s)`. But `log(1+s)/s = 1 − s/2 + s²/3 − …` is **not a
polynomial**, so no such `g` exists.

> **This family cannot reach GMC(2).**

The same obstruction holds in the Gamma family: `U ~ Γ(α)` forces `1 + sg = (1+s)^{1/α}`,
polynomial iff `1/α ∈ ℤ₊`, and only `α = ½` (`d=1`) and `α = 1` (`d=2`) are realizable as
sums of squares of Gaussians.

## (E) What this corrects

- **The dimension table becomes `GMC(N)` false for every `N ≥ 3`** (pad with unused
  coordinates), superseding THM-1480's `N ≥ 4`.
- **HYP-8330 is REFUTED.** I speculated there that the alternating-binomial collapse might be
  the only mechanism, and that its dependence on `k!` "may force an even number of real
  variables, making n=4 sharp." Both halves are wrong: the `d=1` case runs on
  `(1+x)^{−1/2}` composed with a perfect-square substitution, not on `Σ(−1)ʲC(m,j)`, and it
  lands on **odd** `n = 3`. The correct statement is the uniqueness of `f = 1/(1+s)`, of which
  the alternating sum was only the `d = 2` shadow.

## Honest scope

- **(D) is minimality WITHIN THE FAMILY `P = (1+Z)(W − g(Z)U)` with `U` independent of
  `(Z,W)`.** It is *not* a proof that GMC(2) is true. A GMC(2) counterexample would simply
  have to come from a different shape — e.g. `U` correlated with `(Z,W)`, a different outer
  factor than `(1+Z)`, or `Q` not linear. None of those are excluded here.
- The `Φ(−sg) = 1/(1+s)` derivation assumes `U` independent of `(Z,W)`; that independence is
  what factorizes the expectation and is load-bearing throughout.
- The Lagrange-inversion uniqueness (B) is a statement about formal power series with
  `f(0) ≠ 0`; the normalization `f(0) = 1` is implicit in `Φ(0) = 1`.
- **`U = χ²_d/2` is a choice, not a classification.** Other `U` realizable as polynomials in
  Gaussians (not sums of squares) are untested and are the obvious place a `d = 0`-like
  evasion could hide.
- Neither this file nor THM-1480 claims priority on the GMC(3) example — it was supplied by
  the owner from an outside source. What is claimed here is the independent verification, the
  master theorem, the uniqueness, and the minimality analysis.

*Artifacts:* `04-computation/gmc3_master_theorem_macmini_S134.py` (+out).
