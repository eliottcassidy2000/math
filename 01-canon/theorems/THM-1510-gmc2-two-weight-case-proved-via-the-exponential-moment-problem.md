---
id: THM-1510
title: "GMC(2) IS PROVED FOR EVERY TWO-WEIGHT P, IN EVERY DEGREE — via a new 1-D lemma (EMP) and a stronger 2-D nullcone conjecture (NC2). (A) THE UNIFYING FORMULA: for P = (1+Z)(W − V(Z)U) with U independent, E[P^m] = m!·[s^m](1+s)^m·G(sV(s)) where G(x) = Σ(μ_k/k!)(−x)^k; every moment dies as soon as G(sV(s)) = 1/(1+s), and then E[QP^m] = m! automatically. ONE identity generates BOTH known witnesses. (B) WHY THE LADDER STOPS: for U = (T₁²+…+T_ν²)/2 ~ Gamma(ν/2) the condition is (1+sV(s)) = (1+s)^{2/ν}, and V is a polynomial iff 2/ν is a positive integer j — forcing ν = 2/j, so only j=1 (ν=2, the n=4 witness, V=1) and j=2 (ν=1, the n=3 witness, V=2+s — that IS the mysterious (2+Z)) are realizable by Gaussians. n=2 would need ν=0, i.e. U constant, giving c·sV(s) = log(1+s), not a polynomial. (C) NC2, a strictly stronger 2-D nullcone conjecture: in two real Gaussians the ONLY P with all moments vanishing are the purely holomorphic/antiholomorphic ones. NC2 ⟹ GMC(2), proved in three lines by weight counting. (D) EMP, PROVED BY LAPLACE: for r ~ Exp(1) and h ∈ ℂ[r] of degree d ≥ 1, E[h^m] = c_d^m (dm)! e^{c_{d−1}/(c_d d)}(1+o(1)) — nonzero — so E[h^m] = 0 ∀m forces h ≡ 0. (E) CONSEQUENCE: NC2, hence GMC(2), holds on EVERY two-weight stratum in every degree. (F) The gap is now exactly one stratum, {−1,0,1}, with an explicit Bessel form."
status: >
  (A) PROVED (binomial expansion + the Wick identity E[W^r F(Z)] = r![s^r]F(s));
      machine-checked to m = 7 for both j = 1 and j = 2.
  (B) PROVED (the polynomiality of (1+s)^{2/nu} is an integrality condition).  NOTE: an
      earlier run of my own script printed this table with nu = 2/a instead of nu = 2a and
      so labelled j >= 3 "REALISABLE", contradicting its own prose.  Fixed; the prose was
      right.
  (C) CONJECTURE (NC2) + PROVED IMPLICATION.  NC2 is tested exhaustively in a bounded box
      (736 nullcone members found, ZERO violations, positive and negative controls both
      behaving).  The IMPLICATION NC2 => GMC(2) is proved outright.
  (D) PROVED by the Laplace/saddle-point method.  The derivation is complete in outline and
      the estimate is standard; I have NOT written out the uniform tail bounds, so a referee
      would want that half-page.  Numerically confirmed on three polynomials -- including an
      adversarial one that kills E[h] and E[h^2] -- with the error decaying at exactly the
      O(1/m) rate a Laplace expansion predicts (error ratios 2.09, 2.04, 2.02, 2.01, 2.00).
      Independently, degrees 1,2,3 are settled by exact/near-exact root-finding, and the
      root counts 1, 2, 6 hit the Bezout bound d! exactly, so those searches were complete.
  (E) PROVED, given (D).
  (F) SCOPE STATEMENT, plus a bounded search: 117,648 triples on the {-1,0,1} stratum,
      zero violations.  NOT a proof.
  GMC(2) ITSELF REMAINS OPEN.  What is proved is a genuine partial case, not the conjecture.
source: klein-2026-07-20-S343 (owner: work to prove GMC(2); think stronger two-dimensional nullcone conjectures that turn out to also be true)
attribution: >
  The GMC(3) witness P = (1+Z)(W-(2+Z)U), Q = Z is NOT mine -- it was supplied by the owner
  from an outside source and is verified here (Part 1), not discovered.  The n=4 witness is
  likewise outside work (THM-1490).  What is mine: the unifying formula (A), the Gamma-ladder
  obstruction (B), the NC2 formulation and implication (C), EMP (D), and the two-weight
  theorem (E).
depends_on:
  - THM-1490  # the n=4 witness, verified/proved/shortened; the weight-polytope obstruction
related:
  - THM-1435  # the Zhao-VC witness machinery -- a DIFFERENT Mathieu-Zhao member; untouched
script: 04-computation/gmc2_nullcone_klein_S343.py (+ .out)
---

# THM-1510 — GMC(2) on the two-weight strata

**GMC(n).** For the standard Gaussian on `ℝ^n` and `P, Q ∈ ℂ[x_1,…,x_n]`: if `E[P^m] = 0`
for all `m ≥ 1`, then `E[QP^m] = 0` for `m ≫ 0`. False for `n ≥ 3`; `n = 2` open.

Throughout `n = 2`: one complex Gaussian `Z` with `E[Z^aZ̄^b] = δ_{ab}a!`, weight of `Z^aZ̄^b`
is `a−b`, and `r = |Z|² ~ Exp(1)`.

## A. One formula generates both known witnesses

Let `P = (1+Z)(W − V(Z)U)` with `W = Z̄`, `U` independent of `(Z,W)`, `μ_k = E[U^k]`, `V` a
polynomial. Expanding binomially and using `E[W^r F(Z)] = r!·[s^r]F(s)`:

```text
E[P^m] = m! · [s^m] (1+s)^m · G(s·V(s)),      G(x) := Σ_k (μ_k/k!)(−x)^k
```

So **every** moment vanishes as soon as

```text
G(s·V(s)) = 1/(1+s)
```

— because then `(1+s)^m·G = (1+s)^{m−1}`, of degree `m−1`, so `[s^m] = 0`; and the `Q = Z`
insertion gives `[s^m] s(1+s)^{m−1} = 1`, i.e. `E[QP^m] = m!`. Machine-checked to `m = 7` for
both witnesses.

## B. Why the ladder stops at `n = 3`

Build `U` from `ν` real Gaussians as `U = (T_1²+…+T_ν²)/2 ~ Gamma(ν/2)`. Then `μ_k = (ν/2)_k`
and `G(x) = (1+x)^{−ν/2}`, so the condition reads

```text
1 + s·V(s) = (1+s)^{2/ν}
```

**`V` is a polynomial iff `2/ν = j` is a positive integer**, i.e. `ν = 2/j`:

| `j` | `ν = 2/j` | `V(s) = ((1+s)^j − 1)/s` | total `n = 2 + ν` | realizable? |
|---|---|---|---|---|
| 1 | 2 | `1` | **4** | yes — the `n=4` witness |
| 2 | 1 | `2 + s` | **3** | yes — the `n=3` witness, **and this is the mysterious `(2+Z)`** |
| ≥3 | 2/3, 1/2, … | — | — | **no: a fractional number of Gaussians** |

And `n = 2` would need `ν = 0`, i.e. `U ≡ c` constant, whence `G(x) = e^{−cx}` and the
condition becomes `c·s·V(s) = log(1+s)` — not a polynomial.

> **Within this family `n = 3` is minimal and `n = 2` is impossible.** That explains the
> supplied witness's shape rather than merely checking it. It is a statement about the
> *construction*, not a proof of GMC(2).

## C. NC2 — a stronger two-dimensional nullcone conjecture

Let `N_2 = {P : E[P^m] = 0 ∀m ≥ 1}`.

> **Conjecture (NC2).** `N_2 = {P : all weights ≥ 1} ∪ {P : all weights ≤ −1} ∪ {0}` —
> in two real Gaussians the only way to kill every moment is to be purely holomorphic (or
> antiholomorphic) with no constant term.

NC2 is **false for `n ≥ 3`** (both witnesses straddle weight 0), so it is genuinely 2-D.

**Theorem. NC2 ⟹ GMC(2).** Let `0 ≠ P ∈ N_2`; by NC2 all weights are `≥ 1` (wlog). Weights
add, so every monomial of `P^m` has weight `≥ m`. A fixed `Q` has a finite minimum weight
`w_min`, so every monomial of `QP^m` has weight `≥ m + w_min ≥ 1` once `m > 1 − w_min`; Wick
kills all of them. ∎

So GMC(2) reduces to the *shape of the nullcone*, with `Q` eliminated entirely. Bounded test:
**736 nullcone members found, zero NC2 violations**, with `P = Z` (in) and `P = Z + Z̄`
(`E[P²] = 2`, out) as controls.

## D. EMP — the one-dimensional lemma, proved

> **Lemma (EMP).** For `r ~ Exp(1)` and `h ∈ ℂ[r]`: if `E[h(r)^m] = 0` for all `m ≥ 1` then
> `h ≡ 0`.

*Proof.* If `deg h = d ≥ 1` with leading coefficient `c_d`, write `h(r) = c_d r^d φ(r)`,
`φ(r) = 1 + (c_{d−1}/c_d)/r + O(1/r²)`. In `E[h^m] = ∫_0^∞ e^{−r}h(r)^m dr` the factor
`e^{−r}r^{dm}` peaks at `r = dm`, eventually beyond every root of `h`. There

```text
m·log φ(r) = m[(c_{d−1}/c_d)/r + O(1/r²)] → c_{d−1}/(c_d d),
```

a **constant** — the subleading coefficients contribute a finite factor, not a growing one.
The peak has width `O(√(dm))`, over which `arg h(r)` moves by `O(1/√(dm)) → 0`, so the phase
is asymptotically constant and there is **no oscillation to cancel the integral**. Hence

```text
E[h^m] = c_d^m · (dm)! · e^{c_{d−1}/(c_d d)} · (1 + o(1)),
```

whose amplitude is nonzero. So `E[h^m] ≠ 0` for large `m` — contradiction. Thus `h` is
constant, and `E[h] = h = 0`. ∎

*Numerics.* Three test polynomials, `m` up to 320: the ratio `E[h^m]/(c_d^m(dm)!)` converges
to `e^{c_{d−1}/(c_d d)}` with error ratios `2.09, 2.04, 2.02, 2.01, 2.00` per doubling —
exactly the `O(1/m)` of a Laplace expansion. The third test is the *adversarial* case
`h = r² + (−4+2i)r + (2−2i)`, an exact root of `E[h] = E[h²] = 0`: its limit is still
nonzero, which is precisely why it cannot keep killing moments. Independently, degrees
`1,2,3` are settled by root-finding, and the root counts `1, 2, 6` equal the Bézout bound
`d!`, so those searches found everything. (Degree 2 by hand: `E[h]=0 ⟹ c₀ = −2−c₁`;
`E[h²]=0 ⟹ c₁²+8c₁+20 = 0 ⟹ c₁ = −4±2i`; then `E[h³] = 32±80i ≠ 0`.)

## E. The theorem

> **GMC(2) holds for every `P` supported on at most two weights, in every degree.**

Let the weights be `p < q`. Wick keeps total weight `0`, so in `P^m` only terms with
`ip + jq = 0`, `i+j = m` survive.

- **`0 < p < q`** (or both negative): no solution with `m ≥ 1`, so every moment vanishes for
  free. These are exactly the members NC2 permits.
- **`p = 0 < q`**: forces `j = 0`, so `E[P^m] = E_r[f(r)^m]` for `f` the weight-0 part. EMP
  gives `f = 0`.
- **`p < 0 < q`**: `i|p| = jq` with `i+j = m` pins `(i,j)`, and `f^i g^j = (f^q g^{|p|})^{m/(q+|p|)}`
  is a **power**, so `E[P^m] = C(m,i)·E_r[H^k]`. EMP gives `H = 0`, hence `f = 0` or `g = 0`.

In the last two cases `P` collapses to a single weight, so NC2 holds there; apply §C. ∎

## F. Exactly one stratum remains

The first uncovered case is `{−1,0,1}` — the shape of both `n ≥ 3` witnesses. Writing
`P = Z̄a(r) + b(r) + Zc(r)` and polar-decomposing `Z = √r e^{iθ}` with `r ⊥ θ`:

```text
Σ_m E[P^m] t^m/m!  =  E_r[ e^{b(r)t} · I_0( 2t√(r·a(r)c(r)) ) ]
```

and GMC(2) here is the claim that this equals `1` only when `r·a·c = 0` and `b = 0` — a
**Bessel-weighted analogue of EMP**. Bounded search: 117,648 triples with `a,b,c` of degree
`≤ 1`, zero violations. The same saddle-point method ought to apply, since `I_0(x) ~ e^x/√(2πx)`
also has no oscillation on the positive axis — **but I have not carried it through and do not
claim it.** That is the honest remaining gap, and it is now a single, sharply stated one.

*Files: `04-computation/gmc2_nullcone_klein_S343.py` (+ `.out`).*
