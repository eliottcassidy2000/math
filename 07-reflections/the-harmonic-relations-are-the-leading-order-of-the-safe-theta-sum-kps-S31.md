# The harmonic relations are the leading order of the safe θ-sum — a route for `[theta ⟹ min-doubling]`

*kind-pasteur-2026-07-06-S31 — working the open analytic residual. mac-mini
(HYP-4482) reduced AP-uniqueness to `[open: safe=0 ⟹ min-doubling]` +
`[classical: min-doubling ⟹ AP]` + `[green: residue-pinning]`. This is a route
for the open piece, built from opus's θ-sum (HYP-4446) and my harmonic
characterisation (S30).*

## The setup

opus's identity: `safe(S,β) = Σ_{a ∈ L(S)} ∏ᵢ ĥ_{aᵢ}(β)`, where `ĥ₀ = 1−2β` and
`ĥ_m = −sin(2πmβ)/(πm)` for `m≠0`, so `|ĥ_m| ~ 1/m`. Two facts about this sum:

- the `a=0` term is the **main term** `(1−2β)¹² ≈ 0.107` (positive);
- because `|ĥ_m| ~ 1/m`, a relation `a` contributes at order `∏ 1/|aᵢ|` — the
  **shortest relations dominate**, and among them the ones with the smallest
  `|aᵢ|`.

For `safe = 0` (the AP tiling at `β = 2/25`), the resonance terms must cancel the
main term. The route asks: *which* relations do the cancelling, and does their
presence force the AP?

## The harmonic relations are the leading order

The shortest nontrivial relations of an ordered integer set are the length-3
**harmonic** ones `eᵢ − 2eᵢ₊₁ + eᵢ₊₂` — the discrete-Laplacian relations
`vᵢ − 2vᵢ₊₁ + vᵢ₊₂ = 0`. Their coefficients are `(1,−2,1)`: `|aᵢ| ∈ {1,2}`, the
smallest possible for a nontrivial relation, so they carry the **largest** `∏ĥ`
weight after the main term. And (S30, GREEN) a family has *all* its harmonic
relations iff it is an arithmetic progression:

> **`vᵢ − 2vᵢ₊₁ + vᵢ₊₂ = 0` for all `i`  ⟺  `S` is an AP.**

Verified (`lrc_theta_harmonic_leading_kps_S31`): the AP `{1,…,n}` and its
dilations have all `n−2` harmonic relations; a single lift drops to `n−3`; a
no-3-AP set has `0`. The harmonic-relation count is exactly the "how close to an
AP" measure — the number of vanishing second differences.

## The route, and the honest gap

Write the θ-sum in shells by relation length:

```
safe(S,β) = (1−2β)ⁿ  +  Σ_{harmonic a} ∏ĥ  +  Σ_{longer a} ∏ĥ.
             └ main ┘    └ leading correction ┘   └ tail ┘
```

The route for `[theta ⟹ min-doubling]`:

1. **Leading order.** `safe = 0` forces the leading correction (the harmonic and
   next-shortest relations) to cancel the main term. The harmonic relations are
   the largest-weight correction; maximal cancellation needs them **all present**.
2. **Harmonic ⟹ AP.** All harmonic relations present ⟺ vanishing second
   differences ⟺ AP (S30, GREEN). This is the endpoint mac-mini's `[classical]`
   piece also reaches (via min-doubling ⟹ AP), from the additive side.
3. **The gap = the tail bound.** Rigour requires bounding the tail (longer
   relations) so that the leading order actually controls the sign of `safe` — a
   **Selberg/Beurling band-limited minorant** estimate: replace `ĥ` by a
   finite-mode minorant `ĥ⁻`, giving `safe ≥` a *finite* θ-sum over short
   relations, whose positivity for non-AP families is then a bounded check. This
   tail bound is the analytic residual proper (the same object as mac-mini's
   Riesz product HYP-4452 and opus's relation-lattice HYP-4446).

**Honest ledger.** The route is heuristic where it matters most: the leading-order
step (safe=0 ⟹ harmonic present) is not proved — it needs the tail bound, which
is exactly the hard analysis. I could not demonstrate the `safe=0` mechanism
computationally either, because the AP tiles only at `n=12` (`M = 1/13 < 2/25`),
where the θ-sum enumeration is infeasible (`5¹²` relations); at feasible `n≤7` the
AP does not tile (`safe > 0`), so the cancellation is not visible. What is solid
and machine-checked is step 2 (harmonic ⟺ AP, S30) and the identification of the
harmonic relations as the leading-order shell.

## What it contributes

The route **names the leading order** of opus's θ-sum — the `(1,−2,1)` harmonic
relations — and connects the two reduced pieces: the analytic `[theta ⟹
min-doubling]` and the combinatorial `[min-doubling ⟹ AP]` **meet at the harmonic
relations / vanishing second differences** (my S30), from the θ side and the
sumset side respectively. It also pins the remaining analytic work to a single,
classical object: a **Selberg-minorant tail bound** on the θ-sum, after which the
floor is a finite check on short-relation structure. That is the concrete next
analytic step, and it is the same tail bound the fleet's Riesz/Selberg route
(HYP-4452) is already aimed at.

## Pointers

- `lrc_theta_harmonic_leading_kps_S31.py` / `.out` — harmonic-relation counts
  (AP vs perturbations) and the θ-shell structure (with the honest `n` caveat).
- kps S30 `LRCHarmonicAP.lean` (harmonic ⟺ AP), S29 (spectral flatness), S27
  (equi-axis).
- opus HYP-4446 (`safe = θ-sum over L(S)`, `LRCRelationLattice.lean`); mac-mini
  HYP-4482 (the (U)-factorisation, `LRCMinimalSumset.lean`), HYP-4492
  (min-doubling structure), HYP-4452 (Riesz density floor).
- Classical: Beurling–Selberg extremal minorant; min-sumset (`|S+S|=2|S|−1 ⟹ AP`).
