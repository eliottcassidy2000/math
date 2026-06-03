---
source: claude-2026-06-03-S602
status: framework + proved first-moment law + census + falsifiable conjectures
tags: [LRC, covering-depth, arc-covering, tight-config, extremal, entropy, additive-chain, recursive-log, helly]
---

# LRC is a circular-arc covering problem; the depth distribution is the master object

## Repository map

This entry is HYP-2156 after renumbering from collisions with the determinant
Helly logarithm HYP-2152 and S599c's THM-406/HYP-2155.  It sits between HYP-2154, the broad
covering-depth/master-object abstraction, and HYP-2153, the concrete `p_0`
endpoint/shell classifier for the additive-chain collapse family.

## The abstract move

A loneliness certificate is a point `t ∈ ℝ/ℤ` of the clock that avoids every
**forbidden arc** `F_i = {t : ‖v_i t‖ < δ}`. So LRC at gap `δ` is exactly a
**circular-arc covering** question: do the `n` forbidden sets cover the circle?
Each `F_i` is `v_i` equally spaced arcs of total measure `2δ`.

The right object is not "is there a gap?" but the whole **covering-depth
distribution**

```
p_k = meas{ t : depth(t) = k },   depth(t) = #{ i : ‖v_i t‖ < δ }.
```

Lonely times are `{depth = 0}`, of measure `p_0`. Every LRC quantity at gap `δ`
is a functional of this one distribution:

- a lonely time exists a.e. ⇔ `p_0 > 0`;
- **tight / extremal** at `δ = 1/(n+1)` ⇔ `p_0 = 0` (arcs cover up to measure 0);
- loneliness radius `δ_max(V) = sup{ δ : p_0(δ) > 0 }`.

## The one law that is actually a theorem

`∫ depth = Σ_i meas(F_i) = 2nδ`, so

```
Σ_k k·p_k = 2nδ = 2n/(n+1)  at the threshold.
```

The mean depth is **fixed at `2n/(n+1) < 2` for every speed set** — only the
*shape* of `p_k` varies. So the entire LRC problem at the threshold is: given a
fixed mean just under 2, can the distribution avoid putting mass at `k = 0`?

## The order parameter: depth entropy

Among all `V` (same fixed mean), the **depth-distribution entropy `H(p)` is
minimized by tight configurations**, and the AP is the global minimum at each
`n`:

| n | tight `H(p)` | non-tight `H(p)` |
|---|---|---|
| 4 | 0.97 (AP), 1.03 | 1.22 |
| 5 | 1.02 (AP), 1.09 | 1.31 |
| 6 | 1.09 (AP) | 1.31 |

Non-tight configs *pay entropy* for their wasted lonely measure `p_0 > 0`.
Tightness is the **entropy-minimal / most efficient circular cover**. This is the
continuous refinement the recursive-log/Helly frame asks for: the binary Helly
test (do the arcs cover?) is the sign of `p_0`; `H(p)` grades how efficiently.

**Conjecture (clean optimization):** the AP `{1,…,n}` globally minimizes the
depth entropy `H(p)` over all primitive speed sets at `δ = 1/(n+1)`.

## The `p_0 = 0` collapse family (the new sub-problem)

Which `V` are tight at `1/(n+1)`? The AP always is. The prompt's observation —
the family is **larger than the AP**, with sporadic additive chains `(1,3,4,7)`,
`(1,3,4,5,9)` — is correct and sharpens to:

**Census (exhaustive primitive `gcd=1`, `n≤6`, `B≤20`; chain-only to `B≤30`):**

```
n=4: (1,2,3,4)[AP], (1,3,4,7)
n=5: (1,2,3,4,5)[AP], (1,3,4,5,9)
n=6: (1,2,3,4,5,6)[AP]
n=7: (1,2,3,4,5,6,7)[AP], (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)
```

- **Additive-chain is NECESSARY** (every tight set has each element beyond the
  two smallest equal to a sum of two earlier ones) — verified exhaustively for
  `n ≤ 6` — but **far from sufficient** (`n=6`: 257 additive chains, only the AP
  is tight). So "additive chain" is a real but loose necessary condition; the
  true characterization is finer.
- The AP is universal; the sporadic tights are rare and irregular.

## Negative results that save effort

The two natural "explanations" for the sporadics are **coincidences**:

- **Lucas** `(1,3,4,7,11,18,…)` is tight **only at `n=4`** (`(1,3,4,7,11)` is
  not). So `(1,3,4,7)` is not the head of a Lucas-tight family.
- **Paley QR mod `p`** is tight **only at `p=11`** (not `7,13,17,19,23`). So
  `(1,3,4,5,9) = QR₁₁` is an isolated coincidence, not a Paley-tight family.

These matter: they tell the next investigator *not* to chase a Lucas or Paley
classification of the collapse family. It is genuinely sporadic beyond the AP —
the same flavor as the repo's other "sporadic tight" finds (S385, S553).

## Tie to the recursive-log / Helly-entropy frame (HYP-2151)

`p_0 = 0` is precisely the **Helly-infeasible boundary**: the forbidden arcs
cover, so no positive-measure lonely point survives (only measure-zero
witnesses). The depth distribution `p_k` is the covering-**multiplicity** profile
that a binary Helly test collapses to one bit; `H(p)` is its entropy refinement;
`Λ(V) = p_0` is the continuous order parameter whose **zero set is the collapse
family** — the LRC analogue of an exact cover. The Tao-style loglog improvements
to LRC lower bounds are estimates of exactly this `p_k` tail, so `p_k` is the
master object those iterated-log bounds are about.

## Open / next

1. **Characterize the sporadic tights** beyond "additive chain." Candidate finer
   invariants: the witness lattice (S553's half-division points), or the *exact*
   depth profile `(p_1, p_2, …)` — do the sporadics share a profile signature?
2. **Prove the first-moment-plus-tightness rigidity:** with mean fixed at
   `2n/(n+1)` and `p_0 = 0`, how much is `(p_1,…)` constrained? Is there a
   second-moment identity forcing the additive-chain relations?
3. **Prove the AP entropy-minimization conjecture** (a clean variational problem
   on `p_k`).
4. Feed back into the Helly-entropy ledger: the depth-entropy `H(p)` vs the
   obstruction rank of HYP-2151 — is the sporadic-vs-AP entropy gap a rank gap?
