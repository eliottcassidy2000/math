---
source: opus-2026-06-03-S596 (remote-control)
status: SYNTHESIS — Collatz and the LRC residual are the SAME question: a 2-adic ×odd/÷2 orbit hitting a target, easy by measure/density, hard on a measure-zero arithmetic residual that reduces to a bounded TWO-BLOCK (2-power vs 3-power / 2-adic × odd) Diophantine obstruction
tags: [collatz, LRC, two-block, 2-adic, prime-2, prime-3, baker, vitali-wall, rank-1, residue-automaton, n14]
---

# Collatz is the same two-block question as the LRC residual

**Prompt (user):** consider the Collatz conjecture and how it is a similar question to
work already in the repo.

It is the **same shape** as where the LRC arc landed (S589–S595): a `×odd / ÷2` orbit
asking whether it reaches a target, **easy by measure/density**, **hard on a measure-zero
arithmetic residual** that collapses to a **bounded two-block (2-power vs 3-power)
Diophantine obstruction** — and the two problems share their failure mode of measure
(the Vitali wall) and their proof technique (a Baker-type bound on the two-block gap).

## 1. The shared question

| | LRC (this repo) | Collatz |
|---|---|---|
| dynamics | doubling `x↦2x mod n` / runner flow `{v_i t}` | `x ↦ 3x+1`, then `÷2` |
| target | the safe box (loneliness, `M≥1/n`) | reach `1` |
| **easy half** | positive-measure bulk is lonely (S550 measure bound) | density-1 reach 1 (Terras/Korec) |
| **hard residual** | the measure-zero **worry-set** | the density-0 exceptional set (cycles / divergence) |
| **why measure fails** | Vitali-blind at the floor (S551o) | density-blind on the density-0 exceptional set |
| 2-adic structure | doubling = the binary shift (S580); the 2-adic tower (S579) | conjugate to the 2-adic shift (Lagarias) |
| generation / reverse | binary IFS `{x↦2x, x↦2x+1}` (S580) | reverse tree `{x↦2x, the ÷3 branch}` |
| prime split | prime-2 (doubling) vs odd (`2n−1` shells); `n=14`: **2 vs 3** (S589–595) | prime-2 (`÷2`) vs prime-3 (`×3`) |
| residual obstruction | **rank-1 two-block** `det[u_a,r_a;u_b,r_b]=w·n·u_a u_b·ℓ`, bounded slacks (S595) | **two-block** `2^E − 3^k` divides a bounded sum `S` |
| elimination | bounded CRT automaton (verified empty) | linear forms in logs / Baker (no cycles, bounded `k`) |

Both are "**does a 2-adic-flavoured `×odd/÷2` orbit hit the target**", and in both the *only*
obstruction lives at the **prime-2 vs prime-3** seam — the exact place the `n=14` residual
localized to (S589, S592).

## 2. The two-block, concretely (the anchor)

A Collatz cycle on odd `a_1→⋯→a_k` (`a ↦ (3a+1)/2^{e}`, `E=Σe_i`) satisfies the **cycle
equation**
```
a_1 · (2^E − 3^k) = S,    S = Σ_{i=1}^{k} 3^{k−i} 2^{e_1+⋯+e_{i−1}}  (a bounded sum),
```
so a nontrivial cycle needs `(2^E − 3^k) | S` with `a_1` a positive odd integer. The
**two blocks** are the `2`-power `2^E` (the `÷2` / 2-adic side) and the `3`-power `3^k`
(the `×3` / odd side); the obstruction is their **difference dividing a bounded sum**.

This is *exactly* the LRC S595 form: there the cover needs
`det[u_a,r_a;u_b,r_b] = w·n·u_a u_b·ℓ` with **bounded slacks** `|r|<u/n` — a two-block
(owner / slack, i.e. odd × 2-adic) determinant hitting a target it generically cannot.

**Verified** (`collatz_two_block_parallel_s596.py`): the only cycle with `e_i≤5, k≤7` is
the trivial `(1)` (`k=1, E=2, 2²−3=1`); and `|2^E−3^k|/2^E` stays **bounded away from 0**
(min `≈0.0136` at `k=12`) — the Baker / linear-forms-in-logs phenomenon that the two-block
gap is never tiny enough to admit a cycle. Mirror of S595's "bounded CRT automaton empty"
(the bounded slacks never hit the determinant target).

## 3. The two payoffs (each problem lends the other a tool)

**(a) Import to LRC — a Baker-type closure of the residual.** The LRC large-owner residual
(S595) is a *bounded two-block Diophantine* exactly like the Collatz cycle equation. Collatz
*cycles are eliminated for bounded `k`* by **linear forms in logarithms** (`|E log2 − k log3|`
bounded below ⟹ `|2^E−3^k|` not too small ⟹ no division). The same instrument should close
the LRC residual: bound `|det − w n u_a u_b ℓ|` below by a linear-forms estimate on the
owner logs, so the bounded slacks can never hit the target. **This is the concrete
technique transfer the prompt invites:** replace "verified empty automaton" by a Baker
bound on the two-block gap.

**(b) Export to Collatz — the Vitali-wall framing.** The repo's S551o lesson — *measure is
blind on the measure-zero residual; the exceptional set is constructible, not pathological,
and needs arithmetic not density* — is exactly the Collatz situation: density-1 reaches 1
(measure half, done), but the exceptional set is density-0 and density-blind; the proof must
be **arithmetic at the 2-adic / two-block level**, never a density estimate. The repo's
"sidestep resonance energy" (THM-401) = Collatz's "sidestep the density/drift heuristic."

## 4. Why the analogy is tight, not loose

- Both reduce to a **`2`-power-vs-odd-power bounded Diophantine** (`2^E−3^k` / the two-block
  determinant), with the obstruction a **low-rank / divisibility degeneracy** that a
  **bounded-window + linear-forms** argument rules out.
- Both have the **measure/density easy half** and the **measure-zero/density-zero hard
  residual** separated by the **Vitali wall** (S551o).
- Both localize the hard residual to the **prime-2 vs prime-3** interplay — the `÷2` vs `×3`
  of Collatz is the doubling (prime-2, THM-404) vs the `2n−1=3³` shells (prime-3, S592) of
  `n=14`.

## 5. Honest status

- **Established/verified:** the Collatz cycle equation and its two-block form; only the
  trivial cycle for `e≤5,k≤7`; the gap `|2^E−3^k|/2^E` bounded below; the structural
  correspondence table.
- **The transferable claims (directional):** (a) a Baker / linear-forms bound on the LRC
  two-block determinant gap would close the large-owner residual (S595) the way it closes
  Collatz cycles — *untested*, the concrete next instrument; (b) the Vitali-wall framing
  applies to Collatz's exceptional set — a viewpoint, not a theorem.
- **Not claimed:** any new resolution of Collatz or LRC. The contribution is the *identity
  of the questions* and the *shared technique* (bounded two-block Diophantine + linear forms
  in logs), with the Vitali wall explaining why both resist measure.

**Artifacts:** `04-computation/collatz_two_block_parallel_s596.py` (+`.out`). Builds on S595
(rank-1 two-block / bounded CRT automaton), S589/S592 (prime-2/3 localization), S580
(binary IFS / 2-adic shift), S550/S551o (measure bound / Vitali wall), THM-401/404. New:
**HYP-2143**.
