---
source: claude-2026-06-06-S686
status: verified (0 support-violations n≤6; A(T)-refines-H exhaustive n=5) + the orientation-free correction of s685
tags: [arc-flip, delta, hessian, OCF, walsh, odd-cycle, invariant, H-spectrum, arc-interaction-graph, tournaments]
---

# The arc-interaction graph A(T): orientation-free support, and a strictly finer invariant than H

The user's thesis was that "the exact pattern arc flipping has on the deltas of a
tournament is the key to its structure." HYP-2268 (opus-S699i) cast that pattern
as the discrete **Hessian** `M(g,e) = H(T) − H(T^g) − H(T^e) + H(T^{ge})` — the
change in `δ(e)` caused by flipping `g` — and read off its Walsh/OCF form
`M = 4·Σ_{S⊇{g,e}} c_S χ_S` over collections `S` of vertex-disjoint odd cycles.
This session pins that down into two concrete, verified facts about the object I
name the **arc-interaction graph** `A(T)`: vertices = the `C(n,2)` arcs, edge
`g~e` iff `M(g,e)≠0`.

## 1. The support of `M` is orientation-free (correcting s685)

My earlier draft (s685) claimed: disjoint arcs interact **iff a common odd cycle
of `T`** passes through both. That used the *directed* cycles of the specific
tournament — and it is wrong as an invariant, because a **transitive** tournament
has no directed cycles at all yet its sharing-arc Hessians are all nonzero.

The correct statement: `c_S` are the **global** Walsh coefficients; only the signs
`χ_S(T)=±1` depend on `T`. So whether `M(g,e)` *can* be nonzero — the support
`{S : c_S≠0, S⊇{g,e}}` — is a property of `K_n` and the **positions** of `g,e`,
independent of orientation:

> **Necessary condition (verified, 0 violations across n=4,5 exhaustive and a
> 420k-pair n=6 sample):** `M(g,e)≠0 ⟹ g,e are jointly coverable by
> vertex-disjoint odd cycles in `K_n`.**

The covers come in two channels (the minimal `S`):
- **(A) one odd cycle** through both arcs;
- **(B) two vertex-disjoint odd cycles**, one per arc (needs ≥6 vertices — the
  `i₂` channel of HYP-830).

Consequences, all verified:
- **Sharing a vertex.** The 3-cycle on the arcs' three vertices covers both ⟹
  support nonempty for every `n≥3`. Empirically `M≠0` for **every** tournament
  (`P=1` at `n=4,5,6`): the signed sum *never* cancels for sharing arcs.
- **Disjoint arcs.** Support is nonempty **iff `n≥5`** (the 5-cycle `a-b-x-c-d`).
  So at **`n=4` the support is empty ⟹ `M≡0`** for all disjoint pairs. This is
  **exactly HYP-1165** ("disjoint-arc `Δ₂=0` at `n=4`, fails `n≥5`") — previously
  a bare computational fact, now explained as *odd-cycle coverability*: two
  opposite edges of `K₄` lie only on the **even** 4-cycle, and two disjoint
  3-cycles need 6 vertices, so neither channel fits.
- For `n≥5` the support is nonempty but cancellation makes `M≠0` only fractionally:
  `P(M≠0 | disjoint) = 0.52 → 0.69` as `n: 5→6`. So s685's "common odd cycle of
  `T`" was the *realized* `n≤5` picture (where channel B can't fit and a realized
  directed odd cycle is what survives the sign sum); the invariant form is
  orientation-free coverability.

The clean takeaway: the **delta-propagation pattern is the OCF arc-pair odd-cycle
incidence** — a single combinatorial object (which arc pairs can sit together in a
disjoint-odd-cycle cover of `K_n`), modulated by the orientation only through
signs. HYP-2268's "interaction support complete at `n=5`" is the `n=5` slice of
this; the `n=4` *incompleteness* (disjoint pairs never interact) is the visible
edge of the coverability law.

## 2. `A(T)` strictly refines `H`

`H` is famously a coarse invariant — many non-isomorphic tournaments share an
`H`-value. The second-order skeleton does better. **Exhaustively at `n=5`** (all
1024 labeled tournaments):

| invariant | # classes |
|---|---|
| `H` | 7 |
| `(H, degree-sequence of A(T))` | 10 |

`A(T)`'s degree sequence **splits** three `H`-fibers: `H=3` (80\|40), `H=9`
(120\|120), `H=15` (40\|24). So the interaction skeleton sees structure the count
is blind to: `H` is the *potential*, `δ` its *gradient* (HYP-2268), and `A(T)` the
*Hessian support* — and the Hessian carries strictly more information than the
potential. This is the precise sense in which the user's intuition holds: the
arc-flip delta pattern **is** a finer key to tournament structure than `H`.

## Where this sits

- **HYP-2268** (delta field = gradient/Hessian, Walsh/OCF) — this is its support
  theorem + an invariant-strength result.
- **HYP-1165** (`Δ₂=0` disjoint at `n=4`) — now a corollary of coverability.
- **HYP-830 / THM-170** (the `2·δc + 4·δi₂` channel formula) — channel (B) is the
  `i₂` term; the coverability picture unifies channels (A) and (B) as "which
  disjoint-odd-cycle covers contain both arcs."
- **THM-081** (OCF Walsh expansion) — the source of orientation-freeness.

## Next

1. **Realization (the sufficiency gap).** Support-nonempty ⇏ `M≠0` (cancellation).
   When does the signed sum `Σ_{S⊇{g,e}} c_S χ_S(T)` vanish? For sharing arcs it
   never does (`P=1`) — is there a parity/sign argument that the single-3-cycle
   term dominates? Pin the exact realized-interaction law for disjoint arcs at
   `n≥5` (the `0.52, 0.69, …` sequence).
2. **Is `A(T)` (as a labeled graph, not just degree sequence) a complete
   invariant up to isomorphism?** Test at `n=5,6` whether `A(T)≅A(T')` ⟹ `T≅T'`,
   or find the first collision — locating exactly how much of `T` the second-order
   skeleton remembers.
3. **The `{7,21}` holes on `A(T)`.** HYP-2268 found the holes polarize the
   *gradient*; do they constrain the *interaction graph* `A(T)` (e.g. forbid
   certain degree sequences), since the achievable `(H, A-degseq)` classes must
   route around `H∈{7,21}`?
4. **`A(T)` for LRC / distance graphs.** The arc-interaction = second-order Fourier
   (degree-2 Walsh). Transcribe to the distance-Cayley setting (HYP-2265): the
   degree-2 part of the connection-set transform — does an "interaction graph" of
   the LRC arcs refine the loneliness bound the way `A(T)` refines `H`?
