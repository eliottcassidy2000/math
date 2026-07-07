# THM-632 — Even-branch clearance: the canonical mediant family is not a gap member for even N

**Status:** VERIFIED (proof below; N=12 instance formalized sorry-free & kernel-pure in Lean, `LRCEvenBranchWitness.lean`)
**Author:** mac-mini-2026-07-06-S29
**Depends on:** opus-S118 bordered-AP construction (HYP-4506), HYP-4572 trichotomy (mac-mini-S28)
**Relevance:** the LRC(14) second gap (G) — explains the canonical construction's failure at N=12 by PARITY

## Statement

Let `N` be even, `Q = 3N − 1` (odd), and let `F(N) = {1,…,N}\{N−1} ∪ {3(N−1)}` be
opus-S118's canonical mediant family (N speeds). Then

> **M(F(N)) ≥ 3/(3N−1) > 2/(2N+1),**

so `F(N)` is **not** a member of the second gap `(1/(N+1), 2/(2N+1))` — its lonely
value lies strictly **above** the gap.

## Proof

Take the witness time `t = b/Q` with `b ≡ 3·2⁻¹ (mod Q)` (well-defined since `Q`
is odd). For any speed `v`, `‖v t‖ = dist(v·b mod Q, 0)/Q`, and

  `‖v t‖ < 3/Q  ⟺  v·b ≡ {0,±1,±2} (mod Q)  ⟺  3v ≡ {0,±2,±4} (mod Q)`

(multiplying by `2` and using `2·2⁻¹ ≡ 1`; note `v·b = 3v·2⁻¹`). Using the two
identities on `Q = 3N − 1`:

- `3⁻¹ ≡ N (mod Q)`   (since `3N = Q + 1 ≡ 1`),
- `3(N−1) = Q − 2 ≡ −2 (mod Q)`,

the four congruences `3v ≡ ±2, ±4` have the unique solutions

| `3v ≡` | `v ≡ (mod Q)` | value in `[1,Q−1]` |
|---|---|---|
| `2`  | `2·3⁻¹ = 2N` | `2N` |
| `4`  | `4·3⁻¹ = 4N ≡ N+1` | `N+1` |
| `−2` | `−2N ≡ N−1` | `N−1` |
| `−4` | `−4N ≡ 2N−2` | `2N−2` |

and `3v ≡ 0` has no solution with `0 < v < Q` (as `gcd(3,Q)=1`). The four collision
values are `{2N, N+1, N−1, 2N−2}`. Every speed of `F(N)` lies in
`{1,…,N} ∪ {3(N−1)}`, and of the four collisions the **only** one in that range is
`v = N−1` — which is exactly the element removed to form `F(N)`. (Indeed `2N, 2N−2 > N`
and `≠ 3(N−1)` for `N > 3`; `N+1 > N` and `≠ 3(N−1)` for `N > 2`.) Hence **no** speed
of `F(N)` collides: all `‖v t‖ ≥ 3/Q`, with equality at the binders `v = 2`
(`3·2 ≡ 6`, residue `3`) and `v = 3(N−1) ≡ −2` (residue `−3`). Therefore
`M(F(N)) ≥ margin(F(N), t) = 3/Q = 3/(3N−1)`.

Finally `3/(3N−1) > 2/(2N+1)` because `3(2N+1) = 6N+3 > 6N−2 = 2(3N−1)`. ∎

## Why N−1 is removed (and why parity is the operative arithmetic)

The removal of `N−1` is **forced**: it is the unique small speed pinned to residue
`−1` (distance 1) at the binding time — precisely the `3v ≡ −2` collision above.
The trichotomy (HYP-4572) shows the far element `3(N−1)` binds the *smallest
feasible* small speed `s ∈ {2,3,5}` at consecutive denominators `Q ∈ {3N−1, 3N,
3N+2}`. The speed-2 branch (this theorem) is feasible **iff `Q = 3N−1` is odd iff
`N` is even** (`2b ≡ 3` is unsolvable mod an even `Q`). So for even `N` the family
clears at `3/(3N−1)`, above the gap; the tight mediant `3/(3N+2)` is reached only
when both larger branches die, i.e. `N ≡ 1 mod 6`.

**LRC(14) consequence:** `N = 12` is even ⇒ `M(F(12)) ≥ 3/35 > 2/25`. The canonical
mediant construction fails to be a gap member at N=12 **by parity**, not because
`38 = 2·19` is composite. (The "38 = 2·19" is really `3·13 − 1`, the competing
denominator at `N = 13`, whose *evenness* is what lets `N = 13` reach the mediant.)

## Scope

This governs the **canonical** family `F(N)`. It does not by itself prove the full
crux (no 12-speed family whatsoever in the gap) — that remains opus-S118's fleet
sweep (empty at N=12). The reductive next step: show non-canonical gap-member
species obey the same parity gate at N=12.

## Verification

- `04-computation/lrc_mediant_trichotomy_macmini_S28.py` — trichotomy, N=5..30.
- `04-computation/lrc_even_branch_clearance_macmini_S29.py` — the clearance +
  collision-value audit (`inv3=N`, collisions `{2N,N+1,N−1,2N−2}`, N=6..28).
- `04-computation/lean/.../LRCEvenBranchWitness.lean` — N=12 instance, sorry-free,
  axioms `[propext, Classical.choice, Quot.sound]`: `F12_margin` (3/35 ≤ margin at
  19/35) and `F12_reach_above_gap` (2/25 < reach).

→ HYP-4572, HYP-4506/opus-S118, THM-622.
