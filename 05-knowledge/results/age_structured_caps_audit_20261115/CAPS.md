# Age-structured caps: every profile coordinate is at most 30·H, and why linear cones could not see it

Session klein-2026-11-15-age. The "age-structured potential" that THM-5243 proposed for the lemma `Z_ℓ ≤ 38·H` is
not needed as a proof device (AUDIT.md). What survives of the idea is a theorem about the transfer state itself: each of
the 119 profile coordinates is bounded by an explicit multiple of `H`, the multiple depending only on the ages of the
fragment heads and tails recorded in the profile, and the proof is the tail deletion of SHARP_GROWTH applied fragment by
fragment. Tools: `scripts/profile_caps.py` (exact fractions; receipts `profile_caps_frag1to3.out`, `profile_caps_frag4.out`),
`scripts/hratio5.c` (exhaustive maxima over all words of length ≤ 9 and ≤ 10, with history keys; receipts `hratio5_*`),
`scripts/cone_test.py` (history-refined cone tests; receipts `cone_test_K*.out`). The ordered state list with ages is
`data/states3.txt` (regenerated from `ribbon.py`'s `transfer_system(3)` and checked against the exported matrices).

## 1. Statement

Profile states (THM-4880, `ribbon.py`): a sorted tuple of fragments `(h, t)`; `h ∈ {0, 1, 2}` is the age of the
fragment's head (age `a` = the window vertex `n−1−a`), `h = −1` a head that left the window unfed (at most one such
fragment; it must become the global start); `t ∈ {0, 1, 2}` the age of the tail, `t = −1` a tail that left the window
(free: every later arc into it is forward). `x_s(T)` = number of path covers of `T` with profile `s`; `H(T) = Σ_{complete s} x_s`.

**Theorem (PROVED).** For every span-3 ribbon `T` and every profile `s`, `x_s(T) ≤ κ(s)·H(T)`, where `κ(s)` is the exact
rational computed by `profile_caps.py` (`data/kappa_table.txt`). In particular

| fragments in `s` | states | max `κ(s)` (proved cap) | max observed `x_s/H`, words of length ≤ 9 |
|---|---|---|---|
| 1 (complete) | 14 | 1 | 1 |
| 2 | 52 | 68/3 | 14 |
| 3 | 44 | 30 | 14 |
| 4 | 8 | 30 | 7 |

so **every raw coordinate of the span-3 transfer vector is at most `30·H(T)`**. (THM-5252 (B) proved the uniform bound
`K_3·H` with `K_3 = 26,425`; the age-structured count brings the constant to 30, with 14 observed.)
`κ(s) = 1` exactly for the 14 complete states and for the two-fragment singleton states `((0,0),(1,1))`,
`((0,0),(1,1),(2,2))`: 16 states are sharp (κ = observed), 102 are not (receipt `caps_compare.json`).

## 2. Proof

Fix `T` on old sites `0, …, n−1`, a profile `s` with `k ≥ 2` fragments, and a cover `(F_1, …, F_k)` with that profile.
Choose the base `F_1` to be the fragment with the dead head if there is one (otherwise any fragment); every other head
`b_i = head(F_i)` is then a window vertex whose identity the profile fixes (age `h_i` ↦ site `n−1−h_i`). Fix also an
order `F_2, …, F_k` of the remaining fragments (the computation minimises over the admissible bases and orders; any fixed
rule gives a valid bound).

Merge in that order, as in SHARP_GROWTH: let `P` be the path built so far (`P = F_1` at the start). To merge `F_i`,
delete from `P` its maximal terminal block `S_i` of vertices beaten by `b_i` and append `F_i`; if `P − S_i` is nonempty
its last vertex beats `b_i`, so the result is a path of `T`. After the last merge `Q` is a Hamiltonian path of
`T − (S_2 ∪ … ∪ S_k)`, and the `S_i` are pairwise disjoint.

Recovery: given the `b_i` (fixed by the profile) and the sets `S_i`, the cover is recovered from `Q` and the ordered
blocks (Hamiltonian paths of the `T[S_i]`): `F_k` is the suffix of `Q` from `b_k`; the path before `F_k`, with the block
`S_{k−1}` re-appended, is the merged path of the first `k−1` fragments; iterate. Hence, for fixed `(S_2, …, S_k)`,

`#covers ≤ H(T − ∪S_i) · Π_i e_i(S_i) ≤ H(T) · Π_i e_i(S_i)/H(T[S_i])`

by the product inequality for disjoint supports (AUDIT.md §1 step 1), where `e_i(S_i) ≤ H(T[S_i])` counts the admissible
orders of the block. Admissibility: `S_i ⊆ N⁺(b_i)`; `S_i` avoids the heads and tails of the fragments merged later (they
are not yet in `P`); a nonempty `S_i` ends at the tail of `P`, which is `tail(F_{i−1})` — a window vertex fixed by the
profile, or, if that tail is free, some site left of the window that `b_i` beats (then at most three sites, all among the
last six). Only the last six old sites enter: a window vertex of age `a` beats at most `3 + a` vertices, all within
distance three of itself, and every site left of the six beats all of them. Summing the weights over all admissible
tuples `(S_2, …, S_k)` and maximising over the `2^12` orientations of the back arcs among the last six sites gives `κ(s)`.
Every recorded `x_s/H` over the 134,217,728 words of length 9 lies below `κ(s)` (no violation, receipt
`caps_compare.json`), as it must. ∎

Remarks. (i) For two fragments with one live head of age `a`, the crude count gives `κ ≤ 1 + 2^(3+a) ≤ 33`; the weighted
local maximisation gives the exact table (largest 68/3 at `((−1,−1),(2,0))`, `((−1,−1),(2,1))`, `((−1,−1),(2,−1))`,
`((−1,−1),(2,2))`). (ii) The largest proved cap, 30, is attained by the local count for `((−1,−1),(1,−1),(2,−1))` and the
four-fragment `((−1,−1),(0,−1),(1,−1),(2,−1))`; their observed maxima are 7 and 1. (iii) The two-fragment sum is observed
at most `42·H`, the three-fragment sum at most `84·H`, the four-fragment sum at most `24·H` (words of length ≤ 9); summing
the per-state caps gives weaker but proved totals.

## 3. The observed caps are not all attained at finite length

Between lengths 7 and 9 three coordinate maxima grew: state 81 `((−1,−1),(0,−1),(2,2))` from 13/9 to 25/17, state 101
from 8/5 to 44/27, state 114 from 7/5 to 25/17; the other 116 maxima are stable from length 7 on (receipt
`hratio5_9_b*`; the length-10 run `hratio5_10_*` changes none of the 119 maxima: the three values rest at 25/17, 44/27, 25/17 at length 10, so the growth is not monotone in length, and whether they move again at length 11 is the pre-registered P6). So the suprema of `x_s/H` are not all finite-length
maxima (the sequences 7/5, 13/9, 25/17 are `(3·2^k+1)/(2^(k+1)+1) → 3/2`), and the proved caps `κ(81) = 24`,
`κ(101) = 12`, `κ(114) = 12` are far above the observed values. Pre-registered P4 (no growth between lengths 7
and 9; largest cap 14) is **refuted in its first clause and confirmed in its second**; P5 (κ sharp for most two-fragment
states) is **refuted**: only the 16 trivial states are sharp. The sharp per-profile constants remain a finite-exact
table (`data/kappa_table.txt`, columns observed and κ), not a theorem, except for the 16 states where the two coincide.

## 4. History-refined cones (pre-registered P1–P3)

`cone_test.py` tests whether the box cone `C(key) = {x ≥ 0 : f_j·x ≤ cap_j(key)·H(x)}` (features: the 119 coordinates,
the fragment-class sums, the eight `Z_ℓ` and the eight `H(Wℓ)`), with caps the exact maxima over words of length ≤ 9
grouped by the last `K` letters, is letter-invariant (`A_ℓ C(key) ⊆ C(key·ℓ)` by one linear programme per (key, ℓ,
feature)); when it is not, the violated caps are raised to the LP maxima and the test repeated (value iteration).

| `K` | keys | LPs per round | violations, round 1 | after iteration |
|---|---|---|---|---|
| 0 | 1 | 1,080 | 598 | stalls at 550 from round 3: the cone contains directions whose image has no complete state (`H(A_ℓ x) = 0`), so some relative caps are infinite |
| 1 | 8 | 7,272 | 5,328 | diverges: largest `Z` cap 38 → 49,325 and largest coordinate cap 14 → 12,718 in ten rounds, violations stay near 3,900 |
| 2 | 64 | 58,176 | 42,806 | diverges in three rounds: largest `Z` cap 38 → 8,375, largest coordinate cap 14 → 4,095, violations 42,806 → 28,931 |
| 3 | 512 | 8 × 58,176 (one round, eight key slices) | 341,961 (36,774–43,968 per slice of 64 keys), largest excess 198 | not iterated |

P1 (the K = 0 failure reproduces) **confirmed**; P2 (violations decrease with `K`) **refuted** — they grow with the number
of keys and do not vanish; P3 (no convergence within ten rounds for `K ≤ 2`) **confirmed**. No finite-history box cone
with data caps is invariant.

**Why.** The deletion proof bounds a coordinate by `H(T − S)·H(T[S]) ≤ H(T)`: a statement about Hamiltonian-path counts
of sub-tournaments, i.e. about the global orbit, which the 119-dimensional state does not carry (the state records the
window roles of fragments, not how the rest of `T` would count paths with vertices removed). A letter-invariant cone is a
property of the linear dynamics alone, and the linear dynamics admit the dead directions above. This is the conceptual
content of THM-5243 (3)'s two obstructions: span and positivity cannot express the product inequality. The age structure
does enter — through the admissible `S` (which sites a head of age `a` can beat, where the previous tail sits) — but as
a combinatorial injection, not as a potential.

## 5. Consequences and what remains

- Every finite suffix readout of length `m` is bounded by `39^m·H` on physical states (SHARP_GROWTH), and now every single
  coordinate by `30·H`. Both are uniform over the orbit and need no cone.
- The sharp per-profile constants (observed 14 at `((−1,−1),(2,0))` and `((−1,−1),(1,0),(2,−1))`) are a target for a
  sharper injection that uses fragment membership; the growing states 81, 101, 114 show that any proof must allow
  non-attained suprema.
- The lower-growth certificate on rank-constant stretches (HYP-9742 literal, HYP-9802 (iii)) is unaffected by all of
  this: upper caps do not bound growth from below.
