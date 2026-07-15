---
id: THM-866
title: AXIS LEVEL NO-HOLES COMPLETENESS — the populated x-levels are EXACTLY the full step-8 progression from the (near-)regular floor to the transitive ceiling (n³−n)/3; proved by the F3 exchange walk (tie-splitting flips), which also yields the canonical distance-to-transitivity k(T) = ((n³−n)/3 − x)/8
status: PROVED (all n; 5-step proof below) + VERIFIED (levels exact n = 4..14; walk mechanics on 800 random tournaments n = 8..11; floors by explicit circulants n = 3..12)
source: mac-mini-2026-07-15-S109; owner directive 2026-07-15 ("prove the no-holes completeness via the F3 exchange walk"); closes the opus-S316 handoff ("@klein the level-law no-holes completeness wants the F3 exchange-walk proof")
depends_on:
  - THM-855 F3 (per-flip drop law Δx = 4(d_l − d_w) + 8)
related:
  - opus-S316 probe (i) (the mod-8 congruence, two-line parity proofs; exact census n ≤ 8 — both recovered here)
  - THM-790 (the mod-16 selection rule lives on the FLIP layer, not the level lattice)
  - THM-801/855 (the axis as the metagraph's principal coordinate)
  - kind-pasteur cont.21 (I) [filed as THM-854, number collides with opus-S312's THM-854 — kps to renumber] and opus THM-867 — the other two proofs, see banner
script: 04-computation/axis_level_completeness_thm866_macmini_S109.py -> 05-knowledge/results/axis_level_completeness_thm866_macmini_S109.out
---

# THM-866 — axis level no-holes completeness

> **THE TRIPLE CONVERGENCE (2026-07-15, merge point).** Three agents proved this theorem
> independently on the same day, from the same owner directive, at three different levels
> of the same walk:
> - **kind-pasteur cont.21 (I)** (first push): the c₃-shadow — equal-degree pair exists
>   iff non-transitive, THM-833 gives Δc₃ = −1 per tie-split, and x is affine in c₃, so
>   the descent realizes every c₃ level (filed under THM-854, renumber pending).
> - **this file (mac-mini S109)**: the tournament-level walk — the flip happens inside an
>   actual tournament, so Landau-realizability is automatic; floors by explicit
>   circulants; census verified n ≤ 14; the k(T) distance corollary below.
> - **opus THM-867** (written during a network outage; renumbered from 865): the
>   score-sequence-level walk — the LEFTMOST tie-split is always Landau-legal (tight
>   prefix ⟹ contradiction at i+1), ties at 0/n−1 impossible; walk == census referee 7/7.
> Per opus's proposal, THM-866 carries the statement; THM-867 carries the leftmost-tie
> determinization; kps's file carries the c₃-descent form. Triple-certified.

**Theorem.** Let x(T) = Σ_v d_v², d_v = 2s_v − (n−1) (centered doubled scores). Then

> { x(T) : T a tournament on n vertices } = { x_floor, x_floor+8, x_floor+16, …, (n³−n)/3 },
> with x_floor = 0 for odd n and x_floor = n for even n.

No holes: every level in the step-8 progression is populated, and nothing else is.

## Proof

**(1) Parity frame.** d_v ≡ n−1 (mod 2) for every v, and Σ_v d_v = 2C(n,2) − n(n−1) = 0.
So for odd n all d_v are even; for even n all d_v are odd.

**(2) Floor.** Odd n: x = Σd_v² ≥ 0, attained by any regular tournament (circulant
i → i+1, …, i+(n−1)/2 mod n). Even n: each d_v odd forces d_v² ≥ 1, so x ≥ n; attained by
the near-regular circulant (i → i+1, …, i+n/2−1, plus i → i+n/2 for i < n/2), scores
n/2 copies of (n−2)/2 and of n/2. Both constructions machine-checked n = 3..12.

**(3) The F3 step.** Flipping one arc w → l sends d_w ↦ d_w − 2, d_l ↦ d_l + 2, so
Δx = 4(d_l − d_w) + 8 (THM-855 F3). **If d_w = d_l — a tied pair — then Δx = +8 exactly,
whichever direction the arc pointed.** Tie-splitting flips are the exact +8 moves.

**(4) A tie exists off the transitive.** If all n scores are distinct they are n distinct
values in {0,…,n−1}, hence exactly {0,1,…,n−1}; the score-(n−1) vertex beats all, and
downward induction on the remainder gives T transitive. Contrapositive: T not transitive
⟹ two vertices are tied ⟹ the arc between them is a +8 flip.

**(5) The exchange walk.** From any T, repeatedly flip an arc between two tied vertices.
By (3) x strictly increases by exactly 8 per step; tournaments on n vertices are finitely
many, so the walk halts; by (4) it halts only at the transitive order, where
x = Σ_{j=0}^{n−1} (2j−(n−1))² = (n³−n)/3. Hence **every** tournament satisfies
x(T) = (n³−n)/3 − 8k for some integer k ≥ 0 — which proves at once:

- the **ceiling bound** x ≤ (n³−n)/3 with equality iff transitive;
- the **mod-8 congruence** x ≡ 0 (odd n) / n (even n) (mod 8), since 24 | n³−n and
  24 | n³−4n give (n³−n)/3 ≡ 0 resp. ≡ n (mod 8) — opus-S316's two-line parity proofs
  recovered from the walk alone;
- **completeness**: walking from the floor tournament of (2), every level
  x_floor + 8j up to the ceiling is visited, hence populated. ∎

## Corollary — the canonical distance to transitivity

k(T) := ((n³−n)/3 − x(T))/8 is a nonnegative integer for every T, and the tie-splitting
walk is a canonical **transitivization in exactly k(T) moves**: every tournament sits
exactly k(T) tie-splits below the total order. (Sports reading: standings with a tie can
always absorb one more upset among tied teams, raising the spread Σd² by exactly 8;
standings with no ties are already a total order.) Level count: (x_ceil − x_floor)/8 + 1,
= (n³−n)/24 + 1 at odd n — e.g. 6, 15, 31, 56, 92 at n = 5, 7, 9, 11, 13 — matching the
exact censuses.

## Verification record (script above; all checks PASSED)

- Populated levels == predicted progression, **exact, n = 4..14** (Landau enumeration;
  extends the opus-S316 census n ≤ 8).
- Walk mechanics: 200 random tournaments at each n = 8..11 — every step Δx = +8, every
  walk terminates at the transitive with x = (n³−n)/3 and scores {0,…,n−1}.
- distinct-scores ⟹ transitive: exhaustive n ≤ 6 (0 violations in all 2^C(n,2)).
- Ceiling residue arithmetic: n ≤ 60, 0 violations.

## Notes

- The name "F3 exchange walk" is now literal: THM-855 F3's margin structure singles out
  margin-0 (tied) flips as the exact +8 steps; the walk is those flips iterated.
- THM-790's mod-16 alternation is a constraint on the d = m FLIP layer connecting levels,
  not on the level lattice itself (both mod-16 classes of levels are fully populated) —
  as recorded in opus-S316.
