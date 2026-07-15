---
id: THM-865
title: THE LEVEL COMPLETENESS THEOREM — the populated axis levels are exactly the full step-8 progression from the (near-)regular floor to the transitive ceiling (n³−n)/3, with residue x ≡ 0 (mod 8) at odd n and x ≡ n (mod 8) at even n; proof = the leftmost-tie-split walk: any tie splits with Δx = +8 exactly (F3 with zero margin), ties at scores 0 and n−1 are impossible (the mutual arc), and the LEFTMOST tie is always Landau-legal to split (tightness of its prefix forces the bottom-i to be {0,…,i−1} with s_{i+1} = i−1, violating Landau at i+1) — so every non-transitive score sequence climbs by exactly +8, and the walk from the floor visits EVERY level
status: PROVED (one-paragraph argument below) + REFEREED (n = 3..9: the walk from the floor visits exactly the census levels — walk == census, all splits legal, all steps +8, 7/7)
source: opus-2026-07-15-S317 (owner: prove the no-holes completeness via the F3 exchange walk); completes S316's level-law probe (HYP-6935)
depends_on:
  - THM-855 F3 (the per-flip drop law: the +8 tie-split is the zero-margin case)
related: [THM-790 (the d=m line layer's mod-16 selection rule — a FLIP-layer constraint, now cleanly separated from the level lattice), HYP-6935 (the residue laws, proved there)]
verification: 05-knowledge/results/level_completeness_vandermonde_opus_S317.out
---

# THM-865 — the level completeness theorem

**Statement.** Let L(n) = {Σ_v (2s_v − (n−1))² : s a score sequence on n
vertices}. Then L(n) is EXACTLY the arithmetic progression with step 8 from
x_min (= 0 for odd n, realized by the regular sequence; = n for even n, by
the near-regular sequence) to x_max = (n³−n)/3 (the transitive sequence).
In particular: odd n ⟹ every element ≡ 0 (mod 8); even n ⟹ ≡ n (mod 8);
and there are NO holes: |L(n)| = (x_max − x_min)/8 + 1.

**Proof.**
*(Residues.)* Odd n: d_v = 2e_v with Σe_v = 0 forces evenly many odd e_v, so
Σe_v² is even and x = 4Σe² ≡ 0 (mod 8). Even n: every d_v is odd, d² ≡ 1
(mod 8), so x ≡ n (mod 8).
*(The +8 move.)* Splitting a tie — replacing (a, a) by (a−1, a+1) in the
score multiset — changes x by 4(d_u − d_v) + 8 = 8 (THM-855 F3 with zero
margin). A tie at a = 0 or a = n−1 is impossible: the arc between the two
tied vertices gives one of them a win (resp. a loss). Hence any tie has
1 ≤ a ≤ n−2 and the split stays in range.
*(Legality of the leftmost split.)* Sort s₁ ≤ … ≤ s_n and let i be minimal
with s_i = s_{i+1} = a. The split lowers only the i-th prefix sum, by 1; it
is illegal only if Σ_{k≤i} s_k = C(i,2) (tight). Below the leftmost tie the
sequence is strictly increasing, so s_k ≥ s₁ + (k−1) ≥ k−1 for k ≤ i, and
tightness forces equality throughout: s_k = k−1 for k ≤ i, in particular
a = i−1. But then Σ_{k≤i+1} s_k = C(i,2) + (i−1) < C(i,2) + i = C(i+1,2),
violating Landau at i+1 — contradiction. So the leftmost split is always
legal.
*(Completeness.)* A score sequence has no tie iff its n scores are the n
distinct values 0,…,n−1, i.e. iff it is transitive. So from any non-transitive
sequence the walk takes a legal +8 step; starting at the floor it visits
x_min, x_min+8, …, terminating only at the transitive x_max. Every step-8
value in between is therefore realized. ∎

**Referee.** n = 3..9: the walk from the floor visits exactly the exhaustive
census levels (31 levels at n = 9), all splits Landau-legal, all steps +8.

**Reading.** The axis's level lattice is COMPLETE and combinatorially trivial
to walk; all the metagraph's subtlety lives in which levels the various FLIP
LAYERS connect (THM-790's mod-16 selection on d = m lines) and how classes
distribute over levels — the lattice itself has no arithmetic obstructions
beyond the two-line residue laws. The tie-split walk is the score-space
avatar of "majorization toward transitive": x is Schur-convex and the walk
is a maximal chain in the majorization order with constant x-increment.
