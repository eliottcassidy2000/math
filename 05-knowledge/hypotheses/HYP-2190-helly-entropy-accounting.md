---
id: HYP-2190
title: Helly-entropy accounting — the covering-depth distribution as the master ledger for LRC
status: OPEN
source: opus-2026-06-03-S620
depends_on:
  - THM-410   # depth conservation + moment-sieve identity (the rigorous backbone)
related:
  - HYP-2130  # rigidity = orbit-type (perspectives / observer-coupling) — the LRC key
  - HYP-2140  # the 2-adic seam / apex wears every hat
  - THM-404   # doubling-rigidity dichotomy
  - THM-401   # pair-sum sieve (the order-2 layer)
  - lrc_tree_entropy_attack_s543  # the OTHER entropy (time-spread H-matrix entropy) — contrast
---

# HYP-2190 — Helly-entropy accounting

## The frame (what to extend to everything)

A loneliness certificate is a point of the clock avoiding every forbidden arc, so **LRC is a
circular-arc covering problem**, and its master object is the **covering-depth distribution**
`{p_k}` (THM-410). Three independent ledgers act on this one distribution:

| ledger | quantity | status | meaning |
|---|---|---|---|
| **conservation** | mean `= S_1 = 2n/(n+1) < 2` | THM-410, proved | the depth *charge* is arithmetic-independent |
| **Helly** (combinatorial) | the nerve of `{B_i}` has bounded order | partly classical | circular arcs ⇒ Helly number ≤ 3 ⇒ low-order sieve data is structurally complete |
| **entropy** (measure) | `H_depth = -Σ p_k log p_k` | empirical here | spread of the charge; `p_0` is the lonely measure |

**LRC, re-read:** *can a nonnegative-integer depth field with mean pinned just below 2 be forced
to vanish nowhere (`p_0 = 0` a.e.)?* The conservation law says the charge is small; the question
is purely how it is **distributed**, and that is what entropy measures.

## The conjectural content

**(H1) Entropy–tightness dichotomy.** Across speed sets at gap `1/(n+1)`:
generic sets sit at the **independence baseline** `p_0 ≈ (1-2δ)^n → e^{-2} ≈ 0.135` with **high**
`H_depth`; tight extremals collapse `p_0 → 0` with the **lowest** `H_depth`. Verified n=4..7:
AP `H_depth ≈ 1.0`, generic `≈ 1.5`, geometric (most spread) `≈ 1.7`. Conjecture: `p_0` is a
monotone increasing functional of `H_depth` along the resonance axis, and
**`p_0 = 0 ⟺ H_depth` is minimized over primitive sets** at the given `n`.

**(H2) The collapse family is larger than the AP-orbit.** The primitive sets with `p_0 = 0`
("measure-tight") are not just the AP. Computed:
- n=3: `{1,2,3}` only;
- n=4: `{1,2,3,4}` and `{1,3,4,7}` (note `7 = 3+4`);
- n=5: `{1,2,3,4,5}` and `{1,3,4,5,9}` (note `9 = 4+5`).
Conjecture: the collapse sets are exactly the AP together with **near-AP additive chains** — an
AP with its top element replaced by the sum of the two below it (a Fibonacci-type fold). Pin the
characterization and count them.

**(H3) Helly dimension 1 ⇒ order-2 governance.** Forbidden arcs live on a 1-manifold, so the
covering combinatorics has Helly number ≤ 3 and **fractional-Helly dimension 1**: pairwise overlap
data `S_2` should already pin the regime. This is the measure-theoretic home of the **pair-sum
sieve** (THM-401, modulus `2n-1`): order-2 is the fundamental sieve order because the geometry is
1-dimensional. Conjecture: a fixed-`S_1`, fixed-`S_2`, bounded-Helly-order entropy bound forces
`p_0 > 0` off the collapse family — a quantitative measure-LRC.

## How it unifies the thread ("extend to everything")

Every prior LRC invariant is a **functional of the single depth distribution `{p_k}`**:
- **sieve density** `ρ` (S561/HYP-2065) `=` the alternating moment sum `Σ(-1)^k S_k = p_0` (THM-410);
- **certificate-sheaf `H^0`** (S579) `=` the combinatorial type of the witness set `{depth=0}`;
- **corrector arity / single-vs-multi-sieve** (HYP-2075) `=` the order at which the alternating
  sum must be truncated to stay positive (Bonferroni depth);
- **rigidity = orbit-type** (HYP-2130/2135/2140) `=` *where* the depth charge concentrates: the
  apex `t = 1/2` (n even) is where many arcs pile onto one point, the entropy-minimizing pin; the
  observer-coupled (rigid) case localizes the charge to one orbit (single corrector), the
  observer-blind (symmetric) case spreads it across an automorphism orbit (multi-sieve). The
  **2-adic seam** (HYP-2140) is the entropy-minimum of the depth charge.

So the rigidity framework says *which symmetry* controls LRC; Helly-entropy accounting gives the
**quantitative ledger** those symmetries are spending against — a conserved charge `2n/(n+1)`,
distributed under a Helly (order) constraint, with entropy as the dial between safe (generic,
high-entropy) and tight (extremal, low-entropy).

## Contrast (avoid the collision)

This is **not** the S543 "H-entropy" (Shannon spread of the H-matrix over the runner walk). That
entropy is *high* for the AP (the AP spreads in time). The **depth entropy** here is *low* for the
AP (its covering charge piles up, never reaching 0). The two entropies are dual fingerprints of
the same extremal — a duality worth its own note.

## Next (small, chainable — see cluster-session-model)

1. Characterize/count the (H2) collapse family for n ≤ 8; test the "near-AP additive-chain"
   guess and whether each is a single ⟨×2⟩/affine orbit.
2. Try to *prove* (H3): a fixed-`(S_1,S_2)` + Helly-order-3 entropy lower bound on `p_0` off the
   collapse family (would be a genuine quantitative LRC step on the easy side).
3. Locate the depth-charge concentration point for even `n` and confirm it is the apex `1/2`
   (ties THM-404 / HYP-2140 to `H_depth` minimization).
