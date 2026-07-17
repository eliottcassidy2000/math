# THM-967 — The rational-ratio lock and its exact count (death-star-2026-07-17-S48)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCRationalLock.lean`,
standard trio ×5). Source: HYP-7240. Generalizes THM-966 (the case i′ = 1)
and consumes THM-961's mod transport and narrow count. Renumbered from the
colliding THM-963; THM-963 remains the scale-eight owner-nerve theorem.

## Statement

Speeds `g·i′, g·j′` with reduced ratio `i′ : j′` (via the relation
`j′·v_a = i′·v_b`), both failing the safe band at `p/q` with witnesses
`w_a, w_b`:

1. `rational_lock` (i′ + j′ ≤ 13): witnesses lock on the Bézout ray,
   `j′·w_a = i′·w_b` EXACTLY — the exact identity
   `j′(v_a p − w_a q) − i′(v_b p − w_b q) = (i′w_b − j′w_a)·q` plus the strict
   triangle `< (i′+j′)q/14`.
2. `rational_branch_bound` (i′ + j′ ≤ 27): `|j′w_a − i′w_b| ≤ 1` — at most
   THREE witness branches. Every pair of the canonical family {1,…,13} has
   i′ + j′ ≤ 25: the whole pair layer is ≤3-branch.
3. `witness_mod_bridge`: `∃s, 14M|x − sq| < q ⟺ 14M·r < q ∨ 14M(q−r) < q`
   (r = x mod q) — only s ∈ {⌊x/q⌋, ⌊x/q⌋+1} can win. (Also retro-closes the
   witness↔mod gap left open between THM-961's iff and count forms.)
4. `rational_pair_fail_iff` (gcd(i′,j′) = 1, locked): joint failure ⟺ the
   GCD-speed `g` fails the `14·max(i′,j′)`-narrow band.
5. `rational_pair_count`: at `gcd(g,q) = 1`,
   **N_ab = 2·⌊(q−1)/(14·max(i′,j′))⌋** — closed form.

## The discrete–continuous bridge (tandem with boxeph)

`N/(q−1) → 1/(7·max(i′,j′))` reproduces boxeph LEM-044's
`μ(D_k ∩ D_{k+1}) = 1/49 + r(6−r)/(49k(k+1))` exactly for k ≤ 6 (1/21 at k=2,
1/28 at k=3, …), and their zero-excess-iff-7∣k is precisely the locked/sparse
boundary (k + (k+1) ≤ 13 ⟺ k ≤ 6). Same session, both faces went kernel-pure:
this module (discrete, floors at every finite q — what the census consumes) and
boxeph's `LRCTreeHunter.lean` (continuous μ + tree-hunter). Recon:
`rationallock_recon_deathstar_S48.out` — lock 3958/3958, branches 7469/7469,
exact count on all 49 locked pairs of {1..13} (461 cases), bridge table k=2..8.

## Ledger status after this theorem

Pair layer of B5 on the canonical family: 49/78 pairs closed-form (locked);
29 sparse pairs are ≤3-branch with the k = ±1 Bézout branches the named
remainder. Sparse-branch counts + the triple layer (chain iff already covers
divisor chains) are the next rungs toward the a-priori liveCount floor.
