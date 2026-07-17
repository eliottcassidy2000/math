# THM-972 — The relation lock by coefficient weight and the mediant triple count (death-star-2026-07-17-S51)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCRelationLock.lean`,
standard trio ×6). Source: HYP-7260 (the owner's "projective 6-cycle and
triangular prisms" hint, decoded). Subsumes the S46–S50 lock series.

## The principle

**Witnesses inherit every light integer relation of the speeds.** If
`Σαᵢvᵢ = 0` with coefficient weight `Σ|αᵢ| ≤ 14` and every involved runner
fails the safe band at p/q, then `Σαᵢwᵢ = 0` exactly. Proof: the identity
`ΣαᵢXᵢ = −(Σαᵢwᵢ)q` plus the per-term integer gap `14|Xᵢ| ≤ q−1`:
`14|Σαw|·q ≤ (Σ|α|)(q−1) ≤ 14(q−1) < 14q`.

## Statement (Lean ×6)

1. `relation_lock` — the general Finset form (any index set, ℤ weights).
2. `relation_lock3` — the three-term workhorse (direct proof).
3. `sum_triple_lock` — sum-triples {a, b, a+b} have weight 3: ALWAYS locked,
   `w_c = w_a + w_b`, **even when both extreme pairs are sparse** (e.g.
   (5,6,11): pairs at weight 16, 17).
4. `rational_lock_weight14` — the corrected pair boundary: ratio relations
   have weight i′+j′, so pairs lock **through 14** (extends THM-967's ≤ 13;
   the weight-14 pairs of {1..13} — (1,13), (3,11), (5,9) — have EMPTY Bézout
   branches; the true branching regime starts at weight 15).
5. `mediant_triple_fail_iff` — the projective object: the mediant chain
   `(g·i′, g·j′, g·(i′+j′))` (Farey/Stern–Brocot line; the sum-lines of
   {1..7} form the Fano-type incidence whose Levi 6-cycles are
   line-triangles): all three fail ⟺ the gcd-speed fails the
   `14(i′+j′)`-narrow band.
6. `mediant_triple_count` — `N = 2⌊(q−1)/(14(i′+j′))⌋` at coprime moduli:
   **the triple layer's first exact rung.**

## Recon (`relationlock_recon_deathstar_S51.out`)

Relation lock 4471/4471; weight-14 pairs branch-empty 3/3; sum-triple lock
4135/4135 over all 36 sum-triples of {1..13}; mediant counts exact 9/9
(`N = 2⌊(q−1)/(14·c/g)⌋`).

## Consequences for the ledger

The S3 layer of {1..13} opens: the 36 sum-triples are all relation-locked;
those with c ≤ 13 have exact mediant counts. Non-sum triples decompose by
their relation lattice (weight ≤ 14 relations lock; the branch analysis
extends by the same Bézout normal-form pattern). The projective incidence of
shared pairs among sum-lines (the 6-cycles/prisms) organizes which triple
counts are independent — the next structural question.
