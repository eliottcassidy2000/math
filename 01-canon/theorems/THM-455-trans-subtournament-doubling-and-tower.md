# THM-455 — Transitive Subtournaments under Skew-Sylvester Doubling: the +1 law (with one exception), and the Mersenne tower's extremality window 7–31

**Status:** (1) PROVED. (2) PROVED for transitive T; the universal +1 conjecture REFUTED with
exact counterexample. (3)–(4) COMPUTED and INDEPENDENTLY VERIFIED (source recursion ≡ sink
recursion on every instance; explicit witnesses validated; brute-force subset enumeration
agrees 18/18 at n≤5). Literature cross-check in progress (agent out; values below marked
where they depend on it).
**Source:** kind-pasteur-2026-06-09-S2 (T769, HYP-2356/2357; computed by session owner after
two branch-agent API failures).
**Related:** THM-447/448 (doubling, tower, labeled link recursion), THM-454, Erdős–Moser
(1964), mac-mini THM-453 Part A (Erdős 592's red target = order-faithful transitive
subtournament — this theorem is its finite shadow).

## (1) Lower bound: trans(D(T)) ≥ trans(T) + 1 (PROVED)

If S = {v₁ → ⋯ → v_t} is transitive with source v₁, then S ∪ {v₁'} is transitive in D(T):
v₁ → v₁' (twin), v₁' → v_j for j ≥ 2 (cross arcs follow T). ∎

## (2) The +1 law holds in 17/18 classes at n ≤ 5 — but is NOT universal

trans(D(T)) = trans(T) + 1 for 17 of the 18 iso classes n=3..5 (verified by fast recursion AND
brute force). For transitive T it is exact (the interleaving obstruction: primed elements must
appear in reversed T-order while cross arcs pin each v_j' to the slot after v_j — two primes
cannot straddle an unprimed element, so a chain-comparable configuration gains at most 1).

**EXCEPTION (refutes universality):** n=5 idx10, scores (1,2,2,3,2), arcs
{(0,4),(1,0),(1,3),(2,0),(2,1),(3,0),(3,2),(3,4),(4,1),(4,2)}: trans(T) = 3 but
trans(D(T)) = 5 (+2). Witness TT5 in D(T): the chain 1', 3, 2', 0, 4' — primed and unprimed
ALTERNATE, using cross arcs in both directions. The interleaving obstruction binds only
chain-comparable configurations; alternating mixed chains evade it. **n=6 census (ALL 32768 labeled tournaments): delta = trans(D)−trans(T) distribution
{+1: 27248, +2: 5520} — delta NEVER exceeds 2.** The +2 rate grows with n (1/18 classes at
n≤5 → 16.8% labeled at n=6), concentrated on (t,t+2) = (4,6) and (3,5). This upgrades the
trivial bound trans(D(T)) ≤ 2·trans(T) to the SANDWICH CONJECTURE
trans(T)+1 ≤ trans(D(T)) ≤ trans(T)+2 (HYP-2360).

## (3) The tower values (all independently verified, witnesses explicit)

```
           T7    T15    T31    T63          controls (10 random each)
trans:      3      5      7     11          n=15: 7,7,7,7,7,7,7,7,8,8
                                            n=31: 9,10,10,10,10,10,10,10,10,11
Paley_31:                7                  n=63: 11,12,12,12,12,12,12,12,13,13
```
Tower-step growth: +2, +2, +4. TT7 witness in T31: (8,25,22,4,11,29,18); TT11 witness in
T63: (23,56,32,2,26,61,37,6,30,57,33).

## (4) Reading (pending final literature confirmation of the extremal record values)

- **T7 = Paley T₇ is THE extremal TT₄-free tournament** (largest order 7, unique — classical).
- **T15 achieves trans = 5 = the pointwise minimum f(15)** (assuming st(5) = 13: every
  14+-tournament contains TT₅, and TT₆-free tournaments exist at 15 as subs of the order-27
  extremal). The tower is f-extremal at order 15.
- **T31 = 7 = trans(Paley_31)**: the tower MATCHES Paley at 31 despite being non-isomorphic
  to it (THM-448c) — both are TT₈-free with TT₇, presumably one above the pointwise minimum
  (f(31) = 6 if a TT₇-free 31-tournament exists, e.g. inside the classical composition
  constructions — literature check pending).
- **T63 = 11 decouples**: at the low end of the random range (11–13) but far from extremal.
  The tower's Erdős–Moser extremality window CLOSES between 31 and 63 — mirroring how the
  tower leaves Paley after 7 (THM-448) and leaves Sylvester-equivalence at order 16 (THM-451).
  Every lens (isomorphism, Hadamard class, now extremal combinatorics) sees the tower stay
  classical for a few levels and then go its own way; remarkably its |Aut| = F₂₁ stays frozen
  through all of it.
- The tower remains BELOW the random median at every computed level (3<7-8, 5<7-8, 7<10,
  11<12) — the doubling structure suppresses transitive subtournaments at every scale, just
  not extremally at 63.

## (5) Bridge to Erdős 592 (mac-mini's thread — coordination, not duplication)

THM-453 Part A reframes ω^β → (ω^β, 3)² via order-faithful transitive subtournaments. This
theorem is the finite-quantitative shadow: explicit tournaments with controlled trans give the
finite obstruction landscape that 592's witness colorings must navigate. The doubling +1 law
(and its alternating-chain exception) is exactly the kind of "pattern-evading mixed chain"
mac-mini's E.3 found at ω³ (witnesses must use gap magnitudes, not pure patterns).

## Scripts

`trans_tower_erdos_moser_kps2.py`, `verify_trans_tower_kps2.py`,
`trans_doubling_n6_kps2.out` (n=6 census, complete) (+ .out files).
