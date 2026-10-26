# The Fermat five, Pépin reworked for Collatz, and the Fermat line

2026-10-25, klein (worktree `codex/session-tournament-twist-unify`, fourth directive). **PROVED**
elementary statements with proofs in the notes; **FINITE-EXACT** receipts (`fermat-normal.json`,
`merge-normal.json`, byte-identical under `-O`); **CITED / inherited** material credited by file;
the completeness of `{3,5,17,257,65537}` is **OPEN** and no proof is claimed. Theorems F, M, P were
independently audited on 2026-10-26 by kind-pasteur-S5 (confirmed; see THM-4793's `audit` field and
`../fermat_landing_points_pepin_20261026.md`, whose addendum shows the landing points carry Pépin via
`L' = L^2 + (L-1)^2`); Theorem M* (MERGE_TREE.md) is not yet audited.

| Note | Content | Status |
|---|---|---|
| [PEPIN_COLLATZ.md](PEPIN_COLLATZ.md) | Pépin as a cycle-type statement for `x -> 3x+1` on `Z/F_n` (Theorem P); the Fermat line under Collatz: `2^{n-1}-1` valuation-two steps then one of valuation three, landing on `(3^{2^{n-1}}+1)/2` (Theorem F, dual of THM-4556's Mersenne line); consecutive `2^m+1` orbits merge by the `4x+1` rule at exact depth `L+2` or `L+3` (Theorem M, twin of THM-4556 (ii), an explicit THM-4581 pair-chain instance); the Pépin–Collatz dictionary with its loss; seven reformulations of completeness (affine permutation, Sierpiński rows, Paley family, Bernoulli denominators, rigid circulants, near-transitive path counts, heuristic density) and the hostile against a "tricky" congruence route | PROVED / EXACT / HEURISTIC / OPEN |
| [MERGE_TREE.md](MERGE_TREE.md) | The `L ≡ 1 (mod 4)` case (added 2026-10-26): the pair is THM-4581's chain from `(0,-1)` at the consecutive integers `3^{L+1}+1, 3^{L+1}+2`; `d(L) >= L+5` (proved), `d = L+5` iff `L ≡ 81 (mod 256)`, a 2-adic tree of depths with its first rows, no uniform bound; chain enumeration over `2^{18}` classes, 0 mismatches against 375 direct orbits | PROVED / FINITE-EXACT |
| [FIVE_LADDER.md](FIVE_LADDER.md) | `{0..4} -> {1,2,4,8,16} -> {3,5,17,257,65537}` as one tower (index, square-root tower height, prime), with a tournament floor `F_i = H(NT_{2^i+2})`; the odd constructible orders `3,5,15,17,…` as Sierpiński rows and as the orders with free multiplier action; Fermat Paley graphs (`P(5) = C_5` with golden spectrum, `P(17)` the `R(4,4)` witness); the repository's other fives typed as coincidences of cardinality | EXACT / inherited |

Not new: the Pascal-mod-2 identity (`gilbreath_fermat_platonic_20260926`), the valuation-two run count
(THM-4601), the Mersenne line (THM-4556/4605), the Fermat product and path-count identities
(THM-4748/4749), the rigidity theorems (THM-871/4644). New: Theorems F and M as a named family with
exact depths, Theorem P's affine phrasing, the ladder's tournament floor, and the link from THM-4756's
mod-2 forward-step row to the Fermat products.

## Reproduction

```
python3 check_fermat_line.py > fermat-normal.json      # SymPy for n_order/isprime; numpy for Paley spectra; ~2 min
python3 check_merge_points.py > merge-normal.json       # stdlib; asserts Theorem M's exact depths for L = 2..19
python3 check_merge_tree.py 20 > tree-normal.json       # stdlib; THM-4581 chain over all classes mod 2^18 (~3 min)
python3 check_odd_L_depth.py 1500 > depth-normal.json   # stdlib; direct orbits, controls, periodicity test (~1 min)
python3 check_pair_chain.py > chain-normal.json         # stdlib; chain states at 2L+1..2L+3 for L <= 33
```
