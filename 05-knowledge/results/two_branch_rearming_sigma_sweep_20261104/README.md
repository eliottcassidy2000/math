# Rearming the two-branch tree at 12875 and 12877, and the full-digit C2 sweep as the unbounded target selector

klein-2026-11-04 (worktree `codex/session-spectrum-20261104`), continuing the nonzero-preparation lineage ([`nonzero_preparation_20261103`](../nonzero_preparation_20261103/README.md), [`multiport_defect_grammar_20261103/construction/RESULT.md`](../multiport_defect_grammar_20261103/construction/RESULT.md), [`algebra/JOINT_CLASSIFICATION.md`](../nonzero_preparation_20261103/algebra/JOINT_CLASSIFICATION.md), [`algebra/STRONG_COMPILER.md`](../nonzero_preparation_20261103/algebra/STRONG_COMPILER.md)) and klein's spectrum notes (THM-4860, THM-4861, HYP-9453). **Status: FINITE-EXACT rearming witnesses (the +2 branches at both leaves are closed; the +4 branches: not closed at order 15 within this session; digit-2 neighbours of the 12875 witness exist one root-arc flip away (H = 14925 and 15273) and two flips away (H = 12467), and the search is recorded as OPEN); PROVED laws (THM-4862: the one-exit extension law, the digit-equals-cycle-count identity with its Fibonacci bijection, near-join readiness); FINITE-EXACT census theorems (ready floors through order 10; the full-digit C2 sweep realises every odd number in [307, 26551] from strong sole-exit tournaments of order ≤ 10); typed verdicts on the two-branch grammar (its ready spectra R₁, R₂ are sparse with exact small gaps, so unbounded closure of the ±2/±4 digits alone is unlikely, while the full-digit sweep is the natural unbounded selector, HYP-9454). The global spectrum classification (Chan–Pak) remains OPEN.**

Directive (2026-11-04): "use the new strong compiler to rearm 12875 and 12877 and find an unbounded target-selection mechanism that closes both branches toward the global spectrum classification".

## 0. What the lineage means, in the terms used here

A tournament `T` with Hamilton count `h` and a vertex `a` whose only out-neighbour is the root `r` (a **sole exit**) has digit `σ = e_r(T − a)`, the number of Hamiltonian paths of `T − a` ending at `r` (inherited: σ = [P A_a] − h, and the positive joint guard holds iff `a` is a sole exit with `e_r(T − a) > 0`). The two C2 births (`x` with `x → a` and everything else `→ x`; `y` with every old vertex `→ y` and `y → x`) give a strong child with `H′ = 3h + 2σ` (inherited and reproduced here on the prototype: 5 ↦ 17). The two-branch grammar uses digits σ = 1 (+2) and σ = 2 (+4). The depth-2 leaves of the common-state tree from the 11-name H1429 source are, at order 15 with root 8: 12869 (sole exit 14, σ = 1), 12871 (sole exit 12, σ = 2), 12875 and 12877 (no sole-exit vertex with σ ∈ {1, 2}: "none: C2 image"). [`scripts/verify_leaves.py`](scripts/verify_leaves.py) recomputes all four from the packet's own receipt ([`data/direct_children_receipt.json`](data/direct_children_receipt.json)) and reproduces exactly these guards.

**Rearming** a leaf means exhibiting a strong tournament on the same 15 names with the same `H` that carries a sole-exit vertex of the wanted digit; the all-strong native transport theorem (STRONG_COMPILER, PROVED and audited in the lineage) then connects the actual leaf to that target with every intermediate strong and all names retained. So rearming reduces to **target existence**, which is what is computed here; the transport programs themselves are the lineage's compiler and are not regenerated.

## 1. Rearming 12875 and 12877 (FINITE-EXACT)

A constrained annealing walk over the 90 free arcs of a 15-vertex tournament with vertex 1 fixed as a sole exit into vertex 0, started at `E_14 + a` (the transitive order with its extreme arc reversed, `H = 12287`, digit 1), keeping strongness and the digit as hard constraints and minimising `|log(H/target)|` ([`scripts/ready_search.c`](scripts/ready_search.c), mode `csearch`), found:

| target | digit | witness (out-row bit encodings, vertex 1 → vertex 0 only) | check |
|---|---|---|---|
| 12875 | σ = 1 | `14296 1 11219 32758 2050 32727 18 82 2258 10706 31702 195 11222 2514 15319` | strong; `H = 12875`; `e_0(T − 1) = 1` |
| 12877 | σ = 1 | `20116 1 29490 29623 13090 12547 29631 29494 3 12578 31742 29694 8451 259 13106` | strong; `H = 12877`; `e_0(T − 1) = 1` |
| 12875 | σ = 2 | not found (closest digit-2 state at order 15: `H = 12879`) | OPEN |
| 12877 | σ = 2 | not found (closest: `H = 12879`) | OPEN |

(`scripts/ready_search.c check 15 <rows>` re-verifies each line.) Relabel vertex 0 as the old root 8 and vertex 1 as any old name: the strong compiler transports the actual 12875 and 12877 leaves to these targets with all 15 names retained, and the next C2 step produces the level-3 values `38627, 38629, 38633, 38635` as actual children of the tree. The level-3 gap-2 targets listed in the lineage's `surgery/TARGETS.md` (orders 14–15) are consistent with this.

## 2. Laws proved here (THM-4862; proofs in [`PROOFS.md`](PROOFS.md))

- **One-exit extension.** `h(T + x) = h(T) − e_a(T) + [P A_a](T)` for the first birth (`x → a`, all else `→ x`); hence `h(T + x) = 2h + σ − e_a`. Checked on 1,794 cases.
- **The digit is a cycle count.** `C(U + a_exit) = s_r(U)` and `C(U + a_entry) = e_r(U)`: the readiness digit of a sole-exit attachment equals the number of Hamiltonian cycles of the *entry* attachment at the same root. Consequently `(T, a) ↦ ((T − a)^op + a_exit, a)` is a bijection from σ-ready strong tournaments onto strong sole-exit tournaments with exactly σ Hamiltonian cycles, which is why the σ = 1-ready classes are counted by `F(2n − 6)` (1, 3, 8, 21, 55, 144, 377 at orders 4–10; Garey's unique-Hamiltonian-cycle count, J. Combin. Theory Ser. B 13 (1972), CITED via THM-4846) and the σ = 2 classes by 1, 5, 18, 60, 191, 589. The bijection does not preserve `H`.
- **Near-join readiness.** In THM-4860's packet construction `T → M_k → z` with `z → t₀` reversed, the last chain vertex is a sole exit into `z` with digit `h(T)` (k ≥ 2): every attained value occurs as a digit.
- **Sandwich law (verified orders 3–6).** `H(S(X)) = −h + 4K + 2K₃` with digit `h(X)`, in the slot responses of an arbitrary `X`.

## 3. Ready floors through order 10 (FINITE-EXACT) and the name-retention wall

Minimum `H` of a strong tournament with a sole-exit vertex of digit σ ([`data/floor_k.out`](data/), [`data/floor10.out`](data/floor10.out)), with the lineage's cited values at 14–16:

| order | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| σ = 1 | 5 | 9 | 17 | 33 | 61 | 115 | 219 | | | |
| σ = 2 | — | 9 | 23 | 51 | 93 | 177 | 341 | 4841 | 9545 | 18905 |

The σ = 2 floor grows by a factor 1.9–2.0 per vertex (341 → 4841 is 1.94 per vertex over four steps; 9545/4841 = 1.97; 18905/9545 = 1.98), i.e. `floor₂(n) ≈ 0.29 · 2^n` for `n ≥ 14`, while a tree node at depth `ℓ` has `H ≈ 1429 · 3^ℓ` on `11 + 2ℓ` names. The name-retaining all-strong grammar therefore needs `1429 · 3^ℓ ≥ 0.29 · 2^{11+2ℓ}`, i.e. `(3/4)^ℓ ≥ 0.42`, which fails from `ℓ = 4` (order 19, `H ≈ 115,883` against a floor near 150,000) if the doubling continues; depth 3 (order 17, `H ≈ 38,630` against about 37,000) is the last comfortable level. This is the arithmetic form of the lineage's order obstruction and name budget (`3 log₅ 3 − 2 ≈ 0.048` extra names per step): **an unbounded mechanism must drop names**, i.e. select targets of smaller order than the births produce — which is harmless for the integer spectrum, where only existence of a realiser matters.

## 4. The full-digit C2 sweep (FINITE-EXACT) and why it is the right selector

Over every strong tournament of orders 4–10 and every sole-exit vertex, the pairs `(H, σ)` were enumerated (`ready_search.c pairs`, 1.7 million pairs at order 10 alone; per-`H` digit sets in [`data/pairs10_digits.txt`](data/pairs10_digits.txt)). Since each pair realises `3H + 2σ` as the Hamilton count of a strong tournament two vertices larger:

> **Every odd number in [307, 26551] is `3H + 2σ` for a strong sole-exit tournament of order ≤ 10** ([`data/sweep_analysis.out`](data/sweep_analysis.out); with orders ≤ 9 the interval is [307, 5153]). The odd numbers below 307 that are not reached are exactly `1, 3, 5, 7, 9, 11, 13, 15, 19, 21, 23, 25, 27, 33, 37, 39, 43, 45, 47, 49, 55, 61, 63, 67, 79, 83, 103, 107, 121, 125, 127, 167, 203, 305`, and each has the explanation "its unique odd predecessor under `3h, 3h + 2, 3h + 4` lacks that digit" (e.g. `47 = 3·15 + 2` and the digit set of 15 is `{0, 3, 4}`; `63 = 3·21` and 21 has no realiser).

Digit sets are large: for `H ∈ [300, 1000)` a ready value carries on average 86 distinct digits with maximum digit about `0.29 H`; the lines `3H + 2σ = v` therefore meet the ready-pair cloud many times for every target. This is the **unbounded target-selection mechanism** proposed here: for a target `v`, select any ready pair `(H, σ)` on the line `3H + 2σ = v` (for instance the one of least order) and apply C2. It is unbounded in the sense of being defined by the same rule at every scale; its validity for all `v ≥ 307` is the sweep conjecture of HYP-9454, equivalent to a two-dimensional density statement for the ready pairs `(h(U) − e_r(U) + [P A_r](U), e_r(U))` over strong `(U, r)`.

**Why small digits are rare, and why the two-branch grammar lives in the near-transitive corner.** By THM-4862 (b) the digit of a sole-exit vertex is the number of Hamiltonian cycles of the entry attachment at the same root. A digit of 1 or 2 therefore asks for a tournament with one or two Hamiltonian cycles — Douglas's uniquely Hamiltonian structures and their two-cycle relatives, near-transitive by nature — while a generic tournament of order `n` has digits comparable to `H/n`. This is the structural reason for the ready floors (near-transitive tournaments have `H ≈ c·2^n`), for the sparsity of `R₁`, `R₂`, and for the fact that random order-11 realisers of 12875 and 12877 carry digits from 98 to 1648 but never 1 or 2 (band census, §4 below). The full-digit sweep uses the generic realisers; the two-branch grammar uses only the rare ones.

By contrast the two-branch grammar restricts the digit to `{1, 2}`, and the ready spectra `R₁ = {H : digit 1 available}`, `R₂` are sparse: below 321, where orders ≤ 10 are complete, `R₁` holds 108 of the 161 odd numbers and `R₂` 79, with the first gaps `25, 27, 35, 41, 43, …` and `25, 27, 31, 33, 35, 39, …`; the cumulative density of `R₁` falls (0.58, 0.44, 0.40 at 37, 109, 321) and the per-order density of `R₁ ∩ [floor, max]` falls from 1 to 0.22 by order 10. The lineage's reduction (PREPARATION_REDUCTION: unbounded weak-allowed two-branch closure ⟺ every ternary descendant lies in `R₁ ∩ R₂`) thus asks for a sparse set to contain an infinite tree; nothing here refutes it for the particular descendants of 1429, but the full-digit sweep does not need it. Depth-3 evidence at the tree's own scale: random strong tournaments of orders 11 and 12 with a sole-exit vertex were drawn (`ready_search.c sample`, two seeds each) and the pairs with `12000 ≤ H ≤ 14000` kept ([`data/samp_*.txt`](data/), [`band_sampling.py`](scripts/band_sampling.py), [output](data/band_sampling.out)): 43188 pairs, 1000 distinct H of the 1001 odd values in the band, mean digits per H 43.2; C2 children cover 1000 of 1000 odd numbers in [38001, 39999] (1.000); first uncovered: []; all eight level-3 tree values `38609, 38611, 38615, 38617, 38627, 38629, 38633, 38635` are reached. So at the tree's own depth-3 scale the selector works with targets of order 11–12, four to six names fewer than the births produce — exactly the order drop that §3 says an unbounded mechanism must make.

## 5. Verdicts

- Rearming at 12875 and 12877: the +2 branches are closed by explicit order-15 targets (FINITE-EXACT); the +4 branches: open at order 15 (digit-2 states within a few flips exist, e.g. `H = 12467, 14925, 15273` for the 12875 witness); in the weak-allowed language of PREPARATION_REDUCTION a smaller-order digit-2 witness padded by a transitive prefix would also do, and none is known either.
- Unbounded mechanism: the two-digit grammar is walled by the ready floors (name retention) and by the sparsity of `R₁`, `R₂`; the full-digit C2 sweep is defined at every scale and exact on [307, 26551] from orders ≤ 10 (THM-4862 (d)); its completeness is HYP-9454 (1), a density statement, not a theorem.
- Toward the global classification: the sweep replaces "closure of a binary tree" by "every odd number has a ready predecessor pair", which is one operator away from Chan–Pak; the predecessor structure of `3h, 3h + 2, 3h + 4` (every odd number has exactly one odd predecessor, lineage COVERAGE.md) makes digit-0 readiness, which is nearly universal (`R₀` misses only 47 and 67 among the odd numbers from 23 to 599), the free step, and the nonzero digits the scarce ones.

## 6. Reproduction

```bash
cd 05-knowledge/results/two_branch_rearming_sigma_sweep_20261104
cc -O2 -o scripts/ready_search scripts/ready_search.c -lm
python3 scripts/verify_leaves.py                                  # the four leaves from the lineage receipt
scripts/ready_search check 15 14296 1 11219 32758 2050 32727 18 82 2258 10706 31702 195 11222 2514 15319
for k in 4 5 6 7 8 9; do gentourng -q $k | scripts/ready_search floor $k; done   # ready floors (order 10: ~10 min)
for k in 4 5 6 7 8 9; do gentourng -q $k | scripts/ready_search pairs $k > data/pairs_$k.txt; done
python3 scripts/sweep_analysis.py; python3 scripts/one_exit_lemma.py; python3 scripts/lemma2_check.py; python3 scripts/sandwich.py
scripts/ready_search csearch 15 1 12875 1 3000000       # the +2 rearming search (minutes)
shasum -a 256 -c SHA256SUMS
```

Receipts for the Python scripts are byte-identical under `python3` and `python3 -O`.
