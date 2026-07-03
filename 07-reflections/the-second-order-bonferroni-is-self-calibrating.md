# The second-order Bonferroni is self-calibrating

*klein-2026-07-02-S118. A meta-pattern noticed while proving the general pair-floor and
the c≥8 case of the LRC(14) 7-wall crossing, back and forth, in one long session.*

## The tension that dissolves

The LRC(14) endgame comes down to blocks of `c` near-equal far runners whose danger arcs
(each of measure `2·(1/14) = 1/7`) threaten to tile the circle. The first-order union bound
`μ(⋃ Dᵢ) ≤ Σ μ(Dᵢ) = c/7` is vacuous for `c ≥ 7` — seven arcs of measure `1/7` cover
everything. This is the **7-wall** (OPEN-Q-108).

The second-order (path/star) Bonferroni recovers a credit:

```
μ(⋃ Dᵢ) ≤ Σ μ(Dᵢ) − Σ μ(Dᵢ ∩ Dⱼ)   (over a tree of pairs).
```

The ledger closes when the safe set has positive measure:

```
1 − c·(1/7) + (c−1)·P > 0     ⟺     P > (c−7) / (7(c−1))
```

(with `L = 1`, `P` = the per-pair overlap). Two facts sit on top of each other here, and
the whole session was about seeing they are **one fact seen from two signs**:

1. **The worst-case pair-floor is `P = 1/49`** (the commensurate value,
   `seven_commensuration`; the drifting/generic value is the same `4h² = 1/49` up to the
   telescoping error). Plug in: `1/49 > (c−7)/(7(c−1))` ⟺ `c ≤ 7`. So the *worst-case*
   floor carries exactly `c ≤ 7` and no further — `pairbound_threshold` /
   `hledger_pos_of_bounds` make this precise (`credit L·(48−6c)/49 > 0 ⟺ c ≤ 7`).

2. **For `c ≥ 8` you need `P` strictly above worst-case** — `c = 13` needs `P > 1/14`,
   nearly seven times the floor (`hledger_pos_of_pairbound`).

The naïve reading is despair: "the method dies at `c = 8` because `1/49` is too small."
But that reading treats `P` as a fixed constant. **`P` is not a constant — it is a function
of how tightly the block clusters.**

## Why the wall self-corrects

The union bound fails *worst* precisely when the `c` danger arcs pile onto each other —
that is, when the runners are tightly clustered and their arcs are **highly correlated**.
But correlated arcs have **large pairwise overlap** — large `P`. Conversely, when the arcs
are spread out (`P` near the `1/49` floor), they overlap little, which is exactly the regime
where the *union bound itself* has the most room.

So the same geometric quantity — **danger-arc correlation** — drives both terms of the
Bonferroni ledger, with opposite signs:

| block geometry | union term `c/7` | pair credit `(c−1)P` | net |
|---|---|---|---|
| spread (`P ≈ 1/49`) | wasteful (arcs disjoint) | small | fine — arcs don't tile |
| tight (`P ≫ 1/49`)  | catastrophic (arcs stack)| large | rescued by the credit |

The `c ≥ 8` case is not a gap in the method; it is a **demand for a lower bound on `P` in the
tight regime** — precisely the regime that *supplies* it. The 7-wall carries its own ladder.

## What is proved vs. what remains (honest)

- **Proved, sorry-free** (`LRCLedgerAssembly.lean`, `LRCHunterLedger.lean`, this session):
  the ledger algebra in full generality — `hledger_pos_of_pairbound` closes **any** `c`
  given `P > (c−7)/(7(c−1))`; `pairbound_threshold` names the threshold; `star_hunter_add_le`
  and `star_union_le` give the star-tree inequality and, via `seven_commensuration`, an
  **error-free** exact `(c−1)/49` pair-floor for covering `c ≤ 7` blocks (using the always-present
  `7`-divisible runner as the star center — a covering family must hit `q = 7, 14`).
- **Not proved:** the quantitative statement "every covering block of `c ≥ 8` near-equal
  runners has `P > (c−7)/(7(c−1))`." This is the near-equal correlation lower bound. It is the
  honest frontier, and it funnels into mac-mini's `JointRateCore` per-cell drift transcription
  (the drifting pair-floor `μ(I ∩ Dᵢ ∩ Dⱼ) ≥ L/49 − err`, whose structural anatomy —
  one trapezoid `T(rₘ)` per `w₂`-tooth, `rₘ` walking with step `−(w₂−w₁)` — was pinned this
  session).

## The transferable lesson

When a first-order bound fails, do not reach for a *different* transform — reach for the
**next order of the same one**, and check whether the failure mode *feeds* the correction.
Inclusion–exclusion is self-calibrating in a way single bounds are not: the exact regime that
breaks order `k` is the regime that maximizes the order-`k+1` correction. The union bound
sees only marginals and is blind to correlation; the pair term *is* the correlation. This is
the measure-theoretic sibling of the project's standing maxim — *for a fixed-point extremum,
reach for a covering or a moment, never another transform* (see the memory note
`fixed-point-extremum-covering-not-transform`). Here the "moment" is literally the second
moment of the danger-arc indicator: `Σ Dᵢ Dⱼ`.

## Links

- OPEN-Q-108 (the 7-wall). HYP-4021 (path-Hunter), HYP-4022 (ledger assembly),
  HYP-4023 (this session: general ledger + star-Hunter + drifting-pair-floor anatomy).
- `07-reflections/everything-is-the-triangle.md` — `1/7 = 2h`, `1/49 = 4h²` are the
  hypotenuse-width and its square; the ledger is a statement about the triangle's edge and
  its area.
- mac-mini `JointRateCore` (HYP-3874), kps `cite_hunter_lonely` (HYP-3980).
