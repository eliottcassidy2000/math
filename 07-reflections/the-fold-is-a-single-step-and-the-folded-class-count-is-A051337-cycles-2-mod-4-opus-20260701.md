# The fold is a single terminal step (no quarter-tiling), and the folded class-count is A051337 — governed by anti-automorphisms whose cycles are all ≡ 2 mod 4

*opus-2026-07-01-S20. The owner asked (1) whether the half-tiling model `H_n` itself folds to a quarter-tiling,
and (2) for a recursion rule for odd n mirroring the even-n folded-class-count rule. Both resolve cleanly —
and both **correct a claim from S19** (the "dyadic tower" and the "even rule").*

## Task 1 — H_n does NOT fold to a quarter-tiling; the fold is a single terminal step
The blue subgraph = `H_n` = the flip-metagraph of the **complement-fixed (self-complementary) tournaments**
(S19). The folding is driven by **R = complement** (S18), an involution of **order 2**. Its fixed set is
exactly the SC world, i.e. `H_n` itself. So *inside* `H_n` every tournament is already complement-fixed —
applying "fold by complement" a second time fixes **everything**, giving no proper sub-fold.

More sharply: a tiling model folds via its grid reflection, and (S18) the grid reflection **is** the
complement. `H_n` is not a staircase model (`D = ⌊(n−1)²/4⌋` is not `C(k−1,2)` in general), and its
tournaments are complement-trivial, so it carries **no non-trivial grid reflection** — hence no blue/black
2-coloring of its own lines, hence no fold. Verified n=4,5,6: *every* grid-symmetric (half-tiling) tournament
is self-complementary.

> **Correction to S19.** I speculated a "dyadic complement-folding tower" (blue → quarter → …). That is wrong.
> The complement is one involution; folding once exhausts it. **`M_n → H_n` is a single, terminal fold.** The
> SC world does not fold again — it is the fixed point of the only available involution.

(If one wanted a further reduction it would have to come from a *different* structure — e.g. the involutive
anti-automorphism of the SC tournaments, cf. Task 2 — not from iterating the complement.)

## Task 2 — the folded class-count is A051337, and the even rule was a coincidence
The folded class-count is `#SC(n)` = the number of self-complementary (self-converse) tournament classes.
Computing it by Burnside over anti-automorphisms — `#SC(n) = (1/n!) Σ_π #{T : π(T) = complement(T)}`, with the
inner count `= 2^{orbits}` when consistent, `0` otherwise — gives

```
#SC(n) = 2, 2, 8, 12, 88, 176, 2752, 8784, 279968, 1492288, 95458560, 872687552   (n = 3..14)
       = A051337  (self-converse tournaments),  matching the metagraph B+M for n≤7.
```

> **Correction to S19.** I claimed `#SC(2k) = A000568(2k−1)` (the "even rule"). It holds at **n=4** (2=A₃) and
> **n=6** (12=A₅) **only by coincidence** — it **breaks at n=8** (`#SC(8)=176 ≠ A₇=456`) and n=10
> (`8784 ≠ 191536`). The folded class-count is *not* the ordinary count one size down.

### The real (unified) rule, and the even/odd mirror
The inner Burnside count is nonzero **iff every cycle of the anti-automorphism π has length ≡ 2 (mod 4)** —
i.e. π is built from cycles of length `2, 6, 10, 14, …` (= 2 × odd) — **plus exactly one fixed point iff n is
odd**. Verified exhaustively for n=4..14: the contributing cycle types are *precisely* the partitions of `n`
(even) or `n−1` (odd) into parts ≡ 2 mod 4, appended with a single `1` for odd n. Consequences:

- **The even/odd mirror the owner asked for.** The odd-`n` rule **is** the even-`(n−1)` rule with one extra
  **fixed vertex** (the τ-fixed point): the contributing cycle types for `2k+1` are exactly those for `2k` with
  a `1` appended. So odd n mirrors even (n−1), differing only by the anti-automorphism's fixed point — not by a
  clean `A000568` shift.
- **The count of contributing types = `q(⌊n/2⌋)`** — partitions of `⌊n/2⌋` into distinct (equivalently odd)
  parts: `1,1,2,2,3,3,4,4,5` for n=4..14. (Divide each part-≡2-mod-4 by 2 → odd parts of `n/2`; Euler.)
- **The dominant term is the all-transpositions `τ`** (cycle type `2^{⌊n/2⌋}·1^{n mod 2}`), the standard
  fixed-point-free (odd: near-fixed-point-free) anti-automorphism; the correction terms swap blocks of six
  transpositions for 6-cycles (then 10-cycles, …), one "distinct-odd-part" family at a time.

So the "complicated recursion rule" is a **cycle-index over `2 mod 4` cycles** — the same rule for both
parities, with odd `n` = even `n−1` + a fixed point. This is why `#SC` grows like `A051337`, not like a shifted
tournament count.

## Why 2 mod 4 (the mechanism)
Under the complement anti-automorphism `π(T)=T^op`, walking a π-cycle of vertices reverses arcs each step, so
the orientation of each pair-orbit must **alternate**; consistency needs each pair-orbit even, and the pair
structure of a length-`ℓ` vertex-cycle closes up consistently exactly when `ℓ ≡ 2 (mod 4)` (lengths `≡ 0 mod 4`
create a self-paired diameter class that demands `b = ¬b`; odd lengths misalign). Two fixed points force a
fixed pair with `b = ¬b` (impossible), so at most one fixed point — the odd-n vertex. This is the tournament
analogue of the classical `4 | ℓ` rule for self-complementary *graphs*, shifted to `ℓ ≡ 2 mod 4` because the
complement here is arc-*reversal*.

## Status
- **Task 1 (resolved, verified n≤6):** `H_n` does **not** fold to a quarter-tiling; the complement-fold is a
  single terminal step. **Corrects the S19 "dyadic tower."**
- **Task 2 (resolved, verified n≤14):** folded class-count `#SC(n) = A051337` (self-converse tournaments); the
  S19 "even rule `#SC(2k)=A000568(2k−1)`" was a **coincidence (n=4,6), breaks at n=8**. The unified rule is
  Burnside over anti-automorphisms with all cycles `≡ 2 mod 4` (+ one fixed point iff n odd); **odd n mirrors
  even (n−1) plus a fixed vertex**; #types `= q(⌊n/2⌋)`.
- **Honest:** the `≡ 2 mod 4` characterization is verified n≤14 and matches the classical self-complementary
  cycle theory (arc-reversal shift); `#SC = A051337` matches the metagraph for n≤7 exactly.

Related: HYP-3810 (the half-tiling isomorphism — this bounds its recursion: one step, and the class-count is
A051337), S18/S19 (R = complement; the two corrected claims), HYP-3808/3809 (the parity decomposition / atlas).
HYP-3811 (this). Scripts: 04-computation/{sc_count_burnside_A051337, mmg_no_quarterfold,
sc_anti_auto_cycletypes}_opus_20260701.py.
