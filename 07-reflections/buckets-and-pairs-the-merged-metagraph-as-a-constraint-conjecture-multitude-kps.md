# Buckets and pairs: the merged metagraph as a distribution problem — #SC is even, reconstruction is rigid only at n≤4, the blue subgraph IS the half-tiling, and a conjecture multitude

*kind-pasteur-2026-07-01-S13. Taking the owner's frame literally: `2^m` labeled tilings are distributed into buckets (merged nodes) of three kinds — pure-black (NS, even count), mixed (SC, odd), pure-blue (SC, odd) — and the blue/black lines are a perfect matching `t ↔ φ(t)` (complement tiling) that realizes the bucket sizes as degrees. Grounding (n≤6) + mining the repo (half-tiling THM-549/550, THM-281/282/584/582, LEM-003, klein-S13) yields three new results, a settled open question, and a large annotated conjecture list.*

## The process, stated precisely (natural numbers in buckets)

- **Buckets** = merged nodes, of exactly three types: `A` = pure-black (NS), `B` = mixed (SC), `C` = pure-blue (SC).
- **Contents** = tilings; bucket size = tiling count = `H/|Aut|` (SC) or `2H/|Aut|` (NS).
- **The pairs** = the `2^{m-1}` lines, a perfect matching `t ↔ φ(t)` on the `2^m` tilings (`φ` = flip all `m=C(n-1,2)` tiles). A line is **blue** iff its tiling is `σ`-fixed (grid-symmetric), else **black**.
- **Realization**: bucket size = degree in the line multigraph (a self-loop `{t,φ(t)}` inside one bucket = `+2`).

## Three results (two new, one settling an open question)

**R1 — `#SC` (self-converse classes) is EVEN for all `n≥3`.** The mining flagged this as "unknown if a theorem." It is one line: `2^m = Σ_{buckets}(size) = Σ_{NS}(even) + Σ_{SC}(odd) ≡ #SC (mod 2)`, and `2^m` is even for `m≥1` (i.e. `n≥3`). So `#SC ≡ 0 (mod 2)`. (Equivalently `V_merged=(A000568+SC)/2 ∈ ℤ` with `A000568` even.) Verified: `#SC = 2,2,8,12,88` for `n=3..7`.

**R2 — Reconstruction is rigid only at `n≤4`.** The owner's hypothesis "the metagraph *is* its constraints" is **true at `n≤4`, false from `n=5`**: the local data `(category, blue-deg, black-deg per node)` uniquely determines the multigraph at `n=4` (0 legal degree-preserving colour-respecting 2-swaps), but at `n=5` there are **92 black + 19 blue** legal swaps, and `88380` black swaps at `n=6`. So for `n≥5` the bucket/parity/colour constraints are **necessary but not sufficient** — the metagraph carries extra information: the specific `φ`-pairing of tournament classes. The `H`-gradient is *part* of it (the transitive `H=1` tiling always flips to the unique SC class with `H=2^{n-2}+1`, a blue line — the "big SC neighbour"), but `φ` is not cleanly `H`-antitone (H-sums span a range), so `H` alone does not close reconstruction. **The reconstruction defect (number of realizations of the bucket data) is a new complexity measure of the metagraph, zero exactly at `n≤4`.**

**R3 — The blue subgraph IS the half-tiling model.** The grid-symmetric (`σ`-fixed) tilings — on which *all* blue lines live — number `2^{⌊(n-1)²/4⌋}` (exponent `k²` odd `n=2k+1`, `k(k-1)` even `n=2k`; klein-S13's square/pronic), and `⌊(n-1)²/4⌋ = A002620` is exactly the **half-tiling cell count** (THM-549/550, the `σ`-quotient fundamental domain, `σ=` complement `φ(T^op)`). So last turn's "blue-fold recursion" is literally the **half-tiling model**: the blue lines are the pairing process one fold down, and THM-550's parity recurrences (`even: h(n)=2h(n-1)-h(n-2)`; `odd: 4-term`) drive the blue counts. **Blue = the σ-fixed half = the Rédei-witness half; black = the σ-moved half = the LRC-existence half** (mac-mini-S6's two-index split). The whole tripartition is one shadow of the complement involution `σ`.

## Which tilings share a bucket (the fiber, from LEM-003)

Two tilings share a bucket iff their tournaments are isomorphic (or transpose-paired). Concretely: the tilings of a class are the **Hamiltonian paths of a representative, taken as base-path rerootings, modulo `Aut`** — and `Aut` acts **freely** on them (LEM-003, universal), so the count is exactly `H/|Aut|` (integer). The `σ`-fixed (blue) tilings in an SC class number **odd**-many (THM-281), and are the "palindromic" reps — tying blue to THM-582 (`H ≡ #palindromic paths mod 2`). So the symmetry restriction the owner sensed is: **a bucket's contents = free `Aut`-orbits of Hamiltonian paths, with an odd `σ`-fixed core iff the class is SC.**

## The conjecture multitude

**Counts / closed forms**
1. `f =` #anti-diagonal tiles `= ⌊(n-1)/2⌋` (verified n=4..7). *(easy)*
2. grid-sym tilings `= 2^{⌊(n-1)²/4⌋}`; blue lines `= 2^{⌊(n-1)²/4⌋-1}`; black `= 2^{m-1} −` that. *(from R3)*
3. `#A` (pure-black NS) `= (A000568−SC)/2 = 0,1,2,22,184` (n=3..7) `→ A000568/2`. *(identity)*
4. **`#pure-blue` (buckets) `= 1,3,2,?` — find a closed form.** Distinct from the blue *tiling* count (R3); `#C` is small and non-monotone. Target: compute n=7,8; is it `≈ SC(n−?)`?
5. **`#mixed = SC − #pure-blue`**; mixed carries both colours (the bridge) *and* the whole `H`-gradient.
6. Sum of pure-blue sizes `= 1,5,2`; sum of SC sizes `≡ #SC ≡ 0 (mod 2)`; sum of NS sizes even.

**Pairing / H-structure**
7. **Transitive (`H=1`) flips (blue) to the unique SC class with `H=2^{n-2}+1`.** *(provable; the principal-line base)*
8. The `H=2^{n-2}+1` class is the **blue hub** (max blue-degree). *(observed n≤6)*
9. `φ` (complement-tiling) is `H`-*raising on the low end but not globally antitone*; characterize the `H`-pair spectrum of blue vs black lines. *(open)*
10. Blue-line `H`-pairs lie on the SC-spine `H`-values `{1,2^{n-2}+1,…,H_max}`; the spine is a sub-poset under `φ`. *(open)*

**Self-loops / sea (both onset at n=6)**
11. **Blue self-loops occur only on mixed buckets**; count `1,0,2` (n=4,5,6) — parity/formula? *(open; note the `0` at n=5)*
12. **Black self-loops first appear on pure-black buckets at `n=6`** (24 of them). Onset criterion: a class whose complement tiling is iso to itself. *(open)*
13. **Sea-onset** (`pure-black↔pure-black` black lines) at `n=6` (290). Criterion: two distinct NS pairs whose complement tilings meet — first realizable at `n=6`. (klein-S13 has a partial characterization — merge.) *(open)*
14. Sea fraction `#(Bk–Bk)/#black → 1`; rate? *(open)*

**Reconstruction / rigidity**
15. **Reconstruction defect `= 0` iff `n≤4`**; for `n≥5` it grows. Find the *minimal* extra invariant (beyond category+degree) that closes it — candidate: the `φ`-image of each bucket (the complement-tiling node-map), partially graded by `H`. *(open, R2)*
16. The metagraph `=` bucket-degrees `+` the `σ`-quotient (half-tiling) `+` the `φ`-matching; the three data are independent for `n≥5`. *(open)*

**Half-tiling recursion**
17. **Blue subgraph `≅` the pairing process on the half-tiling** (the `σ`-quotient), giving an exact halving recursion via THM-550. *(open, the cleanest structural target)*
18. The half-tiling's fixed diagonal (SC spine) `↔` pure-blue`∪`mixed; off-diagonal free pairs `↔` the sea. *(THM-549 alignment)*

## Concrete next targets (ranked)

- **(a) Prove R3 as an isomorphism** — blue subgraph = merged metagraph of the half-tiling — then blue counts satisfy THM-550's recurrences exactly (closes conjectures 2,17,18 and gives `#pure-blue` a recursive handle).
- **(b) The reconstruction-closing invariant** (R2/15): is `(category, degree, H)` sufficient? Compute the swap count restricted to `H`-preserving swaps at n=5,6 — if it drops to 0, `H` is the missing constraint; if not, identify the next.
- **(c) `#pure-blue` closed form** (conj 4): n=7,8 census + OEIS.
- **(d) Sea/self-loop onset criterion** (conj 12,13): merge with klein-S13; characterize the first `n=6` pure-black self-loop and `Bk–Bk` line arithmetically.

## Honest status

- **Proved:** R1 (`#SC` even, `n≥3`); R2's rigidity computation (n≤6); the fiber = `Aut`-orbits of Ham paths (LEM-003); tiling count = degree, black-deg even, tripartite, SC-odd/NS-even (prior S12).
- **Established by synthesis:** R3 (blue = half-tiling) via `⌊(n-1)²/4⌋=A002620` matching THM-549/550 — the *identification* is clean; the *isomorphism of processes* (target a) is conjectural.
- **Open:** the whole conjecture multitude above; the owner's "metagraph ≡ constraint" holds only at `n≤4` and needs the `φ`-matching as extra data for `n≥5`.

— Related: `the-blueblack-line-pairing-is-a-degree-tiling-count-realization-kps.md` (S12, the tripartite/parity theorem), THM-549/550 (half-tiling = complement quotient), THM-281 (SC sizes odd), THM-282 (blue = SC spine), THM-584 (complement = antipodal, level-parity spectrum), THM-582 (H ≡ #palindromic mod 2), LEM-003 (Aut free on Ham paths), klein-S13 (blue-is-SC-spine, squares/pronic, sea-onset), HYP-3540 (level-multiplicity), `merged-metagraph-invariants.md`, `geometric-alignment-of-merged-metagraph.md`, CLAUDE.md. Scripts: `merged_metagraph_line_pairing_kps.py`, `merged_metagraph_buckets_conjectures_kps.py`, `merged_metagraph_extra_constraint_kps.py` (+ .out). Not a HYP reservation (R1 is a corollary; the rest are conjectures/synthesis).
