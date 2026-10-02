# Shaved tournaments, proved versus conjectured: `u(9) = 14`, the excess law fails from `n = 21`, and what the links to lonely-runner harmful relations and Collatz rigidity are (and are not)

**Thread session `thread-2cjcob`, 2026-10-01** (project thread "Write up shaved tournament cores", requested by the
owner; temporary identity, no registered machine ID). The brief: review
[`shaved_tournaments_unavoidable_cores_20261001.md`](shaved_tournaments_unavoidable_cores_20261001.md) (THM-4526)
and its script, tie the findings to the lonely-runner harmful-relations / Collatz-uniqueness line (the sixteenth to
twentieth notes of the S15 series), and write down cleanly what is proved and what is conjectured.

- Script: `04-computation/experiments/shaved_tournaments_n9_20261001.py` (+ `.out`, ALL CHECKS PASSED), with the C
  helper `shaved_tournaments_n9_20261001.c`. It needs nauty and gcc.
- Builds on: THM-4526 (`01-canon/theorems/THM-4526-shaved-tournaments-path-plus-chord-and-the-largest-unavoidable-graphs.md`,
  independently audited, MISTAKE-555), the designations note
  [`tournament_designations_20261001.md`](tournament_designations_20261001.md), HYP-3817 (the excess law), the
  [sixteenth note](collatz_lrc_relations_20260930.md) (Beck–Everett) and the
  [seventeenth note](collatz_functional_uniqueness_20261001.md) (THM-4523, Collatz rigidity).

**Status.**
- PROVED by exhaustive computation (new): `u(9) = 14`, so `κ(9) = 22`; the maximum 9-shavings are 54 isomorphism
  classes. Three different host tournaments give the same 54 classes, a second argument re-proves `u(9) ≤ 14` from
  the list, and every one of the 54 is re-verified by exhaustive embedding counts.
- PROVED (elementary, new): the excess law of HYP-3817 fails for **every** `n ≥ 21` with `n ≡ 0, 2 (mod 3)`
  (Theorem E); the lazy-caterer formula overshoots `κ(n)` for **every** `n ≥ 38` (and at `n = 34`); the symmetry
  penalty `|Aut S| · 2^u ≤ n!` for any shaving with `u` arcs (Lemma S).
- FINITE-EXACT (new): Theorem A at `n = 10` (all 9,733,056 classes); the excess law holds for `n = 3..9`; every
  maximum shaving for `n ≤ 9` is rigid; symmetric tournaments alone block every one-arc extension of the maximum
  shavings at `n = 7, 8, 9` (and are needed at `n = 7, 8`, but not at `n = 9`); the parity, self-converse and
  Hamiltonian-path census; the lonely-runner tight sets in boxes for `k ≤ 8`.
- Independently re-derived here (no new claim): Theorem A for `n ≤ 9`, and `u(7) = 9`, `u(8) = 11` by an orderly
  host search (at `n = 8` this differs from both earlier derivations, which enumerated forward labellings).
- ANALOGY / NUMEROLOGY: every link to the lonely runner and to Collatz in §5. None of them transfers information
  between the problems.
- Two corrections (MISTAKE-557): a sentence of the first note's §5 on HYP-3802, and "exact only for `n ≤ 6`" for
  the lazy-caterer formula (HYP-3798, THM-4526).
- INDEPENDENTLY AUDITED (2026-10-01, blind re-derivation with the auditor's own code and a different search, §9):
  every claim SOUND; four wording corrections in §5 applied.
- Found independently the same day by the shave4 lane, now THM-4533
  ([`THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md`](../../01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md),
  note [`procgen_shave4_20261001_redei_graphs.md`](procgen_shave4_20261001_redei_graphs.md); its own code,
  independently audited): `u(9) = 14` with the same 54 classes (compared class by class here, §2), and the excess
  law at `n = 9`. THM-4533's audit re-checked `u(9) ≥ 14` with its own code; its `u(9) ≤ 14` rests on the lane's
  search. Here `u(9) ≤ 14` rests on three separate searches: inside three hosts, over all 1739 one-arc extensions,
  and the blind auditor's host-free search. THM-4533 also proves the `n ≤ 6` parity pattern (D69) and finds the
  densest odd-everywhere graphs, `r(n)` for `n ≤ 9`.

## 0. Objects

A *shaving* of order `n` (an *`n`-unavoidable digraph* in the literature) is a spanning oriented graph contained in
every `n`-tournament. Shavings are acyclic (they lie in the transitive tournament) and hereditary (spanning
subgraphs of shavings are shavings). `u(n)` is the largest number of arcs of a shaving and `κ(n) = C(n,2) − u(n)`
is the flip rank of HYP-3798. `T(n)` is the number of `n`-tournament classes (A000568), `H_n` the Hamiltonian path
plus the arc from its first to its last vertex (the bypass `D(n,2)`), and `U(n) = C(n,2) − ⌈log₂ T(n)⌉` the counting
bound. A class is *self-converse* (SC) if it is isomorphic to its converse.

## 1. The ledger

| claim | status | evidence |
|---|---|---|
| **A.** `H_n` lies in every `n`-tournament except `C3` and `C3[1,C3,1]`; the number of copies is `H(T) − n·hc(T)`, odd for even `n` | PROVED; CITED (Grünbaum 1971; the exception list via Havet 2000 and El Zein, arXiv:2204.11211, as THM-4533 reports); independent proof in THM-4526 §1, audited twice | proof re-read line by line here (sound); FINITE-EXACT for all classes with `n ≤ 10` (N2); FORMALIZED in Lean for even `n`, for `n = 4, 5` and for the copy count, not for odd `n > 5` ([formalization note](procgen_lean_20261001_formalization.md) §5) |
| **B1.** `u(n) = 0,1,2,4,6,8` for `n = 1..6`; unique maximiser = path + all span-3 arcs | PROVED by exhaustive computation | THM-4526 (B), audited |
| **B2.** `u(7) = 9` (51 classes), `u(8) = 11` (1617 classes) | PROVED by exhaustive computation | now three independent derivations (first script; S15 auditor; the host search here, N3), and the shave4 lane (THM-4533) ran the same host search with its own code |
| **B3.** `u(9) = 14`, `κ(9) = 22`, 54 classes | PROVED by exhaustive computation (**new**) | §2; N4–N6; §9 (blind audit); THM-4533 (shave4 lane) found the same 54 classes |
| **C.** `u(n) = n log₂ n − O(n log log n)` | CITED (Linial–Saks–Sós 1983) | THM-4526 (C) is superseded (MISTAKE-555) |
| **D.** lazy caterer `κ(n) = 1 + C(n−2,2)` | exact for `n = 3..6` and **`n = 9`**; too small at `n = 7, 8`; too large for **all** `n ≥ 38` (PROVED); "exact only for `n ≤ 6`" was false (MISTAKE-557) | §4; N10 |
| **E.** excess law `κ(n) − ⌈log₂ T(n)⌉ = #{SC classes with \|Aut\| > n}` | holds for `n = 3..9` (FINITE-EXACT); **REFUTED** for every `n ≥ 21`, `n ≡ 0, 2 (mod 3)` (PROVED) | §3; N9 |
| **F.** perfect shavings (each class hit exactly once) exist only for `n ≤ 4` | FINITE-EXACT (`n ≤ 60`); OPEN beyond | THM-4526 §4 |
| **G.** path + span-3 arcs has an odd number of embeddings in every tournament for `n ≤ 6` | PROVED in THM-4533 (D3, independently audited), which also shows it fails for every `n ≥ 7` | THM-4526; §2 adds the census for the maximum shavings |
| **H.** every maximum shaving is rigid | PROVED for `n ≤ 6` (unique linear extension); FINITE-EXACT for `n = 7, 8, 9`; OPEN in general (Q2) | N7 |
| **I.** symmetric classes stop the maximum shavings from growing | FINITE-EXACT for `n = 7, 8, 9`: every one-arc extension of a maximum shaving is avoided by a class with `\|Aut\| > 1`. At `n = 7, 8` symmetry is needed (some extensions are avoided only by symmetric classes); at `n = 9` rigid classes alone also block the extensions (not the 14-arc non-shavings, §9), but the most efficient blockers are symmetric. OPEN in general | THM-4526 (B); §2 here, N6 |
| **J.** "`P_7` is the extremal object of the LRC thread (HYP-3802)" (first note, §5) | **not supported**; corrected (MISTAKE-557) | §5.4 |

## 2. The new case `n = 9`

**Theorem B3.** `u(9) = 14` and `κ(9) = 22`. The 14-arc 9-shavings form exactly 54 isomorphism classes.

*Method.* Every 9-shaving embeds in every 9-tournament, in particular in a fixed host `H`. So it suffices to list the
acyclic 14- and 15-arc subgraphs of `H` up to `Aut(H)` and to test each against all 191,536 classes.
- The listing is orderly: a set of arcs is kept if it is the lexicographically least image of itself under `Aut(H)`.
  Deleting the largest arc of such a set leaves such a set, and acyclicity and containment pass to subsets, so a
  depth-first search that extends only kept sets by larger arcs reaches every orbit once.
- Tournaments that avoided an earlier candidate are tried first ("killers", move-to-front). They only prune. A set
  is accepted only after a full pass over the 191,536 classes.
- Three hosts: `C3[C3] = Cay(Z₉, {1,3,4,7})` (`|Aut| = 81`, the unique most symmetric 9-tournament), `TT3[C3]`
  (`|Aut| = 27`) and `C3[1,P7,1]` (`|Aut| = 21`). They give 64, 104 and 66 host orbits of 14-arc shavings, which
  are the **same 54 classes** in all three cases, and no 15-arc shaving in any of them.
- *Second proof of `u(9) ≤ 14`, given the 54.* A 15-arc shaving would contain a 14-arc shaving, so it would be a
  one-arc extension of one of the 54. Their acyclic one-arc extensions form 1739 classes. For each, a brute-force
  search over all `9!` bijections (no pruning at all) finds a class that avoids it.
- *Positive direction, re-verified.* For each of the 54, a plain depth-first count of all embeddings (no look-ahead)
  into every one of the 191,536 classes is at least 1. The minimum over the classes is 1, 2, 3 or 6 (for 10, 16, 26
  and 2 of the 54).
- *Who blocks.* The class that avoids the most of the 1739 extensions is `C3[C3]` itself: it avoids 1018 of them.
  The 15 most frequent blockers all have `|Aut| ≥ 3`, although 188,337 of the 191,536 classes are rigid. They
  include `C3[C3]` (81), regular classes with `|Aut| = 9`, and the four `P_7`-based classes with `|Aut| = 21`.
  - Every extension is avoided by at least 27 classes, always including a symmetric class (`|Aut| > 1`) and a
    rigid one. So the 3,199 symmetric classes alone block all 1739, and so do the rigid classes alone. The
    symmetric ones are far more efficient.
  - At `n = 7, 8` symmetry is needed. Every one-arc extension of a maximum shaving is again avoided by a symmetric
    class, and some only by symmetric classes: 27 of the 706 extension classes at `n = 7`, and 906 of the 33,384 at
    `n = 8`. One 8-vertex extension is avoided by a single class, `C3 ⇒ R5` (`|Aut| = 15`).
  - A greedy blocking set at `n = 9` has 10 classes (`|Aut| = 81, 21, 21, 21, 9, 15, 9, 9, 7, 1`). The rigid class at
    its end is an artifact of the greedy order.

The same host search re-derives `u(7) = 9` (inside `P_7`: 54 orbits, 51 classes) and `u(8) = 11` (inside
`TT1 ⇒ P_7`: 2290 orbits, 1617 classes) in seconds, with no 10- or 12-arc shaving. At `n = 8` this differs from both
earlier derivations, which enumerated all forward labellings.

**Independent confirmation.** The shave4 lane (THM-4533) ran its own orderly search inside `C3[C3]` and `TT3[C3]`
(the only class with `|Aut| = 27`; [its note, §8](procgen_shave4_20261001_redei_graphs.md)). Its raw list of maximum
9-shavings (the 54 `MAX14` lines of `procgen_shave4_20261001_u9search.out`) gives the same 54 classes as here: the
labelg canonical forms agree one for one. Its two densest odd-everywhere 7-graphs are the converse pair listed below
the census table, and its count of 20 of the 54 with an odd number of linear extensions agrees. These comparisons
were run once here; they are not in the script.

**What the 54 look like.**
- All are **rigid** (`|Aut S| = 1`).
- **None is self-converse**: they form 27 converse pairs. (The set of maximum shavings is always closed under the
  converse, because `T ⊇ S` iff `T^conv ⊇ S^conv`.)
- **None contains a Hamiltonian path.** Path + all six span-3 arcs has exactly 14 arcs, but 81 classes avoid it.
  So the pattern that is exact for `n ≤ 6` and leaves an echo at `n = 7, 8` (path + span-3 arcs minus one span-3
  arc) is gone at `n = 9`.
- None has a bipartite underlying graph; all are weakly connected.
- None has an odd number of embeddings in every class.

**The census through `n = 9`** (maximum shavings; N5, N7):

| `n` | `u` | `κ` | classes | rigid | self-converse | with a Hamiltonian path | odd count in every class |
|---|---|---|---|---|---|---|---|
| 4 | 4 | 2 | 1 | 1 | 1 | 1 | 1 |
| 5 | 6 | 4 | 1 | 1 | 1 | 1 | 1 |
| 6 | 8 | 7 | 1 | 1 | 1 | 1 | 1 |
| 7 | 9 | 12 | 51 | 51 | 1 | 4 | **2** (a converse pair) |
| 8 | 11 | 17 | 1617 | 1617 | 15 | 5 | 0 |
| 9 | 14 | 22 | 54 | 54 | 0 | 0 | 0 |

The two odd-everywhere maximum 7-shavings are `01 03 14 16 23 25 34 45 56` and its converse
`01 03 12 16 23 24 35 45 46` (arcs `ij` = `i → j`).

Whether `u(9) = 14` is in the literature was not checked. The standard names to search are "largest digraphs
contained in all `n`-tournaments" (Linial–Saks–Sós) and "`n`-unavoidable digraphs".

## 3. The excess law: true through `n = 9`, false from `n = 21`

klein's HYP-3817 proposed `excess(n) := κ(n) − ⌈log₂ T(n)⌉ = #{SC classes C with |Aut C| > n}`. It predicted
`excess(8) = 4`, which the first note confirmed.

**At `n = 9` it holds again.** `⌈log₂ 191536⌉ = 18`, so `excess(9) = 22 − 18 = 4`. Of the 14 classes with
`|Aut| > 9` (group orders 15, 21, 27, 81), exactly 4 are self-converse. So the law matches for every `n = 3..9`:
`0, 0, 0, 1, 3, 4, 4`. At `n = 10` it predicts `κ(10) = 24 + 1 = 25`, i.e. `u(10) = 20`, one below the counting
bound `U(10) = 21`. That is untested (Q1).

**Lemma X (substitution).** For an `m`-tournament `X` let `X[C3]` replace every vertex by a 3-cycle.
1. In `X[C3]`, every module with at least two vertices that contains a vertex `v` contains the whole block of `v`.
   *Proof.* Let `v → v⁺ → v⁻ → v` be the block `B` and `M ∋ v` a module with `|M| ≥ 2` and `B ⊄ M`. If
   `M ∩ B = {v}`, pick `w ∈ M ∖ B`; `w` beats all of `B` or loses to all of `B`. Since `v⁺ ∉ M` relates to `v` and
   `w` alike and `v → v⁺`, we get `w → v⁺`, so `w` beats `B`; but `v⁻ → v` forces `v⁻ → w`, a contradiction. If
   `M ∩ B` has two vertices, the third vertex of `B` beats one of them and loses to the other. ∎
   So the blocks are intrinsic (the least modules of size `≥ 2`), and `X[C3] ≅ Y[C3]` implies `X ≅ Y`.
2. `(X[C3])^conv = X^conv[C3^conv] ≅ X^conv[C3]`, so `X[C3]` is self-converse when `X` is.
3. `|Aut X[C3]| ≥ 3^m`, since each block can be rotated independently. (In fact `|Aut X[C3]| = 3^m · |Aut X|`:
   automorphisms permute the intrinsic blocks, the induced permutation lies in `Aut X`, every element of `Aut X`
   lifts, and the kernel rotates each block independently.)
4. The same holds for `TT1 ⇒ X[C3] ⇒ TT1` (`3m + 2` vertices). `X[C3]` has no source and no sink, so the added
   source and sink are intrinsic.

N9 checks 1–4 directly for all `X` with `m ≤ 6` (`m ≤ 5` for item 4).

**Theorem E.** The excess law fails for every `n ≥ 21` with `n ≡ 0` or `2 (mod 3)`.

*Proof.* Put `m = ⌊n/3⌋ ≥ 7`. By Lemma X, there are at least `#SC(m)` self-converse `n`-classes with
`|Aut| ≥ 3^m > n`. On the other side, `excess(n) = U(n) − u(n) ≤ U(n) − L(n)`, where `L(n)` is the proven lower
bound from superadditivity and transitive blocks, seeded with `u(1..9)`.
- For `21 ≤ n ≤ 60` the exact counts `#SC(m) = 88, 176, 2752, 8784, …` (A002785, recomputed by orbit counting)
  exceed `U(n) − L(n)` (N9). At `n = 21`: `#SC(7) = 88`, while `U(21) = 65`, so even `u(21) ≥ 0` suffices.
- For `m ≥ 15`: a fixed involution with at most one fixed point is an anti-automorphism of exactly
  `2^((C(m,2) + ⌊m/2⌋)/2)` labelled tournaments, and each class has at most `m!` labellings, so
  `#SC(m) ≥ 2^((C(m,2) + ⌊m/2⌋)/2) / m!`. This exceeds `log₂((3m+2)!) ≥ U(n)` at `m = 15`. From one `m` to the next
  the left side grows by a factor `≥ 2^(m/2)/(m+1) > 11` and the right side by less than 2. ∎

**What this means.** The law is a small-`n` phenomenon. Asymptotically `excess(n) = O(n log log n)` by
Linial–Saks–Sós (`U(n) = log₂ n! + O(1) = n log₂ n − n log₂ e + O(log n)`), while for `n ≡ 0, 2 (mod 3)` the
number of symmetric self-converse classes is at least `2^(cn²)`. The substitution makes the failure explicit and early. HYP-3817's
mechanism ("high symmetry means few labelled copies, so such classes are hard to cover") survives as a tendency;
what fails is the exact count. It is open whether the law holds for `10 ≤ n ≤ 20` and for `n ≡ 1 (mod 3)` (Q1).

## 4. The lazy-caterer formula

The formula `κ(n) = 1 + C(n−2,2)` (HYP-3798; equivalently `u = 2n − 4`):
- is exact for `n = 3..6` and, newly, for `n = 9` (`22 = 1 + 21`);
- is too small at `n = 7, 8` (`12 > 11`, `17 > 16`);
- is too large for **every** `n ≥ 38`, and at `n = 34`.

*Proof of the last point.* The recursive bound gives `u(n) ≥ L(n) > 2n − 4` for `38 ≤ n ≤ 60`. For larger `n`, every
`n`-tournament contains `TT_k` with `k = 1 + ⌊log₂ n⌋ ≥ 6` (Erdős–Moser), so `u(n) ≥ C(k,2) + u(n − k)`. With
`C(k,2) ≥ 2k` and `n − k ≥ 38`, induction gives `u(n) > 2n − 4`. ∎ (THM-4526 stated this for `38 ≤ n ≤ 60` and all
large `n`.)

So the formula is right at `n = 9` by accident of small numbers, not because the pattern resumes. The behaviour in
`10 ≤ n ≤ 37` (except `34`) is open.

## 5. The links to the lonely runner and to Collatz uniqueness, typed

Each link below gives the source, the target, the map, and what is preserved and lost (AGENTS.md). None has a map
between the objects. Each has an exact statement on each side and a shared logical shape.

### 5.1 Parity certificates (harmful relations ↔ odd embedding counts): ANALOGY

- **Lonely-runner side** (sixteenth note §1, CITED Beck–Everett Theorem 2). If `n` has no strict lonely time, the
  window product vanishes identically, so every Fourier coefficient vanishes. Then even- and odd-sum frequency
  vectors balance in every fibre, and that forces a short relation with odd coordinate sum (harmful).
- **Shaving side.** A shaving whose number of embeddings is odd in every tournament (a *Rédei graph* in
  THM-4533) is unavoidable for a parity reason. Examples: the Hamiltonian path (Rédei), `H_n` for even `n`, path +
  span-3 arcs for `n ≤ 6` (D72).
- **Shared shape.** An existence statement whose proof is a mod-2 obstruction to an exact zero or an exact balance.
  Preserved: the role of parity. Lost: all metric content (relation norms; arc counts).
- **New evidence.** Parity certifies some maximum shaving through `n = 7`, and none at `n = 8, 9`. The unique
  maximiser for `n ≤ 6`, and 2 of the 51 at `n = 7`, have odd counts in every class. None of the 1617 at `n = 8` and
  none of the 54 at `n = 9` does: each is hit an even (nonzero) number of times by some class, so the simple parity
  certificate is unavailable for them.
- **The logical roles are mirrored, not shared.** In the lonely runner a harmful relation is necessary for having
  no lonely time, but not sufficient (the sixteenth note's converse failure). On the shaving side an odd count in
  every class is sufficient for unavoidability, but not necessary.
- **The decisive test** (D72), the largest `e` such that some `e`-arc 8-shaving has an odd count in every class, has
  been run by the shave4 lane: `r(8) = 10`, one below `u(8) = 11`, and `r(9) = 11`, three below `u(9) = 14`
  (THM-4533, by the lane's exhaustive search). So parity certificates fall behind unavoidability from `n = 8` on.

### 5.2 Rigid universal objects, symmetric obstructions (Collatz uniqueness ↔ shavings): ANALOGY with exact cores

- **Collatz side** (seventeenth note, THM-4523, PROVED + AUDITED). The Collatz graph has trivial automorphism group,
  and every vertex prime to 3 is identified by its backward tree. A symmetric fibre could only sit at the 3-adic
  point `−1/5`, which is 2-adically odd.
- **Shaving side** (FINITE-EXACT). Every maximum shaving for `n ≤ 9` is rigid (PROVED for `n ≤ 6`). Symmetric
  classes are the efficient obstructions (§2):
  - at `n = 7, 8, 9`, every one-arc extension of a maximum shaving is avoided by some symmetric class;
  - at `n = 7, 8` symmetry is needed: some extensions are avoided only by symmetric classes (27 extension classes
    at `n = 7`, 906 at `n = 8`, where one is avoided only by `C3 ⇒ R5`, `|Aut| = 15`);
  - at `n = 9` it is not needed for the extensions, since every extension is also avoided by a rigid class. Still,
    `C3[C3]` (`|Aut| = 81`, the unique most symmetric class) alone blocks 1018 of the 1739, and the 15 most efficient
    blockers are all symmetric. One level down symmetry is needed again: the 14-arc graphs contained in every rigid
    9-tournament form 250 classes, and only the 54 shavings survive the symmetric ones (audit, §9).

  The same-day designations note finds two families of singled-out tournaments for `n ≤ 8`: compositions built
  from 3-cycles (`C3`, `TT2[1,C3]`, `TT2[C3,1]`, `C3[1,C3,1]`, `TT2[C3,C3]`) and near-Paley ones (`R5`, `G_par`,
  `P_7`). The blockers here come from the same two families or combine them: `P_7`; `C3[C3]`, a composition of
  3-cycles at `n = 9` (beyond that note's range); and `C3 ⇒ R5 = TT2[C3,R5]`.
- **The exact core on the shaving side** is orbit counting:
  - **Lemma S.** If `S` is a shaving with `u` arcs, then `|Aut S| · 2^u ≤ n!`. *Proof.* `S` has `n!/|Aut S|`
    labelled copies, each lies in `2^(C(n,2)−u)` labelled tournaments, and together they must cover all
    `2^C(n,2)`. ∎
  - Completions that differ by an automorphism of `S` are isomorphic, so the completions meet at most as many
    classes as there are `Aut(S)`-orbits on them.
  - A class `T` receives `emb(S,T)/|Aut T|` completions, so symmetric classes are the scarce targets.
- **This does not prove rigidity.** At `n = 9`, Lemma S allows `|Aut S| ≤ 22`, and a shaving with `|Aut S| = 2`
  would still have at least `2^21 > T(9)` completion orbits. Rigidity of the maximum shavings is an observation (Q2).
- **Shared shape.** The universal object (the Collatz graph; a maximum shaving) has no symmetry. The would-be
  exceptional objects (a second Collatz component; an obstructing tournament; a lonely-runner tight instance) are
  where symmetry or balance would have to live. Preserved: the asymmetry statement. Lost: everything else. No
  information transfers between the problems.

### 5.3 Circulant tournaments and tight lonely-runner sets: NUMEROLOGY (killed by a probe)

An exact enumeration (N11) of the primitive tight sets `{n_1 < … < n_k}` (`gcd = 1`, `lon = 1/(k+1)` exactly)
gives:
- besides `{1, …, k}`, only `{1,3,4,7}` (`k = 4`, speeds `≤ 40`), `{1,3,4,5,9}` (`k = 5`, `≤ 30`), and
  `{1,2,3,4,5,7,12}` and `{1,4,5,6,7,11,13}` (`k = 7`, `≤ 26`);
- none for `k = 6` (speeds `≤ 26`) and `k = 8` (`≤ 22`). Larger boxes were searched once while preparing this note
  (`k = 4` to 120, `k = 5` to 60, `k = 6` to 40, `k = 7` to 32, `k = 8` to 26) with the same result.

For `k ≤ 5` this agrees with the sixteenth note's census.

The tempting pattern:
- `{1, …, k}` is the connection set of the rotational tournament on `Z_{2k+1}`;
- `{1,3,4,7}` is the connection set of `Cay(Z₉, {1,3,4,7}) ≅ C3[C3]`, the host of §2;
- `{1,3,4,5,9} = QR₁₁`, the Paley tournament `P_11`.

The probe kills it at `k = 7`: `3 + 12 = 15` and `4 + 11 = 15`, so neither set is a tournament connection set mod 15.
Recorded so that it is not promoted (the MISTAKE-228 / MISTAKE-554 lineage).

As a DICTIONARY entry for the Tournament Clock thread, with no implication: `C3[C3]`'s connection set contains
`{1,4,7}`, the unit squares mod 9 (the owner's squares-mod-9 cycle `0,1,4,0,7,7,0,4,1`).

### 5.4 Finite checks per parameter, no uniform law (the sixteenth note's typing): ANALOGY

`κ(n)` is a finite check at each `n`, like LRC(`k`) and like the cycle half of Collatz at each period. The uniform
statements available are asymptotic (Linial–Saks–Sós) or cited exception lists (Grünbaum; Benhocine–Wojda;
Havet–Thomassé). The one proposed exact uniform law, the excess law, holds through `n = 9` and fails from `n = 21`.
That is the same lesson as for the other two problems: per-parameter checks work for small parameters, and no
uniform mechanism is in hand.

**Correction to the first note (MISTAKE-557).** Its §5 placed under "STRUCTURAL (proved here)" the sentence "the
same tournament [`P_7`] is the extremal object of the flip-rank thread (HYP-3805) and of the LRC thread
(HYP-3802)".
- The flip-rank half is right: it is the same problem, `κ(7)`.
- HYP-3802 does not make `P_7` an extremal object of the lonely runner. It attaches the Paley orientation to the
  seven roots of `−1` in the lonely measure of the tight core `{1, …, 13}`, as a chosen dictionary.
- The lonely-runner extremal objects are the tight sets of §5.3, and `{1,2,4} = QR₇` is not one of them: its
  loneliness is at least `1/3`, from `t = 1/3`.

## 6. Open questions

Status of the first note's directions:
- **D69** (why path + span-3 has odd counts for `n ≤ 6`): settled in THM-4533 (D3, a proof by deletion–reversal
  steps; independently audited). New data here: among maximum shavings the odd-everywhere ones stop after `n = 7`
  (§2).
- **D70** (`u(9)`): **resolved**, `u(9) = 14` (§2), found independently here and in THM-4533.
- **D71** (the constant): closed by Linial–Saks–Sós; the second-order term is open.
- **D72** (odd-everywhere graphs): partly answered in THM-4533 (infinite families, `r(n)` for `n ≤ 9`). Open
  (HYP-9173): is `r(n)` linear or of order `n log n`?

New:
- **Q1.** `u(10)`. The excess law predicts 20; the counting bound is 21. More generally, does the law hold for
  `10 ≤ n ≤ 20` and for `n ≡ 1 (mod 3)`, and what replaces it?
- **Q2.** Is every maximum shaving rigid? A symmetric maximum shaving at some `n` would settle it.
- **Q3.** For which `n` is the lazy-caterer formula exact? Known: `3..6` and `9`; never for `n ≥ 38` or `n = 34`.
- **Q4.** Is a maximum `n`-shaving ever self-converse again after `n = 8`? The counts are 1, 15, 0 for `n = 7, 8, 9`.

## 7. What this thread checked in the first note

The first note's claims were re-derived before the S15 audit landed on `main`. The two audits agree.
- **Theorem A.** The odd-`n` proof (Steps 1–5, Lemmas 1–2) was re-read line by line and is sound. The theorem was
  confirmed exhaustively for `n ≤ 10`.
- **Theorem B.** The original script was re-run with identical output. The original `u8.c` search was re-run
  against an independent nauty class list. The host search here was added (N3).
- **Literature.** Grünbaum (via Havet 2000), Benhocine–Wojda and Linial–Saks–Sós were found independently. They
  are the corrections that MISTAKE-555 records.
- **New corrections (MISTAKE-557):** the §5 sentence above, and "exact only for `n ≤ 6`" for the lazy-caterer
  formula in HYP-3798's resolution note and THM-4526's `related` line.

## 8. Reproduction

```bash
python 04-computation/experiments/shaved_tournaments_n9_20261001.py          # 57 checks, about 8 minutes on 4 cores
python 04-computation/experiments/shaved_tournaments_n9_20261001.py --full   # 62 checks, about 36 minutes: adds n = 10 for Theorem A and two more hosts
```

The script needs nauty (`gentourng`, `labelg`, `converseg`, `countg`, `pickg`; Debian names `nauty-*` are found
too) and `gcc` with OpenMP. Without them it prints the recorded results. The `.out` is from a `--full` run. The
parallel searches visit candidates in a varying order, so the "full checks" counts and the candidate numbers in
the output can change from run to run; the PASS lines do not.

## 9. Independent audit (2026-10-01)

A blind auditor re-derived every new claim. It was a separate agent, working from a claims file with the
definitions, statements and proofs, with its own code, and without reading this note or its script.
- **A different search for `u(9) = 14`.** The auditor searched level by level over isomorphism classes of all
  acyclic graphs on 9 vertices, with no host tournament. A candidate with `e` arcs is kept only if all its one-arc
  deletions survived level `e − 1`. Levels up to 13 are filtered by small windows of tournaments, which keeps a
  sound superset. Levels 14 and 15 are tested against all 191,536 classes.
  - Survivors: 342,658, 642,806, 447,862 and 56,646 at levels 10 to 13; exactly 54 at level 14; none of the 1739
    candidates at level 15.
  - A second run along a different filter path gave a byte-identical list of 54. Every removal in it (about 13.7
    million deletion witnesses and 328,000 kill certificates) was verified by separate code.
- **Controls.** The unfiltered class counts reproduce A003087 for `n = 5, 6, 7`. The same code reproduces 51
  classes at `n = 7` and 1617 at `n = 8`. Embedding counts agree with brute force over all `n!` bijections. Hostile
  runs behave as they should: with rigid tournaments only, 250 classes survive at level 14; with `C3[C3]` removed,
  still 54. The first of these is also a finding: 196 classes of 14-arc graphs lie in every rigid 9-tournament but
  are avoided by some symmetric one. It comes from the auditor's run only and is not re-checked by this note's
  script.
- **Verdicts.** SOUND: `u(9) = 14`, the 54 and all their stated properties, the minimum counts; the `n = 7, 8`
  censuses; the excess law at `n = 9`; Lemma X, Theorem E and every number behind them (including `#SC(10) = 8784`
  by nauty over all 9,733,056 classes); the lazy-caterer statements; Lemma S; Theorem A at `n = 10` (by an
  independent search); the tight sets, by exact integer arithmetic, including the larger boxes; both MISTAKE-557
  corrections; the blocker counts. SOUND WITH CORRECTION: the typed links of §5.
- **Corrections applied.**
  - §5.1: an odd count is sufficient but not necessary for unavoidability, the mirror image of a harmful relation
    (necessary but not sufficient). "Parity certifies the maximum shavings" became "some maximum shaving", and
    "needs the full search" became "the parity certificate is unavailable".
  - §5.2: the designations note's family is now listed exactly; it does not contain `C3[C3]`. Symmetric blockers
    are not needed at `n = 9`: the auditor found every extension also avoided by a rigid class. N6 now checks both
    directions at `n = 7, 8, 9`.
  - MISTAKE-557: "exact only for `n ≤ 6`" is false, not merely unsupported.
- The auditor's code is not in the repository.
