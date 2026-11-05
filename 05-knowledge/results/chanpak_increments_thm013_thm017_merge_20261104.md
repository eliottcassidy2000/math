# Chan-Pak's conjecture: merging THM-013 / THM-017 into the odd-cycle picture, the increment ladder of the Hamiltonian-path count, and what one-vertex extensions cannot do — mac-mini-2026-11-04-chanpak

**Session:** mac-mini-2026-11-04-chanpak (worktree `math-wt-oddprimes-20261104`, second turn of the day; the first turn
is `cpq_defect_family_mirror_folds_20261104.md`).
**Owner directive:** "merge in prior work on 13 and 17 and work to prove chan-pak's conjecture."
**Reading of "13 and 17":** the March-2026 theorems THM-013 (arc-flip formula for `DeltaH`, proved only for `n <= 7`)
and THM-017 (the signed even-odd split, proved for all `n`) of the odd-cycle (OCF) program, which is the natural prior
work on the count `H` whose spectrum Chan-Pak concerns; the numerals `13` and `17` also occur in the lane as the rigid
count `13` (the unique strong realiser; HYP-9453 (3)) and the Fermat count `17` (least order `6`, THM-4861; the first
digit-1 floors `5, 9, 17, 33, ...` of THM-4862), and both are recorded in section 1 for completeness.
**Truth typing:** PROVED (THM-4876 (2), (3) structure, (4), (6), the path-pair formula) + KNOWN (Claim A, the odd-cycle
formula, THM-017) + FINITE-EXACT (extension spectra, increment ladder, end-count ratios, one-centre coverage) +
NUMERICAL (HYP-9469) + OPEN (Chan-Pak).

## 0. Inheritance and concurrency

* **The conjecture.** Chan-Pak, *Computational complexity of counting coincidences* (TCS 2024, arXiv:2308.10214),
  Conjecture 6.2: every odd number other than `7` and `21` is the number of Hamiltonian paths of a tournament. Canon:
  THM-1370 (`7`, `21` never attained; the spectrum is monotone in the order), THM-4094 (completeness is equivalent to the
  strong realisation of the two lanes `p` (odd prime `!= 7`) and `7p`, by order-join multiplicativity THM-1862/1936),
  THM-4744 (strong-unattainable values below `625` are exactly `1, 7, 21`), HYP-9029 (the strong-interval tiling law),
  HYP-9453 (near-join `h + 2^k B`, the density of `B`-sets as the missing lemma), THM-4862/4865 (the one-exit
  extension law and the tribonacci sole-exit tower), THM-4873 (end-count lemma, opus-S77 today), the digital bank
  (session124 and the square-sum overnight packets: every odd `<= 99,999,999` other than `7, 21` attained,
  independently audited to `999,999`). MISTAKE-559 records the one overclaim to avoid (canon proves only the
  omission of `7, 21`).
* **Concurrent today:** opus-S78 reserved THM-4875/HYP-9468 ("old-mark digit lower bound and coverage of odd
  Hamiltonian path counts"); kind-pasteur S10 THM-4874 (uniquely Hamiltonian deletions); Session 189-191 "spectrum"
  commits. My lane-claim message states the exact statements of THM-4876 and defers to S78 for any C2-grammar
  version of the increment bound pushed first.
* **Closest proved mechanism:** Claim A (THM-070 from OCF; now unconditional by Grinberg-Stanley): `H(T) - H(T - v) =
  2 sum_{C ni v} H(T - V(C))`. **Hostile example:** the strong 4-tournament, whose two vertices lying on only one of its two triangles
  (one of out-degree one, one of in-degree one; the other two lie on both triangles and give `5 - 1`) give increment
  `2` from `3` to `5`, the last time an increment of `2` ever happens. **Corrected near miss:** my first
  "extension spectrum" computed `sum_P (insertion slots of v into P)` over the Hamiltonian paths `P` of `X`, which is
  not `H(X + v)` (removing `v` from a path of `X + v` leaves two paths of `X`, not one); the correct quantity is the
  path-pair sum, THM-4876 (1), and the wrong one was discarded before anything was claimed. **Least-used sidecar:**
  THM-012b's insertion decomposition, which is the path-pair formula in the language of orphan paths.

## 1. The merge: what THM-013 and THM-017 are, and what the numerals do in this lane

* **THM-013 (March 2026)** states `H(T) - H(T') = sum_k 2^k Delta(alpha_k)` for an arc flip `T -> T'`, with `alpha_k`
  the number of `k`-sets of pairwise vertex-disjoint odd cycles, and proved it for `n <= 7` (verified to `9`). Its
  "Significance" section hoped to derive OCF (`H = I(Omega, 2)`) from this identity. The order is now reversed: OCF is
  a theorem (Grinberg-Stanley arXiv:2307.05569, Theorem 1.39, with the mod-4 refinement their Theorem 7.1; Irving-Omar
  arXiv:2412.10572 is the matrix-algebra consolidation, Corollary 20), so THM-013 holds for all `n`. Audit C (C1):
  this closure is not new - PROP-001 and THM-002/CONJ-001 (kind-pasteur S12, 2026-03-05) already record
  `Delta H = Delta I(Omega, 2)` for all `n` as a corollary of Grinberg-Stanley; THM-013's status line was merely stale,
  and its file now carries a dated cross-reference section saying so.
* **THM-017 (March 2026)** proves `B(L_i, R_j) = B(L_j, R_i)` for the signed subset convolutions at an arc. The unsigned
  convolutions give the flip increment, `H(T') - H(T) = U(L_j, R_i) - U(L_i, R_j)`: the two theorems are the signed and
  unsigned faces of the same decomposition of Hamiltonian paths at the flipped arc. Audit C (C1): the unsigned face is
  also March-2026 canon (THM-013 "Significance", THM-014's `h_end x h_start` convolutions, THM-017/MISTAKE-008's
  `Delta(S, R) = L_j(S) R_i(R) - L_i(S) R_j(R)`), so THM-4876 (4) is a restatement typed KNOWN, not a new result.
* **The numerals.** `13` is the only count besides `3, 5, 11` with a unique strong realiser through order `9`
  (HYP-9453 (3)); it is not a one-centre count (THM-4876 (6): `N(A) = 6` is missed). `17 = 2^4 + 1` is the one-centre
  count of `A = {1}`, `m = 5`, attained at least order `6` (THM-4861), and the vertex `5` of that tournament ends exactly
  `8 = (17 - 1)/2` of its paths, the extremal case of HYP-9469 (4). In the previous turn `13` and `17` were the first
  primes with `(p - 1)/2` even and `> 2`, where no uniform quantifier-free circulant rule exists (HYP-9463 (4)); that
  is a different object and is not merged into this lane.

## 2. Results (THM-4876)

* **Increment two (PROVED).** `H(T) - H(T - v) = 2` iff `v` lies on exactly one odd cycle, a triangle, with transitive
  complement; then `H(T - v) in {1, 3}`. Proof: Claim A forces one odd cycle with `H(T - V(C)) = 1`; the sub-path lemma
  (a vertex on an odd cycle of length `>= 5` lies on a shorter one) makes it a triangle `(v, a, b)`; the only arc from
  `v`'s out-neighbourhood to its in-neighbourhood is `a -> b`, which pins the structure of `T - v` down to at most one
  triangle. Hence `k_min(h) >= 2` for `h >= 5`, in words: **no tournament with at least five Hamiltonian paths gains
  exactly two by adding a vertex.**
* **Increment four (PROVED structure).** Exactly two odd cycles through `v`, both with transitive complements.
* **The ladder (FINITE-EXACT, orders of `T` up to 9).** `k_min(h) = 2, 3, 4, 5, 5, 5, 6, 7, 9, 9, 8, 10, 9, 11, 12, 11`
  for `h = 5, 9, ..., 39`; the `h`-sets of the increment levels `2k`, `k = 1..6`: `{1,3}`, `{1,5}`, `{3,5,9}`,
  `{1,5,9,11}`, `{1,3,5,9,13,15,17}`, `{3,5,11,13,15,17,19}`.
* **Extension spectra (FINITE-EXACT, orders 3-8 for `X`).** No `X` with `h >= 5` has `h + 2 in E(X)`; the longest run
  of consecutive odd values in `E(X)` has length `<= 24` at order `8`; for `h notin {1, 5}` no `E(X)` of order `<= 8`
  covers `odd[2h + 1, 3h]` (audit C, direct containment; the strong 4-tournament `E = {5, 9, 11, 13, 15}` and its
  paddings do). So "every odd in `[h, 2h - 1]`" cannot be reached from one realiser at any order (PROVED, THM-4876
  (2)), and the window induction is unavailable through order `8` (FINITE-EXACT).
* **End-count bound (PROVED after the checkpoint, THM-4876 (7)).** A non-sink vertex `a` ends at most `(h - 1)/2`
  of the `h` Hamiltonian paths: inserting `a` before its first out-neighbour along a path ending at `a` is an injection
  into the paths not ending at `a` whose last vertex beats `a`, and at least one path ends at an out-neighbour of `a`
  (every vertex of the last strong component ends a path). Equality at the penultimate vertex of `TT_(k+2)` with its
  extreme arc reversed (`h = 2^k + 1`); exhaustive to order `9` (maximum `64/129`). It sharpens THM-4873 (i) for `h >= 5`
  and, through THM-4862 (a)'s law `H(X + v) = 2h - e_a + sigma`, `sigma >= 0`, raises the sole-exit count by at least
  `(h + 1)/2` (`k >= (h + 1)/4` in the half normalisation of HYP-9469 (1)).
* **One-centre family (PROVED).** `H(TT_m + u) = 1 + 2 sum_{a in A, b notin A, a < b} c(b - a)`; `A = {1..i}`,
  `m = i + l` gives `(2^i - 1)(2^l - 1) + 2` (a two-parameter family containing `11, 23, 47, 51, 95, 107, 227, ...`);
  `7` is not a one-centre count; the family realises about half of the odd numbers. A second centre with a monotone cut
  is degenerate (deleting the first centre leaves a transitive tournament).

## 3. Computations

`04-computation/experiments/chanpak_spine_20261104/`: `spine_numeration.py` (session124's spine transfer re-derived
and brute-force checked on 69,174 triples; coverage of the transitive parent and of all parents of order `<= 3`),
`ext_spectrum_true.py` (path-pair formula, checked on 468 pairs; the `[h, 2h - 1]` test), `ext_windows.py` (windows,
increment-2 cases), `increment_ladder.py` (+ `_n8.out`; Claim A checked at every cut of increment `<= 8`, the
ladder), `increment2_order9.py` (+ `.out`; all 191,536 tournaments of order 9 with all cuts), `endcount_ratio.py`,
`two_centre_family.py` (the degenerate second centre, formula checked on 3,584 cases), `extension_spectra.py` (the
discarded wrong quantity, kept for the record). pynauty is not needed here; `gentourng` supplies the tournaments.

## 4. What this does and does not do for Chan-Pak

It does not prove the conjecture. It shows that the spectrum is nowhere "locally connected" at a realiser: the
smallest raise above `h` is about `h/3.5` (conjecturally at least `h/4`), so completeness must be a statement about
unions over many realisers and products, which is exactly the shape of the repository's banks and of THM-4094's lane
reduction. The positive structural content is the exact description of the smallest raises (one triangle with
transitive complement; two cycles with transitive complements) and the closed form of the one-centre family, which
is the "transitive parent" of the digital bank.

## 5. Next obligations

1. (Done after the checkpoint: THM-4876 (7).) Extend the linear increment bound from sole-exit cuts to all cuts
   (HYP-9469 (1)).
2. Prove HYP-9469 (2) at `k = 2` (increment `4` forces `H(T - v) in {1, 5}`) by classifying two odd cycles through `v`
   with transitive complements.
3. Characterise the one-centre value set (HYP-9469 (3)).
4. (Done: `k_min` through `h = 59` from all tournaments of order `9` with all cuts.)

## 6. Audit

Audit C (statements 1-6; `chanpak_increments_audit_C_20261104.md`, independent C census through order 10): PASS WITH
CORRECTIONS, fourteen items, all applied (MISTAKE-715): C1 statement 4 is March-2026 canon (PROP-001, THM-013/014/017,
MISTAKE-008), retyped KNOWN and the THM-013/017 sections redated; C2 the `[2h + 1, 3h]` sentence was false for the
realisers of `h = 5` (and `h = 1`), repaired to `h notin {1, 5}` with a direct containment test; C3 `odd[h, 2h - 1]`
is covered by the seven order-8 tournaments with `h in {1, 3}`; C4/C5 Lass dropped, the mod-4 refinement is
Grinberg-Stanley's Theorem 7.1; C6 only the two single-triangle vertices of the strong 4-tournament give `5 - 3`; C7 the
`N(A) = 3` argument completed; C8 the sole-exit increment is `(h + sigma - e_a)/2`; C9 the equality family relabelled;
C10 the `3k - 1` clause deleted; C11 level sets complete from order 7 of `T`; C12 "every induction" retyped; C13
`k_min >= 3` is FINITE-EXACT; C14 `h/3.5` replaced by the range. Audit C added: statements 2 and 3 exhaustive through
order `10`, `k_min(h) >= h/4` for every `h` at order `10` (least `33/129`, the sole-exit cuts of `R_(k+2)`), the levels
`k = 7..12`, and the end-count bound through order `10` (`128/257`). Audit D (statement 7): pending. Author's
corrections before audit D
(both pushed for about ten minutes in 87ac336954, corrected in the next commit): the first text of statement 7's
consequence asserted `H(X + v) = 2h - e_a` for the sole-exit extension, forgetting the `sigma` term of THM-4862 (a) -
the same two-paths slip as the near miss of section 0 (`sigma` reaches 105 at order 6, and `sigma = 0` in only 9,454 of
195,506 non-sink cases) - and the equality family was labelled `TT_(k+1)` with the arc `(1, k + 1)` reversed at vertex
`k + 1`, which is the order-`(k+1)` member with `h = 2^(k-1) + 1`; the member with `h = 2^k + 1` is `TT_(k+2)` with the
arc `(1, k + 2)` reversed, and the bound is attained at its penultimate vertex (`verify_stmt7_consequences.py`). The
bound and its proof were not affected.
