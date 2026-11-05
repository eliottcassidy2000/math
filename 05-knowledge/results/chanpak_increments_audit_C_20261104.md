# Audit C (independent, adversarial) of THM-4876, HYP-9469, the chanpak note, the THM-013 / THM-017 dated sections and the `chanpak_spine_20261104` scripts — 2026-11-04

**Auditor:** audit C (independent agent, worktree `math-wt-oddprimes-20261104`, commit `1abbd516be` plus the completed
`increment2_order9.out`).
**Audited:** (1) `01-canon/theorems/THM-4876-deletion-and-flip-increments-...-for-chan-pak.md`; (2)
`05-knowledge/hypotheses/HYP-9469-minimal-deletion-increment-grows-linearly-...-for-chan-pak.md`; (3)
`05-knowledge/results/chanpak_increments_thm013_thm017_merge_20261104.md`; (4) the "Upgrade (2026-11-04 ...)" section of
`THM-013-arc-flip-delta-I.md` and the "Cross-reference (2026-11-04 ...)" section of `THM-017-even-odd-split-proof.md`;
(5) `04-computation/experiments/chanpak_spine_20261104/*` and its two outputs.
**Method:** everything was recomputed with code written for this audit (no code shared with the session):
`04-computation/experiments/chanpak_spine_20261104/audit_C/auditC_census.c` (subset-DP brute force in C; deletion census,
extension spectra by direct DP on `X + v`, random search), `auditC_lib.py`, `auditC_formulas.py`, `auditC_structure.py`,
`auditC_ext_post.py`; all outputs are in the same directory (`del3.out` … `del9.out`, `del10_part*.out` + `del10_merged.out`,
`ext3.out` … `ext8.out`, `rand.out`, `auditC_*.out`, `session_rerun_*.out`). The session's four census scripts were also
re-run at orders ≤ 8 (`session_rerun_*.out`); `increment_ladder.py 8` reproduces the committed `increment_ladder_n8.out`
byte for byte.

## Verdicts

| File | Verdict |
|---|---|
| THM-4876 | **PASS WITH CORRECTIONS.** Statements 1, 2, 3, 6 and all census numbers verified (orders of `T` up to **10**, exhaustive). Corrections: statement 4 is a restatement of March-2026 canon (PROP-001, THM-013/014/017, MISTAKE-008), not new; two sentences of statement 5 are false as written (`h ∈ {1, 5}` resp. `{1, 3}` exceptions); the "Lass" attribution is unsupported; the `N(A) = 3` argument is incomplete; the strong-4-tournament realiser is misdescribed; "impossible" for the window induction is only FINITE-EXACT. |
| HYP-9469 | **PASS WITH CORRECTIONS.** All quoted numbers verified; the three pre-registered tests PASS (independently, at orders 9 and 10). Corrections: the title's "`h <= 3k - 1 (?)`" is contradicted by the hypothesis's own data; the sole-exit increment is `(h + σ − e_a)/2`, not `(h − e_a)/2`; the equality case of (4) names the wrong vertex and the wrong order; "(stable from order 6 on)" refers to the order of `X`; (5) overclaims "every induction". |
| Note | **PASS WITH CORRECTIONS.** Faithful on inheritance, MISTAKE-559, the digital bank, THM-4094/4744/4861/4862/HYP-9453; the "corrected near miss" is real and withdrawn. Corrections: section 1 claims novelty for THM-4876 (4) that belongs to PROP-001/THM-013/014/017 (March 2026); the mod-4 refinement is Grinberg–Stanley's Theorem 7.1, not Irving–Omar's; "Lass"; the hostile-example wording; the `[2h + 1, 3h]` sentence. |
| THM-013 upgrade section | **PASS WITH CORRECTIONS.** Mathematically correct; the closure is not new (PROP-001 and THM-002, kind-pasteur S12, 2026-03-05, already record `ΔH = ΔI` for all `n` as a corollary of Grinberg–Stanley); the status line "PROVED for all n on 2026-11-04" and "also Lass" must be repaired. |
| THM-017 cross-reference section | **PASS WITH CORRECTIONS.** Correct; the same dating repair ("closed by citation" = THM-002/CONJ-001/PROP-001, 2026-03-05). |
| Scripts and outputs | **PASS.** Re-runs reproduce the committed outputs; `increment2_order9.out` is now complete and its table agrees with my exhaustive order-10 deletion census. Two housekeeping items: `endcount_ratio.py` has no committed `.out` although HYP-9469 (2), (4) quote it; the `[2h+1, 3h]` test in `ext_windows.py` depends on the tie-breaking of "longest window" and should be a direct containment test. |

Nothing audited here claims progress on Chan–Pak beyond what is proved; the conjecture stays OPEN (note §4, THM-4876
status). The repository truth types are respected except where noted (statement 4 should be KNOWN-by-repository).

## 1. What was recomputed (all PASS unless stated)

**(a) Statement 1.** Path-pair formula implemented from the bijection (`auditC_formulas.py: path_pair_H`) and compared
with a direct DP on `X + v`: 4,052 `(X, cut)` pairs over all isomorphism classes of orders 2–6 and 1,098 labeled pairs of
orders 1–4, **0 mismatches**. The derivation (delete `v` from a Hamiltonian path of `X + v`; `P_1` ends in `B` or is
empty, `P_2` starts in `A` or is empty; sum over `S = V(P_1)`) is correct. Claim A in the form
`H(T) − H(T − v) = 2 Σ_{C ∋ v} H(T − V(C))` checked at **every** `(T, v)` of orders 3–8 (55,040 pairs at order 8), 0
mismatches; the odd-cycle formula `H = I(Ω, 2)` checked by my own cycle/independent-family enumeration at all tournaments
of orders 3–7, 0 mismatches. Proof 1 is correct: two odd cycles through `v` share `v`, so an independent family contains
at most one of them, and the families containing none are exactly those of `T − v` (this is THM-070 verbatim).

**(b) Statement 2.** Every step of the proof was checked by hand (see §3 for the only textual flaws) and the
conclusion exhaustively: over **all tournaments of orders 3–10** (9,733,056 at order 10) and all vertices `v`, every
pair with `H(T) − H(T − v) = 2` has `H(T − v) ∈ {1, 3}`, `H(T) ∈ {3, 5}`, exactly one odd cycle through `v`, that cycle a
triangle `(v, a, b)` with `T − {v, a, b}` transitive, `a → b` the only arc from `N⁺(v)` to `N⁻(v)`, never both an
in-neighbour of `a` in `A'` and an out-neighbour of `b` in `B'`, and at most one odd cycle in `T − v`
(`auditC_structure.out`: 3, 8, 13, 18, 23, 28, 33, 38 pairs at orders 3–10, profiles `(h, H) = (1, 3)` and `(3, 5)` only,
**0 violations**). The lemma (a vertex on an odd cycle of length ≥ 5 lies on a shorter odd cycle, hence on a triangle)
holds at every `(T, v)` with an odd cycle through `v`, orders 5–8 (53,976 pairs at order 8). Random search at orders 10
(1.8 M samples: `TT_10` with 1–6 random arcs reversed, and uniform), 11 (300 k) and 12 (90 k): **0 violations**
(`rand.out`). Converse symmetry is available: `H(T^op) = H(T)`, `(T − v)^op = T^op − v`, and reversing all arcs swaps
`A ↔ B`, `a ↔ b`, so the hypothesis and the conclusion are invariant. The converse realisations: `C_3` at any vertex
(`3 − 1`) PASS; the strong 4-tournament: see correction C6 (only two of its four vertices give `5 − 3`).

**(c) Statement 3.** Correct by Claim A: a sum of positive odd integers equal to 2 has exactly two terms, both 1, and
`H(T − V(C)) = 1` iff `T − V(C)` is transitive (Rédei). The census is confirmed and extended: at orders 3–10 every
increment-4 pair has exactly two odd cycles through `v` with transitive complements and `H(T − v) ∈ {1, 5}`; the only
profiles are (`h = 1`, two triangles, `H = 5`) and (`h = 5`, a triangle and a 5-cycle, `H = 9`). Increment-6 set
`{3, 5, 9}`, increment-8 set `{1, 5, 9, 11}` through order 10 (`del10_merged.out`).

**(d) Statement 4.** The flip formula was re-derived and implemented (`U_conv`): over all arcs of all tournaments of
orders 3–6 (990 arcs) and all labeled tournaments of orders 3–5 (10,648 arcs), `H(T') − H(T) = U(L_j, R_i) − U(L_i, R_j)`
with **0 mismatches**; `U(L_i, R_j)` equals the brute-force number of Hamiltonian paths of `T` with `i` immediately before
`j` (0 mismatches); `U(L_j, R_i)` is the same whether computed in `T` or `T'` (0 mismatches); THM-017's
`B(L_i, R_j) = B(L_j, R_i)` holds in every case. The index conventions coincide: THM-017 defines
`B(L_i, R_j) = Σ_{S ⊆ W} (−1)^{|S|} L_i(S ∪ {i}) R_j({j} ∪ (W \ S))`, exactly the alternating version of `U`. The deduction
`ΔH = ΔI(Ω, 2)` for all `n` from `H = I(Ω, 2)` is trivially correct. **Priority (see C1):** both halves are already in
canon since March 2026.

**(e) Statement 5.** Extension spectra recomputed by direct DP on `X + v` for every `X` of orders 3–8 and every cut
(no path-pair formula used): the `k_min` table is reproduced exactly (`auditC_ext_post.out`; identical at orders 7 and
8; cross-checked against my deletion census at every order); the `k_min` values of the session's completed
`increment2_order9.out` for `41 ≤ h ≤ 59` (`14, 15, 14, 13, 13, 15, 14, 18, 18, 17`) agree with my **exhaustive order-10
deletion census** (`del10_merged.out`), which also shows `k_min(h) ≥ h/4` for every `h` attained at order 10 (`h ≤ 3357`;
minimum ratio `33/129 = 0.2558`). No `h + 2` for `h ≥ 5` (proved anyway). Longest window of consecutive odd values at
order 8: length 24, `[2763, 2809]` at `h = 657` (matches). Increment-2 cases by `h` at order 8: `{1: 1, 3: 6}` (matches).
**Two claims fail as written** (C2, C3 below): `[2h + 1, 3h] ⊆ E(X)` holds for the realisers of `h = 5` (the strong
4-tournament padded: `E = {5, 9, 11, 13, 15}`; 5 classes at order 8) and trivially for `h = 1`; `odd[h, 2h − 1] ⊆ E(X)`
holds for the 7 order-8 tournaments with `h ∈ {1, 3}`. The session's own `ext_windows.py` output reports nonzero
counts of its own criterion at orders 3–7 (1, 2, 3, 4, 1 of 2, 4, 12, 56, 456) and `0 of 6880` only at order 8.

**(f) Statement 6.** `H(TT_m + u) = 1 + 2 N(A)` by brute force for `m ≤ 7` and all 254 `(m, A)`: 0 mismatches. The
derivation via the odd-cycle formula (every odd cycle passes through `u`, so `I(Ω, 2) = 1 + 2·#cycles`; a cycle
`u → a → S → b → u` has `|S| + 3` arcs; `2^{b−a−2}` even subsets for `b − a ≥ 2`) is correct. Closed form
`(2^i − 1)(2^l − 1) + 2` for `A = {1..i}`, `m = i + l`: numeric for `i, l ≤ 8` and symbolic (SymPy) PASS. `A = {1}`
gives `2^{m−1} + 1` (the reversed-extreme-arc tournament of order `m + 1`). `N(A) = 3` is indeed impossible (attained
`N` below 12 over all words of length ≤ 16: `0, 1, 2, 4, 5, 7, 8, 9, 11`), but the printed argument is incomplete (C7).
Coverage: **50,006 of the 100,000 odd numbers below 200,000**, first misses
`7, 13, 21, 29, 35, 41, 43, 49, 55, 57, 59, 67, 71, 75, 77, 83, 89, 93, 101, 105, 109, 111, 113, 119, 125, …` (matches
both files). The second-centre degeneracy claim is correct (`X − u = TT_m + v` is `TT_{m+1}` when the cut of `v` is
monotone, so `X` is one-centre on `m + 1` vertices).

**(g) HYP-9469.** Numbers verified: `k_min` table; the level sets `k = 1..6` (`{1,3}, {1,5}, {3,5,9}, {1,5,9,11},
{1,3,5,9,13,15,17}, {3,5,11,13,15,17,19}`), unchanged at order 10; end-count maxima over non-sinks `1/3, 2/5, 4/9, 8/17,
16/33, 32/65, 64/129, 128/257` at orders 3–10 (all equal to `(h − 1)/(2h)` at `h = 2^{n−2} + 1`); `e_a ≤ (h − 1)/2`
for every non-sink vertex of every tournament of order ≤ 10 (97,139,015 checks) and in all random samples at orders
11–12. The equality case: in `R_n` (= `TT_n` with the arc `1 → n` reversed, `h = 2^{n−2} + 1`) the vertex ending
`(h − 1)/2 = 2^{n−3}` paths is the **penultimate** vertex `n − 1`; the last vertex `n` (the one carrying the reversed arc)
ends exactly **one** path (`auditC_formulas.out`: end counts `[1, 1, 2, 4, …, 2^{n−3}, 1]`). The note's "vertex 5 of
`TT_5 + u`" is right; HYP-9469 (4)'s parenthetical is not (C9). **Pre-registered tests:** (1) `k_min(h) ≥ h/4` for
`h ≤ 60` with ratio in `[0.25, 0.45]`: PASS (ratios in `[0.264, 0.40]`), and in fact for all `h` at order 10; (2)
increment-14 and -16 levels at orders ≤ 9 (and 10): `{1, 3, 5, 9, 11, 15, 17, 19, 23}` and `{1, 9, 13, 15, 17, 19, 23, 29}`,
finite with `h_max = 23 ≤ 28` and `29 ≤ 32`: PASS (levels `k = 9..12` have `h_max = 33, 33, 39, 37`, all `≤ 4k`); (4)
order 9: maximum `64/129` at the Fermat-type count 129: PASS. **New observation (NUMERICAL, orders 4–10):** the table
entries `k_min = 2, 3, 5, 9, 17, 33` at `h = 5, 9, 17, 33, 65, 129` are exactly the sole-exit cuts at the penultimate
vertex of `R_n`, where `e_a = (h − 1)/2` and the digit is `σ = 1`, so by THM-4862 (a) the increment is
`(h + 1 − e_a)/2 = 2^{n−3} + 1`; hence `k_min(2^k + 1) ≤ 2^{k−2} + 1` whenever `σ(R_{k+2}, penultimate) = 1` (verified
`n ≤ 10`), and the conjectured constant `1/4` in (1) would be asymptotically sharp.

**(h) The note.** The reading of "13 and 17" is declared and the alternatives are recorded; no issue. Inheritance:
THM-4094's reduction (lanes `p ≠ 7` and `7p`, the file adds `p ≠ 3`), HYP-9453 (near-join, `B`-sets, strong singletons
`3, 5, 11, 13`), THM-4744 (3) (`P_d = 0` for odd `d < 625` exactly at `7, 21`), THM-4861 (17 at least order 6),
THM-4862 (c) (digit-1 floors `5, 9, 17, 33, …`), the digital bank (`≤ 99,999,999`, audited to `999,999`, as stated in
HYP-9453), MISTAKE-559: all cited at the strength of the files. "Corrected near miss": real — the discarded quantity
`Σ_P (insertion slots)` is `H(X + v) − #orphan paths` (THM-012b), it differs from `H(X + v)` in 355 of the 468 `(X, cut)`
pairs of orders ≤ 5 (gap up to 27), and nothing in the committed files uses it (`extension_spectra.py` is labelled
discarded). Section 4 is sober ("It does not prove the conjecture"); the overclaims are the novelty statements of section
1 (C1), the mod-4/Lass attributions (C4, C5) and the window sentence (C2).

**(i) Priority.** `git grep` over `01-canon` and `05-knowledge/hypotheses`:
- *Flip decomposition:* THM-013 "Significance" (2026-03-05) states `adj(i,j) − adj'(j,i) = Σ 2^k Δ(α_k)`, i.e. `ΔH` is
  the difference of the numbers of paths with `i, j` consecutive; THM-014 writes those counts as subset convolutions of
  `h_end × h_start` (`nadj_T(x,i,j) = Σ_S h_end(T[S+{x}], x) · h_start(T[{j}+R], j)`); THM-017 "Significance" and
  MISTAKE-008 (2026-03-05) write `Δ(S, R) = L_j(S) R_i(R) − L_i(S) R_j(R)` and `δ_H = Σ_S Δ(S, R)` — which is
  `H(T') − H(T) = U(L_j, R_i) − U(L_i, R_j)` verbatim. *`ΔH = ΔI(Ω, 2)` for all `n`:* PROP-001 (kind-pasteur S12,
  2026-03-05, "PROVED (corollary of OCF / Grinberg-Stanley)") and THM-002/CONJ-001 (same session). So THM-4876 (4) has
  no new content; THM-013's status line was merely stale.
- *Increment-2 lemma, increment-4 structure, extension-spectrum census, `k_min`:* no prior occurrence found (THM-008b
  gives `H(T) − H(T − v) = 2(N_3 + N_5)` at `n ≤ 5` only). New.
- *One-centre closed form, `(2^i − 1)(2^l − 1) + 2`:* no prior occurrence found; the family itself is the session124
  "transitive parent" spine (`spine_numeration.py` re-derives its transfer) and the `A = {1}` member is THM-1370's
  `h(R_n) = 2^{n−2} + 1`. The closed form is new.
- *End-count bound `e_a ≤ (h − 1)/2`:* no prior statement beyond THM-4873 (i). The sole-exit increment formula is
  THM-4862 (a) (klein, 2026-11-04), which HYP-9469 misquotes (C8).

## 2. The proof of statement 2, step by step

All steps are correct: the lemma (three cases; `k ≥ 4` even gives the odd cycle `v, a_2, …, a_{k−1}` of length `k − 1`,
or one of the two triangles); "exactly one odd cycle, a triangle, transitive complement" from Claim A with positive odd
terms summing to 1; "the only arc from `A` to `B` is `a → b`" (any `x → y` closes a second triangle through `v`); the
parity of the closed cycle (`1 + j + j'` vertices); the exclusion of `(j, j') = (2, 2)`; converse symmetry (available,
see (b)); with `a → A'` every `A`-part is `(a)`, so odd `j' ≥ 3` is excluded, so an out-neighbour of `b` in `B'` is the
sink `s` of the transitive `B'`; the three-case enumeration of the cycles of `T − v` (the `A`-side is transitive, the
`B`-side acyclic because `s` has no out-neighbour in `B'`, a mixed cycle re-enters the `A`-side at `a` and its `B`-part is
`(b, s)`); hence `H(T − v) = 1 + 2[b → s]`. The exhaustive data confirm each intermediate claim through order 10.

## 3. Corrections (quoted sentence → repair)

**C1 (major, typing/priority) — THM-4876 statement 4, status line, note §1, THM-013 upgrade.** "Flip increment (PROVED)
and THM-013 for all orders" / status "PROVED … 4" / note: "THM-4876 (4) shows the unsigned convolutions give the flip
increment … The order is now reversed: OCF is a theorem … so THM-013 holds for all `n` (THM-4876 (4))" / THM-013 header:
"**PROVED for all n on 2026-11-04**". → Retype statement 4 as **KNOWN (repository, 2026-03-05)**: the unsigned
decomposition is THM-013 "Significance" (`adj(i,j) − adj'(j,i)`), THM-014 (consecutive-pair counts as `h_end × h_start`
subset convolutions) and THM-017/MISTAKE-008 (`δ_H = Σ_S Δ(S,R)`, `Δ(S,R) = L_j(S)R_i(R) − L_i(S)R_j(R)`); the all-`n`
identity `ΔH = ΔI(Ω, 2)` is PROP-001 (and THM-002/CONJ-001). THM-013's header should read "PROVED for all `n` as a
corollary of OCF (THM-002, PROP-001, 2026-03-05); cross-referenced 2026-11-04", and THM-017's new section should say
"closed by citation in THM-002/CONJ-001 (2026-03-05)". The note's §1 should credit these files instead of THM-4876 (4).

**C2 (major, false as written) — THM-4876 statement 5.** "The longest run of consecutive odd values inside `E(X)` has
length at most `24` at order `8` and is never an interval containing `[2h + 1, 3h]`" (note §2: "no `E(X)` covers
`[2h + 1, 3h]`"). → At every order 4–8 the realisers of `h = 5` (the strong 4-tournament, `E = {5, 9, 11, 13, 15}`, and its
paddings: 1, 2, 3, 4, 5 classes at orders 4–8) satisfy `odd[2h + 1, 3h] = {11, 13, 15} ⊆ E(X)`, and `h = 1` does trivially.
Repair: "for `h ∉ {1, 5}` no tournament of order ≤ 8 has `odd[2h + 1, 3h] ⊆ E(X)`". The session's "longest window"
criterion is tie-break dependent (its own output has nonzero counts at orders 3–7); use direct containment
(`audit_C/auditC_window_direct.out`).

**C3 (minor, false as written) — THM-4876 statement 5.** "among the `6880` tournaments of order `8` none has
`E(X) ⊇ odd[h, 2h − 1]`". → Seven do: `TT_8` (`h = 1`) and the six classes with `h = 3` (`5 ∈ E`). Repair: "none with
`h ≥ 5`" (which is statement 2's `h + 2 ∉ E(X)`, PROVED for all orders).

**C4 (minor, attribution) — THM-4876 Setting/status and THM-013 upgrade.** "The *odd-cycle formula* (Lass;
Grinberg-Stanley)" / "also Lass". → Not supported by CORE-PAPERS (which lists Lass only among the authors consolidated by
Irving–Omar) nor by Grinberg–Stanley's paper, which cites `[Lass02]` only for another proof of Berge's theorem (their
§6) and states the odd-cycle formula as their Theorem 1.39 (with Theorem 7.1 the mod-4 refinement). Repair: cite
Grinberg–Stanley arXiv:2307.05569, Theorem 1.39 (and Irving–Omar Corollary 20 as the restatement, per THM-002); drop
"Lass" or cite Lass only for Berge's theorem.

**C5 (minor, attribution) — note §1.** "Irving-Omar arXiv:2412.10572 for the mod-4 refinement, per CORE-PAPERS". →
CORE-PAPERS attributes the mod-4 refinement to Grinberg–Stanley ("yields the odd-cycle-collection formula … plus its
mod-4 refinement"; Theorem 7.1 of 2307.05569); Irving–Omar is the matrix-algebra consolidation. Repair accordingly.

**C6 (minor) — THM-4876 proof 2 and note §0.** "(realised by … the strong 4-tournament, `v` its vertex of out-degree 1 or
in-degree 1: `5 − 3`)" / "whose two vertices of out-/in-degree one give increment 2". → The strong 4-tournament has score
sequence `(1, 1, 2, 2)`, so every vertex has out-degree 1 or in-degree 1; only the two vertices lying on exactly one of
its two triangles give `5 − 3` (one of out-degree 1, one of in-degree 1); the other two lie on both triangles and give
`5 − 1 = 4` (`auditC_misc.out`). Repair: "`v` either vertex that lies on only one of its two triangles".

**C7 (minor, incomplete argument) — THM-4876 proof 6.** "`N(A) = 3` would be three adjacent `10`-pairs … or a pair at
distance `3` plus an adjacent pair …". → Pairs at distance 2 also have `c = 1`, so the case split is incomplete. Repair:
let `a` be the first `1` and `b` the last `0` with `a < b`; `b − a ≥ 4` gives `c(b − a) ≥ 4`; `b − a = 3` gives `N ≥ 4` in
all four fillings of the two middle letters (`c(3) = 2` plus at least two more pairs, or `c(3) = 2` plus `(a+1, b)` and
one of `(a, a+2)`, `(a+2, b)`); `b − a ≤ 2` confines all pairs to the window `[a, b]`, whose `N` is `1`, `2` or `2`
(`10`, `100`, `110`). So `N ∈ {0, 1, 2} ∪ [4, ∞)`. (Verified: attained `N < 12` are `0, 1, 2, 4, 5, 7, 8, 9, 11`.)

**C8 (minor, wrong formula) — HYP-9469 (1) and (4).** "The sole-exit cut gives the increment `(h − e_a)/2`" / "For the
sole-exit cut the increment is `(h − e_a)/2`". → By THM-4862 (a) the one-exit extension has `h(T + x) = 2h + σ − e_a`,
so the increment is `(h + σ − e_a)/2` with `σ = B_a − h ≥ 0`; equality with `(h − e_a)/2` only when `σ = 0` (53 of the 420
`(T, a)` of orders 2–6; e.g. the strong 4-tournament at its `e_a = 2` vertex has increment 2, not 1.5). The deductions
"(1) would follow from (4) for those cuts" and "`k_min(h) ≥ (h + 1)/4` on sole-exit cuts" survive as inequalities.

**C9 (minor, wrong parametrisation) — HYP-9469 (4).** "equality holds in the reversed-extreme-arc family (`h = 2^k + 1`,
`e_a = 2^(k−1)` at the vertex `k + 1` of `TT_(k+1)` with the arc `(1, k+1)` reversed)". → `TT_{k+1}` with `(1, k+1)`
reversed has `h = 2^{k−1} + 1`, and its vertex `k + 1` ends exactly one path. Repair: "`h = 2^k + 1` in `TT_{k+2}` with
the arc `(1, k+2)` reversed; equality `e = 2^{k−1}` at the penultimate vertex `k + 1`" (the note's "vertex 5 of
`TT_5 + u`, `m = 5`" is the correct instance).

**C10 (minor, self-contradiction) — HYP-9469 title.** "and in general `h <= 3k − 1` (?)". → Contradicted by the listed
levels: `k = 3: 9 > 8`, `k = 5: 17 > 14`, `k = 6: 19 > 17` (and `k = 7..12: 23, 29, 33, 33, 39, 37` versus `20, 23, 26,
29, 32, 35`). The body's "at most `4k`" is consistent with all data through order 10. Delete the `3k − 1` clause.

**C11 (minor, wording) — HYP-9469 (2).** "(stable from order `6` on)" in a hypothesis whose census is described by
"orders of `T` up to `9`". → At `T`-order 6 the `k = 5, 6` sets are still incomplete (`{1,3,5,9,13}`,
`{3,5,11,13,15}`); they are complete from `T`-order 7 (= `X`-order 6, the `n` of `endcount_ratio.py`). Say "from order 7
of `T` (order 6 of `X`) on". Also commit an `endcount_ratio.py` output: (2) and (4) quote numbers that no committed
`.out` supports (reproduced here in `session_rerun_endcount_ratio_8.out`).

**C12 (minor, overclaim) — THM-4876 statement 5 and HYP-9469 (5).** "Hence the two natural single-realiser inductions
… are impossible" / "THM-4876 (5) rules out every induction that extends one realiser by one vertex". → The
`[h, 2h − 1]` induction is impossible for all orders by statement 2 (PROVED: `h + 2 ∉ E(X)` for `h ≥ 5`); the window
induction `[ch, c'h]` is only excluded through order 8 (FINITE-EXACT, and not even there for `h = 5`, C2). "Every
induction" is not established. Repair the typing: "impossible (first kind, PROVED); not available through order 8
(second kind, FINITE-EXACT)".

**C13 (minor, typing) — HYP-9469 status.** "the only proved cases are `k_min(h) ≥ 2` for `h ≥ 5` (THM-4876 (2)) and
`k_min(h) ≥ 3` for `h ∉ {1, 5}` (THM-4876 (3), census)". → The second is FINITE-EXACT (orders of `T` ≤ 10 after this
audit), not proved; move it out of "proved cases".

**C14 (cosmetic) — THM-4876 statement 5.** "about `h/3.5`". → Over the table `h/k_min` ranges from `2.5` to `3.67`,
mean `3.1`; "between `h/3.7` and `h/2.5`" is accurate, and the minimum ratio decreases towards `1/4` along
`h = 2^k + 1` (see (g)).

## 4. What the audit adds (all FINITE-EXACT / NUMERICAL, from `audit_C/`)

- The increment-2 theorem and the increment-4 structure verified exhaustively at **order 10** (`del10_merged.out`,
  `auditC_structure.out`), plus random orders 10–12.
- `k_min(h) ≥ h/4` for **every** `h` attained at order 10 (`h ≤ 3357`), minimum `33/129`; the `k_min` values at
  `h = 2^k + 1` are the sole-exit cuts of `R_{k+2}` with digit `σ = 1`, so `k_min(2^k + 1) ≤ 2^{k−2} + 1` (orders 4–10).
- Levels `k = 7..12` at orders ≤ 10: `h_max = 23, 29, 33, 33, 39, 37` (pre-registered test (2) PASS).
- `e_a ≤ (h − 1)/2` for non-sinks: exhaustive through order 10 (pre-registered test (4) PASS at order 9 with `64/129`;
  order 10 gives `128/257`).

## 5. Files

`04-computation/experiments/chanpak_spine_20261104/audit_C/`: `auditC_census.c` (+ binary `auditC_census`),
`auditC_lib.py`, `auditC_formulas.py` → `auditC_formulas.out`, `auditC_structure.py` → `auditC_structure.out`,
`auditC_ext_post.py` → `auditC_ext_post.out`, `auditC_window_direct.out`, `auditC_misc.out`, `del3.out` … `del9.out`,
`del10_part0..5.out`, `del10.out` (increment-2/4 pairs of order 10), `del10_merged.out`, `ext3.out` … `ext8.out`,
`rand.out`, `session_rerun_{ext_spectrum_true,ext_windows,increment_ladder,endcount_ratio}_8.out`.

*Housekeeping (session, 2026-11-04): `ext8.out` (6.6 MB) and the compiled binary `auditC_census` are not committed; `auditC_census.c` regenerates them. All fourteen corrections were applied in the next commit (MISTAKE-715).*
