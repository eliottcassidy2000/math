# Independent audit: THM-4874 (the uniquely Hamiltonian deletion law) and HYP-9467

**Auditor:** an independent general-purpose subagent of kind-pasteur-2026-11-04-S10. It was briefed only with the
files to audit and instructions to re-derive and break every claim with its own code.
- It audited HEAD 460c918f9a (the pre-audit checkpoint with the order-15 result) plus the uncommitted edits to HYP-9464
  and the fermat note.
- It wrote no `.md` file, so the session transcribed this report from its final message.
- The programs and outputs in this folder are the auditor's own. They share no code with the session or with the
  HYP-9464 auditor.

**Verdict: PASS WITH CORRECTIONS.** 0 major, 11 minor. The mathematics is correct.
- The Setting follows from THM-4869 (C2) with B1-B4 in the stated ranges (k = 3 from n = 6; k = 4, 5 from n = 7).
- The deletion lemma hc(T - v) = e_u(D - v) s_l(R - v) is correct. It needs no strongness, and the closing arc when
  R' = {l} is handled.
- All five cases of (1), and (2)-(5), check out.
- HYP-9467's predictions match the auditor's own census exactly.

All corrections were applied in the commit that adds this file (MISTAKE-714).

## Findings (all MINOR)

1. **Title, last sentence.** "Every core block (AT_m, {m-1, m}) gives exactly three deletions" is false without the
   uniquely Hamiltonian prefix that (5) assumes. With a rooted prefix of hc 2, the join has hc = 2(2^(m-3) + 1) and no
   uniquely Hamiltonian deletion: m = 4, 5, 6 give 0 (`aud_at_title_prefix.out`). The prefix hypothesis was added.
2. **(5), third bullet.** "H(AT_(m-2))" is undefined at m = 4. It now reads 2^(m-4) + 1, = H(AT_(m-2)) for m >= 5.
3. **Title.** "The block also fixes #{v : hc(T - v) = 2}" holds only for n >= 7. At n = 6 the W4-block class gets an
   extra 2 from the l-deletion. "For n >= 7" was added.
4. **Status, dependencies.** (3) also uses THM-4868 (F), and (4) uses THM-4872 (B), which uses THM-4868 (E). Added.
5. **"18 predicted rows".** There are 48 predicted counts in 18 JOINT lines. The blind range is now stated: everything
   is blind at order 15; at orders 10-14 only k2 and the joint correlation are blind.
6. **The census claim's range.** "Agree at orders 6-14" now reads 6-14 for hc = 3 and 7-14 for hc = 4, 5. Order 6
   has 2-strong classes, and (W4, {b, c, d}) with p = 0 gives h0 = 1.
7. **Credit for (5).** The block property and type of the AT_m blocks are the note's section 5 (a), confirmed by the
   HYP-9464 audit (its F6). Only their deletion data are new here.
8. **Credit for the counts.** klein conjectured u_3 and u_4 (HYP-9459); they and the new u_5 law are proved in
   THM-4869 (C4). klein's related hc = 2 segment-deletion remark is now cited as context.
9. **(5) proof.** hc(T) = s_l(R) e_u(D) is B2's product law; the Lemma only covers v ∉ {u, l}. Reworded.
10. **W_5.** The sets A_1 and A_2 were defined only in the program output. An intrinsic description was added: the
    unique cycle x -> y -> t -> s -> r -> x with chords y->s, t->x, t->r, s->x, r->y, A_1 = {x, y} and A_2 = {x, t}.
    The auditor checked it against the program labels.
11. **Companion files.**
    - The fermat note's header said "No new canon theorem".
    - `05-knowledge/results/INDEX.md` had no entry for the note.
    - HYP-9464's edit needed "blind at order 15" and the range "n >= 6 for hc = 3, n >= 7 for hc = 4" in place of
      "all n".

    All are fixed.

## Novelty and credit

- No prior or concurrent statement of the deletion law, the joint (k1, k2, |Aut|) statistic, or the AT_m deletion
  data was found. The search covered 01-canon, 05-knowledge, the S75 and HYP-9464 audits, and commits on all branches
  since 2026-11-04.
- Checked and found unrelated:
  - Session 190/191 "deletion growth", which concerns path counts;
  - S76 and S77;
  - S78's THM-4875;
  - mac-mini's THM-4876 / HYP-9469 (path-count increments).
- There are no numbering collisions.
- The census code is `aud_fewcycle.c`, changed only by NM = 15 and the joint statistic (checked by diff). Its
  completeness argument holds:
  - Moon 1966;
  - insertion at a switch is injective, so hc(T - v) <= hc(T);
  - the min-over-(cycle, start) canonical form is valid, with |Aut| = the number of minimisers.
- Pre-registration integrity: commit 671f7fb87e at 17:01:18; run 17:01:31-17:17:35.

## Reproduced independently

- **`aud_blocks.c`** covers all labelled tournaments of orders 3-7 (1,677,488 strong at order 7). There are exactly
  eight blocks of type <= 5, at orders 3-5. Their delta vectors, #(delta = 1), #(delta = 2), h0 and |Aut(W, A)|
  match the table, and the 3-deletion type-5 block is (AT_5, {4, 5}).
- **`aud_uh5_scores.py`**: W_5 is the unique uniquely Hamiltonian 5-tournament with scores 1, 2, 2, 2, 3.
- **`aud_uhdel_census.c`** uses an individualisation-refinement canonical form.
  - The labelled brute force to order 7 and the vertex extension to order 11 agree.
  - t_k and u_k match canon.
  - The joint counts at orders 7-11 match HYP-9467's formulas exactly, for example order 11: 144/89; 89/89;
    55/34/34/34.
- **`aud_structure_check.py`** covers all 917 in-range classes at orders 6-11.
  - The minimal one-arc-leaving set is unique.
  - Every Setting claim holds.
  - Table (1) agrees on all 9,536 vertex deletions, with 0 mismatches and every case exercised:

    | case | deletions |
    |---|---|
    | W | 3,601 |
    | chain | 1,975 |
    | prefix | 2,126 |
    | u | 917 |
    | l with \|R\| >= 2 | 679 |
    | l with R = {l}, p >= 1 | 232 |
    | l with R = {l}, p = 0 | 6 |

  - (4) has no violations at n = 6-11.
  - An injected wrong prediction was caught 6 times out of 6.
- **`aud_at_family.c`** tests 216 joins for m = 4..9: exactly three uniquely Hamiltonian deletions, interior values
  2^(m-4) + 1, h0 = 2^(m-3), and 0 failures. `aud_at_block_core.py` confirms the block conditions and the terminal
  core for m = 4..10.
- **`aud_gf.py`**: the per-block counts hold from n = w + 2 to n = 40, and HYP-9467's rows come out exactly.
