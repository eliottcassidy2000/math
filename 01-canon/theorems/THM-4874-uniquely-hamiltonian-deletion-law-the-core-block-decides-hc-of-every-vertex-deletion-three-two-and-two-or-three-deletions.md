---
id: THM-4874
title: "The uniquely Hamiltonian deletion law. Let T have no vertex of out-degree 1 and exactly k Hamiltonian cycles, with k = 3 and n >= 6, or k = 4, 5 and n >= 7, and let (W, A) be the block of its terminal core (THM-4869). Then hc(T - v) is fixed by where v sits: a multiple of k in the prefix, 0 at the tail u, 0 at the head l (unless the prefix is the single vertex l), k on the chain, and delta(v) = sum over x in A - v of end_x(W - v) in the block. Hence T - v is uniquely Hamiltonian iff v is a block vertex with delta(v) = 1. There are exactly 3 such v for k = 3 and exactly 2 for k = 4. For k = 5 there are 3 if the block is (AT_5, its two sole exits) and 2 otherwise, so F(2n - 13) of the u_5(n) classes have three and L(2n - 12) have two. The block also fixes #{v : hc(T - v) = 2}; with THM-4872 (B), a three-cycle sole-exit-free tournament has a nontrivial automorphism iff no deletion leaves exactly two cycles. Every core block (AT_m, {m-1, m}) gives exactly three deletions, at 1, m-1, m"
status: >
  PROVED (kind-pasteur-2026-11-04-S10), given THM-4869 (C2), Lemmas B1-B5 and (C4) (opus-S75; its (C1) base, orders
  6-12, is FINITE-EXACT). The proof is the one-arc product law applied to T - v, plus a finite table for the eight
  blocks of type <= 5 (computed, and by hand for C3, W4 and AT_5). Statement (5) holds for every m >= 4 and needs no
  census.
  PRE-REGISTERED (HYP-9467, commit 671f7fb87e, before the run): the joint counts of (#{hc(T - v) = 1},
  #{hc(T - v) = 2}, |Aut|) at orders 10-15, run on the HYP-9464 auditor's independent census. RESULT: CONFIRMED, all
  18 predicted rows exact and no other combination (order 15: 1,090,032 strong classes with hc <= 5).
  Proves HYP-9464 (B), which was FINITE-EXACT at orders 6-14, and the HYP-9464 audit's hc = 5 observation.
  AUDIT: pending.
session: kind-pasteur-2026-11-04-S10
source: 05-knowledge/results/fermat_self_duality_fewcycle_coincidences_20261104.md
scripts:
  - 04-computation/experiments/fewcycle_block_deletions_20261104.py (+ .out)
  - 04-computation/experiments/fewcycle_deletion_positions_20261104.py (+ .out)
  - 04-computation/experiments/fewcycle_deletion_census_20261104.c (+ _explore.out, .out)
  - 04-computation/experiments/at_block_fermat_cycles_20261104.py (+ .out)
related:
  - THM-4869
  - THM-4872
  - THM-4868
  - HYP-9464
  - HYP-9467
  - HYP-9459
---

# THM-4874: the uniquely Hamiltonian deletion law

## Setting (THM-4869, opus-S75)

Let T have order n, no vertex of out-degree 1 ("sole-exit-free"), and exactly k Hamiltonian cycles, where k = 3 and
n >= 6, or k in {4, 5} and n >= 7. THM-4869 (C2), with B1-B4, gives this structure:
- T has arc-strong connectivity 1 and a unique terminal core (D, u). D is the minimal set with exactly one leaving arc
  u -> l, and e_u(D) = k, where e_u(D) is the number of Hamiltonian paths of T[D] ending at u.
- Put R = V(T) - D (the prefix). Either R = {l}, or the contraction Q = T/D is uniquely Hamiltonian with the
  contracted vertex z as a sole exit. In both cases s_l(R) = 1, where s_l counts Hamiltonian paths of T[R] starting at
  l.
- Every arc between R and D points into D, except u -> l.
- D - u = w_1 => ... => w_p => W: a transitive chain (p >= 0) followed by a block (W, A) of type k. Here
  A = N-(u) ∩ W, u beats W - A, u -> w_1 if p >= 1, and p = 0 is allowed only if A != W.

Notation for the statements:
- delta(v) = sum over x in A - v of end_x(W - v), for v in W.
- h0(W, A) = #{Hamiltonian paths of W from a vertex of W - A to a vertex of A}.
- F is Fibonacci and L is Lucas.

## Statement

**(1) The deletion table.** For every vertex v of T:

| v | hc(T - v) |
|---|---|
| v in R, v != l | k * s_l(R - v), a multiple of k |
| v = u | 0 |
| v = l | 0 if R != {l}; if R = {l}: k when p >= 1, h0(W, A) when p = 0 |
| v = w_i (chain) | k |
| v in W | delta(v) |

**(2) The deletion law.** T - v is uniquely Hamiltonian iff v ∈ W and delta(v) = 1. The block data are:

| k | block (W, A) | delta over W | #(delta = 1) | #(delta = 2) | h0 |
|---|---|---|---|---|---|
| 3 | (C3, all) | 1, 1, 1 | 3 | 0 | (p >= 1) |
| 3 | (W4, {c, d}) | a 1, b 2, c 1, d 1 | 3 | 1 | 2 |
| 4 | (W4, {a, c, d}) | a 1, b 3, c 2, d 1 | 2 | 1 | 2 |
| 4 | (W4, {b, c, d}) | a 1, b 2, c 2, d 1 | 2 | 2 | 1 |
| 5 | (W4, all) | a 1, b 3, c 3, d 1 | 2 | 0 | (p >= 1) |
| 5 | (AT_5, {4, 5}) | 1 at 1, 4, 5; 3 at 2, 3 | 3 | 0 | 4 |
| 5 | (W_5, A_1) | 1, 1, 2, 3, 3 | 2 | 1 | 4 |
| 5 | (W_5, A_2) | 1, 1, 2, 3, 2 | 2 | 2 | 4 |

Here W4 = AT_4 is a -> b -> c -> d -> a with a -> c, b -> d. W_5 is the uniquely Hamiltonian 5-tournament with scores
1, 2, 2, 2, 3 and one sole exit (`fewcycle_block_deletions_20261104.out` gives its rows and the two sets A_1, A_2).

So T has exactly:
- 3 uniquely Hamiltonian deletions if k = 3;
- 2 if k = 4;
- 3 if k = 5 and the block is (AT_5, its two sole exits), and 2 for the other three type-5 blocks.

**(3) Counts.** The classes with a given block number:
- F(2n - 2w - 4) for a block of order w with A = W: (C3, all) gives F(2n - 10), and (W4, all) gives F(2n - 12);
- F(2n - 2w - 3) for a block with A != W: F(2n - 11) for each W4 block of type 3 or 4, and F(2n - 13) for each
  5-block.

So for k = 5:
- F(2n - 13) classes have three uniquely Hamiltonian deletions;
- F(2n - 12) + 2F(2n - 13) = L(2n - 12) have two.

Their total is u_5(n) = 4F(2n - 12) - 3F(2n - 14) (THM-4869 (C4)).

**(4) The value 2.** For n >= 7, k2 = #{v : hc(T - v) = 2} is the block's #(delta = 2):
- k = 3: 0 (C3 block) or 1 (W4 block);
- k = 4: 1 or 2;
- k = 5: 0, 0, 1, 2 in the order of the table.

With THM-4872 (B) (opus-S76: Aut(T) = Z_3 iff the block is C3, and trivial otherwise): **a sole-exit-free tournament
with exactly three Hamiltonian cycles has a nontrivial automorphism iff no vertex deletion leaves exactly two
Hamiltonian cycles.** This holds for n >= 6; at n = 6 the l-deletion adds a 2 only to the W4-block class.

**(5) Almost transitive blocks.** For every m >= 4, (AT_m, {m - 1, m}) is a block of type 2^(m-3) + 1. This is the
note's section 5 (a), checked for m = 4..12 by `at_block_fermat_cycles_20261104.py`.

Let T be any one-arc join with this block, any chain, and a uniquely Hamiltonian or one-vertex prefix. Then:
- T has 2^(m-3) + 1 Hamiltonian cycles;
- T has exactly three uniquely Hamiltonian deletions, the block vertices 1, m - 1 and m;
- each other block vertex j gives hc(T - j) = H(AT_(m-2)) = 2^(m-4) + 1, and h0 = 2^(m-3).

The type-3 W4 block (m = 4) and the 3-deletion type-5 block (m = 5) are the first two cases.

## Proof

**Lemma (the product law under deletion).** Let v ∉ {u, l}, and put D' = D - v and R' = R - v. Then

    hc(T - v) = e_u(D') * s_l(R').

*Proof.*
- **Every cycle splits.** In T - v the only arc leaving D' is u -> l. So every Hamiltonian cycle crosses from D' to R'
  exactly once. It is a Hamiltonian path P of T[D'] ending at u, then u -> l, then a Hamiltonian path P' of T[R']
  starting at l, then a closing arc end(P') -> start(P).
- **Every pair closes.** D has at least 3 vertices (B1), so |D'| >= 2 and start(P) != u.
  - If end(P') != l, it beats all of D.
  - If end(P') = l, then R' = {l}, and l beats all of D - u, since u -> l is the only arc from D to l.
- So pairs (P, P') and cycles correspond bijectively. QED.

**The five cases of (1).**
- **v ∈ R - {l}.** D' = D, so hc(T - v) = e_u(D) s_l(R - v) = k s_l(R - v).
- **v = u.** D - u is non-empty and sends no arc to R (u -> l was the only one), while R is non-empty. So T - u is not
  strong.
- **v = l.**
  - If R != {l}: D sends no arc to R - l, which is non-empty, so T - l is not strong.
  - If R = {l}: T - l = T[D]. Cutting each Hamiltonian cycle of D at u gives a Hamiltonian path of D - u from an
    out-neighbour of u to an in-neighbour of u.
  - The condensation D - u = w_1 => ... => w_p => W forces these paths to be w_1 ... w_p P_W, with P_W a Hamiltonian
    path of W.
  - For p >= 1 the start w_1 is an out-neighbour of u, so the count is the number of P_W ending in A, that is
    e_u(D) = k.
  - For p = 0 the count is h0(W, A).
- **v = w_i.** By the lemma, hc(T - w_i) = e_u(D - w_i) s_l(R) = e_u(D - w_i). The tournament D - w_i - u still has
  condensation (chain - w_i) => W. So its Hamiltonian paths end in W, and those ending at in-neighbours of u number
  sum over x in A of end_x(W) = k.
- **v ∈ W.** By the lemma, hc(T - v) = e_u(D - v). Every chain vertex beats all of W, so every Hamiltonian path of
  D - v - u is the chain in order followed by a Hamiltonian path of W - v. It must end in N-(u) ∩ (W - v) = A - v,
  which gives delta(v).

**(2).**
- **Where the values 1 and 2 can occur.** By (1), values in {1, 2} occur only for v ∈ W, or for v = l with R = {l}
  and p = 0. The latter has n = |W| + 2.
- **Ruling out the l-deletion.**
  - k = 3: the C3 block has A = W, so p >= 1; the W4 block gives h0 = 2 at n = 6.
  - k = 4: n >= 7 > |W| + 2 = 6.
  - k = 5: (W4, all) has A = W, and each 5-block gives h0 = 4 at n = 7.
- So in the theorem's range the uniquely Hamiltonian deletions are exactly the v ∈ W with delta(v) = 1.
- **The block table.** The blocks of type <= 5 are THM-4869 B5's eight, re-derived by
  `fewcycle_block_deletions_20261104.py` over all labelled tournaments of order <= 6. The table above is that
  program's output.
  - By hand for W4: W4 - a = (b c d) and W4 - d = (a b c) are transitive with ends in A; W4 - b = (a c d) and
    W4 - c = (a b d) are 3-cycles.
  - For AT_5, see the proof of (5).

**(3).** THM-4869 B3-B4 and (C4) give the classes with a given block as follows.
- **Cores.** For a block of order w there are 2^(d-w-2) cores of order d >= w + 2: chain length p = d - w - 1 >= 1,
  with 2^(p-1) choices of the free arcs. When A != W there is also one core of order w + 1 (p = 0). Distinct choices
  give distinct cores, and classes correspond to (prefix, core) pairs.
- **Generating functions.** With the rooted-prefix factor x(1 - 2x)/(1 - 3x + x^2):
  - A = W: x^(w+2)/(1 - 2x) * x(1 - 2x)/(1 - 3x + x^2) = x^(w+3)/(1 - 3x + x^2), so the count is F(2n - 2w - 4);
  - A != W: x^(w+1)(1 - x)/(1 - 2x) * ... = x^(w+2)(1 - x)/(1 - 3x + x^2), so the count is
    F(2n - 2w - 2) - F(2n - 2w - 4) = F(2n - 2w - 3).
- **Check.** Summing over the blocks of each type recovers THM-4869's u_3, u_4 and u_5. Finally,
  L(2n - 12) = F(2n - 13) + F(2n - 11) = F(2n - 12) + 2F(2n - 13).

**(4).** By the same case analysis, the value 2 occurs only in W, or at v = l in the minimal case. For n >= 7 the
minimal case gives h0 ∈ {4} or is absent. For k = 3 at n = 6 it gives the W4 block one more 2, which does not affect
"k2 = 0 iff C3 block".

**(5).** Label AT_m as 1 < ... < m with m -> 1.
- **It is a block.** This was shown in the note's section 5 (a), together with the type 2^(m-3) + 1.
- **The cycle count.** The lemma and the cases of (1) apply unchanged (they used only the one-arc structure, the
  condensation of D - u and s_l(R) = 1). So hc(T) = e_u(D) = 2^(m-3) + 1.
- **The block values.**
  - W - 1 and W - m are transitive, with ends m ∈ A and m - 1 ∈ A: delta = 1.
  - W - (m - 1) is AT_(m-1) and A - (m - 1) = {m}. Only its transitive path ends at m: delta = 1.
  - For 2 <= j <= m - 2, W - j is AT_(m-1) with extremes 1 and m, in which m - 1 is the second-to-last vertex. So
    delta(j) = end_m + end_(m-1) = 1 + 2^(m-4).
- **The minimal case.** A path from W - A = {1, ..., m - 2} into A is either the transitive path, or A' m 1 B' with
  A' a non-empty subset of {2, ..., m - 2} and m - 1 ∈ B'. So h0 = 1 + (2^(m-3) - 1) = 2^(m-3), which is >= 2 for
  m >= 4.
- **The other positions.** Chain deletions give k, prefix deletions multiples of k, u gives 0, and l gives 0, k or h0.
  None of these is 1. QED.

## Verification

- `fewcycle_block_deletions_20261104.py`: the eight blocks, delta, #(delta = 1), #(delta = 2) and h0, and THM-4869
  B5's list re-derived (type 2: none; types 3, 4, 5: 2, 2 and 4 blocks).
- `fewcycle_deletion_positions_20261104.py`: every vertex deletion of every one-arc join built from the eight blocks
  checked against the table in (1). It uses the chain lengths 0-2 with all free arcs, and every rooted uniquely
  Hamiltonian prefix of order <= 5 up to isomorphism, plus the one-vertex prefix.
- `at_block_fermat_cycles_20261104.py`: the AT_m block conditions and hc for m = 4..12.
- **Census.** The HYP-9464 audit's independent census (`aud_fewcycle.c`, every strong class with hc <= 5) gives the
  distributions of #{v : hc(T - v) = 1}, and they agree with (2) and (3) at orders 6-14. That census plus the joint
  statistic is the pre-registered test HYP-9467 at orders 10-15, CONFIRMED (`fewcycle_deletion_census_20261104.out`).
  At order 15 the joint counts were:
  - hc 3: (3,0,3) 6,765 and (3,1,1) 4,181;
  - hc 4: (2,1,1) and (2,2,1) 4,181 each;
  - hc 5: (2,0,1) 2,584, and (3,0,1), (2,1,1), (2,2,1) 1,597 each.

## Credits

- The structure theory (cores, blocks, prefixes, counts) is opus-S75's THM-4869.
- The Z_3 / asymmetric split is opus-S76's THM-4872 (B).
- The few-cycle counts are klein's HYP-9459.
- The census program and the hc = 5 observation (L(2n - 12) / F(2n - 13)) are the HYP-9464 auditor's.
- New here:
  - the deletion table (1) and the law (2)-(4);
  - the identification of the 3-deletion type-5 block as (AT_5, {4, 5});
  - the AT_m block family (5).

## Audit

Pending.
