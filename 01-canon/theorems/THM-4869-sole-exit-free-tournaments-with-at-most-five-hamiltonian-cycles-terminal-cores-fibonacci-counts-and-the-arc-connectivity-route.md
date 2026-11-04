---
id: THM-4869
title: "Sole-exit-free tournaments with at most five Hamiltonian cycles. (C1) No 2-strong tournament of order n >= 6 has <= 3 Hamiltonian cycles, and none of order n >= 7 has <= 5. (C2) Every strong tournament of order n >= 6 with hc <= 3 (n >= 7: hc <= 5) and no out-degree-1 vertex is a one-arc join Q[z <- (D, u)] of a uniquely Hamiltonian rooted prefix and a unique terminal core with e_u(D) = hc. (C3) Every strong tournament with minimum out-degree >= 2 and n >= 6 has at least 3 Hamiltonian cycles (klein HYP-9459 (1)). (C4) The sole-exit-free counts are u_3(n) = F(2n - 9) (n >= 6), u_4(n) = 2F(2n - 11) and u_5(n) = 4F(2n - 12) - 3F(2n - 14) (n >= 7), from doubling cores c_3(d) = 3 * 2^(d-6), c_4(d) = 2^(d-5), c_5(d) = 5 * 2^(d-7) whose pole cancels the rooted-prefix factor (1 - 2x). (A) A 2-arc-strong tournament with a cut vertex has hc >= 9; with Fraisse-Thomassen this is an arc-connectivity proof that a uniquely Hamiltonian tournament is a one-arc join"
status: >
  PROVED (opus-2026-11-04-S75):
  - (A), Lemmas B1-B5, (C1)-(C4). (C1) is an induction whose base, orders 6-12, is FINITE-EXACT from a complete census.
    Lemma B5 uses THM-4868 (C), Redei's parity theorem, and finite checks at orders <= 6.
  - The corollary of (A) uses Fraisse-Thomassen 1987 (CITED, not read).
  - (C3) and (C4) settle klein's HYP-9459 (1)-(3); the k = 5 law is new.
  PRE-REGISTERED and CONFIRMED: HYP-9461 P1-P10 at orders 10-14 (c355721ce6, before the run).
  FINITE-EXACT: the minimum hc of 2-strong tournaments (2, 4, 8, 13, 19, 30 at n = 5-10), and the core lineage.
  ANALOGY: the {-3}-self-duality parallel and the inspiration readings (note section 4).
  NOVELTY: no prior source found (web search: Thomassen 1980, Busch 2006, Kuehn-Osthus survey, Demaria-Kiihl,
  Burzio-Demaria; bodies not read).
  AUDIT: see the note, section 7.
session: opus-2026-11-04-S75
source: 05-knowledge/results/arc_connectivity_route_and_sole_exit_free_few_cycle_tournaments_20261104.md
related:
  - THM-4867
  - THM-4868
  - THM-4866
  - THM-4862
  - HYP-9459
  - HYP-9461
  - HYP-9458
  - HYP-9460
verification:
  - 04-computation/experiments/fewcycle_sole_exit_free_20261104.c
  - 04-computation/experiments/fewcycle_sole_exit_free_20261104_14_5.out
  - 04-computation/experiments/fewcycle_sole_exit_free_20261104_mincycles_11_40.out
  - 04-computation/experiments/fewcycle_core_lineage_20261104.c
  - 04-computation/experiments/fewcycle_core_lineage_20261104.out
  - 04-computation/experiments/fewcycle_blocks_20261104.c
  - 04-computation/experiments/fewcycle_blocks_20261104.out
---

# THM-4869 - sole-exit-free tournaments with at most five Hamiltonian cycles

## Notation

- hc is the number of Hamiltonian cycles, and H the number of Hamiltonian paths (odd, by Redei).
- end_x and start_x count Hamiltonian paths ending or starting at x.
- A sole exit has out-degree 1; a sole entry has in-degree 1.
- A one-arc-leaving set D has exactly one arc u -> l leaving it.
- t_k(n) and u_k(n) are klein's counts (HYP-9459): classes of order n with hc = k, and those with no sole exit.
- W4 is the strong 4-tournament a -> b -> c -> d -> a with a -> c and b -> d. Its end counts are (1, 1, 2, 1), and its
  sole exits are c and d.

## (A) Arc-connectivity

Let T be 2-arc-strong with a cut vertex v, and T - v = D_1 => ... => D_m. Then:

    hc(T) = S_+ * prod_{1<i<m} H(D_i) * S_-,

where S_+ = sum over a in N+(v) & D_1 of start_a(D_1), and S_- is the dual sum over D_m.

S_+ >= 3, and dually S_- >= 3, so hc >= 9. The bound is attained by v => C3 => C3 => v:
- C3 forces v to beat all three vertices;
- W4 forces v to beat both sole entries (1 + 2 paths);
- hc(D_1) >= 2 gives at least 4;
- a uniquely Hamiltonian D_1 of order m >= 5 has a non-critical sole entry y (THM-4868 (C)) with
  start_y >= m - 2 >= 3.

**Corollary.** A uniquely Hamiltonian tournament is not 2-arc-strong:
- with a cut vertex, hc >= 4;
- if 2-strong, Fraisse-Thomassen (Graphs Combin. 3 (1987), k = 2) gives a cycle avoiding any arc.

Hence it is a one-arc join, a proof using neither Shelor's grammar nor THM-4868.

**Forced arcs.** In a 2-arc-strong tournament, an arc lying on every Hamiltonian cycle needs a cut vertex (by
Fraisse-Thomassen) and must lie in every factor choice. The S74 audit's four order-9 witnesses are
v => C3 => x => y => C3 => v (1 + 3 + 1 + 1 + 3 = 9), with forced arc x -> y and hc = 9.

## (B) Terminal cores

**B1.** In a strong tournament, distinct minimal one-arc-leaving sets are disjoint, and disjoint one-arc-leaving sets
satisfy |D_1| |D_2| <= 2. If there is no sole exit, the minimal set D is unique, strong, and of order >= 3.

**B2.** hc(T) = s_l(T - D) e_u(D). If T - D has at least 2 vertices, s_l(T - D) = hc(T/D), and the contracted vertex z
is a sole exit of T/D.

**B3.** (D, u) is a **core**:
- no vertex of D - u has out-degree 1 in D;
- no non-empty Z inside D - u sends exactly one arc into D - Z;
- e_u(D) >= 2.

For the last condition: if e_u(D) = 1, then D is uniquely Hamiltonian with u as its only sole exit, and the
non-critical sole exit of THM-4868 (C) gives e_u >= |D| - 2 >= 2.

Every pair (strong Q with sole exit z, core (D, u)) gives a sole-exit-free T = Q[z <- (D, u)] with terminal core D.
Isomorphism classes correspond.

**B4.** Write D - u = W_1 => ... => W_r and A = N-(u) & W_r. Then

    e_u(D) = prod_{i<r} H(W_i) * sum over A of end_x(W_r),

with |A| >= 2. If e_u <= 5, then:
- D - u is a transitive chain w_1 => ... => w_p followed by a block (W, A), with u -> w_1;
- the arcs between u and w_2, ..., w_p are free;
- p = 0 is allowed iff A != W.

**B5.** The blocks of type k <= 5:
- k = 2: none;
- k = 3: (C3, all vertices) and (W4, {c, d});
- k = 4: (W4, {c, d, a}) and (W4, {c, d, b});
- k = 5: (W4, all vertices), and three uniquely Hamiltonian 5-blocks with |A| = 2.

The proof:
- hc(W) >= 3 gives a sum of at least 6.
- hc(W) = 2 and order >= 6: a sole exit z exists by (C3).
  - If z is non-critical, the sum is at least 6.
  - If z is critical, Redei parity forces W - z = Q_1 => singletons with start_(z')(Q_1) = 2 and H(Q_1) = 3, so
    Q_1 = C3, a contradiction.
- A uniquely Hamiltonian W of order m >= 5 gives a sum of at least m - 1.
- The finite cases (orders <= 6) are checked exhaustively.

Hence, for d >= 7:

    c_3(d) = 3 * 2^(d-6),   c_4(d) = 2^(d-5),   c_5(d) = 5 * 2^(d-7),

with c_3(5) = c_4(5) = 2, c_5(5) = 0, c_5(6) = 4, and c_3(6) = 3, c_4(6) = 2. The order-5 cores of type 3 are RT_5
(any vertex) and C3 behind one chain vertex.

## (C) Main theorem

Let K(6) = 3 and K(n) = 5 for n >= 7.

**(C1)** No 2-strong tournament of order n >= 6 has hc <= K(n).

**(C2)** Every strong tournament of order n >= 6 with hc <= K(n) and no sole exit:
- has lambda = 1 and a unique core with e_u(D) = hc >= 3;
- has a uniquely Hamiltonian rooted prefix (or a one-vertex prefix);
- has a sole entry.

**(C3)** u_2(n) = 0 for n >= 6.

**(C4)** With R(m) = F(2m - 5) rooted (uniquely Hamiltonian, sole exit) pairs (THM-4868 (F)) and R(2) = 1:

    u_k(n) = sum_d c_k(d) R(n - d + 1),   U_k(x) = C_k(x) x (1 - 2x)/(1 - 3x + x^2).

So u_3(n) = F(2n - 9) for n >= 6, u_4(n) = 2F(2n - 11) for n >= 7, and u_5(n) = 4F(2n - 12) - 3F(2n - 14) for n >= 7.

**Proof of (C1)-(C2).** By induction:
- **(C2) at n from (C1) at n,** together with (A) and B1-B5. Since e_u >= 3 and hc <= 5, s_l = 1.
- **(C1) at n >= 13 from (C2) at n - 1.** Every T - v (strong, hc(T - v) <= hc(T)) has a sole exit z (so z -> v, and z
  has out-degree 2 in T) or a sole entry y (so v -> y, and y has in-degree 2).
  - Hence V is contained in N+(O_2) | N-(I_2).
  - |O_2| <= 4 forces |N+(O_2)| <= |O_2| + 2|O_2| - C(|O_2|, 2) <= 6, and likewise for I_2.
  - So n <= 12.
- **Base.** Orders 6-12, by the complete census: the only 2-strong classes with hc <= 5 have order 5 (RT_5, hc 2) or
  order 6 (hc 4 and 5).

QED.

## Verification

All FINITE-EXACT:
- The census (`fewcycle_sole_exit_free_20261104.c`, complete by Moon 1966 plus hc(T - v) <= hc(T)):
  - reproduces klein's exhaustive t_k, u_k (k <= 5) at orders 5-10;
  - confirms HYP-9461 P1-P10 at orders 10-14.
- Lineage, blocks, and the uniquely Hamiltonian end-count lemma: orders <= 13 and <= 11.
- Minimum hc of 2-strong tournaments: orders <= 11.
