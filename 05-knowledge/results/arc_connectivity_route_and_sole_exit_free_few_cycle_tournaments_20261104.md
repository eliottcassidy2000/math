# The arc-connectivity route to the one-arc join theorem, and sole-exit-free tournaments with few Hamiltonian cycles: terminal cores, doubling, and klein's Fibonacci counts proved

**Session:** opus-2026-11-04-S75. Worktree `C:/tmp/wt-s75`. Canon: THM-4869. Pre-registration: HYP-9461.

## The owner's prompt

The prompt, verbatim:

> think about how the arc-connectivity shortcut is not necessarily ruled out, since four order-9 tournaments lose
> Hamiltonicity when one arc is deleted but A decomposable tournament which has at least 9 vertices is {−3}-self dual if
> and only if it is either transitive or almost transitive. investigate related to a cycle-count table shows three-cycle
> tournaments without an out-degree-1 vertex at orders 6–10 counting 2, 5, 13, 34, 89, the odd-indexed Fibonacci
> numbers. (π^(n/2))/(n/2)! represents the volume of an n-dimensional unit ball and has a maximum integer value around
> 5, like the ratio of 8 to 4 of isomorphism classes in tournaments at n=5. consider past work regarding topics around
> -7/4 and -29/16 and the numerous adjacent topics they have garnered over their many mentions. consider how previous
> microcosm-macrocosm work relates to MRB constant adjacent ideas, the same values realized twice on either side of the
> peak, one instance closer to 0, the other closer to infinity; small changes in the first are big changes in the latter

## Status

| item | typing |
|---|---|
| A. A 2-arc-strong tournament with a cut vertex has hc >= 9 (sharp from order 7). With Fraisse-Thomassen for 2-strong tournaments, this gives a third proof that a uniquely Hamiltonian tournament has arc-strong connectivity 1 (THM-4867's headline). | PROVED; the corollary uses CITED Fraisse-Thomassen |
| A. Why route (b) failed: forced arcs need a cut vertex. The four order-9 witnesses are v => C3 => x => y => C3 => v, with forced arc x -> y. | PROVED structure; FINITE-EXACT minimality of order 9 (the S74 audit's exhaustive order-9 run) |
| A. The {-3}-self-duality parallel | ANALOGY, with the precise statement |
| B. Unique terminal core, product law, the core lemma, and the block list for k <= 5 | PROVED (uses THM-4868 (C), Redei, and finite checks at orders <= 6) |
| C. No 2-strong tournament of order n >= 6 has <= 3 Hamiltonian cycles, and none of order n >= 7 has <= 5 | PROVED (induction, with base orders 6-12 FINITE-EXACT) |
| C. u_2(n) = 0 (n >= 6), u_3(n) = F(2n - 9) (n >= 6), u_4(n) = 2F(2n - 11) (n >= 7), u_5(n) = 4F(2n - 12) - 3F(2n - 14) (n >= 7). This is klein's HYP-9459 (1)-(3), plus a new k = 5 law. | PROVED; pre-registered and confirmed at orders 10-14 (HYP-9461) |
| C. Minimum hc of 2-strong tournaments: 2, 4, 8, 13, 19, 30 at n = 5-10; > 40 at n = 11 | FINITE-EXACT |
| Inspiration threads (section 4) | checked: ANALOGY or NUMEROLOGY, except where marked |

New to the repository. No prior source was found for the sole-exit-free counts, for the bound "minimum out-degree >= 2
and n >= 6 implies at least 3 Hamiltonian cycles", or for the core decomposition. The search covered Thomassen 1980,
Busch 2006, the Kuehn-Osthus survey, Demaria-Kiihl and Burzio-Demaria. Their texts were not read, only abstracts and
citations.

## 0. Notation

- T is a tournament; hc(T) is its number of Hamiltonian cycles, and H(T) its number of Hamiltonian paths (odd, by Redei).
- end_x(T) and start_x(T) count Hamiltonian paths ending or starting at x.
- A **sole exit** has out-degree 1; a **sole entry** has in-degree 1.
- T is **2-strong** if T - v is strong for every v; otherwise a strong T has a **cut vertex**.
- lambda(T) is the arc-strong connectivity.
- A **one-arc-leaving set** D is a non-empty proper vertex set with exactly one arc u -> l leaving it. Then lambda(T) = 1
  iff such a D exists, i.e. iff T is a one-arc join (THM-4867).
- t_k(n) and u_k(n) are klein's counts (HYP-9459): classes of order n with hc = k, and those among them with no sole exit.
- F is Fibonacci and L is Lucas.

## 1. The arc-connectivity shortcut is not ruled out (Theorem A)

**What failed in S74.** S74's route (b) (`uniquely_hamiltonian_one_arc_cut_shelor_20261104.md` section 3) was the
statement "if T is 2-arc-strong then T - e is Hamiltonian for every arc e". It fails at order 9, in exactly four classes
(the S74 audit's exhaustive run).

**What the shortcut actually needs** is much weaker: 2-arc-strong implies hc >= 2. Split by vertex connectivity.

**Theorem A (cut vertex).** Let T be 2-arc-strong with a cut vertex v, and write T - v = D_1 => D_2 => ... => D_m (strong
components, all arcs forward, m >= 2). Then

    hc(T) = S_+ * prod_{1 < i < m} H(D_i) * S_-,   S_+ = sum over a in N+(v) & D_1 of start_a(D_1),
                                                  S_- = sum over b in N-(v) & D_m of end_b(D_m),

and S_+, S_- >= 3. Hence hc(T) >= 9. The bound is sharp: v => C3 => C3 => v (order 7) has hc = 9, and 9 is the
minimum at every order 7-11 (FINITE-EXACT, `fewcycle_sole_exit_free_20261104_mincycles_11_40.out`).

*Proof.*

**The formula.** Every Hamiltonian cycle passes v, leaves it into D_1, and must then cover D_1, ..., D_m in order, each
as one block. All arcs between blocks go forward, so the junctions are free.

**Size of D_1.** Only v sends arcs into D_1. Since T is 2-arc-strong, v -> at least 2 vertices of D_1. So |D_1| >= 2,
hence |D_1| >= 3 (it is strong).

**S_+ >= 3.** Every vertex of D_1 that v does not beat needs in-degree >= 2 inside D_1.
- If hc(D_1) >= 2, every start_a(D_1) >= 2, and |N+(v) & D_1| >= 2 gives S_+ >= 4.
- If D_1 = C3, every vertex has in-degree 1 in D_1, so v beats all three and S_+ = 3.
- If D_1 is the strong 4-tournament W4 (uniquely Hamiltonian), v must beat its two sole entries a, b, and
  start_a + start_b = 1 + 2 = 3.
- If D_1 is uniquely Hamiltonian of order m >= 5, take the non-critical sole entry y at the start of its non-critical
  run (THM-4868 (C)).
  - y has in-degree 1 inside D_1, so v -> y.
  - D_1 - y is strong, so every vertex starts at least one of its Hamiltonian paths. y beats all of D_1 except y-.
  - Hence start_y(D_1) >= m - 2 >= 3, and a second out-neighbour of v adds at least 1.

**S_- >= 3.** Dually. QED.

**Corollary A' (third proof of THM-4867's headline).** A uniquely Hamiltonian tournament is not 2-arc-strong:
- with a cut vertex, hc >= 4 (even S_+, S_- >= 2 suffice, which is elementary);
- if 2-strong, Fraisse-Thomassen (Graphs Combin. 3 (1987) 239-250, k = 2; CITED, not read) gives a Hamiltonian cycle
  avoiding any one arc, so hc >= 2.

Hence lambda = 1, a one-arc join. This proof uses neither Shelor's grammar nor THM-4868. It gives only the one-arc-cut
form, not the degree-1 form (Shelor 1977; THM-4868 (B)).

**Why the strong form fails exactly as observed.** An arc is *forced* if it lies in every Hamiltonian cycle (then
T - e is non-Hamiltonian).
- By Fraisse-Thomassen, a 2-arc-strong T with a forced arc has a cut vertex.
- In the product formula, a forced arc must lie in every choice of every factor. Arcs v -> a and b -> v are never
  forced (|N+(v) & D_1| >= 2). An arc inside a strong middle component is never forced (the rotations of one of its
  Hamiltonian cycles miss each cycle arc once).
- So a forced arc is either:
  - an arc of D_1 lying on every Hamiltonian path of D_1 that starts in N+(v);
  - the dual inside D_m;
  - a junction end(P_i) -> start(P_(i+1)) whose endpoints are forced. The simplest case is two consecutive trivial
    middle components {x}, {y}.
- The smallest instance of the junction case is 1 + 3 + 1 + 1 + 3 = 9 vertices: v => C3 => x => y => C3 => v.
  - v must beat all of the first C3 and be beaten by all of the last; otherwise a C3 vertex has in-degree 1.
  - The arcs between v and x, y are free: 4 classes.
- This matches the S74 audit's exhaustive order-9 count (4 classes, one bad arc each, |Aut| = 9 = 3 * 3, bad arc x -> y),
  and its finding that no witness exists at order <= 8.
- Each witness has hc = 3 * 1 * 1 * 3 = 9.
- Two witnesses (v -> x, y -> v, and the reverse pattern) are self-converse; the other two (v beats both, both beat v)
  are converses of each other.

**The {-3}-self-duality parallel (ANALOGY).** The owner's quote is Theorem 4 of Achour, Boudabbous and Boussairi, "The
{-3}-reconstruction and the {-3}-self duality of tournaments" (arXiv:1204.2513, Ars Combinatoria; abstract and theorem
statement read):
- T is {-k}-self dual if T - X is self-dual for every k-set X;
- a decomposable tournament with at least 9 vertices is {-3}-self dual iff it is transitive or almost transitive;
- an almost transitive tournament is a transitive one with the arc between its two extreme vertices reversed. It is
  the uniquely Hamiltonian A_n, the run-length-2 class of THM-4868.

The parallel is real but only structural:
- The four route-(b) witnesses are decomposable at order 9: both C3's are modules.
- Their 5-vertex quotient on {C3, x, y, C3', v} is the almost transitive A_5 exactly when v beats both x and y, or both
  beat v.
- They are not {-3}-self dual. By the theorem, a decomposable tournament of order 9 would have to be transitive or
  almost transitive, and these are neither.

In both statements, order 9 is where a decomposable structure is first forced into an (almost) transitive skeleton.

## 2. Sole-exit-free tournaments with few Hamiltonian cycles (Theorem B)

### 2.1 Terminal cores

**Lemma B1.** In a strong tournament, distinct minimal one-arc-leaving sets are disjoint. Two disjoint one-arc-leaving
sets D_1, D_2 satisfy |D_1| |D_2| <= 2. So if T has no sole exit, its minimal one-arc-leaving set D is unique, strong,
and of order >= 3.

*Proof.*

**Intersecting minimal sets.** Suppose D_1 and D_2 intersect, with union not V. Submodularity
d+(X & Y) + d+(X | Y) <= d+(X) + d+(Y) = 2 gives d+(D_1 & D_2) = 1, contradicting minimality.

If D_1 | D_2 = V, put A = D_1 - D_2, B = D_2 - D_1 and C = D_1 & D_2.
- D_1 sends exactly one arc into B, and D_2 exactly one into A. So C sends at most 2 arcs out.
- Minimality excludes exactly 1, so C sends one arc into A and one into B.
- Then no arc joins A and B at all, which is impossible.

**Disjoint sets.** Every arc between disjoint D_1 and D_2 leaves one of them.

**Size.** In a sole-exit-free T, a minimal set has order >= 3:
- order 1 would be a sole exit;
- order 2 forces a sole exit or a sink.

**Strong.** The last strong component of T[D] contains u and leaves D only through u -> l. QED.

**Lemma B2 (product law).** For a one-arc-leaving D with arc u -> l,

    hc(T) = s_l(T - D) * e_u(T[D]),

where s_l counts Hamiltonian paths of T - D starting at l (1 if T - D = {l}), and e_u = end_u. Every Hamiltonian cycle
traverses D as one block ending at u, and T - D as one block starting at l. If T - D has at least 2 vertices,
s_l(T - D) = hc(Q) for the contraction Q = T/D. In Q, the contracted vertex z beats only l, so it is a sole exit.

**Lemma B3 (the core).** Let T have no sole exit, lambda(T) = 1, and terminal core D with arc u -> l. Then:
- no vertex of D - u has out-degree 1 in D;
- no non-empty Z inside D - u sends exactly one arc into D - Z;
- e_u(D) >= 2.

Conversely, take any strong Q with a sole exit z (or |Q| = 2) and any pair (D, u) with these properties. Replace z by
D, keeping only u -> l from D. The result T = Q[z <- (D, u)] is strong, has no sole exit, and has D as its terminal core.
Isomorphism classes correspond.

*Proof of e_u(D) >= 2.*
- If e_u(D) = 1, then hc(D) <= e_u(D) = 1 (each cycle, cut after u, is a path ending at u). So D is uniquely
  Hamiltonian.
- Its sole exits lie in {u}, so |D| >= 4, since C3 has three sole exits.
- THM-4868 (C) gives a non-critical sole exit x (the end of the non-critical run). It must be u.
- Then D - u is strong and u is beaten by all of D except u+, so e_u(D) >= |D| - 2 >= 2. Contradiction.

The rest is direct. In T, a vertex of D - u has the same out-degree as in D, and every vertex of Q - z gains at least
|D| - 1 >= 2 out-neighbours. QED.

### 2.2 The core lemma and the blocks

**Lemma B4 (core lemma).** Let D - u = W_1 => ... => W_r be the strong components, and A = N-(u) & W_r. Then:

    e_u(D) = prod_{i < r} H(W_i) * sum over x in A of end_x(W_r),

with |A| >= 2 by minimality (Z = W_r). If e_u(D) <= 5:
- every W_i with i < r is a single vertex, since H(W_i) is 1 or at least 3, and the sum is at least 2;
- the core is a transitive chain w_1 => ... => w_p followed by a **block** (W, A), with W = W_r.

A block (W, A) of type k satisfies:
- W is strong, |A| >= 2, and the sum over A of end_x(W) is k;
- every vertex of W - A has out-degree >= 2 in W;
- no non-empty Z inside W has (arcs Z -> W - Z) + |Z & A| = 1.

Conversely, every such configuration is a core, given that:
- u is beaten by A and beats W - A;
- u -> w_1;
- the arcs between u and w_2, ..., w_p are arbitrary;
- p = 0 is allowed only when A != W (else u is a sink).

Distinct choices give distinct cores, because the chain is the canonical sequence of trivial components. So a block of
order w contributes 2^(d - w - 2) cores of each order d >= w + 2, plus one core of order w + 1 when A != W.

**Lemma B5 (blocks of type k <= 5).** Up to isomorphism, the blocks are:
- k = 2: none;
- k = 3:
  - (C3, all three vertices);
  - (W4, its two sole exits);
- k = 4: (W4, its two sole exits and one of the other two vertices), two blocks;
- k = 5:
  - (W4, all four vertices);
  - three blocks (W, A) with W uniquely Hamiltonian of order 5 and |A| = 2.

Here W4 is the strong 4-tournament a -> b -> c -> d -> a with a -> c and b -> d. Its end counts are
(a, b, c, d) = (1, 1, 2, 1), and its sole exits are c and d.

*Proof.*

**hc(W) >= 3.** Every end_x(W) >= hc(W), so the sum is at least 6.

**hc(W) = 2.** If W has order m <= 5, use the finite check below. Otherwise u_2(m) = 0 (Theorem C, whose k = 2 case uses
only this lemma's k <= 3 part) gives W a sole exit z, and z must lie in A.
- *If z is non-critical,* e_z >= m - 2 >= 4, so the sum is at least 6.
- *If z is critical,* write W - z = Q_1 => ... => Q_s with s >= 2. The only out-neighbour z' of z lies in Q_1 (only z
  enters Q_1). Then:
  - e_z(W) = H(W - z) = prod H(Q_i);
  - hc(W) = start_(z')(Q_1) * prod_{i >= 2} H(Q_i) = 2.
  - Every H is odd (Redei), so start_(z')(Q_1) = 2 and Q_2, ..., Q_s are single vertices.
  - The last of them beats only z, so it is a sole exit and lies in A.
  - So A = {z, q} and the sum is at least H(Q_1) + 2. H(Q_1) is odd and at least start_(z')(Q_1) = 2, so it is at
    least 3; a sum <= 5 forces H(Q_1) = 3, i.e. Q_1 = C3. But then start_(z')(C3) = 1, not 2, a contradiction.

**hc(W) = 1 (uniquely Hamiltonian, order m).** A contains every sole exit.
- **m = 3** (C3): A = W, sum 3.
- **m = 4** (W4): A contains {c, d}, and the sums are as listed.
- **m >= 5:** the non-critical sole exit x of THM-4868 (C) has e_x >= m - 2. So the sum is at least m - 1, and m <= 6
  for k <= 5.

**Order of the proof.** For k <= 3 the case hc(W) = 2 is excluded by the first line alone, since the sum is at least
4. Theorem C (C1)-(C3) uses only the k <= 3 list and the absence of type-2 blocks. The k = 4, 5 lists, which use
u_2(m) = 0 from (C3), are needed only for the counts in (C4). There is no circularity.

**Finite check.** Uniquely Hamiltonian blocks of orders 5 and 6, and blocks with hc(W) = 2 of order <= 5: the exhaustive
check in `fewcycle_blocks_20261104.out` finds exactly the three order-5 blocks for k = 5, and nothing of order >= 6 for
any k <= 5. The same program finds no violation of e_z >= m - 2 (non-critical sole exits) or of its dual for sole
entries, in any uniquely Hamiltonian tournament of order 5-11. QED.

**Core counts.**

    c_3(5) = 2,  c_3(d) = 2^(d-5) + 2^(d-6) = 3 * 2^(d-6)  (d >= 6)
    c_4(5) = 2,  c_4(d) = 2 * 2^(d-6) = 2^(d-5)            (d >= 6)
    c_5(5) = 0,  c_5(6) = 1 + 3 = 4,  c_5(d) = 2^(d-6) + 3 * 2^(d-7) = 5 * 2^(d-7)  (d >= 7)

The two order-5 cores of type 3:
- the C3 block with one chain vertex (scores 1, 2, 2, 2, 3);
- the W4 block with p = 0, which is the regular tournament RT_5 with any vertex as u.

So **RT_5 seeds the three-cycle family.**

## 3. The main theorem (Theorem C)

**Theorem C.**
- **(C1)** No 2-strong tournament of order n >= 6 has at most 3 Hamiltonian cycles. None of order n >= 7 has at most 5.
- **(C2)** Let T be strong of order n >= 6 with hc(T) <= 3 (n >= 7: hc(T) <= 5) and no sole exit. Then T is a one-arc
  join Q[z <- (D, u)]:
  - Q is uniquely Hamiltonian with sole exit z (or |Q| = 2);
  - (D, u) is a core with e_u(D) = hc(T) >= 3;
  - T has a sole entry.
- **(C3)** u_2(n) = 0 for n >= 6. This is HYP-9459 (1), and the open half of THM-4868 (B): every strong tournament with
  minimum out-degree >= 2 and n >= 6 has at least 3 Hamiltonian cycles.
- **(C4)** With R(2) = 1 and R(m) = F(2m - 5) rooted pairs (uniquely Hamiltonian, sole exit) of order m >= 3
  (THM-4868 (F)):

      u_k(n) = sum_d c_k(d) R(n - d + 1),

  and the generating function is U_k(x) = C_k(x) * x(1 - 2x)/(1 - 3x + x^2). The factor (1 - 2x) of the rooted prefixes
  cancels the doubling pole of the cores:

      U_3(x) = x^6 (2 - x)/(1 - 3x + x^2)       u_3(n) = F(2n - 9)                n >= 6
      U_4(x) = 2x^6 (1 - x)/(1 - 3x + x^2)      u_4(n) = 2F(2n - 11)              n >= 7
      U_5(x) = x^7 (4 - 3x)/(1 - 3x + x^2)      u_5(n) = 4F(2n - 12) - 3F(2n - 14) = L(2n - 12) + F(2n - 13)   n >= 7

  At n = 6, u_4 and u_5 also count the 2-strong classes (one with hc 4, four with hc 5).

*Proof.* Let K(6) = 3 and K(n) = 5 for n >= 7. Two statements are proved by induction:
- **P(n):** no 2-strong tournament of order n has hc <= K(n);
- **S(n):** every strong T of order n with hc <= K(n) has a sole exit or a sole entry.

**S(n) from P(n).** Let T have no sole exit.
- By P(n) and Theorem A (hc >= 9), T is not 2-arc-strong, so lambda(T) = 1.
- Lemmas B1-B3 give the core, with e_u(D) >= 2. e_u(D) = 2 is excluded (no blocks of type 2).
- Then s_l = hc/e_u < 2, so s_l = 1 and Q is uniquely Hamiltonian, or |Q| = 2.
- If |Q| = 2, l is a sole entry.
- Otherwise a sole entry of Q other than z (z has in-degree |Q| - 2) keeps its in-degree in T.

**P(n) for n >= 13 from S(n - 1).** Suppose T is 2-strong of order n with hc(T) <= 5.
- Every T - v is strong, with hc(T - v) <= hc(T) (insertion at a switch is injective).
- By S(n - 1), T - v has a sole exit z (then z has out-degree 2 in T and z -> v) or a sole entry y (in-degree 2 in T,
  v -> y).
- Let O_2 and I_2 be the vertices of out-degree 2 and of in-degree 2. Then V is contained in N+(O_2) | N-(I_2).
- Bounding O_2: C(|O_2|, 2) <= 2|O_2| gives |O_2| <= 5, and |O_2| = 5 would make O_2 a sink set, so |O_2| <= 4.
- Hence |N+(O_2)| <= |O_2| + (2|O_2| - C(|O_2|, 2)) <= 6, and likewise |N-(I_2)| <= 6. So n <= 12, a contradiction.

**Base.** P(6), ..., P(12) hold by the complete census `fewcycle_sole_exit_free_20261104_14_5.out`: the 2-strong classes
with hc <= 5 occur only at orders 5 (RT_5) and 6 (hc 4 and 5).

**(C3) and (C4).** These follow from (C2), Lemma B5 and the generating-function algebra:
- R(x) = x^2 (1 - 2x)/(1 - 3x + x^2);
- C_3(x) = x^5 (2 - x)/(1 - 2x);
- C_4(x) = 2x^5 (1 - x)/(1 - 2x);
- C_5(x) = x^6 (4 - 3x)/(1 - 2x).

QED.

**The minimum number of Hamiltonian cycles of a 2-strong tournament** is 2, 4, 8, 13, 19, 30 at n = 5-10, and above 40
at n = 11 (FINITE-EXACT; `..._mincycles_11_40.out`). The general minimum is open (Kuehn-Osthus survey). Exponential lower
bounds are attributed to Thomassen 1980 and Busch 2006 (not read). Theorem C gives the small cases.

### 3.1 Computational record (HYP-9461)

**The census.** `fewcycle_sole_exit_free_20261104.c` generates every strong tournament with hc <= K, by vertex
extension. It is complete because a non-critical deletion never increases hc, and Moon 1966 provides a non-critical
vertex.
- It reproduces klein's exhaustive t_k, u_k (k <= 5) at orders 5-10 exactly.
- Pre-registered in c355721ce6 before the run: P1-P10 at orders 10-14, **all confirmed**. Values:
  - u_3 = 89, 233, 610, 1,597, 4,181;
  - u_4 = 68, 178, 466, 1,220, 3,194;
  - u_5 = 60, 157, 411, 1,076, 2,817;
  - u_1 = u_2 = 0;
  - t_2 = 589, 1,775, 5,257, 15,360, 44,394.
- Every u-class has exactly one minimal one-arc-leaving set, with s_l = 1, and 0 product-law failures.
- The whole core table equals c_k(d) R(n - d + 1). The new cores at d = 9-13 are:
  - c_3 = 24 ... 384;
  - c_4 = 16 ... 256;
  - c_5 = 20 ... 320.
- There are no 2-strong or cut-vertex classes with hc <= 5. The run takes 27 s.

**Lineage** (`fewcycle_core_lineage_20261104.out`):
- Exactly half the cores of each order are *attached*: deleting the chain vertex w_1, a sole entry, leaves a core, and
  every core has exactly one attached child. This happens iff u -> w_2.
- The other half are *orphans*. Their u-out-degrees are the coefficients of z(1 + z)^(d-7)(1 + z + z^2) for k = 3 and
  2z^2 (1 + z)^(d-6) for k = 4: the binomial counts of the free arcs u -- w_i in Lemma B4.

## 4. The inspiration threads, checked

**The unit ball and n = 5.**
- The repository's own facts:
  - V_n = pi^(n/2)/(n/2)! peaks at n = 5, since V_(n+2) = 2 pi V_n/(n+2) and V_6/V_5 = 5 pi/16 < 1 (session127
    `gsma_norm_collapse_20261031/volume/THEORY.md`, VOL-1).
  - The self-converse share of tournament classes (8 of 12 at n = 5; 1, 1/2, 2/3, 12/56, 88/456, ...) has its unique
    interior local maximum at n = 5 (S74 audit m7; `two_selector_boolean_20261102/INSPIRATION.md`).
- What this session adds is a *proved* reason why n = 5 is special for few-cycle tournaments:
  - RT_5 is the only tournament with minimum out-degree >= 2 and two Hamiltonian cycles (Theorem C3);
  - RT_5 with a marked vertex is the order-5 core that seeds the F(2n - 9) family;
  - the 2-strong minimum hc sequence starts 2, 4, 8, ... at n = 5.
- The "peak at 5" shared with V_n is ANALOGY. The RT_5 statements are PROVED. The 8:4 ratio plays no role in the
  mechanism (NUMEROLOGY).

**-7/4 and -29/16.**
- The thread: x -> x^2 - 29/16 has the rational 3-cycle -1/4 -> -7/4 -> 5/4 and eight rational preperiodic points.
  - THM-4139 and THM-4146: the order-six lift and fibre firewall.
  - -7/4 is the period-3 airplane root (THM-4714).
  - Further: the Zsigmondy triad (c in {0, -1, -2, -7/4}), the time reversal L(x) = -x - 1/2 fixing -29/16, and the mod-290
    PRIMEGAME descent F(a) = (a^2 - 29)/4.
- Touchpoints with this session:
  - The directed 3-cycle C3 is one of the two blocks of the three-cycle cores.
  - "Period 3 implies all periods" (Sharkovskii) has the classical tournament counterpart "a strong tournament is
    pancyclic" (Moon; Harary-Moser).
  - Both are ANALOGY.
- The Lucas readings -29/16 = -L(7)/2^4 and -7/4 = -L(4)/2^2 are NUMEROLOGY. The counts here involve L(2m) and F(2m + 1),
  and no mechanism links them to the quadratic map.

**Microcosm-macrocosm and the MRB mirror.**
- The owner's image of one value realized twice across a peak (one near 0, one near infinity, small changes on one side
  large on the other) already has two exact forms in the repository:
  - **Session31's mirror map for x^(1/x)** (`session31_mirror_mrb_20261022/mirror/README.md`): a decreasing involution
    fixing e, with mu e^(-mu) = lambda e^(-lambda), and 2 <-> 4 the only rational pair.
  - **THM-062's palindromic forward-edge polynomial,** with root pairs r, 1/r.
- Here the same shape is exact and PROVED. The denominator 1 - 3x + x^2 is self-reciprocal, with roots phi^(-2) (near 0)
  and phi^2 (near infinity) exchanged by x -> 1/x; and dx' = -dx/x^2 is the "small change near 0, large near infinity"
  of the image.
- The finitely many blocks (C3, W4, three 5-blocks), the microcosm, produce cores with a pole at x = 1/2 (doubling). The
  uniquely Hamiltonian prefixes carry the factor (1 - 2x), which cancels that pole exactly. Only the golden reciprocal
  pair survives in the macrocosm counts F(2n - 9), 2F(2n - 11), ....
- Calling the cancellation microcosm -> macrocosm is ANALOGY. The cancellation itself is the algebra of (C4).

## 5. Typing

| claim | type |
|---|---|
| Theorem A: product formula, S_+, S_- >= 3, hc >= 9 | PROVED (uses THM-4868 (C) for the sharp constant; S_+, S_- >= 2 is elementary) |
| Corollary A': uniquely Hamiltonian implies lambda = 1 | PROVED given Fraisse-Thomassen (CITED) |
| Forced arcs need a cut vertex; the order-9 family | PROVED structure; order 9 minimal FINITE-EXACT (S74 audit) |
| {-3}-self-duality parallel | ANALOGY (theorem statement read; paper body not read) |
| Lemmas B1-B5 | PROVED (B3 and B5 use THM-4868 (C); B5 uses Redei and finite checks at orders <= 6) |
| Theorem C (C1)-(C4) | PROVED; base orders 6-12 FINITE-EXACT; P1-P10 pre-registered and confirmed at 10-14 |
| Minimum hc of 2-strong tournaments, n <= 11 | FINITE-EXACT |
| Lineage halves and Pascal rows | FINITE-EXACT through order 12, explained by Lemma B4 |
| Section 4 | as marked |

## 6. Pre-registration record

HYP-9461 was committed (c355721ce6) with the generator and the order 3-9 evidence before orders 10-14 were run. All of
P1-P10 hold. The proofs of section 2 and of (C1) were found after the run. They explain the doubling that the
pre-registration predicted from the order 5-8 cores.

## 7. Audit

(Filled after the independent audit.)

## 8. Reproduction

    gcc -O3 -march=native -fopenmp -o fewcycle fewcycle_sole_exit_free_20261104.c
    ./fewcycle 14 5            # census, about 30 s, = fewcycle_sole_exit_free_20261104_14_5.out
    gcc -O3 -march=native -fopenmp -DMAXK=60 -o fewcycle60 fewcycle_sole_exit_free_20261104.c
    ./fewcycle60 11 40         # minimum hc of 2-strong / cut-vertex tournaments
    gcc -O3 -march=native -fopenmp -o lineage fewcycle_core_lineage_20261104.c && ./lineage 13
    gcc -O3 -march=native -fopenmp -o blocks fewcycle_blocks_20261104.c && ./blocks 11
