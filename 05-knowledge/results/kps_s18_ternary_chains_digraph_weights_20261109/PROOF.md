# Proofs for THM-5070 and THM-5071 (kind-pasteur-2026-11-09-S18)

**Notation.**
- For a digraph X, I(Omega(X), x) = sum over sets of pairwise vertex-disjoint odd directed cycles of x^(number of
  cycles). This is the odd-cycle packing count; for tournaments, H_x(T) = I(Omega(T), x) (THM-4940).
- H(X) is the number of Hamiltonian paths, H(empty) = 1. OCF (Grinberg-Stanley) gives H(T) = I(Omega(T), 2) for
  tournaments.
- Paths are vertex sequences, and one-vertex paths are allowed. |P| is the number of vertices of P.

---

## THM-5070: port chains

### Setting

- **Ported tournament.** A ported tournament (T, u) is a tournament T with a marked vertex u, the port.
- **State.** Its state is s(T, u) = (h, E, O), where:
  - h = I(Omega(T), x);
  - E = sum over paths P of T starting at u with |P| even of I(Omega(T - V(P)), x);
  - O = the same sum over paths with |P| odd. The path (u) is included in O.
- **Blocks.** A block is a tournament B with two marked vertices w, u' (they may coincide).
- **Gluing.** The glued tournament T (.) B has vertex set V(T) disjoint union V(B). It keeps the arcs inside T and
  inside B. Every t in T beats every b in B, except that the arc between u and w is reversed: w -> u. The new port
  is u'.
- **Port chain.** A port chain is the result of gluing blocks one after another, starting from the one-vertex
  tournament ({u}), whose state is (1, 0, 1).
- **Prior art.** Up to taking converses, this is the seam operation X (.) Y of SEAM.md
  (seven_return_native_recursions_20261029), whose two-block count is Session29 ADDENDUM section 2. At x = 2,
  SEAM.md's (H, S, R) transfer is the case treated here, in other coordinates, by OCF.

### Lemma A (the cycles of T (.) B)

Every directed cycle of T (.) B lies in T, or lies in B, or uses the arc w -> u exactly once.

*Proof.* w -> u is the only arc from V(B) to V(T). A cycle meeting both parts must return from B to T, so it uses
w -> u; a cycle passes through u only once. ∎

### Lemma B (the new count)

Write:
- eta = I(Omega(B), x);
- q_o = sum over paths Q of B ending at w with |Q| odd of I(Omega(B - V(Q)), x);
- q_e = the same sum over paths with |Q| even.

Then I(Omega(T (.) B), x) = eta h + x (q_o E + q_e O).

*Proof.*
- **Shape of the new cycles.** A cycle through w -> u is u = p_1, ..., p_s, q_1, ..., q_r = w, back to u. Here
  P = (p_1, ..., p_s) is a path of T from u, and Q = (q_1, ..., q_r) is a path of B ending at w. The arc p_s -> q_1
  is a cross arc, present because q_1 != w or p_s != u.
  - If s = 1 and r = 1, the "cycle" would be u -> w -> u. That needs u -> w, which was reversed, so it is excluded
    automatically: |P| + |Q| = 2 is even and never counted.
  - Conversely, every such pair (P, Q) with |P| + |Q| odd is an odd cycle.
- **Packings with a new cycle.** All new cycles contain u, so a packing contains at most one of them. With one new
  cycle C, the rest of the packing is a packing of (T - V(P)) together with (B - V(Q)). These two are joined by
  forward arcs only, because w and u are deleted. So no odd cycle crosses, and the packing count factors.
- **Packings without a new cycle.** These are packings of T together with B joined forward, which give eta h.
- **Parity bookkeeping.** |P| + |Q| is odd exactly when (|P| even, |Q| odd) or (|P| odd, |Q| even). This gives the
  terms x E q_o and x O q_e. ∎

### Lemma C (the new port)

E' and O' are integer combinations of h, E, O with nonnegative coefficients that depend only on (B, w, u', x).

*Proof.* Let P' be a path of T (.) B from u'. Because w -> u is the only arc back into T, there are three cases:
- **(i) P' lies in B and avoids w.** The remainder is T together with B - V(P'), still joined with the arc w -> u.
  By Lemma B, applied to (T, u) and the block B - V(P') with the same w, the remainder has packing count
  h I(Omega(B - V(P')), x) + x (r_o E + r_e O). Here r_o and r_e count paths of B - V(P') ending at w, weighted as
  in Lemma B.
- **(ii) P' lies in B and contains w without continuing to u.** The remainder is joined forward only, so its
  packing count is h I(Omega(B - V(P')), x).
- **(iii) P' = Q_1 P_T Q_2.** Here:
  - Q_1 is a path of B from u' to w;
  - P_T is a path of T from u;
  - Q_2 is a possibly empty path of B - V(Q_1). After P_T the path may re-enter B at any unused vertex by a forward
    arc, but it can never return to T.

  The remainder (T - V(P_T)) together with (B - V(Q_1) - V(Q_2)) is joined forward only, so its packing count
  factors. Summing over P_T of even or odd size yields E or O. The parity of P' is that of |P_T| + |Q_1| + |Q_2|.

In every case the coefficient of h, E or O is a sum of packing counts of subtournaments of B, hence a nonnegative
integer when x is a nonnegative integer. ∎

**Theorem (1).** s(T (.) B, u') = N_x(B, w, u') s(T, u), with N_x a nonnegative integer 3 x 3 matrix. It is
computed by `chain_blocks.py`.

**Checks.**
- `chain_check.py` builds 600 random port chains of orders <= 14 from blocks of order <= 5 (both bases).
- For each it compares the predicted h with an independent count of the actual tournament: Hamiltonian-cycle DP
  per vertex subset, then a set-partition DP (`packing_count.c`).
- It also checks OCF, direct H = I(Omega, 2).
- Result: 0 failures (`out/chain_check.out`).

### (2) Apex bridges

Suppose every odd cycle of T passes through u, and every odd cycle of B passes through w.
- **Formula.** Then every remainder T - V(P) (P from u) and B - V(Q) (Q to w) has no odd cycle. So
  I(Omega(T (.) B), x) = I(Omega(T), x) I(Omega(B), x) + x c, where c is the number of odd cycles through w -> u.
- **One-apex blocks.** Take T = TT_k + u, where u beats the chain vertices i with b_i = 1, and
  W = sum_i b_i 2^(k-i). The paths from u are u alone, or u followed by an increasing chain path from a vertex with
  b_i = 1. Counting by parity gives E = ceil(W/2) and O = 1 + floor(W/2).
- **The other side.** Dually, for B = TT_l + w with W' = sum_j [j -> w] 2^(j-1), we get q_e = ceil(W'/2) and
  q_o = 1 + floor(W'/2).
- **The count a.** The number a(W) of odd cycles of the one-apex tournament satisfies a(2W+1) = a(W) and
  a(2W) = a(W) + ceil(W/2). This is the inherited one-apex recurrence of session25 and Session191 one_apex/PROOF.md
  (1), written for a = (H - 1)/2.
- **Result.**
  I = (1 + x a(W))(1 + x a(W')) + x [ceil(W/2)(1 + floor(W'/2)) + (1 + floor(W/2)) ceil(W'/2)].
- **Prior art.** This is the one-arc join law (THM-4866 (B)) specialised to apex-to-apex arcs.

### (3) The tower

- **The block.** For the one-vertex block (w = u' = the new vertex): eta = 1, q_o = 1, q_e = 0.
  - Case (ii) contributes h to O' (the path (w)).
  - Case (iii) with Q_1 = (w), Q_2 empty contributes E to O' and O to E'.
- **The matrix.** So M_T = [[1, x, 0], [0, 0, 1], [1, 1, 0]].
- **The tournament.** Iterating this block builds the tower 1-^n (THM-4865 (c), THM-4952 (C)): each vertex beats
  only its predecessor.
- **Characteristic polynomial.** det(tI - M_T) = (t - 1)(t^2 - 1) - x = t^3 - t^2 - t - (x - 1). This is THM-4952
  (C)'s recurrence f(n) = f(n-1) + f(n-2) + (x-1) f(n-3).
  - At x = 2 it is tribonacci.
  - At x = 3 it is (t - 2)(t^2 + t + 1), with dominant eigenvector v = (3, 1, 2).

### (4), (5): finite computations

**(4) Coverage.**
- `chain_cover.c` takes the breadth-first closure of (1, 0, 1) under the block matrices, keeping h <= LIM.
- **Saturation.** E and O are saturated at CAP = floor(LIM/x) + 1. This is exact for two reasons:
  - Every h-row is (eta, x q_o, x q_e), so a coordinate >= CAP with a positive coefficient forces h > LIM.
  - With nonnegative integer matrices, min(., CAP) commutes with every row.
- **Order <= 5 versus order <= 4.** Blocks of order <= 5 include those of order <= 4. So the order-5 family's miss
  set below LIM is contained in the order-4 family's. That is why the order-5 runs only need to reach 3000.

**(5) Tower directions.**
- `tower_direction.py` multiplies all products of length <= 3 of the 57 order-<= 4 base-3 blocks.
- It keeps the distinct matrices N with N v parallel to v.
- Exactly M_T^1, ..., M_T^12 occur (`out/tower_direction.out`).

**(6) Bases 4-8.**
- **m <= 80.** THM-4940 (6) gives the exact gap sets.
- **80 < m <= 30000.** `bases_gaps.py` runs the order-<= 5 chains in base b, adds the attained alpha-sums <= 80, and
  closes the set under (m, m') -> m + m' + b m m' (order joins multiply I(Omega, b), THM-4940 (1)). The alpha-sums
  left open are:
  - base 4: 149;
  - base 5: 96, 129, 143, 234, 335, 542, 614;
  - base 6: 122, 130, 142, 146, 219, 285;
  - base 7: 131, 146, 150, 200, 278, 404, 466;
  - base 8: 40 values, the largest 8789.
- **Witnesses.** Each open value has an explicit tournament of order 8-11 (`out/bases_4to8_witnesses.txt`, found by
  `hunt.c`). Every witness was rechecked by `verify_witnesses.py` from scratch: tournament check, Hamiltonian-cycle
  DP per odd vertex set, set-partition DP, and OCF against a direct Hamiltonian path count. There were 0 mismatches
  in 61.

---

## THM-5071: the weighted cycle formula for digraphs

### (1) The Grinberg-Stanley formula

**Source.** Grinberg and Stanley, *The Redei-Berge symmetric function of a directed graph*, arXiv:2307.05569.
- For a digraph X on V they define U_X = sum over listings w of V of F_{Des(w, X)}, where Des(w, X) is the set of i
  with (w_i, w_{i+1}) an arc of X.
- Their p-expansion for arbitrary digraphs is
  U_X = sum_sigma (-1)^phi(sigma) p_type(sigma), where:
  - sigma ranges over the permutations of V each of whose nontrivial cycles is a directed cycle of X or of the
    complement Xbar;
  - phi(sigma) is the sum of (length - 1) over the cycles of sigma that are cycles of X.
- **Sign convention, checked on small cases.**
  - X a 2-cycle: U_X = 2 e_2 = p_1^2 - p_2.
  - X empty on two vertices: U_X = 2 h_2 = p_1^2 + p_2.
  - X a directed triangle: U_X = 3 h_3 + 3 e_3 = p_1^3 + 2 p_3.

**Evaluation.** Evaluate at (1, 0, 0, ...):
- F_S(1, 0, ...) = [S empty] and p_k(1, 0, ...) = 1.
- So the left side counts listings with no consecutive X-arc, which are the Hamiltonian paths of Xbar.

**Statement.** Put X = Dbar for a loopless digraph D. Then

    H(D) = sum over sigma, with every nontrivial cycle a cycle of D or of Dbar, of
           prod over cycles [ +1 if a cycle of D ; (-1)^(len - 1) if a cycle of Dbar ].           (1)

**Checks.** `digraph_weights.c` checks (1) against a direct Hamiltonian path DP on:
- every loopless digraph of order <= 4;
- every semicomplete digraph and every oriented graph of order 5;
- 10^5 random digraphs of order 6 and 10^4 of order 7.

That is 232,263 checks with 0 failures (`out/digraph_weights.out`).

### (2) Weighted packing form

**Pair types.** Classify each unordered pair {i, j} of D as:
- a T-pair (exactly one arc);
- a B-pair (both arcs, a reciprocal or "11" pair);
- a Z-pair (no arc, a mirrored "00" pair).

**Which cycles occur.**
- A cycle of D can use no Z-pair, and a cycle of Dbar can use no B-pair.
- A cycle using only T-pairs lies in D exactly when its reverse lies in Dbar.

**Grouping.** Group the two permutation cycles (C from D, its reverse from Dbar) on the same vertex set and pair
order. They contribute 1 + (-1)^(len - 1): 2 for odd length, 0 for even length. Every other cycle of D uses a B-pair
and contributes +1. Every other cycle of Dbar uses a Z-pair and contributes (-1)^(len-1). So

    H(D) = sum over packings of vertex-disjoint "admissible" cycles of prod w(C),                 (2)
    w = 2  : odd cycle of D through T-pairs only,
    w = 0  : even cycle through T-pairs only (may be dropped),
    w = 1  : cycle of D (length >= 2) through at least one B-pair,
    w = (-1)^(len-1) : cycle of Dbar (length >= 2) through at least one Z-pair.

**Special cases.**
- For tournaments only w = 2 survives: OCF.
- For semicomplete digraphs every weight is in {0, 1, 2}.
- For oriented graphs the defect weights are signed. A Z-pair 2-cycle has weight -1.

### (3) Parity

Weights 0 and 2 are even and weights +-1 are odd, so H(D) is congruent, mod 2, to the number of packings (the empty
one included) of defect cycles: B-cycles of D and Z-cycles of Dbar. With no defect pair this is 1, which is
Redei's theorem.

### (4a) Every positive integer is the count of a semicomplete digraph

**Construction.** Take TT_n on 0..n-1 (i -> j for i < j), and add the back arcs from the sink n-1 to a set S of
earlier vertices.

**Count.** A Hamiltonian path either ends at n-1 (one path), or passes n-1 -> x with x in S.
- In the second case the vertices after n-1 form an increasing run starting at x. That run consists of x and any
  subset of {x+1, ..., n-2}.
- The vertices before n-1 are the complement in increasing order.
- So H = 1 + sum_{x in S} 2^(n-2-x).

**Conclusion.** Choosing S by the binary digits of h - 1 gives H = h, with popcount(h - 1) reciprocal pairs.
`one_defect_witness.py` checks h <= 299.

### (4b) Exactly one reciprocal pair excludes 1, 4 and 7

Let D be loopless semicomplete with exactly one B-pair {u, v}.

**Step 1 (factorisation over strong components).**
- The condensation of a semicomplete digraph is a transitive tournament.
- Every Hamiltonian path visits the strong components in that order, so H(D) = prod H(D_i).
- {u, v} lies in one component D_0 (u and v reach each other). The other components are tournaments, so their H
  is odd.

**Step 2 (H(D_0) >= 2).**
- If |D_0| = 2, then H(D_0) = 2.
- If |D_0| = k >= 3, one of the two tournament resolutions of D_0 is strong (Step 3), so H(D_0) >= f(k). Here f is
  the strong floor: f(3) = 3, f(4) = 5, f(5) = 9, and f is nondecreasing (Moon; Busch; THM-1370 (3)).

**Step 3 (a strong resolution exists, self-contained proof).** Suppose both resolutions T (u -> v kept) and
T' (v -> u kept) are non-strong.
- **The two cuts.** T has a cut (S, S-bar) with every arc from S to S-bar. D_0 is strong, so D_0 has an arc back
  from S-bar to S; it can only be v -> u. Hence u is in S and v is in S-bar. Likewise T' has a cut (R, R-bar) with
  u in R-bar and v in R.
- **The four cells.** Write X1 = S∩R, X2 = S∩R-bar (contains u), X3 = S-bar∩R (contains v), X4 = S-bar∩R-bar.
  - Every non-reciprocal pair between X2 and X3 would have to point X2 -> X3 (from the first cut) and X3 -> X2
    (from the second). So X2 = {u} and X3 = {v}.
  - X1 beats every other cell and every other cell beats X4.
- **Contradiction.** Since D_0 is strong, X1 and X4 are empty, so k = 2.

**Step 4 (conclusion).**
- H(D) = 1 is impossible, by Step 2.
- If H(D) is 4 or 7, the odd cofactor forces H(D_0) = 4 or H(D_0) = 7. Then f(k) <= 7 < 9 = f(5) forces k <= 4.
- The exhaustive census of orders <= 4 with one reciprocal pair (`digraph_weights.c`, part 2) contains neither value.
  ∎

### (4c) Census

The exhaustive one-defect census of orders 3..7 is in `out/digraph_weights.out`.
- **One reverse arc added (one B-pair).** The values below 198 are exactly 2..197 without 4 and 7.
- **One arc deleted (one Z-pair).** Every value 0..129 occurs. 7 first occurs at order 5 and 21 at order 6.
- **Explicit matrices** are in `out/one_defect_witness.out`.

### (5) The seven and twenty-one cores

- `core_weights.py` enumerates all semicomplete digraphs of order <= 5. It keeps the strong ones with H in {7, 21},
  reduces them to isomorphism classes, and evaluates (2) term by term.
- It finds 3 classes with H = 7 at order 4 and 16 classes with H = 21 at order 5. Both counts agree with Session234's
  catalogue.
- Every H = 7 core has no odd T-only cycle: H = 1 + 6 (two classes, six pairwise intersecting B-cycles) or
  H = 1 + 5 + 1 (one disjoint pair of B-cycles).
- See `out/core_weights.out` for the 21-cores.
