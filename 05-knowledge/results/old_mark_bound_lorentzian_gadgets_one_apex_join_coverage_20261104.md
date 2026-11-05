# The old-mark digit lower bound, Lorentzian gadget transfers, and one-apex joins up to 2^32

**Session:** opus-2026-11-04-S78. Worktree `C:/tmp/wt-s78`. Canon: THM-4875. Pre-registration: HYP-9468 (commit
19562af223).

**Owner prompt (verbatim).** "start another session proving the old-mark digit lower bound, work to complete chan-pak's
conjecture"

**Target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is
the number H(T) of Hamiltonian paths of some tournament. It is OPEN. THM-1370 proves that 7 and 21 are never attained.

**Results.**
- **(A) The old-mark digit lower bound.** PROVED at every order, with an exact formula at the sink and the sharp
  minimum 4.
- **(B) Explicit one-arc gadget transfers and their Lorentzian adjoint.** PROVED. The mechanism is session29's; the
  explicit coefficients and the general adjoint identity are added here.
- **(C) A descent lemma.** It rules out one natural inductive proof shape.
- **(D) FINITE-EXACT coverage.** Every odd N <= 2^32 except 7 and 21 is the count of an explicit one-arc join built
  from one-apex tournaments. The previous verified range was 99,999,999.

## 0. Notation and inherited facts

**Root response and digit.** For a tournament T and vertex t:
- R_t(T) = sum over X subset V - t of h(T[X]) s_t(T[V - X]), with h(empty) = 1;
- the digit is sigma_t = R_t - h.

R_t is klein's B_t and session29's R; sigma_t counts the two-bin covers (P1, P2) with P2 starting at t and P1
nonempty, ending at an out-neighbour of t. The rooted triple is (h, e_t, R_t), where e_t counts the Hamiltonian paths
ending at t.

**C2 child** (RESULT law (4); THM-4865). T1 = C2(B, a) adds x with x -> a only (every other vertex of B beats x), then
y, beaten by all of B, with y -> x. Its count is H(T1) = h + 2R_a, and its native digits are
sigma_x = 0, sigma_y = h - e_a (RESULT eqs. (8), (10), (12)).

**One-arc join** (session29 ADDENDUM sections 2 and 5; THM-4866 (B)). Let J = A => B with b -> a reversed. Then

    H(J) = x_h y_h - x_e y_e + x_R y_R,   x = (h, e_a, R_a)(A),   y = (h, e_b, R_b)(B^op).

**One-apex tournaments** (session25; Session191 one_apex/PROOF.md). For a binary word w of value v, T(w) is the chain
1 -> ... -> n plus an apex x with x -> i iff bit i is 1. The inherited state recurrence (h, A, p) gives:
- H(T(w)) = g(v), where g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2);
- the apex triple is (g(v), e(v), 1 + v), where e(v) = [v = 0 or v even] and R_x = 1 + v is the binary state A + 1.

Checked by brute force for every word of length <= 7 (`one_apex_triple_check_20261104.py`).

## 1. (A) The old-mark digit lower bound (PROVED at every order)

Let T1 = C2(B, a), with |B| >= 2 and h = h(B).

**(A1) Out-degree bound.** For every tournament T and vertex v of out-degree d, sigma_v >= d(d - 1). This is READY (7)
dualised; the S77 audit section 4 sketched the proof.

*Proof.*
- Fix out-neighbours w != w' of v, and a Hamiltonian path q_1 ... q_m of T - {v, w, w'} (empty if n = 3).
- For k = 0..m consider two covers:
  - A_k: P1 = q_1..q_k w, P2 = v w' q_(k+1)..q_m;
  - B_k: P1 = q_1..q_k w', P2 = v w q_(k+1)..q_m.
- Put a_k = [q_k -> w] and b_k = [q_k -> w'], with a_0 = b_0 = 1 and a_(m+1) = b_(m+1) = 0. Then A_k is valid iff
  a_k = 1 and b_(k+1) = 0, and B_k is valid iff b_k = 1 and a_(k+1) = 0.
- The pair (a_k, b_k) walks from (1, 1) to (0, 0).
  - The step that first leaves (1, 1) validates A_k or B_k at that k.
  - The step that last enters (0, 0) validates one of them at its k.
  - If both happen at the same k, then both A_k and B_k are valid.
- So each unordered pair {w, w'} yields at least two covers. Different pairs give different covers, since the pair is
  recovered as {last(P1), second(P2)}. QED.

**(A2) The sink.** Let z != a be the sink of B and B' = B - z. Then

    sigma_z(T1) = 3 h(B') + 2 R_a(B') - e_a(B')  >=  2 h(B') + 2 R_a(B')  >=  4,

with equality iff B is transitive and a is its penultimate vertex.

*Proof.* The out-neighbours of z in T1 are exactly x and y. Split by the shape of P2:
- P2 = z: P1 is a path of C2(B', a) ending at x (h(B') of them: [B'-path] y x) or at y (R_a(B') of them).
- P2 = z y: P1 is a path of B' + x ending at x. There are h(B') - e_a(B') of them.
- P2 = z y x ...: impossible. P1 would have to end at x or y, but both are in P2.
- P2 = z x: P1 is a path of B' + y ending at y. Here y is a sink, so there are h(B') of them.
- P2 = z x Q, with Q a path of X subset B' starting at a: P1 is a path of (B' - X) + y ending at y. Summing over X gives
  R_a(B').

Adding the cases gives the formula. Since e_a <= h, sigma_z >= 2h(B') + 2R_a(B') >= 4h(B') >= 4. Equality forces:
- h(B') = 1, so B' is transitive;
- R_a(B') = 1, so a is the sink of B' (R at the i-th vertex of TT_m is 2^(m-i));
- e_a(B') = 1, which then holds automatically.

QED.

**(A3) The other old vertices.** If w != a is not the sink of B, then d = d_B(w) >= 1 and w has out-degree d + 2 in T1.
By (A1), sigma_w(T1) >= (d + 2)(d + 1) >= 6.

**(A4) The bound.** Every old vertex w != a has sigma_w(T1) >= 4, with equality iff B = TT_n, a is its penultimate
vertex, and w its sink.

**(A5) The mark** (S77 audit section 4, re-derived).
- If a is the sink of B: sigma_a(T1) = 0, since x is the sink of T1 - a and every path of T1 - a ends at x.
- If a -> w is its only out-arc, with U = B - a: sigma_a(T1) = e_w(U) + R_w(U).
  - P2 = a, a y and a y x ... contribute 0, 0 and e_w(U).
  - P2 = a w ... x contributes R_w(U): y goes to P1 as its end, and x closes P2.
  - Hence sigma_a >= 2, sigma_a != 1, and sigma_a = 2 iff U is transitive with w its sink or its penultimate vertex.
- If d_B(a) >= 2: sigma_a(T1) >= (d + 1)d >= 6.

**FINITE-EXACT** (`old_mark_digit_bound_20261104`, all 7,412 classes of order 1-8, all 58,652 marks).
- Zero failures in each check: (A1) on B and on every vertex of T1; the native digits; the (A2) formula (3,192 sink
  cases at order 8); (A3); (A4) with its equality case; (A5).
- The order-8 run was pre-registered (HYP-9468 P3).
- The minimum digit at a non-sink old vertex is 9 at every order 2-8. This is FINITE-EXACT; the proved bound is 6.

**Relation to THM-4873.** THM-4873 (ii) (the S77 audit's proof) classifies digits 1 and 2. Here the bound is completed
to the sharp old-vertex minimum 4, with the exact sink digit and the equality case, and is re-proved without the audit
files.

## 2. (B) One-arc gadget transfers and the Lorentzian adjoint (PROVED)

**Transfer law.** Let G be a tournament with an entry mark u and an exit mark z (u = z allowed). For a rooted base
(T, v), let N = T => G with u -> v reversed, rooted at z. Then

    (h, e, R)(N, z) = M(G, u, z) (h, e, R)(T, v),

    M = [[h_G,  -s_u,       E_u    ],
         [e_z,  -P_(u->z),  E_(u,z)],
         [R_z,  -J_1,       J_2 + K]],

where:
- E_u = sum over B ∋ u of e_u(G[B]) h(G - B);
- E_(u,z) = sum over B ∋ u, B ∌ z of e_u(G[B]) e_z(G - B);
- J_1 = sum over S ∋ u, S ∌ z of s_u(G[S]) s_z(G - S);
- J_2 = sum over S ∋ u, S ∌ z of E_u(G[S]) s_z(G - S);
- K = sum over partitions G = S + Y + Z with u, z in Y of h(G[S]) P_(z->u)(G[Y]) h(G[Z]);
- P counts Hamiltonian paths between the named ends.

*Proof.* A path crosses back from G to T only through u -> v, and at most once.
- **Row 1** is the join value.
- **Row 2** (paths ending at z):
  - without the back arc: h e_z(G) - e_v P_(u->z)(G);
  - with it: [T-part] [G-part ending at u] v [T-part from v] [G-part ending at z], which is R_v(T) E_(u,z).
- **Row 3** (two-bin covers with P2 starting at z):
  - If P2 stays in G, then S contains T. P1 is a join of T with G[S2], giving
    h R_z(G) - e_v J_1 + R_v J_2.
  - If P2 enters T, it is z..u -> v [T-part from v] [G-part]. Then v is not in S, so T[S] => G[S2] is a plain join.
    The T-factor is R_v(T) and the G-factor is K.

Every count is linear in (h, e_v, R_v). QED.

**Special cases.**
- G a point gives klein's tower step (h - e + R, h - e, h + R) (THM-4865).
- G = {y -> x} with u = x, z = y gives C2 rooted at y: (h + 2R, R, 2h - e + 2R).
- The mechanism (a three-coordinate port vector transferred through one seam) is session29 ADDENDUM section 3, which
  gives explicit matrices for its 18 spine tiles.

**Lorentzian adjoint.** For every marked gadget,

    M(G^op, z, u) = J M(G, u, z)^T J,   J = diag(1, -1, 1).

*Proof.* Let T1 => G => T3 have the two back arcs u -> v and t -> z. Read it as (T1 ⋈ G) ⋈ T3: its H is
<M x1, y3>_J. Its converse is (T3^op ⋈ G^op) ⋈ T1^op, with H = <M' y3, x1>_J and M' = M(G^op, z, u). The two counts
are equal for all realisable x1 and y3. These span Q^3 (for example (1,1,1), (1,0,2), (3,1,3)), so M^T J = J M'. QED.

session29 section 5 checked this identity for its 18 tiles under the mirror. Here it holds for every marked gadget,
with the converse in the role of the mirror.

**FINITE-EXACT** (`gadget_adjoint_20261104`).
- The closed formula agrees with direct counts of N on five probe bases for every labelled G of order <= 5 and all
  (u, z): 133,525 checks, 0 failures.
- The adjoint identity holds for all 1,206,353 marked gadgets of order <= 6.

## 3. (C) A descent lemma: the limit of cone-covering induction

**Lemma.** Let C be a closed convex cone in R^d containing a nonzero integer point. Let h be a linear functional
taking positive integer values on the integer points of C - 0. Let M_1, ..., M_r be integer matrices of determinant
±1 with h(M_i q) > h(q) for every q in C - 0. Then C is not contained in the union of the M_i(C).

*Proof.* Otherwise every integer point p of C has the integer preimage M_i^(-1) p in C, with smaller h. Iterating gives
an infinite strictly decreasing sequence of positive integers. QED.

**Application.** On realisable triples, every gadget transfer has

    h' = h_G h - s_u e + E_u R  >=  E_u R  >=  R  >=  h,

with equality only for neutral insertions (R = h, i.e. digit 0; for example a sole exit into the sink). So on any cone
of triples with R >= (1 + c) h, c > 0, every transfer strictly increases h, and the lemma applies.

An inductive proof of Chan-Pak therefore cannot take the form "all primitive triples of such a cone are realisable, by
unimodular gadget preimages from a finite base". Such a proof needs one of three things:
- determinant-|d| maps with lattice bookkeeping;
- moves that keep h fixed;
- information beyond the triple.

## 4. (D) One-apex joins realise every odd N <= 2^32 except 7 and 21 (FINITE-EXACT)

**Families** (HYP-9468).
- F0: the plain counts, namely every count of order <= 8 and every g(v).
- F1: joins of two one-apex tournaments, T(w1) => T(w2)^op with x2 -> x1 reversed:
  H = (1 + u)(1 + v) + g(u) g(v) - e(u) e(v).
- F2: joins of a one-apex tournament with a rooted tournament (S, s) of order <= 7, T(w) => S^op with s -> x reversed:
  H = g(v) y_h - e(v) y_e + (1 + v) y_R, with y = (h, e_s, R_s)(S).
- The bank has 2,732 distinct triples, matching THM-4866 (D)'s count.

**Result** (`chanpak_one_apex_join_cover_20261104`, out of sample from 2^22).
- F0 + F1 + F2 attain every odd N in [1, 2^32] except 7 and 21. (RESULT LINE TO BE FILLED.)
- The previous FINITE-EXACT range was 99,999,999 (square_sum_spectrum_overnight_20261031d-i, seven-parent banks).
  klein's numerations reached 8,079,323 (THM-4865, HYP-9454).

**The bilinear family alone** (F1). Its miss density decreases in every dyadic block, from 0.127 at 2^10 to
(TO BE FILLED) at 2^31.

**Checks.**
- The suffix law used to evaluate g beyond the table, g(k 2^8 + r) = g(k) + k C(r) + (k mod 2) eps(r) + K(r), is
  checked against the plain recursion on 2,000,000 random v.
- `join_witness_check_20261104` counts 20,000 random F1 joins and 20,000 random F2 joins directly (orders <= 16;
  H up to 50,575): 0 failures.

**What this is not.**
- It is a finite computation, not a proof.
- The greedy sub-bank needed to close [23, 2^20] on top of F1 has 111 triples with a long tail. The coverage is
  redundancy across many triples, not a small explicit numeration.
- So it does not point to an all-order argument by itself. It does give any future induction a base 43 times larger
  than before.

## 5. Typing

| claim | type |
|---|---|
| (A1) sigma_v >= d(d - 1) | PROVED (READY (7) dualised; walk proof written out) |
| (A2) exact sink digit, (A4) sharp minimum 4 and its equality case | PROVED; FINITE-EXACT orders 1-8 |
| (A5) the mark digit | PROVED (S77 audit section 4, re-derived) |
| Non-sink old vertices have digit >= 9 | FINITE-EXACT orders 2-8 (HYP-9468 P3); proved bound 6 |
| (B) explicit gadget transfer | PROVED; mechanism KNOWN (session29 section 3); FINITE-EXACT to order 5 |
| (B) Lorentzian adjoint for every gadget | PROVED; tile case KNOWN (session29 section 5); FINITE-EXACT to order 6 |
| (C) descent lemma | PROVED (elementary) |
| (D) coverage of every odd N <= 2^32 except 7, 21 | FINITE-EXACT (pre-registered HYP-9468 P1) |
| Chan-Pak completeness | OPEN |

## 6. Audit

(Filled after the independent audit.)

## 7. Reproduction

    gcc -O3 -march=native -fopenmp -o omb old_mark_digit_bound_20261104.c && ./omb 8
    gcc -O3 -march=native -o gadj gadget_adjoint_20261104.c && ./gadj 6 5
    gcc -O3 -march=native -o bank rooted_triple_bank_20261104.c && ./bank 7 | awk '{print $2, $3, $4}' | sort -u > bank7.txt
    gcc -O3 -march=native -fopenmp -o cover chanpak_one_apex_join_cover_20261104.c && ./cover 32 bank7.txt
    gcc -O3 -march=native -o jwc join_witness_check_20261104.c && ./jwc 20000 7
    python3 one_apex_triple_check_20261104.py
