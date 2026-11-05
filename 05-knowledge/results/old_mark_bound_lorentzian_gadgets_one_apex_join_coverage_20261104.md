# The old-mark digit lower bound, one-arc gadget transfers, and one-apex joins up to 2^32

**Session:** opus-2026-11-04-S78. Worktree `C:/tmp/wt-s78`. Canon: THM-4875. Pre-registration: HYP-9468 (commit
19562af223). Independent audit: PASS WITH CORRECTIONS (section 6); corrections applied, MISTAKE-716.

**Owner prompt (verbatim).** "start another session proving the old-mark digit lower bound, work to complete chan-pak's
conjecture"

**Target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is
the number H(T) of Hamiltonian paths of some tournament. It is OPEN. THM-1370 proves that 7 and 21 are never attained.

**Results.**
- **(A) Old-mark digits, sharp at every order.**
  - The lower bound itself is KNOWN: THM-4873 (ii), from the S77 audit, proves the sink digit >= 4, every other old
    non-mark digit >= 6, and the mark formula.
  - New here: the exact sink digit and its equality case.
  - New from this session's audit (U1): every non-sink old vertex other than the mark has digit >= 9, with equality
    iff B is transitive, a is its sink and the vertex is its penultimate. This proves HYP-9468 P3's all-order
    conjecture.
- **(B) One-arc gadget transfers.**
  - The general transfer and the Lorentzian adjoint M(G^op, z, u) = J M^T J are KNOWN (SEAM.md, 2026-10-29;
    golden_hamilton RECURRENCE.md (G1)-(G4)).
  - New here: closed subset-sum formulas for the nine matrix entries, and larger finite checks.
- **(C) A descent lemma** (PROVED, elementary). In the audit's corrected form, it excludes only inductions that use
  unimodular gadgets alone to cover a full cone-intersect-coset of triples. Such a target is never the realisable set,
  so this is a weak no-go.
- **(D) FINITE-EXACT coverage** (pre-registered; independently reproduced by the audit). Every odd N <= 2^32 except 7 and
  21 is the count of a one-arc join of a one-apex tournament with either a one-apex tournament or a rooted tournament
  of order <= 7. The previous verified range was 99,999,999.

## 0. Notation and inherited facts

**Root response and digit.** For a tournament T and vertex t:
- R_t(T) = sum over X subset V - t of h(T[X]) s_t(T[V - X]), with h(empty) = 1;
- the digit is sigma_t = R_t - h.

R_t is klein's B_t, session29's R and SEAM's R_end on the converse. The rooted triple is (h, e_t, R_t), where e_t
counts the Hamiltonian paths ending at t.

**Cover reading of the digit.** R_t counts the two-bin covers (P1, P2) with P2 starting at t (P1 may be empty). The
covers with P1 empty, or with last(P1) -> t, are exactly the cuts of the Hamiltonian paths of T just before t. So

    sigma_t = #{covers with P1 nonempty and last(P1) in N+(t)} = E_t + #{such covers with |P2| >= 2},

where E_t = h(T - t) - e_t(T) counts the Hamiltonian paths of T - t that end in N+(t), i.e. the covers with P2 = (t).

**C2 child** (RESULT law (4); THM-4865). T1 = C2(B, a) adds x with x -> a only (every other vertex of B beats x), then
y, beaten by all of B, with y -> x. Its count is H(T1) = h + 2R_a, and its native digits are sigma_x = 0 and
sigma_y = h - e_a (RESULT eqs. (8), (10), (12)).

**One-arc join** (session29 ADDENDUM sections 2 and 5; golden_hamilton (G1); THM-4866 (B)). Let J = A => B with b -> a
reversed. Then

    H(J) = x_h y_h - x_e y_e + x_R y_R,   x = (h, e_a, R_a)(A),   y = (h, e_b, R_b)(B^op).

**One-apex tournaments.** For a binary word w of value v, T(w) is the chain 1 -> ... -> n plus an apex x with x -> i
iff bit i is 1. The apex triple is

    (g(v), e(v), 1 + v),   g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2),   e(v) = [v = 0 or v even].

Credits:
- the state recurrence is session25's;
- session179 spectrum/REPORT.md (1) proves the triple (state map (s, p, h) -> (2s + b, b, h + (1 - b)(s + p)));
- g is session124's digital function f (session124_digital_constructor_inverse_20261031/PROOF.md);
- the direct count is Session191 one_apex/PROOF.md eq. (2);
- THM-4876 (6) gives the same family as the one-centre tournaments TT_m + u.

Checked here by brute force for every word of length <= 7 (`one_apex_triple_check_20261104.py`), and by the audit
for every word of length <= 12.

## 1. (A) Old-mark digits

Let T1 = C2(B, a), with |B| >= 2 and h = h(B).

**(A1) Out-degree bound** (KNOWN: THM-4873 (ii) (D), which is READY (7) dualised). For every vertex v,

    sigma_v >= E_v + 2 * sum over pairs {p, q} subset N+(v) of h(T - {v, p, q}),

and hence sigma_v >= d(d - 1) for out-degree d.

*Proof.*
- Fix a pair {p, q} subset N+(v), and a Hamiltonian path q_1 ... q_m of T - {v, p, q} (empty if n = 3).
- For k = 0..m consider two covers:
  - A_k: P1 = q_1..q_k p, P2 = v q q_(k+1)..q_m;
  - B_k: P1 = q_1..q_k q, P2 = v p q_(k+1)..q_m.
- Put a_k = [q_k -> p] and b_k = [q_k -> q], with a_0 = b_0 = 1 and a_(m+1) = b_(m+1) = 0. Then A_k is valid iff
  a_k = 1 and b_(k+1) = 0, and B_k is valid iff b_k = 1 and a_(k+1) = 0.
- The pair (a_k, b_k) walks from (1, 1) to (0, 0).
  - The step that first leaves (1, 1) validates A_k or B_k at that k.
  - The step that last enters (0, 0) validates one of them at its k.
  - If both happen at the same k, then both A_k and B_k are valid.
- So each (pair, path) yields at least two covers. They have P1 nonempty, last(P1) in N+(v) and |P2| >= 2. The pair
  {last(P1), second(P2)} and the path (P1 minus its end)(P2 minus its first two) can be recovered from the cover. So
  the covers are distinct, and the cover reading of section 0 gives the bound. QED.

**(A2) The sink: the exact digit** (PROVED; the bound >= 4 is KNOWN from THM-4873 (ii)). Let z != a be the sink of B
and B' = B - z. Then

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

QED. (THM-4876 (7) gives e_a(B') <= (h(B') - 1)/2 when a is not the sink of B', which strengthens the middle
inequality there.)

**(A3) The other old vertices: digit >= 9** (PROVED by the audit, U1; checked here). Let w != a be a non-sink vertex of B,
d = d_B(w) >= 1. Then sigma_w(T1) >= 9, with equality iff B is transitive, a is its sink and w its penultimate vertex.

*Proof* (audit section 5).
- If d >= 2, w has out-degree d + 2 in T1, so (A1) gives sigma_w >= (d + 2)(d + 1) >= 12.
- Let d = 1, with w -> c the only out-arc of w in B, and U = B - w. Then N+_T1(w) = {c, x, y} and T1 - w = C2(U, a).
  - E_w = e_c(C2(U, a)) + h(U) + R_a(U), since C2(U, a) has h(U) paths ending at x and R_a(U) ending at y.
  - The three pairs give h(U) (pair {x, y}), h(U - c) (pair {c, x}; y is a sink of the rest) and h((U - c) + x) >= 1.
  - Hence sigma_w >= e_c(C2(U, a)) + R_a(U) + 3h(U) + 4.   (*)
- **Case c = a.** The paths [U - a] y x a give e_a(C2(U, a)) >= 1, so (*) >= 9.
  - Equality forces h(U) = R_a(U) = 1: U is transitive with sink a, so B is transitive with a its sink and w
    penultimate.
  - For that B a direct count gives exactly 9.
- **Case c != a, h(U) >= 3.** R_a(U) >= h(U), so (*) >= 16.
- **Case c != a, U = u_1 -> ... -> u_m transitive, c = u_j, a = u_k.**
  - If k < j: the path [U - {a, c}] y x a c gives e_c >= 1, and R_a(U) = 2^(m-k) >= 2, so (*) >= 10.
  - If k > j: the cover P1 = u_1 ... u_(j-1) y x, P2 = w u_j ... u_m is valid and is not counted in (*). Its
    would-be path u_1 ... u_(j-1) y u_(j+1) ... u_m is not a path, since j < m and y beats nobody in U.
    - If a is not the sink of U, this gives sigma_w >= 10.
    - If a = u_m, the exact count is sigma_w = 4 + 2[j = m - 1] + 3 * 2^(m-j) >= 12 (audit U4).

QED.

**(A4) Summary for the old vertices other than the mark.**
- The sink of B (when it is not a) has the exact digit (A2): at least 4, and 4 iff (TT_n, a penultimate).
- Every other old vertex has digit at least 9, and 9 iff (TT_n, a its sink, w penultimate).
- On TT_n the two minimal configurations are exchanged by swapping the roles of a and w.

**(A5) The mark** (KNOWN: THM-4873 (ii); a corollary of RESULT (12), sigma_a(T1) = 2R_a(B) - h(B) - h(B - a)).
- If a is the sink of B: sigma_a(T1) = 0.
- If a -> w is its only out-arc, with U = B - a: sigma_a(T1) = e_w(U) + R_w(U). This is never 1, and it is 2 iff U is
  transitive with w its sink or its penultimate vertex.
- If d_B(a) >= 2: sigma_a(T1) >= (d + 1)d >= 6, by (A1).

**FINITE-EXACT.**
- `old_mark_digit_bound_20261104`: all 7,411 classes of order 2-8 and their 58,652 marks.
  - Zero failures in each check: (A1) on B and on every vertex of T1; the native digits; the (A2) formula (3,192 sink
    cases at order 8); (A3)'s out-degree bound; (A4) with its digit-4 equality case; (A5).
  - The minimum non-sink old digit is 9 at every order 2-8.
  - The order-8 run was pre-registered (HYP-9468 P3).
- `nonsink_digit_nine_check_20261104.py` (after the audit), on every labelled B of order 2-6, every mark a and every
  non-sink w != a (972,230 cases):
  - minimum 9;
  - the digit-9 cases are exactly the n! labelled transitive B with a the sink and w penultimate;
  - (A3)'s residual formula holds for 2 <= m <= 8 and every j < m.
- The audit's independent class generation (A000568 counts) confirms all of (A) at orders 1-8 and RESULT (12) on all
  58,653 pairs (B, a).

## 2. (B) One-arc gadget transfers: closed entries (transfer and adjoint KNOWN)

**Prior art.** SEAM.md (`seven_return_native_recursions_20261029/SEAM.md`, 2026-10-29) proves the general statement:
- for every block Y with arbitrary entry and exit ports there is an integer 3x3 matrix with v(X ⊙ Y) = M_Y v(X)
  (eqs. (1)-(3));
- composition M_(X⊙Y) = M_Y M_X (eq. (5));
- the adjoint M_(X*) = J M_X^T J, J = diag(1, -1, 1) (eq. (6)), proved by the same two-sided reading and spanning
  vectors;
- the singleton matrix [1 -1 1; 1 -1 0; 1 0 1], which is the tower step below.

SEAM's v = (H, S, R_end) on X is this note's (h, e, R) on X^op, and X ⊙ Y is the converse of T => G. The general
transfer is also golden_hamilton_transfer_20261103/RECURRENCE.md (G1)-(G4). session29 ADDENDUM section 3 had the 18
spine tiles. This section was written without these citations; the audit supplied them (MISTAKE-716).

**Closed entries (new presentation).** Let G have an entry mark u and an exit mark z (u = z allowed). For a rooted base
(T, v), let N = T => G with u -> v reversed, rooted at z. Then (h, e, R)(N, z) = M(G, u, z)(h, e, R)(T, v), with

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
- **Row 3** (covers (P1, P2) with P2 starting at z). Write G[S] for the G-part of P1.
  - **P2 stays in G**, on Y = G - S. Then P1 is a Hamiltonian path of T => G[S]. That is a plain join if u is in Y and a
    one-arc join if u is in S. Summing over Y gives h R_z(G) - e_v J_1 + R_v J_2.
  - **P2 enters T.** It is (z ... u) v [T-part] [G-part G[Z]], and P1 = [T-part T[X]] [G[S]] uses no back arc. The
    T-factor is sum over X of h(T[X]) s_v(T - X) = R_v(T), and the G-factor is K.

Every count is linear in (h, e_v, R_v). QED.

**Special cases.**
- G a point gives klein's tower step (h - e + R, h - e, h + R) (THM-4865; SEAM's singleton matrix).
- G = {y -> x} with u = x, z = y gives C2 rooted at y: (h + 2R, R, 2h - e + 2R).

**FINITE-EXACT** (`gadget_adjoint_20261104`).
- The closed formula agrees with direct counts for every labelled G of order <= 5 and all (u, z): 133,525 checks,
  0 failures.
  - The five probe bases realise only four distinct triples ((1, 0, 2) occurs twice).
  - The audit fitted M from three spanning probes for all 26,705 marked gadgets of order <= 5. The fit equals the
    closed formula, and linearity holds on 456,745 further base checks.
- The adjoint identity holds for all 1,206,353 marked gadgets of order <= 6 (the audit: 0 failures, fitted to order 5).

## 3. (C) A descent lemma, and what it does and does not exclude

**Lemma** (PROVED, elementary).
- Hypotheses:
  - C is a closed convex cone in R^d containing a nonzero integer point;
  - h is a linear functional with nonnegative integer values on C ∩ Z^d;
  - M_1, ..., M_r are integer matrices of determinant ±1 with h(M_i q) > h(q) for every q in C - 0.
- Conclusion: some nonzero integer point of C lies outside the union of the M_i(C).

*Proof.* Otherwise every nonzero integer point p of C has the nonzero integer preimage M_i^(-1) p in C, with smaller h.
Iterating gives an infinite strictly decreasing sequence of nonnegative integers. QED.

**Gadgets raise h on {0 <= e <= h}.** For a triple with 0 <= e <= h,

    h' - h = (h_G h - s_u e) + (E_u - 1) R + (R - h),

and the three terms are nonnegative on realisable triples.
- **Equality.** h' = h iff the base port is the sink of T (e = h = R) and G is transitive entered at its source.
  - By THM-4876 (7), a non-sink port has e <= (h - 1)/2.
  - Earlier wording here ("equality only for neutral insertions, R = h") gave a necessary condition only.
  - The audit's check found exactly 13,661 equality cases among 72,958,060 (gadget of order <= 5, realisable triple)
    pairs.
- **Strict increase.** On a cone inside {0 <= e <= h, R >= (1 + c) h} with c > 0, every transfer strictly increases h.
  This includes the apex (0, 0, 1), where h = 0.

**Corrected consequence** (audit U2).
- Hypotheses:
  - C has nonempty interior and lies inside {0 <= e <= h, R >= (1 + c) h};
  - M_1, ..., M_r are gadget matrices of determinant ±1;
  - Λ is a coset of a full-rank lattice.
- Conclusion: for every finite F, C ∩ Λ is not contained in F ∪ M_1(C ∩ Λ) ∪ ... ∪ M_r(C ∩ Λ).

*Proof.* The lemma gives a point p0 of C outside the closed set ∪ M_i(C). A ball around p0 misses that set and meets the
interior of C. The open cone it generates misses every M_i(C) and contains infinitely many points of Λ. QED.

**What this excludes, and what it does not.**
- The earlier sentence "a finite-base unimodular cone-covering induction for Chan-Pak is impossible" is withdrawn.
  - It handled no finite base.
  - Its target was false: (2, 1, 3) has even h, and (3, 3, 4) has e = h with R > h. So no cone-intersect-coset of
    triples is realisable, and excluding a proof method for such a target says nothing about Chan-Pak.
  - Only 96 of 1,024 order-4 gadgets (9.4%) are unimodular, 1,200 of 25,600 at order 5 (4.7%) and 20,160 of 1,179,648
    at order 6 (1.7%).
- What survives is the heuristic that already motivated THM-4865 (f), HYP-9453 section 2 and klein's root-response
  reflection ("local moves cannot do it"). An induction has to track the realisable set itself, or use moves that keep
  h fixed, maps of determinant |d| > 1, or information beyond the triple.

## 4. (D) One-apex joins realise every odd N <= 2^32 except 7 and 21 (FINITE-EXACT)

**Families** (HYP-9468).
- F1: joins of two one-apex tournaments, T(w1) => T(w2)^op with x2 -> x1 reversed:
  H = (1 + u)(1 + v) + g(u) g(v) - e(u) e(v).
- F2: joins of a one-apex tournament with a rooted tournament (S, s) of order <= 7, T(w) => S^op with s -> x reversed:
  H = g(v) y_h - e(v) y_e + (1 + v) y_R, with y = (h, e_s, R_s)(S).
- The bank has 2,732 distinct triples, matching THM-4866 (D)'s count and the audit's independent rebuild.
- F0, the plain counts (every count of order <= 7 and every g(v)), was in the pre-registered run but is redundant
  (audit U3).
  - g(2t) = F1(0, t) and g(2t + 1) = g(t), so every g value is an F1 value (also Session 193 one_arc_joins/PROOF.md
    section 3, a concurrent post-checkpoint packet).
  - The audit's run without F0 shows that F1 + F2 already cover every odd N <= 2^28 except 7 and 21, which includes all
    bank counts (at most 189).

**Result** (`chanpak_one_apex_join_cover_20261104`, out of sample from 2^22; HYP-9468 P1 CONFIRMED).
- F1 + F2 attain every odd N in [1, 2^32] = [1, 4,294,967,296] except 7 and 21.
- The run took 21.7 minutes on 12 threads.
- The audit reproduced this to 2^32 with its own g, bank and sieve: misses exactly {7, 21}; 47,964,383,242 F1 pairs;
  all 22 F1-only block counts equal.
- The previous FINITE-EXACT range was 99,999,999 (square_sum_spectrum_overnight_20261031d-i, seven-parent banks).
  klein's numerations reached 8,079,323 (THM-4865, HYP-9454).
- Session 193's independent check reaches 2^20. It adds the gcd branch filter and the finite complementarity of F2 on
  the prime and 7p lanes.

**The bilinear family alone** (F1).
- Its miss count per dyadic block grows by a ratio between 1.46 and 1.70 for k = 10..31, and between 1.55 and 1.63 for
  k = 22..31.
- The density falls in every block, from 0.127 in [2^10, 2^11) to 0.00087 in [2^31, 2^32) (HYP-9468 P2, CONFIRMED). So
  F1 misses about X^0.65 of the odd numbers up to X.
- Context. Deleting both apexes leaves a transitive tournament. So F1 lies inside the two-feedback-vertex (two-apex)
  family.
  - That family's completeness is already an OPEN conjecture (square_sum_spectrum_overnight_20261030/README.md, with
    coverage to 99,999).
  - Session191 two_apex/TWO_APEX.md covers it to 230,519.
- New conjecture (a density statement for this subfamily): almost every odd number is the count of a join of two
  one-apex tournaments.

**Residue profile of the F1 misses** (exploratory; `chanpak_exploratory_20261104/f1_miss_residues.c`, on
(2^22, 2^24]; the audit's own g reproduces it).
- **Powers of 2.** The misses are uniform modulo 4, 8, 16 and 32, each class at 0.0066 ± 0.0004.
- **Modulo 3.** N = 1 mod 3 is missed at rate 0.0124, against 0.0036 and 0.0039 for the classes 0 and 2. Modulo 9, the
  classes 1, 4, 7 carry rates near 0.0124.
- **Interpretation.** The one-apex joins are 2-adically balanced; this matches Session 192's dyadic residue
  programming. Their deficit is 3-adic. The order-<=7 bank closes every gap (F2), but which triples do so, and
  whether they act 3-adically, was not analysed.

**Checks.**
- **Suffix law.** g is evaluated beyond the table by g(k 2^8 + r) = g(k) + k C(r) + (k mod 2) eps(r) + K(r).
  - This is a theorem: Session191 one_apex/PROOF.md section 2, and the audit's derivation from the closed form.
  - It was also checked against the plain recursion on 2,000,000 random v.
- **Random joins.** `join_witness_check_20261104` counts 20,000 random F1 joins and 20,000 random F2 joins directly
  (orders <= 16, H up to 50,575): 0 failures.
- **Explicit witnesses.** The audit built 41 witness tournaments (orders 14-26, H up to 40,376,785) and counted their
  Hamiltonian paths directly: 0 mismatches.

**What this is not.**
- It is a finite computation, not a proof.
- The greedy sub-bank needed to close [23, 2^20] on top of F1 has 111 triples with a long tail (the audit's rerun:
  `greedy_repro_20.out`). The coverage is redundancy across many triples, not a small explicit numeration.
- So it does not point to an all-order argument by itself. It does give any future induction a base 43 times larger
  than before.

## 5. Typing (after the audit)

| claim | type |
|---|---|
| (A1) sigma_v >= E_v + 2 sum h(T - {v, p, q}) >= d(d - 1) | KNOWN (THM-4873 (ii) (D); READY (7) dualised); proof written out |
| (A2) exact sink digit 3h(B') + 2R_a(B') - e_a(B') and its digit-4 equality case | PROVED (new); bound >= 4 KNOWN (THM-4873 (ii)) |
| (A3) non-sink old digit >= 9, equality iff (TT_n, a sink, w penultimate) | PROVED (audit U1); FINITE-EXACT orders 2-8 (HYP-9468 P3) |
| (A3) residual digit 4 + 2[j = m - 1] + 3 * 2^(m-j) | PROVED (audit U4); FINITE-EXACT m <= 8 |
| (A5) the mark digit | KNOWN (THM-4873 (ii); RESULT (12)) |
| (B) general transfer and Lorentzian adjoint | KNOWN (SEAM.md eqs. (1)-(6); golden_hamilton (G1)-(G4)) |
| (B) closed subset-sum entries of M(G, u, z) | PROVED (new presentation); FINITE-EXACT to order 5, adjoint to order 6 |
| (C) descent lemma; corrected consequence U2 | PROVED (elementary); the application is a heuristic remark |
| (D) every odd N <= 2^32 except 7, 21 is an F1 or F2 join count | FINITE-EXACT (pre-registered HYP-9468 P1; audit reproduction) |
| F1 miss density -> 0 | conjecture (FINITE-EXACT trend to 2^32) |
| Chan-Pak completeness | OPEN |

## 6. Audit

`04-computation/experiments/old_mark_one_apex_join_20261104_audit_S78/AUDIT.md`: **PASS WITH CORRECTIONS**. The audit
reproduced all four claims with its own code:
- the bank, rebuilt from all 2,131,019 labelled tournaments of order <= 7;
- the one-apex triple to length 12;
- the 2^32 cover;
- 41 explicit witnesses;
- (A) at orders 1-8, from its own class generation;
- (B) by probe fitting.

It confirmed that both pre-registered outputs postdate commit 19562af223.

**Majors (all applied).**
- **M1.** (B)'s general transfer and adjoint were KNOWN (SEAM.md 2026-10-29; golden_hamilton (G1)-(G4)). Only the
  closed entries are new.
- **M2.** (A)'s bounds (>= 4 at the sink, >= 6 elsewhere, the mark) were already THM-4873 (ii). THM-4873 was misdescribed
  as classifying only digits 1 and 2. New: the exact sink digit and the equality cases. (A5) is a corollary of
  RESULT (12).
- **M3.** (C)'s consequence handled no finite base and targeted a false statement. Its hypotheses were unstated, only a
  small minority of gadgets are unimodular, its equality clause was misstated, and its precursors were uncited.
  Replaced by section 3's corrected form.

**Minors (all applied).**
- The class and mark counts mixed order ranges.
- The five probes realise four triples.
- The suffix law is a theorem.
- One-apex credits were missing: session179, session124 and Session191.
- F1 sits inside the two-apex family, whose completeness conjecture and coverage are cited now.
- THM-4876 was uncited.
- (C)'s precursors were uncited.
- The home of (A1) is THM-4873 (ii).
- P2's ratio index: read literally it needs m_32; all nine available ratios lie in the window.
- The runtime was 21.7 minutes, not "about 30".
- The greedy output was uncommitted; the audit reproduced it.
- The row-3 wording reused S.
- F0 is redundant.
- "Proved bound 6" is superseded by U1.

**Upgrades (proved by the audit; adopted).**
- U1: the non-sink >= 9 theorem with its equality case.
- U2: the corrected consequence of (C).
- U3: F0 is redundant.
- U4: the residual digit formula.

## 7. Reproduction

    gcc -O3 -march=native -fopenmp -o omb old_mark_digit_bound_20261104.c && ./omb 8
    python3 nonsink_digit_nine_check_20261104.py
    gcc -O3 -march=native -o gadj gadget_adjoint_20261104.c && ./gadj 6 5
    gcc -O3 -march=native -o bank rooted_triple_bank_20261104.c && ./bank 7 | awk '{print $2, $3, $4}' | sort -u > bank7.txt
    gcc -O3 -march=native -fopenmp -o cover chanpak_one_apex_join_cover_20261104.c && ./cover 32 bank7.txt
    gcc -O3 -march=native -o jwc join_witness_check_20261104.c && ./jwc 20000 7
    python3 one_apex_triple_check_20261104.py

The audit folder's README section (AUDIT.md section 6) gives the independent reproduction, including
`aud_cover -DNOF0` for U3.
