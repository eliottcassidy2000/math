# Independent audit S11: THM-4879, HYP-9471, the seam-gadget thickness note, and the proposed addendum (F)

**Audited files**
- `01-canon/theorems/THM-4879-two-port-transitive-seam-gadgets-closed-transfer-matrices-with-a-three-power-entry-and-the-seam-parity-r-equals-s-mod-2.md`
  (the working-tree version, which adds the (C) proof to commit 4b0b0c740b);
- `05-knowledge/hypotheses/HYP-9471-seam-gadget-thickness-by-piece-order-plane-covering-density.md`;
- `05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md` (called "the note" below);
- the proposed addendum (F) (`scratchpad/THM-4879_F_addendum.md`, not yet in the THM file).

**Method.** Every number below was recomputed with the programs in this folder (`aud_*`). None of the audited
programs was imported or run. They were read only for definitions (seeds, bands, box rule, D_k).

## Verdict: PASS WITH CORRECTIONS

- **THM-4879 (A)-(E).** Every statement is true and every proof is correct, apart from slips in the presentation of
  the (E) proof.
  - **(D) is not new.** Session29 proved it in the same (H, S, R) convention on 2026-10-22. It also follows from an
    exact identity in Session 193, committed at 18:47, before this session's first commit at 20:27. Retype (D) as
    KNOWN (M1).
- **Addendum (F): correct.** I checked the proof step by step and by brute force. Adopt it: it proves the OPEN part
  of (E) and turns P5 into a theorem. One sentence in its computer-check paragraph is wrong (see the (F) section).
- **HYP-9471.** Every pre-registered number reproduces exactly, and git confirms the pre-registration order. All six
  outcomes are reported correctly. There are two small factual slips (m5, m6).
- **The note.** Every computed number reproduces except one after-the-fact statistic (0.8825 against 0.8819).
  - Three of the four "structural" barriers (section 0 item 5, section 4 (b)-(d), section 6) claim more than is shown.
  - One of them, section 4 (c), is a FINITE-EXACT claim that is **false**. The cause is a bug in
    `small_gadget_structure.py` (M2).

## Findings

### MAJOR

**M1. THM-4879 (D) is KNOWN, and the credit is missing.**
- **Where it is claimed.** THM-4879 status ("(D) is a one-line consequence of Redei's theorem and the singleton
  seam"), statement (D), title and Credits. Also the note section 0 item 3 ("Seam parity (THM-4879 (D))") and
  section 7 ("PROVED: THM-4879 (A)-(E) (elementary; (A) KNOWN in substance)").
- **Prior statement 1: Session29.** `05-knowledge/results/session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md`,
  section 5 "Port parity and the Lorentzian seam" (commit 25a4f6ba21, 2026-10-22).
  - It uses the same definitions: H, S, E, Q and R = E + Q.
  - It proves H(X ⊙ X^op) = H² − S² + R² and concludes: "Since every tournament has odd H, reduction modulo two gives
    R≡S modulo two for every X and every a."
- **Prior statement 2: Session 193.** `05-knowledge/results/session193_triangulated_proof_20261105/multiple_triangles/PROOF.md`,
  section 2 "Elementary rooted-path parity identity" (commit 09370ca742, 2026-11-04 18:47). The first commit of this
  session, a98199e8bb, was at 20:27.
  - Section 2 proves O_a + V_a = R_a(A) and O_a − V_a = e_a(A). Here R_a counts covers whose second path *starts* at
    a, which is the converse convention.
  - Applied to X^op this reads R_a(X) − S_a(X) = 2V, an exact strengthening of (D).
- **Correction.**
  - Statement (D): append "(KNOWN: Session29 converse-congruences REPORT section 5, 2026-10-22, in this convention;
    also the mod-2 reduction of Session 193 multiple_triangles PROOF.md section 2, applied to the converse. The
    singleton-seam proof below is a shorter alternative.)"
  - Status: replace "(D) is a one-line consequence of Redei's theorem and the singleton seam" with "(D) is KNOWN
    (Session29, Session 193); the singleton-seam proof given here is new".
  - Add both sources to Credits.
  - In the note, section 0 item 3: "Seam parity (KNOWN: Session29; re-proved in THM-4879 (D))".
  - In the note, section 7: "(A) and (D) KNOWN".

**M2. Note section 4 (c) is false as a FINITE-EXACT statement.**
- **Where.** Note section 4 (c), section 0 item 5 third bullet, section 6 ("closed affine digit grammars, by
  section 4 (c)-(d)") and section 7 ("FINITE-EXACT: ... the eigen-direction and unipotent searches").
- **Quoted claim.** "Among the 281 gadget matrices of order <= 5, the only rational eigen-direction shared by two of
  them with nonzero eigenvalue is (2, 1, -1) on the right and (2, -1, -1) on the left, with eigenvalues ±1. So no
  line of states is preserved by two expanding gadgets."
- **The bug.**
  - `small_gadget_structure.py` takes as candidate eigenvalues the divisors of det. When det = 0 this is only λ = 0
    (`range(1, abs(D) + 1)` is empty).
  - So the nonzero eigenvalues of singular matrices are never examined. 147 of the 281 listed matrices are singular,
    and 68 (matrix, λ) pairs with λ ≠ 0 are missed among the 264 distinct matrices.
  - It also compares kernel *basis vectors* rather than intersecting eigenspaces.
  - My imitation of that procedure (`aud_structure.py`) reproduces its output exactly.
- **Complete search** (`aud_structure.py`; integer roots of the cubic, or of the quadratic cofactor when det = 0, and
  exact eigenspace intersections):
  - 10 right and 10 left rational eigen-directions are shared by two or more distinct matrices with nonzero
    eigenvalue, among the 264 *distinct* matrices of order <= 5.
  - 8 of them are expanding:

    | direction | matrices | eigenvalues |
    |---|---|---|
    | (1, 1, 1) | 12 | 1, 3, 5, 9, 11, 13, 15 |
    | (1, 0, 2) | 11 | 3, 5, 9, 15, 17, 27 |
    | (1, 0, 4) | 8 | 5, 9, 15, 17, 27 |
    | (3, 1, 3) | 5 | 3, 5, 9, 15, 17 |
    | (5, 1, 5) | 4 | 5, 9, 15 |
    | (1, 0, 8) | 4 | 9, 17, 27 |
    | (9, 1, 9) | 2 | 9, 27 |
    | (5, 2, 6) | 2 | 9, 17 |

  - The two non-expanding ones are (2, 1, −1) and (3, 0, −2), both with eigenvalues ±1. The audited search also
    missed (3, 0, −2).
  - All 8 expanding directions are rays through realised states of order <= 5.
  - Example: TT_3 with (u, z) = (1, 0) has M = [1 0 2; 1 0 2; 1 0 2] and M(1,1,1) = 3(1,1,1)
    (`aud_structure_examples.out`).
- **What survives.**
  - Most of these matrices are rank-one "reset" gadgets. On a line through the origin, linear maps act by pure
    scaling, with no digits.
  - `aud_affine_lines.py` finds no pair of distinct gadgets of order <= 5 with a common invariant affine line off
    the origin.
  - So the final sentence of (c) is true, but not for the reason given.
- **Correction.**
  - Replace (c) with: "Among the 264 distinct gadget matrices of order <= 5 (281 counting cross-order repeats), 10
    right and 10 left rational eigen-directions are shared by two or more matrices with nonzero eigenvalue. 8 of them
    are expanding; for example, the singleton ray (1, 1, 1) is fixed by 12 reset gadgets with eigenvalues 1, 3, 5, 9,
    11, 13 and 15. On such a ray the maps act by pure scalings. No two distinct gadgets of order <= 5 share an
    invariant affine line off the origin, so single small gadgets give no 1D affine digit grammar (audit S11)."
  - Fix the eigenvalue candidates in `small_gadget_structure.py`.
  - In section 0 item 5, replace "No two small gadgets share an expanding rational eigen-direction, so there is no
    invariant line for a 1D digit grammar" with "No two gadgets of order <= 5 share an invariant affine line off the
    origin (shared expanding eigen-rays exist, but on them the maps are pure scalings)".
  - In section 7, keep "eigen-direction search" FINITE-EXACT only after the rerun.

**M3. Note section 4 (d): the lemma is correct, but the impossibility claim goes far beyond what it proves.**
- **Where.** Note section 4 (d), section 0 item 5 fourth bullet ("A closed {3N, 3N + 2, 3N + 4} grammar is impossible
  in the triple calculus"), and section 6 ("Ruled out within the triple calculus: closed affine digit grammars").
- **What is correct.**
  - The ordinal sum X => C3 with the exit kept in X sends v to 3v (`aud_ordinal.py`: all 285 rooted tournaments of
    order <= 4).
  - The lemma holds: if the ×3 step is the *scalar* 3I and M_d is linear with H(M_d q) = 3H(q) + d on a set closed
    under 3I, then d = 0.
- **What it does not cover.** The note says "In the triple calculus the step N -> 3N is the ordinal sum X => C3 with
  the exit kept in X". That ignores two other ×3 steps that stay inside the triple calculus (checked by brute force):
  - C3 => X with the exit kept in X: v → (3H, 0, 12R). This is linear and **not scalar**.
  - C3 => X with the exit in C3: v → (3H, H, 3H).
- **Why the lemma does not reach them.**
  - Its proof uses only that M_d commutes with the scalar step. It says nothing about these steps, nor about digit
    maps whose affine law holds on part of D only, nor about state-dependent words.
  - "Seam words cannot supply a pure ×3 step" excludes only the H-row (3, 0, 0). It does not exclude a seam word T
    with H(Tq) = 3H(q) on a sub-cone.
- **Related proved result.** Session 195's `session195/rooted-interfaces/EQUAL_RESPONSE_BARRIER.md` proves a sharp,
  scoped statement in this direction: the counts at most 3h + 4 of the moving-root descendants of a state with R = h
  are exactly {h, 2h − e, 3h}. It should be cited as such.
- **Correction.**
  - Retitle (d) "The scalar ×3 step admits no affine digits".
  - Add: "This covers only the scalar step X => C3 with the exit kept in X. The other ordinal sum, C3 => X, stays
    inside the triple calculus (v -> (3H, 0, 12R) with the exit in X; v -> (3H, H, 3H) with the exit in C3) and is
    not covered."
  - In section 0 item 5: "A closed {3N, 3N+2, 3N+4} grammar built on the scalar step X => C3 is impossible".
  - In section 6, replace "closed affine digit grammars, by section 4 (c)-(d)" with "affine digit grammars on a common
    invariant line of two single gadgets of order <= 5 (section 4 (c), corrected), and ternary grammars built on the
    scalar step X => C3 (section 4 (d))".

**M4. Note section 0 item 5 and section 4 (b): "no orbit contains an arithmetic progression" is not supported.**
- **Quoted claims.** "No short gadget word is unipotent, so no orbit contains an arithmetic progression." and "No
  orbit in this range therefore contains an arithmetic progression, which would have given an easy proof."
- **What the search shows.** No word of length <= 3 over the 56 distinct matrices of order <= 4 has a characteristic
  polynomial with |det| = 1 and |tr|, |c2| <= 3. I reproduced the three hits exactly.
- **What it does not show.**
  - Read literally, with "orbit" meaning the semigroup orbit (as in item 5's first bullet), the claim is unproved.
  - It is also in tension with the session's own coverage data: the H-projection contains the progression 23, 25,
    ..., 19,999.
  - An infinite progression inside it is exactly what the spectrum conjecture asserts for seam words.
- **Extended check.** No such word has (t − 1)² dividing its characteristic polynomial, for any determinant
  (`aud_structure.out`). So no single short word W makes n → H(W^n q) an arithmetic progression.
- **Correction.** Replace with "No word of length <= 3 over the gadgets of order <= 4 is quasi-unipotent, or even has
  a double eigenvalue 1, so iterating one such word never gives an arithmetic progression of H-values."

### MINOR

**m1. THM-4879 (E), third bullet: the S-row-zero hypothesis is missing.**
- The bullet reads "The determinant of the restriction [[h(G), L_u(G)], [R_z(G), K1 + K2]] acting on (H, R) is
  ≡ M_(z,u)(G) (mod 2)."
- Prefix it with "For such G (S-row zero), ...".
- The hypothesis is necessary. Without it the congruence fails for 2, 4, 12, 73, 556, 7,867 and 171,897 of the
  distinct matrices of orders 2-8 (`aud_census.out`, last column).

**m2. THM-4879, proof of (E): two slips.**
- **The s_z step.** "on input (1, 1, 1) the output's S-coordinate is 0 ... this forces s_z(G) = 0". One input gives
  only s_z − P_(z->u) + L_(z,u) = 0. Use the row itself: "the S-row is (s_z(G), −P_(z->u)(G), L_(z,u)(G)) by (A), so
  S_z(G) = s_z(G) = 0".
- **The order of the inputs.** "The inputs (1, 1, 1) and (1, 0, 2) ... give R_z(G) ≡ 0 and R_z(G) - M + (K1 + K2)
  ≡ 0" has them swapped. Write: "give, respectively, R_z(G) − M + (K1 + K2) ≡ 0 and R_z(G) + 2(K1 + K2) ≡ R_z(G)
  ≡ 0".

**m3. THM-4879 "What this does and does not give", Parity: the sentence is inverted.**
- It reads "They forbid odd-determinant maps on P only for gadgets with M_(z,u)(G) odd, and no such gadget has been
  seen."
- Replace with "(E) shows that a plane map of odd determinant can only come from a gadget with M_(z,u)(G) odd; (F)
  shows no such gadget exists."
- Once (F) is adopted, also update:
  - the THM line "Whether M_(z,u)(G) is always even for S-row-zero gadgets is OPEN";
  - the note section 3 ("A proof for all orders is OPEN");
  - the note section 0 item 3 ("So no unimodular map exists on P", which was proved only to order 8 when written);
  - the note section 6 (ii) ("Every plane map of order <= 8 has even determinant").

**m4. THM-4879 cites "(FINITE-EXACT, `ttbfs` in the session's note)", which does not exist.**
- The script is `transitive_words_bfs.py` / `.out`.
- The claim itself holds: `aud_ttwords.py` gives the same 261 misses (first 7, 11, 19, 21, 35, 37) for every k_max
  from 10 to 30, under both rules for recording H (the audited script records H only for kept states; recording
  every image changes nothing).
- Larger k cannot add values <= 2,000. H' >= 2^(k-1-u) R forces k − 1 − u <= 10. A kept state needs R' <= 2,000,
  and R' >= 2^(z-1) H when z > u, while R' >= 2^z when z <= u. So beyond k = 22 every usable matrix already occurs at
  a smaller k.

**m5. HYP-9471: "The realised P-states have R/H >= 2.01 in the runs below" is false.**
- The minimum R/H is exactly 2.0000 in every 3D run, also among states with H >= 801. The seed (1, 0, 2) is an
  example.
- The pure-plane runs give minimum slopes 2.0011-2.0056 on H in [801, 1799].
- Replace with "R/H >= 2, with equality attained".

**m6. Counts that include cross-order repeats are presented as matrix counts.**
- The affected phrases:
  - HYP: "(all 1,941 matrices; ...)";
  - note section 0 item 4 and section 4 (a): "1,941";
  - section 4 (b): "the 61 gadget matrices of order <= 4";
  - section 4 (c): "the 281 gadget matrices of order <= 5";
  - `seam_bfs_p2_p3.out`: "mats=22441".
- The distinct counts are 1,884, 56, 264 and 22,243. The results do not depend on the difference.

**m7. Note section 2: two wrong sentences.**
- **"At order 3 the 12 matrices include the singleton, the transitive gadgets and the directed triangle."**
  - The 12 are the 9 TT_3 port pairs and 3 C3 port classes. The singleton matrix occurs only at order 1.
  - Also, C3 has a unimodular port pair: [3 −1 3; 1 0 1; 3 −1 4] (det 1; it equals the cube of the singleton
    matrix).
- **"Most gadget matrices are singular: 107 of the 220 at order 5 have det 0."**
  - 48.6% is not "most". The singular fractions at orders 2-8 are 75.0, 75.0, 63.6, 48.6, 32.7, 19.8 and 11.3%.
  - Suggest "Many small gadget matrices are singular (107 of 220 at order 5; 11% at order 8)".

**m8. Note section 4 (d): wrong formula.**
- It says "The other ordinal sum, with the exit moved into C3, needs H(X) + 2pc_2(X)".
- For X => C3 with the exit in C3, the state is (3H, 0, 3C2(X)), where C2(X) = sum_U H(U)H(X−U) = 2H + 2pc_2 (pc_2 =
  unordered two-path covers). So R = 6(H(X) + pc_2(X)).
- This is verified by brute force. The state is indeed not a function of v(X, a): v = (1,1,1) gives R = 6, 12, 24
  and 48.

**m9. The after-the-fact statistic 0.8825 does not reproduce.**
- It appears in the HYP Results and status ("0.88") and in note section 5.
- My exact mean of D_8(R/H) over the 54,994 order-<=8 plane states with H in [801, 1799] is **0.8819**. Binned
  variants give 0.8819-0.8820, and the mean over distinct slopes gives 0.8744 (`aud_slope_variants.out`).
- The other slope statistics reproduce exactly: the state count, median 2.64, 90th percentile 4.15, and the
  histogram.
- No committed script computes this number. Replace it with 0.8819, or commit the script.

**m10. The "threshold" framing.**
- HYP title: "D_k(t) >= 1 is the necessary density condition for a finite-alphabet self-covering of a 2D cone of
  states". Note section 0 item 6: "It crosses 1, the necessary condition for a finite alphabet to fill a plane of
  states with positive density, at order 8."
- **The sharp necessary condition.** A counting argument (header of `aud_cone_condition.py`) gives it for
  self-covering the cone t0 <= R/H <= t1:

  (1/(t1 − t0)) Σ_i ∫ (a_i + b_i t)^-2 dt >= 1, where the integral runs over {t in [t0, t1] : A_i.t in [t0, t1]}.

- D_k is decreasing in t, so D_k(t0) >= 1 is necessary but weak.
- **For the order-<=8 alphabet:**
  - the sharp value is at most 0.2634 on all 40 cones tested with 2 <= t0 <= 4;
  - for t0 >= 2.5 even D_8(t0) < 1.
- So order 8 misses the counting condition by about a factor of 4. This agrees with P6 and with the session's own
  untracked `slope_preserving_subalphabet.out`.
- **Suggested wording.** "D_8(2) > 1: the weakest necessary condition now holds at the dense end of the cone. The
  sharp counting condition, restricted to images that stay in the cone, is <= 0.27 on every cone tested (audit S11)."

**m11. The Zaremba framing is presented as a finding.**
- Note section 0 item 5: "The orbit is thin, and the problem is a thin-orbit local-global problem of Zaremba type."
- The thinness is measured: X^1.85 against X^3. I reproduced the exponents 1.852 and 1.848.
- "Of Zaremba type" is an analogy: the semigroup has singular generators, and no expansion property is established.
- Retype as "resembles (analogy)".

**m12. Note section 5 P3: "decays at the 3D exponent (2 - 1.83)" is garbled.**
- Suggest "decays like H^(1.83 − 2), the rate implied by the 3D orbit's growth exponent".

## The proposed addendum (F) (two-endpoint parity)

**Verdict: CORRECT. Adopt it**, with the correction to its computer-check paragraph below.

**Step-by-step check of the proof.**
1. **Reduction.** R_z(T) = Σ_{U∋z} E(T[U], z) H(T − U). Every H(T − U) is odd (Redei, and H(∅) = 1). So
   R_z ≡ pi_z, and with (D), pi_z ≡ S_z. Correct.
2. **(F1).** Pi = Σ_v pi_v ≡ Σ_v S_v = H ≡ 1. Correct.
   - A one-line alternative: Pi(W) = Σ over nonempty U ⊆ W of H(W[U]), a sum of 2^n − 1 odd numbers.
3. **T = G + w, with w → u only.**
   - A Hamiltonian path of T from z either ends at w, preceded by a vertex ≠ u (giving s_z − P_(z->u)), or has w
     followed by u. In the second case the prefix starts at z and the suffix at u (giving L_(z,u)).
   - A path of T ending at z either avoids w (giving pi_z(G)), or is [A] w [C]. Here C is a u→z path of G, and A is
     any path of G − V(C), possibly empty; its last vertex is not u, so it beats w.
   - Both identities are correct.
4. **Parity of the sum.** 1 + Pi(W) is even unless W = ∅, so the sum is ≡ P_(u->z). Combining gives
   L_(z,u) ≡ P_(z->u) + P_(u->z). Correct.
5. **M from L.** M_(z,u)(G) = L_(z,u)(G^op), by reversing both paths, and the right-hand side is symmetric. Correct.
6. **(F3).** By (A) the S-row is (s_z, −P_(z->u), L_(z,u)). If it is zero, then M ≡ L = 0, and (E) gives an even
   plane determinant. If z = u, then M = 0. Correct. Two remarks worth adding:
   - Every seam word is itself a marked gadget, so (F3) covers words.
   - A gadget maps P into P iff s_z = L_(z,u) = 0, iff its S-row is zero (since P_(z->u) <= s_z). So (F3) covers
     every plane-preserving gadget.

**Brute force (own code, 0 violations everywhere).**
- **(F1) and (F2), exhaustive** (`aud_parity.c`): all labelled tournaments of order <= 7.
  - That is 89,084,724 ordered rows z ≠ u: 1,004,340 at orders 2-6 and 88,080,384 at order 7.
  - Plus 300 random tournaments at each of orders 8, 9 and 10 (65,400 rows). (F1) was also checked on random
    tournaments of orders 11 and 12.
- **The proof's two identities, exactly (not mod 2)** (`aud_F_proof_steps.py`), by explicit path enumeration:
  - (i) S_z(T) = s_z − P_(z->u) + L_(z,u);
  - (ii) pi_z(T) = pi_z(G) + Σ_C (1 + Pi(G − V(C)));
  - on every G of order 2-5 and every ordered pair z ≠ u: 21,300 rows. The two reductions mod 2 also hold.
- **(F3)** (`aud_census.out`): all 22,522 distinct S-row-zero gadget matrices of order <= 8 have M even, R_z even and
  an even 2x2 determinant.
- **The addendum's own counts reproduce.**
  - 1,004,340 rows at orders 2-6.
  - 130,364 rows with s_z = P_(z->u) = L = 0, none with M odd.
  - Order 7 adds 7,021,056 such rows, again none with M odd.

**Correction (MINOR): the relation search paragraph.**
- It reads "the exhaustive search over mod-2 linear relations among (...) found exactly the relations M ≡ L and
  M ≡ P_(z->u) + P_(u->z) and their consequences". That is wrong.
- The relation space has **dimension 3** at orders 2-6 and 2-7. There are no affine relations, and 105 (orders 2-6)
  and 128 (orders 2-7) parity patterns occur.
- The third independent relation is **Lmix ≡ Mmix (mod 2)**. It has a two-line proof by the same reduction:
  - (D) on the converse gives e(W, u) ≡ Σ over V ⊆ W with V ∋ u of s(V, u).
  - Substituting this into Lmix = Σ_U s(U, z) e(U^c, u) and into Mmix = Σ_U e(U, z) s(U^c, u) turns both into the
    sum over disjoint (W ∋ z, V ∋ u) of s(W, z) s(V, u).
- Replace the paragraph with "found a relation space of dimension 3, spanned by M ≡ L, M ≡ P_(z->u) + P_(u->z) and
  Lmix ≡ Mmix".

**Optional simplification.** The same substitution proves M ≡ L directly:

  M = Σ_U e(U,z) e(U^c,u) ≡ Σ over disjoint (W ∋ z, V ∋ u) of s(W,z) s(V,u) 2^(n−|W|−|V|) ≡ L,

since only covering pairs survive mod 2. The w-construction is then needed only for L ≡ P_(z->u) + P_(u->z).

**Novelty.**
- A targeted git grep of 01-canon, 05-knowledge and session19* found no earlier statement of
  L ≡ M ≡ P_(z->u) + P_(u->z).
- The reduction step is the mod-2 shadow of Session 193's identity (see M1), so credit Session29 and Session 193
  where (D) is used.

## Reproduced independently (own code; audited values in brackets)

| item | audit result |
|---|---|
| (A) row formulas, permutation brute force | 0 mismatches in 32,045 exhaustive instances (all labelled gadgets of order <= 4, every (u,z), every rooted prefix of order <= 3), 25,600 (all labelled order-5 gadgets, every (u,z), one random prefix) and 6,000 random (gadgets <= 5, prefixes <= 4). A mutation test (wrong K2 index) is detected |
| (A) on every gadget of order <= 7 | route I (prefix fit with (1,1,1), (1,0,2), (1,0,4), different from the audited prefixes) = route II ((A) statistics by subset DP) on all 25,689 (class, u, z) of order <= 7, plus 942 random order-8 triples |
| (B) and (C), TT_k | 285 cases, k <= 9 [204, k <= 8]: closed form = prefix fit = brute-force (A) statistics; a hold-out C3 prefix agrees. The rank rule holds. Both proofs re-derived and complete (slot counts for K1, K2, the 3^(z−u−1) factor, mu) |
| (D) | 14,891,097 rooted tournaments: all labelled ones of order <= 7 [<= 6] plus random ones of orders 8-12; 0 violations |
| (E) | all per-order S-row-zero matrices to order 8 (1,995 to order 7 [1,995], 20,787 at order 8): R_z, M and det2 even, det2 ≡ M |
| Census by exact order | 1, 4, 12, 44, 220, 1,660, 20,500, 427,412 [same]; classes 1, 1, 2, 4, 12, 56, 456 (A000568) |
| Cumulative / S-row-zero / nonsingular restrictions | 1, 5, 16, 56, 264, 1,884, 22,243, 448,915 / 0, 1, 5, 16, 57, 267, 1,919, 22,522 / 0, 0, 1, 5, 29, 191, 1,632, 21,023 [same] |
| `mats1-6.txt` | identical as sets to the audit's `aud_mats1-6.txt` |
| D_k(t), k = 4-8, t = 2, 2.5, 3, 4 | identical to 4 decimals; growth 2.20, 1.89, 1.74, 1.69 [same] |
| Seam BFS, order <= 6 | 240,925 / 4,743,006 / 17,081,577 states at X = 2,000 / 10,000 / 20,000, missing only 7 and 21 [same]; exponents 1.852 and 1.848 [1.85] |
| P2/P3 baseline (order <= 6, X = 5,000) | 1,301,611 states; fills 0.2640 and 0.2389; ratio 0.9051 [same] |
| P2/P3 (order <= 7) | 2,485,464 states, missing 7 and 21; fills 0.4573 and 0.4119; ratio 0.9008 [same] |
| Plane seeds | the 12 seeds are exactly the S = 0 states of order <= 5 (`aud_seeds.out`) |
| P6 | maps / states / fills / ratio: 191 / 20,149 / 0.0256, 0.0222 / 0.8691; 1,632 / 46,521 / 0.0615, 0.0442 / 0.7175; 21,023 / 124,006 / 0.1506, 0.1081 / 0.7174 [same] |
| Slope statistics (order <= 8) | 54,994 states, median 2.64, p90 4.15, histogram 36,235 / 12,363 / 5,221 / 1,150 / 25 [same]; slope-averaged D_8 0.8819 [0.8825] |
| Transitive-only words | 261 misses up to 2,000, first 7, 11, 19, 21, 35, 37 [same]; kept states 11,452 |
| Section 2 determinants | order 5: 107 of 220 singular, range −29 to 288 [same]; order 3: dets 0, 1 and 4 |
| Section 4 (b) | the same three filter hits; no word with a double eigenvalue 1 |
| Section 4 (c) | FALSE (M2) |
| Section 4 (d) | 3v confirmed, plus the omitted steps (3H, 0, 12R) and (3H, H, 3H) |

## Pre-registration integrity (git)

- **Order of commits.**
  - a98199e8bb (20:27:31): P1-P3, under the name HYP-9470.
  - 973c1f8312 (20:29:09): P1 read, P4 and P5 registered.
  - 7297af7bf9 (20:38:05): P2-P5 results; P6 registered together with `plane_bfs.c`.
  - 4b0b0c740b (20:43:30): P6 results.
- **The prediction texts are unchanged** in the current file:
  - P1-P3 as in a98199e8bb;
  - P4 and P5 as in 973c1f8312;
  - P6 as in 7297af7bf9.
- **The `seam_bfs.c` cap.** It was raised from 5,000 to 40,000 matrices in 973c1f8312, which the order-7 run needs.
  My run with all 22,243 distinct matrices reproduces P2/P3 exactly, so the reported run was not truncated.
- **Renumbering.**
  - opus-S79 reserved HYP-9470 in 395a317b0e at 20:04:24, 38 commits before a98199e8bb. That is consistent with
    "outside the 20-commit window".
  - The renumbering commits are 57cb7d51ee and aa81fd0f7b. They are recorded honestly in the HYP status and the note
    section 7.
- **Caveats (not violations).**
  - The order-7 matrix list existed at a98199e8bb (its count is quoted), so P1's blindness rests on the stated "Order
    7 was counted only".
  - The 20:29-20:38 window is consistent with the about-7-minute order-8 run starting after P4/P5 were registered,
    but git cannot prove it.

## Novelty and credit

- **(A).** KNOWN in substance (SEAM.md (1)-(3); opus-S78 THM-4875 (B) closed entries, converse convention).
  Correctly credited.
- **(B).** The single-port case is Session 195 `EQUAL_RESPONSE_BARRIER.md` section 2. I checked the mapping
  j = k−1−u, under the converse convention:

  | Session 195 | THM-4879 (u = z) |
  |---|---|
  | s = [j = 0] | E_u |
  | L = 2^j | L_u |
  | Q = 2^(m−1−j) | R_z |
  | B = 2^(m−1) | c |

  This is credited. For the two-port closed form with 3^(z−u−1) I found no earlier statement (grep for "3^(z-u",
  "two-port transitive" and the slot wording). It is an evaluation of S78's subset-sum entries, and new as a closed
  form.
- **(C) and (E).** New elementary consequences.
- **(D).** KNOWN (M1).
- **(F).** New as far as found.
- **Section 4 (d).** It should cite Session 195's equal-response theorem as the proved, scoped result of this kind.
- **origin/main.** origin/main equals HEAD (4b0b0c740b). No other commit since 2026-11-04 concerns seam transfers,
  this parity or two-port gadgets, besides opus-S78, the Session 193 and 195 packets, and opus-S79's reservation (its
  content is not committed).

## Programs and outputs (this folder)

**Programs.**
- `aud_rows_brute.c`: (A) by permutation brute force.
- `aud_tt.c`: (B) and (C), three routes.
- `aud_parity.c`: (D), (F1), (F2), and the GF(2) relation space.
- `aud_census.c`: classes, the census by two routes, and the (E), (F3) and P5 parities.
- `aud_density.py`: D_k(t).
- `aud_bfs.c`: the 3D and plane BFS.
- `aud_seeds.py`: the plane seeds.
- `aud_ttwords.py`: transitive-only words.
- `aud_structure.py` and `aud_structure_examples.py`: section 4 (b)-(c).
- `aud_affine_lines.py`: common invariant affine lines.
- `aud_ordinal.py`: section 4 (d).
- `aud_cone_condition.py`: m10, and the singular fractions for m7.
- `aud_slope_variants.py`: m9.
- `aud_F_proof_steps.py`: the (F) proof identities.

**Outputs.** The matching `.out` files, plus `aud_plane_seeds.txt` and `aud_mats1-6.txt` (identical to the audited
`mats1-6.txt`).

**Not in this folder.** The order-7 and order-8 matrix lists (617 KB and 14.6 MB) were written to the auditor's
scratchpad. Regenerate them with `gcc -O2 -o aud_census aud_census.c && ./aud_census 8 DIR`, which takes about 20 s.

**Run times.** Every run took under 2 minutes. The longest were the order-7 parity sweep (95 s) and the order-<=7
seam BFS (73 s).

**Not done.**
- The order-8 classes were not canonicalised. The census uses the 58,368 one-vertex extensions of the order-7
  classes, which carry every order-8 matrix. So the count 6,880 rests on A000568.
- Order 9 was not computed.
