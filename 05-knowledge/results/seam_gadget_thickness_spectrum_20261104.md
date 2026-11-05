# Composable seam pieces by order: thickness, thin orbits and what a self-covering proof of the {7, 21} spectrum needs

**Session:** kind-pasteur-2026-11-04-S11 (owner prompt: "look for recursive classifications of pieces of tournament
structure that can compose and grow with size, in order to come to the proof of the full {7, 21} spectrum").

**Canon and hypotheses:**
- [THM-4879](../../01-canon/theorems/THM-4879-two-port-transitive-seam-gadgets-closed-transfer-matrices-with-a-three-power-entry-and-the-seam-parity-r-equals-s-mod-2.md):
  two-port transitive gadgets, the seam parities, and even plane determinants at every order. PROVED, audited.
- [HYP-9471](../hypotheses/HYP-9471-seam-gadget-thickness-by-piece-order-plane-covering-density.md): pre-registered.
  P1, P2, P4 and P5 were CONFIRMED; P3 and P6 were REFUTED.

**Audit:** independent audit S11 PASS WITH CORRECTIONS (4 major, 12 minor, all applied; MISTAKE-717),
`04-computation/experiments/seam_gadget_thickness_20261104_audit_S11/AUDIT.md`.

**Programs:** `04-computation/experiments/seam_gadget_density_20261104/`.

**The {7, 21} spectrum stays OPEN.** No new value is claimed attainable and none is excluded.

## 0. Summary

The composable pieces are SEAM.md's marked gadgets (G, u, z): a tournament with an entry u and an exit z. They are
glued by one reversed cross arc. Every seam word is a product of 3x3 integer transfer matrices acting on
v = (H, S, R), so the attainable H-values of seam words are a matrix-semigroup orbit. This session classifies the
pieces by their matrices and measures how the orbit thickens as the pieces grow.

**Positive findings.**
1. **Census.** The number of distinct transfer matrices of gadgets of order exactly k is 1, 4, 12, 44, 220, 1,660,
   20,500 and 427,412 for k = 1-8. All 6,880 order-8 tournament classes and all 64 port pairs were enumerated.
2. **An explicit all-order family (THM-4879 (B)).** For two-port transitive gadgets TT_k, the R-to-R coefficient is
   2^(k-z+u) 3^(z-u-1) when the entry precedes the exit: each vertex in between has three slots.
   - Seam words of transitive pieces alone are thin. They miss 261 of the 1,000 odd numbers <= 2,000, among them 11
     and 19.
3. **Parities.**
   - R ≡ S (mod 2) for every rooted tournament (KNOWN: Session29; re-proved in THM-4879 (D)). On the plane
     P = {S = 0}, R is even.
   - New, THM-4879 (F): the counts of ordered two-path covers starting at z and u, and of those ending at z and u,
     are both ≡ P(z->u) + P(u->z) (mod 2). Also, every tournament has an odd number of directed paths.
   - Consequence: every gadget, and every seam word, that maps P into P acts there with even determinant, at every
     order. No unimodular map exists on P.
4. **Coverage.** Seam words over all gadgets of order <= 6 (1,884 distinct matrices) attain every odd H <= 20,000
   except 7 and 21. This is consistent with Session 195's grammar and with opus-S78's F1 + F2 to 2^32, which are
   different families.

**The structural finding.**
5. **The orbit is thin.**
   - The order-<=6 orbit has about X^1.85 distinct states with H, R <= X, against about X^3 lattice points in that
     box.
   - Its H-projection is nonetheless full in the tested range. This *resembles* a thin-orbit local-global problem of
     Zaremba type. That is an analogy: the semigroup has singular generators, and no expansion property is
     established.
   - No word of length <= 3 over the gadgets of order <= 4 is quasi-unipotent, or even has a double eigenvalue 1. So
     iterating one such word never gives an arithmetic progression of H-values.
   - No two gadgets of order <= 5 share an invariant affine line off the origin. Shared expanding eigen-rays exist,
     but on them the maps are pure scalings.
   - A closed {3N, 3N + 2, 3N + 4} grammar built on the scalar step X => C3 is impossible (section 4 (d)). The other
     ×3 steps of the triple calculus are not covered.

**The quantitative picture (HYP-9471).**
6. On the plane P, the density D_k(2) = sum (a + 2b)^(-2) over the plane maps is 0.117, 0.257, 0.486, 0.848 and
   1.429 for k = 4-8. It passes 1 at order 8, as predicted.
   - That is only the weakest necessary condition, at the dense end of the cone.
   - The sharp counting condition for a self-covering cone is at most 0.27 at order 8 (audit).
   - The realised plane orbit still thins with scale at order 8: averaged over its own slopes the density is 0.882
     (P6 refuted).
   - So no gadget alphabet of order <= 8 supports a self-covering certificate. Even one that did would also need
     exact coset alignment.

## 1. Inheritance and concurrent work

- **Transfer.** SEAM.md (`seven_return_native_recursions_20261029/SEAM.md`): the universal three-coordinate
  transfer, composition M_(X⊙Y) = M_Y M_X, the adjoint, and the singleton matrix with its tribonacci polynomial.
- **Session29** (`session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md` section 5): R ≡ S
  (mod 2), in this convention.
- **opus-S78 THM-4875.**
  - Closed entries in the converse convention.
  - The descent lemma: unimodular, strictly h-raising cone coverings are impossible. Its listed escapes are h-fixed
    moves, determinants |d| > 1, and information beyond the triple.
  - F1 + F2 attain every odd N <= 2^32 except 7 and 21.
- **Session 193.**
  - The continuation-rank obstruction for nested C3 substitutions.
  - The response rigidity of two-apex macros, and modular synchronisation.
  - The rooted-path parity identity (multiple_triangles section 2), whose mod-2 shadow is R ≡ S.
- **Session 194.** The module barrier, the endpoint floor, and fixed-width windows.
- **Session 195.** The four-coordinate grammar (h, e, R, d) with stationary-root moves; single-port transitive
  formulas; and order-<=5 gadgets plus two order-6 pieces covering to 20,000 except {7, 21}.
  - Its equal-response theorem (`session195/rooted-interfaces/EQUAL_RESPONSE_BARRIER.md`) is the proved, sharply
    scoped barrier of the kind section 4 (d) gestures at: moving-root descendants of a state with R = h have counts
    <= 3h + 4 exactly {h, 2h - e, 3h}.
  - This session found Session 195 mid-run and then did not repeat coverage runs. Its own contribution is the
    thickness analysis and the parities.
- **opus-S79** reserved THM-4878 and HYP-9470 for a "bounded-cutwidth spine grammar", concurrent and not read.
  - This session's hypothesis was first committed as HYP-9470, then renumbered HYP-9471 (see section 7).
- **Hostile control.** Transitive pieces alone (THM-4879 (B)) are an infinite closed-form family that misses 11.

## 2. The classification of pieces by transfer matrix

`gadget_matrices.c` generates every tournament class of order N by vertex extension and canonical form. The class
counts 1, 1, 2, 4, 12, 56, 456, 6,880 agree with A000568. For every class and every port pair it computes the SEAM
matrix from the spanning prefixes (1,1,1), (1,0,2) and (3,1,3).

| order k | classes | distinct matrices (exactly k) | cumulative distinct | S-row zero (cum.) | nonsingular 2x2 on P (cum.) |
|---|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 0 | 0 |
| 2 | 1 | 4 | 5 | 1 | 0 |
| 3 | 2 | 12 | 16 | 5 | 1 |
| 4 | 4 | 44 | 56 | 16 | 5 |
| 5 | 12 | 220 | 264 | 57 | 29 |
| 6 | 56 | 1,660 | 1,884 | 267 | 191 |
| 7 | 456 | 20,500 | 22,243 | 1,919 | 1,632 |
| 8 | 6,880 | 427,412 | 448,915 | 22,522 | 21,023 |

Cross-checks:
- Orders <= 6 agree with an independent labelled enumeration (every labelled tournament, deduplicated by matrix), and
  with the audit's own enumeration.
- Orders <= 5 agree with SEAM's Python `matrix()`.
- The audit reproduced the whole table.

**Order 3.** The 12 matrices are the 9 port pairs of TT_3 and the 3 port classes of the directed triangle C3.
- One C3 port pair is unimodular: [3 -1 3; 1 0 1; 3 -1 4], the cube of the singleton matrix.
- C3 with u = z has determinant 4.

**Singular matrices.** Many small gadget matrices are singular: they are the "reset" pieces whose output depends on
two or fewer coordinates of the input. The singular fractions at orders 2-8 are 75.0, 75.0, 63.6, 48.6, 32.7, 19.8
and 11.3% (audit). At order 5, determinants range from -29 to 288.

## 3. THM-4879 in one paragraph

**Transitive gadgets.**
- For G = TT_k with entry u and exit z, the H-row is (1, -[u = k-1], 2^(k-1-u)). The S-row is explicit.
- The R-row is (2^z, -mu, c), where c = 2^(k-1-u+z) for z <= u and c = 2^(k-z+u) 3^(z-u-1) for z > u. In the second
  case each vertex strictly between u and z has three admissible slots in the two-cover counts K1 and K2.
- (u, z) = (0, k-1) is the only unimodular port pair; all other port pairs with k >= 2 are singular.

**Parities.**
- R ≡ S (mod 2) for every rooted tournament (KNOWN). The singleton-seam proof: H(X ⊙ {z}) = H - S + R, and both
  counts are odd.
- (F): L_(z,u) ≡ M_(z,u) ≡ P(z->u) + P(u->z) (mod 2), proved with one auxiliary vertex w that beats only u.
- So a gadget maps P into P iff its S-row is zero, and then M_(z,u)(G) is even. Its plane determinant is
  ≡ M_(z,u)(G), hence even, at every order.

## 4. Why the orbit is thin, and what the small-gadget structure does and does not exclude

**(a) Growth.** Seam BFS with the order-<=6 matrices (1,941 lines, 1,884 distinct), seeds (1,1,1), (1,0,2), (3,1,3)
and (1,0,2^j):

| X (states with H, R <= X expanded) | distinct states | odd H <= X missed |
|---|---|---|
| 2,000 | 240,925 | 7, 21 |
| 10,000 | 4,743,006 | 7, 21 |
| 20,000 | 17,081,577 | 7, 21 |

- The state count grows like X^1.85 (the audit's exponents are 1.852 and 1.848), against about X^3 lattice points in
  the box. The orbit is a thin subset of Z^3 whose H-projection is nonetheless full in the tested range.
- This resembles thin-orbit local-global problems such as Zaremba's conjecture. There, analytic methods of
  Bourgain-Kontorovich type give "almost all" targets, not "all". The resemblance is an analogy, not a theorem about
  this semigroup.

**(b) No quasi-unipotent short words.** Words of length <= 3 over the 56 distinct gadget matrices of order <= 4 were
searched.
- Only two characteristic polynomials with |det| = 1 and |tr|, |c2| <= 3 occur:
  - t^3 - t^2 - t - 1, the singleton (tribonacci);
  - t^3 - 3t^2 - t - 1, its square, which is the matrix of TT_2 entered at its source and left at its sink.
- No such word has (t - 1)^2 dividing its characteristic polynomial, whatever its determinant.
- So iterating a single short word never gives an arithmetic progression n -> H(W^n q).
- This says nothing about progressions inside the whole orbit's H-values. The coverage data contain 23, 25, ...,
  19,999, and an infinite progression there is exactly what the spectrum conjecture asserts for seam words.

**(c) Shared eigen-directions** (corrected after audit S11). The search is over the 264 distinct gadget matrices of
order <= 5 (281 counting cross-order repeats).
- 10 right and 10 left rational eigen-directions are shared by two or more matrices with nonzero eigenvalue.
- 8 of them are expanding. For example, the singleton ray (1, 1, 1) is fixed by 12 reset gadgets with eigenvalues 1,
  3, 5, 9, 11, 13 and 15. On such a ray the maps act by pure scalings.
- No two distinct gadgets of order <= 5 share an invariant affine line off the origin (audit, `aud_affine_lines.py`).
  So single small gadgets give no 1D affine digit grammar.
- The first version of `small_gadget_structure.py` tried only divisors of det as eigenvalues, so it missed every
  nonzero eigenvalue of a singular matrix. It reported a single shared direction. The corrected script reproduces the
  audit's 10/10/8.

**(d) The scalar ×3 step admits no affine digits** (elementary).
- The ordinal sum X => C3 with the exit kept in X sends v(X, a) to 3 v(X, a), because H, S and R all triple.
- Let D be a set of states closed under that step and under linear transfers M_2, M_4 with H(M_d q) = 3H(q) + d on
  D. Then for q in D,
  3(3H(q) + d) = H(M_d 3q) = 3H(3q) + d = 9H(q) + d,
  so d = 0.
- **This covers only the scalar step.**
  - The other ordinal sum, C3 => X, stays inside the triple calculus but is not scalar: with the exit kept in X it
    sends v to (3H, 0, 12R), and with the exit in C3 it sends v to (3H, H, 3H) (audit, brute force).
  - Neither is covered, nor are digit maps whose law holds only on part of D, nor state-dependent words.
  - X => C3 with the exit in C3 leaves the triple calculus: its state is (3H, 0, 3 C2(X)), where C2(X) = 2H + 2pc_2(X)
    and pc_2 counts unordered two-path covers. That is R = 6(H + pc_2).
  - Session 195's equal-response theorem is the proved, scoped result of this kind.

## 5. The covering density and the order-8 picture (HYP-9471)

**Definition and heuristic.**
- On P every S-row-zero gadget acts as (H, R) -> (aH + bR, cH + dR), and D_k(t) = sum over the distinct nonsingular
  such maps of order <= k of (a + bt)^(-2). For a map A at slope t = R/H, the factor a + bt is its local growth of H.
- Pre-registered heuristic: a point at scale H is covered with expected multiplicity about D·f, where f is the fill
  at smaller scales. The self-similar fill satisfies f = 1 - exp(-D f), which has a positive solution only when D > 1.

**D_k(t):**

| k | D_k(2) | D_k(2.5) | D_k(3) | D_k(4) |
|---|---|---|---|---|
| 4 | 0.1170 | 0.0812 | 0.0596 | 0.0361 |
| 5 | 0.2571 | 0.1787 | 0.1315 | 0.0798 |
| 6 | 0.4861 | 0.3388 | 0.2498 | 0.1521 |
| 7 | 0.8479 | 0.5923 | 0.4374 | 0.2668 |
| 8 | 1.4288 | 0.9999 | 0.7395 | 0.4519 |

**The sharp condition** (audit S11). D_k(t0) >= 1 at the dense end is necessary but weak. Counting lattice points
with H <= X shows that a finite alphabet can self-cover the cone t0 <= R/H <= t1 only if

    (1/(t1 - t0)) Σ_i ∫ (a_i + b_i t)^(-2) dt >= 1,

where the integral runs over the slopes t in [t0, t1] whose image slope also lies in [t0, t1]. For the order-<=8
alphabet this is at most 0.2634 on every one of 40 cones tested with 2 <= t0 <= 4. That is about a factor of 4
short.

**Pre-registered outcomes.**
- P1 (D_7) CONFIRMED.
- P2 CONFIRMED: the order-7 3D-fed plane fill on H in [801, 999], R in [2H, 2.75H] is 0.4573, against 0.2640 at
  order 6.
- P3 REFUTED: the decay ratio stayed at 0.9008 against 0.9051. A slice fed by the whole 3D orbit through rank-2 maps
  decays like H^(1.83 - 2), the rate implied by the 3D orbit's growth exponent. It does not decay at a rate set by
  D_k.
- P4 (D_8) CONFIRMED.
- P5 (no odd determinant to order 8) CONFIRMED. It is now a theorem at all orders (THM-4879 (F)).
- P6 REFUTED. In the pure plane dynamics the fills are 0.0256, 0.0615, 0.1506 for orders <= 6, 7, 8, and the scale
  ratios are 0.8691, 0.7175, 0.7174.

**The reason, after the fact.**
- The realised order-8 plane states with H in [801, 1799] sit at slopes 2-3: median 2.64, 90th percentile 4.15.
  The maps push states toward larger slopes, where D is smaller.
- Averaged over the realised slopes, D_8 = 0.8819 < 1 (the audit's exact mean; the first, subsampled value was
  0.8825). The sharp condition above is far lower still.

**Post-audit exploratory addendum: a slope-preserving sub-alphabet does not help at order 8.** This was added after
the audit, is not audited, and is not pre-registered. Program: `slope_preserving_subalphabet.py`; output:
`slope_preserving_subalphabet.out`.
- Of the 21,023 nonsingular plane maps of order <= 8, 8,815 map the slope band [2, 3] into itself.
- Their density on the band is low: 0.2598 at slope 2, 0.1821 at 2.5, and 0.1348 at 3.
- The pure plane BFS with only these maps thins faster than with the whole alphabet: fills 0.1416 and 0.0928, ratio
  0.6557.
- The maps with small growth a + bt, which carry most of D_8, are exactly those that send slopes upward out of the
  dense band. This is consistent with the audit's sharp condition.

## 6. What a proof along this route needs, and what is ruled out

**(i) An alphabet that satisfies the sharp counting condition on some cone.** At order 8 it is about 4× short. On
the observed growth (D_k roughly ×1.7 per order) that points to order 10 or more, or to pieces outside the plane
calculus.

**(ii) Exact coset alignment.** Positive density is not all points.
- Every plane map, of any order, has even determinant (THM-4879 (E) with (F)). So each image covers at most half of
  the lattice cosets it meets.
- A certificate must combine maps with complementary cosets in every direction sector of a subcone, and cover its
  boundary rays. S78's lemma forbids doing this with unimodular maps, and none exist on P.
- Sectors and residues are finitely many once the alphabet is fixed. A computer-assisted search for such a
  certificate is the concrete next step if (i) is reached.

**(iii) Otherwise an analytic local-global theorem.** Such a theorem would at best give "almost all" odd N, not all.

**Excluded within the triple calculus:**
- affine digit grammars on a common invariant line of two single gadgets of order <= 5 (section 4 (c), corrected);
- ternary grammars built on the scalar step X => C3 (section 4 (d));
- arithmetic progressions from iterating one short word (section 4 (b));
- unimodular cone coverings (S78), and on P any unimodular map at all (THM-4879 (F));
- transitive pieces alone (THM-4879 and the 261 misses).

## 7. Typing, provenance, reproduction

**Typing.**
- PROVED: THM-4879 (B), (C), (E), (F) (elementary). (A) is KNOWN in substance; (D) is KNOWN (Session29, Session 193)
  with a shorter alternative proof.
- FINITE-EXACT:
  - the census;
  - the coverage to 20,000;
  - the corrected eigen-direction search and the quasi-unipotent search;
  - P5 (now subsumed by (F)).
- Pre-registered tests: HYP-9471.
- Heuristic: the covering-density model of section 5. It was used only to make the predictions, two of which failed.
- Audit: PASS WITH CORRECTIONS, MISTAKE-717.

**Numbering.**
- The hypothesis was committed as HYP-9470 in a98199e8bb and 973c1f8312.
- opus-S79 had reserved HYP-9470 at 20:04 (395a317b0e), outside the 20-commit window this session checked before
  reserving. It was renumbered HYP-9471 in 57cb7d51ee and aa81fd0f7b; the predictions were not changed.

**Reproduction** (in `seam_gadget_density_20261104/`):
- **Gadget matrices:**
  `gcc -O2 -o gadget_matrices gadget_matrices.c; for N in 1..8: ./gadget_matrices N | sort -u > matsN.txt`.
  Orders 1-6 are committed; order 8 takes about 7 min.
- **Densities:** `python3 density.py 8 mats1.txt ... mats8.txt` gives `density_orders_1_to_8.out`.
- **Coverage, P2 and P3:** `gcc -O2 -o seam_bfs seam_bfs.c; ./seam_bfs X allmats.txt` (+ `fill.py pstates.txt`)
  gives `seam_bfs_coverage_order6.out` and `seam_bfs_p2_p3.out`.
- **P6:** `plane_bfs.c` gives `plane_bfs_p6.out`.
- **THM-4879 checks:** `transitive_two_port.py 8` and `parity_check.py`.
- **Transitive-only words:** `transitive_words_bfs.py 2000 10`.
- **Section 4 (b)-(c), corrected:** `small_gadget_structure.py`.
- **Post-audit addendum:** `slope_preserving_subalphabet.py`.
- **Audit programs:** `../seam_gadget_thickness_20261104_audit_S11/aud_*`, including the sharp cone condition, the
  exact slope mean and the ordinal-sum checks.
