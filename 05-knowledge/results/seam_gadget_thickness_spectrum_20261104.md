# Composable seam pieces by order: thickness, thin orbits and what a self-covering proof of the {7, 21} spectrum needs

**Session:** kind-pasteur-2026-11-04-S11 (owner prompt: "look for recursive classifications of pieces of tournament
structure that can compose and grow with size, in order to come to the proof of the full {7, 21} spectrum").

**Canon and hypotheses:**
- [THM-4879](../../01-canon/theorems/THM-4879-two-port-transitive-seam-gadgets-closed-transfer-matrices-with-a-three-power-entry-and-the-seam-parity-r-equals-s-mod-2.md):
  two-port transitive gadgets and the seam parity. PROVED, audit pending.
- [HYP-9471](../hypotheses/HYP-9471-seam-gadget-thickness-by-piece-order-plane-covering-density.md): pre-registered.
  P1, P2, P4 and P5 were CONFIRMED; P3 and P6 were REFUTED.

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
3. **Seam parity (THM-4879 (D)).** R ≡ S (mod 2) for every rooted tournament, by Redei and the singleton seam. On the
   plane P = {S = 0}, R is even.
   - Every plane map of order <= 8 has even determinant: 22,522 distinct matrices were checked, and the reduction
     det ≡ M_(z,u)(G) (mod 2) is proved.
   - So no unimodular map exists on P.
4. **Coverage.** Seam words over all gadgets of order <= 6 (1,941 matrices) attain every odd H <= 20,000 except 7 and
   21. This is consistent with Session 195's grammar and with opus-S78's F1 + F2 to 2^32, which are different
   families.

**The structural finding.**
5. **The orbit is thin, and the problem is a thin-orbit local-global problem of Zaremba type.**
   - The order-<=6 orbit has about X^1.85 distinct states with H, R <= X, against about X^3 lattice points in that
     box.
   - No short gadget word is unipotent, so no orbit contains an arithmetic progression.
   - No two small gadgets share an expanding rational eigen-direction, so there is no invariant line for a 1D
     digit grammar.
   - A closed {3N, 3N + 2, 3N + 4} grammar is impossible in the triple calculus (section 4).

**The quantitative threshold (HYP-9471).**
6. On the plane P, the 2D covering density D_k(2) = sum (a + 2b)^(-2) over the plane maps is 0.117, 0.257, 0.486,
   0.848 and 1.429 for k = 4-8.
   - It crosses 1, the necessary condition for a finite alphabet to fill a plane of states with positive density, at
     order 8. This was predicted and confirmed.
   - But the realised plane orbit still thins with scale at order 8: averaged over its own slopes the density is 0.88
     (P6 refuted).
   - So no gadget alphabet of order <= 8 supports a self-covering certificate. Even one that did would also need exact
     coset alignment.

## 1. Inheritance and concurrent work

- **Transfer.** SEAM.md (`seven_return_native_recursions_20261029/SEAM.md`): the universal three-coordinate
  transfer, composition M_(X⊙Y) = M_Y M_X, the adjoint, and the singleton matrix with its tribonacci polynomial.
- **opus-S78 THM-4875.**
  - Closed entries in the converse convention.
  - The descent lemma: unimodular, strictly h-raising cone coverings are impossible. Its listed escapes are h-fixed
    moves, determinants |d| > 1, and information beyond the triple.
  - F1 + F2 attain every odd N <= 2^32 except 7 and 21.
- **Session 193.** The continuation-rank obstruction for nested C3 substitutions, the response rigidity of
  two-apex macros, and modular synchronisation.
- **Session 194.** The module barrier, the endpoint floor, and fixed-width windows.
- **Session 195.** The four-coordinate grammar (h, e, R, d) with stationary-root moves; single-port transitive
  formulas; the equal-response barrier; and order-<=5 gadgets plus two order-6 pieces covering to 20,000 except
  {7, 21}.
  - This session found Session 195 mid-run and then did not repeat coverage runs. Its own contribution is the
    thickness and threshold analysis.
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
- Orders <= 6 agree with an independent labelled enumeration (every labelled tournament, deduplicated by matrix).
- Orders <= 5 agree with SEAM's Python `matrix()`.
- At order 3 the 12 matrices include the singleton, the transitive gadgets and the directed triangle. Its
  (u, z) = (0, 0) gadget has determinant 4, the smallest non-unimodular nonsingular case.

Most gadget matrices are singular: 107 of the 220 at order 5 have det 0. They are the "reset" pieces whose output
depends on two or fewer coordinates of the input. Determinants spread widely; at order 5 alone they range from
-29 to 288.

## 3. THM-4879 in one paragraph

**Transitive gadgets.**
- For G = TT_k with entry u and exit z, the H-row is (1, -[u = k-1], 2^(k-1-u)). The S-row is explicit.
- The R-row is (2^z, -mu, c), where c = 2^(k-1-u+z) for z <= u and c = 2^(k-z+u) 3^(z-u-1) for z > u. In the second
  case each vertex strictly between u and z has three admissible slots in the two-cover counts K1 and K2.
- (u, z) = (0, k-1) is the only unimodular port pair; all other port pairs with k >= 2 are singular.

**Parity.**
- R ≡ S (mod 2) for every rooted tournament. Proof: H(X ⊙ {z}) = H - S + R, and both counts are odd.
- So P-states have H odd and R even.
- Every gadget mapping into P has R_z(G) even, and plane determinant ≡ M_(z,u)(G) (mod 2), where M_(z,u)(G) counts
  ordered two-covers ending at z and at u.
- M_(z,u)(G) is even for every such gadget of order <= 8 (FINITE-EXACT). A proof for all orders is OPEN.

## 4. Why the orbit is thin, and the barriers in the triple calculus

**(a) Growth.** Seam BFS with all 1,941 order-<=6 matrices, seeds (1,1,1), (1,0,2), (3,1,3) and (1,0,2^j):

| X (states with H, R <= X expanded) | distinct states | odd H <= X missed |
|---|---|---|
| 2,000 | 240,925 | 7, 21 |
| 10,000 | 4,743,006 | 7, 21 |
| 20,000 | 17,081,577 | 7, 21 |

- The state count grows like X^1.85, against about X^3 lattice points in the box. The orbit is a thin subset of Z^3
  whose H-projection is nonetheless full in the tested range. That is the regime of thin-orbit local-global problems
  such as Zaremba's conjecture.
- There, analytic methods of Bourgain-Kontorovich type give "almost all" targets, not "all".

**(b) No parabolic elements.** Words of length <= 3 over the 61 gadget matrices of order <= 4 were searched for
characteristic polynomials with |det| = 1 and |tr|, |c2| <= 3. Only two polynomials occur:
- t^3 - t^2 - t - 1, the singleton (tribonacci);
- t^3 - 3t^2 - t - 1, its square, which is the matrix of TT_2 entered at its source and left at its sink.

Neither is unipotent. No orbit in this range therefore contains an arithmetic progression, which would have given
an easy proof.

**(c) No invariant lines.** Among the 281 gadget matrices of order <= 5, the only rational eigen-direction shared by
two of them with nonzero eigenvalue is (2, 1, -1) on the right and (2, -1, -1) on the left, with eigenvalues ±1.
So no line of states is preserved by two expanding gadgets. A 1D affine digit grammar h -> λh + d_i on an
invariant line cannot be built from single small gadgets. This is consistent with Session 193's response rigidity.

**(d) A closed ternary numeration is impossible in the triple calculus** (elementary).
- In the triple calculus the step N -> 3N is the ordinal sum X => C3 with the exit kept in X. It sends
  v(X, a) to 3 v(X, a), because H, S and R all triple.
- Let D be a set of states closed under that step and under linear transfers M_2, M_4 with
  H(M_d q) = 3H(q) + d on D. Then for q in D,
  3(3H(q) + d) = H(M_d 3q) = 3H(3q) + d = 9H(q) + d,
  so d = 0.
- The other ordinal sum, with the exit moved into C3, needs H(X) + 2pc_2(X) and leaves the triple calculus. Seam
  words cannot supply a pure ×3 step either, because their H-row is (H(G*), -S(G*), R(G*)) with R(G*) >= 1.
- This is the triple-calculus form of the obstruction that Sessions 193-195 meet as response rigidity, the module
  barrier and the equal-response barrier.

## 5. The covering density and the order-8 threshold (HYP-9471)

**Definition and heuristic.**
- On P every S-row-zero gadget acts as (H, R) -> (aH + bR, cH + dR), and D_k(t) = sum over the distinct nonsingular
  such maps of order <= k of (a + bt)^(-2). For a map A at slope t = R/H, the factor a + bt is its local growth of H.
- In 2D, a point at scale H is covered with expected multiplicity about D·f, where f is the fill at smaller scales.
  The self-similar fill satisfies f = 1 - exp(-D f), which has a positive solution only when D > 1.
- So D = 1 is the threshold where a finite alphabet starts to fill a plane of states with positive density.

**D_k(t):**

| k | D_k(2) | D_k(2.5) | D_k(3) | D_k(4) |
|---|---|---|---|---|
| 4 | 0.1170 | 0.0812 | 0.0596 | 0.0361 |
| 5 | 0.2571 | 0.1787 | 0.1315 | 0.0798 |
| 6 | 0.4861 | 0.3388 | 0.2498 | 0.1521 |
| 7 | 0.8479 | 0.5923 | 0.4374 | 0.2668 |
| 8 | 1.4288 | 0.9999 | 0.7395 | 0.4519 |

**Pre-registered outcomes.**
- P1 (D_7) CONFIRMED.
- P2 CONFIRMED: the order-7 3D-fed plane fill on H in [801, 999], R in [2H, 2.75H] is 0.4573, against 0.2640 at
  order 6.
- P3 REFUTED: the decay ratio stayed at 0.9008 against 0.9051. A slice fed by the whole 3D orbit through rank-2 maps
  decays at the 3D exponent (2 - 1.83), not at a rate set by D_k.
- P4 (D_8) CONFIRMED.
- P5 (no odd determinant to order 8) CONFIRMED.
- P6 REFUTED. In the pure plane dynamics the fills are 0.0256, 0.0615, 0.1506 for orders <= 6, 7, 8, and the scale
  ratios are 0.8691, 0.7175, 0.7174.

**The reason, after the fact.**
- The realised order-8 plane states with H in [801, 1799] sit at slopes 2-3: median 2.64, 90th percentile 4.15.
  The maps push states toward larger slopes, where D is smaller.
- Averaged over the realised slopes, D_8 = 0.8825 < 1, so the orbit is still subcritical.
- At the observed per-order ratio of about 1.6-1.7, the slope-averaged density would pass 1 at order 9. Order 9
  (191,536 classes, 81 port pairs) was not computed.

## 6. What a proof along this route needs, and what is ruled out

**(i) An alphabet whose realised plane orbit is supercritical.** On present evidence that needs order 9 or a
selected sub-alphabet that keeps slopes near 2. Order 8 is not enough.

**(ii) Exact coset alignment.** Positive density is not all points.
- Every plane map of order <= 8 has even determinant (THM-4879 (E) with P5), so each image covers at most half of
  the lattice cosets it meets.
- A certificate must combine maps with complementary cosets in every direction sector of a subcone, and cover its
  boundary rays. S78's lemma forbids doing this with unimodular maps, and none exist here anyway.
- Sectors and residues are finitely many once the alphabet is fixed. A computer-assisted search for such a
  certificate is the concrete next step if (i) is reached.

**(iii) Otherwise an analytic local-global theorem.** Such a theorem would at best give "almost all" odd N, not all.

**Ruled out within the triple calculus:**
- closed affine digit grammars, by section 4 (c)-(d);
- arithmetic progressions, by section 4 (b), at least from short words;
- unimodular cone coverings (S78);
- transitive pieces alone (THM-4879 and the 261 misses).

## 7. Typing, provenance, reproduction

**Typing.**
- PROVED: THM-4879 (A)-(E) (elementary; (A) KNOWN in substance).
- FINITE-EXACT: the census, the coverage to 20,000, the eigen-direction and unipotent searches, and P5.
- Pre-registered tests: HYP-9471.
- Heuristic: the covering-density model of section 5. It was used only to make the predictions, two of which failed.

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
- **Section 4 (b)-(c):** `small_gadget_structure.py`.
