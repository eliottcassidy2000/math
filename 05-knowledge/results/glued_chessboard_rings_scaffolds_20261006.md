# Pac-man chessboards: the four rings are the bishop's mobility shells, every gluing lives in the outer ring, one parity bit decides whether the two diagonal scaffolds fuse, and sliders only see which lines were joined

2026-10-06, session opus-2026-10-06-S15 (worktree `codex/session-glued-chessboard-20261006`). Owner's seed (verbatim): "Consider an 8 by 8 grid of 64 squares. Their horizontal and vertical connections form a set of concentric rings of size 4, 12,20,28 and their diagonal connections form 2 interwoven but isolated scaffoldings of {1,3,5,5,3,1} by {2,4,6,8,6,4,2} and consider the ability to move along these horizontal/vertical/diagonals one step at a time or unlimited amounts until an obstruction is reached, also consider the outcomes of connecting various edges of the 8 by 8 grid together in various mathematically interesting ways so that they loop together and allow pac man style teleportation"

**PROVED (elementary):** Theorems 1-6 and the mechanisms marked so. **FINITE-EXACT:** every table (exact enumeration on 16 boards; independence counts by two independent algorithms, isometry groups by two independent algorithms). **CITED:** Polya 1918 and Monsky 1989 (toroidal queens, OEIS A085801), Burger-Mynhardt 2003 / Mynhardt 2003 (toroidal queen domination, OEIS A279402), G. Arizmendi Echegaray, "Queens on surfaces" (Bridges 2026 art exhibition: 8x8 torus 6, Mobius band 4, Klein bottle 3 queens). **NUMEROLOGY guard:** no small-number coincidence below is typed as a dictionary without a mechanism. Nothing here touches LRC(14).

Script and saved output: [glued_chessboard_20261006.py](../../04-computation/experiments/glued_chessboard_20261006.py) / [out](glued_chessboard_20261006.out) / [json](glued_chessboard_20261006.json); the SAT-heavy domination and colouring runs are in [glued_chessboard_sat_20261006.json](glued_chessboard_sat_20261006.json).

## 1. Inheritance and portfolio

- **Closest proved mechanisms.** Gauss-Bonnet for flat cone surfaces; wallpaper groups as square-tiled orbifolds (Conway notation o, xx, 22x, 442, 2222); the 45-degree isometry between the max norm and the taxicab norm; Polya's parity argument for toroidal queens.
- **Canonical hostiles.** The diagonal-glide Klein bottle `abba` and the Klein bottle whose flip axis runs through cell centres (non-orientable, yet the two diagonal scaffolds stay separate); the 442 sphere (orientable, yet they fuse); the helical torus (a torus, yet they fuse).
- **Corrected near miss.** The seed's `{1,3,5,5,3,1}` is the inner 6x6 board; see section 2.
- **Least-used sidecars.** The line multigraph of a slider (section 8) and the colour cover (section 7).
- **Observer lens.** The owner's rings are what an observer at the board's centre point sees; every result below is phrased by where that observer's view wraps (ring 4) or where the board's curvature sits (cone points).

Portfolio. **Anchor:** make the rings and scaffolds exact and catalogue what every geometric gluing does to them and to each piece. **Niche:** the colour cover as the home of the bishop (section 7). **Wildcard:** which gluings a slider cannot tell apart (section 8).

## 2. The seed, made exact (and one correction)

Cells `(x,y)`, `0 <= x,y <= 7`. Ring `k` (k = 1..4) is the set of cells at Chebyshev distance `k - 1/2` from the centre point `(4,4)`; it has `8(k - 1/2) = 8k - 4` cells: **4, 12, 20, 28**. The orthogonal-step (wazir) graph is exactly four disjoint cycles `C4, C12, C20, C28` plus `8k` spokes between ring `k` and ring `k+1` (8, 16, 24; total 48 + 64 = 112 edges).

**Correction.** On the 8x8 board each colour class has diagonals of lengths **{1,3,5,7,7,5,3,1}** in one direction and **{2,4,6,8,6,4,2}** in the other (the two colours swap roles). The seed's `{1,3,5,5,3,1}` (with `{2,4,6,4,2}`) is the scaffold of the 6x6 board, i.e. of rings 1-3. In 45-degree coordinates `u=(x+y)/2, v=(x-y)/2` each colour is a diamond-shaped board whose rows and columns are exactly those diagonals; diagonal steps become orthogonal steps there. The two scaffolds are planar duals: every interior cell of one colour is a face of the other colour's diagonal grid (32 vertices, 49 edges, 18 bounded faces).

Moves: wazir W (orthogonal step), ferz F (diagonal step), king K = W + F; rook R, bishop B, queen Q slide along straight lines through any seam until a wall. On a closed line nothing obstructs the slide except the piece's own square, so the slider reaches the whole line.

## 3. Theorem 1: the rings are the bishop's mobility shells (PROVED)

On the plane board a cell of ring `k` has bishop mobility `15 - 2k` (13, 11, 9, 7) and queen mobility `29 - 2k` (27, 25, 23, 21).

*Proof.* With `a = x - 7/2, b = y - 7/2` the two diagonals through the cell have lengths `8 - |a-b|` and `8 - |a+b|`, and `|a-b| + |a+b| = 2 max(|a|,|b|) = 2k - 1`. So the bishop sees `14 - (2k-1) - 1 = 15 - 2k` squares; add the rook's constant 14. QED.

**Mechanism.** The king's (Chebyshev) norm is the taxicab norm in diagonal coordinates. The rings that the orthogonal moves draw are the level sets of the total scaffold length through a cell, `L_diag + L_anti = 17 - 2k`.

## 4. The boards

A gluing is a group `G` of lattice isometries `c -> Ac + t` (A a signed permutation) acting freely on cells with the window as fundamental domain; unglued edges are walls. Sixteen boards were built and checked (`chi` = Euler characteristic of the glued square complex, computed from the identifications; `|Isom|` = automorphisms of the square-tiled surface, by developing maps and independently by VF2 on the W/F edge-coloured graph):

| board | gluing | surface (orbifold) | chi | colour char. | \|Isom\| | cell orbits | cone points |
|---|---|---|---|---|---|---|---|
| plane | none | disk | 1 | - | 8 | 10 | - |
| cylinder | left-right translate | annulus | 0 | 0 | 32 | 4 | - |
| mobius | left-right flipped | Mobius band | 0 | 1 | 32 | 4 | - |
| mobius_diag | top->right by a diagonal glide | Mobius band | 0 | 0 | 2 | 36 | (boundary 90 and 270 degrees) |
| torus | both pairs translate | torus (o) | 0 | 0,0 | 512 | 1 | none |
| torus_k1 | off the right edge one row down | torus (o) | 0 | 1,0 | 128 | 1 | none |
| torus_k2 | two rows down | torus (o) | 0 | 0,0 | 128 | 1 | none |
| torus_k4 | half a board down | torus (o) | 0 | 0,0 | 256 | 1 | none |
| klein | `abab^-1` (flip about a grid line) | Klein bottle (xx) | 0 | 0,1 | 64 | 2 | none |
| klein_cc | as klein, flip axis through cell centres | Klein bottle (xx) | 0 | 0,0 | 64 | 3 | none |
| klein_diag | `abba`, adjacent edges by diagonal glides | Klein bottle (xx) | 0 | 0,0 | 32 | 5 | none |
| rp2 | `abab`, antipodal boundary | projective plane (22x) | 1 | 1,1 | 8 | 10 | two of angle pi |
| sphere442 | adjacent edges by quarter turns about two corners | sphere (442) | 2 | 1,1 | 4 | 20 | pi/2, pi/2, pi |
| pillow | every edge folded at its midpoint | sphere (2222) | 2 | 0,0,0,0 | 16 | 6 | four of angle pi |
| pillow_cyl | translate + top and bottom folded | sphere (2222) | 2 | 0,0,0 | 8 | 8 | four of angle pi |

`rp2_fold` (left-right flipped, top and bottom folded) was also built: its W/F graph is isomorphic to rp2's. It is the same board seen through a window shifted by half a board, so its cone points sit at edge midpoints instead of corners.

## 5. Theorem 2: every gluing happens in ring 4 (PROVED; values FINITE-EXACT)

A seam step leaves the window from a border cell and enters at a border cell, so rings 1-3, their cycles, the 48 spokes and the king shells 4, 12, 20, 28 about the centre point are identical on all 16 boards. Ring 4 gains chords: 8 (cylinder, Mobius), 16 (tori, Klein bottles), 14 (rp2), 13 (sphere442), 12 (pillows). The inner 6x6 is a closed disk meeting ring 4 in a circle, so

`chi(ring-4 region) = chi(S) - 1`:

annulus (plane, chi 0); pair of pants (cylinder) or punctured Mobius band (chi -1); one-holed torus or Klein bottle (-1); **Mobius band** (rp2, 0); **disk** (the spheres, +1). Computed directly on all 16 boards.

**Observer reading.** The centre observer sees three honest rings; the fourth ring is where its view wraps (the cut locus of the centre is the image of the board's boundary). The gluing type is exactly how that ring is stitched.

## 6. Theorem 3: curvature shows up as ring growth (PROVED; FINITE-EXACT check)

Around a lattice point where `m` cell corners meet (cone angle `m * 90` degrees), the king-distance shells from the `m` cells touching it have sizes `m(2j+1)`:

- `4, 12, 20, 28` at a regular point (m = 4, the owner's rings);
- `2, 6, 10, 14` at a half-turn point (m = 2, angle pi);
- `1, 3, 5, 7` at a quarter-turn point (m = 1, angle pi/2).

These hold until the shells reach other cone points.

*Proof.* Near the point, the developing map identifies the board with the plane modulo the rotation group of order `4/m`. That group acts freely on the cells of each square shell of size `4(2j+1)`. QED. Gauss-Bonnet in this setting reads `sum over vertices of (4 - m) = 4 chi`: rp2 has two points with m = 2, giving 2 + 2 = 4 = 4 x 1, and sphere442 has 3 + 3 + 2 = 8 = 4 x 2. Observed full sequences:

- from a cone point of rp2 or the pillow: **2, 6, 10, 14, 14, 10, 6, 2**. The board seen from a cone point is a spindle ending at the antipodal cone point.
- from the quarter-turn points of sphere442: **1, 3, 5, 7, 9, 11, 13, 15**. These are the L-shaped rings `max(x,y) = k`, which the gluing closes into odd cycles `C1, C3, ..., C15`.

## 7. Theorem 4: one parity bit decides whether the scaffolds fuse; the bishop lives on the colour cover (PROVED)

For `g(c) = Ac + t` put `chi(g) = t_x + t_y mod 2`. Since `A` is a signed permutation, `parity(g(c)) = parity(c) + chi(g)`, and `chi` is a homomorphism `G -> Z/2`. The following are equivalent:

- (i) chi vanishes on the gluing maps;
- (ii) the checkerboard colouring survives;
- (iii) the wazir graph is bipartite;
- (iv) the ferz graph has two components;
- (v) the bishop graph has two components.

*Proof.* A step across a seam glued by `g` lands on the window cell `f'` with `g(f') = f + d`. So window parities change by `|d| + chi(g)` mod 2: an orthogonal step changes parity by `1 + chi(g)` and a diagonal step by `chi(g)`.

- If chi vanishes on the generators, it vanishes on `G`. Then every W edge changes colour and every F edge keeps it.
- If a generator has `chi(g) = 1`, a diagonal step across its seam joins the two colour classes of the window, which are each ferz-connected. An orthogonal step across it joins two cells of equal parity; closing with an even path inside the window gives an odd cycle.

QED. Verified on all 16 boards.

For the 8x8 window:

- translations by `(8,k)` have `chi = k mod 2`;
- flips about grid lines and quarter-turns about lattice points have chi = 1;
- half-turns about lattice points, flips about cell-centre lines, and glides with diagonal axes have chi = 0.

So **fusion is not topological**:

| | scaffolds separate | scaffolds fused |
|---|---|---|
| orientable | plane, cylinder, torus, torus_k2, torus_k4, pillow, pillow_cyl | torus_k1, sphere442 |
| non-orientable | mobius_diag, klein_cc, klein_diag | mobius, klein, rp2 |

**Refinement (FINITE-EXACT).** When the scaffolds stay separate they are congruent on every board except `klein_diag` and `mobius_diag`. There the glide axis runs along a diagonal of one colour. That colour's scaffold contains the one-sided core curve and is non-bipartite, while the other scaffold is bipartite (and on klein_diag the two have different ferz diameters, 5 and 6). The two scaffolds are isolated but no longer alike.

**Theorem 5 (double-cover principle).** Let `G0 = ker chi` and `X0 = Z^2/G0`.

- One colour class of `X0` maps bijectively onto `X`, and the map is a local isometry. So `F(X)` is the ferz graph of that colour class, and `B(X)` is its bishop graph (when chi = 0 these are the two components).
- In 45-degree coordinates the colour class is the rotated lattice `L` and `G0` acts on it by lattice isometries. Hence `F(X) = W(L/G0)` and `B(X) = R(L/G0)`: the diagonal world of a board is the orthogonal world of its rotated colour cover.

*Proof.* Path lifting for the double cover, plus the fact that the 45-degree rotation conjugates D4 to itself. QED. FINITE-EXACT checks by graph and line-multigraph isomorphism:

| board | colour cover | rotated lattice | F = W(cover) | B = R(cover) |
|---|---|---|---|---|
| torus | Z^2/<(8,0),(0,8)> | <(4,4),(4,-4)> | yes, yes | yes, yes |
| torus_k2 | Z^2/<(8,2),(0,8)> | <(5,3),(4,-4)> | yes, yes | yes, yes |
| torus_k4 | Z^2/<(8,4),(0,8)> | <(6,2),(4,-4)> | yes, yes | yes, yes |
| torus_k1 | Z^2/<(16,2),(0,8)> | <(9,7),(4,-4)> | yes | yes |
| klein | Z^2/<(8,0),(0,16)> | <(4,4),(8,-8)> | yes | yes |

The klein_cc board, whose colour is preserved, correctly fails against the Klein cover (hostile control).

- **The Klein bottle's bishop lives on a torus**, `Z^2/<(4,4),(8,-8)>`. Every axis line of that torus has length 16, so it has no 8x8 window: it is not one of the window tori above.
- **RP^2 and the 442 sphere share their colour cover.** In both cases `G0` is the p2 group generated by the translations `16Z^2` and the half-turns about `8Z^2`. Its quotient is the 128-cell pillowcase: two chessboards sewn back to back along all four edges. RP^2 is that pillowcase modulo an antipodal glide, and the sphere is the pillowcase modulo a quarter turn; both maps swap the colours. Hence `F(rp2) = F(sphere442)` (isomorphic). Their bishop graphs are even **literally equal** as sets of attacked squares, because in both gluings the diagonal `x - y = d` (d not 0) closes into one line with the same four board segments `{x-y = d, x-y = -d, x+y = d-1, x+y = 15-d}` (traced by hand; FINITE-EXACT). The two main diagonals fold back at cone points in both. A bishop cannot tell the projective plane from the sphere.

## 8. Sliding: lines become loops, and sliders only see which segments were joined

Line census (closed or open; length; folds = the line bounces back at a half-turn cone point):

| board | orthogonal lines | diagonal lines |
|---|---|---|
| plane | 16 open x 8 | lengths 1..8 (4,4,4,4,4,4,4,2) |
| cylinder | 8 closed x 8, 8 open x 8 | 16 open x 8 (helices) |
| mobius | 4 closed x 16, 8 open x 8 | 16 open x 8 |
| torus | 16 closed x 8 | 16 closed x 8 |
| torus_k1 | 8 closed x 8, **1 closed x 64** | **2 closed x 64** |
| torus_k2 | 8 x 8, 2 x 32 | 4 x 32 |
| torus_k4 | 8 x 8, 4 x 16 | 8 x 16 |
| klein | 8 x 8, 4 x 16 | 8 x 16 |
| klein_cc | 10 x 8 (two one-sided columns), 3 x 16 | 8 x 16 (14 cells each) |
| klein_diag | 8 figure-eights (16 steps, 15 cells) | 2 x 4, 15 x 8 |
| rp2 | 8 x 16 | 7 x 16, 2 folds x 8 |
| sphere442 | 8 figure-eights (16 steps, 15 cells) | 7 x 16, 2 folds x 8 |
| pillow | 8 x 16 | 14 x 8, 4 folds x 4 |
| pillow_cyl | 8 x 8, 4 x 16 | 6 x 16, 4 folds x 8 |

Mechanisms (PROVED by tracing the gluing):

- A shifted translation concatenates rows: `gcd(8,k)` lines of length `64/gcd(8,k)`.
- A flip pairs row `y` with row `7-y`; since 8 is even no row is fixed, so pairs make loops of 16.
- A quarter turn continues row `y` as column `y`, giving a figure-eight through the diagonal cell `(y,y)`.
- A half-turn cone point folds the line back.

**Theorem 6 (slider = line graph).** Every cell lies on exactly one line per axis. Build the **line multigraph**: vertices are lines, and each square is an edge (or loop) joining its two lines. The rook (resp. bishop) graph is the line graph of that multigraph. On the plane, each colour's bishop multigraph is bipartite with degree sequences **{1,3,5,7,7,5,3,1}** and **{2,4,6,8,6,4,2}**: the owner's scaffold numbers are exactly the degrees of this multigraph. Consequently a slider's attack sets depend only on **which board segments are joined into one line**, not on how (translation, flip, fold, quarter turn) or whether the line closes. FINITE-EXACT consequences:

- **identical attack sets:**
  - rook: plane = cylinder = torus; rp2 = pillow; klein = pillow_cyl; klein_diag = mobius_diag;
  - bishop: cylinder = torus; rp2 = sphere442;
  - queen: cylinder = torus only.
- **isomorphism classes of line multigraphs:**
  - rook, 7 classes on the 15 distinct boards: {plane, cylinder, torus}, {mobius, torus_k4, klein, pillow_cyl}, {mobius_diag, klein_diag, sphere442}, {rp2, pillow}, and three singletons;
  - bishop, 13 classes.
- **Steppers remember the geometry.** The W and K graphs of the 15 distinct boards are pairwise non-isomorphic. Sliding forgets everything except the line partition.
- On the helical torus (`torus_k1`) one row-line and one line per diagonal direction run through all 64 squares, so **R = B = Q = K64**: every slider sees every square.

(The cylinder and the torus agree for every slider because a helix on the 8x8 cylinder has exactly 8 cells, the same cells as the torus diagonal it covers. Only the walls, i.e. the steppers, tell them apart.)

## 9. Placements, domination, colouring (FINITE-EXACT)

Maximum non-attacking pieces and the number of maximum placements. The DP over blocked future squares is exact. Every count up to 20000 was re-enumerated by SAT, and SAT certified that no larger set exists:

| board | kings | rooks | bishops | queens |
|---|---|---|---|---|
| plane | 16 (281571) | 8 (40320) | 14 (256) | **8 (92)** |
| cylinder | 16 (4460) | 8 (40320) | 8 (147456) | 6 (3072) |
| mobius | 16 (1636) | 4 (26880) | 8 (40320) | 4 (3328) |
| mobius_diag | 16 (2984) | **8 (1)** | 8 (771696) | 4 (1782) |
| torus | 16 (60) | 8 (40320) | 8 (147456) | 6 (3072) |
| torus_k1 | 16 (32) | 1 (64) | 1 (64) | 1 (64) |
| torus_k2 | 16 (32) | 2 (896) | 2 (1024) | 2 (384) |
| torus_k4 | 16 (36) | 4 (26880) | 4 (16384) | 4 (256) |
| klein | 16 (32) | 4 (26880) | 4 (6144) | 3 (1024) |
| klein_cc | 16 (40) | 5 (53760) | 8 (256) | 4 (1792) |
| klein_diag | **16 (2)** | **8 (1)** | 8 (258048) | 4 (1600) |
| rp2 | 16 (4) | 4 (6144) | 4 (10752) | 3 (768) |
| sphere442 | 16 (62) | **8 (1)** | 4 (10752) | 4 (128) |
| pillow | **16 (2)** | 4 (6144) | 8 (451584) | 4 (1088) |
| pillow_cyl | 16 (4) | 4 (26880) | 4 (25600) | 4 (768) |

Proved pieces:

- **Kings: always 16.** The window's sixteen 2x2 blocks are king cliques on every board. The count drops from 281571 to 2, so gluing rigidifies king packings.
- **Rooks.** At most one rook per line. When 4 lines of 16 each meet 8 lines of 8 twice (mobius, klein, torus_k4, pillow_cyl), the count is `C(8,4) 4! 2^4 = 26880`. When 4 row-loops meet 4 column-loops in 4 cells each (rp2, pillow), it is `4! 4^4 = 6144`.
- **Rooks on the figure-eight boards: exactly one maximum placement.** There are 8 lines, a cell on a self-crossing uses one line, every other cell uses two, so 8 rooks force all of them onto the eight crossings (the fold diagonal).
- **Only the plane holds 8 queens.** On the torus (and hence the cylinder) Polya's parity argument applies. If the rows, columns, `x+y` and `x-y` were all permutations mod 8, then `sum(x+y) = 2 sum(x)` would force `sum_{i<8} i = 28 = 0 mod 8`, which is false. The maxima 6 (torus; Monsky, A085801 a(8) = 6), 4 (Mobius) and 3 (Klein) agree with Arizmendi Echegaray's 2026 Bridges piece. The rp2 maximum is also 3.

Domination numbers `gamma` (minimum pieces covering every square) and chromatic numbers `chi_col`, all by SAT (glucose4); the table is in section 9a below.

## 9a. Domination and colouring

`gamma` = minimum pieces attacking or occupying every square, with the number of minimum sets where the enumeration finished (capped at 5000 or 20000). `chi_col` = chromatic number of the move graph. All values are SAT (glucose4 / CaDiCaL).

- Every queen chromatic number below is pinned exactly: SAT colours with `ceil(64/alpha)` (or `omega`) colours, and that number is a proven lower bound.
- For the torus and cylinder, `gamma(B) = 8` is by proof rather than SAT. Per colour, every diagonal meets every anti-diagonal, so an empty diagonal and an empty anti-diagonal would leave their intersection undominated. A dominating set therefore fills all four of one family in each colour.

| board | gamma K | gamma R | gamma B | gamma Q | chi W | chi F | chi K | chi R | chi B | chi Q |
|---|---|---|---|---|---|---|---|---|---|---|
| plane | 9 (3600) | 8 | 8 | 5 (4860) | 2 | 2 | 4 | 8 | 8 | 9 |
| cylinder | 9 (>=5000) | 8 | 8* | 4 (832) | 2 | 2 | 4 | 8 | 8 | 11 |
| mobius | 9 (>=20000) | 4 | 8 | 4 (>=20000) | 3 | 2 | 4 | 16 | 8 | 16 |
| mobius_diag | 9 (11103) | 4 | 8 | 3 (546) | 2 | 3 | **5** | 15 | 8 | 16 |
| torus | **8** (16) | 8 | 8* | 4 (832) | 2 | 2 | 4 | 8 | 8 | 11 |
| torus_k1 | 8 (64) | 1 | 1 | 1 (64) | 3 | 2 | 4 | 64 | 64 | 64 |
| torus_k2 | 8 (160) | 2 | 2 | 2 (1536) | 2 | 2 | 4 | 32 | 32 | 32 |
| torus_k4 | 8 (32) | 4 | 4 | 2 (384) | 2 | 2 | 4 | 16 | 16 | 16 |
| klein | 9 (>=5000) | 4 | 4 | 3 (>=5000) | 3 | 2 | 4 | 16 | 16 | **22** |
| klein_cc | 9 (>=20000) | 5 | 4 | 3 (288) | 2 | 2 | 4 | 16 | 14 | 24 |
| klein_diag | 8 (18) | 4 | 8 | 3 (744) | 2 | 3 | **5** | 15 | 8 | 16 |
| rp2 | 9 (>=5000) | 4 | 4 | **2 (48)** | **4** | 3 | 4 | 16 | 16 | **22** |
| sphere442 | 9 (>=5000) | 4 | 4 | 3 (>=5000) | 3 | 3 | **5** | 15 | 16 | 19 or 20 (pending) |
| pillow | 8 (12) | 4 | 8 | 3 (640) | 2 | 3 | **5** | 16 | 8 | 16 |
| pillow_cyl | 8 (8) | 4 | 4 | 3 (832) | 2 | 3 | 4 | 16 | 16 | 16 |

(* by the proof above.)

**Readings.**

- **Two queens dominate the projective plane** (48 ways). A queen's lines there are 16-loops, so two well-placed queens see everything. The plane needs 5 and the torus 4, matching OEIS A279402 (Burger-Mynhardt).
- **Colouring the queens of the Klein bottle or projective plane takes 22 colours** (plane 9, torus 11). This is forced by `alpha(Q) = 3`, since `ceil(64/3) = 22`, and SAT meets the bound.
- **The orthogonal-step graph of the projective-plane board needs 4 colours, never 3.** This is Youngs' theorem (J. Graph Theory 21 (1996) 219-227): a non-bipartite quadrangulation of the projective plane is 4-chromatic. Our wazir graph is such a quadrangulation once the double edge at each cone point is merged into a neighbouring square. The non-bipartite Klein bottle, Mobius band and helical torus manage with 3.
- **The king needs a fifth colour on pillow, klein_diag, mobius_diag and sphere442.**
  - *Mechanism, boards without cone points (klein_diag):* a proper 4-colouring of the infinite king graph is row-periodic or column-periodic (if a row shows three consecutive distinct colours, the columns alternate). The lift of a board colouring must be G-invariant, and a diagonal glide swaps the two families; the only colourings in both families are the four-colour 2x2 patterns, and the glide breaks those too.
  - *Boards with cone points:* the lift is improper at the cone point, because the diagonal step through it is a loop, so the dichotomy argument does not apply. Here SAT decides: rp2 and pillow_cyl are 4-colourable, while pillow and sphere442 are not (sphere442 has a 5-clique at its quarter-turn points).

## 10. Distances: how much teleportation shrinks the world (FINITE-EXACT)

| board | wazir diameter / mean | king diameter / mean |
|---|---|---|
| plane | 14 / 5.333 | 7 / 3.750 |
| cylinder | 11 / 4.698 | 7 / 3.286 |
| mobius | 7 / 4.318 | 7 / 3.206 |
| torus | 8 / 4.063 | 4 / 2.730 |
| torus_k4 | 6 / 3.873 | 4 / 2.730 |
| klein | 7 / 3.937 | 4 / 2.730 |
| klein_diag | 8 / 4.143 | 6 / 2.873 |
| rp2 | 7 / 4.060 | 7 / 2.913 |
| sphere442 | **14** / 4.630 | 7 / 3.240 |
| pillow | 8 / 4.333 | 7 / 3.095 |

Every slider has diameter 2, or 1 on the helical torus.

- **The 442 sphere keeps the plane's wazir diameter 14.** Its gluings fold each edge onto the adjacent edge at matching distance from the corner, so nothing far away is brought close.
- **King mean exactly 172/63 on torus, torus_k1, torus_k2, torus_k4, klein and klein_cc.** On these boards the king balls of radius 3 embed: 49 cells, then 15 cells at distance 4. klein_diag is the exception (diameter 6).

## 11. Loss ledger (quotient typing)

- **Source:** the infinite board Z^2 with its moves.
- **Target:** `X = Z^2/G`.
- **Map:** projection.
- **Preserved:** every local move away from cone points, and every line as an image.
- **Lost:** the deck coordinate. In particular the parity of a square is lost exactly when `chi` is nonzero.
- **Restoring sidecar:** the colour cover `Z^2/ker chi`, which is where the bishop lives (Theorem 5).
- **Second-level loss:** the scaffold's own checkerboard (column parity), which half-turns and diagonal glides destroy. Ferz components are non-bipartite on rp2, sphere442, both pillows, and on one scaffold of klein_diag and mobius_diag.
- **Third level:** the two families of proper king 4-colourings (row-periodic and column-periodic). Without cone points, the king graph needs a fifth colour exactly when no G-invariant member of either family exists (klein_diag). With cone points the lift is improper and SAT decides (section 9a).

## 12. Verification

- **Positive controls (classical plane values):** 92 eight-queen solutions; 16 kings in 281571 ways; 14 bishops in 256 ways; 5 queens dominate in 4860 ways; chi(Q8) = 9. Toroidal controls: 8 queens impossible and maximum 6 (A085801); queen domination 4 (A279402).
- **Independent paths:**
  - independence counts by DP versus SAT enumeration (all 60 cases up to 20000 agree);
  - isometry groups by developing maps versus VF2 automorphism enumeration (16/16 agree);
  - Euler characteristics from the glued complex versus the known surfaces.
- **Hostile controls:**
  - klein_cc and klein_diag: non-orientable, chi = 0, scaffolds separate;
  - the klein_cc cover test fails as it should;
  - rp2_fold is isomorphic to rp2 (same board, different window).

## 13. Open threads

- **General n.** For odd n the translation by n flips colour, so the torus fuses its scaffolds and Polya allows n queens iff gcd(n,6) = 1. Map the fusion table and the RP^2 / sphere queen numbers as functions of n.
- **Mirror (billiard) edges.** Reflection in a grid line has chi = 1, so a reflecting bishop changes colour. Compare with the fairy-chess reflecting bishop.
- **Cube surface (six boards).** Three cells meet at each corner (cone angle 270 degrees), so the diagonal step through a corner is undefined: discrete curvature that the ferz cannot cross.
- **A count-level explanation** of the king-packing numbers: 60 on the torus, 2 on the pillow and klein_diag, 62 on the 442 sphere.
