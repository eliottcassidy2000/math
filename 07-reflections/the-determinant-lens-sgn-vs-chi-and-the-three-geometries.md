# The Determinant Lens: sgn vs χ, and the Three Geometries of a Tournament

**Source:** mac-mini-2026-06-10-S2 (session theme: Hadamard matrices × tournaments ×
odd functions × simplicial geometry). Companion canon: THM-468 (floor), THM-472
(ceiling), THM-473 (average); HYP-2383..2389; T777.

---

## One identity, three geometries

For a tournament T, the skew matrix S = A - A^T is its **odd part** — literally: a
tournament is an odd ±1 function on ordered pairs (s(j,i) = -s(i,j)), just as a graph
is an even one. The single identity

    det(I + S)  =  Σ_{K even}  Pf(S[K])²

(sum over all even vertex subsets, every Pfaffian odd) carries three different
geometries depending on which question you ask of it:

- **FLOOR (circle geometry).** det = 2^(n-1), every Pfaffian minor ±1 ⟺ T is a
  local order (vortex-free, locally transitive) ⟺ T is a switching of a transitive
  tournament. The switching class of the linear order is the *circular* order —
  Knuth's sweep-line tournaments, Lachlan's dense local order, one of just three
  countable homogeneous tournaments. Attribution: the class equivalences are
  Babai–Cameron (EJC 2000, Lemma 3.3) and Knuth (Axioms & Hulls §4); the principal-
  minor characterization is Boussairi–Ezzahir–Lakhlifi–Mahzoum (Discrete Math 2023,
  class D₁); the det(I+S) floor phrasing and the (2n-2)!! cross-count are the
  project's repackaging (THM-468, with an independent ε-induction proof).

- **CEILING (Hadamard geometry).** det ≤ (n+1)^((n-1)/2) (odd n) with equality
  exactly on switching classes of doubly regular tournaments = skew-Hadamard
  matrices of order n+1 (Reid–Brown); det ≤ n^(n/2) (even n) with equality exactly
  on skew conference matrices. THM-472. The rows of I+S are vertices of the cube
  {±1}^n; the tournament IS a simplex inscribed in the cube; vol = det/n!; the
  ceiling case is the *regular* simplex (Gram = (n+1)I - J, edge √(2n+2)). The
  Hadamard conjecture is literally "which cubes inscribe regular simplices"
  (Grigor'ev), and the skew-Hadamard conjecture is "every n ≡ 3 (mod 4) attains the
  tournament ceiling" — smallest open order 355.

- **AVERAGE (tableau combinatorics).** E[det(I+S)] over uniform random tournaments
  = the involution numbers A000085 = total count of standard Young tableaux;
  E[det(xI+S)] = the signless matching polynomial of K_n = Hermite. (THM-473;
  the char-poly form is Klanderman–Montee–Piotrowski–Rice–Shader 2024 Thm 7.7, with
  Godsil–Gutman 1981 as the symmetric ancestor; the involution evaluation appears
  unstated there.) The staircase's own objects — matchings, involutions, tableaux —
  emerge from the tournament ensemble's spectrum.

Floor = the circle. Ceiling = the regular simplex. Average = the Young lattice.
One Pfaffian expansion away from each other.

## sgn versus χ: the two poles of odd-function space

Every rotational tournament on Z_n IS an odd ±1 function f (connection set J with
J ∪ -J = Z_n \ {0}), and det(I+S) = Π_k (1 + t_k²) where i·t_k is the odd Fourier
spectrum. Two odd functions are canonical:

- **sgn** (the carousel: J = {1, …, (n-1)/2}) — the *smoothest* odd function, the
  circular signum. It is locally transitive: **d = 1, the floor.**
- **χ** (the quadratic character: J = QR, q ≡ 3 mod 4) — the *flattest* odd
  function (|t_k|² ≡ n, perfect Legendre autocorrelation). It is doubly regular:
  **d = max, the ceiling.** χ is odd precisely when q ≡ 3 (mod 4) — the same parity
  of χ(-1) that decides Paley graph vs Paley tournament. The Paley construction is
  an even/odd fork, and the tournament side is the odd side.

So the determinant range of tournament space runs **from sgn to χ** — from the
multiplicatively structured character to the order-theoretically structured signum,
the two most famous odd functions in mathematics.

And the Hamiltonian path count CHANGES ITS ALLEGIANCE between them:

- n ≤ 11: H-max = QR (both crowns on χ: H = 189 at n=7 unique = QR_7; 95095 at
  n=11 = QR_11 = A038375(11)).
- n ≥ 13 (circulants; n ≥ 19 vs Paley directly): **the carousel takes over** —
  H(carousel) > H(QR-type): 3711175 > 3703011 at n=13; H(R_19) > H(both DRT(19)s)
  > H(QR_19). The H-crown migrates from the Hadamard ceiling to the circle floor.

Heuristic worth pursuing: H rewards spectral *concentration* (low-frequency, smooth,
1/f-like odd functions), det rewards spectral *flatness* (white-noise odd functions).
The two extremal problems are Fourier-antipodal, which is *why* d ⊥ H.

## d is a genuinely new coordinate on the metagraph

Attached per iso class (n = 5..9 complete; 191536 classes at n=9):

- d ⊥ H in bulk: Pearson(d,H) = -0.20, -0.06, 0.003, -0.0002, 0.0025 at n = 5..9.
- H is ≈ the second eigenvector of G_n (R² ≈ 0.73–0.83, reconfirmed); **d is NOT**
  (R² ≤ 0.004). The metagraph has at least two nearly orthogonal smooth coordinates:
  H (the score/hierarchy axis) and d (the switching/spectral axis).
- Within a fixed score multiset, d and H correlate *positively* (+0.65/+0.37/+0.33
  pooled at n=5/6/7): once you quotient out the hierarchy, volume and path count
  pull together. Score explains their global decoupling.
- At n=7 the two ceilings touch: argmax H (= QR_7, unique) sits INSIDE argmax d
  (= the 6-class switching family of QR_7, H ∈ {45, 69, 87, 135, 189}). Switching
  pins d and scrambles H. At n=8, 9 they separate completely (H-max has d = 2..9
  vs d-max 32; at n=9 the H-champion is the regular circulant {1,2,3,5}, d = 19).
- d descends to switching classes, hence — via the gauge bijection below — to
  TILINGS.

## Tilings ARE switching classes (gauge-fixing)

Every switching class of tournaments contains EXACTLY ONE tournament containing the
fixed base Hamiltonian path P₀ (any nontrivial cut reverses some base-path arc).
Counts confirm: 2^(C(n,2))/2^(n-1) = 2^(C(n-1,2)) = #tilings. So the project's
staircase tiling hypercube Q_m IS the space of switching classes (oriented two-graphs
= skew two-graphs), and THM-447's "normalizing a skew-Hadamard matrix = fixing the
Hamiltonian path" is the same statement at the matrix level. Three quotient layers
now stack:

    tilings (= switching classes, labeled)  →  G_n (iso)  →  G_n/Z_2 (merged)
    and a NEW coarser layer: S_n = switching classes up to iso (A049313:
    1, 1, 2, 2, 6, 12, 79 for n = 2..8; Babai–Cameron count over LEVEL
    permutations — constant 2-adic valuation across cycles: the dyadic seam again).

Even/odd completion: labeled even graphs (cycle space, 2^(C(n-1,2))) ↔ labeled
switching classes of tournaments (same count); iso-level: two-graphs = A002854 =
V(E_n) (the project's even-graph metagraph counts!) vs oriented two-graphs = A049313.
Mallows–Sloane gives two-graphs ↔ Euler graphs; **the "odd Mallows–Sloane partner"
for A049313 is an open hunt** — the project's tiling/E_n machinery is the natural
place to look.

## The odd map and the hairy ball: Rédei as transversality

x ↦ Sx is an odd tangent vector field on every sphere (xᵀSx = 0 for skew S). At odd
n its zeros on S^(n-1) are exactly ±w/|w|, where adj(S) = wwᵀ — **the Pfaffian vector
is the hairy-ball singularity of the tournament's odd field**. And Σw_i =
√det(I+2A) is ODD (THM-174, Rédei-adjacent parity), so w ∉ 1^⊥: *the singularity
never lies on the sum-zero sphere.* Rédei parity is a transversality guarantee.

New computational law (HYP, proof outstanding): **mod-4 score law**
adj(S)_ij ≡ (-1)^(s_i + s_j) (mod 4) at odd n (exhaustive n ≤ 7, random n ≤ 11);
even n: (adj(S)·1)_i ≡ +(-1)^(s_i) (mod 4) with forced positive sign. The naive
mod-8 lift fails; sign(w_i)·(-1)^(s_i) is not constant. The scores control the
adjugate's 2-adic bottom bits — another dyadic seam.

Tucker/Ky Fan quadrant (confirmed empty in the literature): Fan's alternating
simplices are antidirected monotone chains in the transitive tournament on label
magnitudes; the tournament-side parity results that exist (Rédei; Forcade 1973 at
n = 2^k; El Sahili–Abi Aad 2020: antidirected Hamiltonian paths ≡ 2 mod 4 at even
order, Grünbaum) have no Fan-style topological formulation. OPEN-Q material: a
"tournament Ky Fan" replacing the magnitude order by an arbitrary tournament.
Piquant data point: QR_7 is simultaneously the unique H-maximizer at n=7 AND one of
Grünbaum's three antidirected-path exceptions — maximal forward parity, deficient
antidirected parity.

## The mod-4 stratification and the borrowed eigenvalue

The whole maximal-determinant world splits by n mod 4 (Gram entries ≡ n mod 4 — the
2-adic seam in its classical habitat). Tournament side of the table:

| n mod 4 | ceiling | attained iff | status |
|---------|---------|--------------|--------|
| 3 | (n+1)^((n-1)/2) | DRT exists (skew-Hadamard n+1) | conj. always; A001121 census 1,1,1,2,2,37,722, ≥13330 at 31 |
| 0 | n^(n/2) | skew conference (skew-Hadamard n) | conj. always |
| 2 | EW bound 2(n-1)(n-2)^((n-2)/2) = D(n) | 2n-3 a perfect square (Armario–Frau; Greaves–Suda 2017) | n = 6 ✓ (our 160), 14, 26, 42 |
| 1 | **2(n-1)^((n-1)/2)** (NEW conjecture, HYP-2389) | spectrum (n-2)^((n-3)/2)·(2n-3) | n = 5 (32), 9 (8192 = 2^13, 216 classes, ALL sharing x(x²+7)³(x²+15)) — exhaustively exact |

The n ≡ 1 (mod 4) row appears genuinely untreated in the literature (the even rows
are the Armario-school "skew EW" theory; the odd-3 row is Reid–Brown). The extremal
spectra at the unattainable orders are **two-level with one excited pair at 2n-3**
over a flat base (n-2 or n-3): integrality forbids flatness, and the obstruction
materializes as exactly one borrowed eigenvalue. (At n=9 none of the 216 maximizers
is regular — 69 score sequences, one spectrum. The spectral class is pure; the
combinatorial fibers are wild.)

## What the lens does NOT see

- d is coarse (7 values at n=7) and shape spectra don't separate classes.
- d sees nothing of H's fine structure (the whole point: it's the other axis).
- The d-spectrum has its own forbidden values (d = 7 missing at n=7; {23,25,28,29,31}
  missing at n=9; 20→32 jump at n=8) — a gap phenomenon rhyming with H ∉ {7,21},
  unexplained.

## Engineering hooks

- d(T): cheap (one determinant), switching-invariant, iso-invariant tournament
  feature — add to tournament_tda.py alongside Seidel spectra and Moorhouse
  fingerprints.
- E[det(I+S)] = A000085: an exact null-model statistic — instant bias test for any
  tournament sampler/enumerator.
- OEIS: max-det sequences (2,4,16,32,160,512,4096,8192·…; d-version 1,1,2,2,5,8,32,32)
  are absent — submission candidates; A334123 (labeled max-H counts) extendable from
  our census (a(7) = 240 = 7!/21; a(8), a(9) computable from class data).
- D-optimal designs: tournament matrices are the skew-type D-optimal candidates;
  ratios to the true max D(n): 1, 2/3, 1, 8/9, 1 at n = 4..8.
