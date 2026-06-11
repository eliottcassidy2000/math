# d(T) — the normalized tournament determinant

**Definition:** d(T) = det(I + S) / 2^(n-1), where S = A - A^T is the skew ±1 matrix
of the tournament T on n vertices. Always a positive integer.

**Key equations:**
- det(I+S) = Σ_{K ⊆ [n], |K| even} Pf(S[K])²  (every Pfaffian odd) — THM-468
- d(T) ≥ 1, equality ⟺ T locally transitive ⟺ switching class of transitive
  (THM-468; Boussairi et al. 2023 D₁; Babai–Cameron 2000)
- odd n: d ≤ ((n+1)/4)^((n-1)/2)·… i.e. det ≤ (n+1)^((n-1)/2), equality ⟺
  switching class of a DRT; even n: det ≤ n^(n/2), equality ⟺ skew conference
  (THM-472)
- E[det(I+S)] = involution numbers A000085; E[det(xI+S)] = signless matching
  polynomial of K_n (THM-473; cf. Klanderman et al. 2024)
- circulant T: det(I+S) = Π_k (1 + t_k²), i·t_k the odd Fourier spectrum of the
  connection-set symbol (odd ±1 function on Z_n)
- d(DRT) = ((n+1)/4)^((n-1)/2) — a perfect power of the skew-Hadamard parameter

**Invariance:** isomorphism, switching (S → DSD), reversal (T^op), hence descends to
G_n/Z_2 AND to tilings (= switching classes; one P₀-tournament per class).

**Values:** per-class spectra n=3..9 in 05-knowledge/results/hadamard_det_census_macmini_s2.out
and hadamard_det_n9_macmini_s2.out. Max sequence det: 2,4,16,32,160,512,4096,8192 (n=2..9);
d-max: 1,1,2,2,5,8,32,32. NOT in OEIS (checked 2026-06-11). d ⊥ H in bulk (r ≈ 0 at n=7..9);
NOT the metagraph second eigenvector (R² ≤ 0.004 vs H's 0.73-0.83).

**Geometric meaning:** vol of the tournament simplex (rows of I+S in the cube) = det/n!;
d=1 = circle geometry (sgn pole); d=max = Hadamard geometry (χ pole); squared vertex
distances = 4(1 + #disagreements) — shape IS co-degree combinatorics.

**Related:** [hamiltonian-paths.md](hamiltonian-paths.md) (H ⊥ d), [signed-matrix.md](signed-matrix.md),
THM-174 (det(I+2A) = Pf²), THM-442 (H² − Pf² = 8Q), HYP-2389 (tournament Barba),
T777, reflection the-determinant-lens-sgn-vs-chi-and-the-three-geometries.md.
