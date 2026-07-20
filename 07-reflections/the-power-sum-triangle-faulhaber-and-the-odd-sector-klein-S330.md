# The power-sum triangle, Faulhaber's third perspective, and the odd sector

**Instance:** klein-2026-07-20-S330 (owner: work the S329 next steps + identify the
triangle {1},{2,1},{3,3,1},{4,6,5,1},{5,10,14,9,1},{6,15,30,37,17,1},{7,21,55,101,99,33,1};
relate to Fibonacci, powers of 2, Moser; 2n+1, 2^x+1, the Proth table n·2^x+1; the
rational series; odd degree functions and tournaments. Long session, many pulls.)
All computations exact and frozen (script + out).

## 1. The triangle identified: a power-sum triangle with a Moser break (28/28 exact)

> **T(n,k) = PS_{k−1}(n−k+1) + [k ≥ 4][n−k ≥ 2]**, PS_p(m) = 1^p + 2^p + … + m^p.

Columns: n; triangular T_m (polygonal); square-pyramidal Σm² (polyhedral); Σm³ = T²
(Nicomachus) + break; Σm⁴ + break; … Diagonal k = n−1: 2^m + 1 (the owner's 2^x+1
axis). The +1 break switches on exactly when BOTH indices clear the small cases —
structurally a Moser-circle phenomenon (Moser: R(n) = C(n,4)+C(n,2)+1 tracks 2^{n−1}
until n = 6; here the pure power sums hold until (k, n−k) = (4, 2)).
**Falsifiable prediction: row 8 = {8, 28, 91, 226, 355, 277, 65, 1}.**

## 2. The third perspective on triangular numbers = FAULHABER

Beyond polygonal (2D figurate) and polyhedral (3D figurate): **T is the universal
coordinate of the odd power sums** — verified symbolically: PS₃ = T², PS₅ = T²(4T−1)/3,
PS₇ = T²(6T²−4T+1)/3 — every odd power sum is T²·(polynomial in T) (Faulhaber), with
Nicomachus (PS₃ = T²) as the first instance. The triangle's columns walk exactly this
transition; its Moser break marks where the naive figurate reading fails.

## 3. Fibonacci, powers of 2, Moser, and the Proth table

- Pascal: shallow-diagonal sums = Fibonacci; row sums = 2^n; Moser = the k ≤ 4
  partial row sums (why "2^n breaks"). Our triangle's parallels: shallow sums
  1, 2, 4, 7, 12, 21, 37 (fit: a(n) = a(n−1)+a(n−2)+a(n−4)); row sums 1, 3, 7, 16,
  39, 106, 317 (ratio → 3: the "base-3 Pascal" shadow); the 2^m+1 diagonal.
- **The Proth table n·2^x + 1**: boundary rows (1,n) = 2n+1 (odd numbers) and
  (x,1) = 2^x+1 (Fermat line) — the two "keys to clarification." Repo resonances
  (tagged): 2N+1 = the LRC first-window denominators (K-ladder s/(Ns+1), window top
  2/(2N+1)); 2^x+1 ↔ the binder/gate tower; n·2^x+1 = Proth/Sierpiński arithmetic =
  the shape of death-star's D-rung gate moduli. The triangle is a figurate shadow
  of the table's boundary.

## 4. The rational series

Series 1 = the harmonic numbers H_n EXACTLY (1, 3/2, 11/6, 25/12, 137/60).
Series 2 (1, 5/2, 29/3, 109/12, 1079/60): closest identified companion is the
**row-harmonic of the triangle**, Σ_k T(n,k)/k = 1, 5/2, 29/6, 107/12, … — exact
match at n ≤ 3 under the 29/6-reading of the third term, diverging by small amounts
at n = 4, 5 as given (109/12 vs 107/12). Recorded as OPEN with both readings; the
canonical companion (row-harmonic of the verified triangle) is now computable to any
length. [Owner: clarify 29/3 vs 29/6 and the n ≥ 4 terms if the divergence matters.]

## 5. Odd degree functions and tournaments — the sector unification

Verified this session and armed by the fleet:

- **Tournaments ARE odd pair-functions** (skew ±1 adjacency) — the odd sector of the
  pair-cube; the repo's cut⊕cycle and tournament/even-graph metagraph duality is the
  odd/even decomposition itself.
- **Even score-moments count cycles**: c₃ = C(n,3) − Σ C(sᵢ,2) verified exhaustively
  (all 64 + 1024 tournaments at n = 4, 5) — the 3-cycle count is a PS₂ functional of
  the scores. **Odd score-moments are Faulhaber-polynomials in T = C(n,2) = #arcs**
  on the transitive tournament: odd = the triangle-coordinate, even = the cycle
  detector. (kps's LRC "moment blindness" — S₂,S₄-fibers with ΔM up to 1/12 — is the
  same even-moment blindness in runner clothing.)
- **The Jacobian anatomy is forced into the odd sector**: the fiber cubic is
  trace-zero (odd + constant), its core is the ODD Chebyshev T₃, the counterexample
  carries the sign equivariance F∘σ = τ∘F (graded oddness), and death-star's
  HYP-8160 now reads DC(1)'s A1 weight-triple as an ORIENTED 3-CYCLE — the Weyl
  frontier speaking tournament.

## 6. The S329 next steps, worked via pulls

- **Reconciliation done**: death-star THM-1345 (JC₂ is a THEOREM in the
  ℂ*-equivariant category, every action) + my Euler–Zariski bootstrap compose:
  a 2-variable counterexample must be simultaneously torus-asymmetric AND carry a
  cusped, χ-constrained Jelonek curve with a non-normal S₃ cover. The corner it can
  hide in shrank from two independent sides in one day.
- boxeph-S146 tied the tournament H-template to the KELLER-DEGREE MONOID
  ({7,21}-impossibility ↔ monoid law) — the bridge built from the tournament side.
- Still open and named: the candidate-Jelonek atlas (χ = 1 cusped curves with
  B₃-type groups); the degree-4 bootstrap (A₄/S₄); mac-mini HYP-8155 Hurwitz
  reconciliation (their compactified ledger vs my affine one).

## 7. Cross-links

T1547–T1549 + S323–S329 reflections (the JC arc) · THM-1345, HYP-8160, HYP-8155,
boxeph-S146 (the fleet's concurrent frontier) · everything-is-the-triangle (the
staircase = C(n−1,2) — the repo's own T-coordinate) · kps HYP-7955 (moment
blindness) · Rédei parity / locker parity (the odd sector's arithmetic) ·
CONSTANTS-INDEX (2N+1 windows; binder gates).
