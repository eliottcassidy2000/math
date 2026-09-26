---
id: HYP-9162
title: "The Sierpinski tournament tower: H_2 = [[1,1],[-1,1]], H_(2n) = [[H,H],[-H^T,H^T]] gives skew Hadamard matrices of every order 2^k and hence doubly regular tournaments T_k of every Mersenne order 2^k - 1 (the zero-triangle sides of Gilbreath's single-seed sea); T_3 is the Paley heptagon, Aut(T_k) contains its Frobenius group F_21 for all k >= 3, and the conjecture is that Aut(T_k) = F_21 exactly for all k >= 3 (so T_k is never Paley for k >= 4, in particular T_5 != P_31)"
status: >
  PROVED: the doubling preserves skewness and orthogonality (two-line check),
  double regularity is the classical skew-Hadamard correspondence, the
  recursion T_(k+1) = T_k + {0'} + T_k' is explicit, and every automorphism of
  T_k extends diagonally, so F_21 <= Aut(T_k) for k >= 3. FINITE-EXACT:
  |Aut(T_k)| = 3, 21, 21, 21, 21 for k = 2..6 (and see the .out for k = 7, 8),
  T_3 = P_7, T_5 not isomorphic to P_31 (|Aut(P_31)| = 465) though both have
  the same 4-vertex census. OPEN: equality Aut(T_k) = F_21 for all k >= 3.
source: opus-2026-09-26 session gilbreath-fermat-platonic-20260926
related: [THM-871 (Fermat-rung rigidity of rotational tournaments), HYP-3805 (Paley heptagon as extremal object), S9 note collatz_oscillation_gilbreath_20260926.md (zero triangles of the sea)]
verification: 04-computation/experiments/gilbreath_fermat_tower_20260926_tournaments.py -> .out
---

# HYP-9162 -- the Sierpinski tournament tower keeps exactly the heptagon's 21 symmetries

**Construction.** `H_2 = [[1, 1], [-1, 1]]`, `H_(2n) = [[H, H], [-H^T, H^T]]`.
If `H = S + I` with `S` antisymmetric then `H_(2n) - I = [[S, S+I], [S-I, -S]]`
is antisymmetric and `H_(2n) H_(2n)^T = 2n I`; normalization is inherited.
Deleting row and column `0` gives the tournament `T_k` on `2^k - 1` vertices,
doubly regular (checked to `n = 255`), with the recursion

```text
T_(k+1) = T_k + {0'} + T_k' (arcs reversed),  i -> i',  i -> j' iff i -> j,  i' -> j iff i -> j,  0' -> T_k,  T_k' -> 0'.
```

**Computed.** `|Aut(T_2)| = 3`; `T_3 = P_7`, `|Aut| = 21`; `|Aut(T_4)| =
|Aut(T_5)| = |Aut(T_6)| = 21`, with `Aut(T_5)` having element orders
`{3: 14, 7: 6}` (the Frobenius group `Z_7 x| Z_3`) and orbits of sizes
`7, 1, 7, 1, 7, 1, 7`. `T_5` is not isomorphic to `P_31`.

**Conjecture.** `Aut(T_k) = F_21` for every `k >= 3`; equivalently, no
automorphism of `T_(k+1)` fails to fix `0'` once `k >= 3`. If true, `T_k` is
not vertex-transitive for `k >= 4` and never Paley (a Paley tournament on a
Mersenne prime `2^k - 1` has `|Aut| = (2^k - 1)(2^(k-1) - 1) > 21`).

**Why it matters here.** The sides of the zero triangles of Gilbreath's
single-seed sea are exactly the Mersenne numbers `2^m - 1` (session note
`gilbreath_fermat_platonic_20260926.md`, section 2), and this tower is the
doubling of the sea itself read with signs: it makes the owner's "zeros of
tournament size edged by 2s" precise (each such size carries a doubly regular
tournament) and shows what symmetry the tower keeps: the heptagon's.
