# The Leech lattice seen from its deep holes: a swap theorem, h-cores, divisors of 3n+1, and bad primes (opus-2026-10-22-S53)

**Prompt (2026-10-22, overnight).**
- Pursue the next obligations: prove HYP-9363 and test the other hole types.
- Relate 3^k+1 and k^3+1 to Collatz's 3n+1, and use the pasted Steiner/Sierpinski triangle (x*y = -x-y mod 3; lengths
  3^k+1 collapse to their endpoints).
- Relate {25,1} to {7,21}; 27 in Collatz; the path of least resistance and natural recursion.

**Lanes.**
- mac-mini: the tournament side (THM-4754..4757, HYP-9364..9367). That covers Steiner-quasigroup path signatures,
  H mod 3 = I(Omega,-1), the 7/21 profiles, and H(R_7) = 7*25, H(P_7) = 7*27.
- Session29: the Collatz powers-of-3 / 3^k+1 cylinders.
- This session: the lattice/modular side.
- Pointer sent to mac-mini: in Paley tournaments over F_{3^m}, the Steiner lines {x, y, -x-y} are cyclic triangles
  (MSG to all, a073970a6e).

Scripts are in `04-computation/experiments/deep_holes_theta_20261022/`. Each prints ALL CHECKS PASSED, except
cusp_pairs, whose single failure is a genuine counterexample. All need cypari and numba.
- `swap_theorem_proof.py`: the proof certificates, all 23 types, 8 s.
- `all_holes.py`: exact enumeration of both sides, all 23 types through q^3, plus the eta-quotient checks; 3 min.
- `swap_reciprocity.py`: the off-axis swap symmetry scans.
- `swap_dimension.py`: dimensions of the antisymmetric torus-invariant cusp spaces.
- `cusp_pairs.py`: the generalisation to arbitrary cusp pairs (275 of 276 pass).

## 1. HYP-9363 is a theorem (THM-4758, PROVED for all 23 types)

For every Leech deep hole c of type R, with Niemeier lattice N, Coxeter number h, Weyl vector rho_N and k
components:

**theta(Lambda + c) - theta(N + rho_N/h) = (24 + k) Delta.**

The coefficient 24 + k is the number of hole vertices: 25 for A24 and D24, 26, 27, 28, 29, 30, 32 and 36 for the
others, and 48 for A1^24.

**Proof (sketch; details in THM-4758).**
- In II_{25,1} take the Leech cusp rho and the cusp z of N, from Borcherds' table of norm-0 vectors. Then
  M = {rho, z}^perp has discriminant group (Z/h)^2, with two isotropic lines, one from Lambda and one from N.
- The two theta series in the identity are sums of coset theta series of M over two lines exchanged by a swap sigma.
- The swap is an isometry of the discriminant form, so it commutes with the Weil representation. Hence
  E = Theta_M - sigma^*Theta_M is a vector-valued cusp form of weight 12.
- **Vector-valued Sturm lemma:** a weight-12 form whose components all vanish to order > 1 is zero. (Take the norm
  over cosets and apply the valence formula.)
- So it suffices to check two things using only the vectors of M^* of norm at most 2 (49 to 1149 per type):
  1. E has no terms below q^1;
  2. its q^1 part v is a Weil-invariant vector (checked exactly in Z[zeta_h]).
- Then E = Delta v exactly, and summing over the Leech line gives (24 + k) Delta.

All 23 certificates pass in 8 s. Exact enumeration of both sides through q^3 (all_holes.py) agrees independently.

**The swap theorem, explicitly.** v = sum_{d|h} (T(d) - T(h/d)) 1_{S_d}, where S_d = dZ x (h/d)Z runs over the
maximal isotropic subgroups. Here:
- T(1) = 0 and T(h) = 24 + k;
- for a prime p | h with p < h, T(p) is the number of components with p as a **bad prime** (equivalently, torsion
  prime): D_n gives 2, E6 and E7 give 2 and 3, E8 gives 2, 3 and 5, A_n gives none;
- T vanishes on all other divisors.

This holds for every divisor of every type (exhaustion).

| type | h | c_d |
|---|---|---|
| A1^24 | 2 | {1: -48, 2: 48} |
| A24 | 25 | {1: -25, 5: 0, 25: 25} |
| D4^6 | 6 | {1: -30, 2: 6, 3: -6, 6: 30} |
| E6^4 | 12 | {1: -28, 2: 4, 3: 4, 4: -4, 6: -4, 12: 28} |
| D16E8 | 30 | {1: -26, 2: 2, 3: 1, 5: 1, 6: -1, 10: -1, 15: -2, 30: 26} |
| E8^3 | 30 | {1: -27, 2: 3, 3: 3, 5: 3, 6: -3, 10: -3, 15: -3, 30: 27} |
| D24 | 46 | {1: -25, 2: 1, 23: -1, 46: 25} |

The full table is in swap_theorem_proof.out. Off the isotropic classes the coset theta series are exactly
swap-symmetric (swap_reciprocity.py). In root language:

#{roots r : r.rho = -1, r.z = -j} = #{roots r : r.z = -1, r.rho = -j}  for all j not divisible by h.

This is a reciprocity law for the roots of II_{25,1} between the Leech cusp and every Niemeier cusp.

**Representation theory alone does not force this.**
- For prime h the torus-invariant, sigma-antisymmetric weight-12 cusp space has dimension 1, so the theorem is
  automatic. A two-line proof via W(R) and the sum rule is given in THM-4758.
- For h = 4, 9 and 25 that space has dimension 2, 3 and 5 (swap_dimension.py). The theta function avoids the extra
  directions for arithmetic reasons, which the certificates verify.

**Divisor lattice.** Each S_d defines an intermediate Niemeier lattice L_d between Leech (d = 1) and N (d = h).
- Their Coxeter numbers are often d: A24 gives L_5 with h = 5, which looks like A4^6, and E8^3 gives 2, 3, 5, 6, 10.
- At odd d the D/E types drop by one: D24 gives L_23 with h = 22, D12^2 gives L_11 with h = 10, E8^3 gives L_15 with
  h = 14.
- Presumably these are the orbifolds of V_Lambda by powers of the deep-hole automorphism (Moller-Scheithauer).

## 2. h-cores, 3n+1 and 27 (corollaries; PROVED given Garvan-Kim-Stanton)

For pure A-types A_(h-1)^k:

theta(Lambda + c) = h^(k/2) (eta(tau)^h / eta(tau/h))^k + (24 + k) Delta.

Here eta^h(tau)/eta(tau/h) is the h-core generating function, so Leech points around a deep hole are counted by
k-tuples of h-cores. An h-core is the residue of a partition after all rim h-hooks are slid off the abacus: a
least-resistance reduction, and a natural recursion.

- **A24.** Every partition of fewer than 25 is a 25-core. So N_j = 5 p(j-1) for 1 <= j <= 24, and the first
  departure is N_25 = 5p(24) - 600 = 7275 (from Delta's -24 q^2).
- **A2^12 and "3n+1".** The number of 3-cores of n is a_3(n) = d_{1,3}(3n+1) - d_{2,3}(3n+1) = sum_{d | 3n+1} (d/3),
  i.e. Eisenstein ideals of norm 3n+1 (Granville-Ono). So
  N_j = 729 [x^(j-1)] (sum a_3(n) x^n)^12 + 36 tau(1 + j/3) [3 | j]. Check: N_3 = 729 (12 a_3(2) + 66 a_3(1)^2) - 864
  = 64746, which matches the enumeration. Thus Collatz's expression 3n+1 is literally the Eisenstein norm that counts
  Leech points around the A2^12 hole. 729 = 3^6 = 27^2 is the size of the ternary Golay code, its glue.
  - Exact but trivial bridge: since 2 is inert in Z[omega], a_3(n) > 0 forces v_2(3n+1) to be even. So an odd n with
    a 3-core is followed in Collatz by an even number of halvings.
  - There is no dynamical content (ANALOGY).
- **A8^3 and 27: the balanced hole.** For h = 9 = 3^2 and k = 3, the glue factor 9^(3/2) = 27 equals 24 + k = 27, so
  theta(Lambda + c) = 27 [(eta^9/eta(tau/9))^3 + Delta] = 27 eta^24 [1 + (eta(tau)/eta(tau/9))^3]. It is the only
  pure-A type with glue = vertex count.
- **A1^24.** 2-cores are staircases, so the counts are 4096 x (representations of j-1 as a sum of 24 triangular
  numbers), plus 48 tau.
- **{25, 1}.** The A24 defect lives only on the divisors {1, 25}; c_5 = 0.
- **{7, 21}.** The A6^4 hole (h = 7) has 28 = 7 + 21 vertices, and 28 = 3^3 + 1 is where the families 3^k + 1 and k^3 + 1 meet (k = 3).
  This is NUMEROLOGY: no map from the H-gaps.

## 3. The general cusp-pair swap (tested; REFUTED as stated)

The certificate applies to any two cusps z1, z2 of II_{25,1} (M = {z1, z2}^perp, D = (Z/m)^2 with m = -z1.z2).
Over all 276 pairs from Borcherds' 24 representatives, 275 pass with the line-matching swap. Examples:
- A1^24 x A24 (m = 13): 46 Delta, which is (600 - 48)/12, the prime-m sum rule;
- A2^12 x D24 (m = 16): 70 Delta;
- A24 x E8^3 (m = 6): 13 Delta;
- D4^6 x A6^4 (m = 3): 12 Delta;
- A1^24 x A2^12 (m = 2): 24 Delta.

**The one failure is A4^6 x E8^3 (m = 7).** M^* has 4 vectors of norm 12/7 in the class (3,2), i.e. 4 null vectors y
with (y.z1, y.z2) = (-3, -2), and none in its swap (2,3). Twisted swaps with lambda = 2, 5 are symmetric below norm 2
but do not match the lines. So the Leech cusp is special. Explaining why the Leech cusp always passes is open.

## 4. Prior work and literature

**Repo.**
- S52 HYP-9363 and THM-4752; the A24 cusp nu = (1/2^25 | 5/2) in tournament coordinates.
- Borcherds' table and thesis: z . rho = -h, rho = (rho_N, h, h+1), Venkov's rho_N^2 = 2h(h+1).
- The audit_S52 observation 5p(k-1), k <= 22, is now explained.

**Literature.** The background prior-art search covers Borcherds' thesis, Conway-Sloane, Moller-Scheithauer deep-hole
orbifolds, Shimada and t-core papers; it is recorded in section 6 once complete. Garvan-Kim-Stanton for t-cores and
Granville-Ono for 3-cores are cited for the corollaries.

## 5. Next

1. A conceptual proof of the swap theorem. A natural target is to specialise the fake-monster denominator identity to
   span(rho, z).
2. A conceptual proof of the bad-prime formula. Plausibly T(p) counts the cusps of Coxeter number p between rho and z.
3. Characterise the cusp pairs for which the swap holds.
4. Identify the intermediate lattices L_d by type.
