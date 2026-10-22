# The Leech lattice seen from its deep holes: a swap theorem, h-cores, divisors of 3n+1, bad primes, and umbral and Fricke moonshine (opus-2026-10-22-S53)

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
- For THM-4759 (section 6):
  - `principal_cores.py`: the Kostant/Kac identity component by component, and the glue lemma;
  - `coxeter_frame_twined.py`: Theta_N^[rho/h] = eta^24/eta_pi, all 23 types;
  - `twined_hauptmodul.py`: Leech twined traces;
  - `hauptmodul_cert.py`: lambencies, and the Atkin-Lehner component rule;
  - `multiple_cosets.py` and `deep_hole_multiples.py`: shallow multiples on both sides;
  - `exponent_law.py`: the exponent inequality;
  - `fricke_multiples.py`: Fricke invariance of the shallow-multiple twining;
  - `leech_h_plus_one.py`, `weyl_h_plus_one.py`, `fricke_h_plus_one.py`, `fricke_nudge.py`: the (h+1) Fricke points
    on both sides;
  - `swap_fricke_proof.py`: the Fricke swap certificates;
  - `swap_shell_test.py`: the next shell fails;
  - `fricke_neighbors.py`, `height_descent.py`: Fricke descent;
  - `order3_fricke_census.py`, `order4_fricke_census.py`: HYP-9368 evidence.

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
- Explained in THM-4759 (8): L_d is the orbifold of V_Lambda by g_c^(h/d). Its Coxeter number is d - 1 exactly when
  h/d is a common bad prime || h (Fricke descent, roots of height divisible by h/d), and d otherwise.

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
- **{25, 1}.** The A24 defect lives only on the divisors {1, 25}; c_5 = 0. Also 25.1^-1 is the Coxeter Frame shape of
  A24, the only one that is not a Conway-group Frame shape. Its lambency Gamma_0(25) is the only non-monstrous one,
  but its Fricke point of order 26 = 25 + 1 twines by T_26B (section 6).
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
but do not match the lines. So the Leech cusp is special.

Around the Leech cusp the swap holds for exactly the two nearest cusp types of each N: m = h (THM-4758) and m = h + 1
(the Fricke swap theorem, THM-4759 (11)). It fails for the next shell, m = h + 2 (swap_shell_test.py). Explaining
this conceptually is open.

## 4. Prior work and literature

**Repo.**
- S52: HYP-9363 and THM-4752, including the A24 cusp nu = (1/2^25 | 5/2) in tournament coordinates.
- Borcherds' table and thesis: z . rho = -h, rho = (rho_N, h, h+1), and Venkov's rho_N^2 = 2h(h+1).
- The audit_S52 observation 5p(k-1), k <= 22, is now explained.

**Literature (prior-art agent, 2026-10-22, 113 tool calls).** No published statement was found of:
- the swap theorem;
- the deep-hole theta identity in either the theta or the twined form;
- the A-type eta-quotient formula;
- the A24 counts 5p(j-1).

Sources searched include:
- Borcherds' thesis (Lemma 2.2.7: |rho/h - y|^2 >= 2(1 + 1/h); section 2.7: the holy construction);
- Moller-Scheithauer, Ann. Math. 2023 (Prop. 6.2: Lambda^d with (Z/h)^2 and its two isotropic subgroups; Prop. 6.3:
  deep holes are extremal; weight-1 data only);
- van Ekeren-Lam-Moller-Shimakura;
- Chigira-Lam-Miyamoto (the V_N orbifold by e^{2 pi i rho(0)/h});
- Lam-Shimakura;
- Hohn-Moller;
- Shimada;
- Israel-Melnikov-Proto;
- Cheng-Duncan-Harvey (CDH).

Published numbers agree with the identity, but none of these sources states or uses it:
- OEIS A004034 (the A1^24 hole, from SPLAG p. 511, no formula);
- Dong, arXiv 2608.15431 (A3^8 twined sums computed by dynamic programming: T_4 = -25,520, T_12, T_16);
- Bachoc-Moustrou-Vallentin-Zimmermann, Forum Math. Sigma 2026, Table 3: second radii 2 + 2/h. Their E8^3 entry
  sqrt(32/15) looks like a misprint for sqrt(31/15).

Gaps: SPLAG Ch. 24 p. 511 and Borcherds' thesis Ch. 5 could not be read.

Classical inputs are cited: Garvan-Kim-Stanton (t-cores), Granville-Ono (3-cores), Kostant (1959; principal element
~ Coxeter element), and CDH Prop. 2.1 (T^X is a principal modulus; Remark 2.2 asks for a conceptual proof).

## 5. Next

1. A conceptual proof of the swap theorem. A natural target is to specialise the fake-monster denominator identity to
   span(rho, z).
2. Bad primes. Partly done (section 6): the lower bounds come from the exponent law, but the T(d) values off the
   shallow lines are still case-checked.
3. Characterise the cusp pairs for which the swap holds.
4. Intermediate lattices L_d: DONE (THM-4759 (8)).
5. HYP-9368: genus-zero inner automorphisms of V_Lambda. The prime orders are sketched and orders 3 and 4 are sampled;
   the composite orders are open. They need dim V'_1 = 24e at a Fricke point, or "every Fricke point is a nudged deep
   hole".
6. Prior art for the (h+1)-construction and Fricke twining: a second search is running (see the close-out letter).

## 6. The umbral picture (THM-4759) and Fricke twining

**Coxeter Frame shapes do it for every type.** Let pi^X be the Frame shape of the Coxeter element, i.e.
det(x - c) = prod (x^a - 1)^(m_a), taken over the components. Then:
- **theta(N + rho_N/h) = |N/Q| eta^24 eta_{pi^X}(tau/h).** This is Kostant's conjugacy (principal element ~ Coxeter
  element) applied to V_Q, plus a glue lemma (the special automorphisms fix rho/h). principal_cores.py checks it
  component by component through q^29.
- **Tr_{V_N}(sigma_N q^(L_0 - 1)) = T^X = 1/eta_{pi^X},** the principal modulus of Cheng-Duncan-Harvey's umbral
  moonshine (their Prop. 2.1).
- Combined with THM-4758: **Tr_{V_Lambda}(g_c q^(L_0 - 1)) = T^X + 24 + k,** and
  theta(Lambda + c) = Delta (24 + k + |N/Q| eta_{pi^X}(tau/h)).
- The two traces agree in every weight except weight 1.
- For X != A24 (CDH: pi^X is a Co_0 Frame shape and Gamma^X is monstrous), Carnahan's orbifold duality adds two more
  traces: a Conway lift with trace T^X, and a Monster element of V^nat with trace T^X + k. So four traces agree in all
  weights >= 2.
- A24, with Frame shape 25^1 1^-1, is the unique exception on the deep-hole side: that is the {25, 1}. Its Fricke
  point (order 26) is monstrous, 26B.

**Bad primes = Atkin-Lehner = shallow multiples.** The lambency of X is h + {h/p : p a bad prime of every component,
p || h} (D_n with n even: h + h/2; E6: 12 + 4; E8: 30 + 6, 10, 15). For d | h, the coset N + d rho_N/h:
- has minimal norm 2 - 2d/h, with a unique minimal vector, exactly at those primes;
- has minimal norm >= 2 otherwise.

This is proved from the **exponent law**: sum_i (m_i mod e)(e - (m_i mod e)) >= l e (e - 1)/6. Equality holds iff e = 1
or e = h/p with p bad and p || h. The proof is by hand for A and D; exponent_law.py checks 2815 cases, including E.
The law is the vacuum anomaly of c^(h/e); its case e = h is the strange formula. It is special to ADE: it fails for
H_3, H_4 and I_2(m >= 7), and has no bad-prime equality cases for B_n and C_n.

On the Leech side every multiple jc is one of two things:
- a deep hole of an A-containing type (39 instances);
- at gcd(j, h) = d a common bad prime with d || h, an interior point at squared distance 2 - 2d/h from a unique Leech
  point. Then Lambda + dc also has exactly 24 vectors of norm 2.

**Fricke twining (new, fricke_multiples.py).** For the shallow multiple v = dc (order e = h/d), the trace
Tr_{V_Lambda}(e^{2 pi i v(0)} q^(L_0 - 1)) is invariant under tau -> -1/(e tau). Hence it equals **24 + T_{e+}**, the
normalised Hauptmodul of Gamma_0(e)+:
- D4^6 gives q^-1 + 24 + 783q + 8672q^2 = 24 + T_3A;
- D24 gives q^-1 + 24 + 4q + 7q^2 = 24 + T_23A.

The proof is short. F = Theta^[v]/eta^24 has simple poles only at infinity and 0, because the proper multiples of v
are deep holes. F o W_e has the same principal parts, because the shallow vector is unique. So F - F o W_e is a
holomorphic function on X_0(e) that is odd under W_e, hence 0.

So one deep hole carries both kinds of moonshine: non-Fricke through its Niemeier lattice, and Fricke through its
bad-prime multiples. Level 2 is impossible: an order-2 point v is congruent to -v mod Lambda, so its minimal vectors
come in pairs.

**Every Niemeier lattice has a Fricke point (THM-4759 (10)).** Recall rho = (rho_N, h, h+1) in II_{25,1} = N + U. The
second cusp z'' of type N adjacent to the Leech cusp has -z''.rho = h + 1. Its Leech-side point c'' (order h + 1)
satisfies, for all 23 types:
- Lambda + c'' has a unique vector of norm 2 - 2/(h+1) and exactly 24 vectors of norm 2;
- its proper multiples are deep holes;
- Tr_{V_Lambda}(e^{2 pi i c''(0)} q^(L_0 - 1)) = 24 + T_{(h+1)+}, the Fricke Hauptmodul. For example D24 gives
  order 47 and q^-1 + q + 2q^2 (47AB); E8^3 and D16E8 give 31AB; A24 gives 26+26 = 26B. OEIS confirms
  3A, 4A, 5A, 6B, 7A, 8A, 9A, 10D, 11A, 13A, 14C, 15C, 17A, 19A, 23A and 26B.

On the N side, rho_N/(h+1) is the unique minimum of its coset, and the (h+1)-neighbour of N there is Lambda: a second
holy construction. The prime values of h + 1 are 3, 5, 7, 11, 13, 17, 19, 23, 31 and 47, all supersingular. This is
forced: the twining gives a degree-1 function on X_0(p)/W_p.

**The eleven bad-prime multiples are the (h+1)-points of the descents.** For example, D4^6's 2c is A1^24's
(h+1)-point, which is the M24-fixed point (2^24)/(3 sqrt 8).

The Leech orbifold at pc is the Niemeier lattice of the **Borel-de Siebenthal descent** X_p. Its roots are the roots of
N whose height is divisible by p (height_descent.py); it is proved for prime h/p:
- D24 -> D12^2 -> D6^4 -> A3^8 (D-index halving);
- D16E8 -> D8^3 -> D4^6 -> A1^24;
- D10E7^2 -> A7^2D5^2;
- E6^4 -> A2^12;
- E8^3 -> D8^3, A8^3 and A4^6.

**The 3n+1 the prompt asked about.** With H = h + 1, descent at p is H -> (H-1)/p, i.e. ascent H -> pH + 1:
- p = 3 is literally 3n + 1: A2^12 (H = 4 = 3^1 + 1) -> E6^4 (13), and A8^3 (H = 10 = 3^2 + 1) -> E8^3 (31);
- p = 2 gives the Cunningham chain 5 -> 11 -> 23 -> 47 and the Mersenne chain 3 -> 7 -> 15 -> 31;
- p = 5 gives 6 -> 31.

So 31 (E8^3) has three ascent predecessors: 31 = 2·15 + 1 = 3·10 + 1 = 5·6 + 1. The affine maps are exact; reading
them as Collatz dynamics is ANALOGY.

**27.** The 9+ Hauptmodul of the D10E7^2 Fricke point is x + 3 + 27/x with x = (eta/eta(9 tau))^3. A8^3 is the
balanced hole with glue 27 = 24 + k.

**Collatz remarks (NUMEROLOGY: exact arithmetic, no lattice mechanism known).**
- The accelerated Collatz map T(H) = (3H+1)/2 sends the Mersenne Fricke chain 3, 7, 15, 31 (A1^24, D4^6, D8^3,
  D16E8/E8^3) term by term to the Cunningham chain 5, 11, 23, 47 (A3^8, D6^4, D12^2, D24). Indeed
  T(2^k - 1) = 3·2^(k-1) - 1, which on D-indices is n -> 3n/2.
- The Collatz trajectory of 27 (111 steps, peak 9232) runs 27 -> 82 -> 41 -> 124 -> 62 -> 31 -> 94 -> 47 -> 142 -> 71.
  It contains the Niemeier H-values 31, 47, 23, 10, 5, 8, 4, and the supersingular primes 41, 31, 47, 71, 23, 5, 2.
- The Collatz orbit of 3 (3, 10, 5, 16, 8, 4, 2, 1), with the doublings 6, 12, 24, consists of D-indices of Niemeier
  components. Every halving step in it is a Fricke descent; the odd steps 3 -> 10 and 5 -> 16 have no lattice
  meaning.
- The genuine 3n+1 is the Fricke ascent at the bad prime 3, H -> 3H + 1 (A2^12 -> E6^4, A8^3 -> E8^3), together with
  the Eisenstein norm 3n + 1 of 3-cores at the A2^12 hole.

HYP-9368 conjectures that the genus-zero inner automorphisms of V_Lambda are exactly 0, the 23 deep holes and the 23
Fricke points. A first draft that allowed only the eleven bad-prime multiples was refuted by the (h+1)-points before
registration.

**The Fricke swap theorem (THM-4759 (11), PROVED by the same certificate).**
theta(Lambda + c'') - theta(N + rho_N/(h+1)) = 24 Delta, with the pure defect vector 24(1_{S_{h+1}} - 1_{S_1}). Hence
V_N twined by e^{2 pi i rho(0)/(h+1)} is exactly T_{(h+1)+}. The deep-hole defect 24 + k counts the nodes of the
extended Dynkin diagrams; the Fricke defect 24 counts the finite simple roots.

**Geometry (fricke_nudge.py).** z'' = z + r is the reflection of the deep-hole cusp in an adjacent Leech wall. On Λ,
c'' = c + (lambda - c)/(h+1), the deep hole moved 1/(h+1) of the way towards a second-shell point lambda. Exactly 24
hole vertices stay at sqrt 2, and the k affine nodes move away.

**Levels.** The deep hole is the level-0 Weyl point rho/h (Kostant's principal element). The Fricke point is the
level-1 Weyl point (rho^ + Lambda_0)/(1 + h), of Kac type (2;1,...,1), where the level-1 WZW S-matrix lives. The
swap holds at levels 0 and 1 and fails at the next shell.
