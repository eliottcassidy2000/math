# Selfie tournaments, arc Hamiltonian-path counts, shaving, and the odd Mallows-Sloane partner

**Lane:** procgen selfie lane, 2026-10-01 (owner prompt of 2026-10-01). **Status:** COMPLETE for this pass.
**Promoted:** [THM-4524](../../01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md) (PROVED + INDEPENDENTLY AUDITED, orchestrator 2026-10-01; audit `procgen_selfie_20261001_orchestrator_check.out`); C2 + C5 -> [HYP-9167](../hypotheses/HYP-9167-paley-minus-a-vertex-is-redei-rigid.md), D1 -> [HYP-9168](../hypotheses/HYP-9168-hp-blocking-number-equals-hall-bound.md); OPEN-Q-060 annotated ANSWERED (Burnside form).
**Labels used:** PROVED, CITED, FINITE-EXACT, VERIFIED, EMPIRICAL, CONJECTURE (= OPEN, with evidence), REFUTED.
**Code:** `04-computation/experiments/procgen_selfie_20261001_run.py` (runner, ends ALL CHECKS PASSED),
helpers `procgen_selfie_20261001_lib.py`, C engines `procgen_selfie_20261001_arcs.c` (arc-HP counts,
Aut, shaving, branch and bound), `procgen_selfie_20261001_euler.c` (Euler graphs + twist character),
`procgen_selfie_20261001_paley2.c` (mod-2 DP for Paley minus a vertex up to N = 26).
**Output:** `05-knowledge/results/procgen_selfie_20261001.out` (runner stdout);
`05-knowledge/results/procgen_selfie_20261001_n10_census.out` (exhaustive N = 10 class census, same engine);
`05-knowledge/results/procgen_selfie_20261001_n10_oddfilter.out` (independent mod-2 N = 10 census, few-odd witnesses).

Notation: T a tournament on V, |V| = N; H(T) = #directed Hamiltonian paths (HPs); for an arc e,
c(e) = c_T(e) = #HPs through e; start(v), end(v) = #HPs starting / ending at v; switch_L(T) reverses
every arc with exactly one end in L; tiling model as in THM-022/THM-474 (base path N -> N-1 -> ... -> 1).

## 0. Headlines

1. **Selfie gauge theorem (PROVED, VERIFIED N <= 7).** (tiling t, loop set L) -> switch_L(T_t) is exactly
   2-to-1 onto labeled tournaments (fibre {L, L^c}); the N loops are the switching (gauge) bits and the
   N-1 base-path arc orientations are their discrete derivative: arc (i+1 -> i) is reversed iff exactly
   one of i, i+1 is looped. C(N-1,2) + N = C(N,2) + 1; the +1 is the global complement (H^0 of K_N).
2. **Loop-Walsh degree bound (PROVED).** As a function of the loop spins, H(switch_L T) has Walsh
   support on even sets A with |A| <= 2 floor(N/3); hence sum_L (-1)^|L| H(switch_L T) = 0 for every T,
   the loop enumerator is divisible by (1+y)^(N - 2 floor(N/3)), and for N <= 5 H is an Ising
   Hamiltonian on the switching class with explicit couplings.
3. **Fixed-point OCF (PROVED).** sum_{S indep in Omega(T)} 2^|S| prod_{v not in V(S)} x_v
   = sum_{U subset V} H(T[U]) prod_{v not in U} (x_v - 1). Loops as odd 1-cycles: a selfie with loop set L
   has I(Omega + loops, 2) = sum_{M subset L} 2^|M| H(T - M), always odd and = H(T) + 2|L| (mod 4).
4. **The owner's conjecture is REFUTED in both natural readings**, and the true "5/6" threshold is a
   parity phenomenon: **N = 6 is the smallest order >= 3 with a tournament in which every arc lies on an
   odd number of HPs** (FINITE-EXACT; unique class, = Paley QR_7 minus a vertex, one of the two n = 6 H-maximizers).
   All-odd needs N = 1, 2 (mod 4) (PROVED); none at N = 5, 9 (exhaustive); exactly 2 classes at N = 10
   (exhaustive census of all 9733056 classes); present at N = 6, 10, 18, 22, 26 via Paley QR_q minus a vertex,
   q = 7, 11, 19, 23, 27 (CONJECTURE for all q = 3 mod 4, antipodal arcs PROVED).
   All-even needs N odd (PROVED) and holds for every Cayley tournament of an odd abelian group (PROVED).
   A pattern seen at N = 4, 6, 8 (even N => at least N - 1 odd arcs) is REFUTED at N = 10 (minimum 5).
5. **Shaving.** Single-arc parity-preserving shaving = arcs with c(e) even (PROVED); the parity-break
   number rho(T) is 1, or 2 exactly for all-even T (PROVED). The HP-blocking number beta(T) (fewest arcs
   whose deletion kills every HP) equals a Hall-deficiency bound in every class N <= 9 (CONJECTURE,
   FINITE-EXACT N <= 9); the naive two-source/two-sink bound fails first at N = 7 (explicit).
6. **OPEN-Q-060 answered (PROVED + VERIFIED n <= 10):** A049313(n) = number of unlabeled Euler graphs on
   n vertices all of whose automorphisms reverse an even number of edges of a fixed orientation
   (equivalently: no automorphism has an odd number of even cycles whose antipodal pairs are edges).
   Twisted Mallows-Sloane / twisted Brauer lemma on the tournament torsor.

## 1. (A) The gauge theorem for selfie tournaments

**Definition.** A selfie tournament of size N is a pair (t, L): t in {0,1}^(C(N-1,2)) a tiling (tile bits)
and L subset V the set of looped vertices. Its underlying labeled tournament is switch_L(T_t).

**Theorem A1 (PROVED).** Phi(t, L) = switch_L(T_t) is exactly 2-to-1 onto the 2^C(N,2) labeled tournaments,
with fibres {(t, L), (t, V\L)}. The base-path arc (i+1 -> i) is reversed in Phi(t, L) iff exactly one of
i, i+1 lies in L.
*Proof.* Surjective by THM-474 (every switching class contains a tiling tournament T_t). If
switch_L(T_t) = switch_L'(T_t'), then T_t, T_t' are switching equivalent, so t = t' (THM-474 uniqueness),
and switch_{L xor L'} fixes T_t; a switching by W reverses every arc of the cut (W, V\W), nonempty unless
W in {empty, V}. The derivative statement: arc (i+1, i) is in the cut of L iff [i in L] != [i+1 in L]. QED
*Verified:* N = 2..7 (all 2^C(N-1,2) x 2^N selfies; image multiplicities all 2; fibre = complement;
derivative law on every base-path arc).

**Loops versus path edges (PROVED).** The N loop bits are a 0-cochain l in C^0(P_N; F_2) on the base path
P_N; the N-1 path-arc orientations are its coboundary delta(l) in C^1(P_N; F_2). delta is onto with
kernel the constants: the loops carry exactly one more bit than the path arcs, the global complement
(H^0 = F_2), which acts trivially on the tournament. So the "must exist" path edges are the gauge-fixing
condition (tree/axial gauge on the spanning path), the optional loops are the gauge frame, and the
tiles are the gauge-invariant data (the holonomy around the fundamental cycle of each non-path edge).
Parity corollary: #reversed path arcs = [N looped] + [1 looped] (mod 2).
Equivalently: the tiling model fixes the vertex sequence N, N-1, ..., 1 to be a DIRECTED Hamiltonian path;
the selfie model lets it be an ORIENTED Hamiltonian path of arbitrary type tau = delta(l) in {+,-}^(N-1)
(the objects of Forcade / Grunbaum / OPEN-Q-059), the loops recording the type up to global flip. Each
type occurs for exactly 2 * 2^C(N-1,2) selfies.

**Theorem A2 (PROVED; class-sum CITED THM-1467).** For every tournament T and every vertex ordering pi,
exactly two loop sets (L and V\L) make pi a directed HP of switch_L(T) (the THM-474 argument applied to
the spanning path pi). Hence, for every oriented path type tau in {+,-}^(N-1),
sum_{L subset V} N_tau(switch_L T) = 2 N!, in particular sum_L H(switch_L T) = 2 N! (the switching-class
sum N!, THM-1467 / `switching_class_Hsum_deathstar_S61b.out`), and the forward-edge distribution sums to
sum_L a_k(switch_L T) = 2 N! C(N-1, k). Consequently the switching-class sum of the deformation in the
deformed Eulerian numbers (THM-062) is class-independent:
sum_{T' in [T]} (a_k(T') - A(N,k)) = N! C(N-1,k) - 2^(N-1) A(N,k)
(N = 3: 2,-4,2; N = 4: 16,-16,-16,16; N = 5: 104,64,-336,64,104). *Verified* N <= 5 (every T, every pi).

**Theorem A3 (loop-Walsh expansion; PROVED).** Write sigma_v = (-1)^[v in L] and s_i(pi) = +1/-1 as
pi_i -> pi_(i+1) is / is not an arc of T. Then
  H(switch_L T) = sum_{A subset V, |A| even} h_T(A) prod_{v in A} sigma_v,
  h_T(A) = 2^(1-N) sum_{pi} prod_{i in I_A(pi)} s_i(pi),
where I_A(pi) is the union of the position segments between the 1st and 2nd, 3rd and 4th, ... elements of
A in pi-order.
*Proof.* H(switch_L T) = sum_pi prod_i (1 + s_i sigma_{pi_i} sigma_{pi_(i+1)})/2; expand; on a path the
set of steps with prescribed boundary A is unique. QED
(a) h_T(empty) = N!/2^(N-1).
(b) **Degree bound: h_T(A) = 0 whenever |A| > 2N/3.** *Proof.* Reversing the first odd-length segment of
pi (as a block) keeps all segment positions, multiplies that segment's product by -1 and leaves the others
unchanged: a sign-reversing involution. Surviving pi have every segment of even length >= 2, hence with at
least one interior vertex outside A, all distinct: N - |A| >= |A|/2. QED
(c) **Corollaries (PROVED).** sum_L (-1)^|L| H(switch_L T) = 0 for every T (N odd: L <-> V\L; N even: A = V
violates (b)). The loop enumerator Lambda_T(y) = sum_L y^|L| H(switch_L T)
= sum_A h_T(A) (1-y)^|A| (1+y)^(N-|A|) is divisible by (1+y)^(N - 2 floor(N/3)).
(d) **Ising form for N <= 5 (PROVED):** H(switch_L T) = N!/2^(N-1) + sum_{a<b} J_ab sigma_a sigma_b with
J_ab = 2^(2-N) sum_{k odd} (N-k-1)! sum_{a w_1 ... w_k b} prod s (signed odd-interior walks a -> b).
*Verified* (runner): segment formula for every coefficient (N = 3, 4 all T; N = 5 all tilings), couplings,
odd-degree vanishing, degree bound and alternating-sum identity for all tilings N <= 7. The top degree
2 floor(N/3) is attained by every tiling for N = 3, 5, 6, 7 and by 6 of the 8 tilings at N = 4 (FINITE-EXACT);
the two exceptional N = 4 switching classes carry H = 3 on all eight members (loop enumerator 3(1+y)^4).
(e) **Constant-H classes.** If H is constant on a switching class it equals N!/2^(N-1) (class sum N!), which
must be an odd integer; since v_2(N!) = N - s_2(N), this forces N to be a power of 2 (PROVED). It happens at
N = 1, 2, 4 (2 labeled classes at N = 4) and NOT at N = 8: none of the 30 iso classes with H = 315 = 8!/2^7
has H constant on its switching class (FINITE-EXACT).

## 2. (B) Loops as odd 1-cycles

**Theorem B1 (fixed-point OCF; PROVED).** For indeterminates x_v,
  sum_{S indep in Omega(T)} 2^|S| prod_{v not in V(S)} x_v = sum_{U subset V} H(T[U]) prod_{v not in U} (x_v - 1).
*Proof.* Apply the OCF to every T[U] (odd cycles of T[U] = odd cycles of T inside U), then
sum_U prod_{v not in U}(x_v - 1) sum_{S: V(S) subset U} 2^|S| = sum_S 2^|S| prod_{v not in V(S)} (1 + (x_v - 1)). QED
Univariate H(T;x) = sum_S 2^|S| x^(N-|V(S)|) = sum_U (x-1)^(N-|U|) H(T[U]). In Grinberg-Stanley normalisation
U_T = sum_sigma 2^psi(sigma) p_type(sigma), H(T;x) is U_T at p_1 -> x, p_odd>=3 -> 1 (H(T) = H(T;1) = U_T(1,0,0,..));
the bridge's phrasing "p_odd>=3 -> 2" presupposes U_T written without the 2^psi weight (the 2 is counted once).
Check: U_T(1,1) = sum over ordered 2-block covers of H(T[B1]) H(T[B2]) = sum_S 4^|S| 2^(N-|V(S)|) (VERIFIED N <= 7).
**Meaning (PROVED):** for integer x >= 1, H(T;x) = #{(P, c)}: P a directed path of T (possibly empty), c a
colouring of the off-path vertices by x-1 colours. x = 1: HPs; x = 2: all paths including the empty one;
x = 3: every vertex looped, each off-path vertex visited by its loop in one of 2 ways
(= I(Omega(T) + N loops, 2)); x = 0: sum over partitions of V into odd cycles of 2^#cycles, so
sum_U (-1)^(N-|U|) H(T[U]) = 0 for N in {1,2,4} and = 2 hc(T) (Hamiltonian cycles) for N in {3,5,7};
x = -1: H(T;-1) = (-1)^N I(Omega(T), -2).
**Selfie Redei (PROVED).** For a loop set L, Z_L(T) = I(Omega(T) + loops on L, 2) = sum_{M subset L} 2^|M| H(T - M)
is odd and Z_L(T) = H(T) + 2|L| (mod 4) (each H(T - v) is odd). Meaning: HPs of the selfie in which a looped
vertex may be skipped (visited by its loop, 2 ways). *Verified* N <= 7 for all classes and all L.
Value tables: see .out (H(T;3), I(Omega,-2) ranges and signs, H(T;0)).

## 3. (C) The owner's conjecture

Readings. (i) "For N >= 6 every class has an arc on no HP." (ii) Parity: "classes in which every arc is
on an odd number of HPs exist up to N = 5, and above that every class has an arc on an even number".

**(i) REFUTED.** Tournaments with every arc on some HP exist for every N >= 3 (PROVED: a regular tournament,
N odd, every arc on a Hamiltonian cycle by Alspach 1967 [CITED]; for N even, a regular R_(N-1) followed by a
sink, using the dead-arc formula below). Exhaustive labeled counts (N = 3..9): 2/8, 40/64, 664/1024,
26048/32768, 1934528/2097152, 261710080/268435456, 68219344384/68719476736 (fractions .25 .625 .648 .795
.922 .975 .993, and N = 10: 35109910971392/2^45 = .9979, EMPIRICAL -> 1).
**Dead-arc structure (PROVED).** If C_1 => C_2 => ... => C_k are the strong components, every HP traverses
them in order, every vertex of a strong component of size >= 3 starts and ends a HP of it (Camion), so
z(T) := #arcs on no HP = sum_{j >= i+2} |C_i||C_j| + sum_i z(T[C_i]). In particular all arcs lie on HPs iff
k <= 2 and each component has this property. (Formula VERIFIED on all classes N <= 7.) Strong tournaments
with a dead arc: N = 5: 1 class, 6: 2, 7: 9, 8: 33, 9: 167, 10: 1072.

**No tournament on n >= 3 vertices has every arc on some HP and some arc on every HP (PROVED;
orchestrator's observation for n <= 6).** If c(u->v) = H then no other out-arc of u and no other in-arc of v
is on a HP, so (all arcs being on HPs) u -> v is u's only out-arc and v's only in-arc; then every other w
satisfies w -> u and v -> w, v is a source of T - u, and appending u to any HP of T - u (which starts at v
and ends at some w -> u) gives a HP of T ending at u, i.e. not containing u -> v. Contradiction.

**(ii) REFUTED, and the true threshold.** Parity constraints (PROVED): sum_e c(e) = (N-1)H with H odd,
sum_{w: v->w} c(v->w) = H - end(v), sum_{u: u->v} c(u->v) = H - start(v). Hence
* all arcs odd  =>  N = 1, 2 (mod 4) and end(v) = 1 + outdeg(v), start(v) = 1 + indeg(v) (mod 2);
* all arcs even =>  N odd and every start(v), end(v) odd;
* #odd arcs = N - 1 (mod 2).

| N | classes | all arcs on a HP | all-odd (cls/lab) | all-even (cls/lab) | all-equal | min-max #odd |
|---|---------|------------------|-------------------|--------------------|-----------|--------------|
| 3 | 2 | 1 | 0 | 1 / 2 | 1 (C3) | 0-2 |
| 4 | 4 | 3 | 0 (parity) | 0 (parity) | 0 | 3-5 |
| 5 | 12 | 7 | **0** | 3 / 184 | 1 (RT5) | 0-8 |
| 6 | 56 | 44 | **1 / 240** | 0 (parity) | 0 | 5-15 |
| 7 | 456 | 412 | 0 (parity) | 28 / 96608 | 1 (QR7) | 0-18 |
| 8 | 6880 | 6674 | 0 (parity) | 0 (parity) | 0 | 7-27 |
| 9 | 191536 | 189992 | **0** | 899 / 285894016 | 0 | 0-32 |
| 10 | 9733056 | 9711588 | **2 / 4354560** | 0 (parity) | 0 | **5**-45 |

So the owner's intuition is mirrored: there is NO all-odd tournament for 3 <= N <= 5, and N = 6 is the first
(FINITE-EXACT): the unique class is Paley QR_7 minus a vertex (H = 45 = the n = 6 maximum, |Aut| = 3,
c in {13 (12 arcs), 23 (3 antipodal arcs)}). This is a "Redei-rigid" tournament: no single arc can be
shaved without breaking Redei's parity.

**Theorem C1 (Cayley all-even; PROVED).** Every Cayley tournament on an abelian group of odd order (e.g.
every circulant, every Paley tournament) has every arc on an even number of HPs.
*Proof.* For an arc u -> v, phi(x) = (u+v) - x is an involutive anti-automorphism (T -> T^op) swapping u, v;
P -> reverse(phi(P)) is an involution on the HPs through u -> v; a fixed path would have the self-paired arc
u -> v in the middle position, impossible for N odd. QED (VERIFIED: all circulants N = 3..13, all Cayley
tournaments of Z_3 x Z_3.)

**Conjecture C2 (Paley minus a vertex is Redei-rigid; CONJECTURE).** For every prime power q = 3 (mod 4),
QR_q minus a vertex (N = q - 1 = 2 mod 4) has every arc on an odd number of HPs.
VERIFIED q = 7, 11, 19 (exact counts) and q = 7, 11, 19, 23, 27 (mod-2 DP, N up to 26, all arcs).
**Partial proof (PROVED):** F_q^* acts on HP(QR_q - 0) (multiplication by residues; by non-residues
followed by reversal); arc orbits have size q - 1 except the antipodal orbit {u, -u} of odd size (q-1)/2;
since sum_e c(e) = (q-2)H is odd, the antipodal arcs are odd. Also (PROVED, arc-transitivity of QR_q):
end(u) is odd for u a non-residue (it equals the number of HPs of QR_q with last arc u -> 0, which is
2H(QR_q)/(q(q-1)), odd) and even for u a residue.
Controls: the skew-Hadamard-doubled DRT(15) minus any vertex is not all-odd (#odd in {39, 45, 49} of 91);
non-Paley unions of cubic-residue cosets on Z_19 minus a vertex: 117/153 odd. At N = 10 there are exactly
TWO all-odd classes (exhaustive census): QR11 - v (H = 15745, |Aut| = 5) and a rigid class with H = 3929
(|Aut| = 1, T = 111111001110111111111101110111111111110111110 in gentourng form).

**Observation C3 (EMPIRICAL).** For every n <= 11 with n != 0 (mod 4) some H-maximizer has uniform arc
parity: n = 3, 5, 7, 9, 11 all-even (C3, both n = 5 maximizers, QR7 (unique), the unique n = 9 maximizer,
QR11 with H = 95095, Cayley), n = 6, 10 all-odd (QR7 - v, QR11 - v attaining the known maxima 45, 15745
[A038375, CITED; maximality at n = 10, 11 not re-derived here]). Caveat: n = 6 has TWO maximizing classes
(H = 45); the other one has only 9 of 15 arcs odd. At n = 4, 8 uniform parity is impossible
(#odd = 5/6 at n = 4; 11, 15, 19 of 28 over the six n = 8 maximizers).

**C4 (even N): REFUTED at N = 10.** For N = 4, 6, 8 every tournament has at least N - 1 arcs on an odd
number of HPs, with equality e.g. for transitive and for a source added to an all-even tournament (FINITE-EXACT),
which suggested "even N => at least N - 1 odd arcs". The exhaustive N = 10 census refutes it: the minimum is
**5** (4 classes, each |Aut| = 1; witnesses in `procgen_selfie_20261001_n10_oddfilter.out`, found by an
independent mod-2 DP that reproduces the census histogram exactly, and re-verified with exact counts);
99 classes have 7. So the "N - 1" pattern is a small-N artefact. Still PROVED: #odd arcs = N - 1 (mod 2), so
>= 1 for even N. (The odd-arc graph is connected for N = 4, 6, not always for N >= 8; the GF(2) rank of the
odd-arc matrix can be 1.)

**Corrected statements behind the owner's 5/6.** (1) Dead arcs: no threshold; every N >= 3 has HP-covered
tournaments and their proportion tends to 1 (EMPIRICAL). (2) Parity: the threshold is real and reversed:
"every arc on an odd number of HPs" is impossible for 3 <= N <= 5 and first occurs at N = 6, then at N = 10
(and conjecturally at every N = q - 1, q = 3 mod 4 a prime power); it never occurs at N = 0, 3 (mod 4),
nor at N = 9.

## 4. (D) Shaved tournaments

* **HP-core** K(T) = arcs with c(e) > 0: the smallest spanning subdigraph with the same HP set (PROVED);
  size = C(N,2) - z(T) given by the dead-arc formula above.
* **Redei shaving.** H(T - e) = H(T) - c(e), so a single arc can be shaved keeping H odd iff c(e) is even
  (PROVED). All-odd tournaments are exactly the ones where no single shave keeps parity. Shaving arc by arc
  down to a single HP with H odd at every step is possible for every class N <= 7 except the all-odd one
  (FINITE-EXACT).
* **Parity-break number** rho(T) = min #arcs whose deletion makes H even: rho = 1 iff some c(e) is odd;
  **rho = 2 for every all-even T (PROVED)**: end(v) is odd, so some in-arc e = u -> v has an odd number of
  HPs ending with it; then sum_w c(u->v->w) = c(e) - (that number) is odd, so some pair e, f = v -> w has
  c(e,f) odd and H(T - e - f) = H - c(e) - c(f) + c(e,f) is even. (VERIFIED N <= 7.)
* **HP-blocking number** beta(T) = min #arcs whose deletion leaves no HP. Two sources / two sinks give
  beta <= sigma(T) = min(d-(u)+d-(v), d+(u)+d+(v)) over pairs; beta = sigma for all classes N <= 6, but
  **fails at N = 7** for T = 111000101111111111101 = (0 -> A => B -> 0, A and B cyclic triangles):
  sigma = 4, yet deleting the three arcs of A leaves the three A-vertices with the single common
  in-neighbour 0, so no HP (beta = 3). General **Hall-deficiency bound**
  hall(T) = min over X (|X| >= 2), Y (|Y| <= |X| - 2) of #arcs entering X from outside Y (and the dual with
  out-arcs): deleting them gives the predecessor (successor) bipartite graph deficiency >= 2, so no HP
  (beta <= hall, PROVED). **Conjecture D1:** beta(T) = hall(T) for every tournament. FINITE-EXACT for all
  classes N <= 9 (exact subset tables N <= 7; branch and bound N = 8: 6880 classes, 13 with hall < sigma;
  N = 9: 191536 classes, 180 with hall < sigma, 8.8e8 search nodes). Max beta = N - 1 for odd N (attained
  by the regular classes: 3 at N = 7, 15 at N = 9) and N - 2 for even N, N <= 9: every regular tournament
  survives the deletion of any N - 2 arcs (CONJECTURE beyond N = 9).

## 5. (E) A named open problem: OPEN-Q-060 (the odd Mallows-Sloane partner)

**Theorem E1 (PROVED).** For n >= 1, A049313(n) (switching classes of tournaments up to isomorphism)
equals the number of unlabeled Euler graphs F on n vertices such that every automorphism g of F reverses
an even number of edges of a fixed orientation O of F, i.e.
eps_F(g) = (-1)^#{edges of F whose O-orientation g carries against O} = +1 for all g in Aut(F).
eps_F does not depend on O; in cycle form eps_F(g) = (-1)^(number of even cycles of g of length 2k whose
antipodal pairs {x, g^k x} are edges of F).
*Proof.* Labeled switching classes form the set A = F_2^E / C (C = cut space) and S_n acts affinely,
T_g(x) = g x + c_g, with cocycle c_g = flip vector of g(T_0) relative to a reference tournament T_0
(THM-479 setup). The dual of A is C^perp = cycle space = Euler graphs, with the plain permutation action.
For one g: #Fix(T_g) = |A^g| if c_g in Im(1+g), else 0; |A^g| = |(A^)^g| (Brauer's permutation lemma) and
c_g in Im(1+g) iff chi(c_g) = 1 for all g-fixed characters chi; hence #Fix(T_g) = sum_{F Euler, gF = F} eps_F(g)
with eps_F(g) = chi_F(c_g). Burnside: #orbits = (1/n!) sum_F sum_{g in Aut F} eps_F(g); eps_F restricted to
Aut(F) is a homomorphism (cocycle identity), so the inner sum is |Aut F| or 0, giving the count of orbits of
Euler graphs with eps_F = 1 on Aut F. Independence of O / T_0: changing T_0 by y changes c_g by g y + y,
and chi_F(g y + y) = 1 for g in Aut F. Cycle form: on an edge orbit of length k the product of orientation
signs is the sign by which g^k acts on one of its edges, -1 exactly for diameter edges of even cycles. QED
For graphs (no cocycle) this is Mallows-Sloane's two-graphs = Euler graphs; the tournament torsor twists it.
*Verified:* untwisted Euler graphs counted directly (geng + C backtracking over all automorphisms):
1, 1, 1, 2, 2, 6, 12, 79, 792, 19576 for n = 1..10 = A049313 (Burnside closed form of THM-479 recomputed:
..., 19576, 886288, 75369960 for n = 10..12), out of A002854 = 1, 1, 2, 3, 7, 16, 54, 243, 2038, 33120 Euler
graphs; the per-permutation identity #Fix = sum eps_F(g) checked for every permutation n <= 6 and every
cycle type n = 7, 8; the cycle form agrees on every automorphism enumerated (n <= 10).
Examples: n = 3: only the empty graph (the triangle's transposition reverses one edge); n = 4: empty, C4
(not the triangle); n = 5: empty, C4 (C5, K5, bowtie, K3, K5 - K3 all twisted).
**Relation to repo work.** Refines THM-479's level law (g fixes a switching class iff eps_F(g) = 1 for all
g-invariant Euler graphs F) and is a second incarnation valid for all n, complementing THM-1470 (even
tournaments, n = 1 mod 4 only). Literature: the repo's OPEN-Q-060 record reports no known second incarnation;
no further literature search was done here, so priority is not claimed.

Second lever (new, not yet a named problem): the HP-arc parity results give two clean conjectures for the
backlog: C2 (Paley minus a vertex is Redei-rigid) and D1 (robust Redei: beta = Hall bound).

## 6. OPEN

* C2: prove Paley-minus-vertex all-odd for all q = 3 mod 4 (non-antipodal orbits open).
* What is the minimum number of odd arcs for even N (3, 5, 7, 5 at N = 4, 6, 8, 10)? Is it bounded?
* Conjecture C5: for N >= 3, an all-odd tournament exists iff N = 2 (mod 4). Evidence: none at N = 5, 9
  (exhaustive; N = 0, 3 mod 4 excluded by parity), present at N = 6, 10 (exhaustive/explicit) and
  18, 22, 26 (Paley minus a vertex). First open cases: N = 13 (non-existence), N = 14 (existence; no Paley
  tournament of order 15; none of the 128 circulants on Z_15 minus a vertex is all-odd, max 73 of 91 odd
  arcs; the doubled DRT(15) minus a vertex neither).
* D1: beta(T) = Hall bound (FINITE-EXACT N <= 9); in particular regular tournaments survive any N - 2 arc deletions.
* Redei shaving: is every non-all-odd tournament shavable to a single HP through odd H (FINITE-EXACT N <= 7)?
* C3: does some H-maximizer always have uniform arc parity when N is not 0 mod 4?
* E1 bijective form: an explicit bijection between switching classes of tournaments and
  orientation-parity-preserving Euler graphs (n = 1 mod 4 could go through THM-1470's even tournaments).

## 7. Process notes

* Runner: 133060 checks, ALL CHECKS PASSED, 224 s, peak RSS 322 MB (Python); the mod-2 DP child at q = 27
  uses 256 MB. Exhaustive N = 10 census: 59 min single process, RSS < 2 MB; independent mod-2 N = 10 pass
  about 10 min.
* Incident: one exploratory scratch test (per-permutation check at n = 8 enumerating all 2^28 edge sets)
  reached about 1.7 GB RSS for about 100 s before it was killed; it is not part of any deliverable. The
  runner uses the cycle-space basis instead (peak about 220 MB for that step).
* No git state was changed; files created: the deliverables listed in the header plus scratch under
  `scratch/procgen_selfie/`.
