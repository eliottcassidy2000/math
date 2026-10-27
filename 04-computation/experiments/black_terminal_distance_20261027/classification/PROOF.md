# The KEY LEMMA is false; the corrected classification and the exact black distance

Status:
- **REFUTED** (as stated): the KEY LEMMA fails first at n = 10 (exactly one class), again at n = 11 (exactly one
  class), and for an explicit infinite family E at every order n >= 10.
- **PROVED (all orders)**, modulo the cited repo results listed in Section 2:
  - **Theorem A (classification).** A strong tournament T on n >= 4 vertices has no black descent candidate
    if and only if T is the regular 5-tournament R5 or T belongs to the family E of Section 3.
  - **Theorem B (black distance).** For strong T on n >= 4 vertices, d_black(T) = lambda(T) unless T is R5
    (no black move at all) or T is in E, in which case d_black(T) = 2 = lambda(T) + 1.
- **FINITE-EXACT** controls: complete gentourng censuses n <= 11 (the descent test runs on every strong class for n <= 9 and
  on the nonregular self-converse classes for n = 10, 11, the others being covered by Part I), exhaustive structured
  searches n <= 16
  (Section 13). One sub-case of one lemma (the last-block lemma, case (m, p) = (6, 3)) is closed by a finite
  check rather than by an argument; everything else is a proof.

Conventions are those of the task: Q = (q_1..q_n) is a directed Hamiltonian path; blue iff r_Q
(q_i -> q_(n+1-i)) is an anti-automorphism; f(X) = e(X -> V\X); lambda = min f; a *descent candidate* is a Q with
q_1 in X and q_n not in X for some minimum out-cut X ("black descent" = black descent candidate). s(v) is the
outdegree, d^-(v) = n-1-s(v) the indegree, and

    delta(T) := min( min_v s(v), min_v d^-(v) ).

Singletons and co-singletons give lambda <= delta always. For self-converse T both minima coincide.
Throughout, a "move" is K_Q(T) for a black Q.

---------------------------------------------------------------------------------------------------------------

## 1. The order-10 counterexample, explicitly

Code (gentourng, upper triangle row by row): `111001111101111111111111111111111111000111101`.
Out-neighbourhoods:

    0 -> 1 2 3 6 7 8 9      (s=7)        5 -> 0 6          (s=2)
    1 -> 2 4 5 6 7 8 9      (s=7)        6 -> 7 8 9        (s=3)
    2 -> 3 4 5 6 7 8 9      (s=7)        7 -> 5 8          (s=2)
    3 -> 1 4 5 6 7 8 9      (s=7)        8 -> 5 9          (s=2)
    4 -> 0 5 6 7 8 9        (s=6)        9 -> 5 7          (s=2)

Structure. X = {5,6,7,8,9} carries `5 -> 6 => {7,8,9} => 5` with 7->8->9->7; Y = {0,...,4} carries
`0 => {1,2,3} => 4 -> 0` with 1->2->3->1; every X-Y pair points Y -> X except the single arc 5 -> 0.

- lambda = 1, delta = 2, and X is the **unique** minimum out-cut (its only leaving arc is 5 -> 0).
- Every descent candidate starts in X and ends in Y, so it crosses the cut exactly once, through 5 -> 0:
  it is (6, rotation of 7 8 9, 5, 0, rotation of 1 2 3, 4). These 9 paths are all closed (4 -> 6) and **all blue**.
- Aut(T) = Z3 x Z3 (independent rotations of {7,8,9} and {1,2,3}); there are 9 anti-automorphisms, 3 of them
  involutions, e.g. sigma = (0 5)(4 6)(1 7)(2 9)(3 8). Each candidate's reversal is one of them.
- T is not pure blue: 558 of its 657 Hamiltonian paths are black.
- Labelled BFS with black moves only (`black_distance.cpp`): level 1 has 558 black moves (108 targets with
  lambda = 1, 450 with lambda = 2, none non-strong); level 2 reaches 690 non-strong tournaments.
  **d_black = 2 = lambda + 1.** The "consequence" (d_black = lambda) therefore also fails.

Independent pure-Python verification (`verify_counterexample.py`, shares no code with the C++ checkers):
`results/verify_counterexample_n10.out` (1 minimum cut, 657 Ham paths, 9 candidates, 0 black).
The n = 11 census exception is checked the same way (`results/verify_counterexample_n11.out`: 2 minimum cuts,
1215 Ham paths, 9 candidates, 0 black).

---------------------------------------------------------------------------------------------------------------

## 2. Tools (proved earlier in the repo or below)

- **F1** (given in the task). f_{K_Q T}(V\X) = f_T(X) - [q_1 in X] + [q_n in X]; so
  lambda(K_Q T) >= lambda(T) - 1 with equality iff Q is a descent candidate. Score formulas for K = K_Q T
  (direct count): s_K(q_1) = n - s(q_1), s_K(q_n) = n - 2 - s(q_n), s_K(v) = n - 1 - s(v) otherwise.
- **F2.** A blue Q satisfies s(q_i) + s(q_(n+1-i)) = n - 1 for all i (an anti-automorphism maps out-neighbourhoods
  to in-neighbourhoods). A non-self-converse T has no blue word at all.
- **F3.** If tau is an anti-automorphism and X is a minimum out-cut, then V \ tau(X) is a minimum out-cut:
  f(X) = #{(u,v): u in X, v notin X, tau(v) -> tau(u)} = f(V \ tau(X)).
- **F4** (Camion / Redei). A strong tournament has a Hamiltonian cycle; a vertex v starts a Hamiltonian path
  iff v reaches every vertex (it then lies in the initial strong component; concatenate).
- **F5** (LEM-003). Aut(D) acts freely on Hamiltonian paths of any digraph.
- **F6** (THM-2477 as superseded in `merged_metagraph_finite_pure_atoms_20261020`, Section 3, via
  Benhocine-Wojda's Hamiltonian bypass theorem). The only nontrivial strong tournament with H = |Aut| is C3.
  With F5: if a subgroup G of Aut(R), R strong, acts transitively on Ham(R), then R is K1 or C3 and G = Aut(R).
- **F7** (`blue_quotient_20261020b/structure.md`, Sections 2-3). For a Hamiltonian cycle C with rotation rho, the
  blue cuts form a coset B_C = i + H_C of H_C = {h : rho^(2h) in Aut(T)} (when nonempty); for even n,
  n/2 is in H_C, so |B_C| is even. If all cuts are blue and n is even, there is a Hamiltonian cycle containing an
  arc u1 -> u2 inside one rho^2-orbit U; the orbit scores are alpha (on U) and n-1-alpha. For odd n, all cuts blue
  forces rho in Aut(T), hence T regular. Regular T other than C3, R5 has a black closed Hamiltonian path (for odd
  circulants of order >= 9 this rests on the audited uniform prefix-rerouting theorem,
  `pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`).
- **F8** (`black_descent_20261020c/cycles.md`, Section 3, prime-orbit block normal form). If C is a Hamiltonian
  cycle, n = p*l with p an odd prime and l >= 3, and the rotation by l is an automorphism, then T has a
  Hamiltonian cycle in which each of its l orbits (size p) is one consecutive block.
- **F9** (cycles.md, Sections 4-5). If 1 < |B_C| < n for odd n, or 2 < |B_C| < n for even n, then for an odd prime
  p dividing |B_C| (resp. the odd number |B_C|/2) the rotation by l = n/p is an automorphism, with l >= 3 odd
  (resp. l >= 4 even).
- **F10** (cycles.md, Section 7). If T is strong and self-converse, n >= 4, v has minimum score and T - v is not strong, every
  Hamiltonian-cycle cut starting at v ends at a vertex of score < Delta; it is black by F2.
- **F11.** A self-converse tournament has an **involutive** anti-automorphism sigma (take tau^m, m the odd part of
  ord(tau): its square is an automorphism of 2-power order, and tournament automorphisms have odd order), with at
  most one fixed point (two fixed points u, v would give u->v iff v->u).

---------------------------------------------------------------------------------------------------------------

## 3. The exceptional family E

**Bottom layers.** A *bottom layer* is a pair (S, x0): S a tournament, x0 in S, such that S - x0 is a stack
R_1 => R_2 => ... => R_q (q >= 2; every earlier block dominates every later one) of blocks R_i in {K1, C3},
with R_q = C3, and x0 uniform on every block, x0 => R_1, R_q => x0; the relation of x0 to R_2..R_(q-1) is free.
The number b(a) of bottom layers with |S| = a is

    sum_a b(a) z^a = z^5 (1 + z^2) / (1 - 2z - 2z^3):   b(5..11) = 1, 2, 5, 12, 28, 66, 156.

(Compositions of |S| - 4 into parts 1, 3, times 2^(q-2) free bits.) Census check: `bottom_layers.cpp`
classifies, over all strong S with a <= 9 vertices and all x0, the pairs satisfying Condition A and transitivity
(Section 9) and finds exactly 1, 2, 5, 12, 28 pairs, all of this form.

**The tournaments E_t(S, x0)** (t >= 1). Vertex set L_0 = S, middle vertices c_1..c_(t-1), and a top layer
L_t = sigma(S) carrying the converse of S (sigma(u) -> sigma(v) iff v -> u). The only *upward* arcs are the jumps

    x0 -> c_1 -> c_2 -> ... -> c_(t-1) -> sigma(x0)          (t = 1: x0 -> sigma(x0));

every other pair between different layers points down (L_t => middles => S, and c_j -> c_i for j >= i+2).
Order n = 2|S| + t - 1. sigma extended by c_i <-> c_(t-i) is an involutive anti-automorphism.

The order-10 example is E_1(S5, x0) with S5 = (x0 -> K1 => C3 => x0); the order-11 census exception is
E_2(S5, x0). The number of members at order n is e(n) = sum_{a=5}^{floor(n/2)} b(a):

    n    : 10 11 12 13 14 15 16 17 18 19  20  21
    e(n) :  1  1  3  3  8  8 20 20 48 48 114 114

Relation to the coordinator's family E(M; Z): it is the sub-family t = 1 with x0 beating an initial segment of
the stack. It misses every t >= 2 member (all odd orders, e.g. the n = 11 census class, and E_3(S5) at n = 12)
and the patterns that are not initial segments (first at n = 14: S - x0 = (K1, K1, K1, C3) with
x0 -> R_1, R_2 -> x0, x0 -> R_3, R_4 => x0; S code `101000111111111111101`, x0 = 0).

---------------------------------------------------------------------------------------------------------------

## 4. Part I: non-self-converse and regular tournaments

**I.1.** If T is strong and not self-converse, every Hamiltonian path is black (F2). Take a minimum out-cut X
and a Hamiltonian cycle (F4); the cycle enters X along some arc u -> w, and its cut starting at w is a closed
black descent candidate. If moreover lambda = delta, take X = {v} with s(v) = lambda (or X = V \ {u} with
d^-(u) = lambda): any Hamiltonian path from v (into u) is a *singleton-type* black descent, used in Part VIII.

**I.2.** If T is regular, every singleton is a minimum out-cut (f(X) = |X|(n-|X|)/2 >= (n-1)/2), so every black
word is a descent candidate; by F7 one exists unless T is C3 or R5. R5 (and C3) has no black word at all.

Hence: **every exception (strong, n >= 4, no black descent, not R5) is self-converse and nonregular.**

---------------------------------------------------------------------------------------------------------------

## 5. Part II: lambda = delta implies a singleton-type black descent

**Theorem C.** Let T be strong, n >= 4, T != R5, with lambda(T) = delta(T). Then T has a black descent candidate Q
of *singleton type*: s(q_1) = lambda (so {q_1} is a minimum cut) or d^-(q_n) = lambda (so V \ {q_n} is one).

*Proof.* Non-self-converse and regular T: Part I. So let T be self-converse and nonregular. Let m (resp. M) be
the set of vertices of minimum score delta (resp. maximum score Delta = n-1-delta); m and M are disjoint,
sigma(m) = M, |m| = |M| =: mu. Because lambda = delta, every {v} with v in m and every V \ {u} with u in M is a
minimum out-cut; so every Hamiltonian path starting in m and every one ending in M is a descent candidate.

(c1) If some v in m has T - v non-strong, F10 gives a closed black candidate starting at v. Done.

(c2) Otherwise T - v is strong for every v in m. Suppose, for a contradiction,

    (H)  every Hamiltonian path starting in m, and every Hamiltonian path ending in M, is blue.

By F2 a blue path starting at v in m ends in M.

*Step 1 (cycle cosets give mu <= 2, and mu = 1 for odd n).* Let C = (c_0..c_(n-1)) be any Hamiltonian cycle and
J_m = {j : c_j in m}. The cut starting at c_j (j in J_m) is blue by (H), so J_m is contained in B_C, which is then
a coset i + H_C (F7). We show |B_C| = 1 (n odd), |B_C| = 2 (n even) (B_C contains the nonempty J_m):
- If B_C = Z_n: n odd makes T regular (F7), excluded. For n even, F7 gives a Hamiltonian cycle C' with an arc
  u1 -> u2 inside one rho^2-orbit of common score alpha; the scores of T are alpha and n-1-alpha, so
  alpha is delta or Delta. The cut of C' starting at u2 ends at u1 with score sum 2*alpha != n-1, hence is black;
  it starts in m (alpha = delta) or ends in M (alpha = Delta), contradicting (H).
- If 1 < |B_C| < n (n odd) or 2 < |B_C| < n (n even): F9 gives an odd prime p and l = n/p >= 3 with
  rotation by l in Aut(T); F8 gives a Hamiltonian cycle C' on which every orbit (size p >= 3) is consecutive.
  m is Aut-invariant, hence a union of orbits; take an orbit O in m and consecutive o1 -> o2 in O on C'. The cut
  of C' starting at o2 ends at o1, score sum 2*delta != n-1: black, starts in m. Contradiction.
Hence mu = |J_m| <= |B_C| <= 2, and mu = 1 when n is odd.

*Step 2 (mu = 1 is impossible for n >= 5).* m = {v}, M = {u}, sigma(v) = u. By (H) every Hamiltonian cycle of T
has u as predecessor of v. Let D be a Hamiltonian cycle of T - v (strong by (c2)); A = N+(v), B = N-(v).
For every transition d_(i-1) in B, d_i in A on D, the cycle (v, d_i, ..., d_(i-1)) shows d_(i-1) = u; so A is a
single run (a_1..a_delta) on D preceded by u. If delta >= 2, the path (v, a_2, ..., a_1) along D is a
Hamiltonian path from v ending at a_1 in A, but by (H) it ends in M = {u}, and u is in B: contradiction. So
delta = 1, N+(v) = {a}. As Delta = n-2, u has a unique in-neighbour a*, and sigma(a) = a* (apply sigma to v -> a);
a != u and a* != v. Put T'' = T - {u, v}; sigma maps T'' onto itself as an anti-automorphism.
*Claim:* every Hamiltonian path R of T'' starting at a ends at a*, and such R exist. (A Hamiltonian path of T from
v is (v, a, ..., a*, u); if some R = (a, ..., a*, w, ...) had a* inside, then (v, a, ..., a*, u, w, ...) would be a
Hamiltonian path from v not ending in M, against (H). In particular a != a* as n >= 4.)
Fix such an R = (a = p_1, ..., p_N = a*), N = n - 2, and put W = {p_2, ..., p_(N-1)}.
(i) If p_(j-1) -> a* -> p_j for some 2 <= j <= N-1, then (p_1..p_(j-1), a*, p_j..p_(N-1)) is a Hamiltonian path of
T'' from a ending at p_(N-1) != a*: excluded. So along p_1..p_(N-1) the out-neighbours of a* come first.
(ii) If a* -> p_2 and a -> p_l with 3 <= l <= N, then (a, p_l, ..., p_N, p_2, ..., p_(l-1)) is a Hamiltonian path of
T'' from a ending at p_(l-1) != a*: excluded.
Now a is not in m = {v}, so s(a) >= 2, and all out-neighbours of a lie in T'' (v -> a, and u -> a since a != a*);
so a has an out-neighbour p_l with l >= 3, and (ii) gives p_2 -> a*. Then (i) gives W => a*, and applying sigma
(sigma(a*) = a, sigma(W) = W) gives a => W. If n >= 6, take a Hamiltonian path (w_1..w_k) of T[W], k >= 2: the word
(v, a, w_1, a*, u, w_2, ..., w_k) is a Hamiltonian path of T (a -> w_1, w_1 -> a*, a* -> u, and u beats every
vertex but a*) from v ending at w_k, outside M. If n = 5, W = {w} and s(a) >= 2 forces a -> a*; the word
(v, a, a*, u, w) is a Hamiltonian path ending at w, outside M. Both contradict (H).

*Step 3 (mu = 2 is impossible for n >= 6, n even).* m = {v, v'}; both deletions are strong; name them so that
v -> v'. If delta = 1 then N+(v) = {v'} and v has outdegree 0 in T - v', against (c2); so delta >= 2. Let D be a
Hamiltonian cycle of T - v. For a in A = N+(v), the path (v, a, ..., pred_D(a)) along D is blue by (H), so
pred_D(a) is in M; distinct a give distinct predecessors, so delta = |A| <= |M| = 2: A = {v', a}. The blue path
(v, v', ..., pred_D(v')) has s(pred_D^2(v')) = n-1-s(v') = Delta (positions n-1 and 2, F2); so
M = {pred_D(v'), pred_D^2(v')}, which forces pred_D(a) = pred_D^2(v'), i.e. a = pred_D(v'). The blue path
(v, a, v', ..., pred_D(a)) then has s(pred_D^3(v')) = n-1-s(a) = delta, so pred_D^3(v') is in m \ {v} = {v'}:
D has length 3 and n = 4. Contradiction.

n = 4 falls under (c1) (its only strong class has a min-score vertex with non-strong deletion). So (H) fails:
some Hamiltonian path starting in m or ending in M is black; it is a singleton-type descent candidate. QED

Remark (what (c1)-(c2) prove without lambda = delta). The case analysis never uses lambda = delta (only the last
sentence does), so it proves: **every nonregular strong self-converse T (n >= 4) has a black Hamiltonian path starting at a
minimum-score vertex** (a black path ending in M becomes one starting in m by sigma and reversal). This is weaker than
the "extremal predecessor" statement left OPEN in cycles.md Section 7, both in its closed form and in its open-path form
(a Hamiltonian path from a minimum-score vertex to a vertex of less than maximum score): a black path from m may end in M,
and Steps 1 and 3 use the blueness in (H) (the coset lemma, F2 at inner positions), not only "every path from m ends in
M". The argument settles the open-path form only in case (c1) and when mu = 1 (Step 2 uses (H) only through "every
Hamiltonian path from v ends in M"), and the closed form only in case (c1) and when mu = 1 and delta = 1 (then every
Hamiltonian path from v is closed). Both forms remain OPEN; an independent check finds the closed form true for every
nonregular strong self-converse class of order <= 10.

---------------------------------------------------------------------------------------------------------------

## 6. Part III: structure when lambda < delta

**III.1 (size).** If lambda < delta, every minimum out-cut X has |X| >= 2*delta + 1 and |V\X| >= 2*delta + 1;
hence n >= 4*delta + 2 >= 10. *Proof.* f(X) >= g(|X|) := |X|*delta - C(|X|,2); g is concave with
g(1) = g(2*delta) = delta, so f(X) >= delta > lambda whenever |X| <= 2*delta; similarly with indegrees for V\X.

**III.2 (minimum cuts never cross; chain).** In any digraph
f(X) + f(Z) = f(X cap Z) + f(X cup Z) + e(X\Z, Z\X) + e(Z\X, X\Z), and in a tournament the last two terms add to
|X\Z| |Z\X|. If X, Z are minimum cuts with X cap Z, X cup Z proper and nonempty and X\Z, Z\X nonempty, then
f(X cap Z) + f(X cup Z) <= 2*lambda - 1: impossible. So minimum cuts of a tournament never cross. If they are
disjoint (union proper), f(X cup Z) = 2*lambda - |X||Z| forces |X||Z| <= lambda; if co-disjoint,
f(X cap Z) = 2*lambda - |V\X||V\Z| forces |V\X||V\Z| <= lambda; if Z = V \ X, then
f(X) + f(V\X) = |X||V\X| = 2*lambda. Under lambda < delta, III.1 excludes all three (each product is at least
(2*delta + 1)^2 > 2*lambda).
**So when lambda < delta the minimum cuts form a chain X_1 < X_2 < ... < X_t** (strict inclusions).
Layers: L_0 = X_1, L_k = X_(k+1) \ X_k, L_t = V \ X_t; M_mid := L_1 cup ... cup L_(t-1).
T[L_0] and T[L_t] are strong (L_0 is an inclusion-minimal minimum cut: a non-strong T[L_0] would have a last
component Z with f(Z) <= f(L_0); dually for L_t).

**III.3 (anti-automorphisms reverse the chain).** By F3, X -> V \ tau(X) is an inclusion-reversing permutation
of the chain, hence X_i -> X_(t+1-i): tau(X_i) = V \ X_(t+1-i) and tau(L_k) = L_(t-k).

**III.4 (layer reversal).** If Q is blue, then layer(q_p) + layer(q_(n+1-p)) = t for every p (apply III.3 to
tau = r_Q). In particular a descent candidate whose layer word is not "anti-palindromic" is black, and a
candidate starting in L_0 must end in L_t to be blue.

**III.5 (forward arcs).** Call an arc from L_a to L_b with a < b *forward*; it leaves X_(a+1), ..., X_b and no
other X_i. Since f(X_i) = lambda for all i, every position i in [1, t] is crossed by exactly lambda forward arcs;
all other pairs between different layers point down. Every vertex of U_1 := V \ L_0 is reachable inside T[U_1]
from a head of a forward arc leaving L_0 (follow any path from L_0 and look after its last L_0-vertex).
Consequently (F4) **the initial strong component of T[U_1] contains a head y of a forward arc x -> y (x in L_0),
and T[U_1] has a Hamiltonian path starting at y.**

For an exception, T is self-converse (Part I) and sigma (F11) satisfies sigma(L_k) = L_(t-k) (III.3).

---------------------------------------------------------------------------------------------------------------

## 7. Part IV: 2 <= lambda < delta implies a black descent

**Theorem D.** If T is strong with 2 <= lambda(T) < delta(T), T has a black descent candidate.

*Proof.* Suppose not; T is self-converse. Take x -> y from III.5, any P_0 in Ham(T[L_0] -> x) (T[L_0] strong) and
any P' in Ham(T[U_1] from y). Q = (P_0, P') is a candidate for X_1, hence blue; its first |L_0| positions are L_0,
so by III.4 its last |L_t| = |L_0| positions are exactly L_t, i.e. P' = (M, P_t) with M a Hamiltonian path of
T[M_mid] from y (if t >= 2 and y were in L_t this is already impossible) and P_t a Hamiltonian path of T[L_t]
starting at some y'' with x'' -> y'' a forward arc (x'' = last vertex of M, or x'' = x when t = 1).
Now fix M and consider the whole family Q' = (P_0', M, P_t') with P_0' in Ham(T[L_0] -> x),
P_t' in Ham(T[L_t] from y''). All are candidates for X_1, hence blue. For blue Q', psi := sigma r_Q' is an
automorphism preserving L_0 with psi(x) = sigma(y'') =: z, and P_t' = rev(sigma psi(P_0')). For fixed P_0', distinct
blue P_t' come from distinct restrictions psi|L_0 (every automorphism of T preserves L_0, the smallest minimum
cut); these restrictions lie in one coset g * Stab(x) of the stabiliser of x in Aut(T[L_0]), with g(x) = z. So
x and z lie in one Aut(T[L_0])-orbit, their stabilisers are conjugate, and
h(z) := |Ham(T[L_0] -> z)| = |Ham(T[L_t] from y'')| <= |Stab_{Aut T[L_0]}(z)| (sigma reverses Hamiltonian paths
between L_0 and L_t). By F5 the reverse inequality holds, so **Stab_{Aut T[L_0]}(z) acts simply transitively on
Ham(T[L_0] -> z)**. Note that z is the
tail of the forward arc sigma(y'') -> sigma(x'') leaving L_0, so lambda(z) >= 1, where lambda(w) := number of
forward arcs leaving w; sum over L_0 of lambda(w) = lambda.

Apply Lemma BL (Section 9) to (T[L_0], z): T[L_0] - z ends with a block R_q = C3 with R_q => z. Each r in R_q has
s(r) = 2 + lambda(r) >= delta >= lambda + 1, so lambda(R_q) >= 3(lambda - 1); but lambda(R_q) <= lambda - lambda(z)
<= lambda - 1. Hence lambda <= 1, a contradiction. QED

(Hypotheses of Lemma BL used here: every nonempty W strictly inside L_0 has f(W) >= lambda + 1, since all minimum
cuts contain L_0; and lambda(R_q) <= lambda - 1. They give (H2) of the last-block lemma.)

---------------------------------------------------------------------------------------------------------------

## 8. Part V: lambda = 1 < delta, the exceptions are exactly E

Let T be an exception with lambda = 1 < delta. Notation of Part III; each X_i has a unique leaving arc, the jumps
are e^(k): L_(a_(k-1)) -> L_(a_k) with 0 = a_0 < a_1 < ... < a_s = t (III.5 with lambda = 1). A candidate for X_i
crosses X_i exactly once, so it is (Hamiltonian path of T[X_i] ending at the tail) + (Hamiltonian path of
T[V\X_i] starting at the head). By III.5 (only one head) T[U_1] has a Hamiltonian path from y^(1).

**V.1 (no skipped layer).** Suppose some layer L_c (1 <= c <= t) is not the head layer of a jump, say
a_k < c < a_(k+1). Let P' be a Hamiltonian path of T[U_1] from y^(1) and P_0 one of T[L_0] ending at x^(1).
No forward arc enters L_c, so P' enters L_c from above; the only forward arc crossing positions c and c+1 is
e^(k+1). If k = 0 it is the arc joining P_0 to P'; if k >= 1, P' starts below L_c and must use e^(k+1) to get
above L_c before entering it. Either way, once P' has visited L_c it can never rise above L_c again, so P' ends in
a layer <= c < a_(k+1) <= t, and (P_0, P') is a candidate for X_1 that is black by III.4.
Hence a_k = k: the jumps are e_k = x_k -> y_k with x_k in L_(k-1), y_k in L_k.

**V.2 (middle layers are singletons).** Let k be minimal in [1, t-1] with |L_k| >= 2 (so L_1..L_(k-1) are
singletons c_i = y_i = x_(i+1)). Every simple (y_k, x_(k+1))-path stays inside L_k (leaving L_k downwards needs
e_k to come back, i.e. a second visit of y_k; leaving upwards needs e_(k+1), i.e. x_(k+1)); one exists as T is
strong. Every vertex of U_(k+1) is reachable from y_(k+1) inside T[U_(k+1)] (e_(k+1) is the only arc entering
U_(k+1)), so T[U_(k+1)] has a Hamiltonian path from y_(k+1) (F4). If some (y_k, x_(k+1))-path R misses a
vertex of L_k (e.g. y_k = x_(k+1)), the word
(P_0 ending at x_1, c_1..c_(k-1), R, Hamiltonian path of T[U_(k+1)] from y_(k+1), Hamiltonian path of T[L_k \ R])
is a Hamiltonian path (the step from U_(k+1) down to L_k \ R is a downward arc: the only forward arc between L_k
and U_(k+1) starts at x_(k+1), which lies on R), a candidate for X_1 ending in
L_k, black by III.4. Otherwise every (y_k, x_(k+1))-path is Hamiltonian in L_k; then the BFS levels from y_k in
T[L_k] are singletons, so y_k has exactly one out-neighbour in L_k and f(X_k cup {y_k}) = 1: a minimum cut
strictly between X_k and X_(k+1), contradiction. So all middle layers are singletons {c_k}, and
x_1 -> c_1 -> ... -> c_(t-1) -> y_t are the only upward arcs.

**V.3 (candidates).** c_j's only in-neighbour inside X_(j+1) is its chain predecessor and its only out-neighbour
outside X_j is its successor, so every candidate, for every X_i, is Q = (P_0, c_1, ..., c_(t-1), P_t) with
P_0 in Ham(S -> x_1), S := T[L_0], and P_t in Ham(T[L_t] from y_t). With sigma (sigma(L_0) = L_t,
sigma(c_i) = c_(t-i), sigma(x_1) = y_t), T is exactly E_t(S, x_1) as a tournament (all non-jump inter-layer pairs
point down). As in Part IV, all h^2 candidates are blue iff Stab_{Aut S}(x_1) acts transitively on
Ham(S -> x_1) (the number of blue ones is h * |Stab|; every element of the stabiliser extends to an automorphism
of T: phi on L_0, identity on the c_i, sigma phi^(-1) sigma on L_t). The absence of minimum cuts inside L_0 is
**Condition A**: every nonempty W subset of S - x_1 has e_S(W -> S\W) >= 2.

**V.4.** By Lemma BL (with lambda = 1, z = x_1), (S, x_1) is a bottom layer (Section 3). So T is in E.

---------------------------------------------------------------------------------------------------------------

## 9. Part VI: Lemma BL and the last-block lemma

**Lemma BL.** Let S be strong, z in S, and suppose Stab_{Aut S}(z) acts transitively on Ham(S -> z). Write
S - z = R_1 => ... => R_q (strong components), A_i = N+(z) cap R_i, B_i = R_i \ A_i (B_q is nonempty because
S is strong and R_q has no arc to the other blocks). Assume
(H2) every nonempty W subset of R_q has e_{R_q}(W -> R_q \ W) + |W cap B_q| >= 2, and |R_q| >= 2.
Then every R_i (i < q) is K1 or C3 with z uniform on it, z => R_1, R_q = C3 and R_q => z.

*Proof.* Automorphisms fixing z preserve the condensation, so Stab(z) = prod_i Aut(R_i; A_i), and
Ham(S -> z) = prod_{i<q} Ham(R_i) x Ham(R_q -> B_q); a product group is transitive on a product set iff each
factor is. For i < q, F6 gives R_i in {K1, C3} and Aut(R_i; A_i) = Aut(R_i), so A_i is empty or R_i. For R_q
apply the last-block lemma. z => R_1 because S is strong; q >= 2 because z has an out-neighbour. QED

In Part V, (H2) is Condition A restricted to R_q, and |R_q| >= 2 because a singleton R_q = {r} has outdegree
[r -> x_1] <= 1. In Part IV, (H2) follows from f(W) >= lambda + 1 and lambda(W) <= lambda - 1 (Section 7), and a
singleton R_q = {r} would have s(r) <= 1 + (lambda - 1) < delta.

**Last-block lemma.** Let R be strong on m >= 3 vertices, V(R) = A disjoint-union B, B nonempty,
G = Aut(R; A), and assume (H1) G acts transitively on H' = Ham(R -> B) and (H2). Then A is empty and R = C3.

*Proof.* If A is empty, H' = Ham(R) and F6 gives R = C3. Let A be nonempty.
(a) For a Hamiltonian cycle C = (c_0..c_(m-1)) and b = c_j, b' = c_j' in B, the cuts of C ending at b, b' lie in
H'; the element of G carrying one to the other (positionwise) is the rotation by j' - j. So all rotations by
differences of B-positions lie in G, B is a single coset of the rotation group D_C they generate, and since A is
nonempty, consecutive B-vertices on C are at distance l_C = m/|B| >= 2.
(b) (H2) with W = V(R) gives |B| >= 2; D_C has odd order (tournament automorphisms), so |B| = |D_C| >= 3 is odd.
(c) Let p be a prime dividing |B| and l' = m/p = l_C * |B|/p >= 2; the rotation by l' (a power of the rotation by
l_C, which generates D_C) is in D_C and has order p.
 - l' >= 3: F8 yields a Hamiltonian cycle C' on which each rotation orbit (size p >= 3) is consecutive; B is a
   union of these orbits, so C' has consecutive B-vertices, contradicting (a) for C'.
 - l' = 2: then m = 2p, B and A are the two orbits of rho^2, each a circulant on Z_p, u_i -> w_(i+d) iff d in S_X
   (S_X proper and nonempty: 0 in S_X, -1 not in S_X from C). Choose eps in conn(B), eps' in conn(A) with
   eps + eps' != 0 (possible unless p = 3), then t with t + eps in S_X and t - eps' not in S_X (translation by
   eps + eps' moves the proper subset S_X); the word (u_0, u_eps, ..., u_(-eps), w_t, ..., w_(t-eps')) is a
   Hamiltonian cycle with B consecutive, contradicting (a).
 - (m, p) = (6, 3): **finite check** `lastblock_check.cpp`: over all strong R on m <= 7 vertices and all A, the
   only pair satisfying (H1)+(H2) is C3 with A empty (`results/lastblock_m3-7.out`). QED

---------------------------------------------------------------------------------------------------------------

## 10. Part VII: every member of E is an exception (converse of Part V)

Let (S, x0) be a bottom layer and T = E_t(S, x0).
1. *Condition A holds for S:* for nonempty W in S - x0: if |W cap R_q| >= 2 there are two arcs to x0; if
   W cap R_q = {r}, r beats x0 and its C3-successor; if W misses R_q, W beats all of R_q.
2. T is strong (every vertex of S reaches x0 inside S, x0 climbs the jumps to sigma(x0), T[L_t] is strong, and
   every vertex of L_t beats every vertex below except along the jumps, so L_t reaches all of S and all c_i);
   lambda = 1 (f(X_1) = 1 with X_k = L_0 cup {c_1..c_(k-1)}); delta = 2 (outdegrees: >= 2 on
   S - x0 by Condition A, s(x0) >= 2, the c_i and L_t have outdegree >= |S| - 1; indegrees by sigma; the three
   vertices of R_q have outdegree exactly 2).
3. The minimum cuts are exactly X_1..X_t: they form a chain (III.2); a minimum cut inside L_0 contradicts
   Condition A (or strongness of S, if it contains x0); one above X_t would give, via sigma, one inside L_0.
4. The candidates are the h^2 words (P_0, c_1..c_(t-1), P_t) of V.3, h = 3^(number of C3 blocks).
5. Stab_{Aut S}(x0) contains the independent rotations of the C3 blocks (x0 is uniform on blocks), a group of
   order h acting freely (F5) on the h paths of Ham(S -> x0) = prod Ham(R_i): simply transitive. For any
   candidate, with phi the stabiliser element sending P_0 to rev(sigma(P_t)), psi = (phi on L_0, id on c_i,
   sigma phi^(-1) sigma on L_t) is an automorphism of T and r_Q = sigma psi. So every candidate is blue.

So the members of E are exactly the exceptions with lambda = 1 < delta (Parts V-VII), and with Parts I, II, IV:

**Theorem A.** For strong T on n >= 4 vertices: T has no black descent candidate iff T = R5 or T is in E.

---------------------------------------------------------------------------------------------------------------

## 11. Part VIII: the black distance

d_black(T) >= lambda(T) by F1.

1. **lambda = delta (all T with n <= 9 are of this kind, by III.1).** Theorem C gives a singleton-type black descent
   Q. By F1's score formulas, if s(q_1) = lambda then d^-_K(q_1) = lambda - 1; if d^-(q_n) = lambda then
   s_K(q_n) = lambda - 1; so K = K_Q T satisfies lambda(K) = lambda - 1 = delta(K) (lambda <= delta always). K is
   not in E (members of E have lambda < delta) and is not R5 (that would need lambda(T) = 3 at n = 5), so by
   induction d_black = lambda.
2. **lambda = 1 < delta, T not in E:** a black descent exists (Theorem A), d_black = 1.
3. **2 <= lambda < delta:** Theorem D gives a black descent; K has lambda - 1 and delta(K) >= delta - 1 (score
   formulas), so it stays in the regime lambda < delta. When lambda(T) = 2 (delta >= 3), every vertex of K except
   possibly q_n has outdegree >= 3 (s_K(v) = d^-_T(v) >= 3, s_K(q_1) = d^-_T(q_1) + 1), while every member of E has
   at least three vertices of outdegree 2 (the block R_q). So K is not in E and has a black descent:
   d_black = lambda. For lambda(T) >= 3 the tournament K has 2 <= lambda(K) < delta(K), is not R5 (n >= 18) and
   not in E (lambda(K) >= 2), so Theorem D applies again; induction gives d_black = lambda.
4. **T in E: d_black = 2.** Let r -> r' -> r'' be the block R_q and
   Q = (r, r', r'', x0, c_1, ..., c_(t-1), P_top, P_rest) with P_top a Hamiltonian path of L_t from sigma(x0)
   and P_rest one of R_1 cup ... cup R_(q-1). Q is a Hamiltonian path starting and ending in L_0, hence black
   (III.4, 0 + 0 != t) and not a candidate; so lambda(K_Q T) >= 1, while d^-_K(r) = s(r) - 1 = 1. Thus
   lambda(K_Q T) = 1 = delta(K_Q T), and K_Q T has a black descent by Theorem C: d_black <= 2. Conversely T has no
   black descent candidate (Part VII), so by F1 no single black move reaches lambda = 0: d_black >= 2.

**Theorem B.** d_black(T) = lambda(T) for every strong T on n >= 4 vertices other than R5 and the members of E;
d_black(T) = 2 for T in E; R5 (like C3) has no black move.

---------------------------------------------------------------------------------------------------------------

## 12. The stronger (closed, score-certified) form

It is false beyond n = 9 even outside E: among self-converse classes it fails for exactly 2 classes at n = 10
(E_1(S5), and `110011011110111111111111111111111111100110111`, which is F(R5, x0): black candidates exist but every
start has S-score 2) and exactly 2 at n = 11 (E_2(S5) and one more). The proofs above never need it.
Non-self-converse classes at n >= 10 were not tested for it (irrelevant to the KEY LEMMA).

---------------------------------------------------------------------------------------------------------------

## 13. Computations (all in this directory; `tour.hpp` is the shared C++ toolkit)

| program | what it checks | result |
|---|---|---|
| `keylemma_check.cpp` (gentourng n = 4..9) | C1 min cuts never cross / chain; C2 blue reversals permute min cuts; C3 black descent iff T != R5; C4 lambda = delta propagation (singleton-type black descent with lambda(K) = delta(K)); C5 F10 case; C6 case (c2) of Theorem C; C7 cut identity on every (Q, X) | 0 failures on 1, 6, 35, 353, 6008, 178133 strong classes (C7 at n <= 8); `results/keylemma_check_n*.out` |
| `census_sc.cpp` (gentourng n = 10, all 9,355,949 strong classes) | KEY on the 8,599 nonregular self-converse classes | exactly 1 failure = E_1(S5); strong form fails 2 |
| `census_sc.cpp` (gentourng n = 11, 8 slices, 884,464,590 strong classes) | KEY on 276,971 nonregular self-converse classes | exactly 1 failure = E_2(S5); strong form fails 2 |
| `bottom_layers.cpp` (a = 3..9) | pairs (S, x0) with Condition A and transitive stabiliser | 0,0,1,2,5,12,28 pairs, all bottom layers |
| `lastblock_check.cpp` (m = 3..7) | last-block lemma incl. the (6,3) case | only C3 with A empty |
| `family1.cpp` (all strong S, a = 5,6,7) | every F(S,x0) (all even-order exceptions have this shape) | KEY fails 1, 3, 8 at n = 10, 12, 14 = e(n) |
| `family_odd.cpp` (a = 4..7) | every odd shape with a sigma-fixed middle vertex | KEY fails 0, 1, 3, 8 at n = 9, 11, 13, 15 = e(n) |
| `family2.cpp` (n = 14) | all lambda = 2 < delta configurations (t = 1 is forced at n = 14, so exhaustive in that regime) | 430 configurations, 0 failures |
| `family3.cpp` (n = 15) | lambda = 2 < delta with a fixed middle vertex (t = 2); exhaustive in that regime, since L_0 and L_t have 7 vertices each and the single middle vertex is sigma-fixed | 1131 configurations, 0 failures |
| `e_check.cpp` | all 44 members of E with n <= 16: strong, lambda = 1, delta = 2, cut chain = t, no black descent, explicit two-move certificate of Part VIII.4, >= 3 vertices of outdegree 2 | 44/44 |
| `black_distance.cpp` | labelled BFS from the n = 10 exception; cut identity on all 671,454 (Q, X) pairs | distance 2, 0 identity failures |
| `verify_counterexample.py` | independent Python check of the n = 10 and n = 11 exceptions | confirmed |
| `audit/` (independent auditor, own code, no `tour.hpp`) | last-block lemma for all strong R with m <= 9 and all A; Lemma BL hypotheses for a <= 9; Theorem C on 78,387 random self-converse nonregular strong T with lambda = delta, n = 10..13; the 19 members of E with n <= 15 | 0 counterexamples; BL pairs 1, 2, 5, 12, 28; all 78,387 have a singleton-type black descent; 19/19 |
| second independent audit (2026-10-27, own C + nauty code) | Theorems A and C on all strong classes n <= 9 and all nonregular self-converse classes n = 10; E from its definition (112 members, n <= 18); last-block lemma m <= 9 with (H1) tested directly; Lemma BL a <= 9; exhaustive lambda = 1 < delta searches n <= 17 and lambda = 2 < delta searches n = 14..17; all F/G shapes with a <= 8; exact class-level black BFS n <= 9 | no counterexample; every exception found is R5 or a member of E; d_black = lambda for every strong class n <= 9 except R5 |

Commands: `g++ -O2 -std=c++17 -o X X.cpp`, then e.g. `gentourng -cq 8 | ./keylemma_check 8`,
`gentourng -cq 10 | ./census_sc 10`, `gentourng -cq 11 r/8 | ./census_sc 11` (r = 0..7),
`gentourng -cq 7 | ./bottom_layers 7`, `grep '^S ' bl_full_5.txt | ./e_check 5 16`.
Outputs are in `results/`. Exploration programs (`explore*.cpp`, `describe*.cpp`, `allblue_start.cpp`,
`family_E.cpp`, `tightinfo.cpp`, `antiaut_list.cpp`) are kept for reproducibility only.

---------------------------------------------------------------------------------------------------------------

## 14. Honest notes

- Cited, not re-proved here: THM-2477 (via Benhocine-Wojda 1983, a published theorem), LEM-003, the coset lemma
  and even all-blue construction (structure.md), the regular-case black path (structure.md Section 2, via the
  uniform prefix-rerouting theorem for odd circulants), the block normal form and its consequences (cycles.md
  Sections 3-5, 7; Cauchy-Davenport), Camion. All are marked PROVED in the repo and were read for this note; the n >= 7 step of THM-2477's all-order
  rigidity is the CITED Benhocine-Wojda theorem (publisher abstract), so Theorems A-C are PROVED given that citation.
- An independent audit (fresh reader, own code in `audit/`) found no logical error; its presentation remarks
  (K1 in F6, the complementary case of III.2, why B_q is nonempty, strongness of E, the source of the lower
  bound 2 in Part VIII.4) are incorporated.
- One finite step: the last-block lemma in the case (m, p) = (6, 3) (Section 9), settled by an exhaustive check
  over the 35 strong 6-tournaments and all subsets A (re-done independently for every m <= 9 in `audit/`).
- Theorem C (c2) is by contradiction from (H); Steps 1 and 3 use the blueness of the paths, not only their endpoint
  scores. It shows that a singleton-type black candidate exists but does not exhibit it, and it does not prove the
  "extremal predecessor" statement of cycles.md Section 7 (closed or open-path form) except in the cases listed in the
  Remark after Theorem C; both forms stay OPEN, not only for mu = 2.
- The census data at n = 10, 11 and the structured searches up to n = 16 agree with Theorems A and B in every
  instance; they are controls, not inputs to the proofs.
