# Proofs: Paley doubles and dominating sets of even graphs

Session mac-mini-2026-11-09-paleyeven. Code and outputs: `04-computation/experiments/paleyeven_20261109/`.
Canon: THM-4995 (Part A), THM-4996 (Part B), HYP-9575.

## Part A. Doubles of doubly regular tournaments

### Notation

For a tournament T on vertex set V and a set M with 2 <= |M| <= |V| - 1, the **defect** of M is

    def(M) = sum_{z not in M} min(|N+(z) cap M|, |N-(z) cap M|),

the number of arc reversals that make M a module. The **module defect** md(T) is the minimum of def(M) over all such
M. lambda(T) is the arc-strong connectivity, and delta0(T) is the least of all in- and out-degrees.

A tournament B of order q is **doubly regular** if every pair of vertices has exactly (q - 3)/4 common out-neighbours.
Then q = 3 (mod 4), B is regular of degree mu = (q - 1)/2, and every pair also has (q - 3)/4 common in-neighbours.
Paley tournaments QR_p (p prime, p = 3 mod 4; x -> y iff y - x is a nonzero square) are doubly regular.

The **double** D(B) of a tournament B on [q] has vertices a_x, b_x (x in [q]). A = {a_x} and B' = {b_x} each induce a
copy of B. Between the halves every arc is a_x -> b_y (x != y), and the q matched back arcs are b_x -> a_x. For
B = QR_p we write T_p = D(QR_p).

### Lemma A0 (md <= lambda <= delta0)

For every tournament T on n >= 3 vertices, md(T) <= lambda(T) <= delta0(T).

*Proof.* lambda <= delta0 by taking a single vertex as one side of the cut. Let U be a minimum cut, so that
e(U -> V \ U) = lambda.
- If |V \ U| >= 2, take M = V \ U. Each z in U contributes at most |N+(z) cap M|, so def(M) <= e(U -> M) = lambda.
- If V \ U = {v}, then lambda = d-(v). Take M = V \ {v}, which has n - 1 >= 2 vertices; then
  def(M) = min(d+(v), d-(v)) <= lambda.

(This is the inequality behind THM-4985's dec = md for n >= 3.) QED

### Lemma A1 (doubly regular tournaments have md = (q - 1)/2)

Let B be doubly regular of order q >= 7. Then md(B) = mu = (q - 1)/2, and def(X) >= mu for every X with
2 <= |X| <= q - 1.

*Proof.*
- **Size 2.** The distinguishers of a pair are the vertices beating exactly one of the two. There are
  (q - 2) - 2(q - 3)/4 = (q - 1)/2 of them.
- **Size q - 1.** For X = V \ {v}, def(X) = min(d+(v), d-(v)) = mu.
- **Sizes 3 to q - 2.** Let S be the skew matrix with S_wy = 1 if w -> y, -1 if y -> w, and 0 on the diagonal.
  Double regularity is equivalent to S S^T = qI - J (the off-diagonal entry is
  #agree - #disagree = (q - 3)/2 - (q - 1)/2 = -1). For X with |X| = s, put D_w = sum_{y in X} S_wy, which equals
  |N+(w) cap X| - |N-(w) cap X|. Then

      sum_w D_w^2 = 1_X^T S^T S 1_X = qs - s^2,

  so sum_{w not in X} D_w^2 <= s(q - s). By Cauchy-Schwarz, sum_{w not in X} |D_w| <= (q - s) sqrt(s). Since
  min(a, b) = (s - |a - b|)/2 when a + b = s,

      def(X) >= (q - s)(s - sqrt(s))/2 =: f(s).

  - s = 3: f = 0.634(q - 3) >= (q - 1)/2 iff q >= 10.5.
  - s = q - 2: f = q - 2 - sqrt(q - 2) >= (q - 1)/2 iff q^2 - 10q + 17 >= 0 iff q >= 7.83.
  - s = q - 3: f = (3/2)(q - 3 - sqrt(q - 3)) >= (q - 1)/2 iff 2q - 8 >= 3 sqrt(q - 3), which holds for q >= 7.
  - 4 <= s <= q - 4: s - sqrt(s) >= s/2, so f >= (q - s)s/4 >= 4(q - 4)/4 = q - 4 >= (q - 1)/2 for q >= 7.

  So every q >= 11 is covered. For q = 7 the only doubly regular tournament is QR_7, and md(QR_7) = 3 by exhaustive
  computation (`qr_mdefect.out`). QED

FINITE-EXACT check: md(QR_p) = lambda = delta0 = (p - 1)/2 for p = 7, 11, 19, 23 (`qr_mdefect.out`).

### Theorem A2 (doubles)

Let B be doubly regular of order q >= 7. Then md(D(B)) = lambda(D(B)) = delta0(D(B)) = (q + 1)/2, and D(B) is strong
and prime. Aut(B) acts diagonally on D(B). If B is self-converse (Paley: x -> -x), so is D(B), via a_x -> b_{phi(x)},
b_x -> a_{phi(x)} for an anti-automorphism phi of B.

*Proof.*
- **Upper bound.** a_x has in-degree mu + 1 (its in-neighbours in A, plus b_x) and out-degree mu + q - 1. b_x has
  out-degree mu + 1 and in-degree mu + q - 1. So delta0 = mu + 1. By Lemma A0, md <= lambda <= mu + 1.

- **Lower bound: setup.** Let 2 <= |M| <= 2q - 1. Put S = {x : a_x in M} and T = {x : b_x in M}, with s = |S| and
  t = |T|. For z not in M, write c(z) = min(d+_M(z), d-_M(z)). Then

      for a_w not in M:  d+_M(a_w) = |N+(w) cap S| + |T \ {w}|,   d-_M(a_w) = |N-(w) cap S| + [w in T];
      for b_w not in M:  d+_M(b_w) = |N+(w) cap T| + [w in S],     d-_M(b_w) = |N-(w) cap T| + |S \ {w}|.

  (Neighbourhoods N+-(w) are taken in B.) Write c_X(w) = min(|N+(w) cap X|, |N-(w) cap X|) and
  def_B(X) = sum_{w not in X} c_X(w). Every added term above is >= 0, so c(a_w) >= c_S(w) and c(b_w) >= c_T(w).
  By Lemma A1, def_B(X) >= mu whenever 2 <= |X| <= q - 1.

- **Lower bound: cases (s, t).** In each case def(M) >= mu + 1.
  1. 2 <= s, t <= q - 1: def(M) >= def_B(S) + def_B(T) >= 2mu >= mu + 1.
  2. s = 0, t >= 2: c(a_w) = min(|T \ {w}|, [w in T]) = [w in T], summing to t; the B'-part is exactly def_B(T).
     So def = t + def_B(T). This is q if t = q, and at least 2 + mu if 2 <= t <= q - 1.
  3. t = 0, s >= 2: symmetrically, def = s + def_B(S) >= mu + 1.
  4. s = 1, S = {u}:
     - (4a) t = 1, T = {v}: M is a pair.
       - u = v: its distinguishers are the a_w with u -> w (mu of them) and the b_w with w -> u (mu), so there
         are q - 1.
       - u != v: there are d+(u) + d-(v) + 2([v -> u] - [u -> v]) >= 2mu - 2 = q - 3 distinguishers.
       - Both are at least mu + 1, since q - 3 >= (q + 1)/2 for q >= 7.
     - (4b) 2 <= t <= q - 1: each w in T \ {u} has c(a_w) = min([w -> u] + t - 1, [u -> w] + 1) >= 1, and the
       B'-part is at least def_B(T). So def >= |T \ {u}| + def_B(T) >= 1 + mu.
     - (4c) t = q: for w != u, c(a_w) = min([w -> u] + q - 1, [u -> w] + 1) = [u -> w] + 1. So
       def = d+(u) + q - 1 >= mu + 1.
  5. t = 1, T = {v}, s >= 2:
     - (5a) s = q: for w != v, c(b_w) = min([w -> v] + 1, [v -> w] + q - 1) = [w -> v] + 1. So
       def = d-(v) + q - 1 >= mu + 1.
     - (5b) 2 <= s <= q - 1: the A-part is at least def_B(S). Each w in S \ {v} has
       c(b_w) = min([w -> v] + 1, [v -> w] + s - 1) >= 1. So def >= mu + |S \ {v}| >= mu + 1.
  6. s = q, 1 <= t <= q - 1: for w not in T, c(b_w) = min(|N+(w) cap T| + 1, |N-(w) cap T| + q - 1). This equals
     |N+(w) cap T| + 1, since |N+(w) cap T| <= mu <= q - 3. So def = (q - t) + e_B(T^c -> T).
     - If t <= mu, then q - t >= mu + 1.
     - Otherwise each w not in T has |N+(w) cap T| >= mu - (q - t - 1) = t - mu. With k = q - t in [1, mu], this
       gives def >= k(mu + 2 - k). This is concave in k, so it is at least its endpoint values min(mu + 1, 2mu),
       which is mu + 1.
  7. t = q, 1 <= s <= q - 1: symmetrically, c(a_w) = |N-(w) cap S| + 1 and def = (q - s) + e_B(S -> S^c) >= mu + 1.

  These cases cover every (s, t) with s + t >= 2 and (s, t) != (q, q). So md(D(B)) >= mu + 1.

- **Strong and prime.** These follow from md >= 1 (a module of size >= 2 has defect 0, and so does a non-strong
  split).
- **Diagonal action and self-converse.** Both are checked arc by arc: within A, within B', on the cross arcs and on
  the matched back arcs. QED

Case (4a) is where q >= 7 is needed: D(QR_3) has md = 0 (`paley_mdefect.out`).

FINITE-EXACT: md, twin distance, lambda and delta0 of T_7 and T_11 over all 2^14 and 2^22 vertex sets
(`paley_mdefect.out`) are 4/4/4/4 and 6/7/6/6. T_7 attains md at 35 sets (vertex complements and cross pairs
(a_x, b_y) with y - x a square); T_11 at 22 sets (vertex complements only).

### Proposition A3 (Hamiltonian paths of the doubles)

- **Exact values** (subset DP, `p7_tprop.out` and `p11_tprop.out`):
  - H(T_7) = 4,730,397 and H(T_7 - v) = 686,457 for all v.
  - H(T_11) = 17,667,825,987,845 and H(T_11 - v) = 1,822,800,617,515 for all v.
  - The deletions are equal because Aut(QR_p) is transitive and D(QR_p) is self-converse.
  - H(T_p)/H(QR_p)^2 = 132.4 and 1953.7: the matched back arcs multiply the ordinal-sum count heavily.
- **Brégman bound.** Adding a vertex w with w <-> v for all v turns Hamiltonian paths of a digraph into Hamiltonian
  cycles through w. Brégman's permanent bound then gives

      H(D) <= (n!)^(1/n) * prod_v ((d+(v) + 1)!)^(1/(d+(v) + 1)).

  For D(B) the out-degrees are (3q - 3)/2 on A and (q + 1)/2 on B'. A deletion only lowers out-degrees, and deleting a
  b-vertex leaves the largest product.

### Theorem A4 (conditional escapers of unbounded module defect)

For c > 0 let **(S_c)** be the hypothesis that S(m) contains every odd h <= c m!/2^(m - 1) other than 7 and 21, for
all m >= m_0. (The exact censuses show S(10) and S(11) are solid up to 2.07 and 2.06 times m!/2^(m - 1); klein,
THM-4889.)

Under (S_1), D(B) is an escaper (THM-4985) for every doubly regular B of order q >= 55 with 2q - 2 >= m_0.
Under (S_c) for any fixed c > 0, the same holds for all q >= q_0(c). Hence escapers of unbounded decomposability defect
exist under (S_c): HYP-9565 (A) holds conditionally.

*Proof.* D(B) is strong, prime and has lambda = (q + 1)/2 >= 3 (Theorem A2). It remains to show
H(D(B)) <= (2q - 1)!/2^(2q - 2) and H(D(B) - v) <= (2q - 2)!/2^(2q - 3) for every v.
- For 55 <= q <= 20,000 (q = 3 mod 4), both log-ratios of Brégman's bound to these means are negative and decrease
  (exact log-gamma evaluation, `bregman_table.out`). In log10 they are -1.053 and -0.388 at q = 55, -1.461 and -0.791
  at q = 59, and -2.718 and -2.039 at q = 71. q = 55 is the first q with both negative: at q = 51 the deletion ratio
  is +0.008.
- For large q, Stirling's formula with Robbins' bounds gives log(Brégman / E) = q log(3/4) + O(log q). The same holds
  for deletions. Hence they are negative for all q beyond the verified range.

For (S_c) the margin q log(4/3) - O(log q) eventually exceeds log(1/c). QED

FINITE-EXACT: H(T_11)/E[H(21)] = 0.3626 and H(T_11 - v)/E[H(20)] = 0.3928. So T_11 is an escaper of md 6 already
under (S_0.4) at orders 20 and 21. T_7 is not: H(T_7) = 4,730,397 exceeds 3,719,831, the best known order-13 count
(THM-329), and is 3.11 E[H(13)]. Its deletions, 686,457 = 2.93 E[H(12)], exceed 531,205, the best known order-12 count.
So T_7 is a first creator unless max H(13) >= 4,730,397.

### Theorem A5 (a certified escaper of decomposability defect 5, order 18)

Let C = C_9{1,2,3,5}, the circulant tournament on Z_9 with x -> x + s for s in {1, 2, 3, 5}, and let D = D(C) (order
18).
- **Structure.** md(D) = lambda(D) = delta0(D) = 5, over all 2^18 sets (`d9_mdefect.out`). D is strong and prime.
- **Counts.** H(D) = 6,450,103,269 and H(D - v) = 773,858,997 for all 18 v (`d9_tprop.out`).
- **Certificates.**
  - 773,858,997 is in S(16): `wit_d9_del_m16.out` gives an explicit 16-tournament, re-counted by a direct subset DP
    inside `wit.c`.
  - 6,450,103,269 is in S(17): `wit_d9_H_m17.out` gives an explicit 17-tournament, found by `./wit 17 0.5 4242` and
    re-counted by wit.c's DP and by tprop.c.
  - Hence D is an escaper of order 18 with decomposability defect 5. The previous record was 3 (THM-4985 (2)).
  - Housekeeping: an earlier run (seed 1709) ended, probably with a witness, but its output was lost. The file it was
    writing had been committed, and a `git pull --rebase` replaced it with a new inode.

C is not doubly regular: md(C) = 3 = twin distance. So Theorem A2 does not apply, and md(D) = 5 = 1 + delta0(C) is
the computed value. The other Z_9 circulants with 1 in S give md(D) <= 5. S = {1, 2, 3, 4} gives md(D) = 3, because a
twin distance of 1 in C contributes 1 + 2.

## Part B. Dominating sets of Royle-even graphs

### Notation

- d(G) is the number of vertex sets S that dominate G (every vertex is in S or adjacent to S). It is always odd
  (Brouwer 2009).
- Royle parity: fix the reference orientation u -> v for u < v. An automorphism g reverses {u, v} iff
  g(u) > g(v). sgn_G(g) = (-1)^{#reversed edges} is independent of the labelling and is a homomorphism on Aut(G).
- G is **even** iff sgn_G is trivial (Royle-Praeger-Glasby-Freedman-Devillers 2023; THM-4531).
- E_r denotes the empty graph on r vertices, and v denotes join.

### Lemma B1 (completeness and conciseness of d over all graphs)

Every odd k is d(G) for a graph G on ceil(log2(k + 1)) vertices.

*Proof.* d(K_1 v G) = 2^|G| + d(G), since every set containing the apex dominates and a set without it must dominate
G (and is then nonempty). Padding with isolated vertices does not change d. Induct on N:
- every odd k < 2^N is d of a graph on N vertices;
- if 2^(N - 1) < k < 2^N, write k = 2^(N - 1) + x with x odd and x < 2^(N - 1), and cone over a graph on N - 1 vertices
  with d = x.

These are the threshold graphs. (Elementary; very possibly folklore.) QED

### Lemma B2 (parity calculus)

1. **Twins.** Swapping non-adjacent twins (same neighbourhood) is even. Swapping adjacent twins (same closed
   neighbourhood) is odd: the edge between them reverses, and every other edge pair reverses 0 or 2 times. So a graph
   with adjacent twins is odd.
2. **Unions.** G is even iff every component is even. Aut(G) is generated by component automorphisms and swaps of
   isomorphic components. A swap is even: each edge e and its image phi(e) are exchanged with 0 or 2 reversals in
   total.
3. **Cones.** If G has >= 2 vertices and no universal vertex, then K_1 v G is even iff G is. The apex is then the
   unique universal vertex, so Aut(K_1 v G) = Aut(G). Labelling the apex 0, no apex edge is reversed. In particular
   K_1 v (G u E_j), j >= 1, with G u E_j of order >= 2, is even iff G is.
4. **Joins.** Let A and G be graphs whose complements are connected (for instance A = A' u E_i, G = G' u E_j,
   i, j >= 1) and with A not isomorphic to G. Then A v G is even iff A and G are. The complement of A v G is the
   disjoint union of the two connected complements, so Aut = Aut(A) x Aut(G). Labelling A first, the cross edges keep
   their orientation.
   - A swap of isomorphic join-children of size s reverses all s^2 cross edges between them, so it has sign (-1)^s.
   - Example: K_{r,r} = E_r v E_r is odd iff r is odd.
5. **Complements.** sgn_G(g) sgn_{complement(G)}(g) = sgn(g), the permutation sign: the inversion count of g has the
   parity of g.
6. **Trees.** A tree is odd iff it is bicentral and some automorphism exchanges the two central vertices. Every
   automorphism fixes the centre, and is a product of swaps of isomorphic branches (even by item 2's argument,
   labelling ancestors first) and, in the bicentral case, possibly the flip. The flip reverses the central edge once
   and the other edges in pairs.

### Lemma B3 (exhaustiveness of the census)

If v has at least one neighbour in G, then d(G + v) > d(G). The map S -> S cup {v} injects dominating sets of G into
those of G + v, and V(G) dominates G + v without containing v. Every connected graph has a vertex order whose prefixes
are connected (breadth-first search). So the level-by-level extension in `dom_even.c` reaches every connected graph
with d <= K, provided it keeps exactly the connected children with d <= K. The child counts use

    d(G + v_N) = A(N) + d(G) - C(N),

where:
- A(N) = #{S subset of V : V \ dom(S) subset of N} counts the sets containing v (v dominates N);
- C(N) = #{S dominating G : S cap N = empty} is subtracted from the count of sets without v that must still hit N.

Every stored graph's d is recomputed directly (`count_dom`), and parity is verified to be odd.

### Theorem B4 (exact exceptional set below 512)

Run `dom_even 511 16` (`dom_K511.err`, `dom_K511_summary.txt`):
- **Census.** All connected graphs with d <= 511 number 1, 1, 2, 6, 21, 112, 853, 11,117, 261,080 (every connected
  graph on <= 9 vertices; A001349), then 199,053 on 10 vertices, 1,855 on 11 and none on 12.
- **Check on the parity code.** The Euler transform of the connected-even counts 1, 0, 1, 2, 8, 43, 398, 6,413,
  184,596 (A334335, as in THM-4721) is 1, 1, 2, 4, 12, 56, 456, 6,880, 191,536 = A000568. This is the
  Royle-Praeger-Glasby-Freedman-Devillers equinumerosity, reproduced class by class through order 9. It independently
  validates the parity code: even = all components even (B2.2).
- **Exceptional set.** Since d is multiplicative over components and evenness is componentwise:
  - the even d-values below 512 are all odd numbers except 3, 7 and 13;
  - all graphs give every odd number;
  - connected even graphs miss exactly 3, 7, 13 and 25, and 25 = d(P_3 u P_3).
- **Realisers of the exceptions.** The only graphs, up to isolated vertices, with d in {3, 7, 13} are K_2, K_3 and
  K_4 - e. Each has adjacent twins (B2.1).
- **The value 21** is even: 21 = d(K_1 v (P_3 u K_1)) = 16 + 5. Its only nontrivial automorphism swaps the
  non-adjacent twins.
- **Odd graphs (added during audit P1).** For the complementary class, a graph is odd iff some component is odd.
  Below 512 the d-values of Royle-odd graphs are all odd numbers except 1, 5 and 37.
  - 1 is realised only by empty graphs;
  - 5 only by P_3 (plus isolated vertices);
  - 37 is prime, and its four connected realisers are all even.
  So the even/odd split of the d-spectrum has exceptional sets {3, 7, 13} and {1, 5, 37}.

### Proposition B5 (constructive coverage, evidence only)

`dom_cover.c` starts from the exhaustive base and closes under unions (B2.2), cones K_1 v (G u E_j) (B2.3), E_r-joins
and general joins (A' u E_i) v (G' u E_j) (B2.4). Up to 2^22:
- the general joins reach every odd value except 3, 7, 13 and 12,326 others (`dom_cover_K511_22_genjoin.out`);
- without general joins, 18,766 others are missed (`dom_cover_K511_22_basic.out`).

The unreached values cluster at 2^m + {3, 7, 13} and just below 2^m. Cones cannot produce 2^m + {3, 7, 13}, because
3, 7 and 13 are not even values. These are limits of the chosen families, not counterexamples.

### Proposition B6 (the OCF statistic on even graphs, small range)

I(G, 2) = sum over independent sets of 2^|S| is the OCF value: H(T) = I(Omega(T), 2). It also strictly increases
under adding any vertex. Every connected graph with I(G, 2) <= 31 has at most 15 vertices, and an exhaustive census
(`i2_K31.out`) gives the following.
- I(G, 2) on even graphs takes, below 32, exactly the values 1, 3, 9, 11, 17, 19, 25, 27, 29.
- It misses 5, 7, 13, 15, 21, 23, 31: both 7 and 21, but many more.
- Over all graphs every odd value occurs (I(K_m, 2) = 2m + 1).
