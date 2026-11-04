# The non-critical run of a uniquely Hamiltonian tournament: local roundness, one run from a sole entry to a sole exit, cut arcs + non-critical vertices = n + 1, a Shelor-free proof of the degree-1 theorem, the admissible Shelor parents (Garey's count re-derived), and the run-length recurrence (a second, Shelor-free proof of HYP-9458)

**Session:** kind-pasteur-2026-11-04-S10. Worktree `.claude/worktrees/uh-one-arc-join-7beb92`.

**Task (from the dispatching session).** Prove or settle "every uniquely Hamiltonian tournament is a one-arc join"
(HYP-9457 (b), opus-2026-11-04-S73), check concurrency and the literature, run an exact census, and finish with a canon
theorem, a counterexample, or a partial theorem plus a hypothesis.

**What happened to the target.** While this session worked, the statement was settled elsewhere and pushed first:
- the S73 audit typed HYP-9457 (b) KNOWN (Shelor 1977's in-degree-1 vertex is a single-vertex one-arc cut;
  MISTAKE-705);
- opus-2026-11-04-S74 recorded it as **THM-4867** (KNOWN, with the equivalences, the extension lemma, an independent
  generator through order 16, three failed proof routes, and the cut-arc table T(n, k) of HYP-9458);
- klein cross-cited both.

This note defers to THM-4867 for the statement and does not re-claim it. It adds a structure theorem for the
non-critical vertices, which yields a short proof of the degree-1 theorem that does not use Shelor's insertion grammar,
and it adds several exact consequences.

**Status.**
- **PROVED (elementary) here.** THM-4868 (A)-(H), section 1:
  - local roundness of non-critical vertices and its propagation along the cycle;
  - the degree-1 theorem (KNOWN: Shelor 1977, THM-4867) by a route new to the repository, which succeeds where
    THM-4867's failed route (c) (non-critical vertices) stopped;
  - the non-critical vertices form one run from a sole entry to a sole exit;
  - cut arcs + non-critical vertices = n + 1, which turns HYP-9458's cut-arc table into the run-length table;
  - asymmetry;
  - the refined Fibonacci counts by numbers of sole exits and sole entries;
  - (G) the admissible Shelor parents are y, y- and (when y- is a sole entry) y--. This is a Shelor-free proof of
    Shelor's grammar (KNOWN: Shelor 1977, recorded in `operation_closed_synthesis_20261029` section II), and it re-derives
    Garey's count F(2n - 6) (KNOWN: Garey 1972, Moon 1982);
  - (H) the run transitions r, r + 1, r + 2. They give a Shelor-free proof of the run-length recurrence (HYP-9460 P7)
    and of every formula of opus-S74's HYP-9458.
- **Priority for HYP-9458 and P7.** HYP-9458's closed forms and generating function were PROVED first by opus-S74's
  independent audit (`a0518f53cd`, pushed 12:51), given Shelor's grammar (CITED). Their proof is the count law
  k = 3 + #{b, e} - #{c} along the Shelor word, and it noted that P7 follows from it together with this session's (C).
  - (G)-(H) here are a second, independent proof that does not use Shelor's grammar: (G) proves the needed part of the
    grammar from scratch.
  - Dictionary: Shelor's parent p = pred^j(newest vertex), j = 0, 1, 2, is (G)'s admissible set {y_U, y_U-, y_U--};
    j = 2 occurs only when y_U- is a sole entry.
    The count law's Delta k = +1, 0, -1 is (H)'s Delta r = 0, +1, +2 under k = n + 1 - r.
- **FINITE-EXACT.**
  - On all F(2n - 6) classes at orders 4-18: (A), (C), (D), (E), the degree-1 and (F) counts, and one-arc cuts. Orders
    14-18 were run after the pre-registration commit `a8db6218f9` (HYP-9460 P1-P8).
  - (G) and (H) on every parent of order 4-14.
  - Not checked by these programs: the (F) bijection, the rooted count, the identification of the r = 2 class with A_n,
    and the "equivalently" form of (B). The independent audit checked them at orders 4-16 (the strong form at 5-9).
- **Pre-registered, then proved twice (HYP-9460 P7).** The run-length recurrence
  N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1).
  - It was conjectured from orders 6-13 and pre-registered with blind predictions at orders 17-18.
  - opus-S74's audit then proved it via HYP-9458 and (C), given Shelor's grammar. Independently, this session proved
    it as (H) during the confirmatory run, from the admissible-parent structure (G). (G) and (H) were
    first observed on all parents through order 12, then proved; [`shelor_parents_20261104.c`](../../04-computation/experiments/shelor_parents_20261104.c) checks them on every parent through order 14 ([out](../../04-computation/experiments/shelor_parents_20261104.out)).
- **Novelty.** Not claimed against the literature. Douglas 1970, Shelor 1977 and Moon 1982 characterize or construct
  this class. Their bodies could not be read here: paywall, a ScienceDirect CAPTCHA (not bypassed), and a Cloudflare
  challenge on the Shelor mirror. The arguments below are new to the repository.

Independently audited: section 9.

## 0. Inheritance

- **Closest proved mechanism.** THM-4867's extension lemma: a non-critical vertex of a uniquely Hamiltonian tournament
  has a cyclic-interval out-neighbourhood. klein's block lemma (THM-4862 (f)) is the same insertion fact, and THM-4866's
  cut proof for interval chains is the sibling.
- **Canonical hostile example.** RT_5 = Circ(Z_5, {1, 2}). It is round, has minimum out-degree 2, and has exactly two
  Hamiltonian cycles, the pentagon and the pentagram (HYP-9459). It is the boundary of (B): the second cycle built in
  section 2, proof of (B), step 3, is RT_5's pentagram.
- **Corrected near miss.** HYP-9457 (b) was posed as open and is classical (MISTAKE-705). THM-4867's route (c), "every
  non-critical vertex has in- or out-degree 1", is false from order 5. The true statement is (A) below: non-critical
  vertices are locally round, and only the two ends of their run have degree 1.
- **Least-used sidecar.** The non-critical ("neutral") vertices of strong tournaments: Moon 1966; Las Vergnas 1975, as
  stated in Burzio-Demaria 1988 (read here); Savchenko's papers (arXiv:2311.13377, read here). They have not been used
  in the repository's Hamiltonian-cycle work.

## 1. Statements (THM-4868)

Let T be a tournament of order n >= 4 with exactly one Hamiltonian cycle C. Write v+ and v- for the successor and
predecessor on C, d+ and d- for out- and in-degree, and NC for the set of non-critical vertices (T - v strong). A
**sole exit** has d+ = 1 and a **sole entry** has d- = 1.

- **(A) Local roundness and propagation.**
  - If v is non-critical, then N+(v) = {v+, v++, ..., v^(+d+(v))}: v beats exactly the d+(v) vertices that follow it on
    C, and is beaten by the d-(v) vertices that precede it.
  - If moreover d+(v) >= 2, then v+ is non-critical. If d-(v) >= 2, then v- is non-critical.
- **(B) The degree-1 theorem, by roundness** (KNOWN: Shelor 1977; THM-4867). T has a sole exit and a sole entry.
  Equivalently, every strong tournament with minimum out-degree >= 2 (or minimum in-degree >= 2) has at least two
  Hamiltonian cycles.
  - Two is attained by RT_5 at order 5.
  - klein's exhaustive table (HYP-9459, orders 5-10, FINITE-EXACT) gives minimum 3 at orders 6-10, attained by
    F(2n - 9) tournaments.
  - Raising the bound to 3 for all n >= 6 is HYP-9459 (1)'s open two-cycle half.
- **(C) One run.**
  - NC is a single run y = c_1 -> c_2 -> ... -> c_r = x of consecutive vertices of C, with 2 <= r <= n - 1.
  - Its first vertex y is a sole entry and its last vertex x is a sole exit.
  - For every sole exit z, z+ is critical; for every sole entry z, z- is critical.
  - y is the only non-critical sole entry and x the only non-critical sole exit. So T is Shelor's insertion of a
    vertex (y) into the uniquely Hamiltonian tournament T - y with parent y-, and y is the only vertex whose deletion is
    such a step.
- **(D) Cut arcs.**
  - An arc of C is a cut arc (T minus the arc is not strong) iff it is not internal to the run.
  - The cut arcs form the path x -> x+ -> ... -> y- -> y of C, so #cut arcs + |NC| = n + 1.
  - Hence HYP-9458's T(n, k) equals N(n, n + 1 - k), where N(n, r) counts classes with r non-critical vertices.
  - With Las Vergnas's theorem (the bineutral tournament A_n is the only Hamiltonian tournament with exactly two
    non-critical vertices; CITED via Burzio-Demaria 1988), this gives HYP-9458's top diagonal T(n, n - 1) = 1. (H)
    proves the same independently, since N(n, 2) = N(n - 1, 2).
- **(E) Asymmetry.** T has no non-trivial automorphism.
- **(F) Refined counts** (n >= 5; the total F(2n - 6) is proved in (G)).
  - The classes with 1 or 2 sole exits number F(2n - 8) and F(2n - 7).
  - By (number of sole exits, number of sole entries):

    | pair | (1, 1) | (1, 2) | (2, 1) | (2, 2) |
    |---|---|---|---|---|
    | number of classes | F(2n - 10) + 1 | F(2n - 9) - 1 | F(2n - 9) - 1 | F(2n - 8) + 1 |

  - There are F(2n - 5) rooted pairs (T, sole exit).
  - The classes with two sole exits are in bijection with rooted pairs (T', sole exit) of order n - 1.

- **(G) Admissible parents and Garey's count.**
  - For a tournament U and a vertex p, let U^p be U plus a new vertex y with p -> y and y -> u for every other u. Call
    p admissible in U when U^p is uniquely Hamiltonian, equivalently when e_p(U) = 1 (Lemma 3).
  - For T as above with run start y, the admissible vertices are exactly y and y-, together with y-- when y- is a sole
    entry. So there are 2 or 3 of them, and 3 exactly when T has two sole entries; those are then y- and y.
  - With (C) and (E), T ↦ (T - y, y-) is a bijection from the classes of order m + 1 onto the pairs (U, p) with U of
    order m and p admissible in U.
  - Hence |UH(m)| = F(2m - 6). This is Garey's theorem (KNOWN: Garey 1972, Moon 1982), by a route that is self-contained
    except for Moon 1966 and Camion 1959. So (F) does not need Garey's citation.
  - **Typing.** The statement of (G) is Shelor's grammar (KNOWN: Shelor 1977), as recorded in
    `operation_closed_synthesis_20261029/THEORY.md` section II:
    - type A, one sole entry a = y: parents a (letter b) and pred(a) (letter a);
    - type B, two sole entries a0 = y- and a1 = y: parents a1, a0 and pred(a0) (letters e, d, c);
    - child type B exactly when the parent is a sole entry;
    - counts A' = A + B, B' = A + 2B.

    What is new to the repository is the proof from the non-critical run, which does not use Shelor's paper.
- **(H) Run transitions and the run-length recurrence.**
  - The child U^p has run length r(U), r(U) + 1 or r(U) + 2 according as p = y_U, y_U- or y_U--.
  - Hence N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1) for n >= 6, with the generating function of
    section 4.
  - By (D), HYP-9458's table satisfies T(n, k) = T(n - 1, k - 1) + 2 T(n - 1, k) - T(n - 2, k - 1). All of HYP-9458's
    fitted forms follow: 2^(n-5), n 2^(n-7), (n^2 + n - 8) 2^(n-10), and the diagonals 1, n - 4, C(n - 3, 2).

## 2. Proofs

**Lemma 1 (insertion; THM-4867's extension lemma, klein's block lemma).** Let v be non-critical.
1. T - v has a Hamiltonian cycle (Camion).
2. v has an in-neighbour and an out-neighbour on every Hamiltonian cycle Z' of T - v, since T is strong. So v can be
   inserted between some consecutive u -> w of Z' with u -> v -> w.
3. Distinct (Z', insertion point) give distinct Hamiltonian cycles of T: deleting v and closing u -> w recovers Z'.
4. Since hc(T) = 1, T - v has exactly one Hamiltonian cycle C', and v has exactly one insertion point on it.
5. The insertion points are the places where the cyclic I/O word of v along C' (I = beats v, O = beaten by v) reads
   I then O. So there is exactly one I-run: N-(v) and N+(v) are complementary cyclic intervals of C'.
6. Inserting v gives C. So C' = C - v + (v- -> v+), v- is the last vertex of N-(v), and v+ is the first vertex of
   N+(v).

**Proof of (A).**
- By Lemma 1, N+(v) is the interval of C' = C - v that starts at v+. In the order of C this is
  {v+, ..., v^(+d+(v))}.
- If d+(v) >= 2, then v -> v++. The cycle v -> v++ -> ... -> v- -> v avoids exactly v+, so T - v+ is Hamiltonian,
  hence strong.
- The statement for d-(v) >= 2 is the converse statement (T^op is uniquely Hamiltonian with the reversed cycle). QED.

**Lemma 2.** If z is a sole exit, then z+ is critical: in T - z+ the vertex z has out-degree 0. Dually, if z is a sole
entry, then z- is critical.

**Proof of (B).** Suppose every vertex has out-degree >= 2.
1. By Moon (1966; a strong tournament of order n >= 4 has a strong subtournament of order n - 1), NC is non-empty.
2. By (A), NC is closed under v -> v+, so NC = V. By (A) again, every u beats exactly the k(u) = d+(u) >= 2 vertices
   after it on C. That is, T is the round tournament with step set {1, ..., k(u)} at u. In particular u -> u^(+s)
   whenever s <= k(u).
3. **n odd (n >= 5).** The walk with constant step +2 is a Hamiltonian cycle, because 2 is invertible mod n. It uses
   chords u -> u++, so it is not C.
4. **n even.** First, n >= 6: at n = 4 the out-degrees sum to 6 < 8. The out-degrees sum to n(n - 1)/2 > 2n, so some w
   has k(w) >= 3. Index C from u_0 = w+++. The cycle
   - u_0 -> u_2 -> ... -> u_(n-2) -> u_(n-1) -> u_1 -> u_3 -> ... -> u_(n-3) -> u_0

   uses steps +2, one +1 (u_(n-2) -> u_(n-1)) and one +3 (u_(n-3) = w -> u_0). It visits every vertex once: the even
   indices, then n - 1 and the odd indices 1, ..., n - 3. It contains +2 chords, so it is not C.

Either way hc(T) >= 2, a contradiction. A sole entry exists by the converse. The "equivalently" form: a strong
tournament has hc >= 1 (Camion), and minimum out-degree >= 2 rules out hc = 1. QED.

**Proof of (C).**
1. NC is non-empty (Moon). It is not all of V: a sole exit x exists by (B), and x+ is critical by Lemma 2.
2. Each maximal run of NC is therefore a proper arc of C. Its last vertex z has z+ critical, so d+(z) = 1 by (A), and
   z is a sole exit. Similarly its first vertex is a sole entry.
3. There are at most two sole exits, and two are adjacent: the arc between them must be the unique out-arc of its
   tail, so they are some z and z+. Then z+ is critical by Lemma 2. Three are impossible: they would form a 3-cycle
   with no arc leaving it, and n >= 4.
4. Hence at most one sole exit is non-critical and ends a run, so there is exactly one run.
5. Its ends differ, because no vertex is both a sole exit and a sole entry when n - 1 >= 3. So 2 <= r <= n - 1.
6. A non-critical sole entry y' has y'- critical (Lemma 2), so it starts a run. Hence y' = y, and dually for exits.
7. T - y is uniquely Hamiltonian by Lemma 1, and y beats everything in T - y except y-. That is Shelor's insertion
   with parent y-. QED.

**Proof of (D).**

*Internal arcs of the run are not cut arcs.* Let e = u -> u+ with u, u+ both in the run.
- By Lemma 2, d+(u) >= 2 (otherwise u+ would be critical) and d-(u+) >= 2 (otherwise u would be critical).
- By (A), u -> u++ and u- -> u+.
- T - e then contains the C-path u+ -> ... -> u and the path u -> u++ -> ... -> u- -> u+. So T - e is strong.

*Both C-arcs at a critical vertex w are cut arcs.*
- Write T - w = Q_1 => ... => Q_m (m >= 2). The C-path w+ -> ... -> w- runs through Q_1, ..., Q_m in order, each as a
  contiguous block, so w+ is in Q_1 and w- is in Q_m.
- In T only w can enter Q_1. Suppose w -> z for some z in Q_1 other than w+. Then Q_1 is strong, so it has a Hamiltonian
  path from z to some q. The cycle w -> z => q -> (rest of the C-path from Q_2 on) -> w is Hamiltonian, and w's
  successor on it is not w+, so it is not C. This is impossible.
- So w -> w+ is the only arc entering Q_1, and it is a cut arc. The arc w- -> w is the converse case, using Q_m.

*The count.* The critical vertices are exactly the vertices x+, ..., y- of C between the run's ends (by (C)). The arcs
incident to them are x -> x+, ..., y- -> y, which number n - r + 1. Every other arc is internal to the run.

*Top diagonal.* T(n, n - 1) = N(n, 2), and N(n, 2) <= 1 by Las Vergnas. Equality holds because A_n is uniquely
Hamiltonian:
- A_n has vertices x_1, ..., x_n and arcs x_i -> x_j iff j = i + 1 or j < i - 1.
- The only arc crossing the cut between x_i and x_(i+1) forward is x_i -> x_(i+1). A Hamiltonian cycle crosses every cut
  forward, so it uses the whole path x_1 -> ... -> x_n and closes with x_n -> x_1.
- This is THM-4866's cut argument for interval chains. QED.

**Proof of (E).**
- An automorphism maps C to itself, so it acts as a rotation of C. It permutes the set of sole exits, which is {x} or
  {x, x+} by (B) and step 3 of the proof of (C).
- A non-trivial rotation fixes no vertex. It swaps x and x+ only if it shifts by 1 and by -1 at once, that is, only if
  n <= 2. QED.

**Proof of (F).**

*The bijection.* Let T have two sole exits x and y = x+.
- T - x contains the cycle y -> ... -> x- -> y, because x- is not y+ and so beats the sole exit y. It is therefore
  uniquely Hamiltonian (Lemma 1), and y is still a sole exit of it.
- Conversely, let T' (order >= 3) be uniquely Hamiltonian with a sole exit y. Add x with x -> y and every other vertex
  -> x. Then:
  - hc(T' + x) = s_y(T'), by klein's THM-4862 (b) "C(U + a_exit) = s_r(U)";
  - s_y(T') = s_(y+)(T' - y) = hc(T') = 1.
- Both maps are canonical (x is the sole exit whose out-neighbour is a sole exit), so they are inverse bijections of
  isomorphism classes.

*The counts.* Let a(n) and b(n) be the numbers of classes with one and with two sole exits.
1. Rooted classes number a + 2b at order >= 4: the two sole exits of a class are distinguishable, since only x has a
   sole exit as out-neighbour, and (E) gives the same. With b(4) = 1, b(n) = a(n - 1) + 2 b(n - 1).
2. With a + b = F(2n - 6) (proved in (G), whose proof does not use (F); KNOWN as Garey's theorem),
   b(n) = b(n - 1) + F(2n - 8) = F(2n - 7), and a(n) = F(2n - 8).
3. The bijection preserves the number of sole entries (n >= 5):
   - x has in-degree n - 2 >= 3;
   - y only gains an in-arc, and was not a sole entry, being a sole exit of order >= 4;
   - every other vertex is unchanged.
4. So (2, j)(n) = (1, j)(n - 1) + 2 (2, j)(n - 1). The converse gives (1, 2) + (2, 2) = F(2n - 7).
5. Hence, from the base (2, 2)(4) = 1 (S_4), (2, 2)(n) = (2, 2)(n - 1) + F(2n - 9) = F(2n - 8) + 1. The table follows, and the rooted count is
   a + 2b = F(2n - 5). QED.

**Lemma 3 (entry attachment; klein's THM-4862 (b), C(U + a_entry) = e_r(U)).** hc(U^p) = e_p(U).
- A Hamiltonian cycle enters y from p and leaves to some z other than p.
- The rest of the cycle is a Hamiltonian path of U from z to p.
- Conversely, every Hamiltonian path of U ending at p starts at some z other than p, which y beats.

**Proof of (G).** Write adm(U) = {p : e_p(U) = 1}. Every sole entry z of U lies in adm(U):
- every path ending at z enters it from z-, so e_z(U) = e_(z-)(U - z);
- U = (U - z)^(z-), so Lemma 3 gives e_(z-)(U - z) = hc(U) = 1.

*Base, order 4.* The strong 4-tournament is unique (score sequence (1, 1, 2, 2)) and has exactly one Hamiltonian cycle.
So S_4 is the only uniquely Hamiltonian class of order 4, and u(4) = 1. S_4 has cycle 0 -> 1 -> 2 -> 3 -> 0 and chords
0 -> 2, 1 -> 3; its sole entries are 0 and 1, so b(4) = 1.
- Its Hamiltonian paths end at 3, 0, 2, 1, 2, so e = (1, 1, 2, 1).
- The run is 1 -> 2 and y = 1. The vertex y- = 0 is a sole entry.
- adm = {1, 0, 3} = {y, y-, y--}.

*Step, order m >= 5.* Let y be T's run start, U = T - y and p_0 = y-.
- U is uniquely Hamiltonian of order m - 1 >= 4 (Lemma 1), its cycle is C_U = C_T - y, and T = U^(p_0). So
  e_(p_0)(U) = 1.
- The only in-neighbour of y is p_0. Classify the Hamiltonian paths of T ending at p by the position of y.

The cases:
- **p = y.** The paths arrive from p_0, so e_y(T) = e_(p_0)(U) = 1.
- **p = p_0.** y must come first. So e_(p_0)(T) = e_(p_0)(U) = 1, and y- is admissible.
- **p not in {y, p_0}.** e_p(T) = [e_p(U) - hp_U(p_0, p)] + W(p_0, p), where:
  - hp counts Hamiltonian paths of U with the given ends; the bracket counts the paths with y first;
  - W counts ordered pairs (A, B) of disjoint paths of U covering V(U), with A ending at p_0 and B ending at p. These
    give the paths A, y, B with y not first.

  Let P_p = C_U - (p -> p+) be the Hamiltonian path of U from p+ to p. It passes p_0, and splitting it after p_0 gives a
  pair counted by W. So W >= 1.
  - **p+ is not p_0** (p is not y--). P_p is counted in the bracket, so e_p(T) >= 2.
  - **p = p_0-, p_0 a sole entry** (of U, equivalently of T, since y does not beat p_0). A path ending at p_0- must start
    p_0, y, followed by a Hamiltonian path of U - p_0 ending at p_0-. Since U = (U - p_0)^(p_0-), Lemma 3 gives
    e_p(T) = hc(U) = 1.
  - **p = p_0-, p_0 not a sole entry.**
    - By induction, the admissible non-sole-entries of U are y_U- (if y_U- is not a sole entry) or y_U-- (if it is). So
      p_0 is one of these, and p = p_0- is y_U-- or y_U---. This is not admissible in U (|U| >= 4). Since e_p(U) >= 1 (the path P_p), it follows that e_p(U) >= 2.
    - Let h = hp_U(p_0, p). Here p+ = p_0, so P_p is one of these paths and h >= 1. Splitting each of the h paths
      after p_0 gives distinct pairs for W, so W >= h.
    - Then e_p(T) >= (e_p(U) - h) + h = e_p(U) >= 2.

So adm(T) = {y, y-} together with y-- iff y- is a sole entry. QED for the first claim.

*The bijection.*
- If p is admissible in U (order m >= 4), then T = U^p is uniquely Hamiltonian, and its new vertex y is a non-critical
  sole entry (T - y = U is strong). By (C), y is T's unique non-critical sole entry, so T - y_T = U and y_T- = p.
- Conversely, every T of order m + 1 >= 5 equals (T - y_T)^(y_T-), with y_T- admissible (Lemma 3).
- By (E), pairs (U, p) of a class and a vertex are distinct.

*Sole entries.* T = U^p has two sole entries iff p is a sole entry of U:
- y is a sole entry;
- p keeps its in-arcs (y does not beat p);
- every other vertex gains the in-arc from y, so loses the property.

*The count.* Let u(m) = |UH(m)| and b(m) be the number of classes with two sole entries.
- u(m + 1) = sum over U of |adm(U)| = 2 u(m) + b(m).
- b(m + 1) = sum over U of #(sole entries of U) = u(m) + b(m).
- With u(4) = b(4) = 1, this gives u(m + 1) = 3 u(m) - u(m - 1), so u(m) = F(2m - 6) and b(m) = F(2m - 7). QED.

**Proof of (H).** Let U have run [y_U, ..., x_U] of length r, let p be admissible, and T = U^p.
- C_T is C_U with y inserted after p. The unique Hamiltonian cycle of T is y -> (the unique Hamiltonian path of U ending
  at p) -> y, and that path is C_U - (p -> p+), because e_p(U) = 1.
- By (C), the run of T starts at y.

Three facts:
- Non-critical vertices v of U other than p stay non-critical: in T - v, the vertex y has the in-neighbour p and some
  out-neighbour.
- p is critical in T (Lemma 2).
- x_U is still a sole exit of T with successor x_U+. It is distinct from p: y_U is not x_U since r >= 2, y_U- is
  critical, and y_U-- is critical when y_U- is a sole entry (Lemma 2). So x_U+ is critical in T (Lemma 2).

The run of T is therefore the segment of C_T from y up to x_U:
- **p = y_U:** y, y_U+, ..., x_U, of length r.
- **p = y_U-:** y, y_U, ..., x_U, of length r + 1.
- **p = y_U-- (y_U- a sole entry):** y, y_U-, y_U, ..., x_U, of length r + 2. Here y_U- is non-critical in T, because
  T - y_U- has the Hamiltonian cycle y -> y_U -> (C_U-path) -> y_U-- = p -> y. QED.

*The recurrence.* Let N_2(m, r) count the classes with run length r and two sole entries. By the bijection, (H), and the
sole-entry rule, for m >= 4:

    N(m + 1, r)   = N(m, r) + N(m, r - 1) + N_2(m, r - 2),
    N_2(m + 1, r) = N(m, r) + N_2(m, r - 1).

Use the second line at m = n - 2 and the first at m = n - 2 to get N_2(n - 1, r - 2) = N(n - 1, r - 1) - N(n - 2, r - 1).
Substituting into the first line at m = n - 1 gives the recurrence for n >= 6. The closed forms of HYP-9458 follow by
induction on n from the recurrence and (D):
- T(n, 2) = 2 T(n - 1, 2);
- T(n, 3) = 2^(n-7) + 2 T(n - 1, 3) for n >= 7, from the base T(6, 3) = 3, so T(n, 3) = n 2^(n-7) for n >= 6;
- T(n, 4) = n 2^(n-9) + 2 T(n - 1, 4) for n >= 8, from the base T(7, 4) = 6, so T(n, 4) = (n^2 + n - 8) 2^(n-10) for
  n >= 7 (and T(n, 2) = 2^(n-5) for n >= 5);
- the top diagonals N(n, 2) = N(n - 1, 2), N(n, 3) = N(n - 1, 3) + 1, N(n, 4) = N(n - 1, 4) + n - 4. QED.

## 3. Exhaustive checks (FINITE-EXACT)

**The program.** [`uh_noncritical_run_20261104.c`](../../04-computation/experiments/uh_noncritical_run_20261104.c)
generates every class of order n from those of order n - 1.
- It extends by a vertex whose in-neighbourhood is one cyclic interval of the parent's cycle, and keeps the extension
  when the sum over w in O and u in I of hp(w, u) equals 1. Here hp counts Hamiltonian paths, capped at 2.
- Canonical forms are minima over cycle rotations.
- Generation is complete by Moon's theorem and Lemma 1, so it does not use Shelor's grammar. The counts are checked
  against Garey's F(2n - 6) at every order.
- Per class it checks:
  - (A), (C) and (D), including the cut-arc predicate;
  - the degree-1 counts and asymmetry;
  - one-arc-cut existence;
  - at the cross-check orders: an exact Hamiltonian-cycle count by an independent subset DP, and a brute-force cut search
    over all 2^n vertex subsets.

**Exploratory run** (orders 4-13), [explore.out](../../04-computation/experiments/uh_noncritical_run_20261104_explore.out):
- every check passes;
- hc = 1 is confirmed independently at orders <= 12;
- the cut-arc rows equal HYP-9458's T(n, k).

**(G) and (H)** ([`shelor_parents_20261104.c`](../../04-computation/experiments/shelor_parents_20261104.c), the same
generator):
- for every class U of order 4-14, the vertices p with e_p(U) = 1 (exact subset DP) are exactly {y, y-} plus y-- iff y-
  is a sole entry;
- each child U^p has run length r, r + 1 or r + 2 as stated;
- |adm| = 2 for F(2m - 8) classes and 3 for F(2m - 7) (output: 6,765 and 10,946 at order 14).

**Confirmatory run** (orders 14-18, after pre-registration `a8db6218f9`): section 5 and
[confirm.out](../../04-computation/experiments/uh_noncritical_run_20261104_confirm.out).

## 4. The run-length table and its recurrence (HYP-9460 P7; PROVED as (H))

N(n, r) for r = 2..n - 1. This equals HYP-9458's row read backwards (k = n + 1 - r):

| n | N(n, 2), ..., N(n, n - 1) |
|---:|---|
| 5 | 1, 1, 1 |
| 6 | 1, 2, 3, 2 |
| 7 | 1, 3, 6, 7, 4 |
| 8 | 1, 4, 10, 16, 16, 8 |
| 9 | 1, 5, 15, 30, 41, 36, 16 |
| 10 | 1, 6, 21, 50, 85, 102, 80, 32 |
| 11 | 1, 7, 28, 77, 155, 231, 248, 176, 64 |
| 12 | 1, 8, 36, 112, 259, 456, 608, 592, 384, 128 |
| 13 | 1, 9, 45, 156, 406, 819, 1289, 1560, 1392, 832, 256 |

**Recurrence (HYP-9460 P7; PROVED by (G)-(H); the generating function sums over n >= 4).** For n >= 6,

    N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1).

Equivalently, the generating function is

    sum N(n, r) t^r x^n = t^2 x^4 (1 + (t^2 - t) x) / (1 - (1 + 2t) x + t x^2),

which at t = 1 is Garey's x^4 / (1 - 3x + x^2).

In HYP-9458's variables, T(n, k) = T(n - 1, k - 1) + 2 T(n - 1, k) - T(n - 2, k - 1). This implies all of its
fitted forms (section 2, end of the proof of (H)).

**The mechanism.**
- Under the canonical Shelor deletion of (C), every class has 2 or 3 children:
  - one with the same run length (parent y);
  - one with run length one larger (parent y-);
  - when it has two sole entries, a third with run length two larger (parent y--).
- The "- N(n - 2, r - 1)" term is the count of parents with a second sole entry, expressed one order down.

## 5. Confirmatory results (pre-registered HYP-9460, orders 14-18)

The run `uh_noncritical_run_20261104 18 14 14` was started after the pre-registration commit `a8db6218f9` and took 25 min
45 s on 12 threads. Output:
[confirm.out](../../04-computation/experiments/uh_noncritical_run_20261104_confirm.out).

| n | classes (= F(2n - 6)) | violations of (A), (C), (D), (E), degree counts, one-arc cut | (1,1), (1,2) = (2,1), (2,2) | N(n, r) row = recurrence |
|---:|---:|---|---|---|
| 14 | 17,711 | 0 (and hc = 1 and the brute-force cut search confirmed independently) | 2,585, 4,180, 6,766 | yes (non-blind) |
| 15 | 46,368 | 0 | 6,766, 10,945, 17,712 | yes (non-blind) |
| 16 | 121,393 | 0 | 17,712, 28,656, 46,369 | yes (non-blind; HYP-9458's published row) |
| 17 | 317,811 | 0 | 46,369, 75,024, 121,394 | yes (**blind**) |
| 18 | 832,040 | 0 | 121,394, 196,417, 317,812 | yes (**blind**) |

The blind order-17 row was confirmed before (G) and (H) were proved. The order-18 row was confirmed while they were
being written. All pre-registered predictions P1-P8 are CONFIRMED.

## 6. Relation to the concurrent work

| item | where | relation |
|---|---|---|
| every uniquely Hamiltonian tournament is a one-arc join; in- and out-degree-1 vertices | THM-4867 (opus-S74), KNOWN: Shelor 1977 | deferred to. (B) here is a second, Shelor-free proof |
| extension lemma, independent generator through order 16 | THM-4867 | used as Lemma 1 and re-implemented independently to order 18 |
| failed route (c), non-critical vertices | THM-4867 section 3 | refined: non-critical vertices are locally round (A), and propagation along C proves the theorem |
| T(n, k) cut-arc table and fitted columns | HYP-9458 (opus-S74: conjectured, confirmed at orders 14-16, then PROVED by its audit given Shelor's grammar, 12:51) | (D) gives T(n, k) = N(n, n + 1 - k). (G)-(H) are a second, Shelor-free proof of the recurrence and of every fitted form. Priority for the conjecture, the data and the first proof is opus-S74's |
| cut arcs under one Shelor insertion: nested largest one-arc cuts A(e), each containing the newest vertex | opus-S74 audit (`uh_one_arc_cut_20261104_audit_S74/AUDIT.md`, section A), given Shelor's grammar | (D) is the global form: the cut arcs are exactly the arcs off the non-critical run. (H) is the run-side form of the same insertion step |
| sole exits, digits, Phi bijection, Fibonacci coincidence | THM-4862 (klein) | (F) counts the classes by sole exits; rooted (UH, sole exit) pairs number F(2n - 5). (G) is the Shelor-side grammar behind Garey's count. Its admissible parents are klein's entry attachments with digit 1, via Lemma 3 |
| chord-language rules (end insertion, suffix of the last column) | klein two_branch packet section 5 (FINITE-EXACT, entry-attachment coordinates) | different coordinates. They are necessary conditions on forward chords; (A) and (C) describe the non-critical vertices instead |
| out-degree-1 lemma for one and two cycles | HYP-9459 (1) (klein) | the one-cycle half is (B), KNOWN by Shelor. The two-cycle half, "strong, minimum out-degree >= 2 and n >= 6 implies hc >= 3", stays OPEN; (B)'s round cycle is exactly RT_5's pentagram at n = 5 |

## 7. Inspiration sources (attached by the owner), checked for genuine connections

**Grogan, *An Arithmetic Family of Aperiodic Tilings in Odd Prime Dimensions* (6 Oct 2026).**
- Its tilings have unique hierarchical decompositions ("recognition of complete parents") and finite symmetry groups.
- The structural parallel here is (C) with (E): every uniquely Hamiltonian tournament has one canonical parent, the
  Shelor deletion of y, and no symmetry. So the Shelor tree is a recognizable hierarchy.
- Typed INTERPRETATION. No statement transfers.

**Abbadini-Almeida-Arrieta, *A topos for étale-finite Heyting algebras* (arXiv:2606.03861v2).** No genuine connection was
found. The closest object is the chain of in-closed vertex sets of a tournament (the initial segments of its
condensation). That chain is trivial exactly for strong tournaments, and nothing in the paper's Esakia-dual setting acts
on it.

**The emailed "three polygons tile the plane iff RH" chain.** The chain runs Matiyasevich's criterion, then a 275-state
Lean-checked machine, then Robinson-style Wang tiles, then Demaine-Langerman polygons.
- It is not treated as established, and no RH consequence is claimed or used.
- The only point of contact with this note is the repository's identification of tournaments with binary tilings of a
  staircase (THM-1430). In cycle coordinates, (A) says each non-critical row is a prefix word 1^k 0^(n-1-k): a local,
  tile-like rule. Uniqueness of the Hamiltonian cycle is global, and (A) to (C) say how much of it is local.
- Typed INTERPRETATION.

## 8. Typing

| item | type |
|---|---|
| (A) local roundness and propagation | PROVED (elementary) |
| (B) sole exit and sole entry exist; strong with min out-degree >= 2 implies hc >= 2 | KNOWN (Shelor 1977; THM-4867); a Shelor-free proof, new to the repository, PROVED here |
| (C) one run; the unique Shelor deletion | PROVED; the existence of Shelor's insertion is KNOWN (Shelor 1977), the uniqueness of the deleted vertex is PROVED here |
| (D) cut arcs + non-critical = n + 1; T(n, k) = N(n, n + 1 - k); T(n, n - 1) = 1 | PROVED. The top diagonal also follows from (H); Las Vergnas 1975 (CITED via Burzio-Demaria 1988) corroborates it |
| (E) asymmetry | PROVED |
| (F) refined counts F(2n-10)+1, F(2n-9)-1, F(2n-8)+1; F(2n-7); F(2n-5) | PROVED (using (G) for the total F(2n - 6)) |
| (G) admissible parents {y, y-, (y--)}; the unique Shelor deletion as a bijection; Garey's F(2n - 6) | statements KNOWN: Shelor 1977's grammar and its bijectivity, and Garey 1972 / Moon 1982. A Shelor-free proof, new to the repository, is PROVED here |
| (H) run transitions r, r + 1, r + 2; run-length recurrence; HYP-9458's formulas | PROVED (second, Shelor-free proof). HYP-9458 was conjectured by opus-S74 and first proved by its audit, given Shelor's grammar |
| (A), (C), (D), (E), degree-1 and (F) counts, one-arc cuts at orders 4-18; (G) and (H) on parents of orders 4-14 | FINITE-EXACT (orders 14-18 pre-registered). The audits add (G) and (H) at parent orders 4-15, and the (F) bijection, the rooted count, the A_n identification and the strong 'equivalently' form at orders 4-16 (strong form 5-9) |
| run-length recurrence | PROVED. First by opus-S74's audit via HYP-9458 and (C), given Shelor's grammar; independently by (H). Pre-registered (HYP-9460 P7) and FINITE-EXACT through the confirmatory run |
| tilings and topos readings | INTERPRETATION / no connection |

## 9. Audit

Two independent adversarial audits were run. Each auditor wrote its own code and used none of this session's programs or
outputs as evidence. Both verdicts are **SOUND WITH CORRECTIONS**, with no blocking finding. All corrections are applied,
and MISTAKE-707 logs them.

**Audit 1** ([report](../../04-computation/experiments/uh_noncritical_run_20261104_audit_S10/AUDIT.md)) covers (A)-(F) and
re-read (G)-(H).
- *Proofs.* Every step of Lemmas 1-3 and (A)-(H) holds. This includes the even-order second cycle of (B), the component
  argument of (D), and the base cases of (F): sole-entry preservation genuinely fails at n = 4, which is why the step
  starts at n = 5.
- *Computation.*
  - All labelled tournaments through order 7 (2^21 at order 7), with exact cycle counts and no generation lemma.
  - Plain-extension classes at orders 8 and 9 (55 and 144).
  - Shelor-style insertion at orders 10-16, equal to Garey's F(2n - 6).
  - Zero failures at orders 4-16 for (A)-(F), Lemma 1 and the (D) component structure. The (F) bijection hits every
    rooted pair once, and the r = 2 class is A_n.
  - (G) and (H) on parents of orders 4-15, and P7 at orders 6-16.
  - The quoted Las Vergnas statement over all strong tournaments of orders 4-9.
  - Negative controls fire.
- *Findings.*
  - One correction: the theorem's title stated the "equivalently" form of (B) without "strong", which is false as
    written (a source beating RT_5).
  - Six minor items: the FINITE-EXACT scope, the HYP/THM label map, "new" against the novelty disclaimer, an internal
    pointer, HYP-9460's P4 wording and stale note field, and the unstated base u(4) = 1.
  - Five remarks.

**Audit 2** ([report](../../04-computation/experiments/uh_noncritical_run_20261104_audit2_S10/AUDIT.md)) covers (G), (H)
and everything built on them.
- *Proofs.* No logical gap.
  - Each Hamiltonian path of T ending at p is counted exactly once, including p adjacent to y or p_0.
  - |U| >= 4 is exactly what keeps y_U--- outside {y_U, y_U-, y_U--}.
  - The elimination is valid for n >= 6, and the generating function is confirmed symbolically.
- *Computation.*
  - Labelled tournaments through order 7, and all tournaments on a fixed cycle through order 9 (n F(2n - 6)).
  - Two class generators agree as sets through order 16.
  - On every class of orders 4-15: (G), (H), the decomposition e_p(T) = [e_p(U) - hp] + W with W computed
    independently, and the bijection. Every class of orders 5-16 arises from exactly one admissible pair.
  - The rows N, N_2 and T match the recurrence and every published row. The closed forms hold for n <= 60. S74's audit
    closed form agrees with the recurrence for n <= 40.
- *Findings.* Three minor items: the scope row in section 8, missing ranges on two one-step closed-form recursions (false
  at n = 6 and 7 as written), and missing provisos in the (H) statement. Four remarks.

**Found by the session itself during the audits.** The statement of (G) is Shelor's grammar. It is recorded in
`operation_closed_synthesis_20261029/THEORY.md` section II, so (G) is retyped as a Shelor-free proof of a KNOWN statement
(section 1 (G), typing table).

## 10. Reproduction

```
gcc -O2 -march=native -fopenmp -o uhnc 04-computation/experiments/uh_noncritical_run_20261104.c
./uhnc 13 12 12    # exploratory, orders 4-13 (under 1 s)
./uhnc 18 14 14    # confirmatory, orders 4-18 (12 threads)
gcc -O2 -march=native -fopenmp -o sp 04-computation/experiments/shelor_parents_20261104.c
./sp 14            # (G) and (H) on every parent of order 4-14 (about 1 min)
```
