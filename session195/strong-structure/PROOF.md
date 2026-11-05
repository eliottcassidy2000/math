# Strong-core module hulls and a sharp singleton-growth threshold

Session 195, 2026-11-05. **PROVED elementary structural and counting theorems; FINITE-EXACT controls are separate.** The unrestricted Hamiltonian-path spectrum with proposed exceptions `{7,21}` remains **OPEN**. No claim of literature novelty is made.

## 1. Inheritance and purpose

The closest inherited mechanism is the sharp intact-module threshold in [Session 194, module-barrier/PROOF.md](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/module-barrier/PROOF.md). It says that retaining a nontransitive input as a module cannot give literal `3h+2` or `3h+4`. The new connection is an exact characterization of when a fixed strong induced core can grow by single-vertex **strong** extensions, followed by a sharp numerical sufficient condition for that characterization. The proof below is self-contained.

The canonical hostile is the cyclic extension `S -> a -> b -> S`: it is strong, but neither single-vertex enlargement of S is strong. The corrected near miss is the claim that an arbitrary fixed strong core always has a strong singleton ancestry to its ambient strong tournament. The missing coordinate is the smallest module containing that core. A second hostile below shows that even `(order,H,F2)` does not determine the available one-vertex count responses.

The live concepts are induced strong cores, module hulls, mixed one-vertex additions, two-vertex directed ears, ordered two-path covers, and complete endpoint tables. The aim is a faithful recursive grammar, not an unproved inference from structural completeness to arithmetic coverage.

## 2. Definitions and two monotonicities

All tournaments are finite. A Hamiltonian path is a directed vertex word using every vertex exactly once. Write `H(U)` for its count and `H(empty)=1`. A singleton is strong; there is no strong tournament of order two. A nonempty subset M of a tournament T is a **module** if every vertex outside M either beats all of M or loses to all of M.

Let `F2(U)` count ordered pairs of nonempty directed paths whose disjoint vertex sets partition U. Thus

    F2(U) = sum_(empty != A proper-subset V(U)) H(U[A]) H(U-A).

In particular, `F2(K1)=0`. For order n,

    F2(U) >= 2^n-2,        F2(U) >= (n-1)H(U).            (1)

The first inequality follows because every induced tournament has a Hamiltonian path. The second follows by splitting each Hamiltonian word at each of its n-1 internal gaps; concatenation recovers the word and split.

**Lemma 1 (induced monotonicity).** For nonempty induced `S subseteq U`,

    H(S) <= H(U),           F2(S) <= F2(U).               (2)

**Proof.** A new vertex z can always be inserted into a nonempty directed word: put it immediately before the first vertex it beats, or append it if it beats none. At an internal insertion position the preceding vertex beats z by the definition of the first beaten vertex. Insert all missing named vertices in a fixed order. Deleting them recovers the original word, so this is an injection for H.

For an ordered two-path cover `(P,Q)` of S, insert every missing named vertex into P alone by the same deterministic rule, leaving Q unchanged. Both paths remain nonempty, and deleting the inserted labels recovers `(P,Q)`. This proves the F2 inequality. If S is a singleton its two-cover set is empty and the inequality is immediate. This argument works at all orders and requires no strongness. QED.

## 3. Exact module-hull classification

For a nonempty strong induced U in an arbitrary tournament T, an exterior vertex z is **mixed on U** when it beats at least one U vertex and loses to at least one U vertex.

**Lemma 2.** `T[U union {z}]` is strong if and only if z is mixed on U.

**Proof.** If z is uniform, it is a source or sink relative to U, so the enlargement is not strong. If z is mixed, every U vertex reaches a vertex beating z within U, and z reaches a vertex of U that it beats. Strongness of U supplies all remaining connections. QED.

**Theorem 3 (canonical hull).** Start with any nonempty strong induced core S of an arbitrary tournament T. Repeatedly add any exterior vertex mixed on the current set, stopping when none exists. Then:

1. Every intermediate induced tournament is strong.
2. The terminal set `M_T(S)` is independent of the choices.
3. It is exactly the smallest module of T containing S.
4. There is a chain from S to T adding one vertex at a time and keeping every induced tournament strong **if and only if** S is contained in no proper module of T. When this holds, **every** maximal greedy chain reaches T.

**Proof.** Assertion 1 is Lemma 2. At termination every exterior vertex is uniform on the current set, so the terminal set is a module.

Let M be any module of T containing S. Inductively suppose the current set U lies in M. If an eligible z were outside M, the module condition would make z uniform on M and therefore on U, a contradiction. Thus every added vertex lies in M. Consequently the terminal module lies in every module containing S, proving that it is the unique smallest such module and is independent of all choices.

If this module is T, the construction reaches T. Conversely, a strong singleton chain starting from S must stay inside every module containing S, by the same argument and Lemma 2. A proper such module therefore prevents the chain from reaching T. This proves 4. QED.

Checking only whether S itself is a module is insufficient. For the tournament with outgoing rows

    [22,20,24,19,32,15],

the triangle `S={0,2,3}` is not a module, but its module hull is `{0,1,2,3}`. It can add vertex 1 and then cannot grow further by strong singleton steps. The ambient tournament is strong and has H=37. The verifier checks these facts directly.

This is a criterion for a **specified induced core**, not a claim about arbitrary deletion sequences or an assertion that every strong tournament is prime. In a prime tournament, every strong core of order at least three consequently grows to the whole tournament by strong singleton additions.

## 4. A universal one-or-two-vertex ear grammar

**Lemma 4 (escaping a proper module inducing a strong tournament).** If T is strong and U is a nonempty proper strong induced core with no mixed exterior vertex, there are exterior vertices a,b with

    U -> a -> b -> U.                                  (3)

In particular, adjoining a,b yields a strong induced tournament.

**Proof.** Partition the exterior into A, the vertices beaten by all of U, and B, the vertices beating all of U. Strongness of T forces both classes nonempty. If every B vertex beat every A vertex, B would beat all of `U union A`, contradicting strongness. Hence some `a in A` beats some `b in B`. This gives (3). Every U vertex can reach every other through a,b as necessary, and all new and old vertices communicate. QED.

**Theorem 5 (faithful universal grammar).** Every strong tournament T can be grown from any specified nonempty strong induced core S by these operations, with all induced intermediate tournaments strong:

- add one mixed exterior vertex;
- if no such vertex exists, add a two-vertex ear (3).

Singleton additions terminate canonically at the module hull. Each pair step escapes that proper module, after which singleton closure can resume. Thus the successive saturated cores are a strictly increasing chain of modules of T, each inducing a strong tournament. The procedure ends at T after at most `|T|-|S|` operations. Scanning the remaining vertices for mixed adjacency, and scanning pairs only at a saturated stage, gives an elementary polynomial-time structural algorithm. Counting is a separate task.

**Proof.** Lemmas 2 and 4 show that one operation is available at every proper intermediate core. Each enlarges the core, so termination occurs. Theorem 3 identifies every saturated core as its canonical module hull. As the core sets are nested, their saturated hulls are nested too. QED.

Starting from a singleton is allowed: its first operation must be a pair and produces a directed triangle. Starting from a directed triangle avoids this boundary step. This theorem describes all orders without changing any old arc. It does not say the associated scalar counts cover any specified arithmetic set.

## 5. Exact pair count and the sharp threshold

**Lemma 6 (cyclic pair count).** For every nonempty tournament U, the cyclic extension `C(U)=U -> a -> b -> U` is strong and

    H(C(U)) = 3H(U)+F2(U).                             (4)

**Proof.** If U occurs in one maximal run, the three possible block orders are `U,a,b`, `a,b,U`, and `b,U,a`, each contributing H(U). If U occurs in two runs, the word must be `P,a,b,Q`, and `(P,Q)` is an arbitrary ordered nonempty two-path cover of U. No third U run is possible because only two exterior vertices exist and the reverse edge b->a is absent. These cases are disjoint and exhaustive. Strongness follows from the cyclic block incidences. QED.

**Theorem 7 (sharp anchored-growth gap).** Let S be any nonempty strong induced core in a strong T, and write `h=H(S)`. If

    H(T) < 3h+F2(S),                                   (5)

then the module hull of S is T. Equivalently, S grows to T by singleton additions through strong induced tournaments; every maximal greedy sequence of such additions reaches T.

The strict threshold is sharp **for each S individually**. At equality, `T=C(S)` is a strong example in which there is no first singleton step from S.

**Proof.** If a maximal singleton chain stops at proper U, Theorem 3 makes U a module inducing a strong tournament. Lemma 4 gives an induced `C(U)` in T. The two monotonicities and (4) imply

    H(T) >= H(C(U)) = 3H(U)+F2(U) >= 3h+F2(S),

contradicting (5). For sharpness, `C(S)` has the stated count by (4), while its only exterior vertices are respectively a sink and a source relative to S, so neither can be added alone while retaining strongness. QED.

Equality is not itself a no-growth criterion. The regular cyclic five-tournament (vertex i beats i+1 and i+2 modulo 5) has H=15, contains a C3 core S with `3H(S)+F2(S)=15`, and its module hull is the whole five-tournament. Thus both growth and obstruction occur at equality; only the strict sufficient condition and its sharpness are asserted.

The assumptions allow `S=T`, when the empty chain works, and allow singleton S. In the latter case the threshold is 3 and no proper strong T satisfies the strict inequality: every strong tournament of order greater than one contains a directed triangle and hence has at least three Hamiltonian paths. The equality example is C3. No forbidden order-two intermediate is being asserted.

**Corollary 8 (literal ternary induced ancestry).** If `h=H(S)>1`, S and T are strong, S is induced in T, and

    H(T) in {3h+2,3h+4},

then S has a complete strong singleton ancestry to T, and is contained in no proper module of T.

**Proof.** Here `|S|>=3`, so `F2(S)>=2^3-2=6`. Both target counts are strictly below `3h+F2(S)`. Apply Theorem 7. QED.

This is a **necessary structural consequence of an existing induced realization**, not an existence theorem for either target. It does not assume T is only one vertex larger than S. Conversely, having full module hull is insufficient to force either scalar count.

For the canonical grammar with r pair steps starting from S of order n, (1), (2), and (4) also give the useful depth bound

    H(T) >= H(S) product_(j=0,...,r-1)(n+2j+2).         (6)

Indeed, just before pair step j the core order is at least n+2j, and that step multiplies its count by at least its order plus two. Singleton steps cannot decrease H. This is only a lower bound, not an equality classification.

## 6. Faithful counting at mixed singleton steps

For a labelled U, let `s_x(U)` and `e_y(U)` count its Hamiltonian paths starting at x and ending at y. Define `B_yx(U)` as the number of ordered nonempty two-path covers `(P,Q)` with P ending at y and Q starting at x. Thus `sum_(x,y) B_yx=F2(U)` and `B_xx=0`.

Let a new vertex v beat exactly `X subseteq V(U)`, and put `Y=V(U)-X`. Then the exact response is

    H(U+v) = sum_(x in X) s_x(U)
              + sum_(y in Y) e_y(U)
              + sum_(y in Y,x in X) B_yx(U).            (7)

**Proof.** Partition Hamiltonian words by whether v is first, last, or internal. In the last case delete v to recover the unique cover `P,Q`; its boundary arcs are valid exactly when `y->v->x`, that is, y lies in Y and x in X. Reinsert v to invert the correspondence. QED.

For strong U, the permitted strong singleton steps are precisely the nonempty proper X. Formula (7) is a faithful finite interface for one step. It is **not claimed to be a closed update rule for B**. One explicit exact recursively maintainable state is the complete subset endpoint table:

    E(A,j) = number of Hamiltonian words of U[A] ending at j,
    E({j},j)=1,
    E(A,j)=sum_(i in A-{j}, i->j) E(A-{j},i).           (8)

After adding vertices, retain the old subsets, compute new subsets in increasing size, and use the known incidences. Every H and B above is recovered by subset sums, with starts obtained by reversing arcs. This state grows exponentially; no bounded-size closure or arithmetic coverage is asserted. The two-vertex ear can use (4) for its count and (8) for its full state.

**Finite hostile to scalar compression.** Both strong five-tournaments below have `H=9` and `F2=66`, using outgoing bitmask rows:

    A = [16,1,3,7,14],
    B = [20,1,2,7,14].

Nevertheless their exact sets of strong one-vertex child counts are different:

    A: {15,17,19,23,25,27,29,33,37,41},
    B: {15,17,19,23,29,31,33,37,43}.

Thus even `(order,H,F2)` cannot be substituted for the marked response state. The verifier checks all 30 nonempty proper X for each example using independent endpoint DP. It also checks (7) against direct extension counts throughout the stated finite universe. This hostile refutes scalar sufficiency, not the structural grammar.

## 7. Relation to the previous barriers and the spectrum gap

The exact hull classification explains why the intact-module barrier appears: a proper module is precisely what stops singleton strong growth from the selected core. The sharp H threshold removes that obstruction, but does not create either desired child count.

The [S194 weak-deletion theorem](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/weak-deletion/STRONGER_RATIO.md) concerns an actual one-vertex deletion that is weak. Corollary 8 instead supplies a strong-deletion ancestry for a specified retained strong core and small relative target, irrespective of the order difference. It therefore focuses the unresolved induced ternary mechanism on mixed singleton responses (7), rather than arbitrary weak padding or an intact core module.

The [S194 endpoint-floor theorem](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/endpoint-floor/PROOF.md) still applies to every step: a retained strong core with endpoint floor q forces every positive ambient endpoint count to be at least q. Nothing about a full module hull restores a unit endpoint or permits the known ternary fan to be rearmed.

The existing S194 obstruction to an induced chain `H13 subset H41 subset H127` is not bypassed: Corollary 8 applies only when the specified embedding already exists. Structural singleton interpolation cannot supply an absent embedding.

The remaining arithmetic task is to exhibit a sufficiently rich collection of marked states and constructive transitions, or to reconstruct different witnesses, so that every required prime and 7p atom is covered. The structural grammar, count formulas, sharp obstruction, and exact finite tests here do not establish that coverage. The full `{7,21}` spectrum remains **OPEN**.

## 8. Reproduction

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The standalone verifier uses explicit exception gates, not `assert`. Its receipt specifies the exact exhaustive universe, literal-permutation controls, boundary cases, and hostile checks. Finite checks are evidence for implementation and controls; the all-order conclusions follow from the proofs above.
