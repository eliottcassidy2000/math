# Local reversibility, rigid colors, and hidden sheets

**Status: PROVED elementary constructions and transfer identities; FINITE-EXACT controls.** This packet connects local-to-global inversion, graph coloring, and S165's guard language by explicit maps. Its graph covers are formal observation-state extensions, not claimed tournament realizations. It proves neither a Jacobian conjecture nor any Hamilton-path spectrum conjecture. Literature novelty is not asserted.

## Inheritance and the precise question

The closest inherited mechanisms are [S165's guard-state observation map and forbidden cone](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/guard_geometry/NOTE.md) and [its biclosed guard chamber banks](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/guard_adjunctions/NOTE.md). The hostile is a locally reversible map with several globally different lifts. The corrected near miss is confusing the absence of a symmetry with the absence of a second state. The useful sidecar is transport around cycles, or an explicit seed-language membership certificate.

There are three logically different failures:

1. **Erasure:** a move identifies inputs already at one step. The guard move D does this.
2. **Hidden sheets:** local inverse choices exist and are unique after choosing a sheet, but travel around a cycle changes the sheet. Graph monodromy measures this.
3. **Restricted image:** the coordinate map has a genuine global inverse, but the recovered point does not belong to the allowed source language. S165's three-count observation has this behavior.

Neither a nonsingular Jacobian nor a rigid coloring addresses all three.

## 1. A constant-Jacobian example, with the polynomial hypothesis kept visible

The following is an elementary **real-analytic, non-polynomial** map:

    F(x,y) = (exp(x) cos(y exp(-2x)),
              exp(x) sin(y exp(-2x))).

It is defined on all of R² and has determinant det DF=1. Indeed, with
r=exp(x) and theta=y exp(-2x),

    dF_1 wedge dF_2 = r dr wedge dtheta = dx wedge dy.

Nevertheless F(x,y)=F(x,y+2 pi exp(2x)). Its image is R² minus the origin, and
F(0,2 pi k)=(1,0) for every integer k, so it is not proper. This is a transparent local-to-global failure outside the polynomial category; it is **not** a counterexample to the planar polynomial Jacobian conjecture.

A proper local diffeomorphism R²→R² would be globally invertible: its image is nonempty, open, and closed; properness makes it a covering; and a connected covering of the simply connected plane has one sheet. The graph counterpart below shows exactly why simple connectedness cannot simply be replaced by rigidity.

## 2. Cycle transport is the complete missing datum for graph covers

Let B be a finite connected simple graph, and let p:E→B be a finite covering. Here covering means that for each lifted vertex x, p gives a bijection between the edges incident to x and those incident to p(x). Both graphs retain individual edges. Put d=|p^(-1)(b)|; path lifting makes this cardinality independent of b.

Choose a root b_0 and a spanning tree T. The unique tree paths identify every fiber with a common set F={1,...,d}. Each edge outside T then supplies a permutation of F by traversing its fundamental cycle. There are

    beta = |Edges(B)| - |Vertices(B)| + 1

such permutations. Let G≤Sym(F) be the subgroup they generate.

**Transport criterion.**

- E is connected if and only if G acts transitively on F.
- A chosen root sheet i extends to a path-independent inverse branch on B if and only if every fundamental-cycle permutation fixes i.
- All cycle transports are identity if and only if E is a disjoint union of d copies of B.
- The automorphisms of E that preserve p form the centralizer C_Sym(F)(G).

**Proof.** A path lifts uniquely once its starting sheet is fixed. Its ending sheet is found by composing the edge permutations, with inverse permutations on reversed edges. The sheets reachable over b_0 are exactly its G-orbit, proving the first claim. Two paths with the same endpoints differ by a closed walk; the spanning-tree fundamental cycles generate all such transport, proving the second and third claims. Finally an automorphism preserving p has a permutation sigma on the root fiber and is forced everywhere else by path lifting. It is well-defined precisely when sigma commutes with every closed-walk transport, proving the last claim. ∎

These are standard covering/voltage mechanisms, spelled out to make the test usable. Primary source credit for permutation-voltage constructions and the connectivity criterion: Dalfó, Fiol, Pavlíková and Širáň, [*Spectra and eigenspaces of arbitrary lifts of graphs*, §2](https://doi.org/10.1007/s10801-021-01049-3), which also credits Gross and Tucker's 1977 construction.

**Finite certification rule.** After verifying that p really is locally bijective, it suffices to check the beta fundamental-cycle transports. Checking vertex degrees, counting neighbors, or testing selected local patches without identifying the edge lifts does not establish the premise.

For a connected cover, a path-independent inverse through any root sheet implies d=1. Without connectedness it can select a one-sheet component while other components remain.

## 3. Sharp theorem: rigid colors can conceal arbitrarily many sheets

Color each vertex x of E by the full identity of p(x). This is a proper coloring, since B has no loops. A color-preserving automorphism is exactly an automorphism preserving p. A coloring is **distinguishing** when the only color-preserving automorphism is identity. This is not the same as saying that the uncolored graph has no automorphisms.

**Rigid ambiguity theorem.** A finite connected simple graph B admits a connected nontrivial finite covering whose base-identity coloring is distinguishing if and only if beta≥2. When beta≥2, such a covering exists with every number d≥3 of sheets. Three is the smallest possible nontrivial degree.

**Proof.** Choose two edges outside a spanning tree. Give their oriented transports

    a=(1 2),       b=(1 2 ... d),

give all other selected edge transports identity, and use inverse transports in the reverse directions. Define the lift on vertices (v,i) by replacing an oriented edge v→w with

    (v,i) -- (w, tau(i)).

This is locally bijective. The permutations a and b generate S_d: conjugates of a by powers of b include adjacent transpositions. Thus the lift is connected. For d≥3, the centralizer of S_d inside S_d is its trivial center, so the coloring is distinguishing.

Conversely, beta=0 gives a tree and only one-sheet connected coverings. If beta=1, monodromy is generated by a single permutation; transitivity forces it to be a d-cycle, whose centralizer contains all d powers of that cycle. Consequently a connected d>1 cover is not color-rigid. Finally every connected two-sheet cover has monodromy S_2, whose centralizer has order two. ∎

The theorem separates two questions which a local-to-global argument must not merge:

    Is there a symmetry between two sheets?   Centralizer question.
    Are there two sheets?                    Fiber/transitivity question.

The first answer can be no while the second answer is yes, with arbitrarily large fibers.

### A 12-vertex, four-color witness

Take the diamond graph on {0,1,2,3} with edges

    01, 02, 03, 12, 13.

Use the star 01,02,03 as tree. For d=3 put (1 2) on 12 and (1 2 3) on 13. The lift has 12 vertices and 15 edges. It is connected and 3-to-1, its projection is locally bijective, and its four-color base-identity coloring is both proper and distinguishing. The base has the smallest possible number of vertices for a simple connected graph of cycle rank at least two.

The checker exhausts all 36 pairs of S_3 transports on the two non-tree edges. Exactly 26 lifts are connected: 18 have trivial color-preserving automorphism group and eight have a deck group of order three. An independent test of all 6^4 fiberwise permutations confirms rigidity of the displayed witness.

### Proper coloring as the smallest transport model

For two colors, the constraint along each edge is the swap 0↔1. On a cycle of length ell, transport is the ell-th power of that swap. Every proper path segment can be colored, but the cycle closes if and only if ell is even. This is a global consistency obstruction, not a failure of any individual edge rule. For more than two colors, inequality constraints are not generally permutation transports; one must not apply the covering criterion without additional structure.

## 4. A literal application to S165's guard observation bank

Use the denominators (4,6,9,5,7). S165 proves that the biclosed guard observations, joined by one-bit changes, form

    B = C_6 square C_4.

This graph has 24 vertices, 48 edges, and beta=25. Applying the theorem constructs, for every d≥3, a connected d-sheet graph E_d with:

- 24d vertices and 48d edges;
- the same four labeled guard-flip possibilities as its projected observation;
- exactly d states over each full five-bit observation;
- no nonidentity automorphism preserving those five-bit observations.

For d=3 the checker constructs all 72 vertices and 144 edges and verifies every local neighbor bijection, connectedness, and all possible root-fiber permutations. In its ascending-bitmask BFS tree, the first two non-tree edges are {6,7} and {8,9}; they receive the transposition and three-cycle respectively. The remaining non-tree transports are identity.

**Scope boundary.** E_d is a formal reversible extension of the guard **chamber graph**. No arithmetic FRACTRAN transition law, tournament realization, or blue-metagraph realization of its added sheet coordinate is claimed. S165's actual blue lines connect complementary observations, whereas chamber edges change one guard bit. Those are different graphs.

Even the original unmarked-tournament quotient of the C_6 bank is not automatically a cover. In S165's concrete five-vertex tilings, adjacent observations 100 and 110 have the same canonical tournament code 12. Retaining all six chamber edges, their quotient vertex has degree four (including the loop twice), while each original vertex has degree two. Thus the local-bijection premise already fails. One must first keep enough marked data to establish a covering before interpreting any monodromy.

## 5. The three-count guard observer has integer dynamics

Retain S165's guarded tournament state v=(s,u,d) and temporal moves

    C(s,u,d)=(u,s+u+d,s),
    D(s,u,d)=(s+d,2u,0).

Here s counts Hamilton paths starting at the guard q, u those not starting there, and d counts paths in X−q starting at q's unique in-neighbor. Let

    h=(h_0,h_1,h_2)=(H(X),H(CX),H(C²X)).

The inherited observation matrix and inverse are

        [1 1 0]                     [ 1 -2  1]
    O = [1 2 1],      O^(-1) = 1/2 [ 1  2 -1].
        [3 3 2]                     [-3  0  1]

Direct conjugation gives the following **integer** observed dynamics:

    C_h(h_0,h_1,h_2)=(h_1,h_2,h_0+h_1+h_2),
    D_h(h_0,h_1,h_2)=(h_1,h_0+3h_1−h_2,3h_1).          (1)

Thus the C observer is exactly the three-coordinate tribonacci shift. The integral image of O is the index-two lattice

    O Z³ = {h in Z³ : h_0 = h_2 mod 2}.                 (2)

Both observed moves preserve that lattice. These statements follow either by multiplying the displayed matrices or substituting the inverse formulas. They concern the full integer state space, not just the reachable language.

This is a case where the inverse is globally known and the only integral-coordinate obstruction is explicitly parity. A putative obstruction to seed-language realizability must therefore live beyond the Jacobian and lattice test.

### An actual pair that D erases

Starting at the directed-triangle seed v_0=(1,2,0), temporal words CDC and DDC give

    v_A=(6,9,3),     h_A=(15,27,51),
    v_B=(8,9,1),     h_B=(17,27,53).

Both are literal strong six-vertex guarded tournaments. Applying D gives the same state

    Dv_A=Dv_B=(9,18,0),
    D_h h_A=D_h h_B=(27,45,81).

The difference h_B−h_A=(2,0,2) lies in ker D_h. All subsequent C/D count observations agree. The underlying seven-vertex tournaments are not asserted isomorphic. This is genuine one-step information erasure, so a covering or nonzero-Jacobian argument is unavailable for D.

## 6. A forbidden interval in observable coordinates

Let O_seed denote only the C/D states reached from v_0. The symbol does **not** mean all physically realizable guarded tournaments. S165 proves that an O_seed state with d=0 cannot have

    1/4 < s/u < 3/8.

Since d=0 is equivalent to h_2=3h_0 and

    s/u=(2h_0−h_1)/(h_1−h_0),

the exact observable exclusion is

    h_2=3h_0  and  19/11 < h_1/h_0 < 9/5.              (3)

**Proof of the interval conversion.** On d=0, h_1−h_0=u>0. Multiplying the two strict ratio inequalities by this positive denominator yields, respectively, 5h_1<9h_0 and 11h_1>19h_0. ∎

For the literal strong family Y_n of S165, n≥6, put t=2^(n−4). Its state and observed triple are

    v(Y_n)=(t+1,4t−2,0),
    h(Y_n)=(5t−1,9t−3,15t−3).                          (4)

Because t≥4,

    11(9t−3)−19(5t−1)=4t−14>0,
    9(5t−1)−5(9t−3)=6>0,

so every triple lies in (3). The first is (19,33,57). It satisfies the parity lattice condition, has a perfectly well-defined inverse, and consists of three actual tournament counts. Nevertheless its state is not in O_seed, even up to positive scaling.

There is no failure of invertibility here. What fails is the inference that an invertible coordinate map makes its restricted source language surjective. The distinction is stronger than an arbitrary fake integer example because the excluded points are themselves literal guarded tournaments.

### The obstruction survives every C tail

For all n≥6 and k≥0,

    C^k v(Y_n) is not a positive scalar multiple of any O_seed state.   (5)

**Proof.** It is true for k=0 by the inherited gap. For k>0, all intermediate C^j v(Y_n), 1≤j≤k, have strictly positive d. Any nonempty physical word with positive final d must end in C, since D outputs d=0; the seed also has d=0. Moreover C is invertible, with

    C^(-1)(s,u,d)=(d,s,u−s−d).

Thus a word reaching a positive multiple of C^k v(Y_n) must end in k consecutive C moves. Removing them would reach the same positive multiple of v(Y_n), contradicting the base case. ∎

The observer forms of these excluded triples are C_h^k h(Y_n). This gives an infinite family of honest, locally reconstructible tournament states outside the chosen seed language, not merely a gap confined to the d=0 boundary.

## 7. What transfers, and what does not

- **Analytic map → covering viewpoint.** The map F is locally invertible and becomes an infinite-sheet covering of its punctured image. The missing datum is the sheet/angle lift, and properness fails. This is an analogy with an explicit model, not a reduction of a polynomial problem.
- **Marked graph states → guard observations.** The projection E_d→B preserves every local guard-flip option. It forgets a sheet coordinate. Fundamental-cycle permutations recover the loss; color-rigidity alone does not.
- **Guard state → three counts.** O preserves the entire three-dimensional state over Q, with the exact integral parity lattice. No state coordinate is lost. The forgotten object is instead a construction history and its admissibility in the C3-seeded language.
- **Guard move D → next state.** D loses a rank-one direction. Retaining a suitable extra coordinate, such as the old d together with the output, restores its input over Q: old u=u'/2 and old s=s'−old d. This does not restore the underlying tournament itself.

A valid local-to-global claim must identify which of these problems it solves. For a genuine connected graph cover, trivial cycle transport is sufficient and exact. For the guard observer, the inverse already exists; the remaining task is a semigroup membership or source-construction theorem. Neither statement turns full matrix-algebra generation or Zariski density into positive-word reachability.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. Their JSON outputs are byte-identical. The controls include:

- all 6,859 observation triples in [-9,9]^3 for the exact lattice criterion;
- the transfer conjugations and physical D collision;
- all 65,535 temporal words of length at most 15, including 32,768 ending with d=0;
- independent endpoint Hamilton-path DP for Y_n, n=6,...,10;
- all 36 S_3 diamond voltage assignments and a separate 1,296-map colored-automorphism check;
- connected diamond covers of degrees two through seven, with the expected degree-two failure of rigidity;
- the complete 72-state lift of the 24-state guard bank;
- all two-color assignments on cycles of lengths three through ten.

Finite computations corroborate the proofs. The all-degree cover theorem, its sharp cycle-rank boundary, the observation identities, and the forbidden C tails are proved above independently of the finite search limits.
