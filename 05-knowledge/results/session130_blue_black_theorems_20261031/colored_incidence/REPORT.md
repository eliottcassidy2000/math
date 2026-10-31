# Session 130: blue symmetry, black pairs, and the endpoint-deletion transport

**Scope.** General symbolic theorems for the fixed-Hamiltonian-path cube; an independently rebuilt, **FINITE-EXACT** order-5 classification; exact hostile probes at orders 4–7. No assertion of novelty over the full repository. Session129's `tournaments/REPORT.md` is the inherited source for the incidence span and node gauge; it attributes that mechanism to **THM-796, `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md`**. The present packet does not import that source or its census. The targeted theorem file was absent from the source snapshot available to this worker, so it is not claimed to have been re-audited here.

## Results first

1. **A genuine map preserving complementation:** endpoint deletion sends literal complement lines at order n to literal complement lines at order n−1. Every lower line has `2^(n−2)` lifts.
2. **Exact blue-lift theorem, all n≥4:** a lower line has precisely **two** blue lifts if its marked central `(n−2)`-vertex restriction is reflection-fixed, and no blue lifts otherwise.
3. **Order 5 is a uniform-color special case:** each of the four order-4 literal lines has exactly **2 blue and 6 black lifts**. Therefore upper color is independent of *every* observable of the deleted lower line under uniform sampling. This includes the selected arithmetic endpoint ratios `55/7, 1, 55/17, 17/55`: each appears twice among blue lines and six times among black lines. The map, not the numerical coincidence, explains these multiplicities.
4. **Exact color decomposition of a deletion fingerprint:** the raw vector at a merged node is the sum of its blue-parent vector and the unoriented incidence vector of a canonical multigraph of black reflection pairs. Consequently its raw arithmetic code factors into blue coordinate labels and one product `q_v q_w` for each black parent-pair edge.
5. **A two-number repair of the original Hamiltonian invariant:** `(H(T), number of blue tilings)` distinguishes all ten order-5 nodes. It fails already at order6, even when automorphism count and total fiber size are also retained.
6. **The loss boundaries are explicit:** deletion does not preserve color, does not descend even to the reflection quotient of lines, and primitive-vector normalization loses incidence mass. The order-5 uniformity fails immediately at order 6: only 16 of 32 lower literal lines have blue lifts.

## 1. Carriers and operations

Fix the directed Hamiltonian path `n -> n−1 -> ... -> 1`. Its nonpath arcs are tiles

`I_n = {(x,y): 1 <= y <= n−2, y+2 <= x <= n}`,

listed by increasing y and then decreasing x. Let `X_n={0,1}^{I_n}`. A zero bit directs `x -> y` and a one bit directs `y -> x`.

- `C_n t = t XOR all_ones`: complement only free arcs.
- `R_n t(x,y) = t(n−y+1,n−x+1)`: reverse path positions and take the full tournament converse.
- `p_n:X_n -> M_n`: forget the marked Hamiltonian path and identify isomorphic tournaments and their converses.
- `L_n=X_n/<C_n>`: literal complement lines. Each is an actual two-mask orbit for n≥3. Projecting endpoints to `M_n` permits loops and parallel lines.
- `d_n:X_n -> X_(n−1)`: delete vertex 1 and shift all remaining vertex names down by one. In this bit order, `d_n(t)=t >> (n−2)`.
- `e_n=R_(n−1) d_n R_n`: deletion at the other path endpoint, in the corresponding gauge.

A tiling and its literal line are **blue** if `R_n t=t`, and black otherwise. Since C and R commute, both endpoints of a line have the same color. Since `p_n R_n=p_n`, reflection stays inside every merged-node fiber.

These are three distinct quotients: paths-to-tilings modulo automorphisms, complement lines, and tournament nodes modulo isomorphism and converse. They must not be collapsed into one carrier.

## 2. Klein-four pairing and the Hamiltonian interpretation

### Proposition 1 — black lines have canonical parallel twins (PROVED, n≥3)

The commuting involutions C and R generate a Klein-four action on `X_n`. C has no fixed tilings. Nor can `R_n t=C_n t`: the tile `(n,1)` is fixed by reflection, and its bit would have to equal its own complement.

Consequently:

- a blue tiling has a two-element orbit `{t,Ct}`;
- a black tiling has a four-element orbit `{t,Ct,Rt,CRt}`;
- R fixes a literal line iff that line is blue;
- every black line has a distinct R-partner with exactly the same two merged-node endpoints.

Thus every black edge multiplicity, including a loop multiplicity, is even. This is an actual involution on the edges, not an inference from a parity table. The canonical R-quotient retains one edge per black pair and every blue edge once. It loses which reflected path presentation supplied the edge.

### Hamiltonian-path meaning

For a fixed ordinary tournament T, the automorphism group acts freely on its directed Hamiltonian paths. A path presentation determines exactly one mask; two paths determine the same mask iff they differ by an automorphism. Therefore

`|F_n(u)| = s(T) H(T)/|Aut(T)|`,

where `s(T)=1` for a self-converse class and `s(T)=2` for a non-self-converse pair. A mask is blue exactly when reversal of its path positions is an anti-automorphism of T. Thus blue degree counts symmetric Hamiltonian-path presentations, not all self-converse presentations indiscriminately.

Using the classical oddness of H(T), together with oddness of `|Aut(T)|`, proves that a self-converse node has an odd number of blue presentations. The latter group fact follows because a tournament has no nonidentity involutive automorphism: a transposed vertex pair would have its connecting arc reversed; Cauchy's theorem then rules out an even automorphism-group order. A non-self-converse pair has no blue presentation at all. Hence **a node has some blue incidence iff its tournament is self-converse**, while **pure-blue** says that every Hamiltonian-path presentation is symmetric and is strictly stronger. This corollary uses the classical Hamiltonian parity theorem; the lifting and decomposition theorems below do not depend on it. The checker verifies all these identities independently at orders 3–5.

## 3. Blue extension is a marked two-endpoint condition

### Theorem 2 — exact colored line lifting (PROVED, n≥4)

Deletion commutes with complement:

`d_n C_n = C_(n−1) d_n`.

It therefore induces a well-defined map `D_n:L_n -> L_(n−1)`. Every lower literal line has `2^(n−2)` upper lifts. A lower line has exactly two blue lifts iff the restriction of either endpoint mask to the central `(n−2)` path vertices is reflection-fixed; otherwise it has none.

**Proof.** Deletion forgets precisely the n−2 free bits incident with vertex 1, so each lower mask has `2^(n−2)` lifts. C-equivariance pairs the fibers above the two endpoints of the lower line, giving the same number of upper lines, not twice that number.

Now fix the lower mask, i.e. every arc on vertices `n,...,2`. To make its extension R-fixed, examine the R-orbits of free arcs. Every newly added free arc `(x,1)` with x<n is paired with the retained arc `(n,n−x+1)`; its bit is forced by the already specified retained bit. The only new reflection-fixed arc is `(n,1)`, whose bit is arbitrary. The arcs between central vertices `n−1,...,2` reflect among themselves; compatibility there is exactly reflection-fixity of the central marked restriction. The two forced path arcs at the ends match automatically. Thus there are exactly two blue extensions of an eligible lower mask and none of an ineligible one.

Eligibility is unchanged by C. Pairing blue extensions above the two lower endpoints gives precisely two upper blue lines. ∎

**Why n≥4 is stated.** At order 2 the free cube is a singleton and C is the identity; the order-3-to-order-2 line count does not have the same two-endpoint convention. No false base case is hidden in the formula.

The theorem also gives an explicit recursive parametrization of blue presentations: choose a blue central order-(n−2) presentation and n−2 independent bits, consisting of n−3 reflection-paired outer choices and the corner bit. If `o_n` is the number of R-orbits of free tiles, then

`o_n = o_(n−2)+(n−2)`, and `|Fix(R_n)|=2^o_n`.

This is a two-step recursion on *marked symmetric presentations*, not a node-to-node recursion. Explicitly, `o_n=floor((n−1)^2/4)`, so the fraction of blue tilings is `2^−floor((n−2)^2/4)`.

### Corollary 3 — exact order-5 independence (PROVED)

At n=5 the central restriction has order 3. Its single free tile is reflection-fixed for either bit, so every lower mask and line is eligible. Consequently

`D_5*(blue line counting measure) = 2 * (all order-4 line counting measure)`,

`D_5*(black line counting measure) = 6 * (all order-4 line counting measure)`.

In particular, for every function f on order-4 literal lines,

`sum_(upper blue l) f(D_5 l) = 2 sum_(lower l) f(l)`,

`sum_(upper black l) f(D_5 l) = 6 sum_(lower l) f(l)`.

After normalization by 8 and 24 the two distributions agree exactly. The same assertion holds for endpoint masks: each lower mask has two blue and six black lifts. Thus no statistic of the deleted marked tournament—not even its full mask—can improve the 1/4 prior probability of an upper mask being blue under uniform sampling.

The theorem does **not** say blue lines delete to blue lines. Four upper blue lines delete to lower blue lines, and four to lower black lines. The exact four lower line representatives and all blue lifts are:

| Lower line masks | Lower color | Upper blue line representatives t (line `{t,63−t}`) |
|---|---|---|
| `{0,7}` | blue | 0, 1 |
| `{1,6}` | blue | 10, 11 |
| `{2,5}` | black | 16, 17 |
| `{3,4}` | black | 26, 27 |

### Sharpness at the next order

At n=6 the central restriction has order 4 and need not be reflection-fixed. Exactly 16 of the 32 lower lines are eligible; each has 2 blue and 14 black lifts. The other 16 have 0 blue and 16 black lifts. At n=7 the eligible count is 128 of 512, with 2 blue and 30 black lifts per eligible line. The checker exhausts all upper lines through order 7. More precisely, upper color is independent of the deleted mask or literal line **iff n is 4 or 5** (within n≥4). For every n≥6 the central reflection acts nontrivially on free tiles, so both eligible and ineligible lower masks exist. The conditional blue probability is `2^(3−n)` on eligible masks and zero on ineligible ones. Therefore order-5 independence is a sharply delimited low-order theorem, not an all-order claim.

## 4. A canonical graph behind each deletion vector

For a merged node u, put `F(u)=p_n^−1(u)` and define its raw endpoint-deletion profile by

`r_u(v)=#{t in F(u): p_(n−1)(d_n t)=v}`.

Define the **blue-parent vector**

`b_u(v)=#{t in F(u): R_n t=t and p_(n−1)(d_n t)=v}`.

Define a multigraph `K_u` on the lower merged-node set: every nonfixed R-orbit `{t,R_n t}` in `F(u)` contributes one unordered edge

`{p_(n−1)(d_n t), p_(n−1)(d_n R_n t)}`.

The second endpoint is the merged class after deleting the other path endpoint, because reflection vanishes in the lower merged quotient. Loops and repeated edges are retained. This graph is canonical: swapping the chosen representative of an R-orbit merely swaps its endpoints.

### Theorem 4 — colored deletion-incidence factorization (PROVED)

Let `deg(K_u)` be the unoriented incidence vector, with loops counted twice. Then

`r_u = b_u + deg(K_u)`.

**Proof.** Partition the tilings in `F(u)` into R-fixed singletons and nonfixed pairs. A singleton contributes `e_v` to both r and b. A pair contributes `e_v+e_w` to r and exactly the corresponding edge's incidence vector to `deg(K_u)`. Summing proves the identity. ∎

Consequences include:

- `|F(u)| = |b_u|_1 + 2|E(K_u)|`;
- for every lower-node subset S, `sum_(v in S)(r_u(v)−b_u(v))` has the parity of the number of K-edges crossing the cut `(S,S^c)`;
- if `q_v` are any nonzero elements of a commutative multiplicative group, the **raw** arithmetic code factors as

  `product_v q_v^r_u(v) = (product_v q_v^b_u(v)) * product_({v,w} in E(K_u)) q_v q_w`.

This is the actual map from black reflection pairs to products of two parent labels. It preserves raw deletion incidence and its arithmetic valuation vector; it loses the individual marked masks. The sidecar `K_u` preserves the pairing that the degree vector alone loses. Primitive normalization requires the gcd as an additional sidecar if raw counts are to be recovered.

## 5. Complete order-5 decomposition

Use Session129's canonical node order. Lower coordinates A,B,C have Hamiltonian-path counts 1,3,5 and labels 7,17,55. An entry such as `3CC` denotes three loops at C, each contributing twice to its incidence degree.

| Node u | Raw r_u | Blue b_u | Black parent-pair graph K_u | Raw gcd |
|---:|---|---|---|---:|
| 0 | (1,0,0) | (1,0,0) | empty | 1 |
| 1 | (1,1,0) | (0,0,0) | AB | 1 |
| 2 | (2,3,5) | (0,0,0) | 2AC + 3BC | 1 |
| 3 | (0,1,0) | (0,1,0) | empty | 1 |
| 4 | (2,0,7) | (1,0,0) | AC + 3CC | 1 |
| 5 | (1,4,4) | (0,1,2) | AB + 2BC | 1 |
| 6 | (1,0,10) | (0,0,3) | AC + 3CC | 1 |
| 7 | (0,2,3) | (0,1,0) | BC + CC | 1 |
| 8 | (0,5,8) | (0,1,2) | 4BC + CC | 1 |
| 9 | (0,0,3) | (0,0,3) | empty | 3 |

The blue vectors sum to `(2,4,10)=2(1,2,5)`, exactly twice the lower-node fiber vector, as required by Corollary 3. Blue-parent vectors alone identify only five pairs of nodes: `{0,4}`, `{1,2}`, `{3,7}`, `{5,8}`, `{6,9}`. Black degrees alone also lose information: nodes 4 and 6 have the same black degree, and 0,3,9 all have zero black degree. The complete raw profile and its primitive version distinguish all ten nodes, independently rechecked here.

The three pure-blue nodes have empty K. Their raw arithmetic codes are `7,17,55^3`, whereas Session129's primitive codes are `7,17,55`. The division by gcd is exactly where one loses the threefold path-presentation multiplicity at the regular five-vertex tournament. Thus the selected primitive labels are legitimate static fingerprints but cannot be substituted directly into raw incidence identities.

### The original Hamiltonian count gains exactly the needed color coordinate at order 5

**FINITE-EXACT theorem.** The pair `(H(T), b(u))`, where `b(u)=sum_v b_u(v)` is the **number of reflection-fixed tilings**, is an injective invariant of the ten order-5 merged nodes. The ten pairs in node order are

`(1,1), (3,0), (5,0), (3,1), (9,1), (9,3), (11,3), (15,1), (13,3), (15,3)`.

Thus the three H-collisions are resolved precisely by blue symmetry:

- H=3: a non-self-converse pair versus the pure-blue `1 -> C3 -> 1`;
- H=9: one versus three blue tilings;
- H=15: the mixed node with one blue tiling versus the regular pure-blue node with three.

This combines the repository's original Hamiltonian-path invariant with the new visible color structure. Do not substitute the number of symmetric Hamiltonian paths for b: in a self-converse ordinary class that number is `b*|Aut(T)|`. For instance the regular order-5 tournament has 15 paths and automorphism order5, but three blue tilings.

**Sharp hostile order6 witness.** The fixed-path masks 1 and 77 belong to distinct merged nodes11 and16, but both have

`H=17, |Aut|=1, fiber size=17, blue tilings=5`.

Their raw deletion profiles over the ten order-5 nodes are respectively

`(2,0,0,0,15,0,0,0,0,0)` and `(0,4,6,0,0,7,0,0,0,0)`.

So neither adding automorphism order nor total degree repairs the two-number classifier beyond order5, while the finer deletion vector immediately separates this pair. A second pair, order6 nodes27 and28 represented by masks50 and43, has the common values `H=29, |Aut|=1, fiber size29, blue tilings7`. These are independently recomputed in this packet; all eight collision groups of `(H,b)` at order6 are retained in the JSON output. The successful pair is a finite classification theorem, not a universal sufficiency claim.

### Selected arithmetic line ratios, with an actual preserved operation

Orient each literal line by its smaller mask. For the four order-4 lines, use the ratio of the selected endpoint labels (terminal divided by initial):

| Lower line | Merged endpoints | Label ratio |
|---|---|---|
| `{0,7}` | A -> C | `55/7` |
| `{1,6}` | C -> C | `1` |
| `{2,5}` | B -> C | `55/17` |
| `{3,4}` | C -> B | `17/55` |

Deletion preserves this chosen line orientation because `t<32` implies `d_5(t)<4`, and it commutes with C. Pulling the ratio observable back through D produces each ratio **twice on upper blue lines and six times on upper black lines**. This is a faithful statement about actual complementation and endpoint deletion; it does not claim these ratios are PRIMEGAME transitions, nor that 7,17,55 are forced by tournament structure. Multiplicative independence is needed for decoding node fingerprints, not for the transport identity itself.

## 6. Quotient failures and the necessary sidecars

1. **Color is not preserved:** upper blue line `{16,47}` deletes to lower black line `{2,5}`.
2. **Deletion does not descend to R-orbits of lines:** the upper black masks 2 and 8 are R-partners, and their C-lines are therefore one edge of the R-quotient. Their deleted masks are 0 and 1, lying in distinct lower line R-orbits `{0,7}` and `{1,6}`. This remains a failure even after quotienting the target by R. The correct descended object is the unordered pair of the two endpoint-deletion lower R-orbits, not one chosen parent edge.
3. **Merged-node information is not a recursive state:** the inherited order-4 masks 2 and 4 are one merged node, while the upper pure-blue node with sole mask16 has respectively one and zero lifts over those masks. The extra path-presentation coordinate cannot be suppressed.
4. **The new K sidecar is not the complement-edge graph:** K pairs a tiling with its reflection *inside one node* and records its two deleted parents. The merged metagraph pairs a tiling with its free-arc complement, often crossing nodes. Confusing these two involutions would destroy the theorem's type.

The universal correct carrier is the marked cube with its commuting C action and its two endpoint deletion maps. The two-endpoint object descends farther than either endpoint map separately. At n5 the fortunate three-coordinate injectivity lets one encode the resulting node labels arithmetically; it does not remove the missing path coordinate for dynamics.

## 7. Exact graph totals (supporting the root cycle audit)

All rows use the ten merged vertices, including isolated ones in a monochromatic graph. Loops count as one edge for cycle rank.

| Graph | Edges | Components | Cycle rank |
|---|---:|---:|---:|
| All literal lines | 32 | 1 | 23 |
| Blue | 8 | 3 | 1 |
| Black literal lines | 24 | 4 | 18 |
| One edge per black R-pair | 12 | 4 | 6 |
| Full R-quotient | 20 | 1 | 11 |

The blue cycle is the square `5–6–9–8–5`. Black loops form one R-pair at node6 and one at node8. All black bundles have multiplicity two at order5; the universal theorem only forces an even multiplicity, not necessarily two at higher orders.

## 8. Reproduction and evidence

Run from any working directory, standard-library Python only:

```sh
python3 check_colored_incidence.py > normal.json
python3 -O check_colored_incidence.py > optimized.json
cmp normal.json optimized.json
```

The checker reads no external data and imports no inherited producer. It exhausts literal line lifts through order7; independently canonizes all fixed-path masks through order6 and all labelled tournaments through order5; computes H and automorphism counts, proves finite `(H,b)` completeness at order5 and emits order6 counterexamples; verifies the color-parent decomposition and arithmetic factors node by node; checks every lower ratio's colored multiplicity; checks the explicit quotient obstruction; and computes all graph component/cycle ranks. Failure gates use explicit exceptions and survive `-O`.

**Evidence boundary:** the finite checker supports the low-order classifications and hostile witnesses. The universal statements are justified by the symbolic involution and bit-extension proofs above, not by extrapolating that finite census. No Git or external state was changed by this worker.
