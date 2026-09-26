# Square-sum endpoints, the first Hamiltonian two-edge switch, and a signless lift obstruction

**Status: PROVED (elementary, with no priority claim):** the unicyclic Hamiltonian-path classification in section 2; the first square-sum four-cycle and first Hamiltonian-cycle two-edge switch both occur at `n=46` (section 3); and the exact signless-incidence criterion for designated-square-preserving affine lifts (section 4). **FINITE-EXACT:** an independent path census for every `1<=n<=25`; all 3898 labelled connected unicyclic graphs on 3..6 vertices; triangle/four-cycle birth through 46; two explicit Hamiltonian cycles in `Q_46`; and positive/hostile lift controls. **INHERITED:** the first nontrivial Hamiltonian path at 15, its endpoints 8 and 9, nearby path counts, connectivity, reflection algebra, and the corrected endpoint rule. **OPEN:** any implication for Collatz descent or the occurrence rate of long Collatz trajectories. No tournament is imposed on this undirected relation.

Script: [crossroads_family_20260926_squares.py](../../04-computation/experiments/crossroads_family_20260926_squares.py). Output: [crossroads_family_20260926_squares.out](crossroads_family_20260926_squares.out). All checks execute under ordinary Python and `python -O`.

## 1. Inheritance, the actual question, and the board

Let `Q_n` have vertices `1,...,n`, with an undirected edge `{x,y}` exactly when `x!=y` and `x+y` is a square. This is an intrinsic binary relation. Its edge label is the positive square root of `x+y`; there is no orientation gauge or tie-breaking operation.

The current sources are:

- [Square-sum Hamiltonicity, sections 1–4](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md): the current exact census and the endpoint correction. In particular, 15 is the first **nontrivial** size; the one-vertex path at 1 is the trivial exception.
- [Summand closure and square filter, sections 2–4](collatz_mod6_20260922_w6_summand_closure_square_filter.md): three chains through 12, the merges at 13 and 14, connectivity for all `n>=14`, the degree formula, and reflection composition. The founders are `{1,2,4}`, not the distinct-summand closure's missing set `{1,4,6}`.
- [MISTAKES.md, 2026-09-25 square-sum endpoint correction](../../01-canon/MISTAKES.md): graph degree 2 does **not** force both incident edges in a Hamiltonian path. The vertex may be an endpoint. The repaired endpoint-aware nonexistence argument for `Q_24` is [decoder note, section 4](decoder_prime_square_20260925.md).
- [THM-4497, coprime graph rank compression](../../01-canon/theorems/THM-4497-coprime-graph-rank-compression.md): weighted two-support valuation rows transport a kernel variable across an edge. Odd cycles kill an ungrounded component; even cycles require their exponent ratios. Section 4 below identifies an exact equal-exponent instance of this operator.
- [THM-2521, K13 drift / K14 potential bridge](../../01-canon/theorems/THM-2521-k13-drift-k14-potential-module-bridge.md), equation (4): the signless-incidence map sends a vertex potential to its endpoint sums. Its physical interpretation there is different.

The closest inherited mechanism is the reflection `x -> s^2-x`; the canonical hostile is the degree-2 endpoint 22 in `Q_23`; the corrected near miss is endpoint-blind forcing; and the useful neglected coordinate is whether we are retaining **one selected path** or **all square edges**. Root suggested testing affine lifts against the signless-incidence kernel; that becomes an exact bridge in section 4.

| Live concept | Operation and retained predicate | Lost coordinate / decisive test |
|---|---|---|
| Endpoints at 15 | Delete one edge from a unicyclic graph; retain a spanning path | Deletion must preserve connectivity and maximum degree 2 |
| Local rearrangements | Exchange two edges while retaining all path/cycle degrees | A four-cycle is necessary; cyclic ordering is still needed |
| Square arithmetic | Alternate edge labels around a closed walk | Positive, distinct vertices and a common height bound |
| Affine copies | Square-scale the vertex values and add offsets | Selected-path versus whole-graph edge obligations |
| Valuation graph | Identify the same homogeneous linear operator | Exponent weights, other prime rows, the factor 3, and chronological Collatz realizability |

### Exact nearby census

Paths are counted up to reversal. This independently reproduces the inherited answers through 25:

| n | Number of paths | Possible endpoint pairs |
|---|---:|---|
| 2..14 | 0 | none |
| 15 | 1 | `{8,9}` |
| 16 | 1 | `{8,16}` |
| 17 | 1 | `{16,17}` |
| 18..22 | 0 | none |
| 23 | 3 | `{18,2}`, `{18,9}`, `{18,22}` |
| 24 | 0 | none |
| 25 | 10 | 18 paired with `2,3,4,8,9,10,11,13,22,23` |

The unique path at 15 is

```text
8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.
```

The leaves 8 and 9 force its endpoints. Minimality is elementary from the inherited graph structure: `Q_n` is disconnected for `2<=n<=13`, while `Q_14` is a tree with three leaves 8,9,10. A spanning path cannot have three leaves. This is a finite startup mechanism; it does not itself assert a recurring scale law.

## 2. Why 15, 16, and 17 are rigid: a complete unicyclic rule

**Proposition (PROVED).** Let `G` be a finite connected simple unicyclic graph, with unique cycle `C`. It has a Hamiltonian path exactly when there is a cycle edge `{u,v}` with `deg(u)<=3`, `deg(v)<=3`, and `deg(w)<=2` for every other vertex `w`. Equivalently:

1. Every tree attached to `C` is a simple pendant path, attached by an endpoint.
2. There are at most two nonempty pendant paths, at distinct cycle vertices.
3. If there are two, their attachment vertices are adjacent on `C`.

The number of unoriented Hamiltonian paths is `|C|` when there is no tail, 2 when there is one tail, and 1 when there are two tails at adjacent attachment vertices. It is 0 in every other case.

**Wording repair (root audit, 2026-09-26).** The initial sentence only said the deleted edge's endpoints account for every degree-3 vertex, without explicitly capping those endpoints at degree 3. A cycle vertex with two pendant tails has degree 4 and is a hostile to that incomplete wording. The explicit degree caps above repair the criterion; the equivalent structural conditions, proof, and exact script already imposed the correct caps and are unchanged.

*Proof.* A connected unicyclic graph on `n` vertices has `n` edges. A Hamiltonian path therefore deletes exactly one edge. To keep the graph connected, this edge must be on the unique cycle. The remaining graph is a tree; it is a path if and only if all its degrees are at most 2. Deleting one cycle edge lowers exactly two adjacent cycle degrees by 1. This proves the criterion and the counts. The pendant-path formulation is the same degree condition written structurally. QED.

For `Q_15`, the unique cycle is

```text
1–3–6–10–15–1,
```

with tails `1–8` and `3–13–12–4–5–11–14–2–7–9`. The attachment vertices 1 and 3 are adjacent, so the unique Hamiltonian path deletes `{1,3}`. Thus the omitted square is 4; vertex 4 is an ordinary interior vertex. Adjoining 16 extends the tail ending at 9; adjoining 17 extends the tail ending at 8. Adjoining 18 branches at 7 and produces a third leaf. The same proposition explains all four outcomes without applying the invalid degree-2 rule.

**FINITE-EXACT independent control:** every labelled connected unicyclic simple graph on 3,4,5,6 vertices was checked against a separate exhaustive path enumerator: respectively 1,15,222,3660 graphs, 3898 total. This is validation of the proof, not its logical basis.

## 3. The first local Hamiltonian rearrangement occurs at 46

### 3.1 Square four-cycles require an additive collision

If four distinct positive integers form the cycle `a–b–c–d–a`, write its successive square edge labels as `A,B,C,D`. Then

```text
A+C = (a+b)+(c+d) = (b+c)+(d+a) = B+D.        (1)
```

The unordered pairs `{A,C}` and `{B,D}` must differ. Otherwise some adjacent edge labels agree, forcing two opposite vertices to coincide. Conversely, a nontrivial identity `A+C=B+D` produces the formal reflection cycle

```text
(t, A-t, B-A+t, D-t),                         (2)
```

but it gives a cycle of `Q_n` only when these four entries are distinct integers in `[1,n]`. This last condition is the required height and collision information; an additive identity alone is insufficient.

**Proposition (PROVED).** `Q_n` has a four-cycle if and only if `n>=46`. The unique four-cycle in `Q_46`, up to rotation and reversal, is `1–3–46–35–1`.

*Proof.* For `n<=46`, every edge sum is at most `2n-1<=91`; the possible square labels are exactly a subset of

```text
4,9,16,25,36,49,64,81.
```

The complete unordered pair-sum table (allowing repeated entries) has only one collision between distinct pairs:

```text
4+81 = 36+49 = 85.                           (3)
```

This is a finite eight-element arithmetic check, reproduced independently in the script. Equation (1) therefore requires the edge of square 4, whose distinct positive endpoints must be 1 and 3. The other incident labels are 36 and 49. The two possible unordered choices of the remaining vertices are `{35,46}` and `{33,48}`. Hence the smallest possible maximum vertex is 46, attained uniquely by `1–3–46–35–1`. This cycle remains present for all larger `n`. QED.

The first triangle occurs earlier, at `n=30`, uniquely on `{6,19,30}`. This is recorded here as **FINITE-EXACT** over `Q_1,...,Q_30`, not as a separate analytic classification. Odd cycles begin earlier still: `Q_15` already contains its five-cycle. Triangle birth, odd-cycle birth, and four-cycle birth are distinct events.

### 3.2 A sharp Hamiltonian switch threshold

A **two-edge switch** means that two undirected spanning structures have edge sets whose two set differences each contain exactly two edges. For paths, require the same unordered endpoint set. For cycles, every vertex already has degree 2.

**Lemma (PROVED).** Two Hamiltonian cycles, or two Hamiltonian paths with the same endpoints, differing by a two-edge switch have an alternating four-cycle in their symmetric difference.

*Proof.* Red edges are removed and blue edges are added. Equality of the degree at every vertex gives equal red and blue incidence. With exactly two red and two blue edges in a simple graph, a nonempty balanced symmetric difference can only be an alternating four-cycle. QED.

Consequently every Hamiltonian cycle of `Q_n`, `n<=45`, is isolated under two-edge switches; the same is true in each fixed-endpoint class of Hamiltonian paths. This does not preclude exchanges of three or more edges. The fixed-endpoint qualification is essential for paths.

**Theorem (PROVED, finite explicit boundary witness).** The least `n` admitting two Hamiltonian cycles of `Q_n` related by a two-edge switch is exactly **46**.

The lower bound follows from the four-cycle proposition. For the upper bound, the following is a Hamiltonian cycle of `Q_46` (close the final 15 back to 1):

```text
1,3,33,16,9,40,24,25,39,42,22,27,37,12,13,36,28,21,
43,38,11,14,35,46,18,31,5,44,20,29,7,2,34,30,6,19,
45,4,32,17,8,41,23,26,10,15.
```

Reverse the segment from 3 through 35. This removes `{1,3}` and `{35,46}` and adds `{1,35}` and `{3,46}`; every other edge remains. The result is another Hamiltonian cycle, and the new edge sums are 36 and 49. The script checks the vertex universe and every edge of both cycles, not just the four altered edges. A bounded discovery DFS found this witness in 652 nodes; reproduction verifies the explicit certificate directly.

**Hostile to sufficiency of a four-cycle:** in the cycle `1–2–3–4–5–6–1`, remove `{1,2},{4,5}` and add `{1,5},{2,4}`. The altered edges form an alternating four-cycle, and all vertex degrees stay 2, but the result consists of two triangles. Cyclic order / connectedness must be retained in addition to the additive identity.

The inherited census already has 11 Hamiltonian cycles in `Q_34` ([Hamiltonicity note, section 3](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md)). All 11 are isolated under these local switches. Having many global solutions therefore need not give a useful local deformation process.

## 4. The exact arithmetic bridge: square lifts and the signless kernel

Let `G=(V,E)` be any finite simple graph with distinct positive integer vertex values `x_v`, and designated positive integer edge roots `s_e` satisfying

```text
x_u+x_v=s_e^2                 for e={u,v}.
```

Fix a positive integer `q`. Seek new values

```text
x'_v=q^2 x_v+c_v
```

whose designated edge root is exactly `q s_e`. Then the necessary and sufficient equations are

```text
c_u+c_v=0                    for every e={u,v}.       (4)
```

**Proposition (PROVED).** Over a field of characteristic different from 2, the solution space of (4) has one free parameter for each bipartite connected component of `G` and zero parameters for each nonbipartite component. In a bipartite component the solution is `c_v=b` on one part and `c_v=-b` on the other. An isolated vertex counts as a bipartite component.

*Proof.* Each edge reverses the sign of the parameter. Along a path, the value is fixed by the initial value and path parity. An odd closed walk forces `b=-b`, hence `b=0`; without an odd cycle the bipartition makes sign transport consistent. QED. For integer offsets, the free component parameters must be integers. Positivity and distinctness of the new vertex values are additional inequalities, not consequences of the kernel calculation.

This is precisely the signless-incidence operator `c -> (c_u+c_v)_e`. It is the equal-exponent specialization of the two-support equations in [THM-4497](../../01-canon/theorems/THM-4497-coprime-graph-rank-compression.md), and the same linear operator as [THM-2521, equation (4)](../../01-canon/theorems/THM-2521-k13-drift-k14-potential-module-bridge.md). The meanings of the variables differ; the transport-of-kernel statement is exact.

### 4.1 The first path and the whole graph behave differently

For the whole graph `Q_14`, the tree gives one alternating parameter. `Q_15` is connected and contains the odd cycle `1–3–6–10–15–1`, so it has **no nonzero offset preserving every designated scaled edge label**. Since `Q_n` stays connected and retains this cycle, the same rigidity holds for all `n>=15`.

A selected Hamiltonian path is bipartite. Thus it always permits the alternating assignment `c_(v_i)=(-1)^i b`. For `q=2,b=1`, the path at 15 lifts to

```text
33,3,61,39,25,11,53,47,17,19,45,55,9,27,37.
```

Every consecutive sum is four times its original square. The original unused edge `{1,3}` lifts to values 3 and 11, whose sum 14 is not square. This is a concrete failure of extending the marked-path certificate to the whole graph.

For any integer `q>=2` and `|b|<q^2/2`, the alternating lift of distinct positive integer path labels is positive and injective: different original values differ by at least `q^2`, while their offsets differ by at most `2|b|<q^2`. Hence one obtains infinitely many embedded square-sum paths. Their vertex sets are sparse subsets of a larger interval, not the entire interval `1,...,N`. Dense coverage and junctions remain the work in a true Hamiltonian extension such as the inherited 25-fold construction. Repetition of an embedded motif supplies no occurrence-density estimate.

### 4.2 Correcting the edge roots: the missing compatibility data

More generally prescribe roots `q s_e+r_e`. The exact offset equation becomes

```text
c_u+c_v=d_e,
d_e=2q s_e r_e+r_e^2.                               (5)
```

Choose a root and a spanning tree in each connected component. Tree propagation writes every potential as

```text
c_v=sigma_v t+h_v,       sigma_v in {+1,-1}.
```

Every remaining edge then imposes

```text
(sigma_u+sigma_v)t=d_e-h_u-h_v.                       (6)
```

Equations (6) are a necessary and sufficient compatibility test over the rationals. A coefficient 0 requires its right side to vanish; a coefficient `+2` or `-2` fixes the root parameter, and all such fixes must agree. Integer offsets further require the resulting vertex potentials to be integral. Positive, distinct lifted values must still be checked.

In closed-walk language, even walks require the alternating sum of the `d_e` to vanish; odd walks determine twice the starting potential. It is not sufficient to check simple even cycles alone while forgetting compatibility between multiple odd cycles. The spanning-tree formulation retains that information automatically.

This is the exact form of a possible extension mechanism: assign edge-root corrections and solve their compatibility constraints, then verify injectivity, positive height, dense coverage, and gluing. It does not turn a local family of square paths into a global recurrence theorem by itself.

## 5. Connection ledger and stopping boundary

| Source → target | Map | Preserved predicate | Destroyed information / required sidecar |
|---|---|---|---|
| Square-labelled graph → affine lifted graph | `x_v -> q^2 x_v+c_v` | Every designated edge-square equation iff (4), or (5) for corrected roots | Positive integer values, no collisions, and complete interval coverage |
| Four square edges → additive arithmetic | Alternating sum of edge labels | Identity (1) | Actual four distinct vertices within one height cutoff; recovered by (2) |
| Hamiltonian cycle → two-edge exchange | Reverse a segment using the other two edges of a four-cycle | Vertex degrees and labels | A disconnected 2-factor can result; preserve cyclic ordering |
| Square-lift offset equations → valuation graph kernel | Identify `c_u+c_v=0` with equal-exponent two-support rows | Kernel dimension and odd-cycle grounding | General exponent ratios and all other prime/coordinate rows; THM-4497 retains them |
| Either graph model → actual Collatz orbit | **No trajectory-preserving map supplied** | Only the explicitly shared linear operator | `3m_i+1=2^{k_i}m_{i+1}`, chronological adjacency, sign/carry, absolute height, and descent |

No frequency theorem for long Collatz paths follows from the new square-sum switch threshold or the offset kernel. The useful common mechanism is more specific: an apparently repeatable local path may acquire incompatible equations when its unused edges or arithmetic relations are restored. Odd-cycle grounding gives one exact version of that obstruction. This also explains why the selected motif and its full ambient graph must be kept separate when discussing a recurring family.

## 6. Reproduction and status boundary

```text
python 04-computation/experiments/crossroads_family_20260926_squares.py
python -O 04-computation/experiments/crossroads_family_20260926_squares.py
```

The script does not import the earlier square-sum scripts. It uses exact integer square tests and an independent backtracking census with connectivity and certified endpoint pruning. The general unicyclic theorem is checked over all labelled connected unicyclic graphs on 3..6 vertices. Four-cycles are found by common-neighbour intersections, independently of the eight-square pair-sum proof. The two cycle witnesses at 46 are checked directly. Hostiles cover an actual degree-2 endpoint, a disconnected result of a four-cycle exchange, and failure of a selected-path affine lift on an unselected edge.

The first path at 15 is recovered prior work. The local-switch horizon 46 and the explicit lift/valuation operator connection are the contributions of this lane; no claim of external novelty is made. Literature claims about all sufficiently large square-sum Hamiltonian graphs are not needed here. Collatz remains open.
