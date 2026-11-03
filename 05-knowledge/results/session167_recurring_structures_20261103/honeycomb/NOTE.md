# One honeycomb, two return laws: square reflections and tournament flips

**Status: PROVED elementary covering and constructive-cycle statements; FINITE-EXACT controls.** A single explicitly colored honeycomb covers both an arithmetic three-reflection graph and an actual three-coordinate tournament-flip cube. Their different deck subgroups explain which words return and which retain hidden state. A positive-integer square-sum cycle family realizes the arithmetic winding, with an exact two-lap tournament synchronization law. No Collatz convergence, complete metagraph classification, or literature-priority claim is asserted.

## Inheritance and the new composition

The inherited square-sum mechanism is the reflection `R_s(x)=s−x`, with **s the target sum**, not its square root. The first path at 15 and its extra chord are audited in [THM-4505](../../../../01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md) and [the square-sum crossroads note](../../crossroads_family_20260926_squares.md). The inherited distinction between a selected path and all square edges remains essential.

The fixed-path tiling source is an actual Boolean cube: free arcs are pinned to positions along a prescribed directed Hamilton path. Three chosen free-arc flips generate a literal 3-cube. This source and its difference from the converse-merged tournament node observation are established in the repository's Session 33 full-flip continuation note. [Session 166](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/00-navigation/S166.md) develops graph-cover transport and emphasizes that an unmarked tournament quotient need not be a covering.

What is added here is the explicit common honeycomb carrier, its two exact return subgroups, the shortest winding/synchronization formulas, and the mechanical positive square-cycle family. The Collatz counts below independently recover the earlier growth obstruction; they do not turn that obstruction into a dynamical conjugacy.

## 1. The colored honeycomb and its exact distance

Let `H` have vertices `(u,v,e)` with `u,v∈Z`, `e∈{0,1}`. Every vertex has one edge of each color `a,b,c`. From a vertex in layer zero:

    a: (u,v,0) ↔ (u,v,1),
    b: (u,v,0) ↔ (u−1,v,1),
    c: (u,v,0) ↔ (u,v−1,1).

These are undirected edges; the reverse moves in layer one add the corresponding coordinate. This is the ordinary hexagonal honeycomb. For example, place layer-zero vertices at `u P+v Q` and layer-one vertices at `u P+v Q+R`, where

    P=(sqrt(3)/2,3/2), Q=(−sqrt(3)/2,3/2), R=(0,1).

All three incident edge vectors have length one and the familiar hexagonal faces result.

For two vertices in layer zero,

    d_H((0,0,0),(u,v,0)) = 2 N(u,v),
    N(u,v)=max(|u|,|v|,|u+v|).                    (1)

Indeed two-edge moves on layer zero have displacements `±(1,0)`, `±(0,1)`, `±(1,−1)`, or zero. Each changes each of `u,v,u+v` by at most one, proving the lower bound. If `u,v` have the same sign, use `|u|+|v|` axial moves. If they have opposite signs, first use `min(|u|,|v|)` diagonal moves, then the remaining axial moves. This attains `N` in both cases. The ordinary honeycomb ball has `1+3r(r+1)/2` vertices at graph radius `r`. For completeness, at radius `2t` the black count is `1+3t(t+1)` by (1), and the white count is `3t²`; at radius `2t+1` they are `1+3t(t+1)` and `3(t+1)²`. The latter white region is the integer set `−t−1≤u,v≤t`, `−t−1≤u+v≤t`, counted by subtracting two corner triangles from a `(2t+2)`-square. This also gives shell size `3r`.

Following colors `abcabc` returns to the same honeycomb vertex. Colors `abab` instead reach `(2,0,0)` from the origin. Thus identical local involutions alone do not decide global return.

## 2. Three arithmetic reflections are a cylindrical honeycomb quotient

Fix distinct integers `a<b<c`; write

    p=b−a, q=c−a, g=gcd(p,q), k=(q/g,−p/g).

Let `B_(a,b,c)` have vertices `(n,e)∈Z×{0,1}` and colored edges

    (n,0) ↔ (s−n,1),  s∈{a,b,c}.

This is a three-regular bipartite graph. Even if `2n=s`, its edge joins different layers, so it is not a graph loop. It is the bipartite double cover of the full integer reflection system, before imposing positivity or deleting fixed-point edges.

For any starting integer `x`, define

    π_x(u,v,0)=(x+pu+qv,0),
    π_x(u,v,1)=(a−x−pu−qv,1).                    (2)

**Theorem 1.** The map `π_x` is a color-preserving graph covering onto the component with layer-zero residue `x mod g` and layer-one residue `a−x mod g`. Its fibers are exactly the orbits of translation by `k`. Consequently that component is `H/Zk`.

*Proof.* The endpoint sums of its `a,b,c` edges are respectively `a,b,c`, and the three distinct colors give a bijection of incident edges. Bézout's identity makes every integer with the stated residue reachable. Within a layer two coordinates have the same image exactly when `p Δu+q Δv=0`; coprimality after dividing by `g` gives `(Δu,Δv)=t k`. The layer is retained. These establish the asserted covering, fibers and quotient. ∎

A color word of even length starting at the origin has honeycomb endpoint `(u,v,0)`. It returns to its arithmetic state exactly when

    (u,v)∈Zk.                                    (3)

The zero vector means the return already closed in the honeycomb. A nonzero vector is additional arithmetic winding.

**Shortest-winding corollary.** The shortest additional arithmetic return has length

    L_arith = 2q/g.                              (4)

This follows from (1), because `0<p<q` gives `N(k)=q/g` and `N(tk)=|t|q/g`. A geodesic to `k` attains the bound. The projection is injective on every honeycomb ball of radius at most `q/g−1`, and fails on a ball of radius `q/g`.

For the last statement, two points in a smaller ball have distance less than (4), so cannot differ by a nonzero deck translation. A shortest path from a vertex to its translate has even length (4); its midpoint centers a ball of radius `q/g` containing both. Honeycomb color-preserving symmetries move that midpoint to any chosen vertex, preserving the deck relation up to sign. This proves the sharp radius for balls, not merely a lower bound.

Examples:

- Targets `4,9,16`: primitive period `(12,−5)`, extra-return length 24.
- Targets `9,16,25`: primitive period `(16,−7)`, extra-return length 32; injective balls through radius 15.
- Targets `16,25,36`: primitive period `(20,−9)`, extra-return length 40.

These are distinct roles of the numbers. None explains Hamiltonian thresholds by numerical coincidence.

## 3. The same honeycomb covers an actual tournament-state cube

Choose three distinct pinned free arcs of a fixed-path tournament, keeping all other arcs fixed. There are eight literal tilings, indexed by three bits; the three corresponding arc reversals are the coordinate flips. Such a choice is available from order four onward, since order four already has three free arcs.

Define

    χ(u,v,e)=(u+v+e, u, v) mod 2.                 (5)

An `a` edge flips the first bit, a `b` edge the second, and a `c` edge the third. The map is onto the eight bit states. Its same-layer kernel is `2Z²`, and bit-sum parity equals `e`, so opposite layers are never identified. Therefore

    H/(2Z²) ≅ Q_3                               (6)

as color-preserving graphs. This is a cover of **actual pinned tournament states**, not a claimed cover of tournament isomorphism classes or of the blue/black complement-line metagraph.

For an even color word with honeycomb endpoint `(u,v,0)`, tournament-state return is exactly `(u,v)∈2Z²`. Arithmetic return is (3), and honeycomb return is `(u,v)=0`. These are explicit different equivalence relations on the same move words.

**Theorem 2 (two-lap synchronization).** Among words returning to both the arithmetic state and the three pinned tournament bits, the shortest word with nonzero honeycomb displacement has length

    L_joint = 4q/g = 2 L_arith.                   (7)

*Proof.* The primitive vector `k` has relatively prime coordinates; at least one is odd. Hence

    Zk ∩ 2Z² = 2Zk.

Apply (1) to the smallest nonzero common translation `2k`. A geodesic to it has length `4q/g`. Equivalently any primitive one-winding arithmetic return has nonzero cube holonomy of order two, so repeating it closes the bits. ∎

For `a,b,c=9,16,25`, `abab` returns to the tournament but moves the arithmetic value by 14. The word `(ab)^9(cb)^7` has length 32 and honeycomb displacement `(16,−7)`; it returns arithmetically but flips tournament bits `(1,0,1)`. Repeating it has length 64 and returns to both. The word `abcabc` already returns in all three carriers.

The integer examples in this section may leave the positive integers. Their algebra is exact, but a word belongs to a finite square-sum graph only after every visited value, every distinct-edge condition, and the bound `1≤n≤N` are checked. The next construction does this for an infinite family.

## 4. Positive square-sum cycles with nontrivial tournament holonomy

Here is a general mechanical construction before specializing to squares. Retain `a<b<c`, `p,q,g`, and put

    r=q/g, s=p/g, D=a+b−c+g.

Suppose `D≥2` and choose an integer `1≤L≤D−1`. For `j=0,...,r−1`, set

    x_j=L+g (j s mod r),
    y_j= a−x_j  if x_j < L+q−p,
         c−x_j  otherwise.                       (8)

Read the cyclic sequence `x_0,y_0,x_1,y_1,...,x_(r−1),y_(r−1)`.

**Theorem 3.** This is a positive cycle in the bipartite arithmetic cover, with `2r` distinct layered vertices. Its first edge in each pair has target `a` or `c`, and its second has target `b`. It winds exactly once around the cylinder. If `2L≠a (mod g)`, forgetting the layer gives a simple positive-integer sum-graph cycle. Its largest vertex label is

    N_L=max(L+q−g, b−L).                         (9)

The smallest maximum within this translated mechanical family of layered cycles is

    N_min=ceil((b+c−a−g)/2),                    (10)

attained at `L=floor(D/2)` or `ceil(D/2)`.

*Proof.* Since `gcd(s,r)=1`, the `x_j` run once through `L,L+g,...,L+q−g`. If `x_j<L+q−p`, then `x_(j+1)=x_j+p` and `y_j=a−x_j`; otherwise `x_(j+1)=x_j+p−q` and `y_j=c−x_j`. In both cases `y_j+x_(j+1)=b`.

The `a` case yields the white-label progression from `b−L−q+g` to `a−L`. The `c` case yields the adjacent progression from `a−L+g` to `b−L`. Their union is one progression of length `r`, so all white vertices are distinct. Its minimum is `D−L≥1`; all black values are positive because `L≥1`. Black labels have residue `L mod g`, white labels residue `a−L mod g`; the stated residue condition makes the two sets disjoint after forgetting the layer. This proves the simple-cycle claim and (9).

There are `(q−p)/g` pairs `ab` and `p/g` pairs `cb`. Their honeycomb displacements are respectively `(1,0)` and `(1,−1)`. The total is `(q/g,−p/g)=k`, proving primitive winding. Finally the two affine expressions in (9) cross at `L=D/2`; because `D≥2`, the nearest integer choices lie in the allowed interval. Their common optimized maximum is (10). ∎

**Odd-square specialization.** For every integer `m≥3`, use

    a=(2m+1)², b=(2m+3)², c=(2m+5)².

Then

    g=8, r=2m+3, s=m+1, D=4m²−4m−7≥17.

Because `a` is odd while `g` is even, `2L≠a mod g` automatically. The construction gives a simple square-sum cycle of length

    4m+6

inside the ordinary square-sum graph on

    1,...,2m²+14m+13.                             (11)

It is a shortest additional winding in its three-target honeycomb quotient. Its tournament holonomy is nonzero, and two laps close all three pinned bits. No assertion is made that (11) is the smallest square-sum graph containing a cycle of that length, nor that these cycles are Hamiltonian in the whole interval graph.

For `m=3`, targets are `49,81,121`, `L=8`, and the cycle is

    8,41,40,9,72,49,32,17,64,57,24,25,56,65,16,33,48,73,8.

It has 18 edges, winds by `(9,−4)`, and flips bits `(1,1,0)`. Thirty-six moves, two laps, return to both arithmetic and tournament states. This supplies an actual positive square-sum realization of the otherwise formal deck calculation.

It can be observed even after forgetting the pinned labels. Start with the transitive four-vertex tournament and fixed path `0→1→2→3`. Associate colors `a,b,c` with the free arcs `02,03,13`. Read the displayed square-sum cycle as a word of arc flips. Its five `a`, nine `b` and four `c` moves reverse precisely `02` and `03`. The resulting tournament contains the directed triangle `0→1→2→0`, so it is not even isomorphic, up to converse, to the transitive starting node. A second lap restores the exact labelled tournament. Thus an arithmetic closed loop can induce a visibly nonclosed tournament trace, with no appeal to a purely formal hidden-state model.

## 5. The first square-sum path at 15 is a clipped honeycomb trace

Using only the three targets `9,16,25` on `1,...,15` gives exactly the path

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

Lift it from `(0,0,0)` by its target colors, with initial arithmetic label `x=8` in (2). Its endpoint is `(7,−3,0)`, whose arithmetic label is `8+7·7−3·16=9`. The full square-sum graph adds precisely the fourth-target chord `1–3`, whose sum is 4. That chord creates the odd cycle

    1–3–6–10–15–1.

This is a concrete way to see both the useful three-reflection carrier and the information it omits. The interval graph contains **all** allowed square targets, not merely the chosen three. The common value 15 here and in the injective-radius statement of §2 comes from different equations and is not promoted to a shared mechanism. Endpoint-aware Hamiltonian reasoning remains necessary; a honeycomb lift alone does not prove path existence or uniqueness.

## 6. What the Collatz tree needs that the honeycomb forgets

Use the ordinary Collatz map `C(n)=n/2` for even `n` and `3n+1` for odd `n`. The positive predecessors of `n` are `2n` and, when `n≡4 mod6`, `(n−1)/3`. The backward graph rooted at 5 is a tree: each vertex has a unique forward parent, and 5 is not periodic, as its explicit orbit enters `4→2→1→4` without returning to 5. Thus no vertex can occur at two different backward depths. Every vertex has at most two children, and its undirected degree is at most three.

Let `T_3` be the three-regular tree, represented by reduced words in the involutions `a,b,c`, with immediate repeats canceled. Any rooted tree with at most two children embeds in `T_3`: assign distinct available edge colors to children, avoiding the parent's color. Recursing gives distinct reduced path words and hence an isometric tree embedding. This applies constructively to the backward tree at 5, independently of universal Collatz convergence. Its edge-color choices are auxiliary; they are not the arithmetic parity rule.

The tree `T_3` covers `H` by following the reduced color word. Composing with (2) or (5) produces graph maps, but generally identifies different tree vertices. In particular `abcabc` is nonempty as a reduced tree word while its honeycomb image closes. The diagram therefore gives exact **information loss**, not a conjugacy of Collatz dynamics with square reflections or tournament flips:

    Collatz backward tree ↪ T_3 → H → arithmetic cover,
                                      ↘ pinned 3-flip cube.

An independent integer enumeration confirms the original capacity obstruction:

- Through depth 20, 557 positive integers reach 5; the honeycomb radius-20 ball has 631 vertices.
- Through depth 21, 704 reach 5; the honeycomb radius-21 ball has only 694.
- Through depth 26, 2263 reach 5; the triangular-lattice radius-26 ball has only 2107.

The honeycomb first fails this ball-capacity test at 21, and the triangular lattice first fails at 26. Any injective edge-preserving placement would send a depth-r tree ball inside a graph ball of radius r, so these are rigorous finite obstructions. They apply even without requiring nonedges to be preserved. Non-edge-preserving drawings, stretched edges, or a negatively curved tree carrier are different problems.

The useful constructive replacement is the path-word/tree carrier, followed by explicitly audited quotients. Negative curvature in the informal picture can mean this exponentially growing tree. No specific hyperbolic metric or canonical geometric Collatz embedding is needed or claimed.

## 7. Verification and use

Run `python3 check.py` and `python3 -O check.py`. Their JSON outputs must agree. All tests use explicit exceptions and exact integers, with no dependencies. The checker independently verifies:

- Every local arithmetic/cube edge map in the honeycomb ball of radius 32, and the distance/ball formulas there.
- All 88,573 color words of length at most 10 against direct affine arithmetic and independent bit flips.
- Shortest extra-return lengths and sharp injective-ball radii for four target triples, including a nontrivial gcd control.
- The positive cycle family for every `3≤m≤100`, with every label positive and distinct, every edge a square sum, exact heights, primitive winding and two-lap bit return.
- Every actual Collatz inverse vertex through depth 26, followed forward for its exact depth to recover 5; no depth repeats.
- The three-target path and fourth-target chord on `1,...,15`.

Finite checks validate the stated proofs; they do not extend their scope. The next useful operation is to retain the winding/parity label when transporting a square-sum trace to a pinned tournament experiment. Passing directly to an unmarked tournament node, or dropping the arithmetic layer and positivity filter, would erase the very return data this construction isolates.
