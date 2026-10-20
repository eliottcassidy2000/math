# A uniform three-black detour through every standard Tstar endcap

**PROVED:** the construction below works for every endcap order `r>=4`,
every middle tournament, and its chosen Hamiltonian path. It reproduces
the exact labelled target, not merely an isomorphic or converse state.
**FINITE-EXACT:** an exhaustive three-move search solves the ten-vertex
failed-guard example, and all fixed-path middle presentations through order
six validate the formula for endcap orders four through nine.
All-order black connectivity remains **OPEN**.

This continues the guarded-R4 theorem and its explicit hostile in
[`black.md`](black.md). The closest mechanism is its exact path-edge XOR
ledger; the previously missing boundary arc is the hostile to retain.
The less-used sidecar is the endpoint score charge: path reversal lowers
the first vertex's score by one, raises the last vertex's score by one,
and leaves all other scores unchanged. The new operation enlarges the
terminal ear instead of assuming the failed boundary arc exists.

## 1. The ten-vertex probe and its exact universe

Let `F5=Tstar_5` have rows `(2,4,9,19,7)`, as in the hostile. For
`U=F5^op -> F5`, with the converse-reversed labels, the standard path
`P=(0,...,9)` is blue. The earlier four-vertex ear construction fails
because its first move needs the absent arc `5 -> 7`.

The exact search enumerates all 72 black Hamiltonian paths from `U` and
all 240 black Hamiltonian paths from the literal target `R_P(U)`. For each
of the 17,280 pairs, the symmetric difference of those two path-edge sets
and `E(P)` is forced to be the middle move's edge set. It tests whether
that set is a directed Hamiltonian path in the first intermediate matrix,
then checks blackness and both literal matrix equalities.

This is exhaustive for three black path reversals. Reversing a path
preserves its blue/black predicate: its positional reversal preserves its
path-edge set, so toggling those edges commutes with the anti-automorphism
test. Therefore enumerating the third move backwards from the target
misses no valid third move.

The search finds exactly **18** three-black detours. The first is

```text
P1 = (0,1,2,3,4,6,7,8,9,5)
P2 = (5,9,6,4,3,2,1,0,8,7)
P3 = (7,8,0,1,2,3,4,5,6,9).
```

For this labelled ten-vertex target, three is the minimum possible number
of black path reversals. One move would have to reverse the original
path-edge set, whose only directed ordering is blue. Two moves cannot
work: each reverses nine edges, so their symmetric difference has even
cardinality, whereas the target differs in nine edges. This minimality
claim concerns literal labelled reversals, not distance after arbitrary
isomorphism or converse identifications.

## 2. The all-order construction

Let `Tstar_r` be the fixed-path code-zero tournament: consecutive positions
point forward, and every nonconsecutive pair points backward. Let

    T = Tstar_r^op -> Z -> Tstar_r,       r>=4,

where the left endcap has converse-reversed labels, and let `M` be any
Hamiltonian path of `Z`, possibly empty. The standard order on the left
endcap is a Hamiltonian path. Concatenate it with `M` to obtain `A`.
Write the standard right path as

    a, b, C, z,

where `C=(c,...,y)` is the nonempty list of its interior vertices from
position two through position `r-2`. Define

    P  = (A, a, b, C, z),
    P1 = (A, b, C, z, a),
    P2 = (a, z, b, rev(A), rev(C)),
    P3 = (C, A, a, b, z).                                  (1)

**Theorem.** All three successive moves in (1) are directed black path
reversals, and

    R_P3(R_P2(R_P1(T))) = R_P(T)                            (2)

as adjacency matrices on the original vertices. The target path `P` is
blue exactly when the chosen middle presentation `M` is blue. In
particular every such standard blue path has a three-black replacement.

### Directedness and exact parity

The needed initial arcs in the right endcap are its standard path arcs,
`z -> a`, and `z -> b`. The last two hold because `r>=4`. Every arc from
`A` to the right endcap is forward.

Let `x` and `v` be respectively the last and first vertices of `A`.
The first path traverses `A`, then `b,C,z,a`. After it is reversed:

- `a -> z` is the reversed last edge;
- `z -> b` is unchanged;
- `b -> x` is the reversed new boundary edge;
- `rev(A)` and `rev(C)` follow edges reversed by the first move;
- `v -> y` is an unchanged forward cross-edge.

These are exactly the links in `P2`. After the second move, `C` and `A`
have their original internal path directions, `y -> v` is reversed,
`x -> a` and `a -> b` are unchanged, and `b -> z` has just been reversed.
Thus `P3` is directed.

For the parity ledger, every internal path edge of `A` and `C` occurs
three times. The unintended edges `xb`, `vy`, `za`, and `zb` each occur
twice. The edges `bc`, `yz`, `xa`, and `ab` occur once. The surviving
edges are precisely the original path `P`, proving (2). This ledger
retains every internal arc of the arbitrary middle tournament; no quotient
or uniformity assumption is imposed inside `Z`.

### Blackness from exact endpoint charges

Put `m=|Z|`, so `n=2r+m`. In `Tstar_r` the original scores of `a,c,z`
are `1,2,r-2`. The first vertex of the left path has score one inside
its endcap. The three endpoint-score sums, in the tournaments where their
respective moves are made, are exactly

    m+r+2,       4,       r+1.                              (3)

Indeed P1 starts at the first left vertex and ends at `a`. It raises
`a` by one; P2 then has endpoints `a,c`, of scores two and two. P2 raises
`c` by one, so P3's endpoints `c,z` have scores three and `r-2`.

All three quantities in (3) are strictly below `n-1=m+2r-1` when `r>=4`.
A blue path's positional anti-automorphism would pair its endpoints with
scores summing to `n-1`. Therefore every move is black. The proof works
whether or not the middle path is blue.

Finally, standard positional reversal exchanges the two dual endcaps,
and all cross-arcs are in ordered-join direction. Its only unresolved
anti-automorphism condition is inside `M`; hence `P` is blue iff `M` is.
QED.

The case `r=4` recovers the earlier R4 row-zero construction. The case
`r=5` repairs the failed R4 boundary guard by using the full five-vertex
endcap. The lower boundary is substantive: for `r=3` and empty middle,
the first rotated-endcap path is still blue because rotation is an
automorphism of the directed triangle.

## 3. The precise general endcap criterion

The same macro applies more widely, with an essential symmetry guard.
Let `B` be an endcap on `r>=4` vertices with a chosen directed Hamiltonian
path `Q=(a,b,C,z)`. Require both the closing arc `z -> a` and the chord
`z -> b`. The closing arc already implies that B is strong. Let

    T=B^op -> Z -> B,

use the mirrored left path, and choose a blue middle presentation, so
the original path `P=(A,a,b,C,z)` is blue. Write `rho` for the cyclic
rotation taking the ordered list `Q` to `(b,C,z,a)`.

**General endcap theorem.** The three paths in (1) are successive directed
path reversals and reproduce the same literal target (2). The second
and third moves are always black. The first is black **if and only if**
`rho` is not an automorphism of B. Consequently the macro is a three-black
replacement exactly under this rotation guard.

**Proof.** Directedness and the edge ledger in section 2 use only the
chosen path, the closing arc, the chord, and the ordered-join cross-arcs,
so they apply unchanged. In their respective intermediate tournaments,
the last two endpoint sums are

    d_B(a)+d_B(c)+1,       d_B(c)+d_B(z)+1,

where `c` is the first vertex of C. For any distinct vertices of an
r-vertex tournament, their score sum is at most `2r-3`, because one loses
their mutual edge. Both displayed sums are therefore at most
`2r-2<=n-2`, which proves blackness.

For P1, the original and changed path orders agree on the entire left
endcap and middle, and differ only by rotating the right endcap. The
original blue positional reversal identifies the left endcap with the
converse of B. A second blue reversal for P1 exists exactly when the
changed right-endcap presentation has the same adjacency matrix, which
is precisely the condition that rho is an automorphism of B. If rho is
an automorphism, extending it by the identity on the rest of the ordered
join proves the converse implication as well. QED.

In particular any even r automatically excludes this rotational symmetry:
a full cyclic rotation acting by automorphisms would make B vertex
transitive and regular, which is impossible at even order. Nonregular B
also excludes it. These are sufficient shortcuts; the theorem retains
the exact rotation test and both explicit boundary arcs.

The **regular cyclic five-vertex tournament** is the essential hostile.
Give vertex i arcs to `i+1,i+2` modulo five, and choose its standard cyclic
path. Then both `z -> a` and `z -> b` hold, but rho is an automorphism.
The first move is blue, while the last two are black. Thus the closing
arc and chord alone do not prove a three-black replacement. This regular
tournament is different from the nonregular `F5=Tstar_5` repaired above.

## 4. Scope, checks, and next obstruction

This covers standard endcap Hamiltonian paths for every `Tstar_r`, not
every Hamiltonian presentation of every larger strong endcap. It removes
an infinite family from the remaining blue-edge obstruction, but does
not show that all terminal basins are connected by black moves.

The construction's boundary ledger also makes its remaining input precise:
a path `a,b,C,z` needs `z -> a`, `z -> b`, the appropriate prefix boundary
arcs, and a blackness certificate. For `Tstar_r` all are supplied uniformly.
An arbitrary endcap can fail these requirements, so replacing it by this
family would lose information.

Run from the repository root:

```text
python 05-knowledge/results/deficiency_three_20261020/ear_probe.py
python -O 05-knowledge/results/deficiency_three_20261020/ear_probe.py
```

The checker retains full adjacency matrices, tests every actual path
direction, every positional blue predicate, the unordered edge XOR, and
the final literal labelled target. After the exhaustive F5 search, the
first opaque-middle candidate succeeds. All `1101` fixed-path middle
presentations of orders zero through six pass for each `r=4,...,9`, giving
`6606` family controls and exact agreement with (3). The directed-triangle
boundary is checked as a hostile. These tests support the general proof;
the universal quantifier comes from the edge ledger and endpoint charges.

The additional general-macro replay covers every fixed-path endcap of
orders four, five and six satisfying the closing arc and chord, with
middle orders zero, one and two. There are respectively `2,16,256`
eligible endcaps; `0,1,0` have the rotational automorphism. All 822
expanded cases reproduce the exact target and have colour pattern
`(rotation_is_automorphism, false, false)`, where true means blue.
