# Independent audit of the eight-frame palette/local-system bridge

2026-10-25. **PASS** for [fano_eight/RESULT.md](../fano_eight/RESULT.md),
sections 2--6, in its stated marked finite-palette and graph-constraint
scope. The reviewer did not author that proof or import its checker.
This approves the candidate's elementary propositions; it does not add
a unit-distance realization or a Hadwiger--Nelson conclusion.

Reviewed candidate SHA256 before status promotion:
`b588733745ba40e8b489262cee7b6625d46c2c7b7ce859755c0a3fa48233a662`.
Producer checker SHA256:
`998414c92c62bc4e266f5c7a6512b6dc058eb58ab19e80a841d8b902eb145a05`.

The final inheritance amendment correctly credits concurrent Session62's
twisted-pair packet for the local-coefficient theorem, section counts and
endpoint/conic cube comparison. The reviewer read that incoming source
directly from live origin/main. The additional tetrahedral reduction,
basis-retention index and common-fixed-point witness remain explicitly
separate. That amendment changes provenance, not the audited equations.

## 1. The eight-state carrier and its retained basis

The eight one-endpoint-per-pair choices form F2^3. Each of the six palette
labels is a coordinate facet, so this incidence is faithful only with the
three pair directions retained. The two sets of four parity-transversal
triples, together with the three central pencil lines, satisfy the Fano
pair axiom. Pair permutations preserve parity; endpoint flips change it
by exactly their total flip parity. No identification of the frame points
with the seven Fano *vertices* is made.

All fourteen affine planes give S(3,4,8): for three distinct binary points,
the fourth is their sum. An automorphism fixing an origin preserves this
fourth-point operation and hence addition, so the full affine-design group
is AGL(3,2), order 1344. A palette retains an unordered basis of three
directions; its stabilizer has order 48. The shear hostile really preserves
the block design while losing a designated cube edge. This is an explicit
lost-basis phenomenon, not a count-only warning.

## 2. The tetrahedral split is equivariant

With `c=111`, total parity p has `p(c)=1`. Therefore
`x=p(x)c+(x+p(x)c)` is a direct-sum decomposition into the trivial line
and the even subspace E. Every coordinate permutation fixes c and preserves
E. An affine transport `x -> pi(x)+a` consequently acts as

```
(eta,e) -> (eta+p(a), pi(e)+a+p(a)c).
```

The pair-permutation S3 acts faithfully on E's three nonzero vectors;
adjoining its four translations gives all 24 permutations of four points.
Thus the displayed split is the group isomorphism `C2 x S4`, not merely
an equality of orders. The independent checker constructs the 48 cube
automorphisms directly from all 8! permutations, then verifies this action
on all eight states and **all 2,304 group products**. The tetrahedron is a
four-object permutation carrier; no planar unit K4 is inferred.

## 3. The local-coefficient obstruction dimension is correct

For an oriented edge `u->v`, the affine compatibility equation is
`x_v+pi_e x_u=a_e`. Reversing the edge changes the permutation and gain
as specified. On a graph, every edge assignment is a degree-one cocycle,
because there are no two-cells. A compatible vertex assignment is precisely
a zero-cochain whose differential equals a. This justifies the H1 test.

This target is a section of the endpoint-choice bundle, not a trivialization
of all six transported palette labels. Nontrivial permutation holonomy can
fix an endpoint frame while obstructing the latter. The final amendment
retains this distinction explicitly.

The S3-equivariant split above holds at the level of local systems, so its
H1 sequence splits. In particular the kernel of completion parity is
exactly `H1(Gamma;E_rho)`, with no unmentioned extension obstruction.
For a connected graph with v vertices and e edges, rank-nullity gives

```
dim H1(E_rho) = 2e - rank(d0)
             = 2(e-v) + dim ker(d0)
             = 2(b1-1) + epsilon,
```

where epsilon is the dimension of the even vectors fixed by all linear
monodromies. Its values are two for trivial monodromy, one for a single
transposition subgroup, and zero for C3 or S3. These exhaust subgroups of
S3. The connected-tree case gives zero, as it must. On a unicyclic graph
with three-cycle holonomy the hidden space vanishes; parity is therefore
sufficient there. This positive special case is justified by the formula.

When compatible assignments exist, their differences form H0 of the full
local system. The trivial odd line contributes one dimension and E contributes
epsilon. Hence the number of global assignments is `2^(1+epsilon)`, with
the root value determining every other vertex by tree transport.

The independent probe does not rely only on spanning-tree-normalized
bouquets. It constructs the actual differential matrices on:

- a three-vertex triangle (`b1=1`);
- two vertices joined by three parallel edges (`b1=2`).

For each graph it exhausts all `6^3` linear edge assignments and all `8^3`
translation assignments: **110,592 systems per graph, 221,184 total**.
Direct root propagation counts every section. Independent F2 matrix ranks
give H0 and H1 for both the full and even modules. For every fixed linear
assignment, the ratio of completion-compatible gain systems to fully
compatible systems is exactly `2^[2(b1-1)+epsilon]`. Thus the graph-level
calculation checks the local-system bookkeeping before and after a choice
of spanning tree.

## 4. The common-fixed-point quantifier is retained

For the bit three-cycle pi, the two affine generators `pi(x)` and
`pi(x)+101` have fixed sets `{000,111}` and `{011,100}`. They each have
two fixed frames, one of each completion parity, but their fixed sets
are disjoint. The tetrahedral factor sees two different individually fixed
vertices; it has no common fixed vertex.

This is a hostile to checking each **chosen fundamental loop separately**.
It is not a claim that every word in their generated holonomy group has
some fixed point. The reviewed manuscript explicitly computes the generated
word `g2*g1^-1` as translation by 101, which has no fixed point at all.
The independent checker verifies the word as an eight-point permutation.
The manuscript's simple-graph subdivision of the bouquet preserves this
obstruction, but asserts no Euclidean realization of those seam equations.

## 5. Replays and limits

From this directory:

```
python3 check_palette_independent.py > palette-normal.json
python3 -O check_palette_independent.py > palette-optimized.json
python3 ../fano_eight/check.py > palette-producer-normal.json
python3 -O ../fano_eight/check.py > palette-producer-optimized.json
```

The independent checker gives **52,230 explicit checks**. Its two output
files are byte-identical, SHA256
`a809b957c7ed2055807c550009e6cebed3edf7879053b3fdd46331c920c04442`.
The producer replays are also byte-identical, SHA256
`7842b0e003457faf9538592bb7cf28223887a32c3be32dd39fbed4ef8b2debee`.
Neither checker uses optimization-removable assertions.

No mathematical discrepancy was found. During review the author retained
the explicit “chosen/generating loops” quantifier and the generated
fixed-point-free word; this is a clarification, not a correction of a
published claim. The remaining boundary is substantive: graph-colouring
interfaces must actually produce the specified affine frame equations
before any of these cohomology or tetrahedral tests applies to them.
