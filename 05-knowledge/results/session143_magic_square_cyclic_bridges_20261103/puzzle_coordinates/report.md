# Explicit projective-point / Conway-puzzle / C5 tournament bridge

**Status: PROVED by explicit finite permutations; FINITE-EXACT checker PASS.**
This is a bounded positive construction, not an originality claim, not a projective collineation, and not an identification with the Baby Monster. The earlier sporadic packet remains unchanged.

## 1. Projective coordinates and legal moves

Work in `PG(2,3)`, normalizing each nonzero triple so its first nonzero coordinate is one. Label the thirteen points:

| label | point | label | point | label | point |
|---:|---|---:|---|---:|---|
| 0 | `(1,0,0)` | 1 | `(1,0,1)` | 2 | `(1,0,2)` |
| 3 | `(1,1,0)` | 4 | `(1,1,1)` | 5 | `(1,1,2)` |
| 6 | `(1,2,0)` | 7 | `(1,2,1)` | 8 | `(1,2,2)` |
| 9 | `(0,1,0)` | 10 | `(0,1,1)` | 11 | `(0,1,2)` |
| 12 | `(0,0,1)` | | | | |

A line is the set of normalized projective points in a two-dimensional subspace of `F3^3`. It has four points. Given distinct `a,b`, write their unique line as `{a,b,c,d}` and put `mu(a,b)=(a b)(c d)`. This is a legal puzzle move when the hole is at `a`; it sends the hole to `b`.

This is the basic-game rule in the primary [Conway–Elkies–Martin paper, section 2.2](https://arxiv.org/html/math/0508630v1). Its closed legal sequences based at one hole give the `M12` puzzle group; arbitrary endpoint-compatible sequences give the `M13` pseudogroup. We use our own projective-coordinate labeling, not the paper's point labels.

**Convention:** each move acts on positions. In a chronological word `a0,...,ak`, the resulting permutation is `mu(a{k-1},ak) ... mu(a0,a1)`, with the rightmost factor applied first. The hole starts at point `0`.

## 2. Five-move closed word of order five

Take

`g = [0,1,3,2,4,0]`.

Its five legal elementary moves are:

| step | four-point line | elementary permutation |
|---|---|---|
| `0 -> 1` | `{0,1,2,12}` | `(0 1)(2 12)` |
| `1 -> 3` | `{1,3,8,11}` | `(1 3)(8 11)` |
| `3 -> 2` | `{2,3,7,10}` | `(3 2)(7 10)` |
| `2 -> 4` | `{2,4,6,11}` | `(2 4)(6 11)` |
| `4 -> 0` | `{0,4,8,10}` | `(4 0)(8 10)` |

Multiplication gives

`g=(1 4 2 12 3)(6 11 10 7 8)`, fixing `0,5,9`.

Thus `g` has order five and fixes the hole. Choose its orbit ordered from point `3`:

`V=(v0,v1,v2,v3,v4)=(3,1,4,2,12)`

with projective coordinates

`((1,1,0),(1,0,1),(1,1,1),(1,0,2),(0,0,1))`.

Then `g(vi)=v{i+1}` modulo five.

## 3. The regular five-tournament and its exact spectrum

Define `vi -> vj` precisely when `j-i` modulo five lies in `{1,2}`. In the row-permutation convention `P_ij=1` iff `g(vi)=vj`, its adjacency matrix is `A=P+P^2`:

```
0 1 1 0 0
0 0 1 1 0
0 0 0 1 1
1 0 0 0 1
1 1 0 0 0
```

Exactly one of `j-i` and `i-j` belongs to `{1,2}`, so every unordered pair gets exactly one directed edge. Every vertex has outdegree two. This is the regular five-vertex tournament. Its characteristic polynomial is

`x^5-5x^2-5x-2 = (x-2)(x^4+2x^3+4x^2+3x+1)`.

Equivalently its eigenvalues are `zeta_5^r+zeta_5^(2r)` for `r=0,...,4`. This supplies the requested **explicit projective points -> legal closed puzzle word -> cyclic orbit -> tournament** construction.

## 4. Six-move order-four normalizer: a genuine 5:4 family action

A second closed legal word is

`h=[0,1,7,12,9,4,0]`.

Its successive elementary permutations are

`(0 1)(2 12), (1 7)(4 9), (7 12)(6 8), (12 9)(10 11), (9 4)(1 7), (4 0)(8 10)`.

The corresponding lines are respectively

`{0,1,2,12}, {1,4,7,9}, {6,7,8,12}, {9,10,11,12}, {1,4,7,9}, {0,4,8,10}`.

Their chronological product is

`h=(1 4 12 2)(6 10 11 8)`, fixing `0,3,5,7,9`.

It has order four and satisfies `h g h^-1=g^2` on all thirteen points. On the chosen orbit, `h(vi)=v{2i}`. The twenty distinct permutations `g^a h^b`, `0<=a<5`, `0<=b<4`, form a subgroup isomorphic to `C5 semidirect C4 = 5:4`; their action on `V` is faithful.

**Crucial distinction:** this is not a twenty-element tournament automorphism group. The multiplier `h` sends connection set `{1,2}` through the four-member family

`{1,2} -> {2,4} -> {3,4} -> {1,3} -> {1,2}`.

Translations preserve each tournament; `h` changes the tournament, while `h^2` gives its converse. Thus the larger group acts on the **family of four tournaments**; it does not preserve the chosen orientation. In particular `h^2` explicitly witnesses self-duality of the selected tournament.

## 5. Puzzle permutations are not projective collineations

The line `{0,1,2,12}` is sent by `g` to `{0,3,4,12}`, which is not a line. Indeed the line through `0` and `4` is `{0,4,8,10}`. This gives a direct incidence obstruction.

Independently, `|PGL3(3)|=((27-1)(27-3)(27-9))/2=5616`, not divisible by five. Thus its elements cannot include this order-five permutation. The projective coordinates support the puzzle move rule; they do not turn its entire puzzle group into projective transformations.

## 6. Exact scope and reproducibility

The source is the explicitly labeled projective plane plus legal hole paths. The target is the cyclically ordered five-orbit and its four tournament orientations. The chosen subgroup retains `g` translation and `h` multiplication-by-two actions. Treating its output merely as a tournament discards the other eight positions, line incidence, and the chosen move-word history. No map from all of `M12` or `M13` to this fixed five-set is asserted: arbitrary puzzle permutations need not preserve it.

The Baby Monster connection remains only the already proved generic restriction-to-`C5` operator construction. No Baby Monster subgroup embedding, representation intertwiner, field identification, Hecke map, or canonical correspondence is established here. The generator, orbit, starting point, and connection set are explicit choices.

`python3 check_coordinates.py` executes **609 explicit checks**, including independent span/cross-product constructions for every line, all 156 nontrivial elementary moves, both complete words, the twenty-element subgroup, its faithful action, the tournament, its determinant, and the noncollineation witness. The output is byte-identical under `python3 -O`.

`python3 test_coordinates.py` passes eight positive/negative cases under normal Python and `-O`. It rejects an illegal single-transposition move, the inverse normalizer (which has the wrong conjugation exponent), and a reciprocal-edge adjacency mutation, each with nonzero exit status. The test harness itself also runs unchanged under `-O`. No optimization-sensitive assertions are used.

Optional `python3 search.py` reproduces the deterministic bounded discovery: the order-five search allows at most four internal nonzero hole vertices, and the conjugator search at most five. It finds the normalizer after 25,525 candidate closed words. This is a search over bounded words, not an enumeration of `M12`. Exact verification of the displayed witnesses requires no search and no external software or data.
