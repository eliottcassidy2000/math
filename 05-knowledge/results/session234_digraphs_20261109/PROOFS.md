# Binary matrices, reciprocal arcs and the complete seven/twenty-one families

**Status: PROVED, with explicitly identified FINITE-EXACT inputs; independently
audited.** Local theorem labels DG1–DG8 belong to this file, not the global
theorem-number namespace. Classical ingredients are credited; no novelty
claim is made.

## 1. The object and the first failed implication

Let `A` be an `n by n` binary matrix. Its vertices are `0,...,n-1`, and
`A[i,j]=1` means the arc `i -> j` is present. Define

`H(A) = sum over permutations p of product A[p_i,p_(i+1)]`,

where the product runs over `i=1,...,n-1`. Paths are ordered, visit each
vertex exactly once, and have no prescribed endpoints. A valid reverse
order is another path. We do not count walks, cycles, or path orbits.
Use `H(empty)=1`; the empty-graph convention has no effect on 7 or 21.

For a pair `i<j`, there are four possibilities:

| `(A[i,j],A[j,i])` | meaning |
|---|---|
|`(0,0)`|absent pair|
|`(1,0)` or `(0,1)`|one directed arc|
|`(1,1)`|reciprocal pair, also called a digon|

A loopless **tournament** has one arc at every pair. A loopless
**semicomplete digraph** has one or two arcs at every pair. Thus zero
diagonal and no absent pair imply semicompleteness, **not** tournament
structure. This is the first failed implication behind the conjecture.

**DG1 (diagonal independence).** Changing any diagonal entries leaves `H`
unchanged: consecutive vertices of a permutation are different. In
particular a reflexive node cannot explain either total. Every labelled
off-diagonal matrix has exactly `2^n` diagonal decorations with the same
count. This does not imply a factor `2^n` for *isomorphism classes*.
Under uniform sampling of labelled binary matrices, even after conditioning
on a particular attainable H, the diagonal bits remain independent and
fair. Thus the probability of at least one loop is `1-2^(-n)`. A test
allowing any loop can hide the much smaller all-zero-diagonal sector.

## 2. Minimal witnesses and sharp corrected implications

The four-vertex matrix

```text
0111
0011
1001
1000
```

has exactly paths `0123,1203,1230,1302,2013,2301,3012`. Its two reciprocal
pairs are `{0,2}` and `{0,3}`. All diagonal entries are zero and every
unordered pair is covered. It refutes the seven-path conjecture.

The five-vertex matrix

```text
01111
00110
00011
10001
11000
```

has 21 paths, explicitly listed in [the catalogue](census/results.json).
It has two reciprocal pairs and likewise refutes the conjecture.

These orders are globally minimal even among unrestricted digraphs. An
order-three graph has at most `3!=6` paths. On four vertices a complete
bidirected graph has 24; removing any directed nonloop arc excludes the
`3!=6` permutations containing it consecutively, leaving at most 18.
Thus 21 cannot occur before order five. Loops change neither argument.

One reciprocal pair already suffices for 21 at order six:

```text
011110
001111
000111
100011
000001
100000
```

Its only reciprocal pair is `{0,3}`. The two tournament resolutions have
17 and 19 paths. Their union has 21, because 15 paths use neither direction
of that pair consecutively: `17+19-15=21`.

**DG2 (sharp repair, all orders).** For any binary adjacency matrix,

- `H=7` implies an absent unordered pair **or at least two reciprocal pairs**;
- `H=21` implies an absent unordered pair **or at least one reciprocal pair**.

The first follows from the complete strong-core classification below;
the second also follows directly from the tournament omission of 21.
The displayed matrices attain the respective reciprocal-pair lower
bounds. These are necessary alternatives, not sufficient conditions for
the exact totals. No diagonal condition is needed.

## 3. Why finitely many cores control every order

**DG3 (ordered SCC factorization).** Every semicomplete digraph has a
unique total ordering `C_1,...,C_r` of its strongly connected components,
and every cross-component pair points forward. Its path count is

`H(D)=product_i H(C_i)`.

Proof. Opposing arcs between two components would make them mutually
reachable. Semicompleteness makes the condensation a tournament, and a
condensation is acyclic, hence transitive. A path cannot leave a component
and later return. Therefore it visits the components as consecutive
blocks in their unique order; all component paths concatenate. This gives
a bijection, not only an inequality.

**DG4 (strong semicomplete Hamiltonicity).** A strong semicomplete graph
of order at least two has a directed Hamiltonian cycle; at order two a
two-cycle is allowed.

Proof. Choose a longest directed cycle `C`, initially supplied by strong
connectivity. If an outside vertex `v` cannot be inserted in any cycle
edge, it either strictly dominates all of `C` or is strictly dominated
by all of `C`. Indeed, if `c_i -> v`, noninsertability forbids
`v -> c_(i+1)`, so semicompleteness forces `c_(i+1) -> v`; propagate
around the cycle. If no `c_i -> v` exists, every `v -> c_i` exists.

Let `A` dominate `C` and `B` be dominated by `C`. If vertices remain,
strong connectivity forces both sets nonempty and a directed path from
`B` to `C`. Before first entering `C` it must cross from `B` to `A`.
Choose a crossing arc `b -> a`. Replace a cycle edge `c_i -> c_(i+1)`
by `c_i -> b -> a -> c_(i+1)`, contradicting maximality. Thus the cycle
is spanning.

For `n>=3` its edges occupy distinct unordered pairs. Retain those arcs
and choose one direction at every remaining pair. The resulting spanning
tournament is strong and is a subgraph of `D`. Consequently

`H(D) >= H(T) >= f(n)`,

where `f(n)` is the strong-tournament minimum. The order-two exception
has `H=2`, so it cannot divide either odd target.

**CITED + independently FINITE-EXACT input.** Moon's strong-subtournament
theorem gives a strong induced `(n-1)`-vertex tournament in a strong
`n`-vertex tournament, `n>=4`. Tournament path insertion gives
`H(T)>=H(T-v)`: insert `v` in one feasible position of each old path;
deleting it recovers that path. Thus `f` is nondecreasing. Independent
exhaustive labelled tournament scans through seven verify

`f(3),f(4),f(5),f(6),f(7) = 3,5,9,15,25`.

The order-seven scan covers all `2^21=2,097,152` tournaments. Therefore
strong seven-path semicomplete cores have order at most four, and strong
twenty-one-path cores order at most six. This is the explicit bridge from
finite computation to an all-order conclusion.

Reference: J. W. Moon, *On Subtournaments of a Tournament*, Canadian
Mathematical Bulletin 9(3) (1966), 297–301, Theorem 2,
[DOI 10.4153/CMB-1966-038-7](https://doi.org/10.4153/CMB-1966-038-7).
See the inherited
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
and [independent audit](audit/AUDIT.md).

## 4. Complete all-order construction rules

The exact core census is:

| target | order | strong classes | reciprocal-pair distribution | automorphisms |
|---:|---:|---:|---|---|
|3|3|2|0: one; 1: one|orders 3 and 1|
|7|4|3|2: three|all trivial|
|21|5|16|2: six; 3: ten|14 trivial; two of order 3|
|21|6|16|1: two; 2: eight; 3: six|all trivial|

Both strong three-path cores matter. They can be written as the triangle
`010/001/100` and the triangle plus a reverse arc `010/101/100`.
The added arc `1 -> 0` is unused: candidates `102` and `210` require
absent arcs. Both graphs have exactly `012,120,201`. Strong connectivity
does not make arc addition strictly count-increasing.

Let `X -> Y` denote the **ordered sum**, retaining the induced graphs
and adding all arcs from `X` to `Y`, and none back. Let `T_a` be the
transitive tournament on `a` vertices, with `T_0` empty.

**DG5 (complete classification, up to isomorphism).** Every loopless
semicomplete seven-path graph has exactly the form

`T_a -> G -> T_b`,

where `G` is one of the three order-four seven-path cores and `a,b>=0`.
Every loopless semicomplete twenty-one-path graph has exactly one of:

1. `T_a -> Q -> T_b`, with `Q` one of the 32 strong twenty-one-path cores;
2. `T_a -> G -> T_b -> R -> T_c`;
3. `T_a -> R -> T_b -> G -> T_c`;

where `G` is one of the three seven-path cores, `R` either of the two
three-path cores, and `a,b,c>=0`.

Proof. In an odd target, nontrivial SCC factors are at least three.
For 7 there is only one factor 7. For 21 the possibilities are one
factor 21, or factors 3 and 7. The size cutoffs and complete small census
give precisely the displayed options. Every remaining SCC has count one
and must be a singleton: a nontrivial strong component has a spanning
cycle and hence at least two paths. Uniqueness is DG3.

This gives counterexamples at **every order `n>=4` for 7** and **every
order `n>=5` for 21**, not merely a finite collection. Arbitrary diagonal
bits preserve these counts; arbitrary relabelings give the labelled
matrices. The standalone [generator](census/generator.py) embeds all cores.

**Exact reciprocal-count refinement of DG2.** There are no reciprocal
pairs between SCCs. Therefore a semicomplete graph with `H=7` has
**exactly two** reciprocal pairs. One with `H=21` has **one, two or three**:
the single cores have those counts, and a mixed word has `2+0` or `2+1`.
All listed possibilities occur. For an arbitrary binary matrix, the
stronger necessary alternatives are therefore an absent pair or `b=2`
for seven, and an absent pair or `1<=b<=3` for twenty-one. These remain
necessary rather than sufficient conditions, and loops are unrestricted.

## 5. Exact counts and loop symmetry

Let `u_h(n)` count loopless semicomplete isomorphism classes and `l_h(n)`
labelled matrices. Each term below contributes only when its lower bound
on `n` holds:

`u_7(n)=3(n-3)`, `l_7(n)=3(n-3)n!`, for `n>=4`.

For 21 add:

| source | starts at | isomorphism contribution | labelled contribution |
|---|---:|---:|---:|
|order-five core|5|`16(n-4)`|`(44/3)(n-4)n!`|
|order-six core|6|`16(n-5)`|`16(n-5)n!`|
|three-path and seven-path cores|7|`6(n-5)(n-6)`|`4(n-5)(n-6)n!`|

For example `u_21(7)=92`, `l_21(7)=423,360`; at order eight these are
148 and 5,268,480. These are consequences of the all-order classification,
not claims that all order-seven/eight semicomplete graphs were scanned.

The ordinary generating functions for isomorphism classes are

`U_7(x)=3x^4/(1-x)^2`,

`U_21(x)=(16x^5+16x^6)/(1-x)^2+12x^7/(1-x)^3`.

The exponential generating functions for labelled matrices are

`L_7(x)=U_7(x)`,

`L_21(x)=((44/3)x^5+16x^6)/(1-x)^2+8x^7/(1-x)^3`.

Here `[x^n]L_h=l_h(n)/n!`. A total-order SCC automorphism fixes each
position, so automorphism groups multiply across cores. The coefficient
`44/3=14+2/3` accounts for the two order-five cores with automorphism
group order three; the mixed coefficient is
`2*3*(1+1/3)=8`.

Loops are inert for paths but not for isomorphism. Independent orbit
enumeration of loop subsets gives 12 decorated three-path cores, 48
decorated seven-path cores, 480 decorated order-five twenty-one cores,
and 1024 decorated order-six twenty-one cores. A singleton gap is now an
ordered binary word. Thus the **loop-allowed isomorphism** generating
functions are

`U_7^loop(x)=48x^4/(1-2x)^2`,

`U_21^loop(x)=(480x^5+1024x^6)/(1-2x)^2+1152x^7/(1-2x)^3`.

For labelled matrices simply use `L_h(2x)`. For example each directed
triangle has four loop-subset orbits, not `8/3`; a quotient by symmetry
cannot be replaced by unqualified division.

## 6. Tournament resolution and exact arc response

**DG6 (weighted resolution identity).** If a semicomplete graph has `b`
reciprocal pairs and `T` runs over all `2^b` tournament resolutions, then

`sum_T H(T)=sum_(P Hamiltonian in D) 2^(b-k(P))`,

where `k(P)` is the number of reciprocal pairs used consecutively by `P`.
Proof: every used reciprocal pair forces one resolution choice; every
unused reciprocal pair leaves both choices. This double-counts pairs
`(T,P)` exactly. Consequently `max_T H(T)<=H(D)`, but `H(D)` is neither
the sum nor generally the average of the resolution counts.

For the displayed seven-path matrix the four tournament counts are
`1,3,5,5`, summing to 14. None is seven. The freedom to use different
directions in different Hamiltonian paths is the missing coordinate.

**DG7 (exact increment).** Write `E_D(S,a)` for the number of paths on
`S` ending at `a`, and `B_D(S,b)` for paths starting at `b`. If an absent
arc `a -> b` is added, then

`Delta H = sum_(S: a in S, b notin S) E_D(S,a) B_D(V\S,b)`.

Every new path uses this arc once; splitting there gives the two factors,
and concatenating them is inverse. The increment can vanish, including
in a strong graph. Multiple additions require updated endpoint arrays:
their interactions prevent summing the increments measured only in the
initial tournament. This is a precise construction law, not only parity.

## 7. Arbitrary matrices and a complete, if expensive, generator

The SCC catalogue above classifies the semicomplete sector and therefore
all loopless counterexamples to the proposed implication. It does not
classify the missing-pair sector by those same cores.

At fixed order all matrices can be generated by four pair states, followed
by independent diagonal bits. There are `4^(n choose 2)` loopless digraphs,
`3^(n choose 2)` loopless semicomplete digraphs and `2^(n choose 2)`
tournaments. The endpoint recurrence

`E({v},v)=1`,

`E(S,v)=sum_(u in S\{v}) A[u,v] E(S\{v},u)`,

`H(A)=sum_v E(V,v)`

counts any of them exactly, with no diagonal dependence.

**DG8 (exact path-support generator for unrestricted digraphs).** Fix
`n,h`. Choose a set `P` of exactly `h` vertex permutations and let `U`
be the union of their consecutive directed edges. If `H(U)>h`, discard
`P`; adding arcs cannot repair it. Otherwise `H(U)=h`. For each permutation
`Q` not in `P`, form the nonempty forbidden set `E(Q)\U` of optional
arcs. Choose any set `X` of optional arcs containing none of these
forbidden sets in full. Then `U union X` has exactly the path set `P`.
Conversely every graph with exactly `h` paths arises uniquely from its
actual path set and its optional arc set. Finally choose diagonal bits.

Proof. The paths in `P` are already present in `U`. A new permutation `Q`
becomes a path exactly when every edge in `E(Q)\U` is added. The stated
hypergraph-avoidance condition excludes precisely those events. This is
an exact finite generator, not a claimed efficient classification.

## 8. What remains open

For tournaments the diagonal is zero, every pair has exactly one arc,
and Rédei's parity theorem makes `H` odd. The all-order exclusions 7 and
21 remain valid. Whether every other positive odd integer is attained is
the Chan–Pak conjecture, **OPEN**; the present results neither prove nor
refute it. See Chan–Pak, *Computational complexity of counting coincidences*,
Theoretical Computer Science 1015 (2024), 114776, Section 6.5,
Conjecture 6.2, [arXiv:2308.10214v2](https://arxiv.org/abs/2308.10214v2),
which credits the 2016 bof/Gordon Royle discussion.

The general fixed-count semicomplete grammar, proved separately in
[GENERAL_GRAMMAR.md](GENERAL_GRAMMAR.md), is not a claim that every
positive integer is realized. It classifies each fibre whether empty or
nonempty. No tournament spectrum or unrelated research prize is closed.
