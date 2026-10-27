# A unique minimum cut whose entire geodesic cone is blue

**Session75, 2026-10-27. PROOF CANDIDATE for independent audit.**
No shared repository edits. This packet refutes the proposed all-order equality
`black distance to nonstrong = cut connectivity`. It does not refute the separate
`black distance to singleton + strong core = minimum semidegree` proposal.

## Scope and inheritance

The operation is literal Hamilton-path reversal `R_P`, not single-edge reversal.
A directed Hamilton word is blue exactly when its positional reversal is an
anti-automorphism, and black otherwise. Its literal nonpath-complement `K_P`
version is obtained by alternating converses; the checker retains every named
word and both literal `K` endpoint tournaments.

The pinned inherited source is repository commit
`96e35642110104cc547516a63f37eebcf930f707`. The main dependencies are the signed
cut current and endpoint degree current, already reconstructed in
`/tmp/session75-black-core/BLACK_CORE_GEODESICS.md`, Sections1–2:

`c_(R_P T)(S)-c_T(S)=1_S(last P)-1_S(first P)`.

We read `black_descent_20261020c/cycles.md` but do not use its descending-path
theorem as an endpoint-selection surrogate. No reserved theorem stub is used.
The fixed-shore version of the desired black lemma already fails at order4;
Section7 below records this smaller obstruction to the *method*.

Concept board: minimum-cut family; rooted Hamilton-word orbits; literal reversal
colour; one-back-arc block joins; path-arc rather than score-state changes;
separate nonstrong/core terminal predicates.

## 1. Rooted-converse obstruction criterion

Let `U` be a strong tournament on `m>=3` vertices with distinguished root `r`.
Let `V` be a labelled copy of `U^op` with corresponding root `r_V`. Form `J` by
putting all arcs `V -> U`, except replace `r_V -> r_U` by `r_U -> r_V`.
Assume:

1. `Aut(U,r)` acts transitively on the directed Hamilton words of `U` ending at
   `r` (a rooted word-orbit hypothesis, not an unrooted purity hypothesis).
2. Every nonempty `X subset U\{r}` has at least two outgoing arcs in `U`.

Then `J` is strong, its unique minimum directed cut is `U` with value1, and
**every Hamilton word decreasing this minimum cut is blue**. Therefore no black
move reaches a nonstrong tournament from `J`.

### Strongness and uniqueness of the cut

Strongness follows from strong `U` and `V`, the bridge `r_U -> r_V`, and the
remaining arcs `V -> U`.

A proper nonempty subset of `U` not containing `r` has at least two outgoing
arcs by assumption2. One containing `r` has at least one internal outgoing arc
because `U` is strong, plus the exceptional bridge. The whole shore `U` has
exactly its one bridge outgoing.

A proper nonempty subset of `V` has at least `m-1>=2` outgoing arcs to `U`.
The whole `V` has `m^2-1>=2`. If a shore has nonempty proper intersections with
both blocks, the two internal strong cuts each contribute at least1, hence its
cut is at least2. Finally, the anti-automorphism exchanging the corresponding
copies sends `c_J(S)` to `c_J(V(J)\theta(S))`; this reduces shores containing a
whole block and part of the other to the already checked one-block subsets.
Thus `U` is the unique cut of value1.

### Completeness and colour of the lower-bound cone

A Hamilton word starting in `U` and ending in `V` must use one more `U -> V`
transition than `V -> U` transitions. There is only one possible `U -> V` arc.
Thus it uses that bridge exactly once, no reverse transition, and has the form

`P = alpha_U, beta_V`,

where `alpha` ends at `r_U` and `beta` starts at `r_V`. Identify the reversed
word `beta` with a word in `U` ending at `r`. Assumption1 supplies a root-fixing
automorphism taking the ordered word `alpha` to this reversed word.
Consequently the positional reversal of `P` anti-isomorphically exchanges the
two induced blocks. It also respects every cross pair: the cross orientation
is uniform except for the root-root pair, whose two roots are paired by the
reversal. Thus `P` is blue, with the full anti-automorphism predicate, not merely
an endpoint score equality.

Any nonstrong result after one reversal would have a zero cut. The cut-current
law forces that cut initially to have value1 and the path to start inside it
and end outside. Since `U` is the unique minimum shore, the preceding list is
complete. No black one-step exit exists.

This general criterion gives an obstruction, not by itself a finite black
distance. The next family supplies actual two-black routes.

## 2. Infinite family with exact black distance2

For each integer `a>=1`, let

`U_a = C3(T_a, C3, {r})`,

with module order `A -> B -> r -> A`, where `A=T_a` is transitive on `a`
vertices and `B=C3`. Let `J_a` be its rooted-converse join from Section1.
Then

- `|J_a| = 2a+8`;
- `lambda(J_a)=1`, with unique minimum shore `U_a`;
- `delta(J_a)=2`;
- exactly nine directed Hamilton words decrease the minimum cut, all blue;
- **black distance from `J_a` to nonstrong is exactly2**.

The family begins at order10, not order14. Both blocks are strong; neither is
required to be a strong pure-blue atom.

### Verification of the criterion

Deleting `r` leaves the order join `T_a -> C3`. Every Hamilton word ending at
`r` is the unique increasing word through `A`, then one of the three cyclic
rotations through `B`, then `r`. These three words form a single orbit under
rotation of `B`, fixing `A` and `r`. There are therefore `3*3=9` decreasing
words in the join, with no hidden alternatives.

For a nonempty subset `X` avoiding `r`, if `X` contains at least two vertices
of `B`, their arcs into `r` already contribute two. If it contains exactly one
vertex of `B`, that vertex's arc into `r` and its cyclic successor arc into
`B\X` contribute two. If it contains no vertex of `B`, each of its vertices
in `A` beats all three vertices of `B`. This proves the rooted cut margin.

In the full join the `U`-side outdegrees are

- `a+2-i` for `A` label `i=0,...,a-1`;
- `2` at each of the three `B` vertices;
- `a+1` at `r` (including the bridge).

They have minimum2. The opposite block has complementary full scores under
the block-exchange anti-automorphism, proving minimum semidegree2.

## 3. Two literal black words at every order

Set `m=a+4`. Give `V` the labels `0,...,m-1`; give `U` corresponding labels
`m,...,2m-1`. Internally, `A=0,...,a-1`, `B=(a,a+1,a+2)` is the directed
cycle in that order, and `r=a+3=m-1`. An internal label is shifted by `m` in
`U`. Let `alpha=(0,1,...,m-1)`, a Hamilton word of `U` ending at its root.
Its reverse is a Hamilton word of `V` beginning at its root.

The first word is

`P = (a+1,a,a-1,...,0, m,m+1,...,2m-1, r,a+2)`.

It cuts the Hamilton cycle `alpha_U, reverse(alpha)_V` at the internal `V`
arc `(a+2)->(a+1)`. Its two endpoints lie in `V` and each initially has full
outdegree `2m-3`. Their score sum `4m-6` differs from `2m-1`, so `P` is black.
Both endpoints lie outside `U`, so its unique minimum cut stays1. Every other
initial cut was at least2 and changes by at most1; hence the intermediate
`tournament J'_a` is strong.

The only current `U -> V` arc is now `m -> 0`. The second word is

`Q = (2m-1,2m-2,...,m, 0,1,...,a-1, a,a+2,r,a+1)`.

Its first block is `reverse(alpha)_U`, directed because `P` reversed all arcs
of `alpha_U`. In the second block:

- `0,...,a-1,a` is directed by those reversed path arcs;
- `a -> a+2` is the unchanged third cyclic edge in `V`;
- `(a+2) -> r` is reversed from the original path arc `r -> (a+2)`;
- `r -> (a+1)` is an unchanged original `V` arc.

Together with the new cross arc `m -> 0`, these verify every arc of `Q`.

### Why the second word is genuinely black

Write `U'=R_alpha(U_a)`. If the first move had reversed *every* arc of the
reverse-alpha word inside `V`, its induced result would be `(U')^op`.
Instead it omitted exactly the edge `(a+2)->(a+1)`. Thus the actual `V'`
is obtained from `(U')^op` by flipping one arc between the two `B` vertices
`a+1` and `a+2`.

Both have outdegree2 in `U'` (they are internal vertices of `alpha`), hence
both have outdegree `m-3=a+1` in `(U')^op`. Flipping their arc changes the sum
of internal squared outdegrees by exactly2. Therefore `V'` is not isomorphic
to `(U')^op`.

The word `Q` visits all of `U'` followed by all of `V'`. Its positional
reversal would exchange the two whole blocks and would require exactly that
anti-isomorphism. This is impossible, so `Q` is black. Since it starts in `U`
and ends in `V`, it lowers the sole outgoing cut from1 to0. The terminal is
nonstrong. This proves exact distance2, with a strong intermediate node.

No claim that this particular terminal has a singleton strong core is made.
The separate core-distance question is not settled by this family proof.

## 4. First member and complete exact certificate

At `a=1`, the rows are

`(1008,1001,995,997,494,448,640,768,576,48)`.

Its full scores are

`(6,7,7,7,7,3,2,2,2,2)`.

The unique minimum cut is mask992, labels `{5,6,7,8,9}`. There are657 directed
Hamilton words, of which558 are black. The checker explicitly reverses all558
and checks that all results remain strong. It independently enumerates all
1022 nonempty proper cuts and all root-ending block words. Exactly nine
Hamilton words decrease the unique minimum cut; each full anti-automorphism
equation holds. Thus the failed one-step search is a complete finite
exclusion, not a budget-limited failure to find a witness.

The all-order construction gives these two black words:

`P=(2,1,0,5,6,7,8,9,4,3)`;

`Q=(9,8,7,6,5,0,1,3,4,2)`.

The producer parent has also independently obtained a distinct two-black
route to the stronger singleton-core target at this first member; it is not
a dependency of the present theorem.

## 5. Smallest possible strict cut/semidegree separation at lambda1

For any tournament with minimum semidegree at least2, a shore of size `s`
has outgoing cut at least

`2s - binom(s,2) = s(5-s)/2`.

For `s=1,2,3,4` this is respectively `2,3,3,2`. Applying the same bound to the
complement in the opposite orientation shows that a cut of value1 must have
both shores of size at least5. Hence order10 is the smallest possible order
with `lambda=1<delta`.

This is a structural minimality statement for that parameter pattern. Global
minimality of the black-distance counterexample additionally requires the
parent's complete lower-order checks; this packet alone does not claim it.

## 6. Portable checks and exact scope

Commands:

```sh
python3 -B check_family.py > normal.json
python3 -B -O check_family.py > optimized.json
cmp normal.json optimized.json
```

The checker is stdlib-only and has no producer-module imports. All rejects are
explicit exceptions, not `assert`. It constructs and checks `a=1,...,20`,
orders10 through48, retains all nine decreasing words at every order, verifies
actual anti-equations, both black `R` words, actual intermediate strongness,
terminal nonstrongness, and literal named `K` conversions with colour at both
ends. At order10 it additionally exhausts the complete Hamilton-word universe
and all directed cuts.

For larger orders the minimum cut is checked exactly, not sampled, using
`c(S)=sum_(v in S)d(v)-binom(|S|,2)`: the minimum cut of a fixed size is the
sum of that many smallest scores minus the binomial term. A unique minimizing
size and a strict score gap at its boundary certify the unique shore. Scores
are used only for cuts and a valid obstruction to an anti-isomorphism; they
are never substituted for the actual path-colour test.

The two normal/-O JSON receipts are byte-identical. The proof is all-order;
the finite examples are controls and not its logical basis.

## 7. Why an arbitrary fixed minimum shore was already unsafe

At order4, rows `(8,1,3,6)` have `lambda=1`. Shore `{1}` has exactly one
decreasing Hamilton word, `(1,0,3,2)`, which is blue, whereas another minimum
shore has a black exit. At order6, the one-back-arc join of two `C3` blocks,
rows `(50,60,57,17,32,8)`, has natural minimum shore `{3,4,5}` and its only
decreasing word `(4,5,3,0,1,2)` is blue. Other minimum shores again permit
black exits.

The new family removes that escape hatch: its complete initial minimum-cut
family has just one member. The lost coordinate in the failed all-colour-to-
black transfer was the rooted word orbit at the only legal crossing ports,
not cut capacity or endpoint charges alone.
