# Which ternary tree lives on the fruit elliptic curve?

**PROVED:** the subgroup tree, finite symmetry obstruction, inverse guards,
and real positivity obstruction below. **CITED:** the full rank-one
Mordell--Weil description. **FINITE-EXACT:** the explicit integers, curve
operations and finite controls. No new rank, saturation, leastness, or
Collatz result is claimed.

There is a genuine algebraically defined ternary tree rooted at the large
fruit solution. Its vertices are positive **multiples of a point**, not
positive fruit triples. All three immediate children fail fruit positivity.
The curve's rational inverse-tripling tree is a different object and is
finite at this seed. These distinctions are essential to an isomorphism
claim about the primitive Pythagorean tree.

## 1. Inheritance and exact recovery

The closest proved mechanism is the projective map, torsion and positivity
audit in [the inherited fruit note](catalan_elliptic_20260921_elliptic.md).
The Pythagorean comparison inherits
[THM-3756, odd-square ordinal Berggren affine descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md)
and [the odd-square triangle note](odd_square_20260921_triangles.md).
The canonical hostile is a valid signed fruit solution that is not positive.
The corrected near miss is replacing two decimal-transcribed coordinates
by their intended values without recording the correction. The less-used
sidecar here is the free group coefficient together with its torsion class.

Anchor: an operation-defined tree at the actual seed. Niche: rational
division by three and its lost torsion coordinate. Wildcard: the separate
power-gap curve arising from `2^10=10^3+24`.

| Live object | Operation tested | Required invariant / first failure |
|---|---|---|
| Primitive Pythagorean triples | Three Berggren matrices | Unique parent and primitiveness |
| Fruit cubic | Same three matrices | Fails the cubic equation |
| Fruit coordinate symmetries | Permutations and pair-sum reciprocation | Finite orbit of size 12 |
| Infinite cyclic elliptic subgroup | `P -> 3P+iH` | Free coefficient gives a ternary tree |
| Rational division by three | Solve `[3]P=R` | Free divisibility and torsion congruence |
| Positive fruit solutions | Fixed expanding affine elliptic map | Proper real arcs cannot be globally preserved |

The literal earlier input was

```text
A = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
B = 368751317941299998271978115652254748254929799689719709962831374716372246340555790
C = 43736126779286972578612526023713901528165375581616136186214379933784234677720360.
```

The valid primitive positive solution is `(a,b,c)=(A,B/10,C/10)`, with
coordinate digit counts `(81,80,79)`. The literal `(A,B,C)` does not solve
the equation. Both facts are replayed exactly, using

\[
 F(a,b,c)=(a+b+c)^3-6(a+b+c)(ab+ac+bc)+7abc.
\]

When the three pair sums are nonzero, `F=0` is equivalent to
`a/(b+c)+b/(a+c)+c/(a+b)=4`. The projective linear correspondence is

\[
\begin{aligned}
[X:Y:Z]&=[-28(a+b+2c):364(a-b):6(a+b)-c],\\
[a:b:c]&=[56Z-X+Y:56Z-X-Y:-56Z-12X],\\
 E:&\quad Y^2Z=X^3+109X^2Z+224XZ^2.
\end{aligned}                                                     \tag{1}
\]

The two linear maps compose to `728` times the identity. Substitution
gives the elliptic residual `529984 F`; these identities are inherited
polynomial proofs, not conclusions from point sampling. The repaired
large point is exactly `H=9G`, where `G=(-4,28)`.

**CITED classification.** Bremner--Macleod's [2014 paper](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf),
Remark 2.2 and Section 3, identifies G as generator and the N=4 rank as one.
Together with its Lemma 2.1 this gives
`E(Q)=Z G ⊕ <T>`, with `T=(56,728)` of order six. The
[2025 corrigendum](https://publikacio.uni-eszterhazy.hu/8862/1/AMI_62_from26to27.pdf)
retains the N=4 ninth-multiple entry and the 79/80/81 digit sizes. We use
the full group description only where explicitly indicated; no independent
rank/saturation or global leastness calculation is claimed here.

The weaker facts needed for the forward tree are independently inherited:
the six rational torsion points are
`O,(0,0),(4,±52),(56,±728)`; reductions of orders 12 and 18 at 11 and 17
bound the torsion by six; and G is outside that list, hence has infinite
order. Write `Q=2T=(4,52)` and `U=3T=(0,0)`.

## 2. Two apparent three-branch constructions that fail

The standard Berggren matrices, in odd-leg/even-leg/hypotenuse order, are

\[
L=\begin{pmatrix}1&-2&2\\2&-1&2\\2&-2&3\end{pmatrix},\quad
M=\begin{pmatrix}1&2&2\\2&1&2\\2&2&3\end{pmatrix},\quad
R=\begin{pmatrix}-1&2&2\\-2&1&2\\-2&2&3\end{pmatrix}.
\]

They do not preserve the fruit cubic. Already the valid small signed
solution `(11,4,-1)` is sent by L to `(1,16,11)`, whose residual is
`F=-10920`. The other two matrices fail as well; the companion checks
all three also at the large positive seed. This refutes literal reuse of
the matrices, not the possibility of some different tree morphism.

A coordinate-replacement interpretation also needs a guard. Holding b,c
fixed leaves a cubic in a, not a quadratic. Once one root a is known,
the other roots have discriminant

\[
 -3a^2+6a(b+c)+5b^2+30bc+5c^2.                            \tag{2}
\]

At `(11,4,-1)`, the three coordinate choices give `-200,472,1912`;
none is a rational square. Thus replacing one coordinate by another
rational root does not supply three universally available branches.

The actual coordinate permutations and pair-sum reciprocation do preserve
the cubic. They give the inherited group `S3×C2`: cyclic permutation is
translation by Q, interchange of a,b is negation, and the central
reciprocal involution is translation by U. The orbit of any nontorsion P
under this group is exactly

\[
 \{\,\pm P+kT:0\le k<6\,\},                              \tag{3}
\]

with 12 distinct points. Equality of two of them would either identify
two torsion points or make `2P` torsion. Therefore these symmetries alone
cannot grow an infinite tree. At H, six images are positive coordinate
permutations and the other six have mixed signs. The six torsion points
themselves have forbidden zero pair sums, so cannot be used as fruit seeds.

More generally, finitely many maps `P -> ±P+R_i` have only polynomial
word growth. A length-h word has the form `±P+sum_j c_j R_j`, with each
`|c_j|≤h`, so s generators yield at most `2(2h+1)^s` images. Consequently
three fixed maps of this form cannot generate a full ternary tree with
distinct vertices at every level. This argument does not assume rank one.

## 3. A genuine ternary tree on the infinite cyclic subgroup

**PE1, PROVED independently of the full rank computation.** Let H be any
infinite-order rational point on E. On `{nH:n≥1}`, use the three branches

\[
 \Phi_{-1}(P)=[3]P-H,\qquad \Phi_0(P)=[3]P,\qquad
 \Phi_1(P)=[3]P+H.                                       \tag{4}
\]

Their orbit from H is an exact rooted labelled ternary tree containing
every positive multiple of H once. The free coefficient changes by
`n -> 3n+i`, for `i=-1,0,1`. The three image sets are disjoint: their
coefficients are respectively 2,0,1 modulo three. Every `n≥2` has the
unique parent and label

\[
 p=\lfloor(n+1)/3\rfloor,\qquad i=n-3p\in\{-1,0,1\},
 \qquad 1\le p<n.                                       \tag{5}
\]

Iterating (5) ends at one. Infinite order turns coefficient equality into
point equality, proving injectivity, disjoint branches and exhaustive
coverage of the stated subgroup half. A word of length h has coefficient

\[
 n=3^h+\sum_{j=1}^{h}i_j3^{h-j};
 \quad \frac{3^h+1}{2}\le n\le\frac{3^{h+1}-1}{2}.        \tag{6}
\]

Those intervals are consecutive disjoint levels of sizes `3^h`.

There is now an explicit labelled tree isomorphism to the Berggren tree:
decode a primitive Pythagorean triple by its unique Berggren parent word,
replace `L,M,R` by `-1,0,1`, evaluate (6), and send it to nH. Conversely
(5) reconstructs the word and hence the triangle. This preserves root,
directed adjacency, branch label and depth. It loses side lengths,
odd-square shell, prime factorization, geometric angles, and fruit
positivity. The retained word is the necessary sidecar.

This is more than labelling two arbitrary countable sets: (4) are fixed
algebraic maps supplied by the elliptic group law, and their unique-parent
property was proved. It is still a combinatorial tree isomorphism, not a
rational change of variables between the underlying curves. A nonsingular
conic is a genus-zero curve, and Riemann--Hurwitz forbids a nonconstant
rational map from it to E. The basic curve and group facts, including
`deg[3]=9`, are treated in [Elkies's course, Chapters II and III](https://people.math.harvard.edu/~elkies/M223.24/index.html).

For H=9G this tree is precisely `{9nG:n≥1}`. It does not cover all rational
points of E; it misses even the positive fruit solution at 17G. This last
claim is exact: 17G is positive, and equality with 9nG would imply
`17=9n`. There is no need for a full-group generator theorem to prove
this particular omission.

The canonical-height law gives `hhat(nH)=n² hhat(H)`. Thus the elliptic
canonical heights across depth h are `Theta(9^h)`, whereas logarithmic
Berggren side height is `O(h)`: a fixed matrix-norm bound gives
exponential side-length growth at most. Some boundary paths grow more
slowly. This explains large arithmetic sizes without identifying the two
notions of height.

## 4. The torsion coordinate obstructs extending the same tree to all points

On the entire rational curve, each map in (4) has collisions:

\[
 \Phi_i(P+Q)=\Phi_i(P),\qquad Q\ne O,\quad[3]Q=O.         \tag{7}
\]

It is injective only on the specified infinite cyclic set, which has no
torsion. An exact repair on a fixed torsion coset is available. For a
chosen R define

\[
 \Phi_i^R(P)=3P+iG-2R.
\]

On `{nG+R:n≥1}` it sends `nG+R` to `(3n+i)G+R`. Thus six independent
labelled trees cover the positive-free-coefficient halves of the six
cosets; their negatives cover the negative halves. Under the cited full
group description these twelve trees cover all nontorsion rational
points. The separate torsion coordinate is essential; this is not one
intrinsic tree based solely at the large triple.

## 5. Literal division by three gives seven vertices at 9G

**PE2, PROVED on `<G,T>`; CITED extension to all E(Q).** Write a point
as `(n,k)=nG+kT`, with `k` taken modulo six. Its rational preimages under
tripling, within this subgroup, are exactly

\[
 (m,j):\quad 3m=n,\qquad3j\equiv k\pmod6.                \tag{8}
\]

They exist iff `3|n` and `k=0 or3`. If they exist, there are exactly
three, differing by the order-three translation Q. Among their torsion
labels, exactly one is itself 0 or3, so at most one branch can continue
to another division. This is the guard missed by a picture with three
rational preimages at every vertex.

At the supplied seed the complete subgroup inverse tree is

```text
9G
  3G                 3G+Q          3G−Q
    G, G+Q, G−Q      no children   no children
```

There are seven vertices. The bottom row cannot divide again because
the free coefficient is one. More generally the inverse tree from
`3^rG`, for integer `r≥0`, has exactly `1+3r` vertices, with three vertices
at every nonempty level and a single continuing branch. The cited
Mordell--Weil generation makes these the complete rational inverse trees;
without that external input they remain complete within `<G,T>`.

Over an algebraic closure, [3] has nine preimages, not three. When a
rational preimage exists, the three rational choices form a coset of its
rational kernel; on fruit coordinates they are cyclic permutations.
Neither forgetting those permutations nor
forgetting the free coefficient creates an infinite backward tree.

## 6. Positivity is not hereditary under the forward tree

The inherited real positivity test for an affine curve point `(x,y)` is

\[
 x<-14/3\quad\hbox{and}\quad x^2+112x+784>0.              \tag{9}
\]

It is equivalent to the projective triple in (1) having all coordinates
of one sign. The exact census `1≤n≤40` finds this at n=9 and n=17 only.
In particular **all three children `18G,27G,36G` of H fail positivity**.
The even children already fail at the component level. Keeping only
positive vertices in this particular tree leaves its root with no
children. A different positive recurrence would need an additional
selection or return mechanism.

There is a general obstruction behind the small test.

**PE3, PROVED from the standard real-group description.** No map
`f(P)=[m]P+R`, with integer `|m|≥2` and rational R, maps every positive
rational fruit solution to another positive rational fruit solution.

Proof. The real elliptic curve has two components, each a circle; its
bounded component C contains G and all positive fruit solutions.
The positive locus S is a nonempty open proper subset of C by (9): it
contains H, whereas a whole neighbourhood of G is outside its closure.
The subgroup generated by 2G is dense in the identity circle because it
is an infinite subgroup of a circle. Thus odd multiples of G are dense
in C, and the rational points in every open subarc of S are dense there.

If f sends C into the other component, it already violates the desired
property. Otherwise, in circle coordinates on C its restriction is
`t -> mt+r (mod1)`. Suppose it preserves all rational points of S.
Choose a nonempty small arc I inside S. For every integer j≥1, the dense
set of rational points in I has its image under `f^j` in S; continuity
gives `f^j(I) subset closure(S)`. For large j, `|m|^j` times the arc length
exceeds one, so `f^j(I)=C`. This contradicts `closure(S) != C`.

The proof uses rational density only to pass from a purported rational
invariance to an image containment in the real closure. It does not
assume individual multiplication orbits are dense. It rules out global
positivity preservation by these expanding maps, not every specially
selected orbit, variable first-return rule, or state-dependent tree.
Together with the polynomial word-growth argument for maps `±P+R_i`,
it sharply limits a direct analogue of three globally positivity-preserving
Berggren branches within these affine elliptic operations.

## 7. The gap 24 gives a different elliptic curve

The identity `2^10=10^3+24` is exact. Because ten is even, it gives the
rational point `(x,y)=(10,32)` on

\[
 E_{24}:y^2=x^3+24.
\]

Its double is `(505/256,23053/4096)`, directly checked by the group law.
The gap 24 is outside the gap-one hypothesis of Catalan's theorem; the
[inherited Catalan audit](catalan_elliptic_20260921_catalan.md) explains the
scope of that theorem. Even the equation `2^k=n^3+24` has the smaller
solution `(k,n)=(5,2)`. A finite check through `k=200` finds only these
two pairs; no all-exponent classification follows.

The new curve has `j=0`, whereas the fruit curve has
`j=1408317602329/2153060`, so they are not isomorphic over an algebraic
closure. More strongly, at the common good prime 19 the point counts are
27 and 18, respectively. Isogenous curves have equal good-reduction point
counts, so they are not Q-isogenous. Both counts are independently
computed by enumeration and quadratic characters. The isogeny invariant
used here is explicitly recalled in [Elkies, Chapter VII, Example 3.3.3 commentary](https://people.math.harvard.edu/~elkies/M223.24/index.html).
This retains a legitimate elliptic point from the numerical identity,
while rejecting an identification with the fruit curve.

## 8. Reproduction and finite scope

```text
python 04-computation/experiments/prime_shells_20260921_elliptic.py
python -O 04-computation/experiments/prime_shells_20260921_elliptic.py
```

The [script](../../04-computation/experiments/prime_shells_20260921_elliptic.py)
writes the [JSON certificate](../../04-computation/experiments/prime_shells_20260921_elliptic.json)
beside itself and embeds its LF-normalized source hash. Explicit exception
checks survive optimization. The complete integer-coefficient/Berggren
comparison covers all 9,841 vertices through depth eight. Actual rational
elliptic branch operations at seed G are checked through depth three;
all three children of the positive seed H=9G are also checked. All
`nG` for `1≤n≤40` are computed by sequential and binary addition paths.
Additional controls cover the 12-point symmetry orbit, 972 inverse
coefficient cases, roots `3^rG` for `0≤r≤8`, and the two independent
finite-field counts. Decimal sizes and positive-index lists have exactly
this finite scope; the infinite tree and obstruction statements are
proved in the text.
