**PROVED** for core order at least six; **FINITE-EXACT** for smaller orders.
Companion to [the global cycle spectrum](PROOF.md).

# Sharp cycle maximum with two apex vertices over a transitive core

For a transitive core of order n≥2 and two arbitrary apex vertices, the largest
number of directed Hamilton cycles (up to cyclic rotation) is 1, 3, 5, 14 for
n=2,3,4,5 respectively. For n≥6 it is exactly

    (2^floor(n/2) − 2)(2^ceil(n/2) − 2) + 1 + [n odd].

## Counting and compression

Write each cycle uniquely as x A y B x, with A,B increasing complementary
core words (empty permitted). If one apex does not beat the first core vertex,
then that vertex is forced to the other apex's following word. If one apex
beats the last core vertex, then that vertex is forced to the word not
preceding that apex. Either constraint allows at most 2^(n−1) cycles. But
both apices beating the first floor(n/2) core vertices and losing to the rest
gives (2^floor(n/2)−2)(2^ceil(n/2)−2)+1: each of the two words must meet both
blocks unless a word is empty, with exactly one empty-word cycle. This
exceeds 2^(n−1) for n≥6, as the two factors after division by their powers
of 2 are each at least 3/4. Thus a maximizing pair both beat the first vertex
and both lose to the last.

For interior vertices i=1,…,m=n−2, write x_i,y_i for the outgoing apex bits,
and s_i=x_i+y_i. Under the required boundary orientations the exact count is

    1 + Σ_(i<j) 2^(j−i−1) s_i(2−s_j)
      + Σ_(i<m) [x_(i+1)(1−x_i)+y_(i+1)(1−y_i)]
      + Σ_i [s_i=1].

Indeed if the two extreme core vertices are in the same word, the other
word has an interior minimum i and maximum j. For i<j the in-between
vertices are freely assigned, and the required bits are y_i(1−x_j)
or x_i(1−y_j) according to which word contains the extremes. For i=j
these terms sum to [s_i=1]. If the extremes are in different words,
take the minimum i of the word containing the last core vertex and
maximum j of the word containing the first. For i<j the contributions
are x_i(1−x_j)+y_i(1−y_j), with the same free vertices. The only
other possible relation is i=j+1 (an unbroken cut), giving the adjacent
inversion sum. Counts using extreme indices as i or j give zero by the
boundary orientations. The constant 1 counts the unique empty-word cycle.

Now swap neighboring columns with s_i<s_(i+1). These two binary columns
are componentwise ordered, so each changed row switches 01 to 10, losing
its one internal inversion while only gaining or preserving external
adjacent inversions. The number of changed rows is s_(i+1)−s_i. The
corresponding pair term gains exactly twice that number; all other pair
terms together are nondecreasing, because weights are larger at larger
separation. The singleton sum is invariant. Hence every maximum has
s sorted as 2^a 1^b 0^c. In the 1-run each transition between different
apex bits 01 and 10 contributes one adjacent inversion. Alternating
there gives exactly max(0,b−1), and there are no inversions elsewhere.

## Optimize the compressed state

Put A=2^a, B=2^b, D=2^c. The formula now becomes

    b=0: 2^n − 4(A+D) + 5;
    b≥1: 2^n − (2B+2)(A+D) + B + b + 3.

For the second line the four pair contributions are
4B(A−1)(D−1), 2(A−1)(B−1), 2(B−1)(D−1), and B−b−1;
the other three terms add 1+(b−1)+b. At fixed b, maximizing means
balancing a,c as closely as possible. If b+2≤m, moving from b to b+2
with balanced a,c halves A+D. For b≥1 the change is
−(2B−1)(A+D)+3B+2<0, since A+D≥4 and B≥2. The b=0 to b=2
change is 4−(A+D)≤0. Consequently only b=0 and b=1 need comparison.
When m is even, the b=1 value is strictly smaller; when m is odd it
is exactly one larger. The b=0 value is the basic balanced cut count,
giving the claimed maximum. Every listed maximizing profile is
realizable by the specified orientations.

For core order at least six, equality therefore requires exactly the
balanced identical cuts when the core order is even, and the two cuts
differing at the single center when it is odd, up to exchanging apices.
Either interapex orientation works.

The companion checker enumerates both apex outmasks and chooses the
better empty-run term for n=2,…,7, directly using the rooted words;
it also checks the formulas for the compressed profiles.

Replay: python3 feedback-two-maximum-check.py.
