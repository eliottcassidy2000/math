# The three balanced cuts, Pfaffians, and switch averages

2026-11-27. This is a second construction for the
[four-vertex pair-of-pants question](../pair_of_pants_four_vertex_tilings_20261127/README.md).
The earlier model acted on three free bits of a marked Hamilton
path. This one acts directly on the labelled tournament and
leads to sharp Fourier support bounds in every order.

## 1. A four-holed sphere and the three balanced cuts

Glue two pairs of pants along a cuff to make a four-holed sphere.
The three partitions of the four labelled boundary components
into two pairs index three marked cuts; this is not a classification
of all embedded or isotopy classes of separating curves. On the
four tournament vertices, the corresponding move reverses every
arrow across that cut, the usual [Babai–Cameron switching operation](https://doi.org/10.37236/1516).
The cuts \(01|23,02|13,03|12\) are the three nonidentity elements
of a Klein four group: adding sides by symmetric difference,
modulo taking complements, gives the group law. This geometric
marking does not identify a vertex permutation with arrow switching.

Write \(A_{ij}=1\) if \(i\to j\), and \(-1\) otherwise off the
diagonal. Every switch conjugates the skew matrix \(A\) by a
diagonal sign matrix. A balanced switch preserves each signed
matching term \(A_{01}A_{23},-A_{02}A_{13},A_{03}A_{12}\)
of the Pfaffian, since every term uses every vertex exactly once
and an even switch multiplies it by \(+1\). Each sign triple
leaves eight labelled states because its terms constrain three
disjoint pairs of edges. These eight states split into two
four-element balanced orbits, interchanged by converse.
Converse keeps the matching terms, but negates triangle products,
which cannot change under a cut switch.

| Pfaffian magnitude | Hamilton path counts on the balanced orbit | Labelled orbits |
| --- | --- | ---: |
| 1 | \(1,1,5,5\) | 12 |
| 3 | \(3,3,3,3\) | 4 |

The four-vertex classification proves the link to the Pfaffian:
transitive \(H=1\), absolute Pfaffian one; a triangle with a
source or sink \(H=3\), absolute Pfaffian three; and the strong
class \(H=5\), absolute Pfaffian one (order a directed Hamilton
cycle \(0\to1\to2\to3\to0\), so the first and third terms
cancel). The average is three by Section 2. In a mixed orbit
exactly one nonidentity move is neutral throughout. Read either
transitive state's linear order and split its first two from
its last two. That switch exchanges two adjacent transitive
blocks, pairing the two path-count-one states, and pairs the
remaining two states as well. The other cuts exchange path
counts one and five. The replay exhausts the 64 labelled states.

## 2. Every prescribed path-cover profile has the same switch mean

For \(n\) vertices let \(G=\{\pm1\}^n/\{\pm(1,\ldots,1)\}\)
denote the switch group. Write \(T^\varepsilon\) for a switched
tournament. Fix an ordered composition
\(\ell=(\ell_1,\ldots,\ell_j)\) of \(n\) into positive
path sizes, and let \(F_\ell(\varepsilon)\) count ordered
covers of \(T^\varepsilon\) by directed paths of those sizes.
The exact mean is \(n!/2^{n-j}\), for every \(T\).
Indeed, each permutation and its prescribed cuts specify a
forest of \(j\) paths. Making every requested edge point in
the prescribed direction imposes \(n-j\) independent relative
sign equations, leaving \(2^j\) solutions before quotienting
out the global sign. Summing over profiles gives
\[
\operatorname{mean}_{G} A_j(T^\varepsilon)
   =\frac{n!\binom{n-1}{j-1}}{2^{n-j}}
\]
for ordered \(j\)-path covers. For even \(n\), the same mean
holds on **each** of the two parity cosets, for each profile
separately. It follows in particular that the sum of \(H\)
over either parity coset is \(n!/2\). All means weight labelled
switch states, not unlabelled classes.

### A sharp support bound, and exactly which subgroup means are universal

Expand the path-cover indicator for each permutation as a
product of \( (1+A_{uv}\varepsilon_u\varepsilon_v)/2\) over
required arcs. A choice of expansion edges gives a subforest of
the prescribed runs. If a selected component has an even number
of vertices, reverse that entire contiguous component inside
its run. The selected edge set on labels is unchanged and
the sign product negates, because the component has an odd
number of edges. Boundary edges of the component were not
selected. Choose the leftmost even component for a canonical
sign-reversing involution. Hence only odd components survive.

The switch character is the product of \(\varepsilon_v\)
over the vertices of odd degree in the selected forest.
Every nonsingleton odd component contributes its two
endpoints and uses at least three vertices. Consequently the
only possible nonconstant Walsh characters have even weight
at most
\[
D(\ell)=2\sum_i\left\lfloor\frac{\ell_i}{3}\right\rfloor.
\]
For Hamilton paths this is \(2\lfloor n/3\rfloor\). For even
\(n\) the full parity character has weight \(n>D(\ell)\),
which proves the parity claim.

This uniform bound is sharp character by character: for an even
set \(S\) with \(2\le |S|\le D(\ell)\), place \(|S|/2\)
disjoint three-vertex components into the prescribed runs
and assign exactly \(S\) to their endpoints. The other vertices
are singletons. In the coefficient of the character \(S\),
view the independent signs \(A_{uv}\) (\(u<v\)) as formal
variables subject to Boolean evaluation. The resulting edge
monomial appears with nonzero coefficient: every surviving
copy of that same monomial is merely an arrangement of the
same selected components or a reversal within them, and every
component has an even number of edges. Thus there is no
cancellation of that monomial. Distinct square-free sign
monomials are independent functions on the edge-sign cube.
So some tournament has nonzero coefficient at \(S\).

It follows that a subgroup \(K\le G\) has the random-tournament
mean on **every coset for every starting tournament** exactly
when its annihilator contains no nonzero even set of weight
\(\le D(\ell)\). Sufficiency is character orthogonality. For
necessity choose a tournament for a surviving annihilator
character; the resulting nonzero coefficient on \(G/K\)
prevents the coset means from being constant. This is a binary
code criterion, not a claim that the counts themselves are
constant on cosets.

For the profile \((3,1,\ldots,1)\), the coefficient at the
character \(\{u,v\}\) is exactly
\((n-3)!(A^2)_{uv}/2\): choose the middle vertex of the single
three-vertex run, count both endpoint orders, and permute
the singleton runs. There are no higher characters for this
profile. This gives a direct Gram-matrix measurement via the
switch variance, which vanishes precisely when
\(A^2=-(n-1)I\). For \(n=4\), the Hamilton-path profile itself
has coefficient \((A^2)_{uv}/2\) as well.
The mechanism complements the prior
[Chan–Pak ordered-block work](https://github.com/eliottcassidy2000/math/commit/e56291fd95d0c59ba90c24b64732d8558f78cec1):
individual covers carry more information than their common
switch mean, but their switch response has strictly limited
Fourier degree.

## 3. Closing the path brings in a signed cycle invariant

Let \(C(T^\varepsilon)\) be the number of directed Hamilton cycles
modulo starting position, for \(n\ge3\). Define
\[
S_A=\sum_{\pi\in S_n}
       A_{\pi_1\pi_2}A_{\pi_2\pi_3}\cdots A_{\pi_n\pi_1}.
\]
Then the full-group mean, and the mean on either parity coset
when \(n\) is even, is exactly
\[
\frac{n!+S_A}{n\,2^n}.
\]
To see this, expand the cyclic indicator, then average signs.
The empty edge selection and the whole cycle have empty
boundary; no other selection does. In even parity averaging
the only other eligible boundary is the full vertex set.
Inside one cyclic order that arises from the two alternating
perfect matchings. The signed product for either specified
matching sums to zero over \(\pi\): swap the vertices in one
of its matched positional pairs. Hence the full boundary
contribution cancels on both parity cosets. For odd \(n\),
reversal negates each term in \(S_A\), so the cycle mean is
universal. For even \(n\), \(S_A\) is a switch invariant
rather than a universal constant. At four vertices the
balanced mean is zero in the Pfaffian-three orbit and \(1/2\)
in the mixed orbit, even though both have mean path count three.

The companion Node.js replay checks the balanced orbits
exhaustively and the Fourier and cycle mean formulas on
30 small cases. The general claims follow from the proofs
above, not from the sample.
