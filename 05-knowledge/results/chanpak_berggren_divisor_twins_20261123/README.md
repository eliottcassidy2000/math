# One distinguished branch, a pair of twins, and two different readouts

The owner's three patterns are not interchangeable as scalar predicates.
They do have a shared typed shape: a distinguished choice and two
choices exchanged by a sign or coordinate swap. Here are exact links
that preserve that distinction, together with what each map forgets.

## Signed doubling really is a Berggren branch pair

For Euclid spinors \(m>n>0\) of opposite parity and coprime, put
\(a=m-n,b=m+n\). These are odd, positive, coprime, \(a<b\).
The triangle in odd-leg, even-leg, hypotenuse order is
\[
 (ab,\ (b^2-a^2)/2,\ (b^2+a^2)/2).
\]
The root \((3,4,5)\) has \((a,b)=(1,3)\). In these coordinates the
three [Berggren moves](https://arxiv.org/abs/math/0010281) are
\[
 A(a,b)=(a,b+2a),\quad
 B(a,b)=(b,2b+a),\quad
 C(a,b)=(b,2b-a).
\]
These are the inherited spinor moves in new coordinates. The \(A\)
branch retains the old small factor; the other two retain the old
large factor and choose opposite perturbations around \(2b\).
Their ratio chambers \(b'/a'\) are respectively \(>3\), \((2,3)\),
and \((1,2)\). Each has a unique inverse there:
\((a',b'-2a')\), \((b'-2a',a')\), \((2a'-b',a')\).
The inverse preserves oddness and gcd and decreases the large
coordinate. The only boundary ratio 3 is the root, by coprimality;
ratio 2 is impossible for two odd factors. This proves the
three-way distinction as an actual integer grammar, not a visual
resemblance.

On the distinguished branch from the root the input pairs are
\((1,p)\) for odd \(p\ge3\). They need not be prime. The
[inherited bracket theorem](../procgen_numerology_20260923_odd_square_brackets.md)
selects just \(\{2,3,11\}\) as primes \(p\) for which \(p,2p\)
lie in the same consecutive odd-square bracket. The two odd inputs
are on this Berggren spine; 2 is not an odd factor. Equivalently,
2 already occurs in the seed as its even leg divided by two.
The exact signed readings are:

| Input | Move | Child pair | Child primitive triangle |
|---|---|---|---|
| \((1,3)\) | \(B\) | \((3,7)\) | \((21,20,29)\) |
| \((1,11)=A^4(1,3)\) | \(C\) | \((11,21)\) | \((231,160,281)\) |

Thus the user's \(2\cdot3+1=7\) and \(2\cdot11-1=21\) are genuine
branch **factors**. They are not, in general, the resulting odd
legs or the resulting tournament heights. The first odd leg is
\(3\cdot7=21\); the second is \(11\cdot21=231\).
This identifies an exact nesting of the two numbers in the
triangle tree, including the previously discussed 231.

## A literal 1+2 match on tournament heights, exactly at the seed

Use the inherited [Berggren-to-apex lift](../session256_chan_pak_triples_20261111/spinor/PROOF.md).
The base is a transitive list; O means the apex beats the vertex,
I means the vertex beats the apex. The root word OI is a directed
triangle. On the three root branches the words are
\[
 A:\mathrm{OII},\quad B:\mathrm{OIOIO},\quad C:\mathrm{OOOI}.
\]
Their Hamilton-path counts are **\(5,9,9\)**. In particular 7 is
the midpoint between the distinguished child's count and either
twin. After stripping the trailing singleton from B, the twins are
the strong five-vertex cores OIOI and OOOI. They are not isomorphic,
even though their counts and score multisets agree: their multisets
of pairs (outdegree, triangles through vertex) differ.

This is also a literal pair of arc flips from one on-label low
host. In lexicographic pair order on \(\{0,1,2,3,4\}\), a set bit
means the smaller vertex beats the larger. The lower mask 703
has count 5 and is the four-vertex A-core with one ordinal
singleton. Flip \((0,4)\) to get C's core, or flip \((0,1)\) to
get B's core, in both cases at count 9. Flipping both also gives
9. Thus the count pattern \(5,9,9,9\) lives on a real
two-coordinate square; it is not a relabeling of a missing
7-vertex in the height image.

**The equal twins separate after the seed.** For the spine input
\((1,2j+1)\), \(j\ge1\), the parent apex word is
\(\mathrm{O}\mathrm{I}^j\), with \(H=1+2^j\).
The B and C children use \(\mathrm{OI}(\mathrm{OIO})^j\) and
\(\mathrm{O}(\mathrm{OOI})^j\). The inherited transfer tracks
height \(H\), binary word value \(P\) with O=1, and last bit \(e\):
O changes \((P,e)\) to \((2P+1,1)\) without changing \(H\); I
adds \(P+e\) to \(H\) and sets \((P,e)=(2P,0)\).
An OIO block gives \(P'=8P+5,\ \Delta H=2P+2\), starting at
\((H,P)=(3,2)\); OOI gives \(P'=8P+6,\ \Delta H=4P+4\), starting
at \((1,1)\). Summing the geometric series gives
\[
 H_B=3+\frac{38(8^j-1)}{49}+\frac{4j}{7},\qquad
 H_C=1+\frac{52(8^j-1)}{49}+\frac{4j}{7},\qquad
 H_C-H_B=\frac{2(8^j-8)}7.
\]
The scalar coincidence is therefore **exactly at the root on this
spine**. At \(p=11\) the three heights are \(65,25417,34777\).
Thus the same signed factor read produces 21 but the literal
Berggren tournament lift does not produce a height near 21 here.
The distinction is the order of the apex word, not just its two
letter totals.

## Divisor balance versus Pythagorean leg multiplicity

As in the [inherited classification](../arithmetic_braids_20260917_divisors.md),
\(F\) counts nontrivial proper divisors, \(U\) proper distinct
prime divisors, and \(S\) nontrivial proper squarefree divisors.
Thus \(F-S\) is the number of proper nonsquarefree divisors.
The equality \(F=U+S\) holds for primes, prime cubes, and
\(p^2qr\) with distinct primes.

There is a useful **leg analogue** of the previous [hypotenuse
codec](../session163_ternary_triples_20261103/triple_codec/NOTE.md),
with no assumption that the odd primes are \(1\bmod4\).
A unitary divisor of \(M\) takes all or none of each prime power.
Complementary unitary divisors \(u\leftrightarrow M/u\) form
unordered coprime factor pairs.
- For any odd \(N>1\), the primitive triangles with odd leg
  \(N\) are exactly the pairs \(a<b,\ ab=N,\ \gcd(a,b)=1\),
  by the formula in the first section.
- For \(4\mid N\), use coprime factor pairs of \(M=N/2\) as
  \(\{m,n\}\); one is even, one odd. These parameterize the primitive
  triangles with even leg \(N\). If \(N=2\bmod4\), there are none.

When a positive leg \(N>1\) is admissible, this yields
\(2^{\omega(N)-1}\) triangles (halving a multiple of 4 does not
lose the prime 2). The unit pair \((1,M)\) is distinguished.
The exponent pattern is **not recoverable from the count**:
for composite admissible \(N\), the number of nonunit classes
equals \(U\) precisely when \(\omega(N)=3\), regardless of
prime exponents.

At \(N=p^2qr\) the nonsquarefree divisors including \(N\) are
\(\{p^2,p^2q,p^2r,N\}\). They are exactly the unitary divisors
containing the whole power of \(p\), a canonical transversal
for complementary pairs. For an admissible even leg with this
exponent shape, \(p=2\) necessarily, so divide these four
divisors by 2 to obtain the corresponding transversal for \(N/2\).
The proper nonsquarefree divisors now give exactly the three
nonunit triangles. The \(\{q,r\}\) swap fixes the split
\(p^2\mid qr\), and exchanges \(p^2q\mid r\) with \(p^2r\mid q\).
This is a literal **one fixed point and one two-element orbit**
on the factor-choice addresses. It is not claimed to be an
isometry of the three numerical triangles.

The prime bank itself gives a sharp test. Of
\(2^2\cdot3\cdot11,\ 2\cdot3^2\cdot11,\ 2\cdot3\cdot11^2\)
(all satisfying \(F=U+S\)), only **132** is a primitive leg,
because the others are \(2\bmod4\). Its four triangles have
odd leg and hypotenuse respectively
\((4355,4357),(1085,1093),(475,493),(85,157)\),
all with even leg 132. The corresponding factors of 66 are
\(1\mid66,\,2\mid33,\,3\mid22,\,6\mid11\).
Conversely **231** (the Berggren child above) gives four
odd-leg triangles but does not satisfy divisor balance:
it is squarefree composite, so \(F=S\) and \(U=3\).
A prime cube has only the unit-pair triangle when admissible,
whereas it has one nonsquarefree proper divisor. These cases
separate the two exponent mechanisms.

In particular every odd \(N>1\), including both 7 and 21,
is a primitive odd leg. For 7 its one pair is \(1\mid7\);
for 21 its two pairs are \(1\mid21\) and \(3\mid7\).
The \(3\mid7\) factorization survives in geometry but cannot
build a tournament by strong-component multiplication, since
7 is missing. The \(1\mid21\) geometry imposes no such
component decomposition; the strong 21-core needs the
[separate absolute exclusion](../chanpak_gap_cubical_geometry_20261122/README.md).
This is why the shared branch structure is useful without
being a proof of Chan–Pak.

## Which of the two 9-twins crosses the second hole?

[The literal two-hole square](PHYSICAL-SQUARES.md) puts
\(5,9;19,23\) on one six-vertex host, with **the same horizontal
arc** crossing both gaps. Globally there are exactly two reduced
marked squares, exchanged by converse, and either one has the C
twin OOOI at 9 after a singleton is removed. Neither has the B
twin OIOI at that corner. This is a precise instance of identical
local heights but different global continuation.

There is a single cube that shows both at once. At lower mask 4
on six vertices, use horizontal \(e=(3,4)\), vertical \(f=(0,5)\),
and extra \(g=(2,4)\). With rows indexed by \(f\) and columns by
\(e\), its two opposite faces have heights
\[
 g=0:\begin{pmatrix}5&9\\19&23\end{pmatrix},\qquad
 g=1:\begin{pmatrix}9&9\\29&23\end{pmatrix}.
\]
At \(f=e=0\), the \(g\)-flip selects the B twin. At \(f=g=0\),
the \(e\)-flip selects the C twin. Under the same \(f\)-flip,
the first goes from 9 to 29 and the second from 9 to 23.
The two signed gaps at 7 and 21 therefore do share an exact
geometric socket, but it is **typed by the five-vertex core**,
not determined by a count, a bracket, or a prime factor.

Run python3 replay.py for the factor reads, explicit triangles,
morphic formulas, both 9-cores and all eight cube heights.
