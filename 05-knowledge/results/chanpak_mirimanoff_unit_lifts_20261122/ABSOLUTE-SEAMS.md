# Exact shared-apex seams at further Euclid exponents

Let \(E_p=2^{p-1}(2^p-1)\). These are literal one-arc reversals between
strong tournaments, not just separate realizations of adjacent counts.
For composite \(p\), call \(E_p\) an integer of Euclid's form rather
than a perfect number.

Use outgoing masks for an ordinal deletion \(D\) in bit order with
vertex \(0\) at the low bit. Attach the common apex beating exactly
the vertices in the tag mask. Before \(D\), an optional transitive
left flank uses the standard binary word of \(x\) as apex contacts.
After \(D\), the right flank uses the reverse complement of the
standard binary word of \(y\). An address zero is an empty flank.

| \(p\) | Lower deletion \(D\) | Tag mask | Reversed pair | \(x\) | \(y\) | Order |
|---:|---|---:|---|---:|---:|---:|
| 8 | 2,12,25,129,43,79,31,119 | 51 | 3,5 | 0 | 20 | 14 |
| 9 | 6,20,40,19,37,139,63,95 | 106 | 5,6 | 0 | 74 | 16 |
| 10 | 14,84,56,146,97,75,141,55 | 79 | 3,6 | 0 | 119 | 16 |
| 12 | 4,17,10,19,37,207,31,95 | 4 | 5,6 | 0 | 44964 | 25 |
| 14 | 68,81,74,115,101,71,0 | 65 | 3,5 | 0 | 1191170 | 29 |
| 16 | 4,17,34,7,13,27,63,127 | 86 | 2,3 | 294 | 5699 | 31 |
| 20 | 34,36,49,39,43,0 | 19 | 0,2 | 0 | 6510943581 | 40 |

The two counts in every row are exactly \(E_p-1,E_p+1\). In particular
the rows for \(p=8,9,10\) are independently counted by
literal subset-path dynamic programming on both full tournaments.
For the larger cases the replay computes the exact core profiles by
subset DP and applies the shared-apex ordinal transfer. It also builds
the full graphs, checks strong connectivity and asserts literal
one-arc equality for every row. There is no minimum-order assertion.

## Eliminating the binary read without assuming persistent gaps

Let \(h,c,z,k\) respectively be the deletion's ordinary path count,
unsigned weighted start read, signed weighted start read, and the
count with apex but without flanks. For an empty left flank and fixed
right parity \(e\), write \(K=k-ez\). The exact right-address read is
\(K+cy+hP(y)\), where
\[
 P(0)=0,\quad P(2v+1)=P(v),\quad P(2v)=P(v)+v+[v\text{ odd}]
 \quad(v>0\text{ in the even rule}).
\]
For lower and upper profiles \(h,c,K\) and \(h',c',K'\), set \(N=E_p-1\),
\(R=N-K\), and \(R'=N+2-K'\). The two exact equations for \(y,P(y)\)
have determinant \(ch'-c'h\). In particular, if \(hc'-h'c\ne0\),
\[
 y=\frac{hR'-h'R}{hc'-h'c}.
\]
An integral nonnegative \(y\) of parity \(e\) must still pass the
original two exact reads. If the denominator is zero, unequal
numerator rules out the parity chamber. Otherwise the two equations
are dependent; test the remaining one-variable read in the interval
\(\lceil R/(c+h)\rceil\le y\le\lfloor R/c\rfloor\) for \(c>0\),
using \(0\le P(y)\le y\) and the necessary congruence
\(cy\equiv R\pmod h\). This exact inverse applies even when \(h\) or
\(c\) change under reversal. Thus it is an absolute test for a
specific target rather than a phase correction or a restriction to
one persistent family.

For \(p=8\) the two profiles are
\((h,c,z,k)=(193,1424,-134,685)\) and
\((207,1412,-148,675)\). In the even chamber the determinant
formula gives \(y=20\), and then \(P(20)=18\) verifies both counts.
For \(p=10\) the profiles are \((599,4340,-396,2127)\) and
\((597,4340,-390,2151)\); the odd chamber gives \(y=119\).
The \(p=12,14\) examples have unchanged \(h,c\), but need only the
appropriate parity. The \(p=16\) certificate uses both ordinal
flanks; no scalar change of variable is substituted for them.

## A restriction from positive receiver transmission

For a core reversal with the same \(h,c\) on both sides, let
\(\Delta k=2a,\Delta z=2b\). Append a nonempty marked receiver \(V\)
before the common right flank. Writing \(H=H(V)\) and \(G\) for its
number of full paths starting tagged, the new half-gap is
\[
                         (a-b)H+bG.
\]
This is the shared-apex transfer since the number starting untagged
is \(H-G\). In particular if \(a>b\ge0\), a positive unit gap forces
\(a-b=1\) and \(H=1\); if \(b>0\), it also forces \(G=0\). Thus
these parameters cannot use a growing nontransitive receiver. The
balanced case \(a=b>0\) forces \(a=b=1,G=1\), recovering the positive
unit mechanism. For \(0<a<b\), the necessary arithmetic is
\(bG=(b-a)H+1\), so \(H\) lies in one residue class modulo \(b\) if
and only if \(\gcd(a,b)=1\).

For a receiver that is an ordinal sum \(V_1\to V_2\) of nonempty
parts, both \(H\) and \(G\) factor by \(H(V_2)\), because every full
path starts in \(V_1\). Hence an absolute unit gap forces
\(H(V_2)=1\). In particular, whatever the two half-differences
\(a,b\), all strongly connected components of a unit receiver after
its first must be singletons. The remaining room to develop an
all-exponent construction is in the first component, in an
unbalanced exact start read, or outside this ordinal architecture.

No uniform straddle is asserted.

## A variable-length receiver with closed exact profiles

On vertices \(1,\ldots,n\), let \(U_n\) use the ordinal orientation
from smaller to larger except that every consecutive pair is reversed.
Mark only \(n\). The unique full path starting there is the decreasing
sequence: the newest vertex beats only its immediate predecessor,
and deleting the newest vertex gives \(U_{n-1}\). Starting from the
cycle \(U_3\), each new vertex has a win and a loss against the
preceding strong component, so \(U_n\) is strong for every \(n\ge3\).

For any one-tag receiver \(V\), attach an apex beating only its tag
and give the new apex the sole tag. On the marked core reads the
operation yields
\[
 h_{\rm new}=k,\quad c_{\rm new}=h+c,\quad z_{\rm new}=-h-z,\quad
 k_{\rm new}=k+c-z .
\]
In the start reads, the new singleton contributes \(h\); any longer
selected path continues at the old tag, and gains one in its length.
For a singleton tag the even-path apex identity is \(k=h+c+z\), as
a nonempty even path starting at the sole tag cannot also end there.
These observations give the displayed transformation. For \(U_n\)
it starts at \((h,c,z,k)=(1,1,-1,1)\) and gives the tribonacci
characteristic polynomial \(x^3-x^2-x-1\) on \((h,c,z)\).
In particular with \(A(w)=\sum_{n\ge1} H(U_n)w^n\),
\[
 A(w)=\frac{w+w^3}{1-w-w^2-w^3}.
\]

The other reads also close. Any induced subset decomposes into
ordinally ordered maximal intervals of consecutive original labels,
because all inter-interval edges go forward. Color the \(n\) labels
by the two weak ordered path-cover slots. Each run of one color
contributes the ordinary count \(H(U_j)\) for its length \(j\);
colors alternate between runs. Thus the ordered two-cover mass
has generating function
\[
 Q(w)=1+2\sum_{r\ge1}A(w)^r
     =\frac{1-w^2}{1-2w-w^2-2w^3}.
\]
For the weighted terminal read \(d_n\), require the selected slot
to end at the tagged vertex \(n\), leaving the other slot
unrestricted. If \(g_j\) counts full paths in \(U_j\) ending at \(j\),
then \(g_1=1\) and \(g_{j+1}=h_j-g_j\): the last preceding vertex
can be any except the old tag. The final run has weight \(g_j\),
earlier runs have weight \(h_j\), and the final color is fixed.
Consequently
\[
 \sum_{n\ge1} g_n w^n=\frac{w(1-w)}{1-w-w^2-w^3},
 \qquad d(w)=\frac{w(1-w)}{1-2w-w^2-2w^3}.
\]
The unsigned and signed untagged-end reads are
\(t=Q-h-d\) and \(u=1-h\). All displayed reads can therefore be
computed with linearly many integer arithmetic operations in \(n\).
In particular any fixed balanced core with \(\Delta k=\Delta z=2\)
gives a literal unit seam for every \(U_n\) and every common suffix.

For the fixed published six-vertex core with lower
\((h_D,c_D,z_D,k_D)=(25,28,-10,43)\), the transfer even simplifies to
\[
 L_n(y)=28[Q_n(y+1)-d_n]+25[k_n+c_n y+h_nP(y)-z_n(y\bmod2)]-10
\]
and the single-arc companion is \(L_n(y)+2\). The standalone \(h_n\)
term cancels by \(43=25+28-10\), the singleton-tag identity in
the initial core. Thus both the receiver length and the binary
suffix are freely and independently selectable.

This receiver family by itself has an exact scale obstruction.
Let \(\lambda>2\) be the positive root of \(x^3-2x^2-x-2\);
it is the unique dominant root by the nonnegative companion
matrix, or by the strict triangle inequality in its root
equation. It is simple, and the generating functions give
\(Q_n\sim C\lambda^n\) with \(C>0\), as well as
\(d_n/Q_n\to1/(1+\lambda)\). The other reads \(h,c,z,k\) are
\(O(\rho^n)\) with \(\rho<2\), the positive tribonacci root.
For a fixed balanced marked core, \(c_D>0\).
The ordinary output with receiver \(U_n\) and suffix \(y\ge0\)
is therefore asymptotically
\[
 c_DQ_n\left(y+\frac{\lambda}{1+\lambda}\right)
\]
as \(n\to\infty\), uniformly in \(y\): all other transfer terms
are \(O(\rho^n(y+1))\), and the error in \(t_n/Q_n\) tends to zero.
If \(p=rm+s\), \(n=am+b\) with positive integers \(r,a\), and \(y\)
is a fixed eventually nonnegative dyadic polynomial evaluated at
\(2^m\), its leading height base would have to satisfy
\(\lambda^a2^j=2^{2r}\) for some integer \(j\ge0\).
The irreducible cubic for \(\lambda\) has norm \(2\)
(the rational root test suffices for irreducibility). Taking norms
would give \(a=3(2r-j)\), forcing
\(\lambda=2^{1/3}\), a contradiction. Thus the recursive unit
receiver provides genuine non-dyadic flexibility, but no exact
perfect-number identity on these fixed affine/dyadic width
schedules. Arbitrary target-selected widths and contact words
remain unrestricted by this argument.

## An ordinal factorization cannot hide a larger cross-component seam

Suppose an edge \(a\to b\) goes between two different strongly
connected components \(A,B\) of the lower tournament, in that order.
A reversal \(b\to a\) introduces directed cycles but deletes none;
the old cross edge lay on no directed cycle. The
[odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
therefore says that the increase is twice the sum, over new odd
cycles \(C\) containing \(b\to a\), of \(H(T-V(C))\).
The cycle complement has the same edges in both orientations.
In particular a gap of two forces exactly one new odd cycle
and a transitive complement.

The number of new triangles is
\[
 \deg_A^+(a)+|M|+\deg_B^-(b),
\]
where \(M\) contains all vertices in intermediate components.
Every new directed cycle needs a path from \(a\) to \(b\) after
omitting the new arc, so positive increase implies this number
is positive: if it were zero, \(a\) would be a sink in \(A\),
\(b\) a source in \(B\), and no intervening component could
provide a route. Thus a gap of two requires precisely one triangle.
If it uses the sole intermediate vertex, \(a\) is a sink and \(b\)
a source, so the transitive complement forces both endpoint
components and all exterior components to be transitive.

Otherwise, up to reversing the entire orientation, \(M\) is
empty, \(a\) has a sole outneighbor \(x\) in \(A\), and \(b\)
is a source in \(B\). The new triangle is \((b,a,x)\).
The complement forces \(B-b\) and \(R=A-\{a,x\}\) transitive.
Any new cycle entering \(B\) must return directly to \(b\),
since \(b\) is a source there. If \(x\) beat any vertex of the
linear order on \(R\) other than its final vertex, following
that vertex and its next one would give a second new odd cycle
of length five. Thus \(x\) beats at most the final vertex of
\(R\). All earlier vertices of \(R\) dominate the rest of \(A\).
Its last nontrivial candidate is the triple consisting of \(a,x\)
and the final vertex, which is transitive or a directed triangle.
Hence the lower count is only 1 or 3. This rules out realizing
the large perfect seam by taking the factorization
\(E_p-1=(2^p+1)(2^{p-1}-1)\) and reversing a cross arc between
the factors.

More generally every literal gap-two seam with lower count \(>3\)
has the flipped edge inside one common strong component on both
sides and only singleton components outside it. Indeed on its
lower side the flip must be internal; every untouched component
contributes an odd factor dividing the difference two, hence has
count one and is a singleton. If the changed component split after
reversal, the changed edge would cross two of its new components
(otherwise reversing it back cannot reconnect them). The inverse
flip would then decrease the count by two, impossible for a
cross-component reversal, which only adds cycles in the formula.
