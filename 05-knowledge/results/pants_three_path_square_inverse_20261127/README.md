# Three-path reconstruction and the first hidden Hamiltonian counts

2026-11-27. The [three-run switching response](../balanced_pants_switch_fourier_20261127/README.md)
recovers the labelled square of any skew tournament matrix. Over
the rationals, simple spectral planes leave an exact sign problem.
The same square also determines all four- and five-vertex path
responses, but can already hide a Hamiltonian count at order six.

## 1. Recovering the square and its rational roots

Write \(A_{ij}=+1\) for \(i\to j\), \(A_{ji}=-A_{ij}\) and
\(A_{ii}=0\). For a cut sign vector \(e\), put \(D_e=\mathrm{diag}(e)\)
and let \(P_k(e)\) count directed paths on \(k\) distinct vertices
in \(D_e A D_e\). With \(G=A^2\),
\[
 P_3(e)=\frac{(n)_3}{4}+\frac12\sum_{i<j}G_{ij}e_i e_j.
\]
The second flip difference at coordinates \(i,j\) is \(2G_{ij}\).
The diagonal is \(-(n-1)\). Thus the zero-, one- and two-flip
queries, at most \(1+n+\binom n2\) in all, recover the labelled
square and the full response. This is stronger data than the
eigenvalues alone.

Write \(\det(tI-A)=t^\epsilon h(t^2)\), where \(\epsilon=n\bmod2\).
Even principal subtournaments have nonsingular signed adjacency:
modulo two their matrix is \(I+J\), which is invertible in even
size. Thus the kernel of \(A\) has dimension \(\epsilon\).
If \(h\) is squarefree, the nonzero eigenspaces of \(G\) are real
two-dimensional planes. Every skew \(B\) with \(B^2=G\) commutes
with \(G\); on each plane it is a scalar multiple of \(A\), and
it vanishes on the kernel. Interpolation gives
\(B=A p(G)\) with a unique polynomial \(p\) of degree less than
\(\deg h\). If \(B\) is rational, so is \(p\): the independent
rational matrices \(A,AG,\ldots\) give a consistent rational
linear system for the coefficients. The condition is
\(p^2=1\bmod h\). Hence \(k\) distinct rational irreducible
factors give exactly \(2^k\) rational skew roots by the Chinese
remainder theorem. In particular irreducible \(h\) leaves only
\(A\) and its converse. For reducible squarefree \(h\), a branch
is a tournament exactly when its entries are all integral with
odd off-diagonal entries. The diagonal equation
\(\sum_{j\ne i}B_{ij}^2=n-1\) then forces every odd entry to be
\(\pm1\).

**Decoder from the counts alone.** Recover \(h\) from
\(\det(tI-G)=t^\epsilon h(t)^2\). Solve over \(\mathbb Q\) for
a nonzero skew matrix \(X\) commuting with \(G\), taking it
nonzero on every spectral factor. Such an \(X\) exists (there
is a true tournament among the roots); with a rational basis,
avoid the finitely many proper subspaces where a restriction
vanishes. If \(h\) is irreducible, any nonzero \(X\) works.
Express \(X^2=G R(G)\) on the nonzero eigenspaces by linear
algebra. In the algebra \(\mathbb Q[t]/h\), solve \(z^2=R\)
by factoring \(z^2-R\) over each factor field and assembling
the roots. They are invertible by the choice of \(X\), and
\(X z(G)^{-1}\), with the inverse evaluated modulo \(h\),
gives the candidates. The kernel is killed by \(X\).

For a transitive tournament of odd prime order \(p\),
the coefficients of the reduced polynomial are
\(\binom p{2j}\): each even principal minor equals one
(the transitive Pfaffian is one by its alternating recursion).
This is Eisenstein at \(p\). A nontransitive-switching
example takes the increasing orientation on \(0,\ldots,5\)
and reverses \(1\to4\). Its polynomial
\(h=t^3+15t^2+47t+1\) is irreducible by the rational root
test. The Pfaffian on \(1,3,4,5\) has magnitude three,
so this tournament cannot be switched to a transitive one.

## 2. What the square can hide: an exact five-vertex barrier

Put \(G_e=D_e G D_e\). The next two raw responses are
\[
\begin{aligned}
 P_4(e)&=\frac{(n)_4}{8}+
 \frac{n-3}{2}\sum_{i<j}(G_e)_{ij},\\
 P_5(e)&=\frac{(n)_5}{16}+
 \frac18\sum_{i<j}\bigl(G_e^2+3(n^2-6n+10)G_e\bigr)_{ij}.
\end{aligned}
\]
To prove them, expand the factors \((1+A_{ab})/2\) along
each ordered vertex sequence, working first at \(e=\mathbf1\).
A chosen contiguous run of an odd number of edges cancels by
reversing that run in the sequence. There is no room for two
even-length runs before six vertices. Only a single run of
length two, or for \(P_5\) length four, remains. For length
two the signed walk sum with endpoints \(i,j\) equals
\(G_{ij}\); the unused positions give factors \(n-3\) or
\((n-3)(n-4)\). The run can be placed in two ways for \(P_4\)
and three ways for \(P_5\); the other factor two comes from
both orders of the endpoints.

For the signed length-four walk \(i,a,b,c,j\) with fixed
distinct endpoints, the unrestricted sum is \((G^2)_{ij}\).
The nonadjacent collisions to exclude are \(i=b,i=c,a=c,
a=j,b=j\); adjacent collisions already contribute zero.
The single-collision sums are \(-(n-1)G_{ij}\) for
\(i=b,a=c,b=j\) and zero otherwise (the diagonal of
\(A^3\) is zero). Among pair intersections the nonzero
multipliers of \(G_{ij}\) are, respectively,
\(-1,-1,+1,-1,-1\) at the pairs
\[
 (i=b,a=c),\ (i=b,a=j),\ (i=c,a=j),\
 (i=c,b=j),\ (a=c,b=j).
\]
No triple collision avoids adjacent repeats or \(i=j\).
Thus the injective sum is
\((G^2+3(n-2)G)_{ij}\). Adding the length-two placements
yields the stated coefficient. The same calculation after
switching gives the two formulas for every cut.

This yields an exact obstruction to recovering Hamiltonian
path counts from the shorter responses. Begin with all
increasing arrows on zero-based labels. At order six reverse
the indicated pairs independently in each tournament:

| Tournament | Reversed pairs | \(P_3,P_4,P_5,P_6\) at the empty cut |
|---|---|---|
| \(A\) | 01, 03, 14, 23, 24, 34, 35 | 28, 39, 34, 17 |
| \(B\) | 01, 04, 12, 14, 15, 24, 34, 35, 45 | 28, 39, 34, 13 |

The full labelled squares agree by multiplication. Hence the
three-, four- and five-vertex responses agree on **every** cut,
not only the displayed one, while the Hamiltonian path counts
differ. No pair with the same labelled square can realize such
a discrepancy at order at most five, by the formulas (and
the trivial cases \(n\le2\)).

The phenomenon survives simple spectrum. At order seven,
use the following reverse lists:
\[
\begin{aligned}
 A:\;&02,03,04,13,14,15,24,25,26,35,36,45,46,\\
 B:\;&02,03,04,12,25,26,34,35,36,46.
\end{aligned}
\]
The common square has \(h=(t+7)(t^2+14t+17)\) with distinct
roots. The four rational sign branches are exactly
\(A,-A,B,-B\); thus all are tournament matrices. At the
displayed cut the path sequences from 3 through 7 vertices
are \(57,123,199,211,111\) and \(57,123,199,211,109\).
Only the counts through five are asserted to agree on
every cut.

## 3. Arbitrarily large exact same-square families with simple spectrum

Take distinct odd prime block sizes \(p_i\). In each block
switch a transitive tournament by alternating signs. The
alternating prefix sums make its rows sum to zero. Orient
all interblock arcs from the earlier to the later block.
On zero-sum vectors in a given block the matrix acts only
by that internal block. On the space of block constants,
in the orthonormal block basis it acts by the transitive
matrix weighted by \(\sqrt{p_i p_j}\). These operators
are pairwise orthogonal, with disjoint arrow supports.
Reversing any block or all interblock arrows independently
therefore preserves the labelled square and gives a
tournament. In this particular family the unswitched
Hamiltonian count is the product of the internal counts,
and does not distinguish these branches.

For an odd prime number \(k\) of blocks choose
distinct \(p_i>k+1\) congruent to 1 modulo \(k^2\). Such
primes are infinite: take \(u\) a large multiple of \(k\)
and any proposed finite list; a prime divisor of
\(\Phi_{k^2}(u)=1+u^k+\cdots+u^{k(k-1)}\) is new and the
order of \(u\) modulo it is \(k^2\) (order dividing \(k\)
would make the cyclotomic value equal \(k\) there).
The weighted quotient has reduced coefficients
\(e_{2j}(p_1,\ldots,p_k)\), by the even principal minors.
They are \(\binom{k}{2j}\) modulo \(k\); the last coefficient
is \(k\) modulo \(k^2\). The quotient is Eisenstein at \(k\).
Each internal factor is Eisenstein at \(p_i\), with pairwise
distinct degrees larger than that of the quotient.
Thus all \(2^{k+1}\) rational branches are realized and
there are no more.

The concrete sizes \(3,5,7\) give 16 labelled mates at
order 15. The four factors are
\((t+3),(t^2+10t+5),(t^3+21t^2+35t+7),(t+71)\).
They are distinct irreducibles; the quotient in this
case needs no prime progression.

## 4. Replay and small controls

Run `node replay.js`. It groups squares up to simultaneous
diagonal signs using a deterministic spanning forest of
nonzero off-diagonal entries. For \(a\) labelled
source-normalized representatives in a group, the number
of exact-square mates is \(a\,2^{c-1}\) when the support of
the square has \(c\) components: stabilizing diagonal signs
are constant on each component, modulo a global sign.
The census through order seven is recorded in
[replay.out](replay.out) and checked against all raw labelled
tournaments through order six. The replay also checks the
two unequal Hamiltonian pairs on all cuts through five
vertices, and the 16 exact block squares and sample path
queries for the order-15 example.
