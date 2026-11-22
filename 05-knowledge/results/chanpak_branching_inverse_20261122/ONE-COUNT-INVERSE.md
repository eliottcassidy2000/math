# A branching two-port family with a literal inverse on ordinary counts

**PROVED.** There are actual two-port programs indexed by all finite binary
words, starting from one marked tournament, with distinct ordinary Hamiltonian
path counts. A single count determines the word, including its length.
The output uses only linearly many vertices in the length of the word.
The family contains at least $X^c$ counts through $X$ for some $c>0$,
but still has density zero. Chan–Pak remains **OPEN**.

The inherited input is the physical compilation in [Rooted Resets R1–R2](../two_port_return_induction_20261103/ROOTED_RESETS.md). Its two full-rank return words $041$ and $403$ are the moves below. Here the new work is to separate their complete positive images and to turn the vector distinction into a single unmarked ordinary count. No marked count or extraneous residue is part of the decoder's input.

## 1. Two physical return words with disjoint state images

Use the inherited row state $u=(h-s,q-s,s,F-h-q+s)$ on the root-child
interface. Its entries are nonnegative, the ordinary count is $u c$
with $c=(1,0,1,0)^{\mathsf T}$, and both moves return to the same type
of interface. Start at the actual entry $u_0=(2,0,1,0)$, inherited
entry word $1$. The moves are
$$
 A=\begin{pmatrix}7&1&2&1\\4&2&0&1\\4&2&1&2\\5&1&1&1\end{pmatrix},
\qquad
 B=\begin{pmatrix}6&3&1&4\\3&0&1&2\\4&0&1&2\\5&2&1&2\end{pmatrix}.
$$
Their determinants are $4,-4$. All words act on literal marked
tournaments using the inherited provisional-newborn convention; after the
last move output the actual rooted parent, not its provisional newborn.

Put $v=(-1,5,6,-6)^{\mathsf T}$. Then
$Av=(4,0,0,0)^{\mathsf T}$ and $Bv=(-9,-9,-10,-1)^{\mathsf T}$.
For any nonzero nonnegative row, the last letter is $B$ exactly when
the output paired with $v$ is negative, and is $A$ otherwise.
The chosen inverse matrix then recovers the preceding integer row.
For our seed the $A$ side is strict as well, since the first coordinate
stays positive. Thus the full marked states form an injectively labelled
binary tree at each depth, with an exact one-step inverse.

## 2. A single physical count can recover an entire bounded state

For each $n\ge0$, write $M_n=3\cdot14^n$. Every row after $n$ bits
has nonnegative integer coordinates of sum at most $M_n$, since both
matrices have maximum row sum at most 14. Append a fixed suffix of
$k$ moves $A$. Its output is the **ordinary** count $uA^k c$.
We show that this observation is injective even on all integer rows in
the nonnegative box $\|u\|_1\le M_n$, if $k\ge100+16n$.

The characteristic polynomial and the numerator column of the resolvent
of $A$ are
$$
 f(z)=z^4-11z^3+13z^2-15z+4,\qquad
 g(z)=\operatorname{adj}(zI-A)c
 =\begin{pmatrix}
 z^3-2z^2-3z+2\\
 4z^2+6z-6\\
 z^3-6z^2+19z-8\\
 6z^2-6z+4
 \end{pmatrix}.
$$
The polynomial is irreducible modulo 3: its reduction has no root,
and none of the three monic irreducible quadratics over $\mathbb F_3$
divides it. It is therefore irreducible over $\mathbb Q$ and separable.
There is a real root $\lambda\in(9,10)$ by a sign change. On $|z|=3$
the cubic term has magnitude $297>81+117+45+4$, so the other three
roots have modulus less than 3 by Rouché's theorem.

Given the difference of two rows $d\ne0$, put $D=\|d\|_1$.
The polynomial $d g(z)$ cannot vanish identically: the Krylov matrix
with columns $c,Ac,A^2c,A^3c$ has determinant $496$, and the
resolvent at infinity gives these four coefficients successively.
It follows by irreducibility that $d g(\lambda)\ne0$. Its algebraic
norm is a nonzero integer. At the other three roots $|d g(z)|\le150D$
by the displayed coefficients; hence $|d g(\lambda)|\ge(150D)^{-3}$.
This prevents a tiny nonzero dominant projection from passing as a
collision.

The spectral expansion is
$$
 d A^k c=\sum_{f(z)=0} z^k\,\frac{d g(z)}{f'(z)}.
$$
At $\lambda$, $|f'(\lambda)|<8000$; at the other three roots
$|f'(z)|<500$. The product of the four derivatives is, up to sign, the
nonzero integral discriminant, so at any one of the three small roots
$|f'(z)|^{-1}<2\cdot10^9$. The dominant term has absolute value at
least $\lambda^k/(8000(150D)^3)$; the sum of the others has absolute
value at most $9\cdot10^{11}D\,3^k$. Since $D\le6\cdot14^n$ and
$\lambda/3>3$, the first bound strictly exceeds the second for
$k\ge100+16n$. In particular $dA^k c\ne0$.
This is a uniform height comparison over the whole box, not a
finite training set of observed paths.

## 3. A self-delimiting depth code with linear vertex cost

For a data word $w$ of length $n$, let $m$ be the least power of two
with $m\ge n+1$. The branch word is $w1\,0^{m-n-1}$.
Its last $1$ is a delimiter, and $m<2(n+1)$. Append
$$
 k(m)=100+16m
$$
copies of $A$ to the resulting marked parent, and take its ordinary
count. For allowed widths $m=1,2,4,\ldots$, the possible counts lie in
the integer interval
$$
 [L_m,U_m],\qquad
 L_m=\min_i(A^{k(m)}c)_i,\qquad U_m=M_m\max_i(A^{k(m)}c)_i.
$$
These intervals are disjoint. Put $g_*=g(\lambda)$, a positive eigenvector
whose entries lie in $[1,2300]$ by the explicit polynomials and
$9<\lambda<10$. The vector $Ac=(9,4,5,6)^{\mathsf T}$ gives
$g_*/2300\le Ac\le9g_*$ coordinatewise. Multiplication by $A^{k-1}$
shows
$$
 L_m\ge \lambda^{k(m)-1}/2300,\qquad
 U_m\le 62100\cdot14^m\lambda^{k(m)-1}.
$$
Since $k(2m)-k(m)=16m$ and $9^{16m}>2300\cdot62100\cdot14^m$
for $m\ge1$, we have $L_{2m}>U_m$.

Given a count $N$, find its unique width band, if any. Search the finite
box $\|u\|_1\le M_m$ for its unique nonnegative integral solution to
$uA^{k(m)}c=N$; reject if none. Apply the sign test and exact inverse
$m$ times, rejecting unless the states stay nonnegative integral and
end at $u_0$. Remove the delimiter and following zeroes; reject if
there is no delimiter or it is inconsistent with the least-power-of-two
rule. Conversely any successful reconstruction gives an actual program.
This is exact membership and reconstruction in the absolute count image,
not just a decoder given a profile or a residue. The four-variable box
search can also be solved by [fixed-dimensional integer programming
(Lenstra)](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1983i/art.pdf):
it is one linear equation and linear inequalities in four integer
unknowns. The first width band exceeding $N$ has $m=O(\log(N+2))$
by the lower bound on $L_m$ and doubling; all numbers handled have
$O(\log(N+2))$ bits. The exact sign-and-inverse steps and the literal
compiler then run in polynomial bit time. Thus the restricted image is
also polynomial-time recognizable and decodable, without evaluating
Hamiltonian paths by a general counting oracle.

There are exactly $2^m-2^{m/2}$ admissible branch words at each power
$m\ge2$, corresponding to delimiter positions past the midpoint.
The construction uses $3+3(m+k(m))=O(n+1)$ vertices. A coarse upper
bound on any width-$m$ count is $3\cdot11^{100}(14\cdot11^{16})^m$
using the maximum row sum of $A$. Choosing the largest power of two
below a suitable constant multiple of $\log X$ therefore gives a
lower bound $X^c$ on the number of represented counts through $X$,
for some absolute $c>0$ and all large $X$. On the other hand,
each $A$ or $B$ raises $h$ by at least a factor five, since
$Ac=(9,4,5,6)^{\mathsf T}$ and $Bc=(7,4,5,6)^{\mathsf T}$.
A represented count through $X$ thus has $5^{16m}\le X$.
Even summing over all integral widths gives at most
$O(X^{\log 2/(16\log 5)})$ represented counts. This beats the
polylogarithmic fixed-clock barrier by recurrent branching, but cannot
give dense coverage.

## Transfer of ideas from the linked papers

[Pasten–Vidaux, §§4–7](https://arxiv.org/pdf/2611.24748) first get local
approximation, then use valuation parity and a bound on a nonzero
conjugate difference to turn local information into an exact restriction.
Our slope test is exact already on the marked vector, but a scalar
observation could hide a small difference. The norm bound here plays the
conjugate-difference role: remaining embeddings bound how close two
distinct integral states can be in the dominant one, and the suffix makes
the measurement error smaller than that gap. No elliptic-curve or abc
hypothesis is used. The accompanying odd-adic lane is the direct
analogue of the local approximation step.

[The linked announcement, pp. 1, 3–4](https://alpo.ge/remark.pdf)
describes accumulating multiplicity by iterating lifts and accounting for
its cost, with sign or positivity used to retain information.
The literal positive return matrix is what is iterated here, not a
theta construction. A sign separates the two physical images; a
dominant mode gives the amplification; and counting the padding gives
an explicit cost in vertices and in output density. This uses the
announced strategy as a source of questions, not its asserted RH or GRH
conclusions as a proof input. The next coverage problem is to replace the
costly common suffix by a cheap, regular scalar discriminator while
broadening the set of physically admissible inverse moves.

Reproduction of the new algebraic constants and short literal controls:
`python3 check_twoport_inverse.py`.
