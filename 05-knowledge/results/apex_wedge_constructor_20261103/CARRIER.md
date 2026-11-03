# A minimal linear recurrence for the edge jump

**PROVED scoped.** Let J(w,u)=D(w,u)/2 for the two-apex sources of
[THEORY.md](THEORY.md). The operation language here appends a common
base column (w_i,u_i) and queries J. Arbitrary native moves, re-marking,
and attachment to a different physical port are not included in this quotient.

## C1 Six integral coordinates close under every appended column

Write Q_w,Q_u for the ordinary binary integers represented by w,u.
The empty word has both codes zero. Retain
\[
 z=(1,M,N,a,b,A)^T,\qquad
 M=\lfloor Q_w/2\rfloor,\quad N=\lfloor Q_u/2\rfloor,
\]
where a,b are the last bits, set to zero for the empty word, and
\[
 A=\sum_{i+2\le j}2^{j-i-2}(w_i u_j-u_i w_j).
\]
Appending the bits (c,d) acts by the integer recurrence
\[
 (1,M,N,a,b,A)\longmapsto
 (1,\,2M+a,\,2N+b,\,c,\,d,\,A+dM-cN).                    \tag{C1}
\]
The observation is
\[
 J=A+N+b-M-a.                                            \tag{C2}
\]
The weighted area term A receives exactly the interactions between the new
column and all columns except its immediate predecessor. G1 of
[GEOMETRY.md](GEOMETRY.md) supplies C2, since the linear contribution is
ceil(Q_u/2)-ceil(Q_w/2). This proves the recurrence at every word length.

Each of the four input letters has a 6 by 6 integer matrix. Thus every
fixed repeated block has an exact linear recurrence and rational generating
function. This is a recurrence for the edge derivative, smaller than the
inherited 17-dimensional linear representation of the absolute H value.
It does not replace that representation when absolute counts are requested.

## C2 Six is the sharp linear dimension

Encode a column by c+2d in {0,1,2,3}. Take these prefixes and suffixes,
including the empty word:
\[
 P=(\epsilon,2,1,20,10,201),\qquad
 S=(\epsilon,0,2,00,02,200).
\]
The matrix with row s, column p equal to J(ps) is
\[
 \begin{pmatrix}
 0&1&-1&1&-1&0\\
 0&1&-1&2&-2&2\\
 1&2&0&3&0&3\\
 0&2&-2&4&-4&5\\
 1&3&0&5&-1&7\\
 2&6&-2&10&-5&13
 \end{pmatrix},
 \qquad \det=1.                                        \tag{C3}
\]
Any d-dimensional linear representation factors this matrix through
dimension d. Therefore d>=6 over Q and over every field, including F2.
C1 attains six with integer matrices.

More strongly, the reachable-column matrix R from P and observation-row
matrix O from S each have determinant -1. Their inverses are integral.
Hence the six suffix observations recover all six coordinates over Z and
over every modulus, without a division by a nonunit. The exact matrices and
O inverse are in [jump_carrier.json](jump_carrier.json).

For an even modulus, observe the integer half-jump J before reducing it.
Equivalently compute D modulo twice that modulus. Reducing H(T+) and H(T-)
modulo 2 alone gives no information about J modulo 2; division by the
nonunit 2 is not part of this residue decoder.

This is a minimal **linear dimension**, not a claim that six scalar queries
are optimal among arbitrary nonlinear decoders or that six machine words
are necessary for storing an actual source.

## C3 The edge derivative still remembers the whole source word

Given the six observations, apply O inverse. Check that the first
coordinate is 1, M,N are nonnegative, and a,b are bits. Form
\[
 Q_w=2M+a,\qquad Q_u=2N+b.
\]
Align their binary expansions, pad the shorter with leading zeroes, and
interleave the bits. Replay C1 and check A. This accepts exactly the actual
six-observation profiles and recovers the word modulo common leading 00
columns.

Indeed common leading 00 columns are universal-source padding, invisible
to all future append-and-J observations. Conversely equal future
observations give equal six probes, hence equal codes and the same normalized
word. Thus equality modulo leading 00 is the coarsest future equivalence
for this language.

The code pair is the same actual-word coordinate used by the earlier
[twelve H-probe decoder](../session134_constructive_spectrum_frontier_20261101/two_apex/THEORY.md).
Here it is recovered from six **edge-difference** probes. Each such probe
compares two edge orientations; it is not one old H query. The full integer
lattice is still larger than the actual orbit: z=(1,0,0,0,0,1) is an integral
state, but zero codes force the empty word and area zero, so replay rejects it.

When the original labels, source padding or physical vertex clock matter,
retain the length and named apexes alongside this quotient. A finite
precision residue state does not reconstruct the original unbounded words.

## C4 Dyadic precision is local in separation

Modulo 2^e, interactions at separation >=e+2 vanish in J because their
coefficients contain 2^e. The linear term depends only on the last e+1
columns. Equivalently reduce the six-coordinate matrices of C1 modulo 2^e.

[geometry.py](geometry.py) implements a separate streaming version:
keep the last e+1 columns, accumulate the surviving wedge terms, and keep
the binary difference modulo 2^(e+1) before the final exact division by two.
It agrees with the full contraction formula on every pair of words of
length <=5 at precisions e=1,...,4. C1 proves arbitrary-length closure;
the power-of-two coefficients prove the separation cutoff.

This connects the older digital/FRACTRAN viewpoint to a specific invariant:
binary shifts transport the prefix codes, and a skew area coordinate records
what the two code streams do together. The exact carrier grows in bit size;
finite precision creates finite state spaces without silently losing exact
source identity.
