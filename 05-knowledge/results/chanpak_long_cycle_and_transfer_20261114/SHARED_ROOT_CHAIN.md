# Shared-root splice chains: three exact channels and a binary barrier

## 1. The physical chain and its transfer

Let \(B_1,\ldots,B_m\) be nonempty disjoint tournaments with marked
vertices \(r_i\). Orient every cross-block arc from earlier to later,
except \(r_{i+1}\to r_i\) for consecutive roots. With \(m\ge3\) and
both endpoint blocks singleton, the chain is strong: its descending
root chain closes through the first-to-last edge, every internal
nonroot is reached directly from the first root, and every internal
nonroot reaches the last root directly.

For one block write \(h=H(B)\). Let \(p,s\) be its positive and
signed root-start subset-path responses, and \(q,t\) the same
responses in its converse. Put \(f=F_2(B-r)\), the sum of complementary
induced path-count products, and \(\delta=1\) for a singleton and
\(0\) otherwise. Empty subsets in the unrooted factors contribute 1.

**Transfer theorem.** The Hamiltonian count is the first coordinate
of \(M_m\cdots M_1(1,0,0)^t\), where
\[
 M_i=\begin{pmatrix} h_i&q_i&-t_i\\p_i&f_i&0\\s_i&0&-\delta_i\end{pmatrix}.
\]

**Proof.** Use the tournament expansion over packings of disjoint
directed odd cycles, weighting each nontrivial cycle by two. A cross-block
cycle spans an interval \([i,j]\), \(i<j\). Each internal cut has only
one backward arc. A cycle crosses a cut equally many times in each
direction, so every cut it crosses is traversed exactly once in each
direction. It uses the entire consecutive descending root chain from
\(r_j\) to \(r_i\). The remaining forward leg consists of a starting
path in \(B_i\), an ending path in \(B_j\), and optionally a path on
a subset of each intermediate \(B_k-r_k\). For neighboring extreme
blocks, the forbidden root-to-root forward seam would join two
singletons, which is even and removed by the parity filter. For
non-neighboring extremes, an internal root is already on the
descending leg, so no forbidden forward seam is used. The intervals
of different crossing cycles cannot intersect, since they would
share a root.

The interval weight is thus the unsigned sum minus the signed
sum, including the complementary path-count factors from the
unused vertices:
\[
 p_iq_j\prod_{i<k<j}f_k-s_it_j\prod_{i<k<j}(-\delta_k).
\]
Here the internal root contributes a minus sign of its own, and
the signed two-color cavity sum
\[
 \sum_{U\subseteq D}(-1)^{|U|}H(D[U])H(D-U)
\]
equals \(1\) for an empty cavity \(D\), \(0\) otherwise. To see this,
expand the two factors by the odd-cycle packing formula, then color
each resulting odd cycle or isolated vertex independently by its slot.
The signed color sum of each is \(1-1=0\). An untouched block contributes
\(h_i\); inside the interval all remaining cycles are already included
in the complement counts. Hence the total is a sum of products over
tilings by singleton block weights and disjoint crossing intervals.

Let \(F\) be the full prefix sum, and \(A,B\) the unsigned and signed
open-interval sums. The next left-end contributions are \(Fp_i,Fs_i\);
continuation through an internal block multiplies the respective
channels by \(f_i,-\delta_i\); closing reads \(q_i A-t_i B\).
These are exactly the stated matrix updates. \(\square\)

## 2. Positive coordinates, rank one, and phase controls

Put \(y=-B\), and then \(z=F-y,\ x=A-F\). The singleton implements the
nonnegative Tribonacci rule
\[
 (z,x,y)\longmapsto (z+x+y,y,z).
\]
For any nonsingleton, write \(e=-s,d=-t\) for the end and start path
counts (the known signed endpoint identity). In \((z,x,y)\) its matrix is
\[
 \begin{pmatrix}
 h-e+q&q&h-e+q-d\\
 p-h+f-q&f-q&p-h+f-q+d\\ e&0&e
 \end{pmatrix}.
\]
Every entry is nonnegative: \(p\ge h\) by splitting full paths at the
root, \(h\ge e,\ q\ge d\), and \(f\ge q\) by deleting the root from the
rooted member of the ordered two-cover. After an initial singleton, the
state is \((z,x,y)=(0,0,1)\), so the cone and \(A\ge F\ge y\) are invariant.

A transitive block of order \(m\ge2\) with its source or sink marked,
and \(f=2^{m-1}\), implements respectively
\[
 S_f(F,A,y)=(F+A-y,f(F+A),0),\qquad
 E_f(F,A,y)=(F+fA,F+fA,F).
\]
An interior mark of zero-based rank \(r\) gives
\[
 (F',A',y')=(1,P,0)(F+QA),\quad P=2^{m-r-1},\ Q=2^r .
\]
In fact, any block \(L\to\{r\}\to R\) with both sides nonempty and
\(r\) marked has rank one:
\[
 (F',A',y')=(H(R),F_2(R),0)(H(L)F+F_2(L)A).
\]
This follows from the ordinal factorization of path counts and
of ordered two-covers, with no signed endpoint channels for an
interior singleton component. The prefix scalar is at least three
when the chain starts at a singleton, and persists as a divisor
through every later letter. With a later closing singleton, the
complementary factor is at least three as well: normalized state
\((H(R),F_2(R),0)\) is in the invariant cone and the closing reader is
nondecreasing there. Thus these erasures cannot be used on prime
targets in this chain model. Two consecutive transitive-source
letters also give a rank-one product: \(S_f S_g\) outputs
\((W,fW,0)\) with \(W=F+A-y+g(F+A)\).

## 3. A universal binary obstruction in this physical chain

Keep a final singleton reader separate from the live prefix: its
completed count is \(N=\ell(F,A,y)=F+A-y=2z+x+y\). Assume the prefix
starts with a singleton. Two singleton letters send the reader to
\(6z+3x+5y\ge3N\). A single arbitrary nonsingleton letter sends the
reader to at least \(5N/2\). Indeed, if \(h\ge3\) its output is
\[
 (h+p-e)F+(q+f)A-dy\ \ge\ 3F+3A\ \ge\ 3N,
\]
using \(A\ge y,\ q\ge d,\ f\ge q\ge h\). Otherwise \(h=1\), and the
block is transitive (a directed triangle already contributes two in
the packing expansion). For \(S_f\), the new reader is
\(N+f(F+A)\ge3N\); for \(E_f\), it is
\(F+2fA\ge(5/2)(F+A)\ge(5/2)N\) by \(A\ge F\).
An interior mark gives \((1+P)(F+QA)\ge3N\) with \(P,Q\ge2\).
Singletons are nondecreasing for the reader. Hence these bounds
persist if letters are appended before the common final reader.

For \(N\ge3\), any append word other than the empty word and the one
singleton overshoots \(2N+1\). That singleton is unique, so **a
reusable pair of binary digits cannot live in prefix-preserving
extensions with the same singleton reader**, even allowing arbitrary
cores and arbitrary append lengths. This does not rule out changing
the prefix, changing the final reader, or other physical models.

## 4. An exact ternary defect

The two two-vertex orientations satisfy the precise difference
\[
 \ell(S_2(F,A,y))-\ell(E_2(F,A,y))=z-x,
\]
so \(z-x\) is a phase defect at the absolute target, not a
congruence. After \(E_f\) and \(S_2\), the live state is
\((F+2fA,4F+4fA,0)\), with completed count \(5F+6fA\) and
defect \(-2F\). The subsequent choice \(S_2\) or \(E_2\) yields
\[
 15F+18fA\quad\text{or}\quad17F+18fA.
\]
The signed seam makes the gap \(2F\), independently of
\(f=2^j\) for arbitrary \(j\ge1\).
From the initial singleton \(F=A=1\), both
\(18\cdot2^j+15\) and \(18\cdot2^j+17\) are realized by strong
tournaments. This is not a ternary digit return: after \(S_2\)
the defect becomes exactly \(-2y\) and the signed reservoir becomes
zero; after \(E_2\) the coordinate \(x\) becomes zero. A repeatable
continuation would need to refill a reservoir or change its readout.

## 5. Changing the terminal reader gives a dyadic displacement

Replace the final singleton by a transitive sink-marked block of
size \(r+1\), \(r\ge1\). Its count readout is \(F+2^r A\). So varying
the terminal width doubles the displacement from \(F\), rather than
doubling the total. A starting singleton and this terminal block
still give a strong tournament when there are at least three
blocks: each nonroot of the last block reaches its sink, and the
first root reaches the last-block vertices directly.

The intercept \(F=1\) has a literal obstruction. Any three consecutive
roots form a directed triangle, so a prefix of count 1 has at most two
blocks. In the one-block case it is transitive and \(A=p\) is a power
of two. In the two-block case each block must be transitive. An
outneighbor of the left root or an inneighbor of the right root
would form a triangle with the two roots and the backward seam.
Thus the left root is its block's sink and the right root its
block's source. The left \(A=1\); on the right \(p=f\), so the live
\(A=2f\) is again a power of two. Hence the \(F=1\) slice cannot
choose arbitrary dyadic mantissas; it only gives \(1+2^k\).
