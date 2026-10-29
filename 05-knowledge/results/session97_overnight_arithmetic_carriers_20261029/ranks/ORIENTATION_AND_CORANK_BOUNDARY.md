# A tournament orientation cannot supply the high-rank matrix model by itself

**PROVED elementary linear algebra; exact controls.** This is a boundary on
one proposed bridge, not a rank theorem for elliptic curves. The arithmetic
model in [Park–Poonen–Voight–Wood](https://math.mit.edu/~poonen/papers/bounded-ranks.pdf)
uses integer alternating matrices. Their Remark1.1.4 explicitly distinguishes
a probabilistic model from a naturally attached matrix for each curve.
The session does not provide such an arithmetic attachment.

## 1. Every ordinary tournament dominance matrix has corank zero or one

Let M be an n-by-n integral skew-symmetric matrix with every off-diagonal
entry odd. This includes every ordinary tournament dominance matrix.
Modulo2, M=I+J where J is the all-ones matrix. If n is even then
(I+J)^2=I because J²=nJ=0, so det M is odd and M is invertible over Q.
If n is odd, an even principal(n-1)-minor is invertible, while skew symmetry
in characteristic zero forces det M=0. Therefore, for every n,

    corank_Q M = n mod2.

The n=1 case has rank zero; the empty case has corank zero. This extends the
[session96 odd-square discriminant boundary](../../session96_sha_pairing_growth_20261029/tournaments/PAIRING_GROWTH_MODES.md)
from even matrices to all sizes.

Consequently the raw rational corank of a tournament dominance matrix can
never be2 or larger. It cannot itself model the high-rank tail in PPVW.
This does not prohibit a different tournament-based construction with
specified ports, additional weights, another operator or another observable.
It requires that construction to be supplied and checked, rather than
identifying the word 'rank' across two unrelated carriers.

This rational statement is not a statement in every characteristic. A
four-vertex matrix can have determinant9, remain invertible over Q, and
lose rank modulo3. Changing field changes the observable.

## 2. One fixed orientation supports every positive even matrix rank

The opposite boundary is equally sharp. Fix the transitive orientation on
vertices1,...,n. For any1<=r<=floor(n/2), set

    M^(r)_(ij)=(j-i) sum_(k=0)^(r-1) (ij)^(2k),

for every i,j. This is skew-symmetric and its upper-triangular entries are
strictly positive integers. Thus all these matrices have the **same**
underlying transitive tournament orientation.

Let V be the n-by-2r matrix with columns

    1, i, i², ..., i^(2r-1),

and let J_r be the block diagonal matrix of r copies of
[[0,1],[-1,0]]. Direct expansion gives

    M^(r)=V J_r V^T.

The Vandermonde determinant on any2r distinct labels is nonzero. Hence V
is injective as a map Q^(2r)->Q^n and V^T is surjective in the opposite
direction. Since J_r is invertible, the composition has rank exactly2r.
Therefore this single orientation supports every possible positive even
matrix rank, with coranks

    n-2, n-4, ..., n-2 floor(n/2).

The entry sizes and parity pattern change and are not free. In particular
one cannot retain only signs, or silently require all entries to stay odd,
and preserve these kernels.

Already at n=4, the weights j-i give rank2 while the same orientation with
all upper entries one has rank4. Its Pfaffian cancellation is
1*1-2*2+3*1=0, whereas the unit-weight Pfaffian is1. No edge orientation
changed. The Hamilton-path count of the underlying transitive tournament
remains one in both cases and does not encode the weight cancellation.

## 3. Connection contract

Source: an integral alternating matrix with named entries. Target: its sign
tournament. Preserved: which ordered pair has positive entry. Lost: parity,
magnitude, polynomial cancellations and rational kernel. Required restoration:
the actual weights, or a separately proved sufficient carrier for the target
query. A genus, rank threshold or common integer21 supplies no such map.

The PPVW threshold21 comes from that model's height-count and rank-tail
balance, not from the21 unordered pairs of a seven-vertex tournament. The
present theorem rules out the simplest corank identification and supplies
an explicit all-order family showing exactly where additive information is
lost. It proves neither the boundedness heuristic nor its negation, and
constructs no new high-rank elliptic curve.

## Verification

`check_orientation_coranks.py` uses exact rational elimination, independently
constructed Vandermonde products, every labelled tournament through order5,
weighted examples through order12, and the rational/mod3 hostile. Universal
statements rest on the displayed proofs. Normal and optimized Python execute
the same explicit checks.
