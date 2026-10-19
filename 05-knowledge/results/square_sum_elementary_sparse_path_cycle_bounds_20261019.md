# Elementary polynomial-size sparse square-sum trees and cycles

**Status: PROVED.** This comparison identifies exactly what the AP-paired/common-center restriction adds. Sparse square-sum paths and trees without that restriction already have polynomial *numerical* labels by a greedy construction. A two-vertex closing step gives the same conclusion for cycles. Consequently an elliptic recurrence's polynomial bit bound is structural information about that recurrence, rather than a label-size breakthrough for unconstrained sparse paths.

**Literature status.** Ordinary consecutive-interval square-sum existence is already established: Robert Gerbicz records paths for every n≥25 and cycles for every n≥32 in his January21,2018 [OEIS A090461 contribution](https://oeis.org/A090461), linking his deterministic construction. Any interval-coverage obstacle mentioned here concerns extending this particular structured construction or preserving additional prescribed data, not an open claim of ordinary existence.

## 1. Greedy tree extension

Suppose distinct positive integer labels have current maximum \(M\ge1\), and a new vertex is to be joined to a previously labeled parent with label \(b\le M\). Set
\[
r=\left\lceil\sqrt{M+b+1}\right\rceil,\qquad c=r^2-b.
\]
Then \(c\ge M+1\), so the new label is positive and different from every previous label, and \(b+c=r^2\). Moreover,
\[
c\le M+2+2\sqrt{2M+1}\le(\sqrt M+2)^2.
\]
Thus, starting from any finite labeled seed of maximum \(M_0\), after \(j\) new vertices,
\[
M_j\le(\lceil\sqrt{M_0}\rceil+2j)^2.
\tag{1}
\]
Any finite rooted tree can be labeled in parent-before-child order. Existing seed labels and edges remain unchanged. The construction prescribes a copy of the tree in a square-sum graph; additional square-sum edges are allowed. In particular, any positive integer square-sum path can be extended by arbitrarily many vertices with \(O(j^2)\) numerical labels for a fixed seed, hence \(O(\log j)\) bits per label.

## 2. Close any such path with exactly two new labels

Let a positive integer square-sum path have distinct endpoints \(A,B\), and let \(M\) be its maximum label. Choose
\[
r=\left\lceil\sqrt{32M}\right\rceil
\quad\text{or}\quad
\left\lceil\sqrt{32M}\right\rceil+1
\]
so that \(r\equiv A+B+1\pmod2\). Define
\[
s=\frac{A+B+1-r^2}{2},\qquad
v=r^2-B,\qquad w=s^2-A.
\tag{2}
\]
All quantities are integers, and
\[
B+v=r^2,\qquad v+w=(s-1)^2,\qquad w+A=s^2.
\tag{3}
\]
Thus appending \(v,w\) closes the path to a square-sum Hamiltonian cycle on its enlarged label set.

There are no collisions. Since \(32M\le r^2<60M\),
\[
31M\le v<60M,\qquad
|s|\ge\frac{30M-1}{2}\ge\frac{29M}{2}.
\]
Hence \(w>60M>v>M\). Also \(|s|<r^2/2<30M\), so
\[
\max(v,w)<900M^2.
\tag{4}
\]
Combining (1) and (4), extending a fixed seed path and then closing it produces cycles with \(O(N^4)\) numerical labels and \(O(\log N)\) bits per label. All steps use integer square roots and elementary arithmetic; no factorization or search for primes is needed.

### Sharper bound when a maximum label is an endpoint

This is always the case after greedy path extension. Orient that endpoint as \(A=M\), and set \(K=\lceil\sqrt M\rceil\). Choose
\[
r=\left\lceil\sqrt{M+B+16K+1}\right\rceil
\]
with a possible increment of one to enforce \(r\equiv M+B+1\pmod2\). Define
\[
\delta=r^2-M-B-1,\quad s=-\delta/2,\quad
v=r^2-B=M+1+\delta,\quad w=s^2-M.
\tag{5}
\]
The three square-sum identities (3) still hold. Since
\[
16K\le\delta<16K+4\sqrt{2K^2+16K+1}+4<40K,
\]
we have
\[
M<v<42M<63M\le w<400K^2\le1600M.
\tag{6}
\]
Here \(K\le M\) and \(K^2\le4M\) for integers \(M\ge1\). Thus the closure increases the maximum label by only an absolute constant factor. Reversing or cyclically rotating the resulting cycle restores the seed's original order if desired.

Consequently the greedy-path construction followed by (5) produces square-sum cycles with \(O(N^2)\) numerical labels, and still \(O(\log N)\) bits per label. The more general bound (4) is useful when neither endpoint has the maximum label.

## 3. The exact seed 482, 3362, 6242

The seed has
\[
482+3362=62^2,\quad3362+6242=98^2,\quad482+6242=82^2.
\]
It is already a three-cycle. A four-cycle retaining the same ordered first three labels is
\[
482,\quad3362,\quad6242,\quad2222,
\]
since \(6242+2222=92^2\) and \(2222+482=52^2\). For every \(N\ge5\), greedily extend the seed path to \(N-2\) vertices and apply (5), orienting the maximum endpoint first. This gives an integer square-sum cycle of every size \(N\ge3\) containing the exact, unscaled seed triple. Since \(\lceil\sqrt{6242}\rceil=80\), the generic construction has maximum label less than
\[
1600\bigl(80+2(N-5)\bigr)^2\qquad(N\ge5).
\]

## 4. Boundary of the conclusion

The newly added labels are not required to be AP-paired petals around a common center, and the label set is not an interval of consecutive integers. Therefore this theorem does not contradict the common-center capacity obstruction or the finite rational-contact-atlas height obstruction. Nor does it establish graceful labeling, whose edge-difference requirements are different. It shows why those additional structural requirements, rather than unconstrained sparse path/cycle existence, carry the hard growing-chain problem.

The fixed elliptic interpolations retain value as explicit group-law recurrences, seed-parametrization bridges, and arithmetic mechanisms. Their height bounds should be compared with the elementary bound above when discussing unconstrained sparse label sets.

## Exact controls

`verify_elementary_sparse_square_sum_20261019.py` checks greedy paths and rooted trees, the two-label identity and collision bounds, the exact four-cycle, and cycles retaining the seed triple. It uses only the Python standard library, and normal and `python -O` outputs agree.

From the repository root run `python 05-knowledge/experiments/verify_elementary_sparse_square_sum_20261019.py`. The receipt is [elementary_sparse_square_sum_controls.json](elementary_sparse_square_sum_controls.json).
