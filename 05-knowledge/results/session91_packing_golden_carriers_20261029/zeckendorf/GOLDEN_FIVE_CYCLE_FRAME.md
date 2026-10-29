# An integral golden frame inside the five-cycle

**Status: PROVED by the displayed integer identities; FINITE-EXACT checker.**
Session91 cross-lane bridge. This is a finite operator calculation, not a claim
about physical packing, graph isomorphism, or a solved open problem.

## Source, target, and exact map

Let `C5` have cyclically ordered vertices `0,...,4`, let
\(B=2I+\operatorname{Adj}(C_5)\) on \(\mathbb Z^5\), and let
\(L=\{x:\sum_i x_i=0\}\) be its integer augmentation lattice.
The golden register matrix is
\[
 A=\begin{pmatrix}0&1\\1&1\end{pmatrix},\qquad
 A^2=\begin{pmatrix}1&1\\1&2\end{pmatrix}.
\]
Put \(v=e_0-e_1\), \(u=e_1-e_2\). The column map
\[
 P=(v,Bv-v,u,Bu-u)=
 \begin{pmatrix}
 1&0&0&1\\-1&0&1&0\\0&-1&-1&0\\0&0&0&-1\\0&1&0&0
 \end{pmatrix}
\]
satisfies
\[
 BP=P\operatorname{diag}(A^2,A^2).
\]
All columns sum to zero, and the determinant of the first four rows is
\(-1\). Projection to the first four coordinates identifies \(L\) with
\(\mathbb Z^4\), so **P is an integral basis of L**, not merely a rational
frame. Consequently
\[
 (L,B|_L)\cong (\mathbb Z^2,A^2)\oplus(\mathbb Z^2,A^2)
\]
as integer modules with distinguished operators. In particular
\(B^2-3B+I=0\) on L.

This gives a precise bridge to the Fibonacci register in
[NORM_CONTROLLED_OBSERVERS.md](NORM_CONTROLLED_OBSERVERS.md): the five-cycle
augmentation action is two copies of its two-step clock. No source coordinate
is discarded by P. The map is not a Euclidean isometry or a map of binary-word
legality constraints; those structures are additional data.

## Where the index five actually occurs

Let \(c=(1,1,1,1,1)^t\). Then \(Bc=4c\), and the full frame
\(T=(c,P)\) has determinant \(-5\). Its image is precisely
\[
 \mathbb Zc+L=\{x\in\mathbb Z^5:\sum_i x_i\equiv0\pmod5\}.
\]
Indeed its sums are multiples of five; conversely subtracting
\((\sum x_i/5)c\) puts such an x in L. Thus **index five belongs to
separating the constant mode from the augmentation modes**, not to separating
the two golden blocks inside L. Over \(\mathbb Q\), or modulo any prime
other than five, this frame gives \(B\cong[4]\oplus A^2\oplus A^2\).
Recovering the rational constant coefficient requires the sidecar
\(\sum_i x_i/5\); integral separation is legal only when this is integral.

## Characteristic five: the constant mode enters augmentation

Reduce modulo five. The image \(\bar L\) has dimension four but now contains
c, because its coordinate sum is zero. Explicitly
\[
 P(2,1,3,-1)^t=c\pmod5.
\]
The frame T therefore has rank four. The golden-block decomposition **inside**
\(\bar L\) remains valid: its minor is still a unit. Only the ambient
constant/augmentation direct sum fails.

Write \(D=B+I\) (the shifted operator called \(D_2\) in the operator lane).
For the full five-dimensional space, direct multiplication gives
\[
 D^2\equiv J\pmod5,\qquad D^3\equiv0\pmod5,
\]
where J is the all-ones matrix. For example, writing
\(D=3I+\operatorname{Adj}(C_5)\), the diagonal, adjacent, and distance-two
entries of its square are respectively 11,6,1, all congruent to one.
Also \(Dc=0\), since each row of D sums to five. Exact row reduction gives
\[
 (\operatorname{rank}D,\operatorname{rank}D^2,\operatorname{rank}D^3)
 =(3,1,0).
\]
Its nilpotent Jordan blocks therefore have sizes **3 and 2**. On augmentation,
J vanishes, and \(D|_{\bar L}\) is two copies of \(A^2+I\), each nonzero
of rank one and square zero modulo five. Hence its blocks are **2 and 2**.
The ambient module is not the split sum of a constant size-one block and these
two size-two blocks: that would have rank two rather than three.

The missing split is concretely the length-three chain: for any standard basis
vector \(e_i\), \(D^2e_i=c\ne0\), whereas \(D^3e_i=0\). This is an exact
finite-field extension effect, not evidence that information disappears inside
the integral golden blocks.

## Verification and provenance

Run from this directory:

```sh
python3 check_five_cycle_frame.py
python3 -O check_five_cycle_frame.py
```

The standard-library checker uses active exceptions rather than assertions. It
verifies the integer intertwining identities, determinants -1 and -5, the
constant coordinates, the two rank profiles, the all-ones square identity, and
invertibility of T modulo 2,3,7,11,13. Both modes return identical PASS output.
All fields and dimensions here are explicit; the checker is not a numerical
eigenvalue experiment.

Inherited operator: Session91 operator lane's five-cycle singleton Gram block
\(2I+C_5\). Frame found independently in the parent and operator lanes and
verified here. Related metric boundary:
[packing/GOLDEN_METRIC_AND_PALEY_TORI.md](../packing/GOLDEN_METRIC_AND_PALEY_TORI.md).
No external novelty priority is claimed for this elementary decomposition.
