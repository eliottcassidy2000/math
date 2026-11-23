# A nonnegative collision bound that really can certify a whole window

Let \(I\) be a finite set of \(R\ge2\) admissible target counts.
Choose an actual finite bank of tournaments, with nonnegative
weights, and discard outputs outside \(I\). Write \(c_j\) for
the weight of the exact count \(j\in I\), \(M=\sum_j c_j>0\),
and \(E=\sum_j c_j^2\) for its collision energy. Then
\[
 E<\frac{M^2}{R-1}\quad\Longrightarrow\quad c_j>0\ \text{for all }j\in I .
\]
If one bin were empty, Cauchy–Schwarz on at most \(R-1\) bins
would give the opposite inequality. Equality at uniform positive
weights on \(R-1\) bins shows that the criterion cannot be weakened
without more information. With integer weights it can be sharpened
by distributing \(M\) as evenly as possible among \(R-1\) bins:
if \(M=(R-1)t+s\), \(0\le s<R-1\), the exact obstruction is
\[
 E\ge(R-1-s)t^2+s(t+1)^2.
\]
Thus an **aggregate** estimate can imply pointwise existence,
but only with a proven deficit threshold. The allowable excess
over \(M^2/R\) in the real-weight test is \(M^2/[R(R-1)]\):
a relative excess of order \(1/R\), not an unspecified small
equidistribution error.

## Exact collision count from paired apex contacts

The [weighted apex cut](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md)
has normalized form
\[
 F(x)=a\cdot x-\sum_{i<j}W_{ij}x_ix_j,\qquad x_i\in\{0,1\},
\]
so the absolute height is \(h+2F(x)\). Here \(W\) is half the
weight in the cited unnormalized cut. Write \(d=x-y\in\{-1,0,1\}^n\),
let \(D=\operatorname{supp}d\), and let \(u_j=x_j=y_j\) on its
complement. The product difference for two entries in \(D\)
is \((d_i+d_j)/2\); for \(i\in D,j\notin D\) it is \(d_i u_j\).
Therefore the collision equation becomes a linear subset-sum:
\[
 F(x)-F(y)=L(d)-\sum_{j\notin D}u_j B_j(d),\quad
 L(d)=\sum_{i\in D}a_i d_i
 -\sum_{i<j,\ i,j\in D}W_{ij}\frac{d_i+d_j}{2},\quad
 B_j(d)=\sum_{i\in D}W_{ij}d_i .
\]
For each \(d\), each solution in the free \(u_j\) corresponds to
one ordered pair of masks. The exact energy of the **full** contact
bank is the sum over the \(3^n\) differences of the coefficient of
\(Z^{L(d)}\) in \(\prod_{j\notin D}(1+Z^{B_j(d)})\).
It can also be read modulo a modulus by reducing exponents,
or applied to two simultaneous normalized cuts with two Laurent
variables. No speedup is guaranteed: many different subset sums can
occur, and extra work is needed if outputs are clipped to a window.
A modular collision bound can certify an exact window only if its
residue read is injective on the retained targets and the underlying
bank is restricted to them.

**Physical small test.** The strong four-vertex core with rows
\((8,1,3,6)\) has height 5. Its 16 apex completions have heights
in the five-element set \(I=\{5,9,11,13,15\}\): insertion does not
decrease height, order five has maximum 15, Rédei gives oddness,
and 7 is known absent. The difference formula gives energy
\(58<16^2/4=64\). That proves every target in \(I\) is hit using
only this one contact bank. Counting modulo a small modulus
would not by itself certify the absolute integers. The replay
compares the formula with direct enumeration for three physical
cores over the integers and at four moduli.

This is one precise analogue of a slope-style closing estimate.
It is a local certificate, not a claimed unbounded energy bound.
For the whole conjecture the [SCC-grade test](README.md) is better
at detecting a persistent obstruction from a fixed missing target.
