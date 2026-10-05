# Child-cut repair: exact parameter intervals and the missing integer guard

2026-10-04. **PROVED:** the fixed-word parameter partition and exact
residue-product criterion below. **FINITE-EXACT:** the declared symbolic
universe and the inherited 92 pair diagrams. **REFUTED:** a proposed proof
using only positive rational values, completed routes, and record order.
**OPEN:** the universal integer child-cut repair question and Collatz.
No literature-priority claim.

## 1. Inheritance and precise question

The closest mechanism is the [supplied child-cut law](supplied_state_join_extensions_20261004.md):
a checked diagram \(n\to J\leftarrow m<n\) has lift slope \(\lambda\);
removing a child prefix of coefficient \(\mu\) changes it to
\(\lambda\mu\), while releasing a precisely known ternary source guard.
The canonical expanding hostile is \(233\leftarrow231\), repaired by
cutting at 31. The corrected near miss is \(75\leftarrow73\): cutting
at 71 decreases the child label but leaves the slope greater than one.
The least-used sidecar here is the endpoint's exact ternary integrality
class, together with its least positive odd representative.

The [mediant-tree note](collatz_mediant_tree_K_20260930.md) supplies
prefix affine carries and real threshold interpretations. Its bounded
glide records and named conjectures are not a universal bound on the
products used below. No such bound is imported.

The live board is: original source, child prefix, affine coefficient,
ternary endpoint guard, real size threshold, and first-hit root.
The anchor is supplied-source proof reuse; the niche is a finite
obligation attached to each word pair; the wildcard is the integer
versus rational boundary of the proposed record-minimum argument.

The inherited open question is: for completed positive odd integers
\(m<n\), when their first common endpoint has expanding lift slope,
must some value \(x<n\) on that child prefix have repaired slope at most
one? The present work does **not** prove or refute that integer statement.

## 2. Every compatible pair has a canonical positive parameter

For a positive valuation word \(w\) of length \(r\), total cost \(A\),
and carry \(B_w\), write
\[
 F_w(X)=\frac{3^rX+B_w}{2^A}.
\]
Empty words have \(r=A=B_w=0\). A second word \(v\) has length \(s\),
cost \(D\), and carry \(B_v\).

An odd endpoint \(J\) admits both integer inverse words exactly when
\[
 J\equiv 2^{-A}B_w\pmod{3^r},\qquad
 J\equiv 2^{-D}B_v\pmod{3^s}.                         \tag{1}
\]
The classes are compatible exactly when their residues agree modulo
\(3^{\min(r,s)}\). Let \(H=3^{\max(r,s)}\) and let \(J_0\) be the least
positive odd integer in their common class. Then all positive common
endpoints and their sources are
\[
\begin{aligned}
 J(t)&=J_0+2Ht,\\
 n(t)&=n_0+Pt,\quad &P&=2^{A+1}3^{\max(s-r,0)},\\
 m(t)&=m_0+Qt,\quad &Q&=2^{D+1}3^{\max(r-s,0)},\qquad t\ge0,             \tag{2}\\
 n_0&=(2^AJ_0-B_w)/3^r,\quad&
 m_0&=(2^DJ_0-B_v)/3^s.
\end{aligned}
\]
Moreover \(0<n_0<P\) and \(0<m_0<Q\).

**Proof.** The congruences are the exact integer conditions. Integrality
of the whole inverse word implies integrality of each intermediate
inverse: reduce its numerator modulo three, remove the last inverse,
and induct. An integral inverse of a positive odd target under
\((2^aJ-1)/3\), \(a\ge1\), is positive and odd. Thus the entire inverse
word is an actual valuation word. Also \(1\le J_0<2H\), so the positive
seeds satisfy the strict upper bounds in (2). Every odd common endpoint
has the stated parameter; a negative parameter makes \(J\), and hence
every positive forward source, negative. This proves completeness.

This is an arithmetic statement with formal \(U(1)=1\). A strict
first-hit prefix additionally excludes an interior root. Only the
least parameter \(t=0\) can have this issue: every intermediate value
increases by a positive integer period when \(t\) increases. The
program's first-hit route controls retain that sidecar. A compatible
endpoint alone does not assert that it reaches the root.

## 3. Expanding pairs leave a finite initial obligation

Assume \(Q>P\). The original child is strictly smaller exactly for
\[
 0\le t\le T,\qquad
 T=\left\lfloor\frac{n_0-m_0-1}{Q-P}\right\rfloor.                    \tag{3}
\]
A negative \(T\) means there is no such positive realization. In
particular there are at most
\[
 \left\lceil\frac{P}{Q-P}\right\rceil
 =\left\lceil\frac1{\lambda-1}\right\rceil,\qquad \lambda=Q/P,
                                                                    \tag{4}
\]
smaller-child parameters for this fixed pair.

At child cut \(p\), let its prefix cost be \(C_p\). Along (2) the cut
value is
\[
 x_p(t)=x_{p0}+Q_pt,\qquad Q_p=Q\,3^p/2^{C_p}.
\]
The period \(Q_p\) is an integer, and the repaired slope is \(Q_p/P\).
If \(Q_p>P\), this cut cannot give a nonexpanding family. If \(Q_p<P\),
it gives both the required slope and strict source-size decrease
exactly when
\[
 t\ge L_p=
 \max\!\left(0,\left\lfloor\frac{x_{p0}-n_0}{P-Q_p}\right\rfloor+1\right).
                                                                    \tag{5}
\]
If \(Q_p=P\), it works for every parameter when \(x_{p0}<n_0\), and
for none otherwise.

Consequently, retaining **all** child cuts gives a single upper
parameter ray \(t\ge L_*\), where \(L_*\) is the least valid bound in
(5), or infinity when no cut qualifies. The unresolved parameters for
this pair are exactly
\[
 0\le t\le\min(T,L_*-1).                              \tag{6}
\]
Equations (3)--(6) follow from integer linear inequalities; equality
at a source-size comparison is correctly excluded.

This is an exact finite obligation per fixed pair of words, not a
uniform bound over all pairs. There are infinitely many possible
words, and neither their residual intervals nor universal child
certificates have been discharged. A repaired dependency still
consumes the certificate of its particular child. Cutting preserves
the original source and suffix provenance; it does not substitute a
convenient new input.

## 4. What a record minimum omits

For any checked positive route \(z_0,\ldots,z_j\), define the exact
residue product
\[
 R=\prod_{i=0}^{j-1}\left(1+\frac1{3z_i}\right).
\]
Then its endpoint ratio equals \(3^j2^{-A}R\). For the two routes
\(n\to J\leftarrow x\), with source product \(R_w\) and child-suffix
product \(R_v\), this gives
\[
 \boxed{\lambda_x=\frac{x}{n}\frac{R_v}{R_w}},\qquad
 \lambda_x\le1\iff \frac{R_v}{R_w}\le\frac n x.        \tag{7}
\]
Thus \(x<n\), or even \(x<m<n\), only controls one factor. The missing
comparison is quantitative and path-dependent. Equation (7) is an
identity, not a new upper bound on residue products.

The following exact hostile shows that positivity and record order
alone cannot supply the missing argument:
\[
 n=J=3,\qquad
 m=\frac{77}{27}\xrightarrow{1}\frac{43}{9}
 \xrightarrow{1}\frac{23}{3}\xrightarrow{3}3.
\]
These are positive odd rationals (odd numerator and denominator);
the indicated two-adic valuations are exact. Every value is greater
than one. Both source and child continue to the first root through
\(3\to5\to1\). Yet \(m<n\), the lift slope is \(32/27>1\), the two
interior child states exceed \(n\), and the endpoint is equal to \(n\).
No allowed child cut repairs it.

For this word, \(m=(32J-19)/27\); the interval on which \(1<m<J\) is
\[
 23/16<J<19/5.
\]
But an integral child requires \(J\equiv20\pmod{27}\), whose positive
odd representatives are \(47+54t\). The lattice misses the entire bad
interval. At the first integer endpoint the child is \(55>47\).
This is a concrete repair by the integer guard, not evidence of a
universal guard-versus-interval theorem.

Scaling by 27 produces the exact integer path
\(77\to129\to207\to81\) under \(x\mapsto\operatorname{oddpart}(3x+27)\),
with the same obstruction relative to source 81. It changes the map.
It is therefore a hostile to constant-insensitive reasoning, not a
counterexample for the integer \(3x+1\) system.

## 5. Reproduction, finite scope and stopping reason

Program:
[child_cut_obligation_partition_20261004.py](../../04-computation/experiments/child_cut_obligation_partition_20261004.py).
Output:
[child_cut_obligation_partition_20261004.out](child_cut_obligation_partition_20261004.out).

    python 04-computation/experiments/child_cut_obligation_partition_20261004.py
    python -O 04-computation/experiments/child_cut_obligation_partition_20261004.py

An independent backwards-rational enumeration checks (1) for every
pair of words of cost at most four, including empty words, and every
odd endpoint through 81: 10496 cases. A separate complete symbolic
universe uses the 22 source words of length at most two and cost at
most six, including empty, against all 4096 child words of cost at
most twelve, including empty. Of 44600 expanding shapes, 15116 have
compatible endpoint guards and 11 have any positive smaller-child
seed. All 11 have empty intervals (6); their maximum seed count is
one. Actual replay and threshold-iff checks cover parameters zero
and one. These are finite controls of the formulas.

The inherited source universe is unchanged: all 501 odd integers
through 1001, supplied with completed certificates. Among its
2914 expanding first-coalescence pairs, 2822 already meet below the
source. The remaining 92 each normalize to \(T=L_*=0\): their only
smaller-child parameter is the observed seed, and a repairing cut
already applies there. All 92 complete fixed-word intervals are
therefore discharged, but this adds no new supplied integer coverage.
There are also 2930 exact checks of (7) over their child cuts.
The program uses explicit exceptions rather than removable assertions.

The attempt at a record-minimum proof stops at the unproved comparison
in (7). The rational hostile identifies an indispensable integer
guard, and (6) turns that guard-versus-size question into a finite
obligation for each chosen pair. A uniform integer argument emptying
all such residual obligations remains **OPEN**.
