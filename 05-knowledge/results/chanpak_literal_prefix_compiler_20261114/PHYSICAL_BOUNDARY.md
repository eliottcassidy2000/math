# Changed physical sources on the seven-axis

**Scope.** We give an exact inverse for the near-boundary part of the
source cone, construct an infinite set of real sources at the maximal
relaxed inverse depth, and exhibit a prescribed prime cofactor at which
the old least-positive choice fails but a different point on the
*same affine line* is an actual tournament prefix.
This does not claim all cofactor lines meet the physical set.

## 1. A literal prime example after changing the source

Use the [seven-axis construction](https://github.com/eliottcassidy2000/math/blob/9717032b8bf0afc5371763d0f26b3810d0332f0d/05-knowledge/results/chanpak_terminal_seven_nonuniform_contacts_20261114/AXIS_FRONTIER.md).
Its length-7, width-two terminal pair applies the row \((47,56)\)
and multiplies by seven. The prior maximal relaxed-depth algorithm
uses \(M=\begin{pmatrix}1&2\\4&4\end{pmatrix}\) and
\((u_k,v_k)=(47,56)M^k\).

For the prescribed prime \(p=1137067\), maximal depth is \(k=1\),
so the source line is \(271a+318b=p\). The least positive odd
choice is \((469,3176)\). The dyadic ratio-band test in §4 rejects
it outright, including the reset seams: its only possible width has
exponent 2, but \(3176>4\cdot469\cdot(4/3)\). Move to the next
lattice point on the **same line**:
\[
 (469,3176)+(636,-542)=(1105,2634).
\]
This is the physical gap word \((3,0^8,1)\):
\[
 (1,8)\xrightarrow{\ q=8,\ r=1\ }
 (105+125\cdot8,\ 2(125+149\cdot8))=(1105,2634).
\]
Here a positive gap \(r\) is a width exponent, not a number of
copies of the letter 1. To realize a gap word as a literal tournament,
put anchor vertices at 0 and the successive partial sums of \(g+1\)
over all gaps \(g\). Orient all pairs from the smaller to the larger
vertex except for the reversed edges between consecutive anchors.
The prefix has 15 vertices; independent subset path DP on these
literal arcs yields 3739 Hamiltonian paths. Append one \(q=1,r=1\)
step and the \(q=7,r=2\) terminal. The full 28-vertex tournament has
exactly \(7959469=7p\) Hamiltonian paths, checked by both transfer and
the independent interval-packing recurrence.

Thus choosing the least positive representative is provably not
sound as a physical selection rule, even when the maximal relaxed
depth itself is correct. The coordinate-pair ambiguity in this
example is not a change in the requested cofactor or in the suffix.

## 2. Exact positive geometry near the cone boundary

Call a physical state \((a,b)\), with \(a\) positive odd and
\(b\ge2a\) even. Put \(d=(b-2a)/2\).
If \(d>0\) and \(a>6d\), the following is an equivalence:
\[
 \boxed{
 (a,2a+2d)\ \text{is physically attained}
 \quad\Longleftrightarrow\quad
 \left(d,\frac{a-d}{2}\right)\ \text{is physically attained}.}
\]
Moreover its last corridor and width are uniquely \(q=1,r=1\).
The forward step is \(M(d,B)^t=(d+2B,\,4d+4B)^t\), so the
physical direction is constructive.

For necessity the state is not an initial state (\(a\ne1\)).
If the last width has exponent \(r\ge2\), then \(b\ge4a\)
and \(d\ge a\), impossible. At width exponent 1 let the
predecessor be \((A,B)\). Then
\[
 a=e_qA+f_qB,\qquad d=(f_q-e_q)A+(g_q-f_q)B.
\]
At \(q=0\) we have \(d=0\). At \(q=2\) we have \(d=B\) and
\(a\le(9/2)d\) since \(B\ge2A\). At every \(q\ge3\) the
coefficientwise inequalities
\(e_q\le6(f_q-e_q)\) and \(f_q\le6(g_q-f_q)\)
hold: check \(q=3,4,5\) and propagate by the positive
tribonacci recurrence. Thus only \(q=1\) is possible, which
has \(d=A,\ a=A+2B\). The proof includes either reset seam. The constant 6 is sharp as a slope: a final \(q=3,r=1\) applied to the one-gap predecessor \((1,2^s)\) gives \(d=1+2^s\) and \(a=6d-1\), approaching ratio 6 with a different physical final corridor.
As immediate consequences:
- On the shell with \(d=1\), above \(a=6\) the complete list
  is \(a=1+2^{r+1}\), \(b=2a+2\), \(r\ge2\). The points are
  built by first taking width exponent \(r\) and then \(q=1,r=1\).
- On any positive even-defect shell, the whole range \(a>6d\)
  is nonphysical, since it asks for a predecessor with even first
  coordinate. This is an obstruction beyond real cone inequalities.
- For any fixed odd \(d\), all large points on its shell are
  explicitly obtained from the second-coordinate fiber at \(d\).
  That fiber is a finite union of sequences \(c2^t,\ t\ge1\):
  for a noninitial state the equation \(e_qA+f_qB=d\) bounds
  \(q,A,B\) independently of the final width, and for each physical
  predecessor the new second coordinate is \((f_qA+g_qB)2^t\).
  For \(d=1\) take \(c=1\). Thus the near-boundary inverse is a
  finite set of exact divisibility and power-of-two tests at each
  fixed defect, not a broad cone scan.

In a target-first implementation, at outer row \((u_k,v_k)\)
put \(A_k=u_k+2v_k\). On the shell \(d\) we must have
\[
 B=\frac{p-(u_k+4v_k)d}{2A_k}
\]
equal to the second coordinate of an actual source with first
coordinate \(d\); on success append \(k+1\) copies of \(q=1,r=1\)
and the terminal seven seam. If the defect is not fixed but one
searches \(0<d<A_k\), it is uniquely pinned by
\(p\equiv2v_kd\pmod{A_k}\) since \(\gcd(2v_k,A_k)=1\).
These tests prove exact attainment inside the shell and give a
reasoned rejection when the chosen defect cannot be built.

## 3. Both square-root coordinates can be real at maximal depth

Take any particular physical word with state \((d,c)\). Increasing
only its last positive gap realizes \((d,c2^t)\) for every
\(t\ge0\). The \(q=1,r=1\) step then gives
\[
 a=d+2c2^t,\qquad b=2a+2d,\qquad
 p=(u_k+4v_k)d+2A_kc2^t.
\]
For \(k\) sufficiently large, the \(t=0\) value is less than
\(A_k^2\), since the fixed stem and the outer linear row create a
quantity only \(O(A_k)\). Choose the least \(t\) that reaches
\(A_k^2\). Doubling the second term shows
\[
 A_k^2\le p<2A_k^2<A_{k+1}^2
\]
because \(M(1,2)^t\ge5(1,2)^t\).
Thus the physical square-root source exists *at the maximal
relaxed depth*, not merely after reducing the working depth:
\(a,b\) are both \(\Theta(\sqrt p)\) and the suffix really has
\(k=\Theta(\log p)\) outer steps. Choosing different seed words
changes \(d\); for instance \(d=1,5,9\) give three different
families by the starting words \((1),(1,0,1),(2,0,1)\).
This is an unconditional family of odd cofactors, not a claim
of infinitely many prime terms in those families.

Run physical_seed.py for the literal arc construction,
the full all-seam exact inverse (including 0 and 13), an independent
path DP for the concrete prefix, full count checks, and changed-stem
maximal-depth constructions.


## 4. The ratio determines the last exponent and makes the line inverse smaller

Every state ending with width exponent \(r\) satisfies
\[
                      2^r\le b/a<(4/3)2^r .
\]
For an initial state or the \(q=0\) seam the inner ratio is 1.
For \(q=1\) it is at most \(6/5\) when the predecessor ratio is
at least 2. At \(q=2\) the matrix \(C_2=(3,3;3,4)\) has output
ratio strictly below \(4/3\) as both predecessor coordinates
are positive. For \(q\ge3\), both coefficient ratios \(f_q/e_q\)
and \(g_q/f_q\) are at most \(11/9\): check \(q=3,4,5\) and
propagate the linear inequalities. Thus the last exponent is
uniquely determined by the integer part of \(\log_2(b/a)\),
even before finding a corridor. This includes singular seams.
If the dyadic band or divisibility \(2^r\mid b\) fails, there is no
physical prefix.

On a prescribed cofactor line \(ua+vb=p\) at row \(w^tM^k\),
for every trial \(y=2^r\) the only candidate first coordinates
satisfy
\[
              \frac{3p}{3u+4vy}<a\le\frac{p}{u+vy},
              \qquad a\equiv p u^{-1}\pmod{vy}.
\]
These are necessary ratio-band and integrality filters; the latter
also enforces odd \(a\) for odd \(p\). The full state solver tests
only these thin progressions and reconstructs literal words.
The interval length divided by its stride \(vy\) is exactly
\(p/[(u+vy)(3u+4vy)]<p/(4v^2 y^2)\), so summed over all
widths the branch count is at most \(p/(12v^2)+O(\log p)\).
One can combine this filter with a uniform cap from the
original straight line: allowed odd \(a\)'s are spaced by \(2v\)
and \(a\le p/A_k\). The row obeys \(5v_k/6\le u_k\le v_k\)
and \(A_{k+1}/A_k\le97/17\), by direct substitution in the
two-coordinate recurrence. At maximal depth there are
therefore at most \(1+\lfloor p/(2v_kA_k)\rfloor\le49\)
initial candidates even before imposing the dyadic bands.
This is a complete filter for the seven-axis row, not a claim
that every line clears the remaining corridor test.
