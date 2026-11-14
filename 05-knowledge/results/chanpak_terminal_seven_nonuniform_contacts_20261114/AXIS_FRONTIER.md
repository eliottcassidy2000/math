# Two-sided division by seven at a genuine nonreset seam

The [starting terrace](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/RESETS_AND_DECODER.md) has two-channel transfer \(D_yC_q\), with
positive width powers \(y=2^r\), \(r\ge1\), and symmetric
\(C_q=(e_q,f_q;f_q,g_q)\). A physical prefix ends in state
\(z=(a,b)^t\), where \(a\) is odd, \(b\) is even, and \(b\ge2a\).
A right suffix is read backwards from \(R=(1,1)^t\) via
\(R\leftarrow C_qD_yR\). Thus a seam readout is \(z^tC_qR\)
if the prescribed right width is already inserted into \(R\).
All operations here concern ordinary (unrooted) counts on the
strong minimum-triangle terrace.

## Exact seven-times-norm output and an elastic alternative

Let \(w=(47,56)^t\), the axis whose outer product is \(C_{13}\).
The invertible-over-\(\mathbb Q\) matrix
\(C_7=(57,68;68,81)\), determinant \(-7\), obeys
\[
 C_7(1,4)^t=7w,\qquad C_{11}(1,2)^t=47w.
\]
The shorter first probe gives exactly \(7(47a+56b)\) from **any**
physical prefix. This is not an ordinal product with a seven-path
tournament, nor a rank-one seam. For example \(z=(1,4)^t\) yields
\(7\cdot271\). More generally, at \(q=7\), every
\(y=4\cdot8^j,\ j\ge0\), gives
\[
 H=7(A_y a+B_y b),\qquad
 (A_y,B_y)=((57+68y)/7,(68+81y)/7).
\]
This is an infinite menu of positive primitive coefficient rows:
if \(d\) divides both coefficients, then
\((1,y)^t=-\operatorname{adj}(C_7)(A_y,B_y)^t\) implies \(d=1\).
Here \(A_y\) is odd and \(B_y\) even. The old \(w\)-row is the
smallest member, not the only possible cofactor functional.

In fact modulo 7, \(C_7=(1,5)^t(1,5)\). Exact divisibility
by 7 of any \(q=7\) sandwich means **one side or the other**
lies in \(a+5b\equiv0\pmod7\). In particular, at a terminal
\(q=7\) with right power \(y\not\equiv4\pmod7\), normalize the
left side instead:
\[
  (U,V)^t=C_7z/7,\qquad H/7=U+yV,\qquad
  z=\bigl(68V-81U,\;68U-57V\bigr)^t.
\]
The inverse is integral, and this linear map gives exactly the
lattice \(a+5b\equiv0\pmod7\), not merely one coset criterion.
For example the seed \(z=(1,4)^t\) lies in that lattice; ending
at \(q=7,y=4\) uses both sides, whereas \(y=64\) uses only the left and gives \(7\cdot3631\), with prime cofactor.

## A complete target-directed decoder for terminal \(q=7\)

Given odd prime \(p\ne7\), ask whether \(7p\) has a terrace witness
whose terminal seam is \(q=7\). There are \(O(\log p)\) allowed
width exponents: for \(H=57a+68b+y(68a+81b)\) one has
\(H\ge193+230y\).
* If \(y\equiv4\pmod7\), enumerate \(a\equiv A_y^{-1}p\pmod{2B_y}\)
  in \(1\le a\le p/(A_y+2B_y)\), then set \(b=(p-A_ya)/B_y\).
  The modulus enforces parity and is valid because the row is primitive.
* Otherwise, enumerate **even** \(V\) in
  \[
  \frac{81p}{81y+68}<V\le \frac{230p}{230y+193},
  \]
  then \(a=(81y+68)V-81p\), \(b=68p-(68y+57)V\).
  The interval is exactly \(a>0,\ b\ge2a\); the parity follows
  from \(p\) odd and \(y\) even. The interval width is only
  \(7p/[(81y+68)(230y+193)]\), so large widths give at most one
  candidate rather than scanning a broad set of raw counts.

Check each candidate by inverting the two-state word exactly.
Neither reset \(q=0\) nor \(q=13\) can occur anywhere in a
\(7p\) witness: zero factorizes into ordinary counts on two
nontrivial terrace words (one would need count 7), and 13
factorizes into two glue readouts each at least 159. Hence over
the required prefixes all corridor determinants are nonzero.
At a state \((a,b)\), try \(1\le r\le v_2(b)\) and invert
\(C_q(a',b')=(a,b/2^r)\), testing positive integral preimages
with \(a'\) odd, \(b'\) even, \(b'\ge2a'\). Only \(q\) with
\(e_q+2f_q\le a\) matter. It strictly decreases \(a\), and the
initial state is exactly \(a=1,\ b\) a positive power of two.
This yields a finite exact decision and actual word, not a search
that assumes feasibility of the real cone.

## Complete classification of probes that preserve the glue axis

In the backwards recursion, all coordinates have odds \(v/u\)
in \([1,4/3]\). For \(q=1\) the outgoing odds
\(2(1+t)/(1+2t)\), \(t=yv/u\ge2\), lie in \((1,6/5]\).
At \(q=2\) the odds are \((3+4t)/(3+3t)\in[11/9,4/3)\).
For \(q\ge3\), the coefficient inequalities give the interval
\([7/6,11/9]\) because both \(f_q/e_q\) and \(g_q/f_q\) are
there: check \(q=3,4,5\) and propagate by the recurrence.
Only \(q=0\) emits odds 1, apart from the empty suffix.

To have final odds \(56/47\), the row
\((56,-47)C_q\) must vanish on \((1,t)^t\) with \(t\ge2\).
The finite rows through \(q=16\), followed by positive
Tribonacci propagation, leave exactly \(q=1,t=19/9\);
\(q=7,t=4\); \(q=11,t=2\); or \(q=13\) (any \(t\)).
The first option cannot occur. It requires predecessor odds
\(19/18\) via \(y=2\). Those odds only occur via \(q=1\)
with total incoming \(17/2\), hence require \(y=8\) and still
earlier odds \(17/16\). In turn the latter requires \(q=1\)
with total incoming \(15/2\), which is not a power \(y\)
times any odds in \([1,4/3]\). Both \(q=7\) and \(q=11\)
require preceding odds 1, so a **reset-free** probe proportional
to \(w\) consists of precisely one step, either \(q=7,y=4\)
or \(q=11,y=2\). Longer probes can reach the axis only by
passing through a reset. This is a rigidity of the two-state
geometry, not a statement excluding other tournament constructions.
