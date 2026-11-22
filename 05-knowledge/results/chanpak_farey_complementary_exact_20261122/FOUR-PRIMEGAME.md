# The four-vertex switch square and the priority guard at 13

The [PRIMEGAME return analysis](../metagraph_primegame_20261102/PRIMEGAME_FANO.md)
distinguishes arithmetic feasibility, return words, and
priority-legal executions. The four-vertex tournament gives a
small **literal** analogue in which forgetting the state hides
which arcs actually have to change.

Take the transitive order \(0<1<2<3\), keep four arcs fixed, and
independently reverse \((0,2)\) or \((1,3)\). The two toggles
form a physical square:
\[
\begin{array}{c|c|c}
(a,b)&H&\text{ordinary isomorphism type}\\\hline
(0,0)&1&\text{transitive}\\
(1,0)&3&\text{cycle with an ordinal sink}\\
(0,1)&3&\text{cycle with an ordinal source}\\
(1,1)&5&\text{strong}.
\end{array}
\]
Thus \(H=1+2(a+b)\) on this square. The two height-three
tournaments are not isomorphic, and are **not** single-arc
neighbors even after relabelling: all six flips from either
type go to height one or height five. From either height-one
or height-five type there are three flips to its own ordinary
type, one to each height-three type, and one to the other
odd height. The diagonal \(1\leftrightarrow5\) elsewhere in the
quotient therefore does not make the two height-three states adjacent.

As a toy encoding, store \((a,b)\) as \(11^a13^b\). The
reciprocal fractions \(13/11\) and \(11/13\) take the one-hot
codes into one another. In the square they interchange the
two different height-three states, which requires **two**
literal flips, via height one or five. This does not identify
Conway's whole machine with the four tournaments; it shows
exactly which state a collapsed height would lose.

On the actual [133-vertex lower witnesses](../chanpak_decimal_kernel_20261122/BACKCHAIN-LOWER.md),
the common count has the arithmetic factorization
\[
 H=E_{61}-1=(2^{60}-1)(2^{61}+1).
\]
In particular \(v_{11}(H)=v_{13}(H)=v_7(H)=1\), using orders
10, 12 and 3 and the corresponding elementary prime-power
lift. Both reciprocal clock fractions are integral there.
But in [Conway's ordered PRIMEGAME](https://gwern.net/doc/cs/computable/1987-conway.pdf)
the **first** applicable fraction is \(17/91\) since \(7\cdot13\mid H\).
At the desired upper count \(H+2\), neither 11 nor 13 divides,
and exact evaluation gives the default \(55\) as the first
applicable fraction. So the lower integer's divisibility at
13 is not a priority-legal \(11/13\) step, and does not turn
the known modular upper exits of the 133-vertex host into a
literal \(+2\) edge. The contact-phase lifts across varying
Euclid exponents are recorded separately in
[CONTACT-CLOCKS.md](CONTACT-CLOCKS.md).
