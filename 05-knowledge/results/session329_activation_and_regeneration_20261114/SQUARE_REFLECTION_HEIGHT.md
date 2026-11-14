# When a square identity first fits on the positive integer line

A cycle is not just a formal identity on its edge labels. Each edge on square \(s\) implements the reflection \(R_s(x)=s-x\); the composite may return exactly while every lift hits a nonpositive vertex or fails the height budget. The six-cycle analysis in [the switch note](SQUARE_SUM_SWITCH.md) is an instance. A second useful family comes from the earlier identity \(1^2+7^2=2\cdot5^2\).

For increasing square labels \(a<b\le c<d\) satisfying \(a+d=b+c\), use the cyclic order of edge labels \(a,b,d,c\). Its positive four-cycle candidates are
\[
p,\ a-p,\ b-a+p,\ c-p.
\]
The reverse parameter \(a-p\) instead exchanges the roles of \(b,c\) in the displayed construction when \(b\ne c\); both choices are tested by letting \(1\le p<a\). Remove collisions and minimize \(\max(b-a+p,c-p)\). The two larger vertices balance near \(p=(c-b+a)/2\), except when positivity clips the optimum. For the inherited \(4+81=36+49\) example, this unconstrained optimum exceeds the admissible domain. Taking \(p=3\) gives the first height 46: the successive vertices are \(3,1,35,46\). This recovers the specific witness in the repository's [global order-46 result](../crossroads_family_20260926_squares.md).

For the scale family \((a,b,c,d)=(k^2,25k^2,25k^2,49k^2)\), \(b=c\) makes \(p\) and \(k^2-p\) reverse presentations of the same cycle. Thus the admissible distinct lifts, modulo reversal, are exactly
\[
1\le p\le\left\lfloor\frac{k^2-1}{2}\right\rfloor .
\]
For \(k=1\) there is no lift. For \(k\ge2\), the minimizing \(p\) is the largest allowed one and the first height is \(\lfloor49 k^2/2\rfloor+1\). At \(k=2\), this is 99: \(1,3,97,99,1\) has square sums \(4,100,196,100\). This is a cutoff for that identity, not a second claim about the first four-cycle of the full graph. Nor is the cycle automatically usable as a Hamiltonian repair; the spanning complement must still be checked.

So even a literal number-theoretic symmetry has two separate costs: returning algebraically, and fitting a distinct positive orbit in the window. It is the same gap between an arithmetic candidate for a tournament seam and a physically realizable source.
