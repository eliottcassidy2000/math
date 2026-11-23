# Split and inert congruence phases in signed odd Collatz

For odd $x$ with $3x+\varepsilon\ne0$, write
$$
F_{\varepsilon,A}(x)=(3x+\varepsilon)/2^A,\qquad A=v_2(3x+\varepsilon),\quad \varepsilon=\pm1.
$$
For $\varepsilon=1$, $A=1$ iff $x\equiv3\bmod4$, and $A\ge2$ iff $x\equiv1\bmod4$. Call these the inert and split *congruence phases*, respectively: this is named for the familiar Gaussian-prime criterion when $x$ is a positive odd prime, not an assertion that every odd integer here is prime. For $\varepsilon=-1$ the phases are interchanged; $x\mapsto-x$ conjugates the two maps.

## Centered inert runs and exact finite addresses

Write $u=x+\varepsilon$. On an $A=1$ step,
$$
u'=(3/2)u.
$$
Apart from the fixed point $x=-\varepsilon$, a maximal run starting with $v_2(u)=r$ has exactly $r-1$ inert steps, ending in the other congruence phase. For every odd prime $p\ne3$, the valuation $v_p(u)$ remains unchanged along the run and its centered unit residue rotates by $3/2$ modulo $p^e$; at $p=3$, the valuation increases exactly one per step. If the initial centered residue modulo an odd $q$ has gcd $g$ with $q$ and $3\nmid q$, the multiplication phase repeats after $\operatorname{ord}_{q/g}(3\cdot2^{-1})$ steps, provided one is discussing a long enough inert run. The 2-adic countdown fixes the real run length independently.

More generally fix any finite exponent word $A_1,\ldots,A_m$ with $A_i\ge1$, and fix an odd terminal class modulo four. The inverse branch $x=(2^Ay-\varepsilon)/3$ on $2$-adic odd integers maps a residue class of precision $r\ge1$ to a single class of precision $r+A$, and gives exactly that valuation $A$ because $y$ is odd. Backward iteration yields one starting class modulo $2^{2+\sum A_i}$ prescribing precisely that itinerary and terminal class. For a prescribed starting residue modulo any odd modulus, CRT gives infinitely many positive representatives and infinitely many negative ones. Taking their magnitude sufficiently large keeps the entire finite itinerary in that sign, since all the slopes are positive and the word length is finite.

## Full odd modular control while never entering the inert phase

Use $\varepsilon=1$ (the other sign follows by conjugation). Every exponent word with all $A\ge2$ and terminal class $1\bmod4$ lies entirely in the split phase. Nevertheless it allows any odd modular destination that is a unit at three if the modulus is divisible by three.

If $q$ is odd and $3\nmid q$, the two affine maps $F_2(x)=(3x+1)/4$, $F_3(x)=(3x+1)/8$ are permutations over $\mathbb Z/q\mathbb Z$. Define $D=F_2F_3^{-1}$, using ordinary function composition right to left. Then $D(x)=2x$, and
$$
D F_2 D^{-1} F_2^{-1}(x)=x+1/4.
$$
This unit translation generates the whole additive group, so any prescribed starting and ending residues are connected by a word in $F_2,F_3$ and their inverses. Each inverse can be replaced by a positive power in the finite affine group. The resulting positive word is realized by the exact binary-address argument with prescribed terminal split class. In particular every transfer is attained by genuine finite trajectories in either sign, not necessarily by integer cycles.

For a general odd $q=3^kq_0$, $3\nmid q_0$, let $y$ be the target. For $k\ge1$, necessarily $y$ is a unit at three after any nonempty trajectory. Conversely append $k-1$ uses of $F_2$: since $F_2(x)-1=(3/4)(x-1)$, this forces the intermediate residue to $1\bmod 3^{k-1}$, regardless of the preceding orbit. Choose a final $A\ge2$ with $2^Ay\equiv4\bmod3^k$; this is possible because 2 generates $(\mathbb Z/3^k\mathbb Z)^\times$ (order two modulo 3 and $v_3(2^2-1)=1$). Then $F_A$ yields $y\bmod3^k$. All these appended maps remain invertible modulo $q_0$, so solve for the precursor endpoint modulo $q_0$ by the earlier translation word. Finally use the binary-address CRT with the original prescribed starting residue modulo all of $q$. This proves the exact allowed image even though repeated split steps are contracting in the integer size of large inputs.

## Exact arithmetic for a cycle with one split visit

Again set $\varepsilon=1$. One split exponent $A\ge2$ at $x\equiv1\bmod4$, followed by $k\ge0$ inert steps and a return to $x$, exists precisely when
$$
\Delta=2^{A+k}-3^{k+1}\quad\text{divides}\quad 2^A-2.
$$
Indeed the centered return is obtained by applying $u\mapsto(3u+2^A-2)/2^A$, then $u\mapsto(3/2)^k u$. Solving gives
$$
x+1=\frac{3^k(2^A-2)}{\Delta}.
$$
The denominator is odd and coprime to 3, yielding necessity and integrality under the divisibility condition. Integrality suffices for the exact branch signature: $v_2(x+1)=1$, and the proposed post-split value $y$ satisfies $y+1=2^k(x+1)/3^k$, an integer of valuation $k+1$. Thus $y$ is odd (so the split exponent is exactly $A$), and it has exactly $k$ inert steps before landing at the prescribed split value. The sign is that of $\Delta$, because $x+1$ is nonzero and even. For $A=2$, the only possibilities are $k=0,1$, giving the fixed value 1 or the odd negative cycle $-7\leftrightarrow-5$: for $k\ge2$ the denominator is negative with magnitude increasing from 11. For $A=3$ the cases $k=0,1,2$ have denominators 5,7,5, none dividing 6, and the magnitude thereafter increases from 17.

**Parity obstruction.** No such cycle has both clocks $A+k$ and $k+1$ even. This is the case $A,k$ both odd. If $k\ge A$, the denominator is a nonzero difference of squares and has magnitude at least $2^{(A+k)/2}+3^{(k+1)/2}>2^A-2$. If $k\le A-2$ and $k\ge3$, then
$$
2^A(2^k-1)\ge2^{2k+1}>3^{k+1},
$$
so again $\Delta>2^A-2$. For $k=1$, $2^{A+1}-9>2^A-2$ because $A\ge3$. Thus divisibility fails in every case. By sign conjugation the result also holds for the other numerator sign.

Finally, a cycle entirely in the inert phase has $u=(3/2)^m u$, so is only the corresponding fixed point. In the split phase, $(3x+1)/2^A<x$ for positive $x>1$, and $|(3x+1)/2^A|<|x|$ for every negative $x$. The only all-split cycle is the positive fixed point 1. These facts do not classify mixed cycles with two or more split visits, and congruence reachability is explicitly not a positive-integer-cycle claim.
