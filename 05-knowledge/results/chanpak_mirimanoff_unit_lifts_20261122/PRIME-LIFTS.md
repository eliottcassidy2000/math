# Exact prime-adic branch lifts on the perfect rails

**Proved arithmetic**, not an all-exponent tournament straddle. Put \(n=2r+1\ge3\), \(t=4^r\), and \(E_n=2^{n-1}(2^n-1)=2t^2-t\). All claims below concern the actual target integers \(E_n\pm1\).

## 1. The finite carry is the Collatz branch's Teichmüller defect

Write \([a]\in\mathbb Z_\ell\) for the Teichmüller representative of \(a\in\mathbb F_\ell\), including zero. For an odd prime \(\ell\ne3\), put \(q_\ell(c)=(c^{\ell-1}-1)/\ell\bmod\ell\) for unit \(c\), and \(H_\ell(X)=((1+X)^\ell-1-X^\ell)/\ell\), the classical finite carry. For \(A_j(x)=(3x+1)/2^j\), \(b=\overline A_j(a)\), one has
\[
\frac{2^j(A_j([a])-[b])}{\ell}
\equiv j(3a+1)q_\ell(2)-3a q_\ell(3)-H_\ell(3a)\pmod\ell.
\]
Indeed for an integral representative \(x\) of \(a\), let \(y=(3x+1)/2^j\); then \([a]\equiv x^\ell\), \([b]\equiv y^\ell\bmod\ell^2\), including zero residues. The numerator is \(3x^\ell+1-2^{j(1-\ell)}(3x+1)^\ell\bmod\ell^2\); expanding gives the formula. The \(j q_\ell(2)\) term distinguishes the branch index, and dropping both Fermat-quotient terms is not legitimate in general.

For any composable integral \(\ell\)-adic affine maps \(A,B\) put \(e_A(a)=(A([a])-[\bar A(a)])/\ell\in\mathbb Z_\ell\). Then the exact identity \(e_{B\circ A}(a)=B'e_A(a)+e_B(\bar A(a))\) follows by adding and subtracting \(B([\bar A(a)])\). Thus higher-depth tests for a fixed affine word are transported sums of exact branch defects. Cancellation of the sum need not be branchwise. A formal branch index still needs a dyadic guard to specify a genuine Collatz trajectory.

## 2. Complete odd-prime lifting on the rails

Factor
\[
 E_n-1=(t-1)(2t+1),\qquad E_n+1=Q(t)=2t^2-t+1,\qquad
 8Q(t)=(4t-1)^2+7.
\]
At any odd prime \(\ell\), write \(d=\operatorname{ord}_\ell(4)\), \(s=v_\ell(4^d-1)=v_\ell(2^{\ell-1}-1)\); the equality follows by LTE since \(d\mid\ell-1\) and \(4^{\ell-1}=(2^{\ell-1})^2\).

**Exponent ball.** For a fixed branch \(r\equiv a\bmod d\), put \(\rho=4^a\bmod\ell\). At every \(k\le s\), the image \(4^r\bmod\ell^k\) is the single value \([\rho]\bmod\ell^k\). At every \(k>s\) its images are exactly the residues congruent to \([\rho]\bmod\ell^s\), with each such residue corresponding to a unique exponent class \(r\bmod d\ell^{k-s}\). In fact \(4^a\equiv[\rho]\bmod\ell^s\) by uniqueness of the simple Hensel lift of \(X^{\ell-1}=1\). The order of \(4^d\bmod\ell^k\) is \(\ell^{k-s}\) for \(k>s\) by LTE, so its cyclic set fills the indicated ball by cardinality.

**Simple upper roots (\(\ell\ne7\)).** Only roots \(\rho\) in \(\langle4\rangle\subset\mathbb F_\ell^\times\) are relevant. Let \(\alpha\) be the unique Hensel lift of such a root, and put \(C(\rho)=[\rho]+[1-\rho]-1\). The residue identity \(1-\rho=-2\rho^2\), Teichmüller multiplicativity and \(2\equiv[2]\bmod\ell^s\) give
\[
 Q([\rho])\equiv-C(\rho)\pmod{\ell^s}.
\]
Because \(Q(u)/(u-\alpha)\) is a unit on the root class, \(v_\ell(\alpha-[\rho])\wedge s=v_\ell(C(\rho))\wedge s\). If \(j=v_\ell(C(\rho))<s\), the entire exponent branch has exact valuation \(j\) on \(E_n+1\). Otherwise the whole branch reaches \(s\), and for each \(k>s\) there is exactly one class \(r\bmod d\ell^{k-s}\) reaching \(k\); exactly one of its \(\ell\) children reaches \(k+1\). There is no truncation in \(k\).

Modulo \(\ell^2\), the numerator of \(C(\rho)\) is represented by \(x^\ell+(1-x)^\ell-1\) for any \(x\) representing \(\rho\). Thus \(C(\rho)/\ell\equiv H_\ell(-\rho)\bmod\ell\). **At a base-2 Wieferich prime, an accessible upper root lifts to \(\ell^2\) exactly when this Mirimanoff value vanishes.** Higher carry digits determine the remainder of the plateau. At an ordinary prime, \(s=1\) and every accessible simple root starts its unique infinite tower. Neither statement removes the subgroup condition.

If the roots split, their product is \(1/2\). With \(o=\operatorname{ord}_\ell(2)\), when \(o\) is odd, both roots are accessible or neither is; when \(o\) is even at most one is accessible. In the first case the two plateau depths agree: the Hensel roots have product \(1/2\), and the Teichmüller representatives have product \([1/2]\equiv1/2\bmod\ell^s\), so their distances agree in truncated valuation by division by units.

**Ramified upper root.** At \(\ell=7\) the root is \(t=2\bmod7\), equivalently \(r=2\bmod3\), or \(n=5\bmod6\). The completed square gives valuation exactly one throughout the branch, zero elsewhere. For all other primes, \(-7\) being a square is necessary but not sufficient: the subgroup test remains.

**Lower roots.** Apart from \(\ell=3\), \(1,-1/2\) are distinct roots. Both already agree with their Teichmüller lifts modulo \(\ell^s\), so accessible branches always reach \(s\). The root \(1\) occurs for \(d\mid r\) and has valuation \(s+v_\ell(r/d)\). The root \(-1/2\) is accessible iff \(o\equiv2\bmod4\), because \(2^n=-1\) in the residue field requires the odd exponent \(n\) to be an odd multiple of \(o/2\). On this branch, the odd-multiple version of LTE gives valuation \(v_\ell(2^n+1)=s+v_\ell(n/(o/2))\). At \(\ell=3\), apply LTE to each factor and add:
\[
 v_3(E_n-1)=2+v_3(r)+v_3(n),\qquad 3\nmid E_n+1.
\]
In particular at a prime exponent \(n=p\ge5\), \(v_p(E_p-1)=v_p(2^{p-1}-1)\). The lower rail is always composite, and its multiplicity at its own exponent is exactly that exponent's Wieferich depth. This is not a tournament obstruction.

## 3. The exact first-departure branches in integer Collatz

Let \(M=2^n-1\). In the accelerated odd map \(x\mapsto(3x+1)/2^{v_2(3x+1)}\), the lower rail shadows the branch of \(-1\) (guard 1) for exactly \(2r-1\) steps:
\[
 x_i=2^{2r-i}3^iM-1,\qquad 0\le i\le2r-1.
\]
Before the last displayed state \(3x_i+1\equiv-2\bmod4\); at the last the guard is exactly two because \(3^{2r}M-1\equiv-2\bmod8\). Its first departure, after \(2r\) accelerated steps, lands at \((3^{2r}M-1)/2\).

The upper rail shadows the branch of \(+1\) (guard 2) for exactly \(r-1\) steps:
\[
 y_i=2^{2r-2i}3^iM+1,\qquad 0\le i\le r-1.
\]
Before the last displayed state \(3y_i+1\equiv4\bmod16\); at the last the guard is \(2+u\), where
\[
 u=v_2(3^rM+1)=v_2(3^r-1)=
 \begin{cases}1&r\ \text{odd},\\2+v_2(r)&r\ \text{even}.\end{cases}
\]
The equality holds because \(3^rM+1=2^n3^r-(3^r-1)\) and the displayed LTE valuation is less than \(n=2r+1\). Its first departure after \(r\) steps lands at \((3^rM+1)/2^u\).

For fixed \(n\) and rail, both departure maps are affine in \(M\), with unit slopes at every odd \(\ell\ne3\). Thus a guarded first-departure state transports, and can invert, the corresponding odd-prime residue tower away from 3. This synchronizes departure formulas (step counts in ratio two); it asserts neither coalescence nor anything about eventual stopping.

The finite logarithm is classical: Elbaz-Vincent and Gangl, [*On Poly(ana)logs I*](https://arxiv.org/abs/math/0008089), and the prior [finite-logarithm bridge](../wronskian_frobenius_20261019/FINITE_LOGARITHM_BRIDGE.md) and [physical base-2 phase detector](../chanpak_perfect_wieferich_seams_20261122/PROOF.md#the-first-repeated-width-phase-detects-wieferich-primes). The present addition is the exact branch-specific lifting contract, the valuation behavior at every depth on both perfect rails, and the short guarded integer trajectories. None substitutes for an exact all-exponent tournament switch.

## 4. Explicit nonlifting and a new exact compositeness progression

At \(\ell=1093\), \(d=182\), \(s=2\), and the two upper roots are \(47,500\bmod1093\); just \(47=4^{74}\bmod1093\) is in \(\langle4\rangle\). The carry at 47 has \(C(47)/1093=389\bmod1093\), so
\[
 v_{1093}(E_n+1)=
 \begin{cases}
  1& n\equiv149\bmod364,\\
  0& \text{otherwise}.
 \end{cases}
\]
Thus this base-2 Wieferich prime gives a real plateau failure, rather than an empty conditional criterion. Both roots at 3511 (\(1815,3452\)) are outside \(\langle4\rangle\), so the subgroup gate there rejects both before the carry is tested. These numerical data are obtained by the finite powers and carry replay, without any assumption about how many Wieferich primes exist.
