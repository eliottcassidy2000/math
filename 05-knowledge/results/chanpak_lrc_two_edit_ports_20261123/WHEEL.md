# A finite lonely-runner wheel for the physical factor-seven fork

Write $D=H(M)$, $Q=Q_2(M)$ for a feasible filler. The actual [factor-seven fork and edit route](../chanpak_cofactor_contradiction_20261123/FILLERS.md) have parent $h=7(40Q+D)$, center $3h$, and marked exits $840Q+19D,840Q+23D$. This note uses the forbidden-set covering step of Lonely Runner on finite residue circles. It does not invoke a continuous Lonely Runner bound to claim an exact prime.

## Exact three-form wheel

At each odd prime $\ell$, set $D=1$ in that component and choose the common clock $Q\in\mathbb F_\ell$. The three affine forms
$$
40Q+1,\qquad 840Q+19,\qquad 840Q+23
$$
must avoid zero. At $\ell\ge11$ their slopes are units and their zero points are distinct: the determinants of the three coefficient pairs are $-80,80,3360$, with prime factors only $2,3,5,7$. Thus there are exactly $\ell-3$ permissible clock values. This is the discrete avoided-arc argument, with each forbidden arc a singleton; the exact local density is $1-3/\ell$. At seven only the first form has a forbidden zero, at five there are none, and at three the even-order condition $Q=D=1\bmod3$ makes all forms nonzero. (The cofactor is instead forced to be zero mod three in the odd-order phase.) The parity at two is automatically safe.

The [actual joint odd-prime-power controller](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md) realizes the choices simultaneously instead of merely choosing unrelated residues for stand-alone polynomials. The [height-balanced selector](../chanpak_heegner_absolute_targets_20261123/BALANCED-WINDOWS.md) additionally places the physical fork at any sufficiently high scale. For every fixed wheel, we therefore obtain a literal fork whose cofactor and both marked exits are simultaneously free of the wheel primes.

## The wheel extends to all seven reads of the actual edit route

The physical five-edit route has heights
$$
280Q+7D,\quad780Q+45D,\quad468Q+45D,\quad312Q+67D,\quad840Q+19D,\quad840Q+21D,
$$
and a second marked exit at $840Q+23D$. Strip the fixed scalar factors, retaining the six distinct forms
$$
40Q+D,\quad52Q+3D,\quad(52Q+5D)/3,\quad312Q+67D,\quad840Q+19D,\quad840Q+23D.
$$
Here even filler order is required: $Q\equiv D\bmod3$, so the displayed quotient is an integer. The corresponding fixed factors along the route are $7,15,27,1,1,21$, and 1 for the other exit. No claim is made that these multipliers are universal outside the selected even-order family.

At $\ell\ge7$, each form forbids at most one point after setting $D=1$. A zero slope has a nonzero constant, as is apparent from the coprime slope/constant pairs (the division by three is a unit there). At five the union of forbidden points is exactly $\{0,1,4\}$, so take $Q=2$. Modulo nine use $D=Q=1$, which makes the third form congruent to $(52+5)/3=19\bmod3$ and the other five forms nonzero mod three. Both residues are compatible with the prescribed even order. At two all normalized reads are odd. Thus the scheme applies to all normalized vertices at once, not just the marked exits, while keeping the underlying literal edit route and its height ordering band.

## A growing rigorous wheel in every large scale window

Call a positive integer $y$-rough if it has no prime divisor at most $y$. For real $X$ tending to infinity, set
$$
y=\left\lfloor\frac{\log\log X}{32\log2}\right\rfloor,\qquad
M=9\prod_{5\le p\le y}p .
$$
Then for every sufficiently large $X$ there is a literal fork with parent height
$$
h=X+O\bigl(X^{1-\frac12\log_3(3/2)+o(1)}\bigr)
$$
and all six displayed normalized route reads $y$-rough. The use of a growing modulus is justified by a uniform bound on the literal correction length, as follows.

Use the previously constructive simultaneous controller, with guard $g=2$ because the three-primary component is 9. Each macro expands to at most five contacts. Its common length modulus satisfies $L\le M$ and its inverse exponent bound satisfies $E\le M^2$. The length-class base $\mathcal O^r$ uses at most $2M$ signed macro letters. At a prime $p\ge5$, the shared-conjugation projection used in the earlier proof has at most $4(p^2-1)$ signed letters, a commutator at most $16(p^2-1)$, and at most $p/2$ such directions are used; the pure center is constant size. The three-primary component uses two constant-size special generators. Joint CRT exponents can each be chosen between zero and $M-1$, so the total signed length is $O(M(1+\sum_{p\le y}p^3))=O(My^4)$ once $y\ge5$. Replacing inverse occurrences by positive powers costs at most $E$, and a nonempty identity can be supplied by one macro raised to $E$. Thus the uniform correction bound is $B=O(M^3y^4)$, independent of the chosen suffix, parity-preserving length class, and target values.

To bound $M$ elementarily, $\prod_{m<p\le2m}p\le\binom{2m}{m}\le4^m$; sum over dyadic intervals up to the next power of two above $y$. This gives $\prod_{p\le y}p\le e^{4y\log2}$. Therefore $B=O((\log X)^{3/8}(\log\log X)^4)=o(\log X)$. In the base-$3/2$ native height construction, choose even tournament order with chain length $n\asymp\log_3 X$ (the two real normalized bands tile at either fixed parity), and divide the variable prefix and suffix widths into halves after budgeting $B$. Both widths are at least $n/2-o(n)$. The explicit geometric error of that construction, $O((2/3)^{\min(a,k)}+n(2/3)^n)$ after scaling by $3^n$, is therefore $X^{-\frac12\log_3(3/2)+o(1)}$. This proves the stated absolute error without interchanging a fixed-modulus limit with a growing modulus.

No wheel depth here certifies that a normalized value is prime: $y=O(\log\log X)$ is far below $\sqrt X$, and the selector gives a sublinear height window, not exact equality with $7p$. The useful transfer is the physically simultaneous modular avoidance and its scale control, not a universal prime-to-own-sevenfold flip.

## Third-slot follow-up

The [filler-circle continuation](FILLER-CIRCLES.md) applies same-clock avoidance to a varying third slot while the gain Jacobian stays fixed; the [CRT companion](FILLER-CRT.md) certifies jointly reachable unit phases instead of independent local choices.
