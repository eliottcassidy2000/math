# A clock of order 24 in the actual two-flank reader

**PROVED.** Let \(H_n\) be the canonical lower height at spine \(n=4,8\). Let \(N\) be odd. Over the full padded rectangles \(0\le x<2^m,\ 0\le y<2^l\), the number with \(H_n(x,y)\equiv N\pmod4\) is exactly
\[
 \frac{2^{m+l}+(-1)^{(N-k_n)/2}S_m S_l}{2}.
\]
The sequence satisfies \(S_{j+12}=-64S_j\), with initial values at \(j=0,\ldots,11\)
\[
 (1,0,-2,0,0,-8,-8,0,-16,-32,0,0).
\]
In particular, a free 61-bit address supplies **exactly half** the desired odd residue, whatever the other address is. If both widths are 61, exactly \(2^{121}\) ordered pairs pass the E61 test on each spine. This count is over the full word rectangle, before imposing any bound on the absolute height.

## The finite transfer

The chart coefficients give \(4\mid q,t,e\), \(c\equiv r\equiv2\pmod4\), and \(h\) odd. Thus the nonconstant part of \((H_n-k_n)/2\bmod 2\) is the sum of the independent characters \(F(x)+F(y)\), where \(F(x)=f(x)/2+x_0\bmod2\), with \(x_0=x\bmod2\). From the [binary inversion weights](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md), a 1 at position \(i\) above a 0 at \(j\) contributes \(2^{i-j-1}+[i=j+1]\) to \(f\). After dividing by two, exactly distances 1 and 2 survive modulo two.

Read the bits from least to most significant. The first contributes its value to \(F\), the second contributes \(b_1(1-b_0)\). Every subsequent bit \(v\), with two preceding bits \(a,b\), contributes \(v(a+b)\) in \(\mathbb F_2\). Keeping the signed sums for terminal pairs \(00,01,10,11\) gives the transfer
\[
 M=\begin{pmatrix}1&0&1&0\\1&0&-1&0\\0&1&0&1\\0&-1&0&1\end{pmatrix},\qquad
 M^6=8\begin{pmatrix}0&0&0&1\\0&0&1&0\\0&-1&0&0\\-1&0&0&0\end{pmatrix}.
\]
In particular \(M^{12}=-64I\). The initial two-bit vector is \((1,-1,-1,-1)\). Summing its coordinates after each power gives \(S_j=\sum_{x<2^j}(-1)^{F(x)}\) for \(j\ge2\). Checking \(j=0,1,12,13\) extends the recurrence to the two short words. Character orthogonality gives the displayed two-flank count. The normalized matrix \(M/\sqrt2\) has order exactly 24. Indeed its twelfth power is \(-I\), and its eighth power is not the identity, as seen by multiplying \(M^6\) by \(M^2\).

## What remains visible after both ends of a word are fixed

Leave any consecutive \(s\ge1\) bits free in one word and fix all exterior bits and the opposite word. The same mod-four test is a quadratic polynomial over \(\mathbb F_2^s\). Interactions across either boundary contribute only linear and constant terms. Its polar form is the binary adjacency matrix of the square of the \(s\)-vertex path, with an edge at distance 1 or 2. The radical dimension, in order \(s=0,1,2,3,4,5\pmod6\), is
\[
                   0,1,0,1,2,1 .
\]
To see this, impose zero boundary coordinates \(a_{-2}=a_{-1}=a_s=a_{s+1}=0\). The radical equation is \(a_{i-2}+a_{i-1}+a_{i+1}+a_{i+2}=0\). Its characteristic polynomial is \((t+1)^2(t^2+t+1)\), a divisor of \(t^6-1\) over \(\mathbb F_2\). Starting with \((a_0,a_1)=(u,v)\), the six entries repeat as \((u,v,v,u+v,0,0)\). Imposing \(a_s=a_{s+1}=0\) on that list gives the stated dimensions, including the one-bit case.

For any quadratic phase \(Q\) with radical dimension \(r\), the signed sum has magnitude either zero or \(2^{(s+r)/2}\). Indeed, square the sum and change variables from two inputs to one input and their difference. The inner sum vanishes off the radical, and on it the phase difference is a linear character of the radical. Thus fixing more exterior data cannot make the mod-four sieve exclude every filling of any free block of two or more bits. Its possible bias from half is explicitly bounded, even before evaluating the actual boundary signs.

This period 24 belongs to a signed transfer that counts actual address assignments. The [Ogg-prime connection](../chanpak_profinite_automata_20261122/README.md) concerns characteristics of supersingular elliptic curves and prime divisors of a group order. The period computation supplies no group action of the Monster or Baby Monster on these words.

Run `python3 -B digital_character.py` for the matrix identity, the polar ranks, the one-word sequence against direct integer reads, and small rectangles against both physical chart formulas.
