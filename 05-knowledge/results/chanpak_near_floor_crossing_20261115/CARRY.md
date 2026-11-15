# Golden-ratio carry selection gives more than a two-thirds support exponent

For the physical order-five mask-zero carrier from the [sharp-cut packet](../chanpak_true_cut_fourier_20261115/README.md), let
\[
 F(x,y)=100xy+89x+11\{p(x)+p(y)+1+[x\text{ odd}]\}.
\]
There are \(\gg T^\kappa\) **distinct** values \(F(x,y)\le T\), where
\[
 \kappa=\frac{1+\log_2(1+\sqrt5)}4=0.6735604784\ldots>\frac23.
\]
More strongly, the construction supplies that order of magnitude already with both addresses comparable to \(\sqrt T\) and all its values in \([cT,T]\) for an absolute \(c>0\). The statement is not density one. The exact symbolic growth comes from carry-state selection, not from extrapolating the finite carrier census.

## A physical injective selection from each sum line

Take \(V=2^r\), \(r\ge3\), and \(V\le u<2V\). For \(0\le w<V/8\) write \(s=V/2+w\). Use addresses
\[
 x=Vu+b,\qquad y=Vu+c,\qquad b+c=s,\qquad 0\le b\le w,\quad c=V/2+(w-b).
\]
In the low \(\ell=r-3\) digits of \(b,c\), choose one common carry path for the sum \(w\), with carry zero both in and out. The sum of the two digits at position \(i\) is then fixed: \(d_i=w_i+2C_{i+1}-C_i\). The remaining three digits are fixed as well. In particular \(\operatorname{rev}_r(b)+\operatorname{rev}_r(c)\) is constant within the carry path. By the exact binary split the large \(u\)-dependent part of \(p(Vu+b)+p(Vu+c)\) is constant within this path; the additional parity boundary term is also constant because \(b<V/2\), \(c\ge V/2\).

At fixed \(u,s\), the product increment from \(b_1<b_2\) is \((b_2-b_1)(s-b_1-b_2)\ge V(b_2-b_1)/4\). Multiplication by 100 gives at least \(25V\). The remaining variable checksum contribution has range at most \(11s+11<25V\), by \(0\le p(i)\le i\), and \(89x\) is increasing. Hence the selected outputs are strictly increasing with \(b\). Each free digit \(d_i=1\) supplies two literal choices, and every other digit supplies just one. Thus the largest one-carry-path selection has size \(W_\ell(w)\), the maximum \(2^{\#\{i:d_i=1\}}\) over legal carry paths.

Different lines also have disjoint output intervals. Put \(S=x+y=2Vu+s\) and \(d=y-x\), with \(0<d<V\). Relative to \(25S^2+(89/2)S\), the value has the correction
\[
 -25d^2-(89/2)d+11\{p(x)+p(y)+1+[x\text{ odd}]\}.
\]
Its lower bound is \(-25V^2-(89/2)V\), and its upper bound is \(11S+22\). For successive integer \(S\ge2V^2\) the centers are separated by at least \(50S\), more than these two errors combined when \(V\ge8\). All integer sums arising from the distinct choices \((u,w)\) are distinct. The values are between \(100V^4\) and \(500V^4\) for \(V\ge8\), since both addresses are in \([V^2,2V^2)\). Thus the constructed count of different values is exactly \(V\sum_{w<2^\ell}W_\ell(w)\), with an upper cutoff \(500V^4\).

## The two-state growth that appears after optimizing a carry

Process bits of \(w\) from least to most significant, keeping \((A,B)\) as the maximum path weights ending in carry 0 or 1. Initialize \((A,B)=(1,0)\). The two updates are
\[
0:(A,B)\mapsto(A,\max(A,2B)),\qquad
1:(A,B)\mapsto(\max(2A,B),B),
\]
and the final weight is \(A=W_\ell(w)\). This is a maximum, not a sum: the choices from distinct carry paths on the same line were not certified disjoint.

There are two equivalent cuts through this weighted automaton. For a direct growth proof, require the first bit of \(w\) to be zero, so the state is \((1,1)\); put \(m=\max(A,B)\). An equal state has two children, both unequal and both of weight \(2m\). An unequal state has one unequal child of weight \(2m\) and one equal child of weight \(m\). This uses that \(A,B\) are powers of two, so the larger is at least double the smaller. On the total weights of equal/unequal states the transition is therefore
\[
 \begin{pmatrix}0&1\\4&2\end{pmatrix}.
\]
After the initial zero, exchanging bits \(0,1\) interchanges \(A,B\), so the sum of \(A\)'s is between one half and all of the summed \(m\). This proves the lower bound with growth \(1+\sqrt5\). For the upper bound, allow an initial run of ones before the first zero: a run of length \(k\) weighs \(2^k\); summing \(2^k(1+\sqrt5)^{\ell-k-1}\) gives the same order, and the all-ones word weighs only \(2^\ell\).

In fact, letting \(z\) mark the bit length, the exact sum of \(A\)'s has generating function
\[
 \frac1{1-2z}+
 \frac{z}{1-2z}\,
 \frac{1/(1-z)+2z/(1-2z)}{1-4z^2/(1-2z)}.
\]
To obtain it, use the difference \(d_0=\log_2(B/A)\). After reaching equality, first returns to equality are \(0^j1\) and \(1^j0\) with \(j\ge1\), both multiplying \(A\) by \(2^j\). A final nonreturn can instead be all zeros (no change in \(A\)), all ones (weight \(2^j\)), or empty. The initial approach to equality is \(1^k0\), of weight \(2^k\). This proves the formula, also checking that growth comes from optimizing one carry path and not silently aggregating outputs across carry classes.

Combining \(V=2^r\) and the growth at \(\ell=r-3\) yields \(\gg V(1+\sqrt5)^r\) distinct counts. Choose the largest dyadic \(V\) with \(500V^4\le T\) to obtain the stated all-\(T\) estimate. The [replay](verify_carry.py) tests the optimal path count and injection against the literal carrier formula. This supplies a rigorous lower bound on support in a balanced region. It does not control overlaps between the remaining carry paths, or the unchanged second-moment and absolute-target obligations.
