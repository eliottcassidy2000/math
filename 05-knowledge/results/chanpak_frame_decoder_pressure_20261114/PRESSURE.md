# Why the terrace exponent was 0.98, and how far to trust it

**Proved.** If \(V(X)\) counts the distinct Hamiltonian-path totals
at most \(X\) on the minimum directed-triangle terrace of strong
tournaments, then
\[
       V(X)=\Omega(X^{11/30}) \quad\text{and}\quad V(X)=O(X^{119/125}).
\]
The upper exponent is \(0.952\); the lower exponent is about \(0.367\).
The connected **ordered gap words**, before numerical coincidences
are removed, have a true Dirichlet threshold in \((0.950,0.952)\).
That narrow bracket is not a bracket for the distinct-value exponent.

The argument uses the corridor kernels from
[the terrace factorization](https://github.com/eliottcassidy2000/math/blob/8051b02275faab24bffba7fbed5cc8d3afc29bac/05-knowledge/results/chanpak_target_first_terrace_density_20261114/TERRACE_DENSITY.md).
The finite rational inequalities below are checked by running
pressure_replay.py with Python 3.

## 1. An exact multiplicative cost along a connected corridor

A connected word has positive gap widths \(r_i\ge1\) with \(q_i\ge1\)
zero gaps between successive positive gaps. Its initial two-state vector
is \((a,b)=(1,2^{r_1})\). Across \(q\) zeros and the next positive width
\(r\), put \(y=2^r\). The new state is
\[
  (a',b')=(A_q a+B_q b,\ y(B_q a+D_q b)).
\]
The scalar readout at every positive gap is \(a+b\). The three kernel
sequences start respectively at
\((1,1,3)\), \((1,2,3)\), \((1,2,4)\), and continue by tribonacci.
They satisfy \(1\le A_q\le B_q\le D_q\). This follows, for example,
from the nonnegative coefficients of the generating function
differences \(z(1-z)/(1-z-z^2-z^3)\) and
\(z^2/(1-z-z^2-z^3)\).

Normalize by \(p=b/(a+b)\in[2/3,1]\). This interval is invariant since
\(b'/a'\ge y\ge2\). The exact growth factor and new state are
\[
\begin{aligned}
 G_{q,r}(p)&=A_q+(B_q-A_q)p
                 +2^r\bigl[B_q+(D_q-B_q)p\bigr],\\
 P_{q,r}(p)&=\frac{2^r[B_q+(D_q-B_q)p]}{G_{q,r}(p)} .
\end{aligned}
\]
The first function is nondecreasing in \(p\); the second is monotone
on each interval (not necessarily increasing: it is linear-fractional).
The older \(0.98\) bound treated each corridor via a uniform
Bernoulli-\(2/3\) estimate and then used a convenient geometric-series
convergence test. It was not a fitted exponent, nor a distinguished
constant of tournaments. Keeping the actual evolving \(p\) retains
the consecutive-corridor dependence. The [exact reset analysis](RESETS_AND_DECODER.md#exact-independence-versus-fast-loss-of-information) gives a uniform log-ratio contraction below 0.172; it explains why a coarse projective mesh is effective without making exact integer information disappear.

## 2. Certified weighted branching bounds

For \(s>0\), define
\[
 (\mathcal L_s f)(p)=\sum_{q,r\ge1}
      G_{q,r}(p)^{-s} f(P_{q,r}(p)).
\]
A word's inverse-power weight factors into its initial count raised to
\(-s\) and these exact branch weights. Use the positive test function
\(f(p)=1-\frac{137}{1000}(p-\frac23)\). This choice need not be optimal.
The replay proves the uniform upper inequality
\[
     \mathcal L_{119/125} f < (999/1000)f
\]
and, for the operator truncated to \(q,r\le16\), the lower inequality
\[
     \mathcal L^{[16]}_{19/20} f > (1003/1000)f .
\]
Here is how to reproduce the certificates without floating assumptions.
Divide \([2/3,1]\) into 80 equal closed intervals. On \([l,u]\) an upper
term bound is \(G(l)^{-s} f(\min(P(l),P(u)))\) and a lower term bound is
\(G(u)^{-s} f(\max(P(l),P(u)))\). Divide summed upper bounds by \(f(u)\)
and lower bounds by \(f(l)\). To bound a rational \(x^{-s}\) for
\(x=a/b\) and \(s=n/d\), compare \(k^d a^n\) with \(10^{8d}b^n\)
using integers; round the bounds outward to a grid of size \(10^{-8}\).
The upper certificate includes the omitted terms, using
\(G_{q,r}(p)\ge (1+2^r)A_q\), \(f\le1\), and
\(A_q\ge(9/5)A_{q-1}\) for \(q\ge4\) (check \(q=4,5,6\), then
apply tribonacci to the differences). The two geometric tail bounds
are combined as \(R_{\rm tail}(Q_{\rm fin}+Q_{\rm tail})
+R_{\rm fin}Q_{\rm tail}\). The lower certificate ignores the omitted
positive terms.

Iterate the upper inequality. Since \(f\) is bounded away from zero,
it bounds the total of \(H^{-119/125}\) over all ordered connected
words; the sum of initial weights \(\sum_{r\ge1}(1+2^r)^{-119/125}\)
is finite. Cutting at the zero-length seams expresses every terrace
total as a product of connected totals. For distinct connected
values \(m\ge3\), the product \(\prod_m(1-m^{-119/125})^{-1}\)
converges by the connected sum bound. This product sums multisets,
so it bounds the inverse-power series for distinct terrace totals
even when different multisets give the same integer. Each total at
most \(X\) contributes at least \(X^{-119/125}\), proving the upper
bound. This Euler product does not assert unique factorization of
the integer into corridor totals.

### A genuine pressure, but for words

For each length \(k\), let \(Z_{k,s}(p)=\mathcal L_s^k 1(p)\).
Fix a continuation of positive length; its final count starting from
state \((1-p,p)\) is \(U(1-p)+Vp\) with \(V\ge U>0\).
Indeed the two columns of the first continuation matrix are ordered
entrywise, and later matrices are nonnegative. Its maximum on our
interval is therefore at most \(3/2\) times its minimum, regardless
of length. The inverse-power sums at distinct starting points are
comparable by \((3/2)^s\). Also \(\sup_p Z_{k,s}(p)\) is
submultiplicative in \(k\), and \(\inf_p Z_{k,s}(p)\) is
supermultiplicative. Their positive \(k\)-th roots have the same
limit \(\lambda(s)\), the per-corridor weighted branching rate.
For every \(s>0\) it is finite by the geometric bounds and positive
by using one fixed repeated letter. Hölder's inequality makes
\(\log\lambda\) convex, hence continuous; \(G\ge3\) yields
\(\lambda(s+\delta)\le3^{-\delta}\lambda(s)\) for \(\delta>0\).
The two certificates force exactly one root of \(\lambda(s)=1\)
in \((0.950,0.952)\). At that root, submultiplicativity implies
\(\sup_p Z_{k,s}(p)\ge1\) for every \(k\), so uniform comparability
implies that the length sums still diverge. For larger \(s\) they
converge, and for smaller \(s\) they diverge. The pressure is a
balance between branching entropy and payment in log path count.
Numerical collisions need not respect that balance, which is why a
separate injection is necessary for any lower bound on \(V(X)\).

## 3. Lower bound: combine ternary guard with binary decoding

Start at \((a,b)=(1,2)\). Allow \(1\le r\le20\), with two corridor
lengths and respective maps
\[
\begin{array}{c|c|l}
q&(a',b')&\text{allowed when}\\
1&(a+2b,\ 2^{r+1}(a+b))&a+2b\not\equiv0\pmod3\\
2&(3(a+b),\ 2^r(3a+4b))&b\not\equiv0\pmod3 .
\end{array}
\]
The pair remains odd-even and nonzero modulo \(3\). The last corridor
has length \(2\) precisely when \(3\mid a'\). For length \(1\),
\(r=v_2(b')-1\) and the source pair is
\((2z-a',a'-z)\), \(z=b'/2^{r+1}\). For length \(2\),
\(r=v_2(b')\) and the pair is \((4u-z,z-3u)\),
\(u=a'/3,\ z=b'/2^r\). The known seed is detected as the literal pair
\((1,2)\); growth strictly increases the sum at each step.
Thus the final state encodes the complete ordered history.
Neither a 2-adic nor a 3-adic readout alone is used to claim this;
the ternary guard recovers the corridor, and the binary valuation
recovers the positive gap.

The residue state \(t=a/b\in\mathbb P^1(\mathbf F_3)\) has four
possibilities. For \(q=1\), an even or odd width sends
\(0\to1\) or \(0\to2\), respectively,
\(2\to\infty\) in either case,
and \(\infty\to2\) or \(\infty\to1\), respectively.
State \(1\) cannot use \(q=1\). A \(q=2\) step sends any finite
state to \(0\) and is disallowed at \(\infty\).

For upper growth costs use \(c_{1,r}=2+2^{r+1}\) and
\(c_{2,r}=3+4\cdot2^r\). At \(\alpha=11/15\), let \(E,O\)
be the sums of \(c_{1,r}^{-\alpha}\) for even and odd widths,
and \(J\) the sum of \(c_{2,r}^{-\alpha}\).
The resulting weighted transition matrix, in state order
\((0,1,2,\infty)\), is
\[
 M=\begin{pmatrix}
 J&E&O&0 \\ J&0&0&0\\ J&0&0&E+O \\ 0&O&E&0
 \end{pmatrix}.
\]
The replay verifies \(Mz>(1006/1000)z\) componentwise for
\(z=(1000,475,832,469)^{\mathsf T}\). Assign a transition from
\(t\) to \(t'\) with cost \(c\) probability proportional to
\(c^{-\alpha}z_{t'}\), normalized in row \(t\). Its probability is
at most \(c^{-\alpha}z_{t'}/z_t\).

For any threshold \(Y>1\), stop each history at the first prefix
whose cost product is at least \(Y\). This is a finite cut of
probability one, each leaf has cost at most \(CY\) with
\(C=\max_{q,r} c_{q,r}\), and the probability of each leaf is at
most \((\max z/\min z)Y^{-\alpha}\) by telescoping. There
are therefore \(\Omega(Y^\alpha)\) leaves with \(a+b\le3CY\).

Choose one common power of two \(P\) with \(6CY<P\le12CY\)
and append a terminal \(q=1\) corridor of width \(r\) with
\(2^{r+1}=P\). The residue guard is not imposed on this known
terminal step. Its ordinary Hamiltonian path count is
\[
                    (a+2b)+P(a+b).
\]
The low block \(a+2b<P\) and the high block \(a+b\) uniquely recover
\((a,b)\). Those pairs encode distinct histories, so these ordinary
totals are all distinct. They are bounded by \(6CY+36C^2Y^2\);
choosing \(Y\) proportional to \(\sqrt X\) proves
\(\Omega(X^{11/30})\). This is also a target-first decoder on the
fixed-marker subfamily: split the target at the prescribed marker,
recover both coordinates, and peel by the two residue rules until
reaching the seed or rejecting. The marker is an exact but costly
lift from the physical two-state interface to a scalar. Its
quadratic overhead is one remaining gap in the exponent bounds.

## 4. What selection can and cannot do

A representative can be chosen unambiguously **inside this family**
by decoding, normalizing the gap word under reversal if desired,
and using a fixed total order on the remaining candidate classes.
An exact canonicalizer controls collisions but cannot fill absent
targets: even the entire terrace has density zero. The split between
word pressure and value growth is precisely the distinction between
building many patterns and selecting a *physical pattern with the
prescribed absolute total*. The marker illustrates a useful design
principle for other families: choose a construction that makes a
hidden vector decodable, then prove a separate covering statement
for the image of that decoder.
