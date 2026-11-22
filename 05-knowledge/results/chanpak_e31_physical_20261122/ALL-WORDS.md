# A complete inverse and an upper density bound for the two binary flanks

**PROVED.** The six-vertex construction giving a single arc across \(E_{31}\) admits a complete inverse on its two independent words. For prescribed odd lower count \(N\), all solutions \((v,y)\), \(v\ge1,y\ge0\), can be enumerated in \(O((1+N)^{\log_3 2}\log(2+N))\) arithmetic operations on integers of \(O(\log(2+N))\) bits. The family, including **both** endpoints, has upper relative density below 0.133 among the positive odd integers. Hence other marked cores or other constructions are necessary for Chan–Pak. These bounds concern the binary receiver family in the [exponent-31 proof](README.md); the general unit-transmitting receivers classified there form a larger collection.

## Charging three choices in the middle of a receiver

Write \(b\) for the width of a positive binary receiver address \(v\). The receiver tournament consists of a transitive spine and its tagged apex. Its two-slot count \(Q_v\) counts all assignments of the spine vertices to the part before that apex in its slot, the part after it, or the other slot. Within each part the order is forced. The last bit before the apex must be zero when that part is nonempty; the first bit after the apex must be one when that part is nonempty. Choose the apex slot in two ways.

If \(u\) is positive and even of width \(b\), the leading bit is one and the last bit zero. Put the leading vertex first after the apex and the final vertex last before it. Each of the other \(b-2\) vertices can be assigned independently to any of the three parts. Thus
\[
            2\cdot3^{b-2}\le Q_u\le2\cdot3^b .
\]
Appending a one to the word adds an ordinal sink to the receiver tournament, doubling its two-slot count and preserving its Hamilton count. Every positive word has the unique decomposition
\[
               v=(u+1)2^t-1
\]
with \(u\) positive even, \(t\ge0\), or \(u=0,t\ge1\) for a word of all ones. In the latter case use the single-vertex receiver at \(u=0\) as the seed. In all cases \(Q_v=2^t Q_u\), \(h_v=1+f(u)\), and \(m_v=m_u\) for the dual endpoint mass \(m_v=2^{|v|}-\operatorname{rev}_{|v|}(v)\); use \(m_0=1,Q_0=2\). The last assertion follows by reversing the appended string of ones. Here \(f\) is the binary read called \(p\) in the main proof.

For a fixed right parity \(j\), write the exact lower count as
\[
 L(v,y)=K_{v,j}+C_v y+D_v f(y),\quad C_v=25(v+1)+28Q_v=2^t C_u,\quad D_v=25h_v.
\]
The intercept is \(K_{v,j}=D_v+C_v-28m_u-10-25[v\text{ even and }j=0]\). In particular \(D_v\le C_v\). Since \(Q_v\ge2^{|v|+1}\) from the two-slot subset formula and \(m_v\le2^{|v|}\), the exact word formula gives \(L+10\ge14 Q_v(y+1)\). It also gives \(L+10\ge53(v+1)(y+1)\). The latter follows from \(Q_v(y+1)-m_v\ge 2^{|v|}(2y+1)\) and \(v+1\le2^{|v|}\). All intercepts are positive for positive \(v\).

## Target-first enumeration

Given \(N\), enumerate the all-ones words under \(Q_v\le(N+10)/14\). For every remaining width \(b\) and trailing length \(t\), the necessary budget is \(2^t 3^b\le9(N+10)/28\). For each even \(b\)-bit prefix compute \(Q_u\) by the four-state recurrence, reject if it exceeds the budget, and enumerate its eligible tail lengths. At each retained \(v\), solve the two parity equations using the exact binary suffix inverse. It produces every positive \(y\) in that parity. Test the empty address \(y=0\) separately. [Exact one-dimensional inverse](../chanpak_gap_star_ports_20261122/README.md#one-arm-target-inversion-and-two-exact-landings)

For clarity, splitting \(y=Bz+r\) at a power of two \(B\asymp\max(1,\sqrt{N/C_v})\) costs \(O(1+\sqrt{N/C_v})\) residue and candidate tests. For fixed \(r\) the leading coefficient of \(z\) is at least \(BC_v\). The residual from \(f(z)\) and the boundary bit lies in \([0,D_v(z+1)]\), giving an interval of bounded length because \(D_v\le C_v\) and \(z\le N/(BC_v)\). The implementation uses the equivalent gcd-sieved split. Each test and the prefix calculation cost at most \(O(\log(2+N))\) arithmetic operations. Put \(B_0=\lfloor\log_3(9(N+10)/28)\rfloor\). Summing over the tail lengths and prefixes costs, before the logarithmic factor,
\[
 O\!\left(\sqrt N+\sum_{b=2}^{B_0}2^b\left(1+B_0-b+\sqrt{N/3^b}\right)\right)
       =O((1+N)^{\log_3 2}),
\]
since \(2/\sqrt3>1\). For the all-ones branch the costs form a geometric sum plus \(O(\log N)\).

There are in fact \(\Theta(X^{\log_3 2})\) distinct receiver words whose smallest lower count is at most \(X\). The budget gives the upper bound by the same sum. Conversely at \(y=0\) the upper bounds \(Q_u\le2\cdot3^b\) and \(1+f(u)\le u+1\le2^b\) give \(L(u,0)\le106\cdot3^b\) for every even word of width \(b\). This counts receiver words and does not assert a lower time bound for other algorithms.

**FINITE-EXACT.** For each integer exponent from 2 through 18 inclusive the full two-word inverse is empty at the prescribed lower endpoint \(E_r-1\). This includes exponent 17, where a different core is already known to work. The same inverse equation has explicit solutions at 23 and 31. To reproduce the finite exclusion from this directory, run the following.

```sh
for r in {2..18}
do
    python3 -B all_words.py --exponent "$r"
done
```

## Why this family cannot be an absolute selector

For each \(v\) there are at most \(X/C_v+1\) possible \(y\) with \(L(v,y)\le X\). Since the number of active receivers is \(O(X^{\log_3 2})=o(X)\), the upper relative density of the **union of both endpoints** within the odd integers is bounded by \(4S\), where
\[
 S=\sum_{v\ge1}\frac1{C_v}
   =\frac1{81}+2\sum_{b=2}^{\infty}\ 
                  \sum_{\substack{u\text{ even}\\|u|=b}}
             \frac1{25(u+1)+28 Q_u}.
\]
Here \(1/81\) sums the all-ones branch and the factor two sums \(2^{-t}\) at every nonzero even prefix. Truncating after width \(B\) leaves a positive tail of at most \(\frac3{28}(2/3)^{B-1}\): there are \(2^{b-2}\) even words at width \(b\), and each satisfies the three-choice lower bound. Even the estimate \(S\le1/81+3/28<1/8\) gives density below one half. At \(B=10\), rounding each included term \(2/C_u\) upward to the next multiple of \(10^{-9}\) and adding the rigorous tail gives \(4S<0.133\). This counts all unbounded word lengths; overlaps between different receivers can only decrease the density.

Run `python3 -B all_words.py` in this directory for the upper sum and complete small inverse comparisons. Use `--lower N` to search both words for a prescribed lower count, or `--exponent r` for the center \(E_r\). The optional `--max-width b` restricts the normalized even prefix and leaves all admissible trailing lengths free. For the large confirmed center at 31 the short receiver is already known, and `python3 -B check.py --exponent 31 --reader 6` returns the right address directly.
