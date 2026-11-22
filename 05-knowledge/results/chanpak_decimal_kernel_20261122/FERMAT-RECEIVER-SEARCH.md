# How far can the canonical Farey receiver reach at exponent 61?

Let \(n\) be a positive multiple of four and let \(V_n=C_3[\mathrm{TT}_n,1,1]\) with spine \(\to a\to b\to\) spine. Tag \(a,b\). In spine order the binary expansion, with \(n\) digits, of
\[
 r_n=\frac{3(2^n+1)-1}{5}
\]
marks the **nontags** by ones. Give the four-vertex left block tag mask 12 and the six-vertex middle block tag mask 35 from the [Farey certificate](E15-FAREY.md). Every pair of exterior addresses \(x,y\ge0\) gives an arc of size two. Write \(N=2^{60}(2^{61}-1)-1\) for the prescribed lower endpoint.

**PROVED, using finite exact computation.**

- At every \(n\), if each address is either empty or has just one run of ones, no address pair attains \(N\). Either run may end in arbitrarily many zeros, and the two widths and both shifts are independent.
- For every \(n\ge36\), the exclusion holds for **arbitrary** pairs of binary words. For \(n\ge28\) it holds on each axis, where only one flank is nonempty. At \(n=12\) it also holds on the right-only axis.
- At \(n=72\) the lower count with empty flanks already equals \(115471007792477950669838719687691272735\), exceeding \(N\). This excludes all larger levels. Indeed \(r_{n+4}=9\cdot2^n+r_n\), so the larger receiver contains the smaller with exactly its tags. Its lower host with empty flanks therefore contains the earlier lower host. Any Hamiltonian path in an induced subtournament can be extended by each new vertex, by inserting at an endpoint or at a consecutive pair of the path. Fixing one such insertion per path gives an injection. Nonempty exterior flanks also contain the empty-flank host.

These conclusions hold at every total order in this particular receiver family. The unresolved levels with unrestricted flanks are \(4,8,12,16,20,24,28,32\). The exact bound below makes additional levels individually decidable with more computation.

## Seven readings in time linear in the spine length

A Hamiltonian path of this graph is a permutation of the special vertices interleaved with increasing runs of spine vertices. For \(h,e,r\), use special vertices \(a,b\), imposing the required endpoint tag on the last or first vertex as appropriate. For \(k\), add the physical apex. For \(q,c,t\), introduce a separator joined in both directions to every ordinary vertex. Cutting a path at this separator bijects it with an ordered two-slot cover, including empty slots. Requiring a tagged first vertex computes \(c\). Requiring an untagged last vertex computes \(t\) by interchanging the slots. The separator qualifies for neither endpoint restriction.

For each order of the special vertices, process the spine from first to last. Allocate each new vertex to one run. The state of a run records whether it is empty, occupied with an admissible last vertex, or occupied with an inadmissible last vertex. When entering an empty run, check the incoming edge from its preceding special vertex or the required first contact. On every append, replace the last-vertex state. At termination check all occupied ends, all directly adjacent special vertices separated by empty runs and any first/last restriction. Every path has one allocation. The largest machine has six permutation choices and \(3^4\) states for each. Its arithmetic cost is \(O(n)\) additions on integers of \(O(n)\) bits. Full subset path enumeration independently matches all seven readings at \(n=4,8\).

## Why the searches terminate

Write the lower chart as
\[
 K+t x+c y+qxy+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y)).
\]
For a single-run word \(y=(2^w-1)2^v\) with \(w\ge1\), the address recurrence gives \(f(y)=(2^v-1)(2^w-1)+[v>0]\). Enumerate all single-run \(x\) up to \((N-K)/t\), including zero. For fixed \(x,v\), after subtracting the constants the residual must be divisible by \((c+qx+h)2^v-h\), and the quotient must have the form \(2^w-1\). This tests every width at once.

For unrestricted \(y\) at fixed \(x\) and fixed \(j=y\bmod2\), let \(s=c+qx,\ M=N-K-tx-r[x\text{ odd}]-hf(x)-ej\). Since \(0\le f(y)\le y\),
\[
 \left\lceil M/(s+h)\right\rceil\le y\le \lfloor M/s\rfloor,\qquad
 sy\equiv M\bmod h .
\]
Testing the resulting progression in the original recurrence is exact. On a one-sided word this gives the axis search; exchanging the roles of \(x,y\) tests the other axis. For two nonempty words at least one address is at most \(\lfloor\sqrt{(N-K)/q}\rfloor\). Try that smaller address on each side, precomputing its \(f\) from its half. At \(n=36\) this bound is 29,408,100, and at \(n=40\) it is 3,267,567. It decreases at every tested larger level.

Run `python3 -B fermat_receiver.py --through-36` from this directory for the full exclusion. Omitting the flag starts the two-sided unrestricted search at 40. Planted controls exercise one-run words on the two sides, unrestricted two-sided words in either orientation and unrestricted one-sided words. [Exact machine and inverse searches](fermat_receiver.py)
