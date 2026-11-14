# Approximate readiness forces low outdegree

The [previous ready-root theorem](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/READY_CLASSIFICATION.md)
showed that exact \(R_+=H\) forces root outdegree at most one.
The same split picture yields a quantitative bound for any
rooted tournament, useful when a construction permits a
nonzero but tightly controlled defect.

Let \(d=\deg^+(v)\ge1\), and let \(A\) be all strong
components strictly before the component of \(v\).
Then
\[
             R_+(T,v)-H(T)\ \ge\ H(A)(2^{d-1}-1).
\]
In particular, with defect \(\Delta\), any candidate root has
\(d\le1+\log_2(\Delta/H(A)+1)\).
The exponential base 2 cannot be improved.

## Proof by a gateway that leaves the other choices free

As before, \(R_+-H\) counts exactly the ordered complementary
path pairs for which the first path ends at an outneighbor
of \(v\) and the second starts at \(v\), so they cannot join.

Let \(W\) be the reachable suffix of \(v\), let \(O\) be its
\(d\) outneighbors, and let \(I=W\setminus(O\cup\{v\})\).
If \(I\ne\varnothing\), some \(x\in O\) has an arc into the
initial SCC of \(T[I]\): all of \(I\) is reachable from \(v\)
in \(W\), and a first entry into that component cannot come
from \(v\) or a later component of \(I\).
If \(I=\varnothing\), fix any \(x\in O\).
For each nonempty \(Y\subseteq O\setminus\{x\}\), keep \(x\)
and delete \(Y\). All remaining outneighbors are directly
reachable from \(v\), and \(x\) reaches all of \(I\) through
its initial SCC if there is one. Thus \(W-Y\) has a spanning
path starting at \(v\) (rotate a Hamiltonian cycle in its
first strong component, then join later components).
In the complementary path use any Hamiltonian path in \(A\)
followed by one in \(Y\). Its last vertex lies in \(Y\)
and is beaten by \(v\). These are distinct nonjoinable
splits for distinct \(Y\) or distinct choices of the path
in \(A\). Counting them proves the bound.

For the base-sharp examples take a directed triangle of
modules \(v\to TT_d\to z\to v\), with \(v,z\) singleton.
A complement path containing \(z\) ends there (it is a sink
within the complement without \(v\)), so is not bad.
Every nonempty proper subset \(Y\subset TT_d\) gives one
bad split: its complementary root-start path needs at
least one remaining transitive vertex to reach \(z\).
There are exactly \(2^d-2\) bad splits.
Thus the universal upper and lower examples differ only
by a constant factor in scale. This is not an assertion
of the exact minimum defect for every \(d\).

The gateway lower bound was found on the arbitrary-ready route;
the explicit upper family pins down its exponential scale.
