# Positive square-root reduction in the relaxed seven-axis source problem

In the [terminal-seven decoder](AXIS_FRONTIER.md), fix the right
width power four, so the ordinary total is seven times the
readout of the source by the row \(w=(47,56)\).

The cone and parity alone leave no eventual arithmetic gap for
the smallest right multiplier \(w=(47,56)\): write \(a=1+2i\),
\(b=2a+2j\), with \(i,j\ge0\). Then
\[
                       p=159+2(159i+56j).
\]
Since 159 and 56 are coprime, the exact Frobenius cutoff
shows the largest missed odd \(p\) in this **relaxed cone**
is \(159+2(159\cdot56-159-56)=17537\). This calculation
knows nothing about how to build \((a,b)\) as a prefix.

In fact one can enforce arbitrarily many fixed physical inverse
steps before reaching the remaining relaxed cone. Put
\(M=D_2C_1=(1,2;4,4)\), corresponding to one zero gap and the
next width one, and write \((u_k,v_k)=w^t M^k\) (as a row),
\(A_k=u_k+2v_k\). All these rows remain primitive odd/even:
an odd common divisor of the next row would divide both
entries of the preceding row after multiplying by the integral
adjugate, because \(\det M=-4\). The cone is preserved by \(M\).
Thus for each fixed \(k\), every odd \(p\) greater than
\[
                          2A_kv_k-A_k-2v_k
\]
has an exact representation
\[
                     p=w^t M^k (a,b)^t
                       =A_k+2(A_ki+v_kj)
\]
with \((a,b)\) still in the cone. The cutoff is again the exact
two-generator Frobenius formula since \(\gcd(A_k,v_k)=1\).
It is below \(A_k^2\). Also \(M(1,2)^t\le6(1,2)^t\),
so \(A_k\le159\cdot6^k\).

More quantitatively, for \(p\ge159^2\), take the largest
\(k\) with \(A_k^2\le p\). Then \(\sqrt p/6<A_k\le\sqrt p\),
and \(v_k\ge A_k/3\), because \(v_k\ge u_k\).
Represent \(p\) and take the least positive odd \(a\) in its
required class modulo \(2v_k\); this merely lowers \(a\) and
raises \(b\), preserving the cone. Consequently \(a\le\sqrt p\)
and \(b<18\sqrt p\). So every such target has a target-directed
exact chain of \(k=\Theta(\log p)\) legal inverse transfers
ending at a two-coordinate cone state of size \(O(\sqrt p)\).
To turn it into a physical construction, the unresolved final
condition is that the resulting cone state is the state of an
actual prefix. One may not infer this from the many already
valid inverse stages. The sparse terrace bound shows that this
last condition cannot always hold inside the same family.
