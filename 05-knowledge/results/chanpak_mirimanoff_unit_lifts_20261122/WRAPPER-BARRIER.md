# A speed limit for physical perfect-number lifts

Put \(E_p=2^{p-1}(2^p-1)\). The results concern literal growth of a marked tournament pair at \(E_p-1,E_p+1\), not an arbitrary new pair at the next exponent. The first theorem rules out the common idea of hiding the whole previous tournament in one growing homogeneous module. Its hypotheses leave mixed-contact extensions and re-encoding the old source open.

## 1. Homogeneous substitution cannot give a fixed-step perfect lift

Let \(T\) be an \(n\)-vertex tournament. Substitute it for a marked vertex \(z\) in any quotient tournament \(D\), with every exterior vertex having the same direction to all vertices of \(T\). Let \(Z_T(j)\) count ordered \(j\)-tuples of disjoint directed paths spanning \(T\), with empty slots allowed. Define \(B_r(T)\) in the same way but require all \(r\) slots nonempty.

Each Hamilton path upstairs collapses each maximal run within \(T\) to \(z\). For each compressed path with \(r\) occurrences of \(z\), there are \(B_r(T)\) ways to fill it, since all exterior incidences are uniform. If \(c_r\) is the number of such skeletons, this yields
\[
 H(D[z\mapsto T])=\sum_{r\ge1}c_rB_r(T),\qquad
 c_1=H(D),\quad B_1(T)=H(T),\quad B_2(T)=Z_T(2)-2H(T).
\]

If \(z\) is its own strong component, every path visits it in only one run because a directed path cannot leave and reenter an SCC. Thus the output is exactly \(H(D)H(T)\). If \(z\) lies in a nontrivial SCC, then \(c_1\ge3\) and \(c_2\ge1\). Indeed, a spanning cycle of the SCC gives at least three ordinary paths by cutting the cycle, and following the full cycle from \(z\) back to \(z\) gives a valid skeleton with two marked runs. Glue Hamilton paths in adjacent SCCs to get the full skeleton. The bound is therefore
\[
 H(D[z\mapsto T])\ge H(T)+Z_T(2)\ge(n+2)H(T).
\]
For the second inequality, cut a Hamilton path of \(T\) at any of \(n+1\) boundaries to produce an ordered two-slot cover. All cut results are distinct, since concatenation recovers the path and the first slot recovers the cut. More generally weak cuts to \(j\) slots give \(Z_T(j)\ge\binom{n+j-1}{j-1}H(T)\). The first inequality is attained by the directed-triangle quotient with its other two vertices singleton; its run signature is \((c_1,c_2)=(3,1)\). Hence \(H+Z_T(2)\) is the exact minimum over quotients with \(z\) in a nontrivial SCC.

Apply the two cases to \(T_\pm\) differing on one literal arc with \(H(T_\pm)=E_p\pm1\). If the same quotient is used, a singleton-SCC substitution can preserve the unit odd gap only with \(H(D)=1\), leaving the height unchanged. The other way is too fast: for \(p\ge3\) and \(H(T)\in\{E_p-1,E_p+1\}\), the next endpoint is strictly below \(5H(T)\), and \(n\ge4\) already gives factor at least 6. Indeed with \(t=2^{p-1}\ge4\), \(E_{p+1}=8t^2-2t\) while \(5(E_p-1)-(E_{p+1}+1)=2t^2-3t-6>0\). This excludes even exponent-dependent quotients and does not assume a bounded number of their vertices.

The same gives an asymptotic no-go for every fixed increment \(d>0\), because \(H(T)\le n!\) forces \(n\to\infty\) along the perfect endpoints, whereas the target ratio tends to \(4^d\). More quantitatively the nontrivial wrapper requires \(4^d(1+O(2^{-p}))\ge n+2\), with the error uniform for \(d\ge1\).

## 2. A cyclic wrapper has a separate mod-three obstruction to preserving the gap twice

Consider the particularly natural outer wrapper \(W_k(T)=C_3[T\oplus I_k,1,1]\), where \(I_k\) is a transitive ordinal tail, \(k\ge0\). For any nonempty \(T\), \(j\mid Z_T(j)\): cyclic rotation of the \(j\) slots acts freely, because a repeated nonempty slot would use its vertices twice. This applies equally to a difference of two profiles.

If \(f=Z_T\), the full ordered-slot transfer is
\[
 Z_{W_k(T)}(j)=
 j\big[(j+1)^k f(j+1)+j^{k+1}f(j)-(j-1)^k f(j-1)\big].
\]
For a proof, put \(U=T\oplus I_k\) and orient the two singletons as \(U\to a\to b\to U\). If they occupy one slot in order \(a,b\), then the remaining \(U\)-pieces use \(j+1\) ordered subslots; if they occupy one slot in order \(b,a\), the intervening \(U\)-subslot must be nonempty, giving \(Z_U(j)-Z_U(j-1)\). There are \(j\) choices of that slot. If they occupy different slots, there are \(j(j-1)\) placements and \(j\) independent ordered \(U\)-subslots. Finally ordinal tails multiply \(Z_T(j)\) by \(j^k\). Set \(f(0)=0\) since \(T\) is nonempty.

In particular \(H(W_k(T))=H(T)+2^kZ_T(2)\). Suppose an initial pair has \(\Delta H=2\) and retains this same signed gap after \(W_k\) and after another \(W_\ell\), with arbitrary nonnegative \(k,\ell\). The first equality forces \(\Delta f(2)=0\), and the second forces \(\Delta Z_{W_k}(2)=0\). The displayed transfer then reduces to \(3^k\Delta f(3)=2\), impossible modulo 3, since \(\Delta f(3)\in3\mathbb Z\). This does not use the perfect target and holds at arbitrary heights. In particular unit seams are not invariant under recurrent positive triangle bordering, even where an individual border can preserve one.

**Even allowing the two signs to exchange, four consecutive wraps are impossible.** Let $d_j$ denote the difference of the $j$-slot counts at a given stage, and suppose $d_1=2s$ for $s\in\{1,-1\}$. For any pairs of equal order one has $4\mid d_2$: terms with one empty slot contribute $2d_1$, and each unordered nonempty bipartition contributes twice the difference of two odd products of induced Hamilton counts. Thus a sign flip under $W_k$ requires $k=0$ and $d_2=-4s$; a sign-preserving wrap requires $d_2=0$. In either case the new $d_2'=2[3^k d_3+2^{k+1}d_2-2s]$ is nonzero modulo 3. The next transition, if a unit gap is retained, must therefore flip signs and have zero padding. After the first transition all transitions must flip.

Normalize the initial difference to $2$. If the first transition also flips, its padding is zero. If it keeps the sign, initial $d_2=0$, and the second transition flips, so first-stage $d_2'=-4$ and initial $d_3=0$. If the third transition retains unit magnitude, it must flip again. The two successive flips starting from that first-stage profile require $d_3'=12$ by the two-slot transfer. That profile has $d_3'=3\cdot4^k d_4$, so $4^k d_4=4$. Since $4\mid d_4$, the initial padding is zero in this case as well.

Now all transfers are unpadded. The first five read rows of the transfer, from time zero through time four, are respectively
\[
 (1),\quad (1,1),\quad(-1,5,2),\quad(-11,13,28,6),\quad
 (-37,-43,254,180,24).
\]
They follow by multiplying a row by the operator $f(j)\mapsto j[f(j+1)+jf(j)-f(j-1)]$. If the five heights are $2,-2,2,-2,2$, triangular inversion gives $d_5=199$; if they are $2,2,-2,2,-2$, it gives $d_5=-27$. Both contradict $5\mid d_5$. The conclusion applies to equal-order pairs at arbitrary initial heights.

## 3. Small exact mixed-contact screen

See [the directed-cut extension law and cycle gauge correction](CUT-GAUGE.md) for the coupling that nonhomogeneous attachments must satisfy.

### Small exact extension screen

The companion subset-endpoint DP extends both members by the same labelled exterior contacts and checks both exact counts, keeping the marked old internal edge. From all two lower birth classes at $p=3$ (six eligible arcs as labelled in those representatives), and from every eligible marked reversal on the complete lower birth deck at $p=5$ (eighteen after choosing one representative per class), adding up to two mixed-contact vertices yields no exact $p\mapsto p+1$ seam. From each of the two $p=3$ minimal lower isomorphism classes, three added vertices also give no exact lift; each class's three eligible arcs are a single automorphism orbit, so it is enough to test one each here. For $p=7$, the one named ten-vertex seed in the earlier birth-seam packet has no one-vertex exact lift. The $p=3$ scan also finds no lift from the complete order-seven deck with two new vertices (14 marked arcs in the chosen representatives) or from the order-eight deck with one new vertex (22 marked arcs). Thus there is no such induced lift from any $p=3$ seed to total order at most 9: smaller outputs may be padded by common ordinal sinks, preserving the two Hamilton counts. These are finite architecture tests, not exclusions of independent pairs at the next exponent. [The C++ replay](check_induced_lift.cpp) accepts arguments (number of added vertices, exponent, seed order, optional deck directory); [the Python replay](check_wrappers.py) independently checks the transfer and all tag cuts on three small hosts.
