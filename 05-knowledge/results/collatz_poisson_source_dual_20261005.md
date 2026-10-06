# A source-sensitive Poisson dual and its certificate boundary

**PROVED:** the bounded adjoint-dual criterion, its residual error bill, the finite-support path equivalence, the fixed-price deadline, and the typed lift of the signed atom measurement. **FINITE-EXACT:** the explicit actual and synthetic controls. **OPEN:** independent positive global dual witnesses for every supplied source. A unique infinite solution need not have positive mass at every coordinate.

[Program](../../04-computation/experiments/collatz_poisson_source_dual_20261005.py) and [output](collatz_poisson_source_dual_20261005.out).

## 1. Inheritance and the new obligation

The closest mechanism is the contractive operator in [Green weights, sections 2–3](collatz_green_weight_20261005.md), including the killed ROOT row and its exact column sums. That note and incoming C8 of [critical flow](collatz_three_bits_critical_flow_20261005.md) already prove existence, uniqueness in the specified normed space, computability, and support exactly equal to the ROOT basin. Those facts are inherited here.

The [source refinement note](collatz_source_refinement_floor_20261005.md), SF3–SF4, adds signed-bill cancellation and the finite energy of one fixed evaluation direction. The [refinement/shadow note](collatz_refinement_floor_shadows_20261005.md), T2, identifies the stopping-time content of a positive floor. The sharper [floor deadline](collatz_floor_transport_deadlines_20261005.md) gives a total quantitative extractor once a truthful mixture-weight floor is available. No other broad assertion in the shadow note is needed.

The new object is a **source-sensitive adjoint subsolution**. It is a precise obligation to which one can work backwards from a desired floor. Its inequality direction, norm domain, root boundary, and global residual all matter. The canonical hostile is a contracted component disconnected from ROOT; the corrected near miss is “the fixed point is unique, hence all its coordinates are positive.” The least-used sidecar is the boundary term left when a dual function is unbounded.

The board is **selected source / root forcing / bounded dual / global residual / finite support / independent positivity**. The source of the transfer is the exact labelled base kernel, the target is a floor and literal receipt, the map is adjoint pairing, the preserved predicate is the actual source's ROOT membership, and the data that cannot be discarded are the global inequality and norm bound.

## 2. Explicit operator and adjoint

Let \(U(n)=\operatorname{oddpart}(3n+1)\), \(S(n)=4n+1\), and

\[
 {\cal B}=\{b>0:b\text{ odd},\ v_2(3b+1)\in\{1,2\}\}.
\]

Every positive odd \(n\) has a unique decomposition \(n=S^j(b)\). For a nonroot base \(b\), write \(U(b)=S^{k_b}G(b)\). For \(0<r<1\), define

\[
 (A_r f)(1)=0,\qquad
 (A_r f)(b)=(1-r)r^{k_b}f(G(b)).
 \tag{1}
\]

The root row is deleted. It is not a discounted physical ROOT self-loop.

The possible inverse children of \(c\) are the unique bases \(b_0(S^k c)\), where

\[
 b_0(y)=
 \begin{cases}
  (2y-1)/3,&y\equiv5\pmod6,\\
  (4y-1)/3,&y\equiv1\pmod6.
 \end{cases}
\]

Targets divisible by three and the child 1 are omitted. Since \(S^k c\equiv c+k\pmod3\), one depth class modulo three is forbidden. The inherited exact norm is

\[
 \|A_r\|_{\ell^1\to\ell^1}
 =\|A_r^*\|_{\ell^\infty\to\ell^\infty}
 =\kappa_r=\frac{1+r}{1+r+r^2}<1.
 \tag{2}
\]

The adjoint is therefore the absolutely convergent sum

\[
 (A_r^*\phi)(c)=
 \sum_{\substack{b\ne1\\G(b)=c}}(1-r)r^{k_b}\phi(b).
 \tag{3}
\]

The inherited Green vector is

\[
 g_r=(I-A_r)^{-1}\delta_1
     =\sum_{\ell\ge0}A_r^\ell\delta_1\in\ell^1.
 \tag{4}
\]

It is nonnegative, has root value one, and is zero exactly at bases outside the ROOT component.

Most computations below fix \(r=1/2\), so the edge factor is \(2^{-k_b-1}\), the exact norm is \(6/7\), and the killed root column is \(5/14\). The elementary first-column bound gives \(\|g_{1/2}\|_1\le7/2\). This is a convenient valid bound, not the best current bound: the two-generation estimate in [adaptive mixture flow, P1](collatz_adaptive_mixture_flow_20261005.md) improves it to \(11/4\).

## 3. The source dual and its error bill

**P1 — PROVED.** Let \(q,\phi\in\ell^\infty({\cal B})\). If the coordinatewise global inequality

\[
 (I-A_r^*)\phi\le q
 \tag{5}
\]

holds, then

\[
 \boxed{\ \phi(1)\le\langle q,g_r\rangle.\ }
 \tag{6}
\]

Indeed multiply (5) by the nonnegative \(g_r\), sum, and use

\[
 \langle (I-A_r^*)\phi,g_r\rangle
 =\langle\phi,(I-A_r)g_r\rangle=\phi(1).
\]

Every pairing is justified by \(g_r\in\ell^1\), boundedness of \(\phi,q\), and (2).

In particular take \(q=\delta_b\). A proved rational \(\delta>0\) with \(\phi(1)\ge\delta\) gives the exact selected-source floor \(g_r(b)\ge\delta\). The witness can be analytic or represented by a finite program; finite support is not required. What is required is a proof of (5) on every base and a valid norm domain. A grid of successful inequalities does not supply that proof.

There is a unique bounded exact Poisson solution

\[
 \phi_q=(I-A_r^*)^{-1}q,\qquad
 \|\phi_q\|_\infty\le\frac{\|q\|_\infty}{1-\kappa_r},
 \tag{7}
\]

and \(\phi_q(1)=\langle q,g_r\rangle\). Thus the supremum of the valid lower readouts in (6) is exact. Uniqueness of this bounded solution does not prove that its root coordinate is positive.

At \(r=1/2\), a globally certified approximate inequality

\[
 (I-A^*)\phi\le q+\eta{\bf1},\qquad \eta\ge0
\]

still gives

\[
 \langle q,g\rangle\ge\phi(1)-\frac72\eta.
 \tag{8}
\]

One may substitute any better proved total-mass bound, including \(11/4\). This is the concrete residual bill for a compressed ansatz: prove a uniform residual bound, retain its sign, and make the right side positive. Approximation error is not an independent source floor unless it passes this test.

## 4. Finite support records a finite path

**P2 — PROVED.** Fix \(r=1/2\), a base \(b\), and a finite set \(F\) of bases. There is a function supported on \(F\) satisfying

\[
 (I-A^*)\phi\le\delta_b,\qquad \phi(1)>0
 \tag{9}
\]

if and only if the actual inverse graph restricted to \(F\) contains a finite path from ROOT to \(b\).

For necessity, start at \(c=1\). If \(c\ne b\) and \(\phi(c)>0\), (9) gives

\[
 \phi(c)\le\sum_{G(v)=c}d_v\phi(v).
\]

The total column mass is at most \(6/7\). Hence some child satisfies

\[
 \phi(v)\ge\frac76\phi(c)>0.
 \tag{10}
\]

For finite support choose a child maximizing \(\phi(v)\). Values strictly increase, so no vertex repeats; the process cannot stop except at \(b\). Every visited vertex lies in the positive part of \(F\). Signed entries cannot evade this argument.

For sufficiency, take an actual inverse path \(1=c_0,c_1,\ldots,c_\ell=b\). Put \(\phi(b)=1\), and recursively \(\phi(c_{i-1})=d_{c_i}\phi(c_i)\); put zero elsewhere. Its defect is exactly \(\delta_b\), and its root value is the positive product of path weights. Its unique forward parent path gives the actual receipt.

The global check for a finite packet is itself finite: besides its support and the target, only parents of nonzero support entries can have a nonzero adjoint term. The program checks all those constraints, rather than omitting boundary equations. A finite LP over such a declared set therefore detects precisely its observed rooted paths. It can discover a previously unnoticed path in a union of observations; it cannot claim support outside that graph merely through optimization.

There is also a bounded, possibly infinite version. If \(\|\phi\|_\infty\le M\) and \(\phi(1)\ge\delta>0\), the same argument forces a finite inverse path of length less than any integer \(T\) with

\[
 (6/7)^T M<\delta.
 \tag{11}
\]

For a computable infinite witness, exact maximum comparisons need not be decidable. Choose any rational \(1<a<7/6\), such as \(13/12\). Enumerate children and increasingly precise value intervals until finding \(\phi(v)>a\phi(c)\). Such a child exists by (10), and the search terminates. The resulting path has length less than any \(T\) with \(a^T\delta>M\). No uniform bound on the required sibling depths or numerical precision is asserted.

The bounded exact solution (7) is unique. Arbitrary subsolutions and their finite proof descriptions are not unique. A positive finite Collatz source has one strict actual forward ROOT word; different dual proofs can certify that same word.

## 5. A fixed-price floor gives a logarithmic counter deadline

For a rooted nonroot integer with actual strict ROOT valuations \(a_1,\ldots,a_\tau\), the inherited counters are

\[
 L=\tau-1,\qquad K=\sum_i\lfloor(a_i-1)/2\rfloor,\qquad N=L+K.
\]

The final nonroot valuation is even and at least four, so \(K\ge1\) and \(\tau\le N\). At \(r=1/2\),

\[
 f_{1/2}(n)=2^{-N}.
\]

For \(n=S^j(b)\), a base dual floor \(g_{1/2}(b)\ge\delta_b\) gives \(f_{1/2}(n)\ge\delta=2^{-j}\delta_b\). Thus

\[
 N\le C=\left\lfloor\log_2(1/\delta)\right\rfloor.
 \tag{12}
\]

For exact \(\delta=p/q\), the implementation computes \(C\) as the bit length of \(\lfloor q/p\rfloor\), minus one. ROOT is handled separately. A true floor therefore pays for at most \(C\) actual odd \(U\)-steps. The function receipt_from_green_floor carries out this bounded replay and rejects a contradicted claimed floor. A returned receipt is always literally checked; the routine does not certify an unsupported analytic premise.

The same floor yields a mixture-weight bound. For \(C\ge1\),

\[
 \boxed{
 W(n)\ge
 \frac{2}{(C+2)\binom{C+1}{\lfloor(C+1)/2\rfloor}}.
 \ }
 \tag{13}
\]

For a fixed \(N\), the formal minimum of
\(w(L,K)=2/[(N+2)\binom{N+1}{K}]\), \(1\le K\le N\), occurs at a central binomial coefficient. The minima decrease with \(N\), proving (13). It is a safe minimum over formal counters, without asserting that every minimizing pair is realized by an integer. At \(C=0\), a truthful floor only permits ROOT and its weight is one.

The logarithmic bound is for the fixed-price floor \(\delta\), not the mixture floor \(\epsilon\) of the previous package, whose uniform step bill is \(O(\epsilon^{-1/2})\). Their measurement inputs have different types and may differ greatly in magnitude.

## 6. Exact connection to the signed atom measurement

The preceding direct dual uses one fixed-price base coordinate. It is not already the localized mixture measurement \(H_{m,d}\). The change of carrier is explicit.

Let \(q_m(j)\) be any bounded signed selector on leaf indices, with \(q_m(m)=1\) and \(q_m(j)\le0\) for \(j\ne m\). This includes
\(q_m(j)=(9h_m(j)-8)h_m(j)^d\) from [localized resolvent floors](collatz_localized_resolvent_floor_20261005.md). Extend it to odd sources by

\[
 \sigma(n)=
 \begin{cases}
 q_m((n-3)/6),&3\mid n,\\
 0,&3\nmid n.
 \end{cases}
\]

For each price \(r\), set

\[
 \gamma_r(b)=\sum_{j\ge0}r^j\sigma(S^j b).
 \tag{14}
\]

This observable uses only explicit source arithmetic, not ROOT membership. If \(|\sigma|\le M_\sigma\), the tail after depth \(J\) is at most \(M_\sigma r^{J+1}/(1-r)\), so it is computable at each fixed rational \(r\) provided the selector is computable and a valid bound is known. The concrete localized selector has both properties; boundedness alone does not imply computability.

The unique sibling decomposition and inherited beta mixture give

\[
 \begin{split}
 A_m
 &=\sum_{\ell\ge0}q_m(\ell)\lambda(6\ell+3)\\
 &=\int_0^1 2(1-r)\,\langle\gamma_r,g_r\rangle\,dr.
 \end{split}
 \tag{15}
\]

Absolute convergence follows from
\(\int2(1-r)\sum_b|\gamma_r(b)|g_r(b)\,dr
 \le M_\sigma\sum_{\text{leaves}}\lambda=M_\sigma\).
This uses the inherited probability normalization of the injection measure.

Consequently a measurable family of bounded global subsolutions

\[
 (I-A_r^*)\phi_r\le\gamma_r
 \tag{16}
\]

with root readouts integrable against \(2(1-r)\,dr\) supplies

\[
 A_m\ge\int_0^1 2(1-r)\phi_r(1)\,dr.
 \tag{17}
\]

A proved positive rational lower bound on the last integral is an independent positive measurement bound, hence a floor at the exact leaf \(6m+3\). This is the proposed constructive research obligation. No family meeting it at every unresolved source is produced here.

The sign structure survives the lift: \(\gamma_r\) can be positive only at the single primitive base of \(6m+3\), since every other sibling ray sees nonpositive selector values. Therefore a positive root value in (16) still forces a finite path to that base by the maximum-principle argument. The parameter integration does not erase source identity. No selection of a rational parameter from a merely measurable positive set is assumed; a quantitative bound in (17) can instead be passed directly to the mixture-floor receipt compiler.

## 7. Three boundaries that uniqueness cannot remove

**Disconnected contracted state.** On two labelled states, ROOT and \(z\), take \(A=\operatorname{diag}(0,1/2)\). The unique summable and bounded solution of \(g=\delta_{\rm ROOT}+Ag\) is \((1,0)\). This is a synthetic model, not a proposed positive Collatz cycle. It proves that contraction, a unique solution, and positive root boundary do not alone imply full support.

**Positive regularization with zero limit.** Change the forcing to \(\delta_{\rm ROOT}+\epsilon\delta_z\). The unique solution becomes \((1,2\epsilon)\), positive at both states for every \(\epsilon>0\), but the new coordinate tends to zero. In general,
\(\|g_\epsilon-g\|_1\le
 \epsilon\|h\|_1/(1-\kappa)\)
for forcing \(\delta_{\rm ROOT}+\epsilon h\).
A positive lower bound uniform in \(\epsilon\), after paying this error, would be useful; positivity at each separate regularization is not that bound. The forcing has changed and must be retained.

**Unbounded harmonic boundary.** On the rooted ray \(0,1,2,\ldots\), let
\((Af)(0)=0\), \((Af)(j)=f(j-1)/2\) for \(j\ge1\).
Its norm on \(\ell^1\) is \(1/2\), its Green vector is \(g(j)=2^{-j}\), and the unbounded function \(\phi(j)=2^j\) satisfies \((I-A^*)\phi=0\) while \(\phi(0)=1\). Incorrectly applying (6) would give \(1\le0\). The missing tail is exactly
\(\langle\phi,A^N\delta_0\rangle=1\) for every \(N\).
Boundedness, or an independently proved vanishing boundary term replacing it, is essential. A finite cutoff that drops its last constraint makes the same mistake.

The killed ROOT row is also part of the operator specification. Retaining a ROOT self-loop of weight \(1/2\) changes the root Green value from one to two; it is a different boundary problem.

## 8. Computation and remaining independent step

Run:

~~~text
python 04-computation/experiments/collatz_poisson_source_dual_20261005.py
python -O 04-computation/experiments/collatz_poisson_source_dual_20261005.py
~~~

The verifier checks all 128 positive odd sources below 256 using independently constructed finite control paths, with declared base-depth cap 256. These already completed sources test the dual, exact source retention, inverse-to-ordinary-word decoding, and fixed-price floor compiler; they are not new coverage. There are 124 deleted-target hostiles and 48 exact inverse-column controls.

Independently, all 700 functional parent graphs on two through five labelled vertices are tested with deleted root row and a strict rational contraction. Exact matrix inversion gives 3,412 source-specific Poisson duals; their positive root readout agrees with literal component reachability. These graphs are abstract controls, not invented Collatz transitions.

Further controls include nine regularizations, sixteen harmonic-ray boundaries, 27 leaf/base regroupings, uniform global-residual error, and eleven malformed or false-premise inputs. Every check uses integers or Fraction and survives optimized Python. Normal and optimized output agree, with **26,112 exact checks**.

The next useful input is a bounded analytic or compressed \(\phi\) with a globally proved residual small enough to make (8) positive, or a price-indexed family satisfying (16)-(17). A finite exact dual can organize and verify known edges; its positivity supplies a real path. An infinite uniquely defined dual is also a legitimate object, but neither uniqueness nor convergence of its approximations supplies the remaining strictly positive root readout.
