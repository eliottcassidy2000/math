# Dyadic guards in the ternary descent tree: a window and controller obstruction

Let \(h_s=7(40\cdot2^s+1)\), \(s\ge1\). The [literal seven-ray fork](../chanpak_ternary_buckets_20261122/PHYSICAL-FORKS.md) attains \(h_s\) and supplies a common center of height \(3h_s\) with reversals to \(3h_s\pm2\). Study an optimistic relaxation of **this source ray alone**: after entering one of these guards, grant every later ternary correction, whether or not physical.

## Exact inverse and missed windows

Put \(q=(H-1)/2\). The three moves \(H\mapsto3H-2,3H,3H+2\) become \(q\mapsto3q,3q+1,3q+2\). The ray inputs become \(q_s=140\cdot2^s+3\). So the relaxed set \(\mathcal R\) consists exactly of positive integers whose standard base-three expansion begins with that of some \(q_s\). If \(N\) were a first additional missing odd height, membership of \((N-1)/2\) identifies an explicit point on its inverse orbit where the ray could start a proof, but does not supply the missing physical transitions. Nonmembership rules out a proof using this inverse and only these guard stops, even after granting all suffix digits. It does not exclude another construction for \(N\).

Set \(X=3^{L-1}\), so the length-\(L\) shell \([X,3X)\) contains \(2X\) integer targets. A guard of length \(\ell_s\le L\) gives exactly the descendants in
\[
 [q_s3^{L-\ell_s},(q_s+1)3^{L-\ell_s}).
\]
It contributes nothing if \(\ell_s>L\). The covered fraction, even before removing overlaps, is
\[
 \frac32 \sum_{s:\ell_s\le L}3^{-\ell_s}
 < \frac32\sum_{s\ge1}\frac1{140\cdot2^s}=\frac3{280}.
\]
There are at most \(\log_2(3X/140)\) contributing intervals. Their complement has at most one more component than there are intervals, and hence contains a consecutive gap of length \(\Omega(X/\log X)\) in **every** sufficiently large shell. The result transfers to odd heights \(H=2q+1\). For any fixed \(0\le\varepsilon<1/2\), arbitrarily large whole windows of width \(X^{1/2+\varepsilon}\) are therefore untouched even by the optimistic relaxation.

More generally, for any anchor set, a shell of size \(\asymp X\) with \(K\) contributing anchors and uncovered proportion at least \(\eta>0\) has a contiguous uncovered gap of length \(\gg\eta X/(K+1)\). With uncovered proportion bounded below, a worst-case window bound \(O(X^{1/2+\delta})\) needs \(\Omega(X^{1/2-\delta})\) active anchors rather than \(O(\log X)\). This is a scale mismatch with the square-root balance mentioned in the RH argument: a contradiction would require a proved signal on the **actual least-hole orbit** with strictly smaller error, not simply an envelope valid for another window or on average. No zeta estimate is used here.

## Dense leading digits do not imply a hit; no finite exact phase screen

The guard cylinders are nevertheless dense in the space of infinite ternary strings with nonzero first digit. Indeed
\[
 \log_3q_s=s\log_3 2+\log_3 140+O(2^{-s}).
\]
Since \(\log_3 2\) is irrational, its irrational-rotation orbit visits every open subinterval of the circle infinitely often. The leading-digit interval for any specified valid prefix has nonempty interior. Shrink it slightly and choose arbitrarily large \(s\) with the unperturbed normalized mantissa in the smaller interval; then the vanishing perturbation preserves membership. For large enough \(s\) the guard itself extends the specified prefix. The complement of the union of guard cylinders is thus closed and nowhere dense but has measure at least \(277/280\) under uniform ternary cylinder measure on leading-nonzero strings.

**Every regular superlanguage of \(\mathcal R\) among canonical ternary strings accepts a fraction tending exponentially to one in complete length shells.** To prove this, fix a complete deterministic automaton accepting the superlanguage. For any state reached after a valid nonempty prefix \(p\), choose a guard extending \(p\). From the state reached on that guard *all* continuations accept, since they remain in \(\mathcal R\). Such a permanently accepting state remains permanently accepting under each digit, and a shortest path to the set of such states has length at most the number \(D\) of states. On independent uniform suffix digits, each block of \(D\) steps has conditional probability at least \(3^{-D}\) of reaching that set if it has not been reached already. The rejection fraction at length \(L\) is therefore at most \((1-3^{-D})^{\lfloor(L-1)/D\rfloor}\). The exact guard language has acceptance fraction below \(3/280\), so it is **not regular**. A fixed finite-state screen with no false negatives accepts asymptotically at least \(277/280\) of every shell as false positives relative to actual relaxed coverage.

This excludes an exact finite-state absolute-height stopping rule only for this one ray and inverse. It does not exclude unbounded integer registers, different guard sets, or new physical transition rules.

## The exceptional set is large even on the two first-hole candidate types

The [inherited first-hole classification](../least_hole_near_twins_20261117/README.md) restricts a first additional hole to a prime or seven times a prime. This restriction does not rescue the seven-ray inverse. In each sufficiently large height shell \([2X+1,6X+1)\) the fraction of prime heights in \(\mathcal R\), and separately the fraction of seven-times-prime heights in \(\mathcal R\), is at most \(3/280+o(1)\), where membership refers to the half-index.

Indeed the union of covered intervals has at most \(O(\log X)\) disjoint components and total real length at most \(4X(3/280)\) after mapping \(q\mapsto2q+1\). On each component \([a,b)\) the [unconditional quantitative prime number theorem](https://doi.org/10.1016/j.jmaa.2023.127460) gives
\[
 \pi(b)-\pi(a)=\int_a^b\frac{dt}{\log t}
             +O\bigl(Xe^{-c\sqrt{\log X}}\bigr)
\]
with harmless endpoint adjustments. The main term is \((b-a)(1+O(1/\log X))/\log X\), uniformly in this shell, while the sum of the error terms is \(o(X/\log X)\). Dividing the endpoints by seven gives exactly the same fraction bound for seven-times-prime heights. Moreover the single largest uncovered interval already contains \(\Omega(X/\log^2 X)\) heights of each type, by applying the same estimate there and on its dilation by \(1/7\). These are missed by the relaxed inverse, not proved missing from the tournament spectrum.

There is an abstract first-counterexample model, not just a counting objection. Fix any sufficiently large missed prime or missed \(7p\) from these intervals, and call it \(N\). Let \(S_N\) be the multiplicative closure of (i) all allowed odd values below \(N\) and (ii) all heights of the relaxed guard language. It obeys ordinal-product closure, contains every certified small height, and contains the ray with **all** its arithmetic descendants. Yet \(N\notin S_N\): it is not itself a generator, and has no factorization into nontrivial available factors (a prime cannot factor, and \(7p\) would require 7). Thus \(N\) is exactly the first additional hole in this abstract set. This is not a model of all tournament operations, but proves that the least-counterexample principle, small coverage, ordinal multiplication and arbitrary suffix moves from this one ray together do not force a contradiction.

The regular-superlanguage obstruction also persists on these candidate types. Its exponentially small rejection fraction on shells gives \(O(X^{1-\gamma})=o(X/\log X)\) rejected integers for some \(\gamma>0\) depending on the fixed automaton. Therefore a finite regular overapproximation misses a vanishing fraction of either prime candidate type, while admitting at least \(277/280-o(1)\) of each as false positives relative to the relaxed guard language.
