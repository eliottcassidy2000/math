# An effective Newton-circuit law for separated narrow root blocks

Status: **PROVED + FINITE-EXACT; [independent audit accepted](continuing8_20260906_newton_ballot_audit.md)**. The all-profile assertion is proved below by an explicit finite threshold construction. The 703-gate exact bank is corroboration and a reproduction of the symbolic identities, not an extrapolation from finitely many multiplicities.

## 1. Inheritance, object, and new scope

The closest proved mechanism is the boundary spike calculation in [THM-3004, circuit sign-change cluster law and classifier refutation](../../01-canon/theorems/THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation.md), §3b: equal-multiplicity separated clusters force at least (2K-3) changes, and doubled roots attain the maximum. Its general-profile exact count and its arbitrary-separation upper bound have finite evidence rather than a universal proof. The canonical hostile is the degree-five positive-root-parameter row ((1,1,3,3,8)), with circuit signs (-,+,-), which defeats a two-end curvature classifier.

The corrected near miss is the finite-edge interpretation of [THM-3000, fixed-edge cumulant-curvature universality](../../01-canon/theorems/THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer.md): its limit is degree tending to infinity with bounded jet invoices; it is not a global circuit classification at fixed degree. The least-used sidecars here are exact reversal in [THM-3001](../../01-canon/theorems/THM-3001-newton-circuit-reversal-involution-and-two-end-curvature-law.md) and factorial PF rows in [THM-3079, Laguerre PF row transform](../../01-canon/theorems/THM-3079-laguerre-pf-row-transform-and-strict-integer-mesh-terminal-minus-one.md).

The live concept board is: boundary spikes; interior factorial logarithmic curvature; reversal; root spread and coefficient-error budgets; factorial PF factors; global no-return. The new step is to join spikes to strictly ordered interior circuits. This both controls zeros and removes the equal-multiplicity restriction, while permitting a quantified nonzero width in every block.

For (N(n)=\prod_{i=1}^d(n+r_i)), (r_i>0), let (e_k) be the elementary symmetric root parameters, (h_k=e_k/\binom dk), and

\[
R_k=\frac{h_k^2}{h_{k-1}h_{k+1}},\quad 1\le k\le d-1,
\qquad C_k=R_k/R_{k-1},\quad 2\le k\le d-1.
\]

Circuit signs mean (operatorname{sgn}(C_k-1)=operatorname{sgn}(\log C_k)). Sign changes are counted after deleting zeros.

**Theorem.** Fix (K\ge2) and integers (m_1,\ldots,m_K\ge2), with (d=\sum m_j). The construction in §3 supplies explicit rational (\varepsilon_0>0) and integer (T_0>1), depending only on this profile, such that every multiset

\[
r_{j,i}=\rho_j w_{j,i},\quad 1\le w_{j,i}\le1+\varepsilon_0,
\quad \rho_j/\rho_{j+1}\ge T_0
\tag{1}
\]

has exactly (2K-3) circuit sign changes. Write (b_j=\sum_{i\le j}m_i). At every boundary (b_j), (C_{b_j}>1) and (C_{b_j+1}<1). The whole first circuit block is positive and the whole last block is negative. Each intervening block has exactly one negative-to-positive change, allowing a zero at the change.

In particular, for (m_j=2) the number (2K-3=d-3) is maximal. This holds on open sets of distinct root parameters inside (1), as well as on coalesced roots. For any fixed profile, coalesced or sufficiently narrow blocks whose adjacent separations tend to infinity eventually satisfy the conclusion.

This theorem does not assert the arbitrary-separation upper bound for a prescribed number of distinct roots. It does not classify the untouched wall-stripped first-gap cores. Since it has both circuit signs, it obstructs both possible single-sign versions of no-return on the constructed class; it does not infer no-return properties from only a resultant norm or a finite edge jet.

## 2. Interior factorial curvature

First put all weights (w=1). If (k=A+l) fills all (A) roots before a block of size (m) and then (l) roots of this block, the dominant coefficient is

\[
t_k D_k,\qquad t_k=\binom ml,\qquad
D_k=\Bigl(\prod_{j\text{ before}}\rho_j^{m_j}\Bigr)\rho_{\text{current}}^l.
\]

At a boundary either description gives (t_k=1) and the same (D_k). Set (t_0=t_d=1), and define finite rational constants

\[
\kappa_k=\frac{t_k^2}{t_{k-1}t_{k+1}}
\frac{\binom d{k-1}\binom d{k+1}}{\binom dk^2},
\qquad q_k=\kappa_k/\kappa_{k-1}.
\tag{2}
\]

The dominant (R_k) equals (kappa_k) except at (k=b_j), where it equals (kappa_{b_j}\rho_j/\rho_{j+1}). No two boundaries are adjacent because (m_j\ge2).

For (1\le l\le m-1), put (B=d-A-m). Direct factorial cancellation gives

\[
\kappa_{A+l}=f_A(l)f_B(m-l),\qquad
f_A(l)=\frac{(l+1)(A+l)}{l(A+l+1)}=1+\frac A{l(A+l+1)}.
\tag{3}
\]

For (A>0), writing (phi_A=\log f_A),

\[
\phi_A'(l)=-\frac1{l(l+1)}+\frac1{(l+A)(l+A+1)}<0,
\]

\[
\phi_A''(l)=\frac1{l^2}-\frac1{(l+1)^2}
-\frac1{(l+A)^2}+\frac1{(l+A+1)^2}>0.
\tag{4}
\]

The latter is the difference of a strictly decreasing positive function (u\mapsto u^{-2}-(u+1)^{-2}) at (l) and (l+A). For (A=0), (f_0=1). Consequently (log\kappa_{A+l}) is strictly convex for every block, since (K\ge2) implies (A+B>0). It increases through the first block and decreases through the final block. Thus all bounded interior (q_k) in the first block exceed one, those in the final block are below one, and the bounded interior (q_k)'s in every middle block are strictly increasing. These are finite strict rational inequalities, even when the dimensions have no uniform bound.

## 3. Explicit thresholds and uniform errors

Form the finite set of strict margins

\[
\mathcal G=\{q_k:2\le k<b_1\}
\cup\{q_k^{-1}:b_{K-1}+2\le k\le d-1\}
\cup\bigcup_{j=1}^{K-2}\{q_{k+1}/q_k:b_j+2\le k\le b_{j+1}-2\}.
\]

Every member exceeds one by (3)–(4). Let (eta=\min\mathcal G), using (eta=2) if the set is empty. Choose a positive dyadic integer (M) large enough that

\[
\tau=(M+1)/M,\qquad \tau^8<\eta.
\]

Let

\[
S=\max_j\{q_{b_j}^{-1},q_{b_j+1}\},\quad
T_0=\max\{3M2^d,\lceil\tau^4S\rceil+1\},\quad
\varepsilon_0=(6Md)^{-1}.
\tag{5}
\]

This is a finite exact algorithm using only the profile. No numerical root isolation occurs.

To prove a coefficient bound, consider a size-k subset with block counts a_j, and let a_j^* denote its top-filling counts. The ratio of base monomials is

\[
\prod_j\rho_j^{a_j-a_j^*}
=\prod_{j<K}(\rho_j/\rho_{j+1})^{\sum_{i\le j}(a_i-a_i^*)}.
\]

Every exponent is nonpositive, because top filling maximizes each prefix count; unless all counts agree, at least one exponent is a negative integer. Thus every nondominant base monomial is at most (D_k/T_0). There are at most (2^d) subsets and all products of weights lie between (1) and ((1+\varepsilon_0)^d). For every (k), including the endpoints,

\[
1\le E_k:=\frac{e_k}{t_kD_k}
\le(1+\varepsilon_0)^d(1+2^d/T_0)<\tau.
\tag{6}
\]

For the last inequality, ((1+\varepsilon_0)^d\le(1-d\varepsilon_0)^{-1}\le1+(3M)^{-1}), and (2^d/T_0\le(3M)^{-1}); the product of these two bounds is less than (1+M^{-1}).

The circuit error is

\[
\frac{C_k}{C_k^{\rm dom}}=
\frac{E_k^3 E_{k-2}}{E_{k-1}^3 E_{k+1}}\in(\tau^{-4},\tau^4).
\tag{7}
\]

The ratio of errors in two adjacent circuits lies in ((\tau^{-8},\tau^8)). Thus the strict finite margins preserve the signs in the first and last bounded interiors, and preserve the increasing order of all middle bounded-interior circuits. At a boundary,

\[
C_{b_j}^{\rm dom}=(\rho_j/\rho_{j+1})q_{b_j},\qquad
C_{b_j+1}^{\rm dom}=q_{b_j+1}/(\rho_j/\rho_{j+1}).
\]

The choice (T_0>\tau^4S) makes their actual signs strictly positive and strictly negative respectively.

There are therefore (K-1) positive-to-negative changes across the boundary pairs. Between successive boundary pairs the circuit begins negative, has a strictly increasing bounded-interior portion, and ends positive, producing exactly one further change. This remains true for empty or one-element interiors and when an interior circuit is zero. The first and last portions have no changes. The total is ((K-1)+(K-2)=2K-3).

## 4. A typed factorial-row consumer

Let (F_j) be any monic polynomial of degree (m_j\ge2) with strictly negative real roots, written (F_j(n)=\prod_i(n+r_{j,i})). Let

\[
B_j=1+\max\{\lvert\text{nonleading coefficient of }F_j\rvert\}.
\]

The elementary Cauchy bound gives (0<r_{j,i}<B_j): beyond (1+A), where (A) is the coefficient maximum, the leading power dominates the geometric sum of the other powers. Choose (L_j\ge B_j/\varepsilon_0), any bases as in (1), and (c_j=\rho_j/L_j). Then

\[
G_j(n)=c_j^{m_j}F_j(n/c_j+L_j)
\]

has root parameters (\rho_j(1+r_{j,i}/L_j)\in[\rho_j,(1+\varepsilon_0)\rho_j]). Consequently (prod_jG_j) satisfies the theorem. Rational inputs give an explicit rational target polynomial, using coefficient arithmetic alone. A positive leading coefficient other than one is removed first.

This applies to the actual negative-root factorial PF rows supplied by THM-3079. The source checks the literal monic factors

\[
F_m(n)=\sum_{j=0}^m\frac{m!\binom mj}{j!}n^j,
\qquad m=2,3,4,
\]

independently certifies their negative roots by exact Sturm counting, and constructs the target by polynomial substitution, without numerical roots. The resulting degree-nine circuit word is (+,-,+,+,-,-,-), with exactly three changes.

The source-target map is **translate each factor, dilate each factor independently, then multiply**. It preserves strictly negative real roots and positive coefficients. It destroys the original coefficient row, its normalization, and its original circuit ratios; the sidecar retained is each degree and Cauchy coefficient bound. Thus this is a genuine explicit new class built from factorial rows, not a theorem that an unchanged factorial row, Laurent return, or first-gap core has the claimed sign pattern. The cheapest decisive test is the literal degree-nine substitution above.

## 5. Hostiles, controls, and reproduction

Three exact hostiles delimit the theorem. Removing the multiplicity restriction allows the separated roots ((100,10,10,1)), interpreted as profile ((1,2,1)), whose circuit has only one change instead of three. Removing (K\ge2) admits a narrow single block ((1,1,1001/1000,1001/1000)) with a change. The inherited classifier hostile ((1,1,3,3,8)) has signs (-,+,-). These do not contradict (1); they identify the first assumptions that fail. An unspecified, possibly wide notion of “cluster” supplies none of the estimates (6)–(7).

The [source](../../04-computation/continuing8_20260906_newton_clusters.py) checks 15 declared profiles from two to six blocks, each with a coalesced geometric control and a distinct-root, unequal-gap control. It reconstructs exact elementary coefficients, all errors in (6), every boundary sign and interior order, and the exact reversal identity (C_k^*=C_{d+1-k}^{-1}). The reciprocal profile ((3,3,3)) exercises the zero convention: its coalesced word is (++-0+--), whereas a spread control has (++-++--); both have three changes. The factorial transport is additional, not a substitution for a root-bank control.

There are **703 always-active gates**. [Certificate](continuing8_20260906_newton_clusters_certificate.json), [normal output](continuing8_20260906_newton_clusters.out), and [optimized output](continuing8_20260906_newton_clusters_optimized.out) are frozen; the two outputs agree as raw LF bytes.

```text
python continuing8_20260906_newton_clusters.py
python -O continuing8_20260906_newton_clusters.py
```

SHA-256 source `29c4fdc5614380b11e1e69770fc141e831e27f35618b681aae0ced3bb439775e`; output `3d980a902fd23299521e6a3d51b3878d787760853ebcefbe1cf09d7668282281`; certificate `1d35470321fbc24289b62fb49955cafe9e10e7518f26b9ae17f59ee12145e17e`.

The next unresolved mathematical boundary is arbitrary separation or wide internal root spread. The present mechanism requires the quantified narrowness needed to preserve the factorial curvature margins. It supplies no proof of the larger finite-evidence conjecture by relabelling a wide factor as one cluster.
