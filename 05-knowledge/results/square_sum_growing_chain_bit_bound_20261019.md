# Boundary-start recurrence and singly exponential-bit AP contact chains

**Status: PROVED intermediate, pending independent audit.** The companion controls are **FINITE-EXACT**. Its all-mixed-sign construction remains valid; for compatible arbitrary-sign chains the stronger [fixed-chart cubic theorem](square_sum_cubic_contact_chains_20261019.md) supersedes its quantitative bound. The boundary-start Lemma 2 below remains essential to that improvement. This improves the dependence on a growing chain length; it does not give a polynomial-in-length chain algorithm or a polynomial numerical label bound.

## 1. Statement

Let `T_0=(a,b,C)` be an admissible primitive Pythagorean petal, let `epsilon in {-1,+1}` be its prescribed outgoing port, and let `k>=1`. Put `L=1+ceil(log_2 C)`.

**Theorem 1.** A deterministic exact algorithm constructs distinct admissible primitive petals `T_1,...,T_k`, all different from `T_0`, such that

\[
 1+\varepsilon f_{T_i}-\varepsilon f_{T_{i+1}}\in\mathbb Q^2,
 \qquad 0\le i<k.                                      \tag{1}
\]

At a common integer square scaling the resulting `2k+3` positive integer labels have the inherited AP-paired square-sum Hamiltonian cycle. Every label has

\[
 \boxed{\quad O\!\left(k^4\,176^k\,L^2(L+k)\right)\quad} \tag{2}
\]

bits, with an absolute effective implied constant. The complete algorithm takes `exp(O(k)) poly(L)` bit operations. Numerical labels are consequently at most `exp(exp(O(k)) poly(L))`; no fixed-power numerical bound in `k` is asserted.

The improvement is that the seed-to-small-root transition is performed **once** on the seed's fixed elliptic curve. The rest of the chain iterates one fixed degree-11 map inside a controlled cusp chart. Collision avoidance consumes only `k+1` candidates at the first step, rather than a growing atlas at each new seed.

## 2. Inheritance and the boundary point that changes the construction

We inherit the quartic/cubic isomorphism and port transfer from [square_sum_ppt_bridges_20261019.md](square_sum_ppt_bridges_20261019.md), the real-period bound and signed circle pigeonhole argument from [square_sum_contact_recurrence_bound_20261019.md](square_sum_contact_recurrence_bound_20261019.md), the `r`-point amplification from [square_sum_fixed_count_contacts_20261019.md](square_sum_fixed_count_contacts_20261019.md), and division-polynomial height estimates from [square_sum_contact_output_height_20261019.md](square_sum_contact_output_height_20261019.md). The two degree-11 maps are the audited plus and minus bounded cusp sections.

Write `f_0=2ab/C²`, `k_0=|a+epsilon b|/C`, and `K=k_0²=1+epsilon f_0`. (The subscript on `k_0` distinguishes this rational square root from the requested chain length.) The fixed seed curve is

\[
 z^2=K(u^2+1)^2-4u(u^2-1),\qquad
 E:Y^2=(X-K)(X^2-1).                                    \tag{3}
\]

Instead of starting the translating orbit at the original seed's quartic point, start at the rational **zero-spread boundary**

\[
 (u,z)=(\varepsilon,2k_0).
\]

Under the inherited birational map this becomes

\[
 B_\varepsilon=
 (2K+\varepsilon,\ 2\varepsilon k_0(K+\varepsilon)).       \tag{4}
\]

Substitution verifies (4) directly. It is an ordinary smooth point: `K in (1/2,3/2)`, `K!=1`, and its inverse-map factors `X` and `X-K` are nonzero. Its coordinates have `O(L)` bits. The input petal is not used as a boundary point, and no assertion that its old translate coset approaches the cusp in bounded index is made.

## 3. A quantitative supply of roots in an arbitrarily small cusp arc

**Lemma 2.** For rational `0<eta<=1/64` and integer `r>=1`, one can find `r` pairwise distinct admissible primitive petals `U_j` of spreads

\[
 0<f_{U_j}<\eta,
 \qquad 1+\varepsilon f_0-\varepsilon f_{U_j}\in\mathbb Q^2,
                                                               \tag{5}
\]

by testing even signed translates `B_epsilon+qP`, `P=(0,k_0)`, with

\[
 |q|\le O(r^2L/\eta).                                  \tag{6}
\]

Each resulting rational parameter and certificate has

\[
 O\!\left(Q^2[L+\log(Q+1)]\right),\qquad
 Q=O(r^2L/\eta),                                       \tag{7}
\]

bits. The algorithm is polynomial in `r,L,1/eta`.

**Proof.** On the positive-z low-parameter arc starting at `(epsilon,2k_0)`, increasing target spread `s` from zero to `eta` has invariant length

\[
 \ell=\frac14\int_0^\eta
 \frac{ds}{\sqrt{1-s^2}\sqrt{1+\varepsilon(f_0-s)}}
       >\frac\eta5.                                    \tag{8}
\]

This is a genuine open one-sided arc: the target low parameter starts at `1`, spread is strictly increasing there, and both square-root factors stay nonzero. It is a proper embedded arc in the real component containing the boundary point. Translation by `A=2P` stays in that component and is non-torsion. The inherited `r`-point signed pigeonhole lemma produces `r` points at displacements `alpha,2alpha,...,r alpha` strictly inside the arc, with `0<alpha<ell/r`, by some signed translate `nA` satisfying `n<=floor(r Omega/ell)`. Hence the corresponding indices are `q=2 sigma j n`, with `1<=j<=r`.

For a rational search budget, choose `J>=0` least with `2^J f_0>=5/2`, and put `W=5(4+J)/7`. The inherited period estimate gives `Omega<W`, and `J=O(L)` because `f_0>=2/C²`. Thus one may test

\[
 q\in2\mathbb Z,\quad
 0<|q|\le 2r\left\lfloor\frac{5rW}{\eta}\right\rfloor.   \tag{9}
\]

Filter exactly for ordinary inverse coordinates, `z>0`, `1<epsilon u<4/3`, and distinct spreads in `(0,eta)`. The points supplied by the low monotone arc satisfy these filters. It is not necessary to know their indices or compute a real period.

For completeness, only `ell<Omega` is needed here; no extra period lower bound is imported. Set `N=floor(r Omega/ell)`. The `N+1` distinct circle points `0,A,...,NA` have a circular gap at most `Omega/(N+1)<ell/r`. The difference of its endpoints is `sigma n A`, with `1<=n<=N`, at signed positive displacement `alpha<ell/r`. Its first `r` multiples remain in the chosen embedded arc because `r alpha<ell<Omega`. This proves exactly the form used above, including strictness and absence of wraparound.

The height proof is unchanged when the `O(L)`-height point `B_epsilon` replaces the former `O(L)`-height anchor. Division polynomials bound the coordinates of `qP` by (7); one ordinary addition with `B_epsilon` preserves that order. Identity, inverse, and tangent cases use their usual separate formulas. The inverse quartic map and primitive normalization have bounded degrees. This proves (7) and the stated exact complexity. QED.

The possible presence of the original petal among these candidates is intentionally allowed. It will be removed together with every future collision with that petal in Section 5.

## 4. One fixed map at each port, on a uniformly contracting input chart

Let `F_+(t)` be the parameter of the plus cusp section `3Q+2P`; explicitly,

\[
 F_+(t)=
 \frac{t(t^4-8t^3+2t^2+1)
       (5t^6-12t^5+23t^4-17t^2+12t-3)}
 {(2t^5-t^4-4t^3-2t^2+2t-1)
  (t^6+20t^5-29t^4+16t^3+3t^2-4t+1)}.                  \tag{10}
\]

Let `F_-(t)=1/F_+(1/t)`, the inherited minus section `3Q_-`. Both are rational maps of degree exactly `11`, and each supplies the contact

\[
 1+\varepsilon f(t)-\varepsilon f(F_\varepsilon(t))
          \in\mathbb Q^2.                              \tag{11}
\]

**Lemma 3.** On the entire rational chart `0<h<=1/64`,

\[
 h<F_\varepsilon(1+h)-1<4h,\qquad F_\varepsilon'(1+h)>0.
                                                               \tag{12}
\]

In particular all outputs lie in `(1,17/16)`, where the inherited positivity and contact certificates apply.

**Exact polynomial proof.** Write `F_epsilon(1+h)=N(h)/D(h)` with integer polynomials and `D(0)=32`. Define

\[
 P_1=(N-(1+h)D)/h,\quad
 P_2=((1+4h)D-N)/h,\quad P_3=N'D-ND'.
\]

Their constant terms are respectively `64,32,3072` for both signs. The following upper bounds on the absolute negative-coefficient tails at `h=1/64` are strict:

- plus: `D:1`, `P_1:0` (identically no negative tail), `P_2:1`, `P_3:0`;
- minus: `D:1`, `P_1:1`, `P_2:3`, `P_3:1`.

Each bound is below the corresponding positive constant term; nonnegative terms can be dropped. This proves `D,P_1,P_2,P_3>0` throughout the chart. The companion checker displays the exact rational gaps, verifies divisibility by `h`, and checks the polynomial identities. This is a uniform sign proof, not a sample. QED.

The same fixed maps satisfy the usual rational height inequality

\[
 h_{\rm rat}(F_\varepsilon(t))\le11h_{\rm rat}(t)+O(1), \tag{13}
\]

where `h_rat(a/b)=log max(|a|,|b|)` for coprime integers. It follows simply by evaluating the homogeneous numerator and denominator degree-11 forms and reducing the fraction.

## 5. Build the whole chain without a new atlas at every step

Set

\[
 \eta=\frac{4^{-k}}{64},\qquad r=k+1.                   \tag{14}
\]

Apply Lemma 2 once, obtaining `r` distinct candidate roots, with low parameters `t_j=1+h_j`. On the whole low cone, `f(1+h)>=h`; for instance this follows immediately from the rational formula on `0<h<1/3`. Hence `h_j<eta`.

For each candidate root form the `k` parameters

\[
 t_j,F_\varepsilon(t_j),...,F_\varepsilon^{\,k-1}(t_j).  \tag{15}
\]

By (12), every iterated displacement is at most `4^(k-1)eta<1/64`; every transition uses a valid cusp chart. The parameters and spreads strictly increase down this chain, so these `k` petals are mutually distinct.

For each fixed index `i`, the function `f(F_epsilon^i(t))` is strictly increasing on this root interval. Thus at most one candidate root can make its `i`th petal equal to the original seed's spread `f_0`. There are `k` possible indices and `k+1` distinct candidates. At least one candidate avoids the seed at every index. Test the candidates exactly and retain one that passes. Primitive-petal uniqueness from spread proves all required distinctions.

The seed is joined to this root by (5); every following pair is joined by (11). Every edge uses outgoing port `epsilon` and incoming port `-epsilon`, so consecutive contacts at interior petals use opposite ports. The inherited common-center transfer supplies the asserted cycle. The construction is non-induced: extra square-sum edges are harmless.

## 6. Height and time

From (9) and (14), the largest tested index satisfies

\[
 Q=O(k^2 4^k L).
\]

By (7), a chosen root has rational height

\[
 O(k^4 16^k L^2(L+k)).                                 \tag{16}
\]

Here `log(Q+1)=O(k+L)`; `k,L>=1`. Iterating (13) at most `k-1` times gives per-petal and per-rational-certificate bit lengths

\[
 O(k^4 176^k L^2(L+k)).                                \tag{17}
\]

Primitive PPT normalization has bounded degree and can only change the absolute implied constant. The divisibility proposition below shows that the `k` cusp-orbit hypotenuses are nested, each dividing the next. Therefore their common multiple is the last one. Taking `H=lcm(C,C_last)<=C C_last` adds only `O(L)` bits to (17); no extra factor `k` is necessary. Every final graph label is less than `4H²`, proving (2). The integer edge-sum roots have no larger asymptotic bit length.

Testing all `O(Q)` translates and all `O(k²)` root/iterate pairs uses a number of rational operations bounded by `exp(O(k)) poly(L)`, on operands of the same form. Schoolbook arithmetic therefore proves the running-time claim. This deliberately conservative estimate does not claim a small practical constant.

## 7. Exact prime-valuation transport and an orbit-specific lower obstruction

**Proposition 4.** For either degree-11 cusp map, let `C` and `C'` be the primitive input and output hypotenuses. Whenever both parameters define positive nondegenerate PPTs,

\[
 C'=C^2 R,\qquad R\in\mathbb Z_{>0},\qquad\gcd(C,R)=1. \tag{18}
\]

Thus every already-present hypotenuse-prime valuation doubles exactly; all new prime factors are disjoint from the input hypotenuse. In particular the cusp-chain hypotenuses are nested, and `log C_j>=2^j log C_0`.

**Proof.** For the numerator and denominator `N,D` in (10), exact integer-polynomial identities give

\[
 N^2+D^2=(t^2+1)^2P(t),\qquad
 \operatorname{Res}(N,D)=-2^{102},                     \tag{19}
\]

where `P` has degree 18 and

\[
 P(t)\equiv-49152(t+1)\pmod{t^2+1}.                 \tag{20}
\]

For reference the factorization is

\[
\begin{aligned}
P={}&(29t^8-132t^7+216t^6-108t^5-6t^4+20t^3-4t+1)\
 &\cdot(t^{10}-8t^9+133t^8-64t^7-150t^6+144t^5
          +10t^4-64t^3+37t^2-8t+1).
\end{aligned}
\]

Take coprime integers `n,m` with `t=n/m`. The input primitive hypotenuse is `(n²+m²)/delta`, with `delta in {1,2}`. The homogeneous numerator and denominator have no common odd prime divisor: finite projective roots are excluded by the resultant, and at infinity their leading coefficients `5,2` are relatively prime. Thus fraction reduction and subsequent PPT normalization remove only powers of two from their sum of squares.

Every prime `p|C` is `1 mod 4`. In particular `p!=2,3`, and `m` is a unit modulo `p`. The residue of `t` satisfies `t²=-1`, so `t+1!=0` modulo `p`. Equation (20) says the homogeneous residual factor is a `p`-adic unit. Equation (19) therefore gives exactly `v_p(C')=2v_p(C)`. This proves (18) for plus. Replacing `t` by `1/t`, then reciprocating the output, preserves primitive hypotenuses; the exact duality `F_-=1/F_+(1/t)` proves the minus case. QED.

This is an exact lower obstruction for this particular iterative construction, not for all contact chains. Even with optimal cancellation and lcm sharing, a long orbit of the degree-11 map cannot have polynomial-in-length bit size. The one-step polynomial bound was therefore never enough to remove the long-chain difficulty.

A useful contrast is the odd section `F_5` used in the balanced-tree companion note. Its degree-33 numerator and denominator have resultant `2^968`, and their sum of squares has exactly one factor `t²+1`, with residual congruent to `-2^52` modulo that factor. Thus that section preserves every old hypotenuse-prime valuation exactly, while its remaining factors are new. This does not by itself bound the sizes of the new factors or prove polynomial-size long chains.

## 8. Controls, connection contract, and frontier

The companion script verifies both fixed maps, the strict chart bounds, monotonicity, the boundary point on the cubic and quartic, exact finite chains, and the final integer square-sum cycles. The analytic existence of all boundary-arc candidates comes from the inherited real recurrence proof, not from finite enumeration.

Source: an arbitrary admissible primitive petal and a target chain length. Target: a port-compatible distinct contact chain and its common-center integer cycle. Map: one fixed-seed elliptic boundary recurrence, then iterates of one fixed cusp rational map. Preserved: the prescribed port, rational square contacts, strict distinctness, and AP pairing after common scaling. Lost: economical numerical labels and interval coverage. Sidecars: the cusp radius, monotonicity, a finite set of root exclusions, and the lcm denominator.

The canonical failed route was to repeatedly apply a many-neighbor atlas to ever-growing seeds. The replacement makes all exclusions at the first parameter and uses a monotone one-map orbit thereafter. For a chain, the needed initial arc still shrinks exponentially and the rational-map degree grows as `11^k`. Thus polynomial-in-length bit complexity remains open. The separate common-center capacity theorem already rules out polynomial numerical labels in `k`, even without contacts; it does not rule out polynomial bit complexity.

No graceful labeling of a prescribed tree, Hamiltonian cycle on a consecutive interval, or Collatz termination statement follows.


## Reproduction from the repository root

The exact control requires Python 3 and SymPy. Run `python3 05-knowledge/experiments/square_sum_growing_chain_bit_bound_20261019.py`. The matching receipt is `05-knowledge/results/square_sum_growing_chain_bit_bound_20261019.json`; the script prints a fresh receipt.
