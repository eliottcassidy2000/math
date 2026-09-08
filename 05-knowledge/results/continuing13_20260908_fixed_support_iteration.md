# A fixed critical orbit keeps the pointed unit fixed while multiplicities grow

**Status: PROVED analytic corollary + FINITE-EXACT controls; INDEPENDENTLY AUDITED.** For one fixed surface `W2`, one fixed source coordinate change, and the fixed shear parameter `lambda=1`, this gives polynomial global submersions whose complete pointed unit Weyl modules are isomorphic while their ambient torsion multiplicities and intrinsic rational-pair degrees tend to infinity. The generic original affine fibres are non-isotrivial punctured rational curves. All mates here are rational with poles; there is no polynomial mate or general Jacobian-conjecture conclusion.

The iteration connection was derived by the parent researcher and independently checked by this report's producer. The fixed-source-coordinate improvement uses a second fixed point of the same cubic. The independent referee supplied the algebraic-integer argument allowing `lambda=1` uniformly.

## 1. Supplier, retained data, and the new operation

The direct supplier is the [global polynomial mutation theorem](continuing13_20260908_quintic_mutation.md), Sections 2--5. For a squarefree polynomial `f` of degree `q>=2`, with `f(0)=0`, write `Q'=f^2, Q(0)=0`. If `a0=f(-2h)!=0` and `lambda*(lambda+a0)!=0`, it gives

```
u=x-h, z=u-2h+u^3 t, A=u f(z),
T=lambda*A+Q(z), G=1/(2A^2), J_(x,t)(T,G)=1.
```

The supplier proves both actual charts on fixed `W2`, every original-source fibre and rational constant field, the complete pure order-two principal parts, all ambient arms, and the intrinsic degree of every rational mate. Its original component/connection suppliers are THM-3412, [Hamiltonian principal-part differential and Pruefer torsion arms](../../01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md), and THM-3770, [vertical principal-part equalizer and log-canonical dressing gate](../../01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md). The [previous hidden-arm example](continuing12_20260908_fixed_unit_hidden_arms.md) and [pair-degree corollary](continuing12_20260908_pair_degree.md) establish the closest fixed-pointed-module mechanism.

The current source is a cubic whose critical point has a finite forward orbit and never returns to itself. The map sends iterates to squarefree square roots of derivatives, and then to the supplier's global first functions. It preserves the full target support while multiplying the number of components mapped to each value. The lost coordinate is critical-point multiplicity within a fixed critical **value**, not an unproved change of target coordinates. The fixed-point choice below also keeps `h` and `lambda` unchanged through the whole sequence.

The live comparison board is: the critical orbit; disjoint preimage levels; squarefree derivative square roots; target-support collisions; pointed versus ambient Weyl modules; rational-pair degree; and affine puncture moduli. The cheap hostile is a returning critical point: it keeps the critical-value set small but makes the derivative square root nonsquarefree, invalidating the submersion supplier.

## 2. One exact cubic and fixed source parameters

Choose complex numbers `zeta,c` with

\[
\zeta^2+\zeta+1=0,\qquad c^2=\zeta-1,
\qquad P(z)=z^3+c. \tag{1}
\]

Thus `c!=0`, `zeta!=1`, and

\[
0\longmapsto c\longmapsto\zeta c\longmapsto\zeta c. \tag{2}
\]

Indeed `P(c)=c(c^2+1)=zeta*c` and `P(zeta*c)=c^3+c=zeta*c`. Neither nonzero orbit value is zero, so zero never returns to zero.

The fixed-point polynomial factors as

\[
P(z)-z=(z-\zeta c)(z^2+\zeta c z-\zeta^2). \tag{3}
\]

Choose either root `eta` of the displayed quadratic, and fix

\[
h=-\eta/2,\qquad\lambda=1. \tag{4}
\]

Both choices work. Its constant coefficient is nonzero and its discriminant `1+3*zeta^2` is nonzero. Substitution of `zeta*c` in the quadratic gives `2-3*zeta^2!=0`. Also `eta!=c`, because `P(c)=zeta*c!=c`. Hence

\[
\eta\ne0,c,\zeta c,\qquad P(\eta)=\eta. \tag{5}
\]

No source parameter in (4) depends on the iteration length.

## 3. All derivative square roots are squarefree

Choose a square root `sigma` of three once. For every integer `k>=2`, put

\[
Q_k(z)=P^k(z)-\zeta c,\qquad
f_k(z)=\sigma^k\prod_{j=0}^{k-1}P^j(z),\qquad
q_k=(3^k-1)/2. \tag{6}
\]

Here `P^0(z)=z`, and `P^j` means compositional iteration. The chain rule gives

\[
Q'_k(z)=3^k\prod_{j=0}^{k-1}P^j(z)^2=f_k(z)^2. \tag{7}
\]

Moreover `deg f_k=q_k`, `f_k(0)=0`, and (2) gives `Q_k(0)=0` for every `k>=2`.

For completeness, set `E_j={rho:P^j(rho)=0}`. These finite sets are pairwise disjoint. If `i<j` and `rho` belonged to both, then `P^(j-i)(0)=0`, contradicting (2). Every `E_j` consists of exactly `3^j` simple points: a multiple zero of `P^j` would, by the chain rule, have some earlier iterate equal to zero, and would again force a positive iterate of zero to be zero. Consequently the factors in (6) are squarefree and pairwise coprime, so `f_k` is squarefree.

Its zero at zero has nonzero first derivative

\[
f'_k(0)=\sigma^k c(\zeta c)^{k-2}\ne0. \tag{8}
\]

At the fixed source point `eta`,

\[
a_{0,k}=f_k(\eta)=(\sigma\eta)^k\ne0. \tag{9}
\]

It also never equals minus one, so the fixed choice `lambda=1` is admissible for every `k`. To prove this without a numerical bound, note that `c` satisfies the monic integer polynomial `c^4+3c^2+3=0`. Thus `c` is an algebraic integer, and `eta`, satisfying `eta^3-eta+c=0`, is an algebraic integer as well. If `(sigma*eta)^k=-1`, then `eta^(2k)=1/3^k`. The left side is an algebraic integer, while the rational number on the right is not an integer. A rational algebraic integer is an integer, a contradiction.

This pays all supplier hypotheses uniformly, including the global source-boundary condition that a rational Jacobian identity alone would not guarantee.

## 4. The three supports and every component multiplicity

For `rho in E_j` with `j<k`,

`P^k(rho)=P^(k-j)(0)`.

If `j=k-1` this is `c`; otherwise it is `zeta*c`. Therefore the complete critical-value set of `Q_k` is exactly

\[
\{0,\,a\},\qquad a=(1-\zeta)c. \tag{10}
\]

The corresponding numbers of distinct roots of `f_k` are

\[
e_0=\sum_{j=0}^{k-2}3^j=(3^{k-1}-1)/2,
\qquad e_a=3^{k-1}. \tag{11}
\]

The source component `Eu={u=0}` has the fixed target value

\[
b=Q_k(\eta)=\eta-\zeta c. \tag{12}
\]

The three values `0,a,b` are pairwise distinct by (5): `a!=0`, `b!=0`, and `a-b=c-eta!=0`. There is exactly one additional `A=0` component at `b`, namely `Eu`. Thus for every `k` the **same three target values** support ambient torsion, with respective arm counts

\[
\boxed{(e_0,e_a,e_b)=((3^{k-1}-1)/2,\,3^{k-1},\,1).} \tag{13}
\]

These are `A=0` component counts and ambient torsion arm counts. Each complete special fibre has one further regular component, so its total irreducible-component count is `e_c+1`. The regular component, on which the primitive has principal part zero, is part of the supplier's exactness proof and must not be omitted.

There are useful literal polynomial checks behind these counts. Writing `R_j=product_(i=0)^(j-1) P^i`,

`Q_k-a=(P^(k-1))^3`,

and `Q_k` is divisible by `R_(k-1)^3`, with a quotient of degree `q_k+2` coprime to `R_(k-1)`. These are exact cubic critical multiplicities, consistent with `Q'_k=f_k^2` and squarefree `f_k`. The critical **points** grow exponentially while their two critical **values** stay fixed.

## 5. Global maps, identical pointed unit, and unbounded geometry

Use the fixed coordinates from (4), and for each `k>=2` set

\[
u=x-h,\quad z=u-2h+u^3t,\quad A_k=u f_k(z),
\quad T_k=A_k+Q_k(z),\quad G_k=1/(2A_k^2). \tag{14}
\]

By the supplier and the verified hypotheses, `T_k` is a polynomial global submersion on the same `W2`, and `J_(x,t)(T_k,G_k)=1` in the original rational function field. Its actual source `t`-degree is exactly `3^k`, with leading coefficient `u^(3*3^k)`, because `Q_k` is monic of degree `3^k` and the `A_k` term has smaller `t`-degree.

At every root component of `f_k`, the complete scalar principal part of `G_k` is `1/[2(T_k-c)^2]`; at `Eu` it is `(1+a_(0,k))^2/[2(T_k-b)^2]`. There is no simple-pole coefficient. At the regular component the principal part is zero. Both displayed coefficients are nonzero by Section 3.

Let `tau` denote one common abstract target coordinate, with the action `tau=T_k` in each module. The canonical pointed unit Weyl module is, for **every** `k`, one principal-part tower at each of `tau=0,a,b`, with distinguished vector of pure order two at each support. Rescaling each nonzero coefficient vector identifies these pointed modules while preserving both Weyl actions. This is an isomorphism of the unit-generated pointed modules, not an identification of their differing ambient embeddings.

In particular the exact scalar annihilator is the same ideal for the whole sequence:

\[
\boxed{\operatorname{Ann}_{\mathbf C[\tau]}([1])
=\bigl([\tau(\tau-a)(\tau-b)]^2\bigr).} \tag{15}
\]

Its degree is six, with local order two at each of three supports. The unit is nonzero and generates exactly three full arms. The ambient torsion has

\[
\boxed{q_k+1=(3^k+1)/2\text{ arms},} \tag{16}
\]

so its excess over the unit-generated part tends to infinity despite completely fixed target supports, source parameter `h`, and shear parameter `lambda`.

The supplier also gives

\[
\boxed{[\mathbf C(x,t):\mathbf C(T_k,G_k)]=2\cdot3^k.} \tag{17}
\]

Its paid constant field is `ker D_(T_k)=C(T_k)`. Every rational mate with nonzero constant Jacobian is `alpha*G_k+R(T_k)` with `alpha!=0`, and has the same embedded pair field. Thus (17) is intrinsic to `T_k` among all rational mates; it cannot be removed by choosing a different rational primitive.

For target values outside the same fixed set `{0,a,b}`, the complete original-source fibre is `P1` minus

\[
\boxed{3q_k+2=(3^{k+1}+1)/2\text{ points}.} \tag{18}
\]

For each fixed `k`, these affine curves form a non-isotrivial family as the target value varies, by the supplier's three-fixed-puncture argument. Their smooth projective completions all have genus zero. Across `k`, even their number of punctures is unbounded. The same pointed unit therefore misses ambient arm multiplicities, intrinsic rational-pair degree, and this increasing affine complexity.

## 6. Boundaries, alternate source choice, and exact controls

The no-return condition in (2) is load-bearing. With `c=0`, the cubic is `P=z^3` and already the product for `k=2` is `z*P(z)=z^4`; it is not squarefree. A finite critical-value set by itself does not authorize the submersion or module conclusion. For `k=1`, the particular shift in (6) gives `Q_1(0)=(1-zeta)c!=0`, so it fails the displayed integral normalization. The stated family begins at `k=2`.

One can instead choose `z0` from the next level `E_k` and set `h=-z0/2`. This gives the parent's original fixed support `{0,(1-zeta)c,-zeta*c}`, with the same counts and invariants, but with `h` depending on `k`. All such `z0` are admissible by disjointness of the preimage levels. The fixed choice `lambda=1` remains valid: the product `B=product_(j<k)P^j(z0)` is an algebraic integer, and `sigma^k B=-1` would force `B^2=1/3^k`. The fixed-point version (4) is stronger because it keeps all source parameters unchanged.

The [standalone exact producer](../../04-computation/continuing13_20260908_fixed_support_iteration.py), [transcript](continuing13_20260908_fixed_support_iteration.out), and [certificate](continuing13_20260908_fixed_support_iteration_certificate.json) pass **115 always-active exact gates per mode** using the exact field `Q[c]/(c^4+3c^2+3)`, whose defining polynomial is Eisenstein at three. They test `k=2,3,4`, with source degrees `9,27,81`, ambient arms `5,14,41`, rational-pair degrees `18,54,162`, and puncture counts `14,41,122`. The square identity is checked as `Q'_k=3^k R_k^2`, so no numerical choice of `sqrt(3)` or `eta` is used. Normal and optimized Python runs have byte-identical raw LF stdout and regenerated certificates. A source under `04-computation` writes its certificate to `05-knowledge/results`; an outside source writes beside itself.

The controls include the complete orbit, simple and pairwise-disjoint levels, both critical groups with their exact cubic multiplicity, the original `E_k` variant, fixed-point remainders for the uniform `h`, and the `k=1` and critical-return hostiles. They are finite checks of the corollary's polynomial input; the all-`k` module and geometry claims follow from the analytic proof and the independently paid global mutation theorem. No source-engine import or numerical-root inference is used.

## Independent acceptance

The [independent referee](continuing13_20260908_fixed_support_iteration_audit.md) accepts the all-`k` fixed-`h`, `lambda=1` corollary without mathematical repair. It separately checks the fixed-point factorization, the uniform algebraic-integer exclusion, simple and disjoint preimage levels, all three fixed supports and their complete component counts, the pointed-unit identification, every-rational-mate degree, and non-isotriviality of the punctured fibres. Its analytic proof supplies an independent path; the finite engine remains a small exact control bank rather than a substitute for that proof.
