# Parabolic clocks and a sharp cotangent integrability barrier

19 October 2026. **PROVED:** the statements proved below. **CITED:** the explicitly named one-dimensional analytic inputs. **FINITE-EXACT:** the companion controls. There is no historical-novelty claim, no polynomial Keller counterexample, and no resolution of a Jacobian or coloring conjecture.

## 0. Inheritance and the actual new question

The closest earlier mechanism is [quadratic_clock_return_20261012a.md](../quadratic_clock_return_20261012a.md), especially Sections 2 and 5: polynomial cotangent lifts retain covering sheets and critical divisors, and have no nonconstant polynomial invariant. Its critical-value compiler concerns polynomial regularity of observables, not invariant functions. The new Section 3 below extends that last invariant statement to the entire rational function field by retaining the leading Laurent exponent at infinity.

The newer [planar_autonomous_lie_jets_20261019](../planar_autonomous_lie_jets_20261019/README.md) checkpoint separates autonomous polynomial Hamiltonians from arbitrary formal logarithms. The [metabelian_global_invariants_20261019](../metabelian_global_invariants_20261019/README.md) checkpoint detects rational logarithms by translation-orbit sums of principal parts. The present question is different: when does the formal first integral of a cotangent lift exist analytically, and what does a scalar bottleneck clock retain about it?

The canonical hostile is a rational determinant-one covering with a removed critical divisor. A corrected near miss would identify its formal Hamiltonian, a sectorial Fatou clock, and a single-valued analytic Hamiltonian near the fixed point. These are three different objects.

The exact connection is:

    scalar germ f → cotangent germ C_f(q,p)=(f(q),p/f′(q));
    iterative logarithm v → Hamiltonian H=pv(q);
    Fatou coordinate Phi → symplectic coordinates (Phi(q),p/Phi′(q)).

Composition and the canonical symplectic form are preserved. Global inverse sheets, affine poles, and sector-to-sector transition data are not removed by this lift. The necessary sidecars are the domain f′≠0, the function-field degree, and the Fatou transition derivatives.

## 1. Formal logarithm and complete formal invariant ring

Let K have characteristic zero and let

    f(q)=q+a q^r+O(q^(r+1)),   r≥2, a≠0,

be a formal one-variable germ. Let v=itlog(f) be its normalized iterative logarithm: the unique order-r formal vector-field coefficient for which the time-one map of v(q) partial_q is f. Equivalently, v=a q^r+… and

    v(f(q))=f′(q)v(q).                                      (1)

Set C_f=(f,p/f′), a formal map tangent to identity in (q,p). Use the standard bracket {A,B}=A_qB_p−A_pB_q and D_H={−,H}.

**Theorem 1.** The unique normalized formal Hamiltonian logarithm of C_f is

    H(q,p)=p v(q).                                          (2)

Moreover its entire formal invariant ring is

    K[[q,p]]^(C_f)=K[[H]].                                  (3)

Here K[[H]] means all formal compositions h(H), h∈K[[T]], not just polynomials in H.

**Proof.** Cotangent lifting respects composition by the chain rule. If f_t is the formal flow of v, then C_(f_t) has infinitesimal generator

    v(q) partial_q − p v′(q) partial_p = D_(pv).

Therefore its time-one map is C_f. Filtration-raising exponential and logarithm are inverse, giving uniqueness of (2). This also proves preservation of dq∧dp; directly, C_f preserves the canonical one-form p dq.

Write I=Σ_(m≥0) A_m(q)p^m. Since C_f is tangent to identity, I is fixed precisely when D_H I=0. Its coefficient equations are

    v A_m′−m v′ A_m=0.

In K((q)) this says (A_m/v^m)′=0, hence A_m=c_m v^m. Characteristic zero makes the constants of formal differentiation exactly K. Conversely every h(H) is fixed. Since ord(H)=r+1≥3, all such compositions are defined. This proves (3). ∎

A useful related lemma is that the only B∈K((q)) satisfying B(f)=B are constants. After removing a constant leading term, the first term b q^s has s≠0, and B(f)−B has first term s a b q^(s+r−1), which is nonzero. Consequently every nonzero Laurent solution of

    A(f)=f′^m A,   m∈Z,

is exactly A=c v^m.

## 2. Meromorphic first integral if and only if analytic flow embedding

Now take K=C and f a nonidentity holomorphic germ tangent to identity.

**Theorem 2.** The following are equivalent:

1. C_f has a nonconstant meromorphic first-integral germ at (0,0).
2. C_f has a nonconstant holomorphic first-integral germ there.
3. v=itlog(f) converges near zero.
4. f is the time-one map of a holomorphic vector field vanishing to order at least two.
5. The normalized formal Hamiltonian logarithm pv converges near (0,0).

**Proof.** Conditions 3–5 are equivalent by the flow construction and logarithm uniqueness. Under 3, H=pv is a holomorphic nonconstant invariant, so 3⇒2⇒1.

For 1⇒3, expand a meromorphic invariant I in p at zero, over the field of meromorphic germs in q. A quotient of two convergent two-variable power series has such a Laurent expansion: factor the first nonzero p-power from the denominator, then invert its nonzero meromorphic coefficient. There are finitely many negative powers, and every coefficient is meromorphic in q.

Its leading nonzero term A_m(q)p^m satisfies A_m(f)=f′^m A_m. If m=0, the preceding lemma makes A_0 constant; subtract that constant. Since I is nonconstant, a new nonzero coefficient then appears. We thus obtain some nonzero integer m and a nonzero meromorphic germ A_m=c v^m. If m<0, invert A_m. Hence a positive integer power of v is meromorphic. Its known positive formal order makes it holomorphic, and after factoring q^(r|m|), an analytic root of its nonvanishing unit exists. The root matching the formal leading coefficient is v, proving convergence. ∎

This theorem is an exact transfer between one-dimensional embedding and two-dimensional meromorphic integrability. It does not assume the first integral is linear in p. A nonlinear or quotient-valued proposed invariant cannot evade the leading-coefficient argument.

**CITED analytic consequence.** Aschenbrenner–Bergweiler, *Julia’s equation and differential transcendence*, Theorem 1 and the discussion before it, show that iterative logarithms of nonlinear entire maps are differentially transcendental over entire coefficients; the paper also records the classical divergence result for globally meromorphic maps except parabolic Möbius maps. Their Theorem 2 explicitly covers f(q)=q+q² with the stronger coefficient field of convergent germs. Thus a nonlinear polynomial f tangent to identity has divergent v. By Theorem 2, its cotangent lift has no nonconstant meromorphic first-integral germ, although (3) supplies a full nontrivial formal invariant ring. Primary text: https://www.mat.univie.ac.at/~maschenbrenner/pdf/itlog-final.pdf , pages 3–4.

The rational and meromorphic obstructions persist over R: a real-analytic or real-meromorphic germ complexifies. Merely adjoining complex constants does not repair the obstruction. This differs from the field-dependent allowed-factor phenomena in [metabelian_field_sensitive_jets_20261019](../metabelian_field_sensitive_jets_20261019/README.md).

## 3. A simpler global theorem: no rational invariant for any nonlinear polynomial lift

The following statement needs neither a parabolic fixed point nor analytic dynamics.

**Theorem 3.** For every characteristic-zero field K and polynomial f∈K[q] of degree d≥2,

    K(q,p)^(C_f)=K.                                         (4)

The notation denotes rational functions I for which I(f,p/f′)=I as a rational identity; C_f need not be birational.

**Proof.** Expand I at p=0 over K(q). If its leading exponent is zero, the leading coefficient A_0 satisfies A_0(f)=A_0. A nonconstant rational map A_0 has degree deg(A_0∘f)=d deg A_0, so A_0 must be constant. Subtract it. Unless I is constant, one obtains a nonzero integer m and nonzero A∈K(q) satisfying

    A(f)=f′^m A.                                           (5)

At infinity write f(q)~c q^d and A(q)~B q^s, where c,B≠0 and s∈Z may be negative. Comparing leading exponents in (5) gives

    ds=m(d−1)+s,   hence s=m.

Comparing leading coefficients now gives Bc^m=B(dc)^m, or d^m=1. This is impossible for a nonzero integer m in characteristic zero. ∎

This is the all-rational extension of the earlier polynomial-invariant calculation in [quadratic_clock_return_20261012a.md](../quadratic_clock_return_20261012a.md), not a claim that the cotangent construction itself is new. The earlier critical-value observable ring concerns a different question and remains valid.

There are sharp positive and hostile controls. For the parabolic Möbius map M(q)=q/(1−q), C_M is exactly the repository source flow

    (q,p)↦(q/(1−q),p(1−q)²),

with polynomial Hamiltonian and invariant pq². Conversely, for f(q)=q^d the rational function pq is multiplied by 1/d, not preserved. For a Chebyshev polynomial T_d, (q²−1)p² is multiplied by 1/d². Rational eigenobservables are not rational first integrals.

## 4. The cusp clock: formal residue and coalescing complex poles

For f(q)=q+q², exact formal calculation gives

    v(q)=q²−q³+(3/2)q⁴−(8/3)q⁵+(31/6)q⁶+…,
    1/v(q)=q^(−2)+q^(−1)−1/2+(2/3)q−(13/12)q²+… .

Consequently a formal Abel coordinate, normalized by its constant, is

    Phi(q)=−1/q+log q−q/2+q²/3−13q³/36+… .                 (6)

Its logarithmic coefficient is the residue 1 of the time differential dq/v. Equation Phi(f)=Phi+1 follows by integrating the time differential along the formal flow. The logarithm is sectorial/multivalued; (6) is not a globally meromorphic function at zero.

For epsilon=delta²>0 the corrected real clock used in the companion escape analysis is

    Psi_delta(q)=atan(q/delta)/delta + (1/2)log(q²+delta²),
    omega_delta=dPsi_delta=(1+q)dq/(q²+delta²).               (7)

It is the exact time form of the rational model vector field (q²+delta²)/(1+q), not the exact global logarithm of the discrete recurrence. The model agrees with the exact cusp logarithm through q³ at delta=0; their q⁴ coefficients already differ.

The complex residues at the two poles are

    Res_(i delta) omega_delta = 1/2+1/(2i delta),
    Res_(−i delta) omega_delta = 1/2−1/(2i delta).

Thus the positively oriented loop periods are

    pi/delta+pi i,     −pi/delta+pi i.                       (8)

The poles merge to the double pole with residue 1 in (6). Their summed residue is exactly 1. Along a fixed symmetric real interval [−a,a], the logarithms in (7) cancel and the arctangent passage is pi/delta−2/a+O(delta²). This explains the leading bottleneck constant as a degenerating complex period; an integer escape count still needs endpoint and discrete-error control.

There is also an exact local discrete period. The fixed points of f_delta(q)=q+q²+delta² are ±i delta, with multipliers 1±2i delta. They are repelling when delta>0. Koenigs linearization near each fixed point makes the local map z↦lambda z. Choosing the branch log(lambda)→0 as delta→0, the local time differential has residue 1/log(lambda). At the upper fixed point its loop period is

    2pi i/log(1+2i delta)
      =pi/delta+pi i+(pi/3)delta−(pi i/3)delta²+O(delta³).    (9)

This exact local period and the model period agree in the divergent and constant terms, but not in the next term. Equation (9) does not assert a single exact Fatou coordinate across the whole real passage, or identify its O(delta) term with the entire real escape correction.

For a general simple parabolic germ g(u)=u+a u²+b u³+…, a≠0,

    itlog(g)=a u²+(b−a²)u³+…,
    Res_0(du/itlog(g))=1−b/a².                             (10)

This residue is invariant under formal changes of the scalar coordinate because the time differential transforms by pullback. It can distinguish germs with the same leading quadratic bottleneck type.

## 5. Exact new period-three calculation: residue 51/49

Let f_c(x)=x²+c with c=−7/4 and

    P(x)=x³+x²/2−9x/4−1/8.

The earlier period-three and cotangent-return notes establish

    f_c³(x)−x=(f_c(x)−x)P(x)²,

with P having three distinct real roots, none fixed by f_c. For any one of these roots r, set g_r(u)=f_c³(r+u)−r. Its derivative is one. Reduction modulo P gives

    g_r(u)=u+a_r u²+b_r u³+O(u⁴),
    a_r=14r²−28r−21/2,
    b_r=−56r²+56r,
    a_r²=1372r(r−1).                                      (11)

Neither r=0 nor r=1 is a root of P, so a_r≠0. Therefore

    b_r/a_r²=−2/49,
    Res_0(du/itlog(g_r))=51/49.                            (12)

All identities hold in Q[r]/(P), so (12) holds at all three orbit points without choosing a numerical root. In the normalized coordinate z=a_r u,

    g_r:z↦z+z²−(2/49)z³+…,
    itlog(g_r)=z²−(51/49)z³+…,
    Phi_r(z)=−1/z+(51/49)log z+O(z).                       (13)

The clock advances by one per return, which is three original quadratic steps. Mixing those two units would introduce an erroneous factor of three.

The rational number 51/49 survives the cubic irrational orbit coordinates. It differs from the cusp residue 1, so the period-three return and q+q² are not formally conjugate as time-one scalar germs. This is stronger than observing that both have multiplier one, but it neither decides their full analytic conjugacy moduli nor asserts a new arithmetic classification of three-cycles. Theorem 3 applies to g_r over any characteristic-zero field containing r; Theorem 2 and the cited divergence theorem show that its cotangent return has no meromorphic invariant germ.

## 6. Arbitrarily long agreement with the integrable source flow

Let N≥2 and define the polynomial

    f_N(q)=q+q²+…+q^N.

Compare it with the Möbius map M(q)=q/(1−q). Their cotangent lifts have the same total-degree N-jet at (0,0), and their first difference is

    C_(f_N)−C_M
      = (−q^(N+1), (N+1)pq^N) + O((q,p)^(N+2)).           (14)

Indeed f_N−M=−q^(N+1)/(1−q); differentiating and reciprocating gives the second component. C_M has the polynomial logarithm pq² and its analytic invariant. Every C_(f_N) has no nonconstant rational invariant by Theorem 3 and no nonconstant meromorphic invariant germ by Theorem 2 plus the cited polynomial divergence theorem.

Thus no fixed finite jet certifies even meromorphic first-integral existence in this explicit family. This extends the nature of the obstruction visible in the metabelian finite-jet examples. It is a different family: C_(f_N) has generic covering degree N, whereas every element of the metabelian group is birational. Dropping that distinction would be an invalid transfer.

The already-proved symplectic-jet completion theorem also supplies a polynomial symplectic automorphism matching every finite jet of C_(q+q²). If A_N is such an automorphism and N≥2, its first coordinate must have degree at least N+1: otherwise that coordinate would equal q+q² exactly, which cannot be a polynomial coordinate because its differential vanishes along q=−1/2. In particular these completions cannot have a common degree bound. No optimal degree estimate is claimed.

The formal inverse makes the hidden branch explicit:

    q=(sqrt(1+4Q)−1)/2
      =Q−Q²+2Q³−5Q⁴+14Q⁵−… .                            (15)

The coefficients are nonzero alternating Catalan numbers. The other global branch is lost by choosing the germ at q=0. The source involution (q,p)↦(−1−q,−p) has the same cotangent image. The map is undefined on q=−1/2, so this lift does not cover the initial cotangent fibre of a scalar orbit starting at x=0, q=−1/2. The local bottleneck construction begins after restricting to its regular germ.

## 7. Sectorial gluing is the exact missing condition

For a simple parabolic analytic germ, take the canonical univalent attracting and repelling Fatou coordinates Phi_i on overlapping sectors of opening greater than pi, with their common formal asymptotic expansion and Phi_i∘f=Phi_i+1. Such coordinates are unique up to additive constants. Their existence and asymptotic expansions are classical; a primary proof and the horn-map formulation are in Dudko–Sauzin, *The resurgent character of the Fatou coordinates of a simple parabolic germ*, Definition 1, Theorem 2 and Section 2: https://arxiv.org/pdf/1307.8093 . Arbitrary Abel solutions on narrow sectors, which may admit periodic perturbations, are not the coordinate choices in this statement.

Define on each corresponding cotangent sector

    t_i=Phi_i(q),   E_i=p/Phi_i′(q).

Then dt_i∧dE_i=dq∧dp and C_f becomes (t_i,E_i)↦(t_i+1,E_i). If h_ji=Phi_j∘Phi_i^(-1), the overlap map is exactly

    (t_i,E_i)↦(h_ji(t_i), E_i/h_ji′(t_i)),
    h_ji(t+1)=h_ji(t)+1.                                  (16)

The normalized sectorial Hamiltonians E_i agree on overlaps precisely when all h_ji′=1, that is, the Fatou transitions are translations. In that case the sectorial v_i=1/Phi_i′ glue to a punctured-neighborhood analytic vector field with v_i=O(q²); the singularity is removable and the glued vector field is the convergent v. Conversely a convergent v provides its primitive time coordinates, whose overlaps differ by constants. Thus (16) identifies the gluing obstruction in Theorem 2 concretely.

For the cusp the +log q term creates a constant additive monodromy 2pi i. Such a constant is allowed by the derivative condition. The obstruction to analytic Hamiltonian gluing is the nonconstant part of the horn transitions, not the logarithmic residue alone. In particular matching the residue is necessary for formal conjugacy but does not establish analytic integrability.

The metabelian orbit-sum test and the [unit-distance stress-period test](../unit_distance_stress_periods_20261019/README.md) likewise concern solvability of specified difference equations after local conditions. But their coefficient objects and operators differ. There is no established map transporting a finite-field stress syndrome into a horn map, and no unit-distance or coloring consequence is inferred here. The exact transferable lesson is to state the actual operator, its obstruction space, and the category in which a primitive must exist; the proofs above do that independently.

## 8. Verification, provenance, and stopping point

Run `python3 verify_parabolic_geometry.py`, optionally with `-O`. The script uses exact SymPy arithmetic, always-active failure exceptions, two independent formal-log coefficient constructions, Julia equations, Hamiltonian invariance identities, the period-three quotient calculation, the jet-hostile family, model/discrete residues, and Catalan inverse coefficients. The JSON records the finite universe and the checked coefficients. Neither a finite coefficient list nor a growth plot is offered as a proof of divergence.

Prior content is inherited from:

- [quadratic_clock_return_20261012a.md](../quadratic_clock_return_20261012a.md), cotangent return and the polynomial-invariant obstruction;
- [period3_graph_square_bridge_20261019.md](../period3_graph_square_bridge_20261019.md), the parabolic three-cycle and its real/rational boundary;
- [mobius_affine_source_return_20261011.md](../mobius_affine_source_return_20261011.md), Möbius cotangent source flows;
- [planar_autonomous_lie_jets_20261019/finite_symplectic_jets_and_boundary_control.md](../planar_autonomous_lie_jets_20261019/finite_symplectic_jets_and_boundary_control.md), finite symplectic-jet completion;
- the three metabelian/stress checkpoints named in Sections 0 and 7.

The exact advance is the meromorphic embedding iff criterion, the all-rational field extension of the earlier invariant obstruction, and the period-three residue calculation. The analytic divergence theorem and classical Fatou/horn-map framework are attributed, not reclaimed. The polynomial cotangent family remains globally ramified and outside the polynomial Keller-map class. No claim about arbitrary planar polynomial maps follows.
