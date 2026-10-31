# Session128: triadic equivalences, Pisot attractors, and exact feedback

**Status:** proved representation and feedback statements, exact certificates,
and explicitly cited Diophantine inputs. The attachment supplies inspiration
and a concrete loop, not a blanket equivalence or convergence theorem.

The strongest new connection is exact: the exceptional fraction56/103 from
Session127 is one of **only three rational gains** that make a specified
Tribonacci feedback resolvent integral. Its arithmetic explains both the
continued-fraction305 phenomenon and a nonperturbative closed-loop operator.

## 1. The feedback law and its complete cubic classification

Let C be multiplication by the positive Tribonacci root beta in the basis
(1,beta,beta^2), where beta^3=beta^2+beta+1. Consider the implicit linear step

    y=Cx+gamma Cy,
    C_gamma=(I−gamma C)^(-1)C.

This is a specified algebraic feedback loop, not the attachment's nonlinear
threshold update. For a reduced rational gain gamma=a/b with b>0, its field
representative is

    lambda_gamma=b beta/(b−a beta).

The [uniform norm gate](resolvent_norm/THEORY.md) proves, for orders k=3,4,5,
that algebraic integrality of lambda_gamma, integer entries of C_gamma and
of the resolvent, and the equation F_k(b,a)=+/-1 are equivalent. Necessity
comes from the exact norm and coprimality:

    Norm(lambda_gamma)=(-1)^(k+1)b^k/F_k(b,a),
    gcd(F_k(b,a),b)=1.

Sufficiency uses the integral adjugate of bI−aC. The cited complete cubic
Thue classification therefore gives **all** rational cubic gains:

| Gain gamma | Exact closed-loop operator | Principal generator |
|---|---|---|
| 0 | C | beta |
| 1/2 | 2C^4 | 2beta^4 |
| 56/103 | −103C^17 | −103beta^17 |

The identities beta^(-3)=2−beta and beta^(-16)=56beta−103 prove the two
nonzero reductions. The last gain lies beyond the Neumann-series radius:
(56/103)beta>1. Nevertheless its algebraic resolvent is exactly the displayed
integer matrix. Static solvability, convergence of a matrix series, and
stability of a time update are different questions.

There is an additional exact observation map:

    constant coefficient of lambda_gamma
      = U_k(gamma)
      = gamma^(k−1)/(1−gamma−...−gamma^k).

At56/103 the value is−323008, the scalar recovered in Session127. It is now
identified as one coordinate of the closed-loop generator. Within this
coprime rational-gain family, integrality of that single coordinate detects
integrality of the entire generator; it does not reconstruct the generator.

For a fixed finite prime set S, the same proof identifies S-integral gains
with the corresponding Thue–Mahler condition on F_k(b,a). Classical
finiteness is cited; no explicit complete quartic/quintic gain list is claimed.

## 2. Three genuine equivalent representations

The [triadic packet](triadic/THEORY.md) constructs explicit equivalences among
fixed-field modules, rational companion dynamics, and bilateral recurrence
streams. Objects, all intertwining morphisms, inverse maps and naturality
are specified. These are fixed-degree additive categories, not a claim that
three vaguely similar mathematical constructions are equivalent.

The observation stream of a field element x is the constant coefficient of
beta^n x. Its first k entries form a unimodular coordinate system. In this
system the rank-two sparse-unit condition becomes the vanishing of the
middle k−2 stream coordinates. The source, the complete observation window,
and the named companion representation therefore retain the same information.
A scalar norm or a single recurrence term does not.

The substitution basis also makes the earlier47 relation precise:

    beta^(-16): power coordinates (−103,56,0),
    substitution coordinates (−47,56,0),
    reversed recurrence-window coordinates (56,0,−103).

These are different, explicitly related conventions. The negative entry
means the vector is not the letter-count vector of a finite positive word.
This is a valid linear representation change, not automatic word realization.

The rational resolvent transformation similarly gives an equivalence between
its corresponding operator categories by rational functional calculus.
Intertwiner spaces are preserved, but a time step is redefined. This does
not make C and C_gamma similar: their spectra differ.

## 3. The arithmetic cost of absorbing feedback

The rational field remains the same because the Möbius transformation is
invertible. The integral order generally changes. At an integral gain,

    [Z[beta] : Z[lambda_gamma]]=b^(k(k−1)),
    [Z[beta] : C_gamma Z[beta]]=b^k.

The cubic order indices are64 and103^6 for the two nonzero gains. These are
not the image indices8 and103^3. After inverting b the two orders agree;
without this localization, integral reversibility and arbitrary integral
categorical equivalence must not be inferred.

The principal generator at56/103 is negative, so it is **not** a positive
Pisot number under the standard definition. Its negative103beta^17 is
positive Pisot, as is2beta^4. Their contracting conjugates have magnitudes
103beta^(-17/2)<1 and2beta^(-2)<1, respectively. The sign change is explicit.

## 4. Pisot attractors and what spectral transport preserves

The [Pisot packet](pisot/README.md) retains the standard substitution
1→12, 2→13, 3→1, its incidence matrix, and its ordered prefix digits.
Writing alpha for a nonreal conjugate of beta, |alpha|=beta^(-1/2)<1, the
three compact internal sets satisfy

    R1=alpha(R1 union R2 union R3),
    R2=alpha R1+1,
    R3=alpha R2+1.

Contraction proves the unique graph-directed compact attractor. This alone
is not a tiling or pure-spectrum theorem.

In real coordinates, a shear S_h fixing the digit1 gives the exact operator
B_h=S_h A S_h^(-1) and attractors R'_i=S_h R_i for every real h. No small-gain
restriction is needed; the contraction metric is transported. The ordered
edge Fourier cocycle has an exact change-of-coordinate and endpoint-gauge
formula, verified independently including nonzero component translations.

The sidecars are essential. No nonzero shear in this family preserves even
the original embedded number field, let alone its order. The arithmetic
module and multiplication must be transported. Component-dependent
translations preserve address equality within each labelled component,
not arbitrary spatial coincidences between different components.

An explicit hostile is the twisted substitution1→21,2→13,3→1. It has the
same incidence matrix and eigenvalues, but different prefix digits obstruct
the stated edge-preserving affine conjugacy. Spectrum alone is not enough.
Also the sparse unit beta^(-16), tiny in the principal embedding, has star
image of modulus beta^8>100, outside the displayed attractor's radius<4.
Norm collapse is not an attractor-membership test.

## 5. Nonunit feedback reveals missing finite places

The positive matrices2M^4 and103M^17 are valid incidence matrices, realized
by repeated substitution words. Their expanding roots are nonunit Pisot
numbers. For N=qM^n, complex contraction q alpha^n has area factor

    delta=q^2 beta^(-n),     Perron root lambda=q beta^n,
    delta lambda=q^3.

Thus [essentially disjoint planar graph pieces](nonunit_feedback/README.md)
cannot have a nonzero finite area vector when q>1: the Perron-weighted area
equation would require delta lambda=1. This is an obstruction, not a computed
area or a proof that a particular attractor tiles.

The standard nonunit Pisot representation space retains the finite places
above2 or103 as well. Their total Haar factor q^(-3) restores the necessary
measure scaling1/lambda. Restoring this balance is necessary arithmetic data,
not by itself a proof of disjointness, tiling, or pure spectrum.

## 6. What the attached Amanda loop actually does

The attached *The Konrad backbone in an Amanda loop* (Charmley and Konrad,
version1.0, December22,2025) uses a different feedback law:

    b_(t+1)=(1−lambda)b_t+lambda Pi_t,
    sigma_j(t+1)=sgn_+(sin(theta_j)−b_(t+1)).

The [paper audit](paper_audit/PAPER_AUDIT.md) separates its correct algebra
from normalization, domain and enforcement claims. Equation60 has a visually
confirmed factor-two mismatch; the later cosine normalization is used here.
That common factor cannot affect the parity counterexample.

The [native feedback theorem](feedback/README.md) proves:

- With relative event data fixed, common positive gains do not change the
  memory/chirality trajectory. A positive Fibonacci/Tribonacci amplitude
  adapter changes scale, not this feedback's sign dynamics.
- The two harmonics3 and6 have zero sum on each three-point rotation orbit.
  With ordinary sign(0)=0 this gives the sharp full-carrier bound |Pi|<=1/3.
- Starting b0=0, memory remains in[-1/3,1/3]; generated signs whose seed
  lies outside that interval cannot change.
- Fixed schedules reduce exactly to a one-dimensional piecewise-affine map
  on threshold cells. Fixed points, boundary conventions and periodic
  itineraries have explicit tests; global convergence is not automatic.

A two-event example using the canonical golden geometry and tier weights
has certified persistent bias

    −0.208981210545197037 < Pi < −0.208981210545196940.

It works for every0<lambda<=1 and positive common gain. The explicit update
therefore does not guarantee vanishing parity. An unspecified hard acceptance
filter could reject this trajectory, but is not a supplied repair mechanism.
A proposed odd-harmonic kernel, such as3 and9, instead guarantees zero parity
by an exact sign-reversing symmetry; it changes the model and its audits.

## 7. Synthesis and limits

There are three distinct operations here: quotienting a positive amplitude
in a nonlinear sign loop; conjugating a graph-directed system with its metric,
digits and arithmetic transported; and solving a specified linear algebraic
feedback loop that changes the time generator. Their hypotheses cannot be
swapped. Exact functors explain what is preserved, and the counterexamples
identify which additional data a claimed connection needs.

The complete cubic arithmetic gain classification is the strongest closed
family in this session. Higher-order gain lists, nonlinear global convergence,
Rauzy tiling, and pure-spectrum conclusions remain outside the proved claims.
Classical sources and prior repository mechanisms are credited; no historical
novelty or general open-conjecture solution is asserted.

[Verification and full-session replay](VERIFICATION.md).
