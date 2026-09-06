# A conserved amplitude separates global root geometry from return signs

**Status: PROVED ANALYTICALLY + FINITE-EXACT;
[independent analytic and exact audit passed](long_frontier2_sep06_duality_audit.md).**
The regular beta step preserves every coefficient ratio to its genuine row.
Already at height one, an entire exact parameter orbit can keep the genuine
first equations, constant and leading response coefficients, paired
characteristic factors, and strictly negative roots at every parameter,
yet reverse the same-root response sign. The sharp uniform thresholds for
the single interior attenuation are different:

```text
all-parameter strictly negative response roots: lambda >= 1/sqrt(10),
negative response at every genuine first root: lambda >= 13/40.
```

This is a coefficient-model obstruction and a calibrated positive invariant
for that model. It is not an actual Laurent counterexample or an all-height
sign theorem.

## 1. Inheritance and the selected operation

Read first the audited
[regular reflection and paired factors](continuing4_20260906_regular_duality.md),
its [independent audit](continuing4_20260906_regular_duality_audit.md),
the [degree-eight geometry decoder](continuing4_20260906_moments_packet.md),
and the [continuing synthesis](continuing4_20260906_synthesis.md).
The decoder makes lower-moment relaxations an obsolete target: the two
ordinary fifth Hankels exactly recover the prescribed model. Its fixed
coefficient fibres are compact intervals. That separate finite-phase
problem is not explored here.

The closest mechanism is reflection of the discarded carried complement
into an actual regular factorial family. The canonical hostile is the
[carried beta-step model](long_frontier_sep06_allh.md), which keeps its
two endpoint amplitudes but changes the interior ones. The corrected near
miss is blaming that failure on a lost carry, or on root reality failing
after a few further steps. The least-used sidecar is a coefficient ratio
to the exact factorial orbit. The live concepts are reflection, beta
integration, conserved amplitudes, global root geometry, paired factors,
and the original first-root predicate.

The source-to-target map is the exact regular parameter step. It preserves
the genuine coefficient recurrence but conserves every initial amplitude
defect as well. It does not identify an orbit with the actual factorial
one. A bounded search of the relevant Laurent/duality reports and canon
recovered no prior statement of the thresholds or full-orbit obstruction
below; this is a retrieval statement, not a literature-priority claim.

## 2. The regular beta step conserves all amplitude ratios

For H>=1 and z>=1, use the incoming monic rows

```text
Phi_(H,z)(t)=sum_(j=0)^H H! (z+2j+1)^[2H-2j]/[j!(3H-3j)!] t^j,
Psi_(H,z)(t)=sum_(j=0)^(2H) (2H)! (2z+2j+1)^[4H-2j]/[j!(6H-3j)!] t^j.
```

Products are rising, and Theta=t d/dt. Direct coefficient ratios give

```text
(z+1+2Theta) Phi_(H,z+1)=(z+2H+1) Phi_(H,z),
(2z+1+2Theta)(2z+2+2Theta) Psi_(H,z+1)
  =(2z+4H+1)(2z+4H+2) Psi_(H,z).                       (1)
```

Equivalently the first step integrates `u^z Phi_(H,z)(u^2 t)` on [0,1]
with prefactor z+2H+1. The second integrates
`u^(2z)(1-u) Psi_(H,z)(u^2 t)` with prefactor
(2z+4H+1)(2z+4H+2). There is no inverse carry in these complete regular
fibres; its absence follows from the literal count equations, not deletion
from the old carried family. All denominators are positive in scope.

Suppose a real polynomial response family M_z, of degree at most2H and
with no terms outside the genuine support0,...,2H, satisfies the second
recurrence on all integer z>=1. Divide its coefficient of t^j by the nonzero genuine
coefficient of Psi_(H,z). The two recurrence multipliers cancel, so

```text
[t^j]M_z / [t^j]Psi_(H,z) = lambda_j,
```

independent of z. Conversely, any fixed diagonal multiplier lambda_j gives
an exact orbit. This proves an all-height structural statement: the beta
step commutes with every fixed coefficient attenuation. Fixing the constant
and top coefficients sets only lambda_0=lambda_(2H)=1. It leaves 2H-1
conserved interior amplitudes. The genuine orbit requires all of them to
equal one.

## 3. The height-one model and two sharp thresholds

Set H=1 and define, for real z>=1,

```text
a_z=(z+1)(z+2)/6,
A_z=(2z+3)(2z+4),
R_z=(2z+1)(2z+2),
c_z=R_z A_z/360,
Phi_z(t)=t+a_z,
M_z^lambda(t)=t^2+lambda A_z t/3+c_z,  lambda>=0.       (2)
```

The genuine response is lambda=1. Every fixed lambda satisfies the exact
step (1) for all z, and both endpoint coefficients remain genuine.
The displayed real-parameter formula also extends the discrete orbit.

**Global root threshold.** M_z^lambda has two distinct strictly negative
roots for every real z>=1 if and only if lambda>=1/sqrt(10). Indeed its
discriminant is

```text
Delta_z=(A_z^2/9)(lambda^2-R_z/(10A_z)).                (3)
```

Here 0<R_z/A_z<1 because A_z-R_z=8z+10, and the ratio tends to one.
At lambda>=1/sqrt(10), (3) is strictly positive at every finite z. All
three polynomial coefficients are positive, so both roots are negative.
Below that threshold, the discriminant is negative for sufficiently large
z. The same iff holds if z is restricted to positive integers. Equality
at the threshold retains distinct roots at every finite parameter.

**Uniform same-root sign threshold.** Write

```text
ell_z=(13z^2+31z+16)/(10A_z).
```

Direct evaluation at the genuine original root gives

```text
M_z^lambda(-a_z)=(a_z A_z/3)(ell_z-lambda),
13/40-ell_z=(29z+46)/(20A_z)>0,
ell_(z+1)-ell_z=(29z^2+121z+120)/(5A_z A_(z+1))>0.     (4)
```

The supremum of ell_z is 13/40, approached at infinity. Thus
M_z^lambda(-a_z)<0 for every real z>=1 if and only if lambda>=13/40;
again the same statement holds on integer parameters. The boundary value
13/40 gives strict negativity at every finite z. The gap is nonempty since
(13/40)^2-1/10=9/1600>0.

Consequently the exact beta step does preserve a calibrated sign region:
lambda>=13/40 is invariant and sufficient, with a sharp uniform constant.
The weaker globally real-rooted region has the strictly smaller threshold
1/sqrt(10). The genuine lambda=1 satisfies the stronger condition; its
actual height-one sign was already known.

## 4. One global rational orbit reverses sign at primitive adjacent masses

Choose lambda=8/25. Since lambda^2>1/10 and lambda<13/40, every polynomial
M_z has two distinct negative roots throughout the entire positive orbit.
Equation (4) simplifies to

```text
M_z(-a_z)=a_z (z^2-69z-112)/150.                       (5)
```

The quadratic factor has one positive root, strictly between70 and71,
and one negative root. Therefore the response is negative for integers
1<=z<=70 and positive for every integer z>=71. In particular

```text
z=70: a_z=852, M_70(-852)=-5964/25,
z=71: a_z=876, M_71(-876)=876/5.                       (6)
```

The reflected actual supports at these two parameter addresses are
`(-3,70,216)` and `(-3,71,219)`. Their designated first masses are73 and74;
both z values are coprime to3, so both are genuine first support-return
addresses. The model has the genuine first equations at these addresses,
not merely phases chosen independently of the first row.

An exact double-cancellation model is available without losing global root
reality: set lambda_*=ell_70=10981/34320. Then

```text
1/sqrt(10)<lambda_*<13/40,
Phi_70(-852)=M_70^(lambda_*)(-852)=0.                  (7)
```

Its entire parameter orbit is still strictly negatively rooted by (3).
Equations (6)-(7) concern model response polynomials. Changing alpha,beta,
gamma of an actual Laurent polynomial changes its phase and monomial
anchor, but does not change the normalized interior factorial amplitude
from one to8/25 or lambda_*. No actual first/doubled moment cancellation
or sign failure is asserted here.

## 5. Paired characteristic factors survive the obstruction

At H=1, quotient multiplication has one eigenvalue M_z^lambda(-a_z).
Its nonleading characteristic coefficient is the negative of this value.
The incoming paired factor is D_(1,1)=(z+1)(z+2)=6a_z. It divides that
coefficient for every lambda, with exact residual

```text
c_(1,1)^lambda / D_(1,1)
  =[10lambda A_z-(13z^2+31z+16)]/180.                 (8)
```

For lambda=8/25 the residual is `(-z^2+69z+112)/900`. The paired factor
and all beta-step identities therefore survive while residual positivity
fails. At lambda>=13/40, every coefficient of (8) is nonnegative and its
constant and linear coefficients are positive; at equality the quadratic
coefficient vanishes. Thus the sharp sign threshold is also an exact
coefficient-positivity threshold for this height-one residual.

The all-height divisor theorem is unchanged. The obstruction says that its
factors, even together with global root reality and the complete parameter
step, do not supply the missing positivity. The actual interior amplitude
is a conserved coordinate, not an error that more beta steps can remove.

The reflection map itself still excludes the old singular zero block.
At old h=1,x=-1, the generic carried quotient response is -1/90720, while
the specialized raw row is t^2/720 and is zero at t=0. The controls retain
this inherited carry hostile. No model argument above identifies those
two values or removes the original carry before reflection.

## 6. Exact controls and stopping boundary

The [standalone source](../../04-computation/long_frontier2_sep06_duality.py)
and [frozen output](long_frontier2_sep06_duality.out) pass **93 explicit exact
gates**, with byte-identical normal and optimized output. The exact
universe contains formal step and threshold identities, the complete
height-one global orbit, the two primitive sign controls, the exact model
double cancellation, complete literal genuine rows at named parameters,
one nonzero complementary reflection, and the inherited singular-carry
control. Universal statements rest on (1)-(8), not on parameter samples.

The source uses standard-library Fraction arithmetic and imports no research
producer. Symbolic coefficient-ratio checks at heights1,2,6 are controls for
the all-height analytic derivation. Literal first/doubled charge enumeration
at z=1,5,70,71 independently checks all genuine coefficients and positive
multinomial normalization factors. These actual controls preserve the
distinction between lambda=1 and the model amplitudes.

```bash
python3 04-computation/long_frontier2_sep06_duality.py > /tmp/duality-normal.out
python3 -O 04-computation/long_frontier2_sep06_duality.py > /tmp/duality-optimized.out
cmp /tmp/duality-normal.out /tmp/duality-optimized.out
cmp /tmp/duality-normal.out 05-knowledge/results/long_frontier2_sep06_duality.out
```

Source and output are frozen; the root integrator owns their checkpoint.
Independent analytic/source review passed; see the linked audit.

The next structural question at higher height is which conserved amplitude
region guarantees the actual same-root sign, or which operation couples
those amplitudes instead of commuting with every diagonal multiplier.
Another negative-root table or another paired-divisor calculation cannot
settle that missing implication. No higher-height amplitude cone is claimed
by this bounded result.
