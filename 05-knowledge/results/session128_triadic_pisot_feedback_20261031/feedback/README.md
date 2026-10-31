# Exact gain absorption and the native Amanda feedback carrier

Status: **PROVED** statements below, with finite exact regression controls.
Session 128, 2026-10-31. No claim of general nonlinear-feedback convergence.

## Source and inheritance

The supplied Charmley--Konrad paper, *The Konrad backbone in an Amanda loop*,
v1.0, 22 December 2025, DOI 10.6084/m9.figshare.30934796, pp.15--20,
Eqs.(57)--(83), defines the source, Poisson solve, sign audit, and threshold
memory. Its native feedback is **not linear**: only its memory filter is linear.
The paper's Eq.(60) has a factor-two inconsistency between its exponential
and cosine expressions. This report chooses

    K(theta) = (1/2) cos(3 theta) - (1/4) cos(6 theta),

which agrees with Eq.(69) and Proposition 3.2. The alternative global factor
2 changes the field amplitude, not any sign statistic below.

Closest inherited mechanism: THM-3305,
`01-canon/theorems/THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md`,
read at repository commit `03a9a2c8839e62880cfe84c74ccfb0a2413d2f8d`
through the GitHub file connector. Its Section 4 uses an exact Woodbury
correction in a static graph. We transfer only that classical algebra, not
its graph, observer, FC/GMC/LRC conclusions. Closest hostile: unique static
solution does not imply convergent fixed-point iteration. Corrected near miss:
the source calls its two-harmonic Galerkin solve an approximation, but correctly
notes on p.17 that it is exact for its prescribed source. Least-used sidecar:
positive amplitude scale is indispensable for reconstructing a field after
passing to the sign-feedback quotient.

Concept board: (1) mean-zero Poisson inversion, (2) positive-scale quotient,
(3) three-shift sign census, (4) threshold cells with endpoint ownership,
(5) exact linear resolvent, (6) commuting-field realization from the companion
triad lane. None is identified with another without its indicated map.

## 1. Positive-gain absorption of the actual nonlinear feedback

Write `q_j=w_(j-1) p_j^2`, fix event geometry theta_j, and suppose the gain
`g_t=kappa r^2 rho*_t` is strictly positive. From Eq.(69), the native field is

    phi_t(theta) = g_t F_sigma,t(theta),
    F_sigma,t(theta) = sum_(j in S_t) q_j sigma_j(t)
       [cos(3(theta-theta_j))/18 - cos(6(theta-theta_j))/144].       (1)

Here the event schedule and weights may depend on time, but must be identical
in the two systems being compared. The audit uses ordinary `sgn(0)=0`:

    Pi_t = (1/(2 pi)) integral sgn(phi'_t(theta)) dtheta.

Since `sgn(g_t x)=sgn(x)` for every `g_t>0`, Pi_t is exactly independent of the
entire positive gain sequence. If two engines have identical event schedules,
weights, geometry, memory parameter lambda, initial memory, and initial signs,
but arbitrary different positive gains, induction in Eqs.(82)--(83) proves
that **their complete memory and chirality trajectories are identical**.
There is no small-gain assumption and no limit argument.

In particular, positive Fibonacci, positive Tribonacci, and positive blended
ratios produce the same parity-controller trajectory when those ratios affect
only Eq.(57). They may produce different field amplitudes and different
backbone/adapter audits. This is a quotient/semiconjugacy of the feedback
observable, not an equivalence of full fields or full engine states.

Failure boundary: zero gain collapses a nonzero field to zero; negative gain
reverses Pi; nonuniform event rescaling changes shape; changing the placement
angle or firing schedule changes the controller. A positive gain tending to
zero still has the same Pi at each finite time, so convergence of fields to
zero does not imply convergence of sign audits to zero.

## 2. Sharp universal parity bound and frozen events

**Theorem.** For every real nonzero trigonometric polynomial g supported on
Fourier frequencies `{+3,-3,+6,-6}`,

    abs((1/(2 pi)) integral sgn(g(theta)) dtheta) <= 1/3.          (2)

**Proof.** The three values at theta, theta+2pi/9, theta+4pi/9 sum to zero,
because each supported frequency has a nontrivial third root of unity as its
shift multiplier. Outside a finite set of theta, none of the three values
vanishes; their signs therefore contain both + and -, so their sign sum is
+1 or -1. Integrating and using translation invariance gives (2). If g is
identically zero, its ordinary sign integral is zero and (2) remains true.
The bound is sharp: `g(theta)=cos(3theta)-cos(6theta)` is positive for two thirds
of the circle and negative for one third, since with x=3theta and u=cos(x),
`g=(1-u)(2u+1)`. Its sign mean is exactly 1/3. Negation gives -1/3. QED.

This sharpness is for the complete four-dimensional harmonic carrier. It does
not assert that the particular default prime/golden-weight event configuration
attains equality.

**Corollary.** With `b_0=0`, `0<lambda<=1`, and the native two-harmonic kernel,
all memories lie in `[-1/3,1/3]`, for arbitrary finite firing schedules. This
follows by induction from `b_(t+1)=(1-lambda)b_t+lambda Pi_t`.
Every event whose seed h_j is greater than 1/3 is always assigned +; every seed
less than -1/3 is always assigned -. These conclusions apply to either
`h_j=sin(theta_j)` or the optional `h_j=sin(3theta_j)`.
At equality, the source's one-sided tie convention must be retained.

For an irrational placement angle, the natural density of indices that can
possibly change sign is at most `2 arcsin(1/3)/pi`; this is the length fraction
of `|sin(theta)|<=1/3`. Irrational rotations are equidistributed: every nonzero
Fourier average is a geometric sum divided by N and tends to zero, then
trigonometric approximation bounds interval indicators. The same fraction
holds for the triadic seeds. This is a necessary condition for flipping, not
a claim that every eligible event flips.

### Cyclic cancellation and a deliberately modified symmetry guarantee

More generally, let a measure-preserving transformation tau satisfy tau^q=id,
q>=2, and suppose `sum_(j=0)^(q-1) g(tau^j theta)=0` almost everywhere. With
ordinary `sgn(0)=0`, each orbit sample either consists entirely of zeros or
contains both a positive and a negative value. Thus its sign sum has absolute
value at most q-2. Integration gives

    abs(mean sign(g)) <= (q-2)/q.                               (2a)

The proof does not require a finite zero set. The q=2 case forces exact zero
parity. Consequently a **proposed change**, not the paper's original kernel,
to any kernel supported only on odd multiples of frequency 3, for example
modes 3 and 9 instead of 3 and 6, gives

    phi'(theta+pi/3)=-phi'(theta),   hence Pi=0

for every finite set of translated impulses and all real amplitudes. Poisson
inversion preserves those modes. The native memory then remains b_t=0 from
its prescribed initial state, without any corrective work. This guarantees
the specific J2 sign-parity statistic by imposed half-turn antisymmetry; it
is a different engine, not a proof that the original 3-and-6 feedback works.

## 3. Fixed-schedule feedback is an exact finite piecewise-affine map

Assume now that S, q_j, theta_j, and the choice of seed h_j are fixed. List the
**distinct** thresholds as `h_1<...<h_m`; events with equal threshold flip as
one group. Use the source's `sgn_+(h-b)` convention. The correct cells are

    I_0=(-infinity,h_1],
    I_i=(h_i,h_(i+1)]  for 1<=i<m,
    I_m=(h_m,infinity).                                         (3)

Let sigma^(i) be the sign pattern in cell I_i, and let c_i be the exact sign
integral of the corresponding normalized derivative in (1). Then
`|c_i|<=1/3`. For all t>=1, regardless of the original initialization of the
signs, the native feedback obeys exactly

    b_(t+1)=T(b_t),
    T(b)=a b + lambda c_i for b in I_i,   a=1-lambda.             (4)

Thus at most m+1 threshold patterns occur, not 2^|S|. Only cells meeting the
invariant memory interval `[-1/3,1/3]` can be visited after the first step.
The reduction discards the common positive gain but preserves every future
memory, chirality assignment, and parity audit. To reconstruct the potential,
retain gain, geometry, event weights and the field formula (1).

### Complete branch-fixed classification

On a cell I_i, a fixed point exists **iff** `c_i in I_i`, and then it is exactly
`b=c_i`. For any starting b in that cell, the whole forward orbit remains in
the same cell and

    T^n(b)=c_i+a^n(b-c_i).                                      (5)

Indeed (4) is a convex combination of b and c_i, and the half-open interval
is convex. Interior fixed points are locally attracting from both sides;
a fixed point at the included right endpoint is guaranteed attracting only
from its owning cell. For lambda=1, (5) is understood for n>=1.

If `c_i` lies strictly outside the closure of I_i, every orbit starting in
that cell exits it after finitely many steps. Otherwise (5), valid for as
long as the branch is followed, would converge outside its closure.

There is one crucial excluded-endpoint case: if `c_i=h_i` is the **excluded
left endpoint** of a finite-lower-end cell, and 0<lambda<1, every branch orbit
stays in I_i while converging to h_i. That limit need not be a fixed point of
the actual discontinuous map, because h_i belongs to the previous cell.
For lambda=1 the branch jumps to h_i and leaves immediately. Thus bounded
branch contraction does not justify dropping the endpoint ownership.

### All proposed periodic words have an exact test

For a proposed cyclic word `(i_0,...,i_(p-1))`, p>=1, the only possible
periodic orbit with this itinerary starts at

    b_0 = lambda * sum_(k=0)^(p-1) a^(p-1-k) c_(i_k)
                         / (1-a^p).                            (6)

Generate the next p states by (4). The word is realizable **iff** each b_k
belongs to its stated half-open cell I_(i_k). This is an exact finite test,
including tie ownership; no perturbative expansion is involved. The word
may have smaller minimal period. A realizable orbit avoiding all thresholds
is locally attracting, with p-step multiplier a^p.

For lambda=1, T has finite image `{c_i}`. Every orbit is therefore eventually
periodic after entering this finite set; its classification is the functional
graph `i -> cell(c_i)`. For 0<lambda<1 we make **no claim** that all orbits are
periodic or converge to fixed points.

A hostile for the generic piecewise-affine template (not asserted to be a
realizable native Amanda instance): on b<=0 take c=1/3 and on b>0 take c=-1/3.
There is no fixed point and, for every 0<lambda<=1, there is the attracting
2-cycle `+-lambda/[3(2-lambda)]`. A second abstract hostile takes c=1/4 on
b<=0 and c=0 on b>0, lambda=1/2: positive starts tend to zero, but T(0)=1/8.
The companion paper-audit lane supplies a separate actual native two-event
persistent-bias witness, not merely this abstract control example.

### Actual native persistent-bias witness (independent paper-audit lane)

The companion `../paper_audit/PAPER_AUDIT.md` and exact-rational interval
certificate `../paper_audit/verify_persistent_bias.py` use the source's
canonical golden angle, weights `w_k=varphi^(-k-2)`, S={1,2}, and initial
signs (+,-). They certify h_1>1/3, h_2<-1/3 and

    -0.20898121054519704 < Pi < -0.20898121054519694.

Our invariant-memory theorem freezes these signs for every 0<lambda<=1.
Consequently Pi is constant and nonzero forever, while
`b_t=Pi*(1-(1-lambda)^t)`. This is a genuine native instance, unlike the
abstract two-cycle controls above. It refutes an unconditional reading of
the paper's claim that the explicit feedback updates enforce long-run J2
symmetry. The source calls J1--J3 required but specifies no hard
reject/retry/stop rule. Adding such a rule could reject this trajectory, not
correct it or prove eventual arrival at an accepted state.
The conclusion holds for every strictly positive gain sequence, including
the positive backbone/adapter choices; it does not require a gain limit.

## 4. Explicit proposed linear extension: exact finite-rank absorption

This section is a **newly specified linear extension**, not Step F of the
paper. Let H be the real mean-zero L2 space on the circle, let
`L=r^(-2)d^2/dtheta^2` on the mean-zero H2 domain, and let `R=L^(-1)`.
For bounded maps `U:R^d->H`, `V:H->R^d`, and a real d-by-d matrix C, put
`B=UCV`. For arbitrary real gamma, consider

    L phi = f + gamma B phi.                                   (7)

Set `S=I_d-gamma V R U C`. Then (7) is equivalent, without a smallness
assumption, to

    S y=V R f,      phi=R f+gamma R U C y,      y=V phi.         (8)

If det S is nonzero, the unique solution is

    phi=[R+gamma R U C S^(-1) V R]f.                            (9)

If S is singular, solutions exist exactly when VRf is in im S (equivalently
all left nullvectors annihilate VRf), and all solutions are reconstructed from
all solutions of (8). The reconstruction is injective on ker S: if its phi
is zero, applying V gives y=0. At gamma=0 the matrix S is the identity.
Hence the exceptional gains form a finite set, zeros of det S, unless empty.
The dimension of the homogeneous solution space is exactly dim ker S.
These are elementary finite-rank resolvent/Woodbury identities, proved here
by substitution, not claimed as new identities.

**Spectral and stability boundaries.** An invertible static equation is not
convergent naive iteration. On the triadic carrier E with r=1 and B=P_E,
there is a unique solution for every gamma other than -9 and -36. However
`phi_(n+1)=Rf+gamma R P_E phi_n` converges from every initial state iff
`|gamma|<9`. For example gamma=18 is a perfectly nonsingular static problem
but has error multiplier -2 on the third harmonic. At gamma=-9 the third
harmonic forcing must vanish; when it does, arbitrary third-harmonic
homogeneous solutions remain. Neutrality of forcing is only the zero-mode
compatibility condition, not compatibility at feedback resonance.

**Carrier boundary.** If feedback maps E into E and the complementary problem
is uniquely solvable, triadic forcing stays in E. Without that invariant
subspace, the old projection is no longer exact. A minimal leak: on a third
harmonic e and a first harmonic z, take `Le=-9e`, `Lz=-z`,
`Be=z`, `Bz=0`, and f=e. The exact solution for every gamma is
`phi=-e/9+gamma z/9`; naive projection to E loses the z term.
For general blocks with an invertible complementary block D, eliminating it
gives the exact Schur complement, not merely the E-to-E compression.

**Field-linear bridge.** The companion-triad lane proves the additional,
separate statement: if B is K-linear for a fixed multinacci field K, then
`det_Q(I-gamma B)=Norm_(K/Q)(det_K(I-gamma B))` for gamma in K, using its
Q-linear realization. This handles its exact linear absorption via a norm.
It requires the commuting/field-linear sidecar and is not a property of
arbitrary feedback or the nonlinear Step-F thresholds.

## Verification and limitations

Reproduction (from this directory):

    python3 verify_feedback.py > primary.normal.json
    python3 -O verify_feedback.py > primary.optimized.json
    python3 independent_feedback.py > independent.normal.json
    python3 -O independent_feedback.py > independent.optimized.json
    cmp primary.normal.json primary.optimized.json
    cmp independent.normal.json independent.optimized.json

The two independent paths agree on all 12,960 proposed words in the declared
27-template, four-gain abstract universe: 600 admissible words (not minimal
cycle counts) and 108 branch-fixed instances. The independent path solves
simultaneous cyclic equations rather than reusing formula (6), and compares
81 direct Gaussian inverses to their rank-two update. These finite templates
are regression controls for the general proof, not a native event census.

`verify_feedback.py` checks exact symbolic carrier shifts, sharpness
factorization, gain scaling, branch/periodic formulae, hostile controls,
resonance, leakage, and an exact rank-two Woodbury instance.
`independent_feedback.py` separately checks branch and resolvent identities
using only rational arithmetic and elementary matrix elimination, importing
neither the primary checker nor its results. Both run in normal and `-O`
modes with identical outputs. No regression check is a substitute for the
all-parameter proofs above. None proves general nonlinear feedback stability,
zero parity in the native engine, an LRC/JC statement, or a literal equivalence
of the unrelated harmonic and algebraic-recurrence carriers.
