# Paying ternary depth: a uniform regular-density contraction and the atom cost

**Continuation:** [finite signed atomic packets](collatz_reset_debt_20261015f.md#5-conductor-moments-cannot-absorb-a-finite-signed-source-packet-for-free)
have the same obstruction even with vanishing moments: their weighted
Wiener norm diverges exactly for `R>=1/3`. W1/W2 below remain valid on their
stated regular-density domains; this does not yield an unpriced native
source bound.

2026-10-15. **PROVED:** a weighted Fourier norm bound at infinite depth, an
explicit marked-atom cost, and the floor-error repair. **FINITE-EXACT:** grading
controls through depth 8. **OPEN:** an applicable source-preserving regularity
bound for actual Collatz rows. This does not prove an energy or ROOT bound.

[Program](../../04-computation/experiments/collatz_conductor_wiener_20261015d.py),
[output](collatz_conductor_wiener_20261015d.json).

## 1. Inherit the stronger finite-depth result, then retain its price

Incoming commit a405c50b8a adds statement 7 to
[THM-4658, two-octave-memory-energy-descent-critical-lines](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md):
the odd branch lowers Fourier conductor, so every finite-depth mean-zero
excess block has only -2^(-rho) times roots of unity as eigenvalues. Its
[HYP-9284](../hypotheses/HYP-9284-collatz-odd-step-suppression-wavelet-energies.md)
keeps all-scale odd suppression open. The main theorem's unweighted operator
norm can exceed one, so finite-depth spectral radii cannot be used as a
uniform norm contraction. This is the closest mechanism and hostile.

The least-used sidecar is the *cost of the initial fine digits*. The
[odd-ray note](collatz_odd_ray_overlap_20261015d.md) shows both a paid first
block and a cyclic potential whose price grows with ternary depth. Here a
weighted norm makes that price explicit. The operator is the incoming density
operator, not an asserted conjugate of a native source controller.

## 2. A weighted Wiener norm at all depths

Let f be a mean-zero function on Z_3 with an absolutely summable Fourier
series over characters of 3-power order. For a nontrivial character xi let
q(xi)>=1 be its exact conductor exponent. For R>=1 define

    ||f||_R = sum_(xi!=0) |fhat(xi)| R^q(xi).

The completed space is a Banach space of continuous mean-zero densities.
For rho with sigma=Re(rho)>0 use THM-4658's excess operator

    Ptilde_rho f(a) = -2^(-rho) f(a/2)
       -(3/2)^rho f((3a+1)/2)/3
       +(3/2)^rho integral_(2+3Z_3) f.                        (1)

**Proposition W1 (PROVED).** Uniformly in depth, including infinite depth,

    ||Ptilde_rho f||_R <= gamma_R(sigma) ||f||_R,
    gamma_R(sigma)=2^(-sigma)+(3/2)^sigma/(3R).              (2)

The operator norm is exactly gamma_R(sigma) on the full space. It contracts
whenever

    R > 3^(sigma-1)/(2^sigma-1).                            (3)

Proof: the halving branch permutes characters within each conductor level.
A character of conductor q>=2 is sent by the odd branch to a single character
of conductor q-1, with phase of modulus one and coefficient (3/2)^rho/3.
Its norm therefore pays a factor 1/R. For q=1 that branch is constant and the
rank-one correction cancels it exactly. Triangle inequalities prove (2) on
finite Fourier polynomials; continuity gives the completion. A single
character of conductor at least 2 has its two output frequencies at different
levels, so their absolute masses add and attain (2). Rearranging gamma<1 gives
(3). There is no exchange of an unbounded-depth limit with a finite matrix
spectral estimate.

At sigma=1/2, R=sqrt(2),

    gamma=1/sqrt(2)+1/sqrt(12)<1,
    gamma^2=7/12+1/sqrt(6)<1,                               (4)

where the last comparison is exactly 24<25 after squaring positive sides.
The threshold R is (1+sqrt(2))/sqrt(3). Thus a mean-zero excess fixed point
at Re(rho)=1/2 cannot belong to this weighted Wiener class, even with
infinitely many nonzero Fourier coefficients. More generally (3) excludes
fixed points with that norm finite. The Neumann series gives
||(I-Ptilde_rho)^(-1)||<=1/(1-gamma) in the same norm.

This extends the incoming finite-depth exclusion to a specified regular
infinite-depth class. It does not exclude arbitrary distributions, nor does
it assert that actual Collatz fluctuations have this regularity. In
particular it does not turn the density operator into a native integer-clock
controller. Its sqrt(2) is a paid norm weight, not an imported Q--Pell
covariance law.

## 3. A marked source does not have bounded cost in that norm

For a marked residue a modulo 3^Q, put

    f_(Q,a)=3^Q 1_(a+3^Q Z_3)-1.

This is the mean-zero excess of a normalized source atom. Every nontrivial
Fourier coefficient through conductor Q has modulus 1. Exactly 2*3^(q-1)
characters have conductor q, hence

    ||f_(Q,a)||_R = 2R ((3R)^Q-1)/(3R-1).                   (5)

The price is independent of the location a but grows as (3R)^Q. A point atom
has infinite norm. Replacing it by a bounded-norm density discards exactly
the fine source address that fixed-source ROOT arguments must retain.
With a uniformly bounded forcing g_j in this norm, (2) does give the usual
geometric forced-recurrence estimate. An unpriced atom cannot be inserted
as such a forcing.

The genuine map here is Fourier expansion of an excess density; it preserves
the linear transfer equation and exposes the conductor filtration. The
lost coordinate in a finite truncation is the fine Fourier tail. Its sidecar
is (5), or a stronger bound on the *actual signed weighted tail*. The cheapest
decisive test is to apply the proposed regularity claim to this marked atom.
This is directly relevant to the residual invoice in the odd-ray note: the
new question is whether source-weight cancellation pays the tail before it
is smoothed, not whether smoothing makes an operator contract.

## 4. Integer floors leave a modal residue

Incoming HYP-9285 and its companion claimed that the exact integer identity

    H(N)=-[S(N)+3S(floor(2N/3))+3S(floor(N/3))]              (6)

annihilates every power N^rho with D(rho)=0 exactly. That inference drops the
floors. For Re(rho)>0, define 0^rho=0. At N=1 the bracket evaluated on the
power is 1, not 0. This is a minimal exact hostile, independent of numerical
root evaluation. At multiples of 3 the cancellation is exact.

The strongest general survivor is, for D(rho)=0 and sigma=Re(rho)>0,

    N^rho+3 floor(2N/3)^rho+3 floor(N/3)^rho
      =3[floor(2N/3)^rho-(2N/3)^rho]
        +3[floor(N/3)^rho-(N/3)^rho]
      =O_rho(N^(sigma-1)).                                (7)

For 0<sigma<=1 and N>=3 an explicit bound is
3|rho|(6^(1-sigma)+3^(1-sigma))N^(sigma-1), by integrating the derivative
on each interval between a real argument and its floor. Thus the leading
power cancels, but the integer-clock remainder survives. A finite packet
of modes inherits the sum of these explicit bounds. An infinite expansion
needs a summability or justified truncation estimate; one cannot sum the
individual O_rho terms without it. The exact identity (6), conditional pole
certificates, finite-depth spectral theorem and numerical fits are unchanged.
This correction is recorded as MISTAKE-620.

## 5. The critical boundary pays one conductor moment

2026-10-15 continuation, clock_duality_oct15. **PROVED.** Incoming commit
e62673a9f7 constructs unbounded-depth eigen-distributions at every
Re(rho)!=0 by successively choosing Fourier-block preimages. This does not
contradict W1: those distributions need not have its summable weighted
Fourier tail. There is also a useful strengthening exactly at the boundary
where W1's strict contraction disappears.

Let 0<sigma=Re(rho)<=1, and set

    R_c=3^(sigma-1)/(2^sigma-1),
    a=2^(-sigma),             c=1-a.

Here R_c>=1: the convex function 2^(-sigma)+(3/2)^sigma/3 lies above its
value 1 at sigma=1 on (0,1), since its derivative at 1 is negative.
Write F_q=R_c^q sum_(conductor(xi)=q) |fhat(xi)|. Thus ||f||_(R_c)=sum F_q.

**Proposition W2 (critical strong stability and moment resolvent).** On the
mean-zero weighted Wiener space at R_c,

    ||Ptilde_rho^m f||_(R_c)
      <= sum_(q>=1) F_q Pr{Binomial(m,c)<q}.              (8)

In particular Ptilde_rho^m f tends to zero in this norm for every f in the
space, even when it has infinitely many conductor levels. There is no
nonzero fixed point there. If

    M_1(g)=sum_(q>=1) q R_c^q sum_(conductor(xi)=q)|ghat(xi)| < infinity,

then the equation (I-Ptilde_rho)f=g has a unique solution in this space and

    ||f||_(R_c) <= M_1(g)/c.                             (9)

Also ||Ptilde_rho^m g||_(R_c)<=M_1(g)/(c*m) for m>=1.

Proof. Dominate each Fourier branch by its absolute coefficient. A halving
letter keeps the conductor and has weighted mass a; an odd letter lowers
the conductor by one and has weighted mass
(3/2)^sigma/(3R_c)=c. At conductor one that odd letter is killed by the
rank-one correction. Thus survival from conductor q for m letters is exactly
the displayed binomial tail in this positive *majorant*. Fourier phases can
only decrease the actual norm; no independence of Collatz trajectories has
been asserted. The binomial variable counts branch letters in the expansion
of a linear operator, not sampled integer orbits.

For each fixed q the tail tends to zero. Dominated convergence applies
because sum F_q is finite, proving strong stability and excluding fixed
points. The expected waiting time for q lowering letters is q/c, so

    sum_(m>=0) Pr{Binomial(m,c)<q}=q/c.

Tonelli's theorem applied to nonnegative majorants yields
sum_m ||Ptilde_rho^m g||<=M_1(g)/c. The absolutely convergent Neumann sum
therefore solves the equation; strong stability proves uniqueness. Markov's
inequality for the same waiting time gives the last rate. QED.

This is an endpoint regularity statement, not a uniform geometric norm
bound or an exclusion of the incoming arbitrary distributions. Its source
price remains explicit. For the marked atom of section 3, the additional
moment is

    M_1(f_(Q,a))=sum_(q=1..Q) 2q*3^(q-1)*R_c^q,

which is unbounded in Q. The native counterpart is the survivor residual
in [P1--P2, native packets](collatz_survivor_packets_20261015e.md): both retain
the mass left over after an allowed reduction, but no source/clock-preserving
intertwiner between these two reductions has been proved.

The [critical-boundary control](../../04-computation/experiments/collatz_critical_conductor_20261015e.py)
checks exact finite branch survival at rho=1,R=1 and exact waiting-time
identities at three rational lowering probabilities. These are checks of
the proof's majorant; they do not certify the unknown regularity of an actual
Collatz main term.

Normal/optimized W2 control output SHA256:
`42364748e11673c6ef9b3e6ab69dd90097ccb675f4b7a372c4bb9a6ee6c337bc`.

## Reproduction

```sh
python3 04-computation/experiments/collatz_conductor_wiener_20261015d.py
python3 -O 04-computation/experiments/collatz_conductor_wiener_20261015d.py
```

The script checks every nonzero Fourier mode at Q=1..8 for the exact conductor
map, the rank-one cancellation phase, and the source-atom formula at R=3/2.
The complex-parameter norm assertion has the proof above; no floating
spectral search is presented as its proof. The normal and optimized outputs
are deterministic.

Normal/optimized output SHA256:
`ec94abe78f95c4cc4b4c15b94b28885902417b96e43328879de51ca76f481f6a`.
