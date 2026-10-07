# Overlap depth: integration, a sharper bound, and conditioning repairs

**Status:** PROVED elementary formulas and counterexamples; FINITE-EXACT
controls. The distribution of overlap depths under a useful coarse history,
debt recurrence, and integer Collatz coverage remain OPEN.

## Incoming mechanism and notation

The incoming checkpoints `9dc16269c` and `27d1d417b` supplied
[THM-4564, aligned tape readers](../../01-canon/theorems/THM-4564-two-orbits-read-one-tape-alignment-coupling-of-collatz-exponent-streams.md)
and [THM-4565, readers at any offset](../../01-canon/theorems/THM-4565-two-readers-any-offset-cross-covariance-of-coupled-collatz-exponents-is-the-saturation-kernel.md).
Their commit messages said audit pending. This note independently checks
the elementary saturation kernel and joins it to our
[adaptive-precision law](collatz_lazy_precision_innovations_20261007.md).
It does not promote their empirical depth distributions to theorems.

Use D for their saturation depth, to distinguish it from our M, the number
of revealed source bits. At a pair of current windows, let
lambda=v+A_pre-B_pre. The endpoint precision deficits are
qY=(-lambda)^+ and qX=lambda^+. At least one reader is fresh.

On an overlap, let F be the fresh exponent and R the earlier reader's
exponent minus the offset. Conditional on the complete joint past, D is
known, F is positive geometric, and

\[
R=\min(F,D)\quad(F\ne D),\qquad
R=D+G\quad(F=D),
\]

where the extra G in the tie is fresh geometric. Both F and R have
geometric marginals conditional on this past. This agrees with the lazy
rule: the earlier query either reuses the fresh observation or asks for
additional bits at its critical boundary.

Direct summation gives

\[
\operatorname{Cov}(F,R\mid D)=2-6\,2^{-D},\qquad
\Pr(F=R\mid D)=1-2\,2^{-D},\qquad
I(F;R\mid D)=2(1-2^{-D}).                         \tag{1}
\]

The formulas extend to D=infinity by taking 2^(-D)=0. Nonoverlapping
windows have conditional covariance zero. At every joint past, at least
one of the two raw conditional exponent means is 2. Since the Y exponent
has unconditional mean 2, the conditional-mean term in total covariance
vanishes, including when which reader is fresh depends on the past.
Thus, for the actual overlap event O,

\[
\operatorname{Cov}(a_s,b_t)
=\mathbb E[(2-6\,2^{-D})\mathbf1_O],\qquad
-\Pr(O)\le\operatorname{Cov}(a_s,b_t)\le2\Pr(O).       \tag{2}
\]

## A sharper executable long-gap bound

An overlap of these two windows requires A_s>B_(t-1)-v. Hence our
negative-binomial overlap bound combines with (2) to give

\[
\boxed{\quad
|\operatorname{Cov}(a_s,b_t)|\le
2\inf_{H\in\mathbb Z_{\ge0}}\min\{1,
\Pr(\mathrm{NB}(s)>H)+2^v\Pr(\mathrm{NB}(t-1)<v+H)\}.
\quad}                                                     \tag{3}
\]

Any fixed H gives an exact rational upper bound through the companion
`overlap_bound` function. This improves the earlier Holder square-root
bound: the cost is linear in the overlap probability. It retains the
proved exponential decay for t>=(1+epsilon)s with fixed v, while leaving
the near-diagonal depth law and return problem open.

## Three matching scalar statistics still lose the depth law

First consider the following hidden-depth law, independent of the fresh input,

\[
\Pr(D=1)=1/3,\qquad\Pr(D=2)=2/3.
\]

It has E[2^(-D)]=1/3, exactly the positive-geometric depth value.
Equation (1) therefore gives zero covariance, agreement probability 1/3,
and average **conditional** mutual information 4/3 bits. All three are
functions of the same moment; they are not three independent tests.
Yet P(F=R=1)=1/3 rather than 1/4, and P(F=R=2)=0 rather than 1/16.
The pair is dependent. The next section realizes this very law on actual
Mersenne initial alignments at arbitrarily large gaps.

The full diagonal profile restores the lost information:

\[
\Pr(F=R=j)=2^{-j}\Pr(D>j),\qquad j\ge1.                \tag{4}
\]

Differences recover each finite depth probability, and the limiting tail
recovers any infinity atom. In particular full independence is equivalent
to geometric depth in this kernel family. Unconditional I(F;R) must not
be confused with the average conditional quantity in (1): it is zero
under geometric-depth mixing and positive for the two-depth hostile.

## Corrections carried into the incoming truth surface

**A measurable depth cannot remain Haar after conditioning on itself.**
The depth D is a function of the full joint past G. Therefore
P(D>=2 | G) is a zero-one indicator, differing from the geometric value
1/2 by exactly 1/2. It cannot have an error tending to zero along any
unbounded sequence of gap values. The original full-past formulation of
[HYP-9218](../hypotheses/HYP-9218-asymptotic-haar-alignment-of-two-collatz-orbits.md)
is repaired to a refuted formulation with an OPEN coarse-conditioning
research direction. A replacement must specify a smaller sigma-field
which has not already recorded D. Averaged offset/gap classes and selected
finite older observations are possible objects to test, not automatically
valid hypotheses. The existing numerical tables concern such averages.

Unbounded aligned gaps do occur on explicit positive-probability cylinders.
Take raw y=3 mod4, v=4k, and the paired source U(y), for k>=1. The X orbit
has its first 2k-1 valuations equal to2, then valuation3. Its total is
4k+1, exactly the initial paired scale. At the alignment (s,t)=(0,2k),
kappa=(1-3^(2k))/2 and D=2+v2(k). Thus the depth is fixed on this entire
cylinder at arbitrarily large gaps; its conditional probability given
v=4k is1/2. Under the incoming geometric v law its probability is2^(-4k).
This makes the full-past refutation nonvacuous. It also rules out a uniform
claim of geometric depth at every such initial alignment, even if only
the pair index and its gap are retained. A viable averaged replacement
must specify its sampling and transient exclusions. The family varies v,
so it does not contradict our fixed-v off-diagonal estimate.

In fact the complete coarse law is explicit. Keep v=4k but allow every
odd raw y, and put a0=v2(3y+1), nu=3+v2(k). The first2k-1 X valuations
are2, and the next is2+c with c=v2(1+3^(2k+1)y). The two numerators
defining a0 and c differ, after multiplication by the unit3^(2k), by
1-3^(2k), of valuation nu. Thus c=a0 iff a0<nu, exactly the alignment
condition. On it D=nu-a0. Consequently

\[
\Pr(D=d\mid v=4k,\text{initial alignment})
=\frac{2^d}{2^\nu-2},\qquad 1\le d<\nu.                \tag{5}
\]

For every odd k, nu=3, and (5) is exactly the two-depth law1/3,2/3
above. The fresh endpoint U(y) is Haar independently of a0, so the
saturation mixture applies without additional conditioning losses.
This is an actual infinite family of exponent pairs whose covariance
vanishes although they remain dependent at initial alignment separations2k. It retains the scale
v=4k; discarding that growing scale would hide the failure boundary.

**The lockstep depth has a cancellation exception.** For aligned states
x=3^k y+kappa, put E=3kappa+1-3^k and delta=v2(E). If the shared exponent
a<delta, then

\[
E'=3E/2^a-(3^k-1).
\]

With nu=v2(3^k-1), unequal valuations give delta'=min(delta-a,nu).
At equality, cancellation can give delta'>nu, even arbitrarily large.
There is no unconditional post-first-step cap by nu.

A lawful positive-integer witness in the raw v=2 Mersenne family is
raw Y=19 and X=12*19+1=229. Compare U(19)=29 with X:

\[
29\to11\to17\to13,\qquad
229\to43\to65\to49\to37.
\]

The aligned states y=17, x=49 have k=1, kappa=-2, E=-8, delta=3,
and a=2. The next alignment again has E'=-8 and delta'=3>nu=1.
Every displayed step precedes ROOT. The exact transition survives; the
blanket cap does not. The source note and hypothesis also keep numerical
near-independence and variance estimates separate from universal claims.

The normalized-debt recurrence retains a forcing term 2^(-L_next).
Removing it gives a comparison perpetuity, not automatically the exact
law of the coupled debt. The relevant tail transfer needs a bound on
that discarded forcing.

## Reproduction and next obligation

Run `python 04-computation/experiments/collatz_overlap_kernel_integration_20261007.py`
and repeat with `python -O`. The checker makes **2,685 exact checks**:
kernel moments with full analytic geometric tails, eight literal affine
maps on every odd residue with ten extra bits, the non-Haar mixture,
the full diagonal profile, geometric mixing, and the actual lockstep
counterexample, plus the unbounded-gap family for1<=k<=16 and four source
representatives per cylinder, and complete odd residue classes proving
the coarse-law controls through k16. No empirical table is used to prove a limiting law.

The focused next task is to choose a coarse conditioning that retains the
debt-return obligation while leaving the comparison depth genuinely
unrevealed, then prove quantitative depth-tail bounds under that choice.
The full diagonal profile in (4) is a stronger diagnostic than the three
redundant scalar summaries. Recurrence additionally needs control after
failed merge attempts; neither unconditional mixing nor (3) supplies it.
