# Scoped corrections to the incoming Moran-family theorem

**PROVED elementary corrections / FINITE-EXACT controls.** Read-only review of
incoming `fdad326e9`, followed by authorized repairs of THM-4504 and its
`procgen_family27_20260926_long_orbit_families.md` proof note. The incoming
2^32 scan and literature review were not repeated. This audit does not
blanket-retract the fair-word or finite-window theorems.

## 1. The conditional mean needs a nonempty conditioning event

For k=1,W=2, both possible values M_1=1/2,3/2 are below W. Thus
P(tau<=k)=0 and eps_k=E[M_k]=1. The conditional expectation
E[M_tau|tau<=k] is undefined, and a claimed strict lower bound0<P is false.
The universally correct identity is

    E[M_tau;tau<=k]=1-eps_k.

Only on positive hit probability may it be divided by the conditional mean.
The infinite-horizon strict bound survives: eps_k tends to zero, and the
first overshoot lies in [W,3W/2). All finite thresholds, including those too
high to hit, are retained by the undivided identity. The note also retains
the n=1 finite exception when translating to trajectories stopped at1.

## 2. Supremum is not pointwise attainment

Start M_0=1. Choose the next odd letter if M_j<2, and the next even letter
otherwise. Then 1<=M_j<3 forever. Setting S_j=log M_j gives exactly rotation
by log(3/2) on the circle of length log3: below log2, add log(3/2);
at or above log2, subtract log2. The rotation ratio is irrational, since a
rational ratio would imply 3^u=2^v for nonzero integers u,v.

An irrational circle rotation has dense forward orbit (pigeonhole gives
arbitrarily small nonzero return steps; their multiples meet every interval).
Consequently sup_j M_j=3, but M_j never equals3. Every finite prefix extends
by all-even letters to a sequence whose entire supremum is strictly below3.
Thus {sup_j M_j>=3} is not open. The parity isometry transfers this example
to a dyadic address; no ordinary positive-integer realization is asserted.

The repair is the open set R_W^hit={exists j:M_j>=W}. Under the fair-word
law, the mean logarithmic increment is log(sqrt3/2)<0, so M_j tends to zero
almost surely. Its supremum is therefore attained almost surely. The hit
and supremum events have the same Haar probability, though they are not
equal sets. The stated infinite-riser dimension argument is unchanged.

## 3. Two scales give two different exponents

The exact window identity is count=2^(k-m)W_m. If
log_2 W_m=h m+o(m) and m/k tends to c>0, then

    log_2(count)/k -> 1-c(1-h).

The covering dimension is h in prefix scale m. At the stated exact-window
edge c=log_3 2, the integer-height exponent is larger than h. Only c=1
makes these numerical exponents equal, and that case needs the separate
THM-4495 result rather than Theorem G's short-window argument. The formal
Theorem S has 1<beta<log_2 3; summaries have been restricted to that proved
domain instead of silently including its upper endpoint.

## 4. What to integrate with the poset and atomic-carrier work

The fair-word martingale and finite complete-residue transfers are valid
probabilistic tools. The branch exponent identity remains conditional on
regular variation and a limiting residue-class fraction. The model delay
constant41.677648 and independent-model record rate2/(n+1) are not actual
Collatz orbit bounds. Fractions near0.3927 refer to the incoming finite
scan; existence of that branch's natural density remains open.

These distinctions sharpen the current source-identity barrier: a Haar-null
set can contain a prescribed integer, whereas vanishing mass under a measure
with positive atoms at every positive integer would exclude each integer.
No such mass-decay estimate follows from the fair-word martingale here.
The threshold3 example supplies a second boundary failure: finite prefixes
can approach a threshold forever without any finite prefix hitting it.

Reproduce the small independent audit with
`python 04-computation/experiments/crossroads_poset_20260926_moran_audit.py`.
Checks use exact fractions and remain active under `-O`. The finite greedy
probe checks the invariant, not the infinite density theorem by simulation.
