# First descent and first slope descent agree through a polynomial horizon

Use the ordinary shortcut Collatz map on positive integers: even x goes to x/2 and odd x goes to (3x+1)/2. Write x_0=y, and S_k for the number of odd steps among the first k steps. Let sigma(y) be the first k≥1 with x_k<y, and kappa(y) the first k≥1 with 3^S_k<2^k (infinity if absent).

**Proposition.** There are absolute effective c>0 and A>0 such that if kappa(y)=k is finite and y>k^(A+1)/(3c), then sigma(y)=k. In particular a positive value future minimum y has no coefficient descent for its first floor((c y)^(1/(A+1))) steps (decreasing c if necessary).

The constants may be taken from any effective two-logarithm bound: for 0≤s≤k and k≥1,
```
if k log 2 − s log 3 > 0, then k log 2 − s log 3 ≥ c k^(−A).
```
This follows from [Matveev's linear-form bound](https://www.mathnet.ru/eng/im314) for the fixed positive rationals 2,3; the case s=0 is immediate. This application only uses the existence of such constants, not an optimized irrationality measure.

**Proof.** At every j<k, the coefficient is at least 1 and the affine carry is nonnegative. Thus x_j≥y. At the first coefficient drop the exact logarithmic telescoping identity gives
```
log(x_k/y)
 = −(k log 2 − S_k log 3) + sum_{i<k, x_i odd} log(1 + 1/(3x_i))
 ≤ −c k^(−A) + k/(3y).
```
Under the claimed inequality the right side is strictly negative. All previous x_j≥y, so both stopping times are k. A value future minimum has x_k≥y for every k, giving the horizon claim.

**Scope and contrast.** This is not equality of the two stopping times at all sources. It forces a polynomially long interval of exact slope protection starting at *every* value future minimum, independently of reciprocal summability or the behavior later in the orbit. Compare the [September spine theorem](../collatz_posets_dags_20260927_spine_descent_tree.md), which distinguishes actual and coefficient minima, and [Session 306's logarithmic buffer](../session306_relative_isolators_20261113/ESCAPE_BUFFER.md). Under the inherited THM-4476 summability and any tail minimum y of a divergent orbit, the full future slope deficit is O(y^(h−1+eps)), with h=h₂(log₃2); this separate bound is uniform over all future lengths but does not give an exact slope inequality at every length. The elementary argument above exploits the *first* contraction and costs only k/y. No finite-depth parity certificate or unconditional Collatz conclusion follows.

## Uniform affine-threshold interpretation

Within any fixed length-k parity cylinder whose first contraction occurs at k, write the prefix as x_k=(3^s y+B)/2^k with B≥0. The first descent occurs exactly at k if and only if y>B/(2^k−3^s). The proposition bounds this transition height by k^(A+1)/(3c) uniformly over all such cylinders. Indeed for this stronger assertion one can apply the same real affine trajectory to all real y>0: earlier prefixes have nonnegative carries and slopes ≥1, and each odd letter is by definition the affine operation (3x+1)/2, so the same logarithmic identity applies without integer integrality or parity on those real arguments. This gives a polynomial source cutoff, not a finite residue test for the conjecture.

For comparison, replacing the linear-form bound with the elementary inequality 2^k−3^s≥1 only produces an exponential threshold. The real logarithmic gap and the positive carry must both be controlled; sign and parity frequency alone do not measure when the actual map descends.

**Provenance refinement.** [THM-4512](../../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md) (see also its [source paper](../collatz_boundary_compiler_20261004.md)) already gives the affine cylinder threshold, the estimate B≤j·3^(j−1) before a first odd-accelerated drop, and an explicit application of Matveev's theorem to 2^A/3^j−1. These combine to an analogous polynomial threshold in the odd-accelerated word length j. The point here is the uniform horizon at *every actual value future minimum* in the shortcut coding, with the simple telescoping proof and the following refinement.

## A sharper joint horizon from orbit thinness

Put h=h₂(log₃2). For every a>−0.986211… let M(t) be the count of visits at most t by the fixed divergent orbit. The repository's [THM-4499](../../../01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md) gives M(t)≤C_a t^h (log t)^a with an orbit-independent constant. At a value future minimum y, all future values are ≥y. Integration by parts (or dyadic shells) therefore gives
```
sum_{i≥0} 1/x_i ≤ C'_a y^(h−1) (log y)^a,  y≥2.
```
For dyadic shells the ratio is bounded by a geometrically summable sequence times a polynomial in the shell index. At the first coefficient drop, the logarithmic additive error is thus at most the smaller of k/(3y) and a constant times y^(h−1)(log y)^a. Combining with the two-logarithm separation proves protection through any k satisfying
```
c k^(−A) > min( k/(3y), C''_a y^(h−1)(log y)^a ).
```
The second branch uses inherited thinness, not mere convergence. It also bounds how far below zero any future slope discrepancy can reach, uniformly over all time horizons. This does not exclude a first drop at longer times.
