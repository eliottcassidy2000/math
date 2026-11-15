# Complete residue reach does not rescue a switching guard grammar

**Mechanism.** The dyadic clock theorem controls the depth of one guard
tower. Perhaps a succession of different starting hosts can spend its
digits as absolute choices. In the following grammar we can change the
host between *every* two wraps, with no restriction on the number of
small strong components used at any stage. A growth-cost theorem shows
why that still cannot suffice.

Fix a finite collection of initial tournaments. The allowed steps are
clean two-guard wrapping $G$, as in [the clock](GUARD_DEPTH.md), and
ordinally adjoining a singleton or any strong tournament on at most six
vertices, on either side. Equivalently, one may insert arbitrary ordinal
sums of these small components at each stage. The resulting spectrum
contains only $O(X^s)$ distinct counts at most $X$, where
$$
 s=-\log_2(253/500)=0.9827907099\ldots<1.
$$
This is an upper bound, not a claim about the optimal exponent.

## Proof by a weighted address budget

For an ordinal sum, $P_{A\to B}(t)=P_A(t)P_B(t)$, since each slot
has its A-portion then its B-portion. Thus order and side of placement
can be forgotten between two wraps. In any such interval only a
multiplicity multiset of slot polynomials matters. By odd-cycle packs,
for orders 3 and 4 the polynomial has leading coefficient 1 and at
most one other coefficient; for order 5 or 6 it is
$t^n+c_1t^{n-2}+c_2t^{n-4}$. Thus $H=P(1)$ and $Q=P(2)$
uniquely identify it at any of these four fixed orders (at most two
coefficients with distinct weights in the second equation).
The [literal profile replay](replay_menu.cpp) gives numbers
1, 1, 5 and 16 at orders 3, 4, 5 and 6 respectively.
There is also exactly one polynomial of order 1. Put
$$
 A(z)=\frac{1}{(1-z)(1-z^3)(1-z^4)(1-z^5)^5(1-z^6)^{16}}.
$$
This enumerates available multisets by order, possibly counting the
same polynomial product more than once.

Canonicalize the addresses by stripping the terminal singleton
components applied *after* the last wrap; they do not change $H$.
A prefix ending at that last wrap has weighted order cost $w$, where
each wrap costs 2 and each core insertion costs its order; the fixed
initial host's order $n$ is excluded. The last wrap alone gives
$H\ge1+2^{n+w-2}$, regardless of what was placed before it. A prefix
weight $2^{-sw}$ therefore bounds its inverse count to power $s$,
up to a constant depending on (n,s). One wrap with its preceding
multiset has weight generating function $z^2A(z)$. A whole nonzero
wrap prefix consequently has finite total inverse-count mass if
$z=2^{-s}$ and $z^2A(z)<1$. At the rational point (253/500), direct
integer comparison gives $z^2 A(z)<989/1000$.

It remains to restore the terminal nontrivial ordinal components,
and also allow an address with no wrap. Each nontrivial terminal
component multiplies $H$ by its own $H_i\ge3$. Summing over every
terminal multiplicity gives the finite factor
$\prod_i(1-H_i^{-s})^{-1}$, by the finite menu. This includes
arbitrarily many terminal components. Therefore the total normalized
address Dirichlet mass is finite at $s$, and the number of distinct
values at most $X$ is at most this mass times $X^s$.

With just singleton padding instead, the same test is
$z^2/(1-z)<1$. Elementary Fibonacci address counting then gives
the simpler bound $O(X^{\log_2\varphi})$, including the boundary
exponent: prefix addresses contain at most $\log_2 X+O(1)$ steps
when padding costs 1 and wrapping costs 2.

## What has to change

This is the same kind of bookkeeping as a Kraft bound, but the weights
here are *charged by a physical count formula*, not chosen for
convenience. Even allowing different hosts at every wrap leaves a
sublinear number of absolute baselines. And this grammar contains the
six towers that hit every jointly compatible dyadic/ternary residue
in [the clock construction](GUARD_DEPTH.md). Local control and
absolute scarcity hold in the **same** switching family. The desired
replacement rule needs to improve the growth cost, beat the profile
entropy bound by changing the seam, or both. No density conclusion for
all tournaments follows.
