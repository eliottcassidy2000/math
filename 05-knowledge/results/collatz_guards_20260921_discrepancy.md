# Collatz scaling discrepancy: no bounded strip and a genuine Sturmian clock

**Date:** 2026-09-21. **Status:** PROVED elementary statements with FINITE-EXACT
controls. No novelty claim. These are restrictions on possible infinite words,
not a proof that all positive Collatz trajectories converge.

## Inheritance and the coordinates being retained

The closest mechanism is the ordered-carry identity in
[arithmetic braids](arithmetic_braids_20260917_collatz.md). The hostile is that
**every finite positive halving word has a positive integer realization**;
see [inverse completion](arithmetic_braids2_20260917_inverse_completion.md).
The corrected near miss is confusing an infinite compatible sequence of
2-adic cylinders with a single positive ordinary integer. The underused
sidecar is the number of distinct ordinary integers available below a height.

The historical [rapidity-defect reflection](../../07-reflections/collatz-rapidity-defect.md)
already retains the accumulated additive correction in a logarithmic identity.
We use that identity, not its empirical bound or its invalid replacement of
an inverse tree by its immediate predecessor star. Independently,
[HYP-2456, Beatty-Pell crossover](../hypotheses/HYP-2456-beatty-pell-crossover-word.md)
distinguishes a genuine Sturmian *address clock* from a carry-decorated
visible word. That file remains an OPEN synthesis; no unproved classifier
is a dependency here. The transferable operation is to separate the clock
from the additional arithmetic coordinates and prove each formula directly.

The live board has five objects: halving words, ordered carry, ordinary
height, residue-unit density, and irrational mechanical clocks. The anchor
is descent; the niche is the capacity bound; the wildcard is the signed
`an+b` extension. Passing to the clock preserves the multiplicative scaling
and discards the carry and source guard. They must be restored before any
integer-orbit conclusion.

## 1. An exact additive form of the discrepancy identity

Let `n_j=T^j(n_0)>0` be odd, `k_j=v_2(3*n_(j-1)+1)` for `j>=1`, and

```text
K_0=0, K_j=k_1+...+k_j,
alpha=log_2(3), Delta_j=K_j-j*alpha,
q_j=2^(K_j)/3^j=2^(Delta_j).
```

The elementary step relation gives

```text
q_(j+1)*n_(j+1) = q_j*n_j + q_j/3.
```

Therefore, for every finite j,

```text
n_j = q_j^(-1) * (n_0 + (1/3)*sum_(i=0)^(j-1) q_i).       (D1)
```

This is the same ordered carry in another normalization: `B_j/3^j` is
`(1/3)*sum_(i<j)q_i`. In the logarithmic representation the same information
is `log(n_j/n_0)=-Delta_j*log(2)+sum_(i<j)log(1+1/(3*n_i))`.
Deleting the sum changes the dynamics. In particular the leading-factor
condition `2^K_j>3^j` is necessary, but alone not sufficient, for descent.

## 2. Direct capacity bound, before using stopping-time density

**PROVED.** If an entire positive odd Collatz orbit satisfies
`A<=q_j<=B` for all `j>=0`, with `0<A<=B<infinity`, then

```text
B/A >= 9.                                                  (D2)
```

Equivalently, if `a<=Delta_j<=b` for all j, then `b-a>=log_2(9)`.
No assertion of sharpness or existence at equality is made.

First the orbit cannot repeat a state. If it had a cycle of r odd steps
with total halving count s, then along repetitions q would be multiplied
by `2^s/3^r`. This factor cannot equal one, by unique prime factorization.
Its repeated powers contradict either the positive lower bound or the
finite upper bound. Thus every n_j is distinct.

By (D1), for `1<=j<=N`,

```text
n_j <= n_0/A + (B/(3A))*N = M_N.                           (D3)
```

Every image n_j, j>=1, is odd and prime to three. The count of such
positive integers at most an integer M is

```text
floor((M+5)/6)+floor((M+1)/6) <= M/3+1.
```

There must be room for the N distinct images below M_N. Consequently

```text
N <= n_0/(3A) + (B/(9A))*N + 1.
```

Divide by N and let N tend to infinity to obtain (D2). Both coordinates
matter: ordinary growth is bounded using the carry sum, and the congruence
restriction limits available states. A real-valued surrogate without
integer spacing would not satisfy this counting argument.

For a **finite** prefix the same counting inequality holds if distinctness
of `n_1,...,n_N` is assumed separately. If `B/A<9`, it gives
`N<=(3*n_0+9*A)/(9*A-B)`. Bounded q on a finite prefix alone does not
imply distinctness; finitely many repetitions of the state 1 are a hostile.

The argument applies also to a bounded tail, after reindexing and dividing
all q's by their value at the tail start. It does **not** determine the
direction of escape or provide a uniform finite-time descent bound. The
next section strengthens it to exclude every finite strip for `3n+1`;
the direct counting argument remains useful for the more general maps in
section 4. A convergent orbit
repeating the fixed odd state 1 has `q_j` eventually proportional to
`(4/3)^j`, and so is consistent with the theorem.

## 2a. A source-density argument excludes every bounded strip

**PROVED.** No positive integer Collatz orbit satisfies
`0<A<=q_j<=B<infinity` for all j. Equivalently, its discrepancy
`Delta_j=K_j-j*log_2(3)` cannot stay in any bounded real interval,
even after discarding a finite prefix.

The missing ingredient in the first capacity argument is an elementary
fixed-factor stopping-time density fact. For any fixed `epsilon>0`, let

```text
E_epsilon={positive odd n: T^ell(n)>=epsilon*n for every ell>=0}.
```

This set has relative natural density zero among the odd integers.
For epsilon>1 it is empty already at ell=0; the substantive case is
`0<epsilon<=1`.
Here is a self-contained proof; it is a version of the classical
finite-stopping-time method, not a new almost-everywhere convergence
claim. For background, see [Terras's 1976 paper, bibliographic record](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/30/3/101028/a-stopping-time-problem-on-the-positive-integers)
and [Lagarias's exposition of the stopping-time argument](https://www.cecm.sfu.ca/organics/papers/lagarias/paper/html/node4.html).

A specified word `(k_1,...,k_L)` has odd-source density `2^(-sum k_i)`.
Consequently the set `K_L<7L/4` has density equal to the probability of
that event for L independent auxiliary geometric variables with
`Pr(k=t)=2^(-t)`, t>=1. This is an identity of **source counts**:
the event is a finite union of disjoint word cylinders. No probabilistic
assumption is placed on any one orbit. The auxiliary sum has mean 2L
and variance 2L, so Chebyshev gives

```text
density_odd{K_L<7L/4} <= 32/L.                            (D2a)
```

Uniformly over every actual word of length L, its affine intercept is
at most `(3/2)^L-1`. Indeed each individual affine step has slope at
most 3/2 and intercept at most 1/2, and composing these inequalities
gives the stated bound. Take L=4m. On the complement of the bad event,

```text
T^L(n) <= (81/128)^m*n + (3/2)^L-1.                     (D2b)
```

For arbitrarily large m the leading factor is less than epsilon/2.
For each such fixed m, all sufficiently large good sources therefore
have `T^L(n)<epsilon*n`. Thus the upper density of E_epsilon is at
most 32/L. Let L tend to infinity through these multiples of four.
This proves density zero, without exchanging the limits in source
height and word length.

Now suppose `A<=q_j<=B` throughout an orbit. The no-repeat argument
from section 2 still applies, and (D1) gives the two-sided linear bounds

```text
(A/(3B))*j <= n_j <= n_0/A+(B/(3A))*j.                   (D2c)
```

For any ell>=0, divide the lower bound at j+ell by the upper bound at j.
For all sufficiently large j this proves

```text
n_(j+ell)/n_j >= (A/B)^2/2 = epsilon > 0,  for all ell>=0.
```

Every sufficiently late n_j belongs to E_epsilon. But the n_j are
distinct and their upper bound is linear in j, so the set of these
tail values has **positive lower natural density**: for height X take
all indices up to `(3A/B)*(X-n_0/A)`. This contradicts the density-zero
lemma. If the displayed density lower bound would exceed the available
integers, the counting contradiction is already immediate.

This is a legitimate use of spatial density on a particular orbit:
the hypothetical strip first forces that orbit to occupy a positive
density of sources, all of which have the exceptional property.
Without the linear bound and the uniform future/source ratio, this
inference is unavailable. The result supplies no bound on arbitrary
unbounded excursions and does not establish global descent.

## 3. A precise Sturmian connection, and its integer obstruction

The threshold for a contracting *leading multiplier* is exactly

```text
2^K_L > 3^L  iff  K_L >= floor(L*alpha)+1.                 (D4)
```

For any real intercept rho define the mechanical word

```text
K_j = floor(j*alpha+rho)-floor(rho),
k_j = K_j-K_(j-1) in {1,2}.
```

Its discrepancy lies in an interval of width one. To check balance, write
`beta=alpha-1` in `(0,1)`; the number of symbols 2 in a block of length m
is either `floor(m*beta)` or `ceil(m*beta)`, by telescoping the floors.
The slope is irrational since a rational `log_2(3)` would give `2^u=3^v`
for positive integers u,v. Thus this is the usual irrational mechanical,
or Sturmian, coding with binary slope `log_2(3/2)`. A golden-ratio slope
does not enter this derivation.

**PROVED.** No such infinite mechanical word is the exact halving word
of a positive integer Collatz trajectory. Its q-values have supremum to
infimum ratio at most two, contrary to (D2). Nor can a positive integer
trajectory have an eventually mechanical tail of this type.

Nevertheless every finite prefix is legal on a nonempty positive integer
progression. The exact source class for a word of length L is

```text
n = (2^K_L-B_L)*(3^L)^(-1) mod 2^(K_L+1).                 (D5)
```

The source cylinders for successive prefixes are nested and determine a
2-adic integer. The least positive representatives are nondecreasing and
unbounded: if bounded, they would eventually be constant, giving one
positive integer realizing the forbidden infinite word. This is an
explicit example in this workspace of **finite realizability at every
length without positive-integer realizability of the infinite limit**.
It is an obstruction theorem, not an infinite divergent Collatz example.

## 4. Signed and higher-multiplier extension

Let a>1 and b be odd integers, with a positive and `gcd(a,b)=1`. Consider
an infinite orbit of positive odd integers under

```text
n_(j+1)=(a*n_j+b)/2^k_(j+1),
k_(j+1)=v_2(a*n_j+b)>=1.
```

The domain hypothesis includes `a*n_j+b>0` at every step. Set
`q_j=2^K_j/a^j`. Exactly as above,

```text
n_j*q_j = n_0 + (b/a)*sum_(i<j)q_i.                      (D6)
```

**PROVED, b>0.** If `0<A<=q_j<=B` throughout the orbit, then

```text
B/A >= 2*a^2/(b*phi(a)).                                  (D7)
```

The no-repeat proof still holds since odd a>1 cannot have a positive power
equal to a power of two. Every image is a unit modulo 2a, since b is a
unit modulo a. Such integers have density `phi(a)/(2a)` with bounded
counting error for fixed a. Apply that density to
`n_j<=n_0/A+(b/a)*(B/A)*N` and divide by N. When the right side of (D7)
is at most one it supplies no improvement over the tautology `B/A>=1`.

The arithmetic sidecar can be sharpened. For b>0 let `d=gcd(n_0,b)`.
Since a is a unit modulo b and division by two is invertible modulo odd b,
`gcd(n_j,b)=d` is invariant. The possible images therefore have density
`phi(a)*phi(b/d)/(2*a*b)` rather than just `phi(a)/(2*a)`. The same proof gives

```text
B/A >= 2*a^2/(phi(a)*phi(b/d)).                             (D7a)
```

Here the integer-count error is bounded for each fixed a,b. On the
invariant sector d=b this recovers the b=1 bound, consistently with the
scaling conjugacy `T_(a,b)(b*m)=b*T_(a,1)(m)`. This refinement illustrates
exactly when restricting the available integers improves a capacity proof.

**PROVED, b<0.** Positivity alone forces

```text
sum_(j>=0)q_j <= a*n_0/abs(b),
q_j -> 0, and K_j-j*log_2(a) -> -infinity.                (D8)
```

Indeed every partial sum in (D6) is strictly below `a*n_0/abs(b)`.
The positive-term series therefore converges, and its terms tend to zero.
No bounded-strip assumption is needed. For `3n-1`, this gives a genuine
signed distinction from `3n+1`, while remaining compatible with all three
known positive `3n-1` cycles. It neither excludes further cycles nor proves
that every positive orbit enters one of them. The same calculation applies
to `3n-5` on any orbit remaining positive; at n=1 that map leaves this
domain, so no global-positive-domain claim is made.

There is also a precise escape test. The nonnegative limit
`ell=n_0-(abs(b)/a)*sum_(j>=0)q_j` equals `lim(n_j*q_j)`.
If ell>0, then n_j tends to infinity because q_j tends to zero.
Every bounded positive orbit thus has ell=0. The converse is not claimed.

The stronger no-bounded-strip theorem from section 2a also extends to
every fixed **positive** odd b in `3n+b`, even if b is divisible by three.
Finite halving-word cylinders still have density `2^(-K_L)` because
three is invertible modulo every power of two. The uniform intercept
bound becomes `b*((3/2)^L-1)`, and both linear growth constants acquire
the factor b; the limiting future/source ratio is still `A^2/B^2`.
This extension does not apply unchanged to `5n+b`: its source-average
halving count 2 is smaller than `log_2(5)`, so the stopping-time density
argument has the wrong drift. The direct bounds (D7)-(D7a) still apply.

## 5. What this changes in the research program

The [first-reset calculation](collatz_guards_20260921_valves.md) uses the
same threshold alpha and obtains an exact density whose binary positions
are `floor(L*alpha)`. This is a proved connection between two derived
objects, rather than an analogy between visually similar sequences.
Its source-counting density does not prove a sampling law on each orbit.

The [prime-square reset obstruction](collatz_guards_20260921_squarefree.md)
retains another lost coordinate: prime-square positions along repeated
affine branches. Our capacity argument constrains words nearly balanced
between powers of two and three; that result constrains squarefree runs
of a strongly contracting word. Neither supplies coverage of all words.

**OPEN research question:** can a counterexample orbit be divided into
finite excursions with a uniform, demonstrable restriction on *ordinary
integer availability* strong enough to force a descent? The present proof
rules out every fixed bounded discrepancy strip. Unbounded excursions
and the carry-dependent return inequalities remain uncontrolled. A
possible extension would replace all units modulo six
by a rigorously narrower set occupied by a specified excursion family;
its density and its compatibility with the height bound must both be
proved. Merely drawing a smaller residue graph would discard the needed
lift and would not suffice.

## Reproduction

Run `python 04-computation/experiments/collatz_guards_20260921_discrepancy.py`
and the same command with `python -O`. The [script](../../04-computation/experiments/collatz_guards_20260921_discrepancy.py)
uses exact integers and fractions. It checks all mechanical prefixes
through length 128, direct odd orbits for all odd starts 1..1999 through
100 steps or first repetition, the exact unit count through height 2999,
and specified `3n-1`, `5n+1`, `5n-1`, `3n+5`, and positive `3n-5` controls.
The `3n-5` source 3 is also checked as a hostile that leaves the positive
domain via 1; it is not silently included in the infinite-positive theorem.
The source-density calculation also checks the negative-binomial formula
`Pr(K_L=s)=binomial(s-1,L-1)/2^s` against an independent convolution,
and exhaustively checks the bad-event cylinders for L=4 and L=8.
The [JSON](../../04-computation/experiments/collatz_guards_20260921_discrepancy.json)
is a finite certificate of these controls, not a numerical proof of the
infinite statements above. The proofs of (D2), the stronger section 2a
theorem, and (D7)-(D8) were also
independently audited by the formalization lane; they have not been
translated into Lean in this checkpoint.
