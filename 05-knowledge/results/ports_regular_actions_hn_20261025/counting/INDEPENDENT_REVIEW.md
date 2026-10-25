# Independent review of the fixed-degree response quotient

2026-10-25. Reviewer: `nonconsecutive_ports`, independently of the author of
[README.md](README.md). **VERIFIED:** no blocking mathematical or implementation
defect found in the stated ordinary fixed-target theorem. **FINITE-EXACT:** the
independent checks below; they do not replace its all-order proof.

## Scope and algebra

The claimed quotient keeps the stage, used degree, and every affordable
coordinate `D(d)=d*m-f(d)`. It removes the past cycle count only after multiplying
the accumulated weight by `2^(R*m)`. This is the correct direction of scaling:
for a continuation of actual degree R, its interaction exponent is

    sum_d x_d*f(d) = R*m - sum_d x_d*D(d).

The equality applies to each individual continuation, with its original
denominator unchanged. It therefore justifies summing scaled weights in the
same key, including histories with different numbers of cycles.

For k new ell-cycles and a smaller continuation nu of degree R'=R-ell*k,
the internal exponent is `(ell*k*k-k)/2`. Updating
`D'(d)=D(d)+k*(d-gcd(ell,d))` gives

    -sum_(d in nu) D'(d)
      = -sum_(d in nu) D(d) - k*R' + k*sum_(d in nu)gcd(ell,d).

Consequently the multiplier that belongs outside the smaller F is

    (ell*k*k-k)/2 - k*D(ell) + k*R'
      = k*(R-D(ell)) - (ell*k*k+k)/2.

This independently reproduces equation (5), including the sign of D, the
quadratic term, and the factor `ell^k*k!`. Negative exponents are legitimate
exact rational factors. No hidden integrality claim is required at a transition.

The boundary R=0 has one empty continuation. At ell=1, D(1)=0 and there is one
all-ones partition, giving `2^binom(R,2)/R!`. Skipping ell>R does not require
reading an absent coordinate. For positive N and a nonempty affordable
coordinate set, past lengths at least ell+2 imply

    0 <= D(d) <= m*(d-1) <= u*(d-1)/(ell+2) < N.

The producer implements these same rules and uses explicit exceptions, so its
validity checks survive Python optimization.

## Coarseness and what it does not prove

At fixed R and a decreasing stage, the all-ones continuation has normalized
interaction 1. Replacing d ones by a single d-cycle changes it to `2^(-D(d))`.
These are valid tests for every affordable odd d. Thus agreement on all
individual normalized continuations forces equality of all the retained D
coordinates, and that equality is also sufficient.

Equivalently the vectors `e_d-d*e_1`, for affordable d>1, span the kernel of
the actual-degree functional. Its annihilator is exactly the degree direction
that the scalar gauge removes. This argument would fail if the all-ones probe
were forbidden or had already been processed. It is not minimality among
representations that preserve only the final sum, and does not give a bounded
state recurrence, a runtime improvement, or a Burnside-free classification.
The README makes all these restrictions explicit.

The strict mergers `(23)` versus `(5,9,9)` and `(27)` versus `(9,9,9)` have equal
used degree and unequal old responses, but the stated equal D values. Their
positive weights prevent cancellation from being the explanation for the
ordinary state-graph reduction.

## The excess-degree hostile is an actual term

The clone-removal substitution in
[the inherited boundary packet, section 3](../../defect_boundary_20261020/boundary.md)
is `p -> p/(1+p)`. In particular

    [p^4] p/(1+p) = -1.

Thus a single actual future fixed point and three excess units is a nonzero
available term, not a hypothetical assignment outside that substitution.
For prefixes `(27)` and `(9,9,9)`, with R=4, normalizing its actual interaction
by `2^(4*m)` gives `2^(-3*m)`: respectively 1/8 and 1/512, although D agrees.
The sign and substitution coefficient are common and cannot repair this
missing distinction. Actual degree or a suitable excess grading must be kept
to extend the quotient. Literal ports, clone identities, and canonical owners
are already absent before this numerical quotient; no conclusion about their
transport follows from its exact counting property.

## Independent coefficient route and finite checks

For an odd partition lambda, let m be its number of parts and let b_q count
its parts divisible by q. The divisor identity
`gcd(a,b)=sum_(q|a,q|b) phi(q)` gives

    e(lambda) = (sum_(odd q) phi(q)*b_q^2 - m)/2.

Indeed the diagonal terms of the double gcd sum are the part lengths; removing
m and dividing by two supplies exactly the within-cycle contributions and the
unordered cross-cycle contributions. Summing `2^e/z_lambda` by multiplicity
vectors is an alternative coefficient computation, without the new response
state or its transition formula. The code below agrees with the producer for
every N=0..32.

A second control uses the literal action of every permutation on **ordered**
vertex pairs, without a gcd formula. An orbit equal to its reversed orbit
forbids an invariant tournament. Otherwise reversal pairs the ordered-pair
orbits, with one independent orientation choice per pair. Averaging these
counts agrees for N=0..7, yielding `1,1,1,2,4,12,56,456`.

The producer was also rerun in normal and `-O` modes. Both passed the 2,904
normalized-future comparisons, independent partitions through 27, literal
edge actions through 7, invalid inputs, both strict mergers, excess hostile,
and the 30030-offset boundary. Its old/new graph counts at 30,40,60,80,100
were reproduced exactly. The independent script below also passes both modes.

The 30030 example retains equal degree 60060 and equal cycle count. Its first
different odd response is 59 because `30031=59*509`, whereas 30013,30047,30029
are prime. It is correctly labeled an old-style finite-horizon merger rather
than evidence specific to the new gauge.

To reproduce the independent check, run this code from the repository root.
It imports only the candidate answer function; the two controls are separate.

```python
from fractions import Fraction
from itertools import permutations
from math import factorial, gcd
import importlib.util

path = '05-knowledge/results/ports_regular_actions_hn_20261025/counting/check.py'
spec = importlib.util.spec_from_file_location('candidate', path)
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)

def phi(q):
    return sum(gcd(i, q) == 1 for i in range(1, q + 1))

def multiplicities(n, d=1):
    if d > n:
        if n == 0:
            yield {}
        return
    for k in range(n // d + 1):
        for tail in multiplicities(n - k*d, d + 2):
            if k:
                tail = dict(tail)
                tail[d] = k
            yield tail

def divisor_count(n):
    total = Fraction(0)
    for mult in multiplicities(n):
        m = sum(mult.values())
        top = max(mult, default=0)
        square_sum = sum(
            phi(q) * sum(k for d, k in mult.items() if d % q == 0)**2
            for q in range(1, top + 1, 2))
        if (square_sum - m) % 2:
            raise RuntimeError('nonintegral divisor exponent')
        denominator = 1
        for d, k in mult.items():
            denominator *= d**k * factorial(k)
        total += Fraction(2**((square_sum-m)//2), denominator)
    return total

def ordered_pair_count(n):
    total = 0
    for p in permutations(range(n)):
        remaining = {(i, j) for i in range(n) for j in range(n) if i != j}
        orbits = []
        while remaining:
            cur = next(iter(remaining))
            orbit = set()
            while cur not in orbit:
                orbit.add(cur)
                cur = (p[cur[0]], p[cur[1]])
            remaining -= orbit
            orbits.append(frozenset(orbit))
        seen, pairs, valid = set(), 0, True
        for orbit in orbits:
            reverse = frozenset((v, u) for u, v in orbit)
            if orbit == reverse:
                valid = False
                break
            if orbit not in seen:
                seen.add(orbit)
                seen.add(reverse)
                pairs += 1
        if valid:
            total += 2**pairs
    if total % factorial(n):
        raise RuntimeError('literal divisibility')
    return total // factorial(n)

small = []
for n in range(33):
    answer = candidate.count_target(n)[0]
    if answer != divisor_count(n):
        raise RuntimeError(('divisor mismatch', n))
    if n <= 7:
        literal = ordered_pair_count(n)
        if answer != literal:
            raise RuntimeError(('ordered-pair mismatch', n))
        small.append(literal)
print('PASS: divisor coefficients 0..32; ordered-pair actions 0..7;', small)
```

## Addendum: independent audit of the excess extension

2026-10-25. **VERIFIED after the two explicit wording repairs below:**
[EXCESS_EXTENSION.md](EXCESS_EXTENSION.md), including X1--X9. The exact
bivariate recurrence and both individual-future equivalence claims hold.
This promotes no state-count, runtime, bit-complexity, or structural
clone-allocation claim.

For a fixed coefficient b_(M,E), only actual future cycle degree R belongs
in the scalar exponent `R*m`. Excess Q contributes the independent factors
`beta(k,e)` and never changes a gcd interaction. Because the choice of
excess at a length is completed before that length is discarded, no
unresolved earlier allocation is omitted by the state (R,Q,D). The
ordinary transition derivation above therefore gives X4 unchanged, times
beta. In particular ell>R cannot support excess by itself. The R=0 and
ell=1 boundaries in X5 are correct.

**X7, exact live mask.** If d<R, one actual d-cycle and R-d singletons
support every Q by assigning all excess to the singleton block. Its beta
coefficient is nonzero. If d=R, one d-cycle exhausts actual degree, so it
supports Q exactly when R divides Q. The all-singleton continuation is
always available for R>=1. Dividing each continuation by its known
nonzero intrinsic coefficient reads `2^(-D(d))` for every coordinate in
X7. Coordinates outside the mask cannot occur. Both directions of the
stated semantic equivalence follow. The wording now distinguishes this
readout property from independent coordinate freedom or storage necessity
on a restricted reachable-state set.

The checker's support test aggregates excess assignments within an actual
partition. This does not accidentally erase support by cancellation:
all lengths are odd, so `sum e_d` has the parity of `sum d*e_d=Q`. Every
nonzero beta product has the same sign (-1)^Q. Its aggregate is nonzero
exactly when an individual assignment exists. X8 also checks directly:
beta(3,2)*8/6=8; beta(3,3)*8/6=-40/3; and the sole actual 3-cycle at Q=3
adds `-(2/3)*2^(-D(3))`.

**X9, projective equivalence.** At total ungraded remainder S>=2 the
one-actual-singleton and two-actual-singleton probes both have nonzero
coefficients. Their quotient, after removing known intrinsic factors,
determines `2^f(1)`. For 1<d<S the d-cycle plus one singleton probe,
relative to the one-singleton probe, then determines `2^f(d)`. At d=S the
one-cycle no-excess probe instead determines `2^(f(S)-f(1))`; f(1) is
already known. Thus proportional individual-future vectors force exact
equality of the old affordable f profile. Equality is sufficient term by
term. The S=0,1 and no-singleton exclusions are correctly retained.

**Correction and stronger meaning of the bivariate coefficients.** The
draft's assertion that individual b_(M,E) need not be integral was false.
They have the signed combinatorial interpretation

    (-1)^E b_(M,E)
      = number of isomorphism classes of M-vertex tournaments with
        nonnegative integer vertex colors summing to E.

To prove this, first bound the allowed colors by E, so ordinary finite
weighted Burnside applies coefficientwise. A coloring fixed by a
permutation is constant on each cycle. Only permutations with odd cycles
can fix a tournament. On an odd cycle of length d, summing its color j
with weight (-z)^(d*j) gives 1/(1+z^d), up to the harmless truncation.
Multiplying these factors and averaging fixed tournaments gives precisely
`Z_T(s^d/(1+z^d))`. Every coloring in the coefficient has sign (-1)^E.
This proves integrality and the fixed sign at every M,E; it also explains
the integer coefficients despite rational intermediate cycle-index terms.
The section 1 text has been repaired explicitly, preserving the draft
correction lineage. M and E remain auxiliary to the core substitution:
they describe a colored source tournament, not two vertex counts of a
particular resulting core. This adds no literal owner or port recovery.

**FINITE-EXACT controls rerun.** Both normal and `-O` executions of
`check_excess.py` pass: 120 bivariate coefficients with M+E<=14; literal
pair-free classes at orders 0..5; 288 formal-profile comparisons; the
scalar hostile; visibility and invalid-input controls. Separately, direct
enumeration and simultaneous permutation canonicalization of tournament
arcs and vertex colors agrees with the signed interpretation on all 30
pairs M=0..4, E=0..5. Integrality and the prescribed sign were additionally
checked for M+E<=20. These checks corroborate the arguments above; the
unbounded statements are proved by the displayed continuations and
weighted Burnside, not by extrapolation from this finite universe.

The independent colored test is preserved as `check_colored_inventory.py`;
run it normally and with `-O`. The main `check_excess.py` now also rejects
nonintegrality or the wrong sign for every coefficient in its own universe.
