# The same precision problem, two different observers

Status: **PROVED CONSEQUENCE of inherited fixed-fiber classification and this
session's granularity bound**, not a new Chan–Pak coverage theorem.

## Why a fully accurate leading census still cannot select a target

Write A_n(h) for the number of ordinary n-vertex tournament isomorphism classes
with exactly h directed Hamiltonian paths. Fix a positive integer h. The
inherited strong-core decomposition (Session266, `fibers/PROOF.md`, crediting
Session24 and opus S50) yields a finite sum

    A_n(h) = Σ_W binom(n-s(W)+k(W), k(W)),

where W is an ordered word of nontrivial strong components, their path-count
product is h, and s and k are its total size and word length. The summand is
zero before n>=s. A nontrivial strong tournament has at least three directed
Hamiltonian paths, so k<=floor(log_3 h). Hence for fixed h

    A_n(h) = O_h(n^floor(log_3 h)).                       (1)

This statement includes the identically zero polynomial if no realization
exists. It does **not** assume Chan–Pak. For h=1 the only class is transitive.
For a finite fixed set F of target values the number of classes whose counts
lie in F is also bounded by a polynomial in n.

Compare (1) with the new positive-sector lower bound

    t(λ)>=2^((n-(1+log_2 n)^2)/2)

for every odd-cycle partition λ. The ratio of that lower bound to the right
side of (1) tends to infinity. Thus even ONE omitted Burnside sector eventually
has much larger absolute mass than ALL classes realizing ANY chosen finite
collection of prescribed H targets.

This gives a precise limitation of a particular information channel. For a proper positive-sector partial sum L and any certified upper bound
U>=A000568(n), the granularity bound implies A000568(n)-L exceeds that
fiber mass eventually. Both the true census and the census with those rare
fibers removed therefore lie in [L,U]; this interval alone cannot distinguish
them. This claim is about this specified lower-sum interval, not arbitrary
intervals of the same width.
Those removed collections are abstract comparison populations, not alternative
universes of tournaments. We do not claim an impossibility of all proofs using
enumeration; we identify what a scalar approximation fails to retain.

For uniform ordinary classes the probability of a fixed attained target is
polynomially bounded divided by A000568(n); under uniform labeled tournaments,
its numerator is at most n! times that polynomial and its denominator is
2^binom(n,2). Either way it vanishes at the quadratic-exponential scale. A
relative-distribution approximation can be excellent while being completely
uninformative about a fixed target such as 23, 7 or 21.

## Two exact completions, with different missing information

For A000568 the missing state is arithmetic: the interval and an exact residue
with modulus larger than its width isolate a unique integer. The engine gives
that residue through the actual symmetry-sector sum. It does not invent an
arbitrary residue compatible with the estimate.

For Chan–Pak the missing state is constructive: a path-count target must be an
actual value of a realizable endpoint/cover interface. A residue condition, a
range estimate, a gap between two values or a large counting exponent does not
supply that interface. The families lane retains a multivariate coefficient
selector before taking a univariate aggregate. Its regular-quotient collapse
is a controlled example of exactly the target information that averaging loses.

The common principle is **observer-relative completion**: attach the data the
consumer actually needs, and prove the reconstruction or lifting map. This is
not a numerical analogy between unrelated constants; it specifies both source
objects, both forgetful maps and the distinct data required to reverse them.

## Incoming improvement, not our claim

While this session ran, concurrent S93 published commit
`f90b22ab5daa626c21d92da4fac37d9d5052cbd8`, reporting the reviewed bound
A(Y)>=c Y^(3/4)/(log Y)^4 for the count of attained Hamiltonian-path values at
most Y (THM-5141). Here A(Y) is a spectrum-counting function, NOT A000568(n).
We inherit its stated scope and do not re-audit or rename it. It strengthens
global abundance, while the prescribed-target and eventual-threshold steps
remain open. A positive count of many targets and a complete specification of
all targets are different statements even when the exponent improves.
