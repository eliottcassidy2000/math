# Finite polynomial valuations cannot repair a globally decreasing logarithmic Collatz potential

**Status: PROVED elementary theorem; independently audited proof; FINITE-EXACT
controls. Collatz remains OPEN. No literature-priority claim.**

This strengthens the single-resource reset obstruction in
[the crossing note, section 5B](crossroads_crossing_20260926_potential.md).
It addresses arbitrary nonlinear dependence on finitely many unbounded
valuation features, not merely bounded colour or residue corrections.

## 1. Statement and scope

Let S be any finite set of primes and let P_1,...,P_m be nonzero polynomials
in Z[X]. Let f be any real-valued function of their finitely many valuation
coordinates, let g be any bounded function on the positive integers, and
let c>0. At integers where all P_i(n) are nonzero, define

    V(n)=c log n+f((v_p(P_i(n)))_(p in S, 1<=i<=m))+g(n).

**Theorem.** There is no cutoff H such that V(T(n))<=V(n) for every n>H
under shortcut Collatz, with the potential defined at the endpoints.
In fact, for each H there are arbitrarily large complete positive orbit
segments staying above H, whose endpoint valuation vectors are exactly
equal and whose endpoint height ratios are arbitrarily large.

The coordinates need not be bounded; f need not be linear, continuous,
computable, or bounded. The finite set of positive integer roots of the
P_i can be absorbed into H. The same conclusion holds if nonincrease is
required at every sufficiently large accelerated odd-to-odd transition:
each repeated word has odd-step valuation word 1^(a-1),(B+1), and its
endpoints are odd. The result includes valuations at {2,3,11},
at the mod30 primes, and of forms n+1,n-1,3n+5,2n+1 together. Bounded
Fibonacci, log-periodic, finite-state and residue corrections can enter g.

This does not exclude a ranking function with an unbounded history, an
unbounded number of arithmetic features, general digit transducers,
selected return times, or a fundamentally different leading height term.
In particular it is not a theorem against all possible Collatz potentials.

## 2. A rational periodic shadow avoiding every selected polynomial zero

For p>=5 in S, let ord_p(b) denote multiplicative order modulo p. Choose
B>=1 divisible by every ord_p(2), and L>=1 divisible by every ord_p(2)
and ord_p(3). Empty lists impose no condition. For arbitrarily large
a=1 mod L, we have 3^a>2^(a+B). Put

    D_a=3^a-2^(a+B),
    r_a=-(3^a-2^a)/D_a,
    M_a=3^a/2^(a+B)>1.

The denominator D_a is odd and prime to3. At every p>=5 in S it is
3-2=1 modulo p. Thus r_a is p-adically integral at every selected prime,
and also at2 even when2 was not selected.

These rational points are all distinct:

    r_a=-1-(2^B-1)/[(3/2)^a-2^B],

so they approach -1 strictly from below. The finitely many polynomial
roots exclude only finitely many a. Fix a single a for which every
P_i(r_a) is nonzero, and abbreviate r=r_a, M=M_a, ell=a+B.

The word w=1^a0^B is the actual parity word of a rational Collatz cycle
through r, with the usual parity of a rational having odd denominator.
This can be checked without an abstract conjugacy theorem:

    r+1=-2^a(2^B-1)/D_a,       v_2(r+1)=a,
    T^a(r)=2^B r,              v_2(r)=0.

The first a sources are odd, and the next B sources are even, returning
to r. The affine word map is

    T_w(x)=r+M(x-r).

The cycle lives at negative rational values. It is a shadow used to build
positive finite segments, not a positive integer counterexample.

## 3. Arithmetic precision preserves all endpoint features

For p in S union {2}, choose

    e_p > max_i v_p(P_i(r)).

If n-r is divisible by p^e_p, polynomial continuity gives
v_p(P_i(n))=v_p(P_i(r)) for every i. This follows directly by factoring
P_i(n)-P_i(r)=(n-r)Q_i(n,r), where Q_i(n,r) is p-adically integral.

For any integer k>=1 use the Chinese remainder theorem to choose an
arbitrarily large positive integer n satisfying

    n=r modulo 2^(k ell+e_2),
    n=r modulo p^e_p for every odd p in S.

Rational residues are well-defined because the denominator of r is a
unit at these primes. The first k ell actual parity bits of n are w^k:
each shortcut step loses at most one bit of agreement with the rational
cycle. At the endpoint,

    y=T^(k ell)(n)=r+M^k(n-r).

For p=2 the initial k ell extra bits pay for all divisions. At p=3 the
agreement gains ak digits; at the other odd selected primes it loses
none. Therefore n and y have exactly the same full valuation vector.

Choose n>2^(k ell)H, and larger if necessary. Every shortcut step is at
least half its positive source, so every intermediate value remains
above H. As r<0,

    y/n=M^k+(1-M^k)r/n > M^k.

The endpoint value of f cancels exactly, while |g(y)-g(n)|<=2 sup|g|.
Thus

    V(y)-V(n)>ck log M-2 sup|g|,

which is positive for large k. Summing the hypothetical one-step
nonincrease inequalities along this actual high positive block gives
the opposite inequality. Contradiction.

## 4. What changed in the research direction

The inherited all-odd-prefix obstruction excluded bounded corrections.
The 41->31 reset then excluded arbitrary functions of v_2(n+1).
Rational periodic shadows now exclude arbitrary dependence on any fixed
finite collection of polynomial valuations, even if those valuations can
be arbitrarily large across the integers. The mechanism is local constancy
away from polynomial zeros, combined with expanding positive shadows.

The argument's quantifiers matter. One fixed negative rational cycle and
one fixed valuation vector produce positive witnesses for every block
length, but the positive starting integer changes with that length. This
is enough to refute a universal one-step potential. It says nothing about
whether one fixed positive integer follows the shadow forever.

Source: a negative rational periodic orbit. Target: arbitrarily long
positive integer prefixes, using simultaneous p-adic approximation.
Preserved: the word, endpoint polynomial valuations, and actual affine
carry. Lost: a common positive starting integer across lengths. Necessary
sidecar for a convergence argument: a coordinate sensitive to that changing
precision or an independently well-founded stopping construction.

The concrete next question is therefore narrower: can a source-preserving
crossing potential use variable-depth information with an amortized cost
that is not locally constant near every expanding rational cycle? A finite
list of additional valuation features cannot answer this by itself.

## 5. Inheritance, controls, and audit

Closest proved mechanism: [the blueprint energy audit, section 3](collatz_blueprint_20260921_energy.md).
Canonical hostile: the all-odd 2-adic point -1. Corrected near miss:
finitely many unbounded valuation features are stronger than a bounded
phase, so the old argument alone did not exclude them. Least-used sidecar:
all places at which the rational cycle denominator must be invertible.
The live board is rational cycles, positive shadows, polynomial roots,
valuation precision, crossing order, and fixed-integer quantifiers.

Run `python -B 04-computation/experiments/crossroads_crossing_20260926_resource.py`.
The [output](crossroads_crossing_20260926_resource.out) tests five prime sets
up to {2,3,5,7,11,13}, eight polynomial forms, and 1,2,5,10,20 repeats of
each selected word. It verifies the rational cycle directly, constructs
actual positive CRT witnesses, checks every parity and intermediate core
bound, checks the unwrapped affine endpoint, and compares every endpoint
valuation exactly. All checks remain enabled under Python -O.

The form X+5 is a hostile control: it invalidates the naive shadow -5,
so the constructor must select a different rational cycle. The proof was
independently audited by the colour lane, including the direct rational
word certificate, polynomial-root exclusion and all p-adic precision losses.
The crossing lane independently probes the {2,3,11} instance. No density,
random-word law, or assumption of Collatz divergence is used.
