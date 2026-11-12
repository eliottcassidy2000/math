# Four exact corrections, one principle: retain interactions before estimating error

**Session278, 12 November 2026. Chan–Pak and planar JC remain OPEN.**

Session277 showed how to recover an integer from an absolute interval and an
exact residue. This session improves what enters that recovery: rather than
correcting a scalar approximation blindly, retain the interaction that produced
its error, resum an exactly manageable part, and certify the remaining sign.
The gains below concern distinct mathematical quantities and are not conflated.

## 1. A000568: turn arithmetic symmetry into an interaction gas

In the odd-cycle Burnside formula, replace each gcd(a,b) interaction by
min(a,b). Unlike gcd, the latter admits a size-only dynamic program. This is
a certified upper envelope, not an exact count. Its missing information is
exactly a graph of incompatible cycle lengths, with edge weight

    q(a,b)=1-2^(gcd(a,b)-min(a,b)).

Each exact symmetry contribution is its envelope times product(1-q).
Elementary-symmetric corrections alternate with rigorously signed, two-sided
remainders. Resumming **every interaction involving a 3-cycle** before expanding
the rest substantially tightens the absolute certificate.

At order64 the answer has1721bits. The resulting interval requires only738 sufficient low
residue bits to identify it uniquely, versus1544 for Session277's support-cut3
interval. This is a stronger certificate, **not a speedup claim**: the current
correction-coefficient implementation still enumerates partitions, and these
are different workloads. Rational controls and certificates run through64;
independent permutation Burnside controls run through7.

There is a sharp structural boundary: four signed terms in the resummed
expansion are exact through order23 and first fail at24. Explicit symmetry
families prove that **every fixed correction depth eventually loses integer
accuracy**. Recovery by this expansion alone therefore needs adaptive depth or further exact
resummation; an independently computed residue is another route.

The unlikely source connection is an old LRC Brownian min-kernel. Both min and
gcd have Gram representations, but their difference is only positive on
nonnegative occupancies, not positive semidefinite. Keeping that distinction
produces the useful defect graph without a false spectral analogy.

[Proof, algorithm, exact error tables](enumeration/README.md).

## 2. Chan–Pak constructions: lost root memory is a signed matching polynomial

For the inherited moving-root chain of one-apex tournament blocks, discarding
one root coordinate gives a convenient two-state product and the wrong count.
Eliminating that coordinate exactly yields

    x_i=B_i x_(i-1)-D_i x_(i-2).

The correction is a signed sum over disjoint eligible boundaries of the block
path. **Its exact degree is the maximum matching number of those boundaries.**
A positive causal operator proves every truncation has an alternating signed
remainder, with an explicit bound. This is stronger than merely observing
alternating coefficients.

Repeated directed triangles give a particularly revealing infinite family.
At30vertices the blind product gives59,932,665 rather than53,277,425. The
relative error grows exponentially with chain length. Exact generating
functions for individual correction layers prove that **no fixed layer depth
repairs even relative accuracy**: full resummation changes the dominant growth
pole. Very small local memory losses can therefore alter the infinite system's
growth law.

This is physically constructive: each block word produces an actual tournament,
and its signed bounds can reject an incorrect proposed target before building
it. It still does not prove that every prescribed target has such a word.
Tests include1554 algebraic chains,258 independently counted physical chains,
and correction-layer generating functions through length30.

[All-length proof, sharp termination, examples and code](physical/PROOF.md).

## 3. Simultaneous arc edits: cancel shared error before taking absolute values

For a tournament and a chosen set of reversals, form the union digraph and
retain jointly how many original and reversed edited arcs each Hamilton path
uses. This paired ledger yields a finite exact signed correction at every
order, with forced-path contractions computing its moments.

A four-vertex example changes5paths to1. Its first correction predicts-8
instead of the true-4. The second correction is exact even though both separate
count estimates still have error1: their tails cancel. Independent error bars
miss that cancellation; the joint response retains it. The underlying forest
insertion and single-edge laws are inherited/classical, not claimed new.

[Paired certificate and255 subset-DP controls](edits/PROOF.md).

## 4. Gaussian/JC ideas: correct the polynomial representative first

For a polynomial potential W and a chosen vector field xi, retain the exact
Euler defect R=xi(W)-W. We derive a finite noncommutative recurrence reducing
W^m h modulo polynomial divergences, with explicit vector-field certificates.
It never treats the Stein quotient as an algebra. A simultaneous
homogeneous-defect generating formula retains mixed corrections; the quartic
case specializes to finite Laguerre polynomials.

Then apply a second correction axis: rigorous signed Gaussian Taylor bounds.
For the specified quartic tilt g=1/1000 and radial power m=6, rewriting first
improves the absolute interval width by a factor5062 at the same expansion
order. Every bound is rational; no numerical quadrature supplies the proof.
An explicit divergent-series boundary prevents a false arbitrary-precision
promise. This is a tool for controlled representatives, not a JC proof.

[Defect compiler, mixed formula and exact interval experiment](cross/PROOF.md).

## 5. The unifying rule and the remaining proof obligations

The common pattern is **interaction-first correction**:

1. Name what the approximation forgot: divisibility pairs, root memory,
   shared edited paths, or defect multiplication.
2. Preserve that information in a finite object with a precise operation law.
3. Resum a manageable part exactly, then certify the residual sign and size.
4. Intersect certificates from different representations only when they refer
   to the same quantity. Use an integer residue only for an integer target.

These mechanisms extend signed correction across four systems. They do not
imply that their correction objects are interchangeable: a negative Gaussian
coefficient is not a negative number of tournament witnesses.

The sharp remaining Chan–Pak question is still a **realizable absolute-target
selector**, not merely an accurate census, a gap, or a formal expansion. For
A000568, the next algorithmic opportunity is computing low defect coefficients
without full partition enumeration. For planar JC, it is a uniformly controlled
representative of the relevant obstruction with the needed source algebra.

An independent new-claim review covers all four lanes, with small notation,
provenance and finite-moment repairs, not a re-audit of the inherited corpus.
[Review and exact-source receipts](audit/REVIEW.md).
Concurrent S94 work at2b9e27578588 is preserved; its affine-cube/anchoring
obstructions are compatible with these deliberately scoped local certificates.
