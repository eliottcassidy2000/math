# Every fixed residue observer loses the twin-prime predicate

**PROVED elementary CRT obstruction; FINITE-EXACT implementation controls.**
This does not rule out unbounded sieves, algorithms with growing moduli, or
observers that retain the integer itself. It is a precise failure of a fixed
modulus to recognize the prime conditions in the new tree correspondence.

## Statement and construction

Let M>=1 and let p,p+2 be odd twin primes. There are infinitely many odd
composite x such that x+2 is composite and x=p modulo M. Consequently the
primitive Pythagorean triples

    T(x)=(x(x+2),2(x+1),(x+1)^2+1)

and T(p) are coordinatewise equal modulo M, but only T(p) has both its
inradius and its even-leg exradius prime. Both triples lie on the same
parabolic Berggren ray, so a branch-determinant parity bit cannot repair
the loss along this ray.

**Proof.** Choose distinct odd primes q,r not dividing2M. The simultaneous
congruences

    x=p mod2M,       x=0 modq,       x=-2 modr

have a solution by the Chinese remainder theorem. Every sufficiently large
positive integer in that residue class exceeds q and has x+2>r, making both
numbers composite. The class is an infinite arithmetic progression, and
x is odd. Since T is an integral polynomial map, equality modulo M passes
to all three coordinates. The formula is the Euclid lift of `(x+1,1)`;
its parameters are coprime and have opposite parity, so it is primitive
for every positive odd x, independently of primality. Its two distinguished
radii are x,x+2 by direct substitution. QED.

The theorem also works with x=p modulo2M, which preserves the ray depth
`(x-1)/2` modulo M. Thus finite residue data on depth cannot separate the
two cases either. For a finite collection of moduli, take their least
common multiple. A finite-valued observer that does not factor through a
residue modulus is outside this statement; the truth value of primality
itself is an obvious reason not to assert a broader impossibility.

## The explicit F7 / F8 hostile

At p3,M7, choose q3,r5. One such x is213:

    213=3*71,        215=5*43,        213=3 mod14,
    T(3)  =(15,8,17),
    T(213)=(45795,428,45797),
    both triples reduce to (1,1,3) mod7.

In the residual slope chart r=n/m, both have `r=1/4=2 mod7`. The finite-field
note uses the reciprocal coordinate t=m/n=4 for its exponential marking.
Any fixed marking of
that conic point into F8, and the resulting Fano/Steiner block-incidence
data, therefore agrees for these two triples. Yet their radius pairs are
respectively twin primes and two composites. Reduction to F8 here means
the explicitly marked finite-set correspondence, not a field embedding
F7 into F8. No such field embedding exists because the characteristics differ.

This explains how the finite-field map can preserve a real geometric
predicate while losing the prime condition defining the selected tree
vertices. It provides no reduction of twin-prime infinitude to a finite
design question or to the ordinary tournament H gaps7,21.

## Inheritance, proof boundary and replay

The source mechanism is the new exact twin/Euclid chart, together with the
classical elementary CRT. The inherited warning is the [previous packet's
marked F8 rank and map distinction](../../aliquot_power_124_20261026/CONNECTION_BOARD.md).
The least-used sidecar is the complete integer lift, including factorization;
the finite residue only retains its class. This uses the repository's
"controlled forgetting and unlabeled quotients require a sidecar" and
"type every analogy and every implication" cards.

Run `python3 -B check.py` or `python3 -B -O check.py` from this directory.
The checker constructs the CRT classes for all twin lower primes through500
and every modulus1..120. It checks literal factors, depth congruences,
primitive triples, and the explicit seven-field witness. These are controls
for the all-parameter proof, not an extrapolation to infinitely many primes.
