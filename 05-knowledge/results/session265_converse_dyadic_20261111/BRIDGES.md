# A physical dyadic amplifier of the Fano-hidden response

**PROVED composite of the new Fano certificate and the inherited retained-root
endpoint consumer.** No scalar length projector or Chan–Pak closure claimed.

Let A,B be the labelled converse twins in `doubling/PROOF.md`. Their weighted
root-start path polynomials differ by the positive two-length packet

    P_A(w)-P_B(w)=2w²+w⁷.

For any rooted T, append an m-vertex transitive tail D_m by an ordinal join,
leaving the original root marked. Then

    P_(T=>D_m)(w)=P_T(w)(1+w)^m.                         (B1)

Proof: a rooted directed path first follows a rooted path inside T, then
visits an arbitrary subset of the tail in its unique increasing order. The
complement is the old complement followed by the remaining tail. Choosing
j tail vertices supplies binomial(m,j), with no extra completion multiplier.
Summing proves (B1), including m=0.

Now wrap each T=>D_m with two new vertices u,v, oriented u->(T=>D_m)->v->u.
The Session264 endpoint consumer identifies P_(T=>D_m)'s length-l coefficient
with the number of Hamiltonian paths starting at the old root and visiting v
at position l+1. Therefore the two actual (9+m)-vertex tournaments have

    identical complete unmarked path-cover profiles,
    identical complete root-deleted path-cover profiles,
    start counts 260·2^m and257·2^m,
    gate-position difference w(2w²+w⁷)(1+w)^m.          (B2)

The profile equality follows from converse invariance before the tail, ordinal
profile multiplication by r^m, and the inherited wrap operator
Q_W(r)=r²Q(r)+rQ(r+1)-rQ(r-1). It is not deduced from equal H alone.

All coefficients in the difference are nonnegative. Hence EVERY nonnegative
length-dependent observation of this gate sees a weakly larger response in
the first construction, and the total difference is exactly3·2^m. This is a
physical dyadic amplifier, not formal multiplication of a scalar path count.
For m>=1 the odd-minus-even root-leg difference is zero, since (1+w)^m
vanishes at w=-1. Thus one legal tail erases the parity difference while
doubling the total difference, without changing the input scalar H=171.

This gives a concrete consumer for the Fano character sector lost by the
converse fold. It also distinguishes two uses of powers of two in this
session: positive binomial channel amplification here, and residue-orbit
control in `collatz/PROOF.md`. Neither operation selects a single residue
channel as an unconditioned scalar tournament count.

`bridge.py` checks m=0,1,2,3 using literal subset Hamiltonian counts and a
gate-position dynamic program on all eight resulting tournaments. The
all-m statement is the bijection and inherited wrapper law above.

## The exact arithmetic distinction at63

Accept the inherited all-order exclusions of7 and21 from
`THM-1370-h-spectrum-omits-7-21-all-n.md`. Then every tournament with H=63 has
exactly ONE nontrivial SCC, whose count is63; all other SCCs are singletons.
Indeed SCC counts multiply, and every proper factorization of63 into odd
factors greater than1 includes7 or21. Neither can be a component count. Thus
63 is multiplicatively irreducible in the realizable-count monoid, although
it is composite as an integer. Its realization is inherited from
`THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`,
not newly discovered here.

The new shadow theorem locates that unique63 core and its singleton tails
using H-valued cuts, without arc directions. It does not obtain63 by multiplying
the unrealizable7 by two triangles. This is a useful meeting point of the
arithmetic and the literal weighted graph, and a concrete reason not to treat
strong components as ordinary integer primes.
