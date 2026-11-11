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

## A two-sided physical reader extracts a step of2

There is a sharper cross-lane consequence than the primitive difference gcd.
Let L be the two-vertex transitive tournament whose root is its source, and
let R be a rooted triangle. Form X=L⋆A⋆R and Y=L⋆B⋆R, using the SAME labelled
contexts on both sides. The Fano data have

    delta(H,h,Q,k,U,V,d,e)=(0,0,0,0,3,-3,-1,1).

The fusion matrix law therefore gives, for arbitrary pointed L,R,

    H(L⋆A⋆R)-H(L⋆B⋆R)
      = h_L(3V_R-e_R)-h_R(3U_L-d_L).                    (B3)

For this L, (h,U,d)=(1,2,0); for R, (h,V,e)=(1,3,1).
Thus the difference is (9-1)-(6-0)=2, with the actual counts

    H(X)=1771,    H(Y)=1769,    |X|=|Y|=10.

Both outputs are strong. The root-sharing A/R portion is strong. The extra
nonroot vertex of L receives an arc from the root and sends arcs to all
nonroot A vertices, so it lies in that same strong component. The argument
applies equally to B. The executable provides both literal adjacency lists,
subset-DP counts and independent Hamiltonian-path enumeration.

This is not a formal subtraction construction: X and Y are two explicitly
built positive tournaments. The subtraction compares their counts. It shows
that the converse-hidden response can produce the smallest possible nonzero
scalar difference after a fixed ordinary context. It does NOT supply such a
step at an arbitrary starting count, a reusable +2 operation, or an interval
of attainable counts. The values1769 and1771 themselves were already within
the inherited covered range; the new information is this context and its
information-loss explanation.

One-sided fusion cannot make the same extraction. With only a right context
R, write a_e,a_o for its nontrivial root-ending even/odd weighted path masses.
The difference is4a_e+2a_o. If a_o>0, a nontrivial odd root-ending path has an
incoming final edge, and the two-vertex root-ending path contributes at least1
to a_e. Thus a nonzero one-sided response is at least4. The two-sided reader
is genuinely using a cancellation between two physical response channels,
not a hidden one-sided selector.

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
