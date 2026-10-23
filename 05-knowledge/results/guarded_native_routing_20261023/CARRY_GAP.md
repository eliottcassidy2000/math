# Paying the missing kernel direction in an actual guard-preserving family

23 October 2026. **PROVED scoped**, extending the native order-20 example
in incoming [Session44 lift holonomy](../session44_guarded_native_moves_20261023/lift-holonomy/README.md),
sections 3 and 4. The inherited example authenticates the two generators
and shows why the ternary quotient loses connected components and slow
modes. Here the construction is extended to cell sizes M=3^k and the
repaired chain's full gap and return costs are computed. It stays within
one authenticated source-automorphism orbit. It does not connect arbitrary
path orbits or declare a new elementary metagraph edge.

## 1. The complete marked carrier and its retained kernel

Let M=3^k, k>=1, and K=M/3. Take a source singleton s, two regular cyclic
M-cells A,B, and a sink singleton t. Internally, i points to j iff
j-i mod M belongs to {1,...,(M-1)/2}. Across the cells put A_i -> B_j iff
i=j mod3, and reverse the other cross arcs. The singleton s beats all
other vertices and all vertices beat t.

In fact every argument below works for ANY odd positive multiple M of
three. The M=15 control explicitly checks this extension beyond powers of
three; the power family is a convenient recursively growing subfamily.

For a,b modulo M with a=b modulo3, the actual path is

    P_(a,b)=(s,A_a,A_(a+1),...,A_(a-1),
               B_(b+2),B_(b+3),...,B_(b+1),t).

The cross junction is legal because a-1=b+2 modulo3. Each path visits
every owner exactly once. All M^2/3 paths have the same literal endpoints
s,t. They are exactly this specified cyclic-cut language; no claim is
made that these are all Hamiltonian paths of the source.

Every pair (a,b) in this language is the image of P_(0,0) under the actual
source automorphism rotating A by a and B by b. Thus the literal complement
partners are also carried to one another. Their actual target class is
constant. EVERY supplied class potential, hence every actual harmonic
character at any modulus, is preserved. This uses isomorphisms rather than
inventing a large-order numerical electrical atlas.

The accepted phase group has a useful exact decomposition:

    Z/M x Z/K -> {(a,b):a=b mod3},
    (u,v) -> (u+3v,u).                                  (1)

Its inverse is u=b and v=(a-b)/3 modulo K; the latter quotient is defined
by the congruence modulo M and retains all K choices. In these coordinates
the two native moves are independent generators:

    g:(a,b)->(a+1,b+1)       = (u,v)->(u+1,v),
    h:(a,b)->(a+3,b)         = (u,v)->(u,v+1).             (2)

Both are full-domain actual source automorphisms, preserve s,t, and
transport all interior owners explicitly. They do not fix arbitrary
interior owner pins. Such extra pins must be imposed and the available
generators recomputed. The g move changes two cells, and h changes one
when K>1; respectively 2M and M actual owner clocks change. For K=1,
h is the identity. No kernel coordinate is discarded.

## 2. Exact spectrum: connectivity is now paid quantitatively

Give g and its inverse probability p each, h and its inverse probability
e each, and the identity probability 1-2p-2e, where p,e>0 and p+e<=1/4.
Repeated identity proposals when K=1 are counted as holding events. This
is a lazy, uniform-reversible native chain. Characters of (1) give all
eigenvalues, with multiplicities:

    lambda_(j,l) = 1 - 2p[1-cos(2pi j/M)]
                       - 2e[1-cos(2pi l/K)].            (3)

Each term is nonnegative in the energy 1-lambda. The smallest positive
one-coordinate energy occurs at frequency +/-1. Consequently, for K>1,

    gap = min(2p[1-cos(2pi/M)], 2e[1-cos(2pi/K)]).       (4)

For K=1 omit the second term. Formula (3) follows by direct substitution
of the characters into the four translations, and the characters form
a complete orthogonal basis of this finite group. It is not a quotient
spectrum or a floating-point eigenvalue estimate.

With p=e=1/8 the smaller term is the M-cycle term, so

    gap = (1/4)[1-cos(2pi/M)] >= 2/M^2.                 (5)

The inequality follows from sin(pi/M)>=2/M for M>=2. It retains the full
kernel: without h, there are K separate components and the global gap
is zero for k>=2. If e decreases with scale, (4) records precisely the
new slow direction; full generator rank alone would miss it.

This also shows why the fixed ternary quotient does not supply a
scale-independent gap. The positive bound (5) decreases as M grows,
despite every native generator moving at most two cells. In owner-clock
units those generators are not bounded-size as M grows.

## 3. Paid return routing and distinct owners

For an ordered pair of states, use the unique shortest signed u displacement
on the odd M-cycle, then the unique shortest signed v displacement on the
odd K-cycle. Its length is at most

    D=floor(M/2)+floor(K/2).

For k>=2 this is 2M/3-1; for k=1 it is one. These are actual generator
events, not quotient returns. The mean number of positive-direction steps
for a uniformly chosen displacement on an odd c-cycle is

    (c^2-1)/(8c).

Translation invariance makes the resulting directed-edge load independent
of its starting state. Majorizing each route length by D therefore gives
all-pair length-weighted congestion at most

    A = D max((M^2-1)/(8Mp), (K^2-1)/(8Ke)),             (6)

with the second term zero when K=1. This can also be derived simply by
counting positive and negative signed displacements separately. It uses
the actual N=M K states, whose factors cancel against edge capacity p/N
or e/N; it does not put uniform measure on a ternary quotient.

For ANY uniform-stationary return kernel on this same family, choose an
independent intermediate uniformly in Z/M x Z/K and concatenate the route
to it with the reverse route from the desired endpoint. The argument in
PROOF.md section 6 applies literally: at most 2D events and congestion
at most 4A. For p=e=1/8, (6) is D(M^2-1)/M. Actual changed-owner work
is at most 4MD per compiled return; the original excursion clock is
retained and the new clock charged separately. For positive actual weights
with ratio kappa, Metropolis and the same intermediate sampler give the
bound 4kappa A and gap at least (2/M^2)/kappa.

Complete unit-capacity families also meet the permutation compiler's
hypotheses, now with this smaller D and N=M^2/3. Every slot can retain its
actual target class throughout, not merely a frozen mod-eight value.
Heterogeneous slot eligibility and selected subsets still need the exact
matching and shared-variable tests; the incoming required three-slot
cycle is not replaced by unrestricted pair swaps.

## 4. Why this is useful, and the surviving boundary

The final incoming [stable-padding theorem](../session44_guarded_native_moves_20261023/stable-padding/THEOREM.md)
extends endpoint-fixed partner isomorphisms through transitive prefixes
and suffixes. Applying that proved closure here gives a quantitative
corollary at EVERY source order n>=8. Choose

    M=3+6 floor((n-8)/12),

and pad the 2M+2-owner source with the remaining zero through eleven
transitive owners. The generators fix both old endpoints and every new
owner; the full path and partner transports extend. The same M^2/3-state
rotation family and spectrum survive, and

    gap >= 2/M^2 >= 8/(n-2)^2.

The actual clock shift from prefix padding is retained. This does not
bound changed-owner work independently of n. The incoming padding packet
also already provides genuine TWO-SLOT CROSS-ORBIT moves at every n>=8;
we do not claim first existence of such moves. Our additional statement
is a quantitative growing rotation family, while universal connectors
across different accepted orbits remain unproved.

The source-to-product map (1) retains the actual kernel missing from the
ternary quotient. The target predicate is legal native return with all
class-based guards preserved. The map preserves that predicate because
the full source/partner isomorphisms certify both generators. The lost
information under reduction modulo3 is v and the higher digits of u;
they reappear as the slow Fourier modes in (3). The cheapest hostile is
to remove h: the quotient remains connected but K components survive.

Thus this is an all-size quantitative certificate for an ACTUAL
guard-preserving automorphism family, complementing the count-guard
cross-layer routing theorem. It does not solve the cross-orbit connector
problem, prove a uniform gap as the source grows, or show that every
accepted original electrical path enters this carrier.

`carry.py` checks literal paths, both source-automorphism generators, partner
transport and converse; separately checks the product coordinates and
retained components; and computes exact rational route loads for the
M=3,9,15,27 cases. M=81 checks all 2187 paths and their product chart, with
explicitly selected partner transports and route pairs. The receipt marks
that pair universe as incomplete. The all-k spectrum is proved by (3),
not inferred from this finite census.
