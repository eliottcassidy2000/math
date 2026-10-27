# A uniform prefix-rerouting obstruction for pure-blue strong tournaments

Status: PROVED, independently audited, 2026-10-20.
The proof is elementary once the audited odd-circulant reduction is imported.
No large-order search or Hamiltonian-cycle counting bound is used.

## Precise target and inherited input

A tournament is **pure-blue** if the positional reversal of every directed
Hamiltonian path is an anti-automorphism. The audited note
`PURE_BLUE_STRONG_CYCLE_REDUCTION.md` in the companion foundations packet
proves that a nontrivial strong pure-blue tournament has odd order and is
a circulant along every directed Hamiltonian cycle. That reduction uses
Camion Hamiltonicity and Lu's crossing theorem; those cited inputs are
not reproved here.

We prove that **no odd circulant tournament of order at least seven is
pure-blue**. Combined with that reduction and the verified positive
controls C3 and the regular C5, this classifies the nontrivial strong
pure-blue atoms as exactly C3 and C5.

Write n=2m+1. On Z/nZ let x -> y mean y-x is in the skew connection set S:
exactly one of d,-d lies in S for each nonzero d. The Hamiltonian-cycle
normalization means 1 is in S. Put s_d=+ if d is in S, for 1<=d<=m.
An arbitrary circulant can also be put in this normalization by
replacing every label x by -x if necessary, since exactly one of
1 and -1 belongs to S.
All examples here already have the directed Hamiltonian path
P0=(0,1,...,2m); its reversal is r(i)=2m-i.

The theorem does not identify purity with circulancy, vertex transitivity,
regularity, or path rigidity. C5 is pure-blue but nonrigid.

## 1. Disjoint-support prefix lemma

Suppose q is a permutation of {0,...,m} fixing m. Extend q by the identity
on {m+1,...,2m}, and suppose

Pq=(q(0),q(1),...,q(m),m+1,...,2m)

is a directed Hamiltonian path. If T were pure-blue, then both r and
q r q^{-1} would be anti-automorphisms. Consequently

g=(q r q^{-1})r=q (r q^{-1}r)

would be an automorphism.

The support of q is contained in {0,...,m-1}; the support of
r q^{-1}r is contained in {m+1,...,2m}. These supports are disjoint.
Therefore the order of g equals the order of q.

No tournament has an even-order automorphism: taking the odd part of
such an element's order and then a suitable power gives a nontrivial
involution, and any transposed pair of vertices has its unique arc
reversed by that involution. Thus:

**Every such directed prefix rerouting q has odd order in a pure-blue
tournament.**

It suffices below to display one even-order q. A displayed finite word
is extended by the unchanged increasing tail. All displayed changed
positions lie strictly below m.

## 2. Three explicit forbidden sign configurations

Assume m>=4.

### F1: s2=- and s3=+

Use

q=(2,0,3,1,4,5,...,m).

The changed path increments are -2,3,-2,3. They are allowed by the
hypotheses. As a permutation of positions, q has the 4-cycle
(0,2,3,1), so its order is four.

### F2: k>=4, k<=m, s_k=+ and s_(k-2)=-

Use

q=(k-1,1,2,...,k-2,0,k,k+1,...,m).

Only 0 and k-1 are exchanged, so q is a transposition and fixes m.
Its non-unit path increments are -(k-2),-(k-2),k. Hence it is directed
under precisely the stated sign hypotheses, with no requirement on
s_(k-1).

### F3: the first negative sign is k, where 3<=k<=m-1

Thus s1,...,s_(k-1) are positive and s_k is negative.

If k is odd, use

q=(1,k,0,2,3,...,k-1,k+1,k+2,...,m).

The non-unit increments are k-1,-k,2,2, so the path is directed.
The changed positions form the single cycle
(0,1,k,k-1,...,2), of even length k+1.

If k is even, use

q=(1,3,4,...,k,0,2,k+1,k+2,...,m).

The non-unit increments are 2,-k,2,k-1. The permutation has exactly
the two nontrivial cycles

(0,1,3,5,...,k-1), of length k/2+1,

(2,4,6,...,k), of length k/2.

One of these consecutive integers is even; their least common multiple
is therefore even. Again q fixes m because k<=m-1.

### Direct violated-arc certificates, without group theory

There is an independent, more direct check of every family. On the low
half {0,...,m-1}, the reversal of Pq has the explicit formula

r_Pq(x)=2m-q^{-1}(x),

because q fixes the reflected high half. For F1 and both forms of F3,
the original arc 0 -> 1 fails under anti-reversal:

| Family | r_Pq(0) | r_Pq(1) | Required image increment r_Pq(0)-r_Pq(1) |
|---|---|---|---|
| F1 | 2m-1 | 2m-3 | 2, absent since s2=- |
| F3, k odd | 2m-2 | 2m | -2, absent since s2=+ |
| F3, k even | 2m-k+1 | 2m | -(k-1), absent since s_(k-1)=+ |

For F2, take whichever of 0 -> k-1 or k-1 -> 0 is the original
arc. The permutation q exchanges these two vertices, so their
required reversed image has exactly the opposite increment and is
absent from S. Explicitly r_Pq(0)=2m-k+1 and r_Pq(k-1)=2m.

These four certificates prove the same forbidden configurations
without using automorphism groups or permutation orders at all.
The disjoint-support lemma explains the discovery and generalizes
the obstruction, but is not logically necessary for this particular
classification once these direct arc checks are supplied.

## 3. The three-pattern reduction

Every pure-blue normalized circulant of order 2m+1, m>=4, must satisfy
exactly one of:

A. S intersect {1,...,m} = {1};

B. S intersect {1,...,m} = {1,...,m-1};

C. S intersect {1,...,m} = {1,...,m}.

Proof. If s2 is negative, F1 makes s3 negative. F2 says that
s_k positive implies s_(k-2) positive for every 4<=k<=m.
Induction separately on the two parities makes every s2,...,s_m
negative. This is A.

If s2 is positive, a first negative at any index 3<=k<=m-1 would
give F3. Thus s1,...,s_(m-1) are all positive; the last sign is free,
giving B or C. These cases are necessary, not yet sufficient.

## 4. Unit-step reindexing eliminates all three patterns

For any unit a modulo n with a in S, the order
(0,a,2a,...,(n-1)a) is a directed Hamiltonian cycle. In these coordinates
the new connection set is

S_a={d : ad mod n is in S}.

Purity is invariant under this relabelling. In particular the three-pattern
reduction applies again. Because n is odd, 2 and -2 are always units;
there is no primeness assumption or hidden failure at composite order.

### Case A

Here S={1} union {m+1,...,2m-1}. Take a=-2=2m-1, which belongs to S.
The new sign s'_2 is positive: -4 mod n=2m-3 is in S for m>=4.
Set k=floor(m/2)+1. Then 3<=k<=m-1, and

  - if m is even, -2k mod n=m-1;
  - if m is odd, -2k mod n=m.

These residues are neither 1 nor in {m+1,...,2m-1}, so s'_k is
negative. This contradicts the three-pattern reduction, whose only
patterns with s'_2 positive have every sign through m-1 positive.

### Case C

Here S={1,...,m}. Take a=2 in S. The new s'_2 is positive because
4<=m. Set k=floor(m/2)+1, again 3<=k<=m-1. Then
m<2k<=m+2<n, so s'_k is negative. This is the same contradiction.

### Case B

Here S={1,...,m-1} union {m+1}. Take a=2 in S.

For m>=5, s'_2 is positive because 4<=m-1. If m is even,
put k=m/2; then 3<=k<=m-1 and 2k=m is absent from S.
If m is odd, put k=(m+3)/2; then 4<=k<=m-1 and 2k=m+3 is
absent from S. Thus a negative sign occurs before m, again impossible.

For m=4, S={1,2,3,5} in Z/9Z. Reindexing by a=2 gives
S_2={1,5,6,7}, which is precisely case A at m=4. Case A has already
been ruled out. Explicitly the composite reindexing used to reach a
forbidden prefix is by 2*(-2)=-4 modulo nine.

Therefore no circulant tournament of odd order n>=9 is pure-blue.

## 5. The order-seven boundary

The normalized connection sets (1 in S) are exactly the four below.
Each displayed directed Hamiltonian path has a positional reversal
that fails to reverse at least one arc. All can be checked directly.
Here the listed pair (u,v) satisfies u->v but r_P(v) does not point
to r_P(u).

| S | Directed Hamiltonian path P | Violated arc (u,v) | Absent required image arc |
|---|---|---|---|
| {1,4,5} | (0,1,2,6,4,5,3) | (0,5) | 1 -> 3 |
| {1,2,4} | (0,1,2,4,5,6,3) | (0,2) | 5 -> 3 |
| {1,3,5} | (0,1,2,5,3,6,4) | (2,0) | 4 -> 3 |
| {1,2,3} | (0,1,3,4,6,2,5) | (0,2) | 1 -> 5 |

These words are imported from the independent Hamiltonian-cycle certificate
packet `cycle_certificates.json`, then
checked again by the present script. They are a four-case hand proof,
not an extrapolation from a larger census.

## 6. Positive controls and final scope

At order three the unique strong tournament C3 is pure-blue: its three
Hamiltonian paths are rotations of (0,1,2), with reflections as their
positional reversals.

At order five the two normalized skew connection sets {1,2} and {1,3}
give isomorphic regular C5 tournaments, by multiplication by 3.
For S={1,2}, translation lets us restrict to Hamiltonian paths starting
at zero. A direct two-choice branching gives exactly

(0,1,2,3,4), (0,2,3,4,1), (0,2,4,1,3).

Their positional reversals are respectively x -> 4-x, x -> 1-x,
and x -> 3-x modulo five. All three are anti-automorphisms of the
circulant. Translating recovers all 15 Hamiltonian paths and proves
purity. The accompanying script independently enumerates these small
path sets, rather than assuming that circulancy alone gives purity.

Together with the inherited odd-circulant reduction:

**The nontrivial strong pure-blue tournaments, up to isomorphism,
are exactly C3 and the regular C5.** A singleton may be included when
the term atom allows the trivial strong component.

The proof does not invoke the new path-rigidity classification, a
Hamiltonian-cycle counting lower bound, or the finite sweep through31.
Its distinctive coordinate is the order of a relative prefix
permutation, preserved by the disjoint-support product of two
path-reversal anti-automorphisms. The direct violated-arc table then
removes even that group-theoretic dependency from the explicit proof.

## Reproduction and independent checks

Run `python3 verify_prefix_classification.py --max-exhaustive 31`.
The script constructs the explicit prefix witness (with at most two
unit reindexings) and checks directedness and an actual failed
anti-automorphism pair. It also independently checks permutation
orders and the disjoint-support identity, the four order-seven words,
and all Hamiltonian paths in the order-three and order-five controls.
Random larger examples and hostile composite structured sets are
included only as checks of the implementation, not proof dependencies.

The pre-proof probes `probe.py` and `probe_prefix.py` explain the
discovery route; they are not needed to reproduce the theorem.

The separate audit implementation `independent_prefix_obstruction.py`
imports none of this script. It verified all 32,760 normalized presentations
of orders 9 through 31, the four order-seven certificates, all small
positive-control Hamiltonian paths, and the three residual patterns at
every m from 4 through 200. The symbolic audit checked the full parameter
ranges and the direct violated-arc table, not merely the finite tests.
