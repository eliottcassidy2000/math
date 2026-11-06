# A two-retained-port, one-state binary digit obstruction

Status: PROVED by exact finite certificate; independently audited.
This is not a disproof of F(8,137), and it does not produce a transducer.

## Statement

Retain two distinct active vertices, labeled 0 and 1, in an arbitrary actual
tournament. A digit operation is any one of the following 24 legal updates:

1. Introduce a vertex labeled 2. It beats an arbitrary subset A of the two
   old ports and loses to every other old vertex.
2. Forget any one of the three ports.
3. Label the surviving two ports 0 and 1 in either order.

Thus there are two retained ports between digits and peak register width
three during an update. Choose one fixed update M0 for digit 0 and one
fixed update M1 for digit 1. The updates need not be distinct.

Choose one fixed closing context from the following six:

- mark old port 0 or old port 1, forgetting the other port;
- introduce a final root beating any chosen subset of the two old ports
  and losing to every bulk vertex, then forget both old ports.

For the resulting rooted tournament let B=h+L-s and put D=3B+2s.

**Theorem.** No actual starting tournament and no real constant kappa make
both identities

    D(M0 U) = 2D(U) - kappa,
    D(M1 U) = 2D(U) + 6 - kappa

hold at every prefix U of every binary word.

No bound on the root-start count s is needed. The starting tournament may
have arbitrarily many vertices, and kappa is unrestricted. In particular,
this excludes the elementary fixed-context, one-state binary architecture
for the desired invariant D=6N+kappa, including large-prefix initializations
and affine shifts of the integer being encoded.

The theorem does not cover multiple control states, digit blocks containing
more than one update, different closing contexts selected by finite control,
more than two retained ports, or larger closing contexts.

## Exact forest presentation

Use the published p=2 directed path-forest presentation. A component records
its ordered pair of endpoints; -1 denotes a bulk endpoint and 0,1 denote
active ports. Components form an unordered multiset. Retain at most two
components whose start is bulk. Every active port is an endpoint of at most
one component, except that a singleton has endpoint pair (u,u).

The 46 formal patterns are generated in verify.py. Some formal patterns are
unreachable; this causes no problem because their physical counts are zero.
Introduction of a new vertex enumerates its possible predecessor and
successor components, forbids choosing the same component on both sides,
and merges the selected components with the new vertex. Forgetting changes
the forgotten endpoint label to -1 and discards patterns with more than two
bulk starts. A discarded cover can never contribute again: no future vertex
can enter a bulk-start component. These are exact nonnegative integer
linear transfers, not an arbitrary nonnegative-vector realization claim.

The matrices M0,M1 act on column vectors x of counts. For each closing
context, B and s are exact final linear readouts. A one-component cover
contributes to h; it contributes to s precisely when it starts at the root.
A one- or two-component cover contributes to L precisely when one component
ends at the root. The designated root-ending component is unique, including
in a two-component unordered cover. Consequently the stated B and D rows
are physical readouts.

The checker regenerates every update and closing row from these rules. Its
additional literal controls use a separate subset-partition path-forest DP.
Those finite controls are regression checks; the transfer argument above is
the all-order justification.

## Two necessary physical facts

We may assume the bulk is nonempty. If the initial tournament consists only
of the two active ports, pass to either one-digit child. It has three
vertices and two active ports, and inherits all required identities.

For nonempty bulk, the four patterns consisting only of active singletons,
including the empty pattern, have count zero. In the checker's ordering they
are coordinates 0,5,9,27.

Let e41 select the pattern

    {(-1,-1), (0,0), (1,1)}.

Its count is exactly H(bulk): the two ports are singleton components and the
remaining vertices form one spanning directed path. Every nonempty
tournament has a Hamilton path, so

    e41 x >= 1.

These facts apply to every actual seed of every order. They are the physical
constraints that eliminate spurious formal nonnegative vectors.

## Necessary invariant rows

For a chosen closing row ell=D, define

    r0 = ell M0 - 2ell,
    delta = ell(M1-M0).

The proposed digit recursion requires r0 x=-kappa and delta x=6 at every
reachable prefix x. Hence, for every word w,

    r0(Mw-I)x=0,
    delta(Mw-I)x=0.

Only words of lengths one through five are needed in the supplied
certificate. Stack their 124 rows, followed by the four impossible-pattern
coordinate rows, into a matrix H with 128 rows. Then Hx=0.

Let

    G = 6e41-delta.

Since delta x=6 and e41 x>=1, we have Gx>=0. Also x>=0 coordinatewise.
No bound, congruence condition, or assumed regularity of s is used.

## Exact spectral reduction

Each digit matrix must have eigenvalue 2. Indeed, let x0=M0 x and x1=M1 x.
The difference of their scalar readouts is 6. Applying either fixed digit i
repeatedly doubles this difference, because the same affine constant
cancels at every step. Therefore

    ell Mi^n (x1-x0) = 6*2^n   for every n>=0.

The characteristic-polynomial recurrence implies pi(2)=0. The checker uses
exact fraction-free determinants to test det(Mi-2I).

Precisely 14 of the 24 updates have eigenvalue 2, with indices

    0,1,2,3,4,5,8,10,11,13,16,17,22,23.

The other ten updates are excluded independently of the starting vector and
closing context. Operation indices are lexicographic triples
(mask, forgotten port, swap), where mask ranges 0..3, the forgotten port
ranges 0..2, and swap ranges 0..1.

## Exact Farkas certificate

For each of the six closing rows and each ordered pair of the 14 remaining
updates, CERTIFICATE.json supplies a rational row multiplier y and a
nonnegative rational scalar z satisfying the coordinatewise inequality

    H^T y - z G^T - delta^T >= 0.

There are exactly 6*14*14=1,176 cases. Multiplication by the actual
nonnegative column x gives

    0 <= -z Gx - delta x <= -6,

a contradiction. This proves the theorem.

The generator used floating-point linear programming only to discover
multipliers. The proof and verifier do not trust an LP solver, tolerances,
rank estimates, or rounded numbers. The released checker uses exact Python
integer and Fraction arithmetic, regenerates every constraint, verifies
all coordinate inequalities, and verifies that no case is missing. Sparse
rational multipliers are the complete mathematical certificate.

## Reproduction and provenance

Run:

    python verify.py
    python -O verify.py

Both runs must report PASS. The bundle contains the proof, exact certificate,
self-contained transfer/certificate checker, and independent literal-forest
control routine. The latter was supplied by the concurrent two-port seed
construction and uses direct subset DP rather than the transfer matrices.
The matrix semantics originate in the published finite-profile interface.

No repository checkout was modified by this scratch proof packet. Earlier
searches on bounded seed banks motivated the result, but their negative
outcomes are not used in the proof. The current obstruction is genuinely
all-seed within its explicitly limited architecture.
