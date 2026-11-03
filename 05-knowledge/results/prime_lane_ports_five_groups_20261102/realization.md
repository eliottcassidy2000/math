# Exact component relocation and constructive prime-lane exchanges

**PROVED all-size identities and monotonicity; FINITE-EXACT independent
implementation controls. Prime output is an additional arithmetic test.
Uniform realization of every prime and every required 7p remains OPEN.**
No claim of external novelty or a new finite coverage record is made.

## Inheritance, research board, and the changed operation

Source snapshot: `48298050de`. The closest proved mechanism is the actual
marked-parent block action in
[`prime_inverse_ports_20261102/BRANCHES.md` §1](../prime_inverse_ports_20261102/BRANCHES.md),
with the uniform-endpoint specialization in its §3. The canonical hostile
is the sharp order-five pair
[`U=1->1->C3`, `V=1->C3->1`](../tournament_joint_parity_20261102/growth_obstruction.md):
they have identical H, automorphism order and full spanning-cover profile,
but different one-vertex attachment distributions. The corrected near miss
is inferring an actual source from an algebraically legal inverse state;
[`spectrum_descent_jacobian_20261102`](../spectrum_descent_jacobian_20261102/README.md)
exhibits the impossible state (9,7,3,5) despite attainability of scalar 9.
The least-used sidecar is **the order of one marked singleton relative to
a vertex-transitive SCC**, with all apex ports held fixed.

Anchor: the two realization lanes in
[THM-4094, prime/7p completeness reduction](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
Niche: a constructive operation on the exact information lost by the cover
profile. Wildcard: the quadratic form controlling a relocation commutator.
The live board is actual marked sources / SCC order / endpoint response /
component symmetry / dyadic output image / prime divisibility.

New incoming inputs were recovered before deriving this result:

* [`uniform_prime_transport_20261103/THEORY.md`](../uniform_prime_transport_20261103/THEORY.md)
  gives homogeneous-module run-content obstructions and native six-edge
  exchanges, as well as exponential prime-candidate families. It does not
  prove primality for every parameter or universal prime entry.
* [Session150's triangle-wrapper transfer](../session150_hamiltonian_spectrum_20261103/transfer/TRANSFER.md)
  gives H(C3[T,1,1])=3C_1(T)+C_2(T) and the full cover transfer. That is a
  homogeneous cyclic wrapper. The present operation moves an existing
  component across a singleton under a **nonuniform apex cut**.
* [Session150's square-root family](../session150_hamiltonian_spectrum_20261103/induction/RESULT.md)
  is an injective lower bound on distinct attainable values, not all-prime
  coverage. Its prime counts and the older finite banks remain different
  scopes.
* The current cited unrestricted positive bank reaches 99,999,999, as routed
  by the descent and uniform-transport packets. It is inherited and was not
  rerun here; the present result does not extend it.

The general sorting theorem R0 below was supplied by the coordinator after
the vertex-transitive calculation and independently audited here.

The new graph operation keeps every old named vertex, every internal block
edge and every apex-to-old-vertex edge. It reverses only the complete cut
between the adjacent singleton and block. This is an explicit construction
on actual tournaments. It is **not asserted to be a legal word in the
native complement/re-marking language**; those separate guards must not be
erased by calling both operations “exchange.”

## 1. The actual four-coordinate state

A parent tournament A and a marked vertex x carry a fixed subset S of A,
with x->v exactly for v in S. Write

```
F=H(A+x),
h=H(A),
g=#{Hamilton paths of A ending in S},
q=#{ordered spanning two-slot covers (P,Q) of A:
      Q nonempty and first(Q) in S}.
```

The first slot may be empty. The empty parent has state (1,0,0,1).
All coordinates are literal nonnegative integers, with h>=1 and 0<=g<=h.
Append a nonempty block B after A in the dominance order, retaining its
actual subset S_B beaten by x. The inherited action is

```
F'=HF+Pq+(H-s)g+(f-H)h,
q'=Eq+Qh,
g'=(H-t)h,
h'=Hh.                                                   (1)
```

Here H=H(B), E counts ordered two-slot covers of B allowing empty slots;
P counts covers whose first slot is nonempty and ends outside S_B;
Q counts covers whose second slot is nonempty and starts in S_B;
s counts Hamilton paths starting in S_B; t counts those ending outside S_B;
and f=H(B+x). The scalar state is used to evaluate the specified actual
blocks; it is not a license to substitute arbitrary nonnegative states or
profile values.

For a singleton with cut bit e in {0,1}, denote its matrix by S_e:

```
S_e = [[1,1-e,1-e,0],
       [0,2,  0,  e],
       [0,0,  0,  e],
       [0,0,  0,  1]].                                   (2)
```

Matrix products act on column states, with the rightmost factor applied
first. Thus S_e M_B means append B and then the singleton.

## 2. Strict sorting for every nontrivial actual block

The coordinator generalized the first vertex-transitive result to arbitrary
blocks. The following stronger theorem and proof were independently checked
here; the more precise rank-two increment remains useful in the vertex-transitive
case below.

**Theorem R0 (PROVED).** Let B be any tournament with n>=2 vertices and an
actual cut S_B. In any actual prefix and common suffix context:

* moving a zero-cut singleton from before B to after B strictly increases H
  if S_B is nonempty;
* moving a one-cut singleton from after B to before B strictly increases H
  if S_B is not all of B.

In the complementary empty/full boundary, the relevant entire state matrices
commute. No regularity, transitivity or automorphism hypothesis is needed.

**Proof.** Use H,E,P,Q,s,t,f from (1), and put D=E-H-P. Direct multiplication
with the actual block matrix M gives

```
S_0 M-M S_0 = [[0,D,-s,H+Q-t],
               [0,0, 0,    Q],
               [0,0, 0,  t-H],
               [0,0, 0,    0]],

M S_1-S_1 M = [[0,P,s-H,H+P-s],
               [0,0,  0,E-H-Q],
               [0,0,  0,   -t],
               [0,0,  0,    0]].                       (R0)
```

D counts covers whose first slot is nonempty and ends in S_B; E-H-Q counts
covers whose second slot is nonempty and starts outside S_B. Both are
nonnegative. Also Q-s counts proper two-path covers whose second path starts
in S_B. For each v in S_B take that second path to be (v), and choose a
Hamilton path of B-v as the first. Thus Q-s>=|S_B|. Similarly P-t>=n-|S_B|
by taking the first path to be a singleton outside S_B.

The common suffix's first matrix row (alpha,beta,gamma,delta) satisfies
alpha>0, beta>=0 and 0<=gamma<=alpha, for the reason proved in (6) below.
Using 0<=g<=h, contraction of the first matrix in (R0) is bounded below by

```
alpha*D*q + alpha*(Q-s)*h + beta*Q*h,
```

strictly positive when S_B is nonempty. The second is bounded below by

```
alpha*P*q + alpha*(P-t)*h + beta*(E-H-Q)*h,
```

strictly positive when S_B is not full. Empty/full cuts give zero matrices
in the matching cases. Every equality uses the actual shared cut and source
state. QED.

This explains why loss of SCC order can affect arbitrary attachment even for
asymmetric blocks. It supplies extremal directions of controlled relocation,
not an assertion that all intermediate counts occur.

## 3. An all-size relocation law for vertex-transitive blocks

Let B be a nontrivial vertex-transitive tournament of order n, and put

```
u=H(B)/n,           d=(E(B)-H(B))/n.
```

Uniform endpoint counts make u,d positive integers. Since
E(B)=C_2(B)+2H(B), where C_2 counts nonempty ordered two-path covers, we have
d>u for n>=2. Let a=|S_B|. The literal cut pattern may affect f=H(B+x), but
transitivity gives

```
P=(n-a)d, Q=ad, s=au, H-t=au.
```

Consequently

```
M_B = [[nu,(n-a)d,(n-a)u,f-nu],
       [ 0,n(u+d),   0,    ad],
       [ 0,     0,   0,    au],
       [ 0,     0,   0,    nu]].                          (3)
```

There is no assertion that arbitrary cuts of the same size give the same
f. The important fact is that f cancels in the following exact difference.

**Theorem R1 (PROVED).** Define

```
K(u,d) = [[0,d,-u,d+u],
          [0,0, 0,  d],
          [0,0, 0, -u],
          [0,0, 0,  0]].
```

Then for e in {0,1},

```
S_e M_B - M_B S_e = (a-ne) K(u,d).                       (4)
```

In particular, after any actual prefix state (F,q,g,h), the two orders
have difference

```
Delta(F,q,g,h)
  = (a-ne) (d q-u g+(d+u)h, d h, -u h, 0),             (5)
```

where Delta means “B then singleton” minus “singleton then B.”
This is direct multiplication of (2)–(3), so it holds for every actual
cut pattern, including non-arc-transitive B.

The source parent's H and full path-cover profile are unchanged by this
reordering: the source SCC multiset is unchanged, and
[THM-3166, falling-factorial order-join transform](../../../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md)
makes that profile commutative under order join. The marked child H changes
by the explicitly retained response in (5). This identifies what the
commutative quotient lost.

### Arbitrary common suffixes preserve the direction

Let W be any actual common nonempty block word appended after the two
orders, with fixed named cuts; the empty suffix is allowed. Write the first
row of its product as (alpha,beta,gamma,delta). Then

```
alpha>0, beta>=0, 0<=gamma<=alpha.                       (6)
```

Indeed alpha is the product of block H values. The input g coordinate can
contribute only at the suffix's first block: q,h do not depend on g, and
that block resets g independently of the input. Hence, if its H and s are
H_1,s_1, then gamma/alpha=(H_1-s_1)/H_1 in [0,1]. All relevant matrix
entries are nonnegative (f>=H by insertion monotonicity), proving beta>=0.
For the empty word alpha=1,beta=gamma=0.

The final H difference is

```
(a-ne) [alpha(dq-ug+(d+u)h)+(beta*d-gamma*u)h].           (7)
```

The bracket is strictly positive because g<=h, d>u, q>=0, h>=1 and (6):

```
bracket >= alpha*d*q + alpha*(d-u)*h > 0.
```

**Corollary (PROVED).** In every actual prefix and suffix context:

* moving a zero-cut singleton from before B to after B strictly increases
  H when a>0;
* moving a one-cut singleton from after B to before B strictly increases
  H when a<n;
* when a=ne, the entire four-coordinate state agrees in both orders.

These are actual graph inequalities with fixed cuts, not merely average
or modular statements. They do not say the values between the endpoints
are attained.

## 4. The quadratic connection is an exact matrix identity

With E_14 denoting the matrix having one at row1,column4 and zero elsewhere,

```
K(u,d)^2=(u^2+d^2)E_14,       K(u,d)^3=0.                (8)
```

This is the same quadratic form appearing in the inherited two consecutive
prime-symmetric-event residue law in
[`BRANCHES.md` §3](../prime_inverse_ports_20261102/BRANCHES.md).
There it is the coefficient of the two unpaid boundary transitions modulo
the prime module size. Here it is the second-order coefficient of the
relocation-difference operator. The shared map is the concrete matrix K,
not the numerical equality of a dimension or a prime.

The difference operator K is signed and need not itself map actual states
to actual states. Equation (8) is therefore an algebraic connection, not
an additional positive constructor. The two words on the two sides of (4)
are the actual constructions. No signed linear combination is accepted as
a tournament.

For a triangle, (n,u,d)=(3,1,3), giving

```
S_e C_a-C_a S_e=(a-3e)K(1,3),    K(1,3)^2=10 E_14.
```

This specializes the previously finite order-five obstruction to an
all-context operation law.

## 5. Four explicit strong dyadic families and a lane exchange

Fix r>=1. On the same named vertices take a transitive prefix R of r
vertices, a singleton z, and a directed triangle C. Let the marked vertex
x beat every vertex of R and exactly a vertices of C, where a is 1 or 2;
let z beat x. Fix the triangle orientation and the actual selected subset
of C once and for all. Define

```
U_(r,a): parent R -> z -> C, with the retained apex cut;
V_(r,a): parent R -> C -> z, with the same retained apex cut.
```

Both have r+5 vertices. Only the three z-to-C edges change. Every old
vertex and every x-port retains its identity. Their source parents have
H=3, automorphism order3, and the same full cover profile.

**Theorem R2 (PROVED).** These marked children are strong, and

| a | H(U_(r,a)) | H(V_(r,a)) | difference |
|---:|---:|---:|---:|
|1|15*2^r-7|18*2^r-7|3*2^r|
|2|9*2^r-1|15*2^r-1|6*2^r|

**Proof.** After the r prefix singletons, the state is
(1,2^r-1,1,1): the parent is transitive, x is a source, and every nonempty
choice for the second slot contributes one to q. A triangle with cut size
a has H3,E12,P9-3a,Q3a,s=a,t=3-a,f=3+a(3-a). Apply its matrix and S_0
in the two orders to obtain the displayed formulas. Strongness is literal:
x reaches the first parent SCC, every parent SCC reaches the last, and the
last SCC reaches x. In U the last SCC is the mixed-cut triangle; in V it is
z. Every vertex is therefore mutually reachable with x. QED.

Small witnesses illustrating both directions between the lanes are

```
r=1,a=1:   23 <-> 29,
r=2,a=2:   35=7*5 <-> 59,
r=3,a=2:   71 <-> 119=7*17,
r=5,a=2:  287=7*41 <-> 479.
```

The displayed primes were checked exactly. These are not new attainable
values beyond the existing bank. Their contribution is a common all-size
construction and explicit named-edge transport between the two endpoints.
Primality for general r is neither assumed nor proved.

### Exact support image, not just a forward formula

Membership in the union of these four supplied families is equivalent to
at least one of

```
(h+7)/15, (h+7)/18, (h+1)/9, (h+1)/15
```

being an integer power 2^r with r>=1, in the corresponding row and column.
The successful test recovers r, a, U/V, and the actual named tournament.
Thus this restricted image has an exact terminating inverse without a
separate source-entry debt. It has only O(log X) values below X and already
misses the attainable prime 13. This decisive hostile prevents promoting
it to a full prime or 7p compiler.

### Exact elementary divisibility branches

All four family values are 2 mod3. Both a=1 families avoid divisibility
by7, so they can contribute only to the prime lane among the two targets.
For a=2, the following are exact for every r>=1:

```
7 | H(U_(r,2)) iff r=2 mod3,
7 | H(V_(r,2)) iff r=0 mod3,
5 | H(U_(r,2)) iff r=2 mod4.
```

Hence the two 7p candidate branches are

```
U: r=3k+2, k>=0,  H/7=(36*8^k-1)/7;
V: r=3k,   k>=1,  H/7=(15*8^k-1)/7.                    (9)
```

On the U branch, r=2 mod12 forces a factor5 in H/7; only the first r=2
can have prime quotient there (it gives5). This removes an infinite
subbranch without losing a potentially prime exception. Also
5 | H(V_(r,1)) iff r=2 mod4, so that entire prime-candidate subbranch fails.
These are necessary arithmetic tests, not sufficient primality tests.
No assumption on the infinitude of primes in exponential sequences enters
any construction or claim.

## 6. Controls and remaining interface

`realization_probe.py` imports no repository implementation. It builds
literal ordered blocks and named apex cuts, counts H on all subsets by
last-vertex DP, and obtains q and g independently from actual endpoint
counts. It compares the complete four-coordinate state with the transfer
on all 723 words of at most four blocks and at most eight parent vertices,
using the alphabet S_0,S_1,C_0,C_1,C_2,C_3. It tests 2408 contextual swaps,
1098 actual block/cut profiles through order four with no symmetry filter,
156 general vertex-transitive matrix identities (including different f values to attack the
cancellation), and the nilpotent-square identity. The dyadic families are
checked directly as graphs through r=5, including strongness, and exactly
by their transfer through r=30. All family membership inversions and the
stated divisibility branches are tested, with13 as a negative control.
The general conclusions follow from the proofs, not the finite ranges.

Run from this directory:

```sh
python3 -I -B realization_probe.py > /tmp/realization.json
python3 -I -B -O realization_probe.py > /tmp/realization-optimized.json
cmp /tmp/realization.json /tmp/realization-optimized.json
```

Both runs agree; checks use explicit exceptions, not removable assertions.
No random weights, isomorphism-class resampling, reset of source ports, or
implicit clock change occurs.

The strongest remaining realization interface is now concrete: find a
positive family of actual marked blocks/cuts for which the exact inverse
image contains every required p and7p. The present operation supplies
controlled variation in the previously lost SCC-order coordinate, an
all-context sign theorem, and fully decoded sparse families. It supplies
neither an interval-filling theorem nor a choice of parameter for each prime.
