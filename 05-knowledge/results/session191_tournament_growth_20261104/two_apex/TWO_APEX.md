# A four-letter two-apex grammar and its exact continuation dimension

**PROVED transfer; finite exact linear-algebra certificate; finite coverage only.**
This is a different construction architecture from the same-clock rooted
wrapper. It permits an arbitrary attachment word rather than a fixed core and
one varying block length. The global spectrum remains open.

## 1. Physical objects and the inherited boundary

Fix two named vertices a,b with a->b. A word w=c_1...c_m over {0,1,2,3}
defines T(w) by adding the transitive chain x_1->...->x_m. Bit 0 of c_i says
that a->x_i, and bit 1 says that b->x_i. Thus

    0: x_i -> a,b       1: a -> x_i -> b
    2: b -> x_i -> a    3: a,b -> x_i.

Every tournament with a specified two-vertex transitive deletion has such a
presentation after naming the two deleted vertices in their arc order.
One-apex tournaments are contained by using a redundant apex. This statement
concerns presentations, not unique isomorphism representatives.

There is also a short exact strongness test. For a nonempty word, T(w) is
strong if and only if

    c_1 != 0, c_m != 3, and
    (some letter is 2, or a letter 3 occurs before a later letter 0).   (1)

Indeed a->b is fixed. A path from b to a either uses a letter 2 by itself,
or starts at a chain name beaten by b and ends at a later chain name beating
a. If no letter 2 exists, that is precisely a 3-before-0 pair. Once a,b lie
in the same strong component, that component reaches the first chain name
exactly when it is not 0, and the last chain name reaches the component
exactly when it is not 3. The increasing chain then brings every name into
that component. If the first is 0 it is a global source; if the last is 3 it
is a global sink. The empty word is the nonstrong edge a->b.

The [one-apex result](../one_apex/PROOF.md) proves that H=13 is impossible when one
vertex deletion is transitive, but exhibits a five-vertex two-apex witness.
This grammar crosses that genuine all-order one-apex obstruction. In the
present naming, the word 120 has H=13, with outgoing rows (6,8,26,17,3) on five vertices.

## 2. The exact 36-state transfer

Every Hamilton word contains the two apexes in exactly one order (p,q), either
(a,b) or (b,a). It has three intervening chain runs:

    L_0, p, L_1, q, L_2.

Each chain name is assigned to exactly one run; within a run its order is
forced. Empty runs are allowed. This is a bijection, subject precisely to the
boundary edges below.

For each of the two apex orders retain 18 states (x,y,z), where:

- x in {0,1,2}: run L_0 is empty, or nonempty with its last name pointing to p,
  or nonempty with its last name not pointing to p;
- y in {0,1,2}: run L_1 is empty, or nonempty with its last name pointing to q,
  or nonempty with its last name not pointing to q;
- z in {0,1}: run L_2 is empty or nonempty.

A partial assignment never violates a first-name boundary. Its current last
boundary may be bad, since a later chain name can repair it. This is why bad
last states must not be deleted while reading the word.

Initialize one empty state in each apex order. Read the chain from left to
right. For the current letter put P=[p->x_i] and Q=[q->x_i]. Every old state
contributes its count to each allowed assignment below, adding multiplicities
when targets coincide:

    assign to L_0: (x,y,z) -> (2 if P else 1, y, z);
    assign to L_1: allowed iff y!=0 or P=1;
                   target (x, 2 if Q else 1, z);
    assign to L_2: allowed iff z!=0 or Q=1;
                   target (x,y,1).

At the end accept exactly x!=2 and either y=1 or (y=0 and p->q).
There is no final condition on z. In particular empty middle runs are accepted
only for the apex order (a,b).

Summing these accepting states is **exactly H(T(w))**, for every finite word,
including the empty word. The proof is the run-assignment bijection: assigning
a name to a nonempty run preserves every old first boundary and replaces its
last boundary; assigning its first name checks the incoming apex edge; final
acceptance checks the two remaining last boundaries and the empty middle
edge. No other adjacency can fail because each chain run is increasing.

Let M_c be the resulting 36 by 36 nonnegative integer transition matrix,
v_empty the two empty states, and f the accepting row. With column states,

    H(T(c_1...c_m)) = f M_(c_m)...M_(c_1) v_empty.

These formulas themselves define the matrices; `verify.py` generates them
rather than reading an opaque fitted matrix file.

## 3. Exactly 17 rational linear continuation coordinates are needed

The reachable vector space over Q has dimension 28. After quotienting vectors
indistinguishable by **all suffix words**, its dimension is exactly 17.
Thus H on this full word language has a 17-dimensional rational linear
representation and no representation of smaller rational linear dimension.
This does not forbid nonlinear encodings, restricted sublanguages, or a
finite-state grammar whose separate states carry other coordinates.

The certificate is finite linear algebra with an all-word closure check:

1. Starting from v_empty, `verify.py` adds a vector only when it is independent
   of the existing span, and tests every M_c image of every added vector.
   Its 28 vectors therefore span an invariant reachable space, not merely a
   bounded-depth sample.
2. On that space, start from f and prepend each letter to each suffix row,
   again continuing until every image is in the existing span. Exactly 17
   independent observable rows remain.
3. Their independent columns give the 17 by 17 matrix
   [H(T(uv))], with prefixes u and suffixes v generated in `receipt.json`.
   Its determinant is 150994944 = 2^24 * 9, nonzero. Every one of the 289
   entries is independently recomputed by vertex-subset Hamilton-path DP.

For completeness, a d-dimensional linear representation factors every such
Hankel matrix through Q^d, so its rank is at most d. The displayed minor gives
d>=17. The invariant observable quotient gives d<=17. This proves the exact
claim, without extrapolating a finite word census.

## 4. Variable clocks and the boundary of the grammar

On the entire reachable space every letter matrix satisfies

    M_c^4 - 6 M_c^3 + 11 M_c^2 - 6 M_c = 0.

The verifier checks this exact identity on each of the 28 spanning vectors.
Consequently, for every fixed prefix u, suffix v and letter c, there are
rational A,B,C such that, for all k>=1,

    H(T(u c^k v)) = A + B*2^k + C*3^k.

The identity has distinct roots 0,1,2,3, so the positive-clock solutions have
no polynomial-in-k factor. For k=0 the zero eigenvalue can contribute an
additional boundary term. This is an all-clock statement, not a fit to a few
samples. Multiple independently variable runs give a multivariate sum of
products of these bases, separately in each positive run length. For example,

    H(T(12 0^k)) = -1 + 4*2^k + 2*3^k,       k>=1.

The first three exact counts are 13,33,85, which determine the three
coefficients; the proven recurrence then establishes the formula for every
k. These tournaments are strong: a reaches the first chain name; every chain
name reaches the final 0-name and then a,b; and b reaches the second chain
name and hence a. Thus this is an explicit infinite strong two-apex family
starting at the all-order one-apex obstruction 13. It is not an infinite
arithmetic progression or arbitrary-target construction.

A fixed finite set of run patterns nevertheless cannot be a universal
spectrum construction, even after arbitrary order joins: see the proved
multiplicative sparsity theorem in `SPARSITY.md`. The candidate route left
open here is unbounded changes of attachment letters. No `3H+2 / 3H+4`
closure, target solver, or infinite interval theorem is asserted.

## 5. Reproduction and finite scope

    python3 verify.py
    python3 -O verify.py

No external packages or input bank are needed. Each run checks every one of
the 1365 presentations of chain length at most five independently by both
literal run assignment and full vertex-subset DP; 40 deterministic larger
samples; all 289 Hankel entries by DP; invariant-space closure; and every
clock-polynomial identity. The JSON receipt is deterministic and all checks
use explicit exceptions.

The standalone C++ scan `scan.cpp`, driven by `run_scan.py`, enumerates all
22,369,621 words of chain length at most 12 with no isomorphism filtering.
It retains 200,235 distinct values, maximum 744,405, and an unbroken odd
prefix through 230,519 except 7 and 21. The first additional missing value in
**that finite scan** is 230,521. This is not an all-order absence.

Run `python3 run_scan.py` to regenerate the complete witness TSV and scan
receipt. The driver independently checks all 499 witnesses through 1001 by
full vertex-subset DP, plus 20 deterministic larger witnesses. Its complete
witness TSV SHA-256 is pinned in `scan_receipt.json`; the smaller independently
checked witnesses are generated in `prefix_witnesses.json`. The full finite prefix was also checked independently with a literal-slot subset-zeta calculation, not vertex-subset DP of all 200,235 witnesses. Every integer is below the 64-bit safety bound: at maximum
chain length 12 there are at most 2*3^12 complete run assignments, and a state
at any prefix counts no more assignments than that.

The larger prefix is finite computational evidence for a possible two-apex
route, not a uniform arithmetic induction. The proof of the exact transfer,
minimal response dimension, and clock identities does not depend on this
large scan.
The inherited fixed-feedback presentation and run-slot ideas are in S26's
`session26-fixed-apex/proof.md`; its random-attachment results count average
mass rather than this individual-word transfer. No literature-priority claim
is made. The new object here is an explicit exact continuation language and
its minimal rational linear response dimension.
