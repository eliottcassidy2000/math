# Every finite closure system is a denominator bank

**Status: PROVED representation and optimality statements; CITED only for the
standard Dilworth equality; FINITE-EXACT small certificates.** This extends
[the guard adjunction and dihedral banks](NOTE.md). It does not make guard
observations into a dynamical state, does not make arbitrary root-system
symmetries into tournament isomorphisms, and makes no literature-priority claim.

## 1. Universal squarefree realization

Let E be a finite labeled ground set and c a closure operator on E. Let C be
its closed sets, and M the proper meet-irreducible closed sets. Here a proper
closed set M is meet-irreducible when M=A intersection B for closed A,B forces
M=A or M=B. Equivalently it has exactly one immediate larger closed set.

Assign a different prime p_M to each M in this family, and set

    b_i = product_(M : i not in M) p_M,       i in E.          (1)

Empty products are 1. Then for every subset S of E,

    {j : b_j divides lcm_(i in S) b_i} = c(S).               (2)

**Proof.** Every closed set is the intersection of the meet-irreducible closed
sets containing it. To see this without invoking a lattice representation
result, for j outside a closed set X choose a maximal closed M containing X
but excluding j. Every strict closed superset of M contains j, so M cannot
be the intersection of two such supersets. This M is meet-irreducible and
separates X from j.

The prime p_M is absent from lcm_(i in S) b_i exactly when S is contained in M.
Thus b_j divides that lcm exactly when j belongs to every M containing S.
By the separation argument this is exactly j in c(S). The same reasoning
covers S empty and closures with nonempty c(empty). ∎

So every finite closure system, not just a special divisibility geometry,
has an exact static arithmetic guard realization. Choosing ordered fractions
1/b_i makes these actual reduced FRACTRAN denominators, but their transitions
are not prescribed by the closure system.

## 2. The exact minimum number of primes

The squarefree construction uses |M| primes, and **this is minimum among
squarefree representations** of c. More generally, if arbitrary exponents
are permitted, the minimum number of distinct primes is

    minimum number of inclusion chains covering M
      = width(M).                                           (3)

The last equality is the standard finite
[Dilworth theorem](https://mathweb.ucsd.edu/~asuk/dilworth.pdf).
The arithmetic representation and lower bound are proved here.

**Lower bound.** In any denominator representation define, for every prime p
and positive integer t up to its maximum occurring exponent,

    X_(p,t) = {i : v_p(b_i) < t}.                            (4)

Each X_(p,t) is c-closed. Indeed the lcm of its denominators still has p
exponent less than t, so no index outside it can enter its closure. For a
fixed p these sets form an inclusion chain. Moreover every c-closed S is an
intersection of the X_(p,t) that contain S: any j outside S has some p with
v_p(b_j) greater than the maximum exponent of p on S, so a threshold in (4)
separates it. This remains valid for S empty, using maximum 0.

Every proper meet-irreducible M must therefore equal some X_(p,t). If an
intersection representation used only strict supersets, finite iteration of
meet-irreducibility would give a contradiction. Consequently one prime can
account only for a chain of meet-irreducibles, proving the chain-cover lower
bound. In the squarefree case it has only the nontrivial threshold t=1 and
can account for at most one M, giving |M|.

**Upper bound.** Partition M into k chains. For a chain

    M_1 strictly contained in ... strictly contained in M_h

choose one prime p, and give index i exponent

    e_i = min{t-1 : i in M_t},  or h if no such t exists.      (5)

Then {i:e_i<t}=M_t for 1<=t<=h. Multiplying these prime powers across all
chains gives denominators whose threshold-closed sets include every member
of M. They generate exactly c by intersections, proving (2) for the new code.
This uses k primes. If M is empty, all b_i=1 and zero primes suffice. ∎

A further exact cost bound is

    sum_p max_i v_p(b_i) >= |M|.

Every meet-irreducible needs its own prime threshold. The chain construction
attains equality: a chain of length h has maximum exponent exactly h because
its largest closed set is proper. Thus compressing to width(M) primes need
not increase the total number of exponent thresholds.

This is an exact coordinate-cost statement: primes are independent exponent
registers, and exponent thresholds in one register can encode only a chain.
It is not a bound on integer bit-length, numerator size, execution time or
minimal tournament order.

For interval closure on m>=2 ordered indices, M consists of the m-1 proper
nonempty prefixes and the m-1 proper nonempty suffixes. Its width is 2.
The familiar code b_i=p^(m-i) q^(i-1) from NOTE.md is therefore prime-optimal,
while any squarefree code needs 2m-2 distinct primes.

## 3. Type A gives exact permutohedral arithmetic banks

Take E={(i,j):1<=i<j<=r}. Define c by repeatedly applying

    (i,j),(j,k) in S  implies  (i,k) in S,       i<j<k.        (6)

This is the transitive closure on naturally ordered pairs. Under
(i,j) -> e_i-e_j it is the usual two-root positive closure for type A_(r-1).
No root-system classification is needed for the next statement.

A subset S and its complement are both closed under (6) exactly when S is
the inversion set of a total order of {1,...,r}. To prove this, orient the
pair i<j backward when (i,j) belongs to S and forward otherwise. On each
triple i<j<k, the only directed cycles correspond to

    (i,j),(j,k) in S but (i,k) not in S,
    (i,j),(j,k) not in S but (i,k) in S.

These are precisely the two forbidden closure patterns. Thus biclosure is
equivalent to having no directed triangle. A tournament with no directed
triangle is transitive: a shortest directed cycle of greater length has a
chord producing a shorter directed cycle. Its vertex order is unique.
Conversely every total order satisfies both implications.

Applying the universal code to (6) therefore gives exactly r! arithmetic
observations whose complements are also realizable. Two such states differ
in one guard bit exactly when the total orders differ by interchanging two
consecutive vertices: reversing a nonconsecutive pair would violate its
comparison with an intervening vertex. Thus the guard-flip graph is the
permutohedron graph on S_r, not merely an abstract r!-element set. Complement
reverses the total order.

Now apply NOTE.md's reflection-orbit embedding with g=r(r-1)/2 guards. The
bank produces exactly r!/2 individual blue lines, for every r>=2, using any
n for which floor((n-1)^2/4)>=g. This is a statement about line identities and
marked presentations. It is not a count of distinct endpoint classes, nor
a formula for the endpoints' Hamilton-path counts.

The intermediate r-vertex comparison tournaments are all transitive and
hence all isomorphic when labels are forgotten. This gives an extreme
quotient warning: all r! different arithmetic chamber states disappear into
one unmarked comparison-tournament class. The larger actual blue embedding
is a separate construction that retains the labeled guard bits.

## 4. Explicit prime-optimal A_3 code

For r=4 order the guards as

    12,13,14,23,24,34.

One prime-optimal denominator list is

    1089, 4235, 330, 875, 100, 24.                            (7)

The only primes are 2,3,5,7,11. Its realizable observations are exactly the 40
transitively closed subsets of those six pairs. Precisely 24 are biclosed,
forming the S_4 permutohedron. In particular the reduced ordered program
whose fractions are the reciprocals of (7) has that exact static bank.
No transition correspondence with adjacent swaps is asserted.

The 11 proper meet-irreducibles, as masks in the stated pair order, admit the
five-chain partition

    (11,15,31), (26,62), (33,39,55), (53), (56,60).

The five masks

    15,26,39,53,60

are pairwise incomparable. The equal chain-partition and antichain sizes are
an independently checkable certificate that five primes are necessary and
sufficient, without relying on a matching algorithm's optimality claim.
A squarefree code needs 11 primes. The general representation from (1)
provides it, but (7) spends exponent depth to reduce the prime count.

Exact small results are:

| r | Guards | Closed observations | Biclosed observations | Squarefree minimum primes | General minimum primes |
|---|---:|---:|---:|---:|---:|
| 2 | 1 | 2 | 2 | 1 | 1 |
| 3 | 3 | 7 | 6 | 4 | 2 |
| 4 | 6 | 40 | 24 | 11 | 5 |
| 5 | 10 | 357 | 120 | 26 | 10 |

The final two columns beyond these rows are not extrapolated.

## 5. Scope and finite controls

The exact connection has source a finite closure system, target a denominator
bank, map (1) or (5), and preserved predicates closure/biclosure and all guard
bits. Prime names and exponent labels are additional implementation choices.
The operation preserves neither arbitrary dynamics nor preferred priorities.
The powers-of-four obstruction in NOTE.md still rules out a finite exact
observer for its stated countdown program.

Run `python3 universal_check.py` and `python3 -O universal_check.py` here.
Both produce byte-identical output. The checker enumerates all 61 closure
systems on a labeled three-element ground set, including those with nonempty
closure of empty, and tests both code constructions on every subset. It then
checks every subset for the type-A systems r=2,...,5, compares biclosed sets
with independently generated permutation inversion sets, and compares every
one-bit edge with an adjacent-swap edge. Each minimum-prime claim has an
explicit equal-sized chain partition and antichain, both independently
verified. `universal_certificates.json` records those certificates and codes.

The root-system connection is fully proved here in type A and in NOTE.md's
rank-two case. More general finite Coxeter inversion-set banks can be obtained
by combining the universal representation with the standard finite-biclosed
inversion-set theorem, as recorded by
[Barkley--Speyer](https://arxiv.org/abs/2207.05998). That broader corollary is
**CITED**, not required for any proof or finite claim above, and no general
root-system result is inferred for arbitrary denominator families.
