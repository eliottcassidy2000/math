# A one-arc strong size repair with the whole rooted response retained

**PROVED:** arbitrary-module neutral repair, its physical doubling square,
exact primitive guard and finite target selector. **FINITE-EXACT:** the
small quotient coefficients and displayed prime fixtures. **OPEN:** an
all-target selector, a general nonzero-digit return, and Chan–Pak.
The scalar frame15p−h is inherited; the result here is an actual
root-star-preserving one-arc/one-vertex repair between different strong
orders, with equal full rooted response for every nonempty module.

## 1. Recovery and the change of question

[S158's square-sum Zeckendorf code](../../session158_creative_20261103/square_overlap/ZECKENDORF.md)
and [S165's strong code](../../session165_geometry_representation_20261103/geometric_carry/STRONG_ZECKENDORF.md)
encode additive Fibonacci addresses into products of(1+2^j). They prove
injectivity and actual strong guard returns, not additive completeness of
the attained counts. The S165 code even has an explicit sparse2-adic image.

The scalar frame15P_M(2)−H(M) and its exact ternary children are already
in [prime_ternary_lift CORE_REPLACEMENTS §2](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md).
That packet explicitly does not preserve endpoint profiles. Our new claim
is the full rooted equality, named size repair, and commuting growth below.

The useful inherited physical interface is
[S164 GUARD_ALGEBRA §§5–6](../../session164_emergent_20261103/profile_algebras/GUARD_ALGEBRA.md).
The current methodological comparison is
[S312 MIXED_SHORTCUTS](../../session312_cluster_pythagorean_symmetry_20261113/MIXED_SHORTCUTS.md):
its Pythagorean short/long relation is legal on every seed, rather than
just an equality of unguarded matrices. We ask for that same distinction
at a tournament port.

The nearest response calculus is
[THM-4865 / sole-exit tower](../../sole_exit_tower_tribonacci_20261104/PROOF.md),
whose missing variable is a controlled rooted response at a chosen target.
The exact consumer used below is
[ROOTED_RESETS R1](../../two_port_return_induction_20261103/ROOTED_RESETS.md):
equal rooted four-profiles become equal **all-future two-port** responses
after a common fresh seal. Old arbitrary marks are not silently retained.

The live concepts were: additive addresses versus multiplicative codes;
physical polar ports; exact short/long words; uniform-module run bounds;
and target-dependent core replacement. The two first exploratory routes
and their stopping boundaries are recorded after the main result.

## 2. The literal repair and its arbitrary module

Use out-neighbor bitmasks with names0,...,5:

    S=(46,12,56,48,35,2), root r=4.

Reverse precisely the old arc1→2. Then add a new named vertex6 which beats
only1. The result on0,...,6 is

    T=(110,8,122,112,99,66,2), root r=4.              (1)

No old arc at r changes. Its incoming neighbors are2 and3 on both sides;
it gains only an outgoing arc to the new vertex6. The intermediate flip
is strong too. A Hamilton cycle common to S and the flipped six-vertex
skeleton is

    0,2,5,1,3,4,0;

in T it becomes0,2,5,6,1,3,4,0.

Replace vertex0 in **both** skeletons by the same arbitrary nonempty
tournament M, with all incidences inherited uniformly. Write

    h=H(M), p=P_M(2), c=C_2(M)=p−2h,

where P_M(2) counts covers in two labelled slots allowing empty slots,
and C_2 counts two nonempty ordered paths. Call the resulting graphs
S[M],T[M]. Their orders are |M|+5 and |M|+6. A Hamilton path of M expands
the displayed skeleton cycles, so both graphs and the intermediate flip
are strong even when M itself is not strong.

The named module has only one incoming outside neighbor, r. It beats all
other outside vertices in either frame. The repair still changes only
old arc1→2 and adds6→1. It touches no arc inside M.

### Universal rooted response identity

Use z=(H,s,q,F), with s counting HPs starting at r, q=H after deleting r,
and F counting ordered spanning covers(P,Q) with P nonempty ending at r
and Q possibly empty. Then

    z(S[M],r)=z(T[M],r)
      =(29h+15c,5h,5h,35h+15c)
      =(15p−h,5h,5h,15p+5h).                       (2)

This is not an extrapolation from transitive inflation. Here is a complete
finite coefficient certificate for the arbitrary-module statement.

A single Hamilton path can visit M in at most two runs: every run except
a possible initial run must enter through r, which occurs only once.
A path starting at r uses exactly one M run, immediately after r. Deleting
r makes M a global source module, again forcing one run. For F, slot1
ends at r, so the only possible external incoming arc r→M cannot be used.
Consequently each of the two slots contains at most one M run, at its
beginning. Thus **all four observations have at most two M runs**.

Collapse the ordered nonempty M runs into occurrences of0. For k runs,
their fillings are counted by C_k(M), ordered by global slot and then by
position in that slot. The complete quotient-word counts are:

| observation | one M run | two M runs |
|---|---:|---:|
|H|29|15|
|s|5|0|
|q|5|0|
|F|35|15|

The table is identical for S and T. `lift_check.py` independently counts
these finite quotient words, including the empty second-slot case, without
using a tournament inflation or fitting unknown coefficients. Its universe
is all permutations of the five/six outside names and k identical zeros,
with k=1,2. For F the cut is forced immediately after r; every other seam
must be an actual directed arc. The run bound proves there are no omitted
coefficients. This establishes(2) for every tournament M.

The module position is essential. Inflating old vertex1 by TT2 instead
of vertex0 gives source/target states(95,11,11,105) and(71,11,11,85).
Scalar equality of the singleton pair alone would not justify inflation.

The intermediate flipped skeleton has the equally complete coefficient
rows(17,9),(3,0),(3,0),(21,9). Consequently its full rooted state is

    (9p−h,3h,3h,9p+3h).                             (2a)

The single arc reversal therefore implements the decreasing target map
H→(3H−2h)/5 on this actual frame. The sole-exit helper then restores
15p−h and the entire initial four-profile. Both stages are strong, and
under(3) both intermediate and terminal rooted states are primitive.
For h=1 the examples are29→17→29,59→35→59,119→71→119.
This is a frame-guarded rational descent; no arbitrary target is declared
to have that frame or that module.

## 3. A complete future interface, and its exact content

After either repair endpoint, add the same new named seal vertex z with

    r→z→all other vertices,

and retain the ordered pair(r,z). By ROOTED_RESETS R1, equation(2) implies
equality of **every** future H response under the six-command two-port
language. This includes changing the active root between available ports
and future newborns. It excludes returning to discarded older names.
The extra helper6 is not silently made an available third port.
The exact content statement uses ROOTED_RESETS R3, formula R5a: the
first four sealed query columns have determinant−4. Their dyadic inverse
and the oddness of all tournament counts give equality of the sealed
all-future gcd with the original rooted gcd; rank alone would not suffice.

The original rooted content is exactly

    gcd(H,s,q,F)=gcd(h,15p).                         (3)

Indeed gcd(s,q,F−H)=gcd(5h,6h)=h. Intersect with H=15p−h.
Thus the sealed interface is primitive if and only if gcd(h,15p)=1.
This is an iff statement, not an inference from rigidity or strongness.
For M=C3 the common content is3; for the order4 H5 core it is5.

The rooted vertex r has indegree two. In particular the identity
F=H+q−s valid at an indegree-one root does **not** apply here:
F−H=6h, while q−s=0. The actual four-profile and fresh seal are the
retained resource.

## 4. The neutral repair commutes with an actual doubling digit

A root-child birth at r, retaining r, produces a new vertex x with
r→x→M and x beating every other old vertex. Hence x and M together form
the homogeneous module TT1→M in exactly the same skeleton. Literally,
up to the declared new name at the start of that module,

    R(S[M],r)=S[TT1→M],
    R(T[M],r)=T[TT1→M].                             (4)

The one-arc repair and this birth commute: they use the same old names1,2,4,
and the repair's vertex6 lies after the whole enlarged module.
Ordinal joins multiply every labelled-slot evaluation, so h stays fixed
and p doubles. Therefore(4) has the exact marked law

    (H,5h,5h,H+6h) → (2H+h,5h,5h,2H+7h).           (5)

For transitive M, h=1 and this is an actual2H+1 digit. The known scalar
ray15·2^m−1 is not claimed as new. What is added is a commuting physical
size repair at **every** stage, preserving all rooted responses and the
old root star while changing the strong order by one.

More generally M itself may grow by actual ordinal insertions. All old
vertices and arcs of M persist, as do all external port incidences. The
neutral repair can be applied or undone in either growing frame. Undoing
it means deleting its specifically named helper6 and reversing the same
arc; it is not an inferred deletion of an arbitrary graph vertex.

The forward neutral step S[M]→T[M] cannot simply be repeated to create
unbounded strong order at a fixed H. Its inverse is available as the
explicit helper deletion and reverse flip. Indefinite growth uses(4),
or changes M itself, and changes H accordingly.

This is a source-closed family with exact digits h. It does not permit
choosing h independently of the actual module whose H is h and whose
second response is p. That joint realizability is still the target gate.

## 5. Unbounded cores, a finite inverse, and exact seven phases

Let

    D=(6,28,24,17,1), H(D)=11, P_D(2)=100.

Use the actual module

    M_(n,m)=D^(ordinal n)→TT_m,
    n,m≥0, n+m>0.

The profiles multiply slotwise on ordinal joins, giving h=11^n and
p=100^n2^m. Equation(3) is always1. Thus the two strong orders
5n+m+5 and5n+m+6 have equal primitive rooted state and count

    H_(n,m)=15·100^n·2^m−11^n.                      (6)

Letting n grow changes the nontrivial core itself; letting m grow gives
the literal doubling digit(5). Neither operation alters the old module
internally. No statement of prime infinitude is made.

This language has an exact finite target inverse. For a requested Y,
only n with14·100^n≤Y are possible, since H_(n,m)≥14·100^n2^m.
For each such n, test whether

    (Y+11^n)/(15·100^n)

is an integral nonnegative power2^m, with n+m>0. Every passing address
constructs both physical witnesses, and every address in(6) is tested.
Uniqueness of the address is not required or asserted.

The following seven-adic calculation was supplied by the parent lane.
Modulo7,

    7|H_(n,m) iff m≡n (mod3).

Write m=n+3k, allowing integer k subject to m≥0. Modulo49,
200/11≡1+35, 15≡1+14, 8≡1+7. Factoring the unit11^n gives

    v_7(H_(n,m))=1 iff 2+5n+k is nonzero modulo7.     (7)

Negative k uses the inverse of8 modulo49; the same first-order expansion
holds. This is a divisor selector for an already realized family, not a
prime-cofactor conclusion.

Small exact fixtures, each with both strong orders and identical rooted
response, are:

|n,m|H|factor status|orders|
|---|---:|---|---|
|0,3|119|7·17|8,9|
|0,6|959|7·137|11,12|
|1,0|1489|prime|10,11|
|1,7|191989|7·27427|17,18|
|2,5|4799879|7·685697|20,21|

The prime/cofactor checks use trial division through the integer square
root. Counts in the larger fixtures use(2) and ordinal multiplication,
not a claimed20-vertex path census. All values in(6) are4 modulo5, so
this family alone cannot select every target.

## 6. The first physical short/long signal and its missing port

For an indegree-one guard, S164 uses nonnegative state(s,u,d), with H=s+u,
q=s+d and F=H+d. Its actual one-vertex operations are

    C(s,u,d)=(u,s+u+d,s),
    D(s,u,d)=(s+d,2u,0),

where C selects the newborn and D retains the old guard. Each preserves
strongness from a nontrivial strong input. On every canonical state
v=t(1,2,0), the two legal words

    S0=D^8 C D,             length10,
    L0=C^4 D(CD)^3,         length11,

both return513v. For S0 this is1+2^9. For L0, C^4 sends(1,2,0) to
(11,20,6), D gives(17,40,0), and the successive CD blocks give57,171,513
times(1,2,0). Prefixing DCD sends any guarded input to
(s+2u+d)(1,2,0), so DCD S0 and DCD L0 have **identical universal
three-state matrices**. This is a legal all-seed identity, not a formal
matrix relation outside the guard.

After k choices of S0 or L0, a canonical seed has the same rooted
profile at every one of the k+1 orders n+10k+j,0≤j≤k. We do not infer
2^k distinct graph classes from2^k histories. A common fresh seal then
makes the complete future two-port interface equal by ROOTED_RESETS R1.

Reusing the old last child without that seal is false. At scale t, the
last child's deletion count is1025t after S0 and855t after L0, despite
the equality at the active guard. The receipt checks the full unequal
child profiles. Moreover the returned rooted content is513t: this
size-adjustment relation by itself traps every future sealed count in a
large divisor. It cannot be a primitive prime-return theorem. The
arbitrary-module repair(2) is the positive escape from that limitation.

## 7. Why simply adding the opposite guard did not supply digits

A cheap orthogonal test retained an indegree-one guard r and an
outdegree-one guard t, r≠t. Right births at r and converse left births at
t physically commute and preserve the opposite guard; their cross arcs
always point from the right newborn to the left newborn. Each side may
retain its old guard or select its newborn. This is a genuine enlarged
physical language, not a forced numerical analogy.

It still cannot attain11 or13 from C3. Every unique-entry birth satisfies
H'=H+u+d≥3H/2, because u≥s+d: u is the ending-cover response at the old
predecessor, which is at least the deletion path count s+d by splitting
an HP at that predecessor. The exit case is the converse. The first
birth from C3 gives5; either next guarded choice gives9 (the two possible
order4 guard states are(1,4,0),(2,3,1)). A third birth gives at least15,
and all later counts grow. Thus11 and13 are absent at every depth.

`polar_probe.py/json` records the bounded discovery control through order9;
the all-depth stopping statement is the preceding inequality, not that
finite list. A new core or a new operation is required even after both
polarities are allowed. This was the reason to look for the physical
core replacement(1), rather than extend the same list.

## 8. Reproduction and scope

    python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/selection/lift_check.py

`lift_receipt.json` contains the complete quotient coefficient table,
selected literal transitive and nontransitive modules, strongness and
old-root-star controls, actual named doubling embeddings, primitive and
nonprimitive controls, the arithmetic phase tests and prime fixtures,
and the short/long child-port hostile.

`response_collision.py/json` preserves the discovery from the inherited
ordinary strong bank through order7; it is not a new class census.
`module_probe.py/json` preserves the five possible nonroot inflation
positions and exposes why the module slot matters. The main theorem uses
the run proof and complete tiny quotient table, not the scan's absence
of other hits.

The remaining selection gate is explicit: given Y, find an actual module
with15P_M(2)−H(M)=Y (or a different controlled frame), then use the repair
without losing its marked continuation. The construction solves the
physical size/port compatibility for every such module. It does not solve
the arbitrary-Y module equation.
