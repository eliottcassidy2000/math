# Root-child interfaces, a four-coordinate grammar and amplifying queries

**PROVED scoped factorization and observation theorems; FINITE-EXACT bounded
grammar coverage. Prime/7p realization at all sizes remains OPEN.**
The closest mechanisms are the inherited seventeen-query observer and rooted
four-profile transfer. The hostile is the equal-profile H37 pair when an
older third name is used. The corrected near miss is to mistake a low-rank
individual command for a small observer under every command. The least-used
sidecar here is the moment at which the other old name is deliberately
discarded.

## R1. A universal root-child interface really has four coordinates

Let (T,v) be any nonempty rooted tournament, and write

    z=(h,s,q,F),
    u=(a,t,s,delta)=(h-s,q-s,s,F-h-q+s).

The definitions and nonnegative physical coordinates are inherited from
[ADAPTIVE sections2-3](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md).
Add x with v->x->T-v and retain the ordered ports(v,x). Call this operation
S(T,v). Let X be the seventeen continuation queries, in precisely the order
of [live_port.json](../adaptive_state_induction_20261103/live_port.json).
Then

    X(S(T,v);v,x)=u F_*,                              (R1)

where F_* is the explicit nonnegative integer4-by-17 matrix in
`rooted_resets.json`. It has rank4; its first four columns, corresponding to
queries empty,0,3,00, are

    [[2,4,3,6],
     [1,1,2,3],
     [1,1,3,3],
     [1,3,2,5]],  determinant -4.                    (R2)

Thus equality of the original rooted profiles is equivalent to equality
of **all future two-port H responses after this interface is formed**.
This is stronger than closure under keeping the root or marking the newest
vertex: subsequent choices may use any of the six two-port commands. It
does not allow returning to discarded older names.

**Proof.** For |T|>=2, temporarily forget all names except v. Put e equal
to the number of Hamilton paths ending at v. The only relevant endpoint
cover signatures, with forgotten labels denoted*, are

    (*,*),  (*,v),  (v,*),
    {(*,*),(v,v)},  {(*,*),(*,v)}.

Their respective coefficients are

    h-s-e,  e,  s,  q,  F-e-q.                      (R3)

At most one forgotten end is allowed, since every future newborn dominates
every forgotten old vertex. Join the new x to at most one incoming component
ending at v and one outgoing component starting elsewhere; exclude choosing
the same component twice. This is the inherited endpoint-cover bijection.
The resulting five-by24 birth map has rank5. Composing with the exact
seventeen-query columns annihilates the e coefficient identically. Substituting
z in terms of u gives F_*. These are exact coefficient identities, not fits
to a finite source sample. The prior closed observer propagates them through
every future word; (R2) gives the converse and sharp four-dimensional bound.
Actual rooted source profiles span four dimensions, as the small controls
also certify, so the lower bound is not based on unattainable formal inputs.

The singleton boundary follows by adding a forgotten global sink to T before
the operation. This preserves h,s,q,F, including H(empty)=1. The sink stays
a global sink under every subsequent two-port command and hence preserves
all H observations. The |T|>=2 argument therefore applies. Literal singleton
controls are also included in the receipt.

In the six-letter encoding `3*root_slot+discard`, the two interface commands
are1 and3: keep the chosen root and the newborn, discard the other old
name. Their17-by17 maps have rank4. Commands0,2,4,5 have rank13. After3 we
rename the controller slots so the old root is slot0; this renaming is
compiled into later letters and is not an unrecorded graph operation.

## R2. An explicit finite positive grammar for a new induction attempt

Starting from a root-child interface, take zero, one or two commands from
{0,2,4,5}, then finish with1 or3 and normalize the controller slots. For a
word w of this form, exact factorization gives

    F_* A_w = B_w F_*,                              (R4)

including the slot swap after a terminal3. The resulting42 words have
**28 distinct nonnegative integer4-by4 maps B_w**. Every map is supplied
in the receipt with its word, rank and determinant. Thus their action is
on actual named sources, not just a cone of numerical profiles.

The two nonsingular maps are

    w=041: [[7,1,2,1],[4,2,0,1],[4,2,1,2],[5,1,1,1]], det4;
    w=403: [[6,3,1,4],[3,0,1,2],[4,0,1,2],[5,2,1,2]], det-4.

The other26 maps have rank at most3. These full-rank excursions recover
possibilities that fresh-block-only marking loses. No assertion that the
two invertible maps alone cover targets is made.

The boundary state u describes the rooted parent **before** its final
interface newborn. To compile a sequence of maps, create that provisional
newborn, perform the literal words, and omit the last newborn when outputting
the final parent. Equivalently perform the complete raw program except its
last birth. The retained old root is already named in that parent. The
compiler records every raw letter and slot permutation; no arbitrary
deletion or graph substitution is needed. Each macro adds |w| vertices to
the output parent, preserves strongness from the C3 seed, and strictly
increases its H count.

Initial root-child interfaces are obtained from C3 by a prefix of length
at most2 over{0,2,4,5}, then a terminal1 or3. Their parent profiles give
**21 distinct actual entry states**, each with its source program. This
entry bank and the28 maps are a concrete finite positive system proposed
for induction. Its exact all-word source semantics are proved by(R4) and
the compiler. Its full numerical image is OPEN.

One complete bounded saturation gives:

| Initial prefix limit | Entry profiles | States at parent H<=999 | Missing allowed odd values |
|---|---:|---:|---|
|1|5|32499|55,67|
|2|21|39329|none|

Every accepted scalar witness in the second row is decoded and independently
recounted on its literal parent graph. Rejection above999 is permanent by
strict H growth; there is no hidden vertex-order cutoff. This does not
extend the inherited numerical record. Its significance is the finite
21-entry/28-map **four-coordinate** language, substantially smaller than
the unrestricted two-port state saturation at this same bound.

The next constructive theorem must find an invariant actual subfamily of
this grammar and show that each target has a smaller admissible predecessor
or a proved reset route to one. Inverting only the two nonsingular maps
would omit26 source-producing branches. An auxiliary inverse parent may
be composite even when the target is prime or7p.

## R3. The exact unavoidable continuation divisor

For any actual two-port source in the observer's domain, set

    D(T;ports)=gcd(X_1,...,X_17).

Every X_j is odd. Every coefficient of the six query-update matrices has
denominator1 or2, as checked over Q in the inherited exact certificate.
Consequently D divides H after **every** legal future two-port word: the
dyadic denominators are units modulo this odd D. Conversely each defining
query is a permitted future word. Therefore

    D = gcd{H(after w): w any six-command word}.       (R5)

It follows that D(current) divides D(next). For a strong source of order at
least three, if D>1, no later nonempty adaptive continuation can have prime
count: a prime future count would have to equal D, but strict H growth and
D<=H(current) prevent that. The strongness guard is necessary: TT2->C3,
with both ports in TT2, has D=H=3 and a birth at its universal source can
preserve H3. The universal divisor theorem still applies to that weak source.
A7p target must
retain the necessary divisibility D|7p; this alone is not sufficient.

The static seventeen-query basis used by the native-reset probe has the
same divisor. [The basis certificate](observer_induction_basis.json) gives
`X_static=X_adaptive B`, det B=-1/4, and denominator LCM2 for both B and
B inverse, with identities on all25 universal carrier coordinates. Each
odd gcd therefore divides the other. No odd-index ambiguity is hidden in
the choice of seventeen probes.

At a root-child interface there is a simpler exact formula:

    D(S(T,v))=gcd(a,t,s,delta)=gcd(h,s,q,F).           (R5a)

The integer factor F_* proves one divisibility; the dyadic inverse of its
first four columns proves the reverse for the odd divisor D. The profile
gcd is odd because it divides h=a+s. The two coordinate lists are related
by an integral invertible change of basis. Thus (R5a) is valid for every
actual nonempty rooted source, without a strongness assumption.

For the native module-reversal family V_k in N7, forming the next interface
at root3 gives D=3 again, while forming it at the oldest leaf z gives D=1.
Indeed the former profile is(h,3,3,h); for the latter, gcd(h,q_z)=3 but
F_z=3*2^k+2*3^k-1 is not divisible by3. The choice of root is an exact
arithmetic guard. Primitivity removes this obstruction without supplying
a target-coverage theorem. Spectrum independently audited this corollary.

When a rooted profile is b times another actual rooted profile, (R1)
transports that literal scaling to every two-port continuation after sealing.
This extends the earlier keep/new factor obstruction to the whole future
six-command language. The native resets in
[NATIVE_RESET_INDUCTION](NATIVE_RESET_INDUCTION.md) deliberately change X
and can decrease D while preserving the current H. They are outside the
six-command divisor-monotonicity theorem.

## R4. Every nonzero observer difference can be amplified

Let A_0,...,A_5 be the inherited17-by17 query matrices and f the column
reading H. For the two commands retaining both old ports put

    g_r=(A_(3r+2)^2-A_(3r+2))f,  r=0,1.

They read twice the rooted coefficient F_r-q_r. Exact closure finds
seventeen columns A_w g_r, with |w|<=4, forming a square matrix of
determinant128. The receipt supplies all words and every rational entry;
spectrum independently reconstructed the columns and determinant.

For any two unequal observer rows X,Y, their difference d has nonzero
pairing with at least one of these columns. Append that prefix w and then
keep adding at its indicated root. The inherited fixed-root formula gives,
for k>=1,

    H_X(w r^k)-H_Y(w r^k)=alpha*2^k+beta,
    2alpha=d A_w g_r !=0.                            (R6)

Thus every unequal observer pair has unbounded discrepancy under a short
prefix followed by a constant-root tail. This is an all-word consequence
of the exact finite basis, not a sampled growth claim.

In particular no nonzero finite set of error vectors can be closed under
all six right actions while retaining bounded H corrections. A uniform
ternary replacement whose error is required to stay in a finite bounded
register under **every** continuation therefore fails. Selective command
languages, unbounded integer carries, source-dependent choices and native
resets are not excluded. These are precisely the additional options the
constructive proof must retain.

## Verification and connection boundary

`rooted_resets.py/json` certifies the five-signature factor, nonnegative
four-coordinate matrix, all28 macro factorizations, dyadic coefficients,
the17-column amplification basis, literal small-source continuation controls
and the two explicitly bounded saturations. It uses the previous positive
kernel and a separate literal adjacency DP path for the scalar controls.
There are646 literal small-source continuation checks,609 direct rooted
entry/macro compiler checks, and497 decoded scalar witnesses in the complete
second saturation. A separate `observer_induction_reset_audit.py/json`
recovers the factor from four actual rooted sources, checks all42 words and
replays44 controller programs, without importing this producer.
Normal and optimized runs must agree. No unrestricted spectrum, minimal
positive realization or arbitrary native-operation closure is inferred.

The connection is rooted actual graph -> paired covers -> sealed queries ->
finite macro programs. It preserves the chosen future H response and its
integer divisor, loses unused old names, and restores a graph with a literal
word and slot-transport sidecar. The cheap hostile is a required inverse
state outside the actual entry language, followed by the marked H5 example
in [OBSERVER_INDUCTION](OBSERVER_INDUCTION.md). The native reset family
separately tests whether a source can escape a permanent divisor while
returning to a legal marked interface.
