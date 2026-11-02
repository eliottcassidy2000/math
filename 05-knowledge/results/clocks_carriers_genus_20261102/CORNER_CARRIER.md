# Exact clocked corner carriers and the re-marking boundary

Status: PROVED, with exact consequence checks and same-agent proof review;
no independent-agent or formal audit. This is a restricted native endpoint-edit grammar,
not a carrier for all metagraph operations or for tournament isomorphism.

## C1. The smallest carrier for the stated observation

Use THM-796's fixed directed path n -> ... -> 1 and free bits t_(i,j),
1<=j<i<=n, i-j>=2. Bit 0 means i->j. Fix an integer clock
0<=h<=n-3. The current output is the endpoint chord t_(n,1).

Operations are:

- complement all free bits, and path reflection, at zero clock cost;
- delete either path endpoint, at cost 1;
- attach a new endpoint with an arbitrary specified orientation row and
  the required path edge, at cost 1.
Every cost 1 operation changes h to h-1. Inputs retain their ordered ports.
All state counts are conditional on the shared current ordered vertex list
and clock: this external type information is not included in the bit count.
Relabelling to positions1..n is only a notation for that retained path.
This clock is a query budget, not a claim that one operation is a Collatz step.

Define the triangular address

    c_h(t)=(c_(a,b): a,b>=0, a+b<=h),
    c_(a,b)=t_(n-b,1+a).                              (C1)

**Theorem C1.** For fixed n,h, the coarsest exact deterministic carrier for
all endpoint-chord observations under these operations is c_h. It has
exactly 2^((h+1)(h+2)/2) states, independent of n once n>=h+3.

**Sufficiency and actual attachment compiler.** Complement flips every c;
reflection exchanges a,b. Here reflection means reversing the ordered
labels together with taking the tournament converse, so the forced path
remains directed. After deleting the low endpoint,

    c'_(a,b)=c_(a+1,b);

after deleting the high endpoint, c'_(a,b)=c_(a,b+1).
After attaching a high endpoint, let u_a be the specified new chord to
old vertex1+a. Then, at the new clock h-1,

    c'_(a,0)=u_a,       0<=a<=h-1;
    c'_(a,b)=c_(a,b-1), b>=1, a+b<=h-1.               (C2)

Low-end attachment is the reflected formula: c'_(0,b)=v_b and
c'_(a,b)=c_(a-1,b) for a>=1. Thus only h bits of a full attachment row
can affect the remaining observations. The path edge is forced separately;
these h queried chords never coincide with it. All updates preserve source
and named ports, and all needed old coordinates are already in (C1).
Induction on the remaining budget proves sufficiency, including arbitrary
interleavings of zero-cost complement/reflection.

**Necessity.** If two corners differ at (a,b), delete the low endpoint a
times and the high endpoint b times, then read the chord. This uses at most
h steps and distinguishes the states. Every corner bit assignment extends
to a full fixed-path tournament independently. Hence the displayed state
count is both sufficient and necessary. At h=n-3 all free bits are recovered.

This is an explicit all-order extension for an endpoint-chord query, not an
extension of THM-796's much richer node/line/colour observation census.
For n100,h2 it retains6 bits instead of 4,851 free bits. Most information
is intentionally unobserved; no H value or node class is asserted preserved.

## C2. Where compression stops

If checking the legality of an arbitrary proposed Hamilton-path re-marking
is an allowed one-step query, THM-4807's reconstruction theorem makes the
response injective on all labelled tournaments. Even one fresh query forces
all 2^binom(n-1,2) fixed-path states to remain distinct. This is a direct
corollary of the inherited theorem, not a new reconstruction claim.

There is also a small explicit hostile even if both re-markings are legal.
At n5 set every free bit to 0 except t_(5,2)=1, and vary t_(3,1)
between 0 and 1. Both sources have the same h1 corner. Both admit the path
(3,2,5,4,1). One re-marking step within that same clock makes their endpoint
chord (3,1) different.
Thus the new carrier cannot silently add arbitrary re-marking to its
operation alphabet. Reading H or blue/black colour also exceeds its scope.

## Connection contract

Source: marked tournament with budget h. Target: binary triangular corner.
Map: (C1). Preserved: every permitted endpoint-chord response. Lost:
interior arcs outside the corner, class/H and arbitrary re-marking legality.
Sidecars: n,h, port order, attachment row's queried h bits, complement phase.
Decisive tests: full endpoint-update replay and the n5 same-clock re-marking hostile.

FRACTRAN's h+2-state observer and this triangular carrier are different
objects. They share a proved operation: delete information that no legal
remaining query can reach, and decrement the query clock honestly.
