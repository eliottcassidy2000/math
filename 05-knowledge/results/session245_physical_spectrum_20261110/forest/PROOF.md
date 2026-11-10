# Constant-letter spectral calculus for the inherited apex tape

**Status:** PROVED elementary deductions from an inherited positive transfer;
FINITE-EXACT independent controls and a restricted forest census. All-target
Chan–Pak coverage remains OPEN. No literature-priority claim.

## Inheritance and scope

The positive root-permutation/gap automaton is **inherited**, not a new
compiler: see Session221 `automata/PROOF.md`, which credits Sessions215/216.
Its dimension is k!·2·3^k and its alphabet is the possible neighborhoods of a
chain vertex against k fixed apexes. Its exact finite inverse and target-order
compression are also inherited. Session199's four-block source compiler is a
different, smaller physical language; neither modular surjectivity nor scalar
state compression licenses the required physical inverse.

The two-star formula below is **inherited THM-4748**, file
`01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`,
crossing case p=q=0. The older proof already counts assignments to c+1 runs.
THM-3181, `01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`,
gives polynomial-exponential repeated ordinal-join laws for path-colour jets.
Those are not claimed new here and are not the operation considered below.

The useful deduction isolated here is that each *individual tape letter* is
exactly a sum of commuting idempotents, hence has a square-free universal
annihilator with explicit eigenvalues 0,...,k+1. This gives a short all-context
run interpolation law without polynomial prefactors. The independent-star
family then gives a transparent nonvanishing leading coefficient k!.
A bounded search of the routed material above found no statement of this
specific projection identity. This is not an exhaustive novelty search.

Concept board: ordered forests; fixed apex tape; independent star interleaving;
commuting projections; target inverse. The source-to-target map sends a physical
leaf word to its exact gap-assignment state, preserving H and all continuations.
It loses the written leaf word unless a predecessor is retained. Projecting to H
loses continuation data, as the explicit hostile below shows.

## 1. The inherited automaton, restated exactly

Fix a tournament A on k named roots. All other vertices form a transitive chain,
in tape order. Fix a root permutation p=(p_0,...,p_(k-1)). A Hamilton path has a
unique form: an increasing chain run in gap0, root p_0, run in gap1, ..., root
p_(k-1), run in gap k. Assigning every chain vertex to its gap specifies its
position because each gap run is increasing.

A gap state is 0 if empty, 1 if occupied and its latest vertex cannot enter the
right root, 2 if occupied and it can. The final gap has states {0,2}, since it
has no right root. Each other gap has {0,1,2}. Thus each p has 2·3^k states.
The initial state has all gaps empty.

Reading a chain vertex of signature a, assign it to a gap g. If that gap is
empty and has a left root, the left root must beat this vertex. Otherwise
this assignment is disallowed. Once a gap is occupied, the previous chain
vertex beats the new one because tape order is transitive. Replace its state
by 1 or 2 according to whether the new last vertex beats the right root,
using 2 in the final gap. All other coordinates stay unchanged.

A terminal state is accepted iff every occupied nonfinal gap has state2,
and every empty internal gap has the required root arc p_(g-1)->p_g.
Summing over all p counts each Hamilton path exactly once. This is a direct
bijection, with no OCF assumption and no extrapolation from examples.

## 2. Commuting-projection theorem

Work over the rational vector space with basis the states for one p, adjoining
an implicit zero vector for a disallowed assignment. Let E_(a,g) be the linear
map that assigns the next signature-a vertex to gap g, or gives zero when the
empty-gap left condition fails. Then

    E_(a,g)^2 = E_(a,g),
    E_(a,g) E_(a,h) = E_(a,h) E_(a,g) for g != h.

For idempotence, after a permitted assignment the gap is occupied and has the
fixed signature-a status; a second assignment fixes that status. A forbidden
assignment remains zero. For commutation, the eligibility condition and update
for g read only coordinate g, whereas those for h read only coordinate h.

The full letter transfer is A_a=sum_(g=0)^k E_(a,g), with multiplicity: different
gap assignments can produce the same resulting state and must be added.
Commuting idempotents are simultaneously diagonalizable over Q: split the
space successively into the kernel and image of each projection. On each joint
eigenspace, A_a acts by the sum of k+1 zero/one eigenvalues. Consequently

    product_(j=0)^(k+1) (A_a-jI) = 0.                         (1)

This holds on every root-permutation block and hence on the whole automaton.
Not every listed eigenvalue must occur in a particular observable quotient.
The theorem does not assert that different letter transfers commute.

Define the rational spectral projection

    P_j(A)=product_(0<=i<=k+1, i!=j) (A-iI)/(j-i).

For every m>=0, A^m=sum_j j^m P_j(A), with 0^0=1. For fixed physical prefix u
and suffix v, therefore,

    H(u a^m v)=sum_(j=0)^(k+1) c_j j^m,                    (2)
    c_j=initial * A_u * P_j(A_a) * A_v * terminal.

All coefficients are rational, and all resulting values are nonnegative
integers. Coefficients need not be positive. For m>=1 only bases1,...,k+1
remain, so k+1 exact values determine all positive-length values by a Vandermonde
inverse. The value at m=0 needs its own zero-eigenvalue correction.

For multiple runs, insert these projections in their actual order. This gives
an exact multivariate exponential polynomial, but the mixed projection products
retain that order. Reordering runs is not justified.

For k=2, every fixed-prefix/suffix constant run obeys

    h_(m+4)=6 h_(m+3)-11 h_(m+2)+6 h_(m+1), m>=0,

and the corresponding order-three recurrence starts at m=1. There are no
polynomial-in-m prefactors because the annihilator is square-free.

## 3. Independent-star specialization and nonvanishing top coefficient

Take transitive roots r0->...->r_(k-1), followed by a tape over {-1,0,...,k-1}.
A -1 vertex loses to every root; symbol i beats exactly root i and loses to
the others. Chain vertices still point forward to later chain vertices.
In this chronological order every reversal is a leaf-to-root edge. Thus every
word is a depth-one forest tournament, with arbitrarily interleaved stars.
This is narrower than all depth-one histories, since all roots are initial.

Let F_k(a) have tape (-1)^a followed by 0,1,...,k-1. By (2), for a>=1,

    H(F_k(a))=sum_(j=1)^(k+1) c_(k,j) j^a,
    c_(k,k+1)=k!.                                        (3)

For the coefficient proof, the top spectral projection of the -1 letter is
product_g E_(-1,g): its joint eigenspace requires every gap projection to have
eigenvalue1. Applied to the empty state, this fills every gap, producing bad
last status in each nonfinal gap and good status in the final one. The suffix
has exactly one vertex that can finish the gap before each root. Each of those
k vertices is therefore forced into its matching gap. All gaps were already
occupied, so no first-vertex condition obstructs this unique assignment.
Each of the k! root permutations contributes exactly1 to the top coefficient.
Thus the (k+1)^a term is genuinely present, not merely an upper eigenvalue bound.

Exact interpolation, using the proved law, gives:

- k=1: 1+2^a for a>=1; a=0 has H=1.
- k=2: -1+4·2^a+2·3^a for every a>=0. This is inherited THM-4748.
- k=3: -3+5·2^a+20·3^a+6·4^a for a>=1; a=0 has H=29.
- k=4: 3-16·2^a+102·3^a+120·4^a+24·5^a for every a>=0.

The verifier independently computes the interpolation values and extra values;
(1)-(3), not the extra finite values, justify the universal formulas.
No assertion of prior novelty is attached to these particular families.

## 4. Exact finite evidence and a continuation hostile

`verify.py` contains only stdlib Python. It independently compares the tape
count with ordinary endpoint-subset DP on all words through six leaves for
k=1,2 and four leaves for k=3: 1561 separate words total. It checks 27288
idempotence/commutation instances. Optimized and ordinary receipts agree.

The restricted two-initial-root census exhausts all 3^m words at each
0<=m<=10. At order12 it gives 10283 distinct counts, maximum30061, and first
allowed miss3091. This is a restricted word census, not the full forest
census, not a count of isomorphism classes, and not a global exclusion.

In this depth-one subclass, H=173 already occurs at order8 with tape

    (-1,0,0,0,-1,1),

and outgoing masks

    (198,124,248,241,225,193,128,2).

This is a positive check that the subclass is not limited to the old fixed-root
unit-start source image. It does not establish a new global coverage cutoff.

A scalar continuation hostile is immediate even at the same order: words
(-1,1) and (0,-1) both have H=3, but appending symbols(-1,0,1) respectively gives
menus (3,15,5) and (3,9,9). The whole transfer state, not just H, is needed.

## 5. Honest inverse boundary

The exact spectrum statement remaining open for this restricted family is:
for which odd h does there exist k and a word over {-1,0,...,k-1} whose
root-permutation transfer outputs h? This family is physically valid but does
not exhaust freely rooted forests. Neither the uniform annihilator, the k!
leading term, nor the small census supplies a predecessor for an arbitrary h.

For a fixed prefix/suffix and one variable run, (2) makes testing long runs
cheap. It does not imply that every target is hit, and a fixed finite collection
of fixed-run families is not promoted to a covering argument. The main task's
all-target construction has not been obtained in this lane.

## Reproduction

    python3 -I -B verify.py > REPLAY.json
    cmp RECEIPT.json REPLAY.json
    python3 -I -B -O verify.py > REPLAY_OPT.json
    cmp RECEIPT.json REPLAY_OPT.json

Every check uses explicit exceptions, not Python assertions. Scratch discovery
files are not required to reproduce the proof-facing receipt.
