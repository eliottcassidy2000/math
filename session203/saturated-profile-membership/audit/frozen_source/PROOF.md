# Exact bounded-profile membership by saturated forest states

**PROOF CANDIDATE, not independently audited.** This is a corollary of the
independently audited finite-profile-interface theorem. It proves a
terminating exact algorithm, with an explicit but enormous finite state
bound. It does not prove F(8,137), numerical tail coverage, or polynomial
running time in the binary length of a requested count.

## 1. Fixed physical interface and inputs

Let K and X be positive integers. We ask for every actual rooted profile
(B,s) in the rectangle1<=B<=X,1<=s<=K. A single membership query (B0,s0)
uses X=B0,K=s0. The interface theorem supplies W=W(K) registers and finite
nonnegative integer transitions on spanning path-forest endpoint patterns
with at most two bulk-start components. Every actual profile in this
rectangle has a W-register representative, and every legal word constructs
an actual tournament.

Retain the current active-label set A as finite control. Root marking can
be postponed until the final readout: marking does not alter any graph or
forest coordinate. At a readout, choose any r in A. The B,s readouts are
valid even if other ports are still active. Indeed forgetting those other
ports preserves every one- or two-component forest, since it has at most
two bulk starts, and does not change whether its first or last vertex is r.
Thus the scalar readouts agree before and after those forget operations.

The initial control is A empty, with one empty forest. Introduce u only
when u is a free label, choosing its out-neighbor subset among A. Forget u
only when u is active. A final root is chosen only among active vertices;
no root is silently forgotten. Equivalently append a legal root-mark
operation and then forget the other ports after a witness is found.

## 2. Saturation is an exact congruence for bounded outputs

Set c=max(X,K)+1 and define sigma(z)=min(c,z), coordinatewise on vectors.
For nonnegative integers a_i,x_i,

    sigma(sum_i a_i*x_i)=sigma(sum_i a_i*sigma(x_i)).    (1)

If every contributing x_i<c, the sums are unchanged. If some x_i>=c has
positive coefficient, both sides are c. Zero coefficients simply ignore
that coordinate. This also proves commutation with every nonnegative
integer transition and scalar readout.

The exact B readout must be used in its nonnegative form: count every
one-component forest not starting at r, and every one- or two-component
forest having a component ending at r. Its coefficients are0,1,2. The s
readout counts one-component forests starting at r. Using h+L-s as a signed
readout would not justify (1); there is no cancellation argument here.

Maintain only the saturated forest vector. For every legal instruction
word w, induction using (1) shows that its stored vector equals sigma of
its actual forest vector. Consequently a saturated readout equal to a
requested positive value below c is exactly that actual value. A saturated
readout c means only that the true value is at least c.

Two words with the same active control and saturated vector may be merged.
Every legal continuation has the same saturated observation from either,
by (1). Large counts can be discarded by a zero coefficient, but saturation
handles that exactly; a positive coefficient cannot turn an unknown count
above c into a false small output.

## 3. Finite reachable-state algorithm and physical witnesses

Run breadth-first search from the initial state. A state is its active-label
set and saturated forest vector. For each newly reached state:

1. Read B,s at each active root r and record any pair in the requested
   rectangle, together with the instruction-word predecessor trace and r.
2. Apply each legal introduction or forget operation, saturate, and queue
   the result if its complete state has not previously been visited.
3. Stop after the queue is empty. Return the recorded rectangle exactly.

Only reachable states are used. Arbitrary nonnegative or saturated vectors
are never treated as physical. Each first-discovery trace is a legal graph
construction, and appending a root mark at its recorded active r gives a
literal witness. Every finite word has a representative reachable state,
so interface completeness proves that no actual bounded profile is missed.

An implementation stopped by a time, memory, or state-count limit must
return INCOMPLETE for unrecorded pairs. It may return recorded witnesses
as VERIFIED. It may return ABSENT only after the queue has been exhausted.
An early successful query may stop with its witness.

## 4. Explicit size and complexity bounds

There are W+1 endpoint types (the W ports and bulk), hence at most
T=(W+1)^2 component endpoint-pair types. A retained pattern has at most
W+2 components. Counting all multisets, including invalid patterns as an
upper bound, gives

    D <= binomial(T+W+2, W+2).                         (2)

The exact admissible basis is a finite subset of this set. Unlisted
coordinates are zero. There are at most2^W active controls and(c+1)^D
saturated vectors, so breadth-first search has at most

    N=2^W*(c+1)^D                                    (3)

states. There are at most W*2^W introductions and W forget actions per
state. Root readouts add at most W checks. Formula (3) is an upper bound,
not a claim that all these formal states are reachable.

For fixed K, W and D are constants. Sparse introduction maps choose at most
(W+3)^2 pairs of component slots per input pattern; forget maps regroup a
pattern once. Thus each transition can be evaluated in a number of
saturated arithmetic operations bounded in terms of K times D. Every
coordinate uses O(log(c+1)) bits. Straightforward stored-state search uses
O_K(N*D*log(c+1)) bits plus graph-search bookkeeping and
O_K(N*D*poly(log(c+1))) bit operations, with finite precomputation depending
only on K. Conservative dense-matrix evaluation gives the same fixed-K
conclusion with D^2 in place of D.

This is polynomial in the NUMERICAL target bound X when K is fixed, with
an enormous exponent D. It is generally exponential in the binary length
of X, and the constants and degree grow very rapidly with K. For K=26,
W=15 but even the pattern bound (2) is enormous. No practical full search
at that budget is claimed.

The algorithm proves decidability and supplies physical witnesses or
all-order bounded absences. It does not supply one finite verification
that works uniformly for every X, nor does it decide F(8,137) by any stated
finite cutoff. Applying it to larger and larger finite rectangles remains
a sequence of finite questions, not a tail-coverage proof.
