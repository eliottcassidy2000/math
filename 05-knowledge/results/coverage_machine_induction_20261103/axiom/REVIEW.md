# Independent proof and implementation review

**Status: VERIFIED scoped mathematical audit.** The axiom lane reviewed
the machine, gluing, and junction arguments independently of their authors.
This is an internal agent audit, not a formalization or external review.
The parent separately replayed the five frozen jobs normally and with `-O`.
No broader census was requested or used for this review.

## Machine lane

Reviewed `machines/NOTE.md`, `machines/check.py`, `machines/hankel.py`,
and the scopes of their receipts.

* The twenty-five fractions have disjoint data/control primes and retain
  exactly one control prime. Every source and target control is different
  on an instruction, so fraction reduction cannot erase a control guard.
  The main nonzero test restores its consumed register unit. The division
  two-unit rule precedes the one-unit and fallback rules. Quotient restore,
  prefix scratch, and final constants give precisely
  P'=4P+2b, h'=h+(3-b)P+2b, using the old P on the right. Each loop has a
  decreasing finite register; unread input strictly halves at digit returns.
* The input convention is least-significant-first, with empty input for
  zero. The initialized integer has exponent N at prime2; no claim of a
  compact binary integer encoding of the FRACTRAN input is made. Exact
  marked decoding uses the unique transitive order after deleting the apex.
  The noncanonical same-H9 hostile is correctly identified as a legal
  construction prefix, not a collision between two claimed canonical inputs.
* The six-symbol CA is the radius-one data shift with a stationary marker.
  It is a finite-word reader, not a universal CA. The two-tape machine is
  the specified streaming transducer, with one writing phase per output
  symbol and an explicit final blank test. These machines preserve the
  emitted word; their microstep clocks are not conflated.
* The source always has the fixed Hamilton spine. Native complement is
  consequently legal in both directions. Conjugating the two births by
  that involution preserves the apex name and explicitly replaces the
  terminal port. The native base is strong through overlapping consecutive
  directed triangles; an input1 supplies the missing arc into the apex.
  Empty and all-zero-prefix exceptions are retained.
* The reflection-defect weight is twice the number of nonfinal backward
  apex arcs: 2(2l-popcount(N)). Thus only the empty-input code is blue among
  the canonical inputs, with the inherited literal-line convention.
* The original code's supply bound uses its final1 digit and counts all
  canonical inputs of eligible lengths, including zero. For native exits,
  floor((2+2l)/3)=ceil(2l/3) disjoint base triangles give the lower bound
  3^ceil(2l/3). The resulting exponent
  3 log(2)/(2 log(3)) is strictly below one. The arbitrary-word bound differs
  by at most the stated factor two. Other native paths and gluing are
  explicitly outside these image bounds.
* The 7-by-7 matrix is a genuine concatenation-response minor with anchored
  prefix1. A fixed d-dimensional linear word realization factors it through
  dimension d, so its nonzero determinant proves d>=7. This does not show
  seven suffices, and does not constrain nonlinear or state-dependent
  encodings. The odd-cycle implementation chooses each directed cycle once
  by its least vertex and each packing once by its first unused vertex;
  it is independent of the subset Hamilton DP.

One prose clarification was requested: in the inherited one-apex pair
formula, state “earlier1 at i and later0 at j>i, weight 2^(j-i-1).” The
initial draft wrote 2^(i-j-1) without defining which index was earlier.
The recurrence, code, and receipts already use the correct positive-gap
weight. The parent applied that prose-only correction during integration;
the final note explicitly names i as the earlier position and j as the later.

## Gluing lane

Reviewed G1--G4b of `gluing/NOTE.md` and the corresponding operator,
history, cellular-automaton, Taylor-jet, and parity code.

The component-selection proof retains the unique old cover after deleting
the two wrapper vertices. Selecting the same old component on both sides
is valid only without the joining a->b edge, exactly as the code checks.
At most two distinct anonymous components are chosen, explaining the
second derivative and its positive coefficients. The temporary third-name
claim and the immediate two-port compiler failure are properly scoped.

Uniform-context divisibility uses maximal runs in slot/path order, so its
coefficients count actual exterior patterns. It does not apply to pointed
births which break the homogeneous module. Full coefficient gcd one from
the all-singleton cover is correctly distinguished from continuation
content and fixed slot divisor.

For G3, weighted samples B_r=(r-1)!F(r) turn the wrapper recurrence into a
monic tridiagonal integral update. Hence the first n histories and first n
B samples are related by a unitriangular matrix. The finite-difference
recurrence of the degree-n polynomial, multiplied by (m-1)!, expresses
each later B_m integrally in earlier samples; F(0)=0 handles the initial
boundary. This proves the exact finite gcd formula. Modulo a prime, the
first nonzero sample at r<ell gives a nonzero history by step r-1 because
(r-1)! is a unit. Combining the prime and n bounds gives the stated
min(n-1,ell-2). No unproved equality of prime valuations is used.

For G4b, integer-polynomial Taylor expansion gives
F(7a+j)=V_j+7a D_j modulo49, also for negative a. Differentiating the exact
wrapper law gives the displayed derivative update. Its cyclic indexing is
valid modulo7; the value update instead uses F(7)=7D_0 and
F(-1)=V_6-7D_6 modulo49. The Leibniz product rule is exact. Differentiating
F(-r)=epsilon F(r) gives D_(7-j)=-epsilon D_j; Taylor expansion then gives
V_(7-j)=epsilon(V_j-7D_j). Thus thirteen registers reduce to at most seven
plus the parity sign, without a minimality claim. These registers retain
all modular wrapper responses but do not recover magnitude or source entry.

No correction was requested for these gluing proofs or their code.

## Junction and controller entry

Reviewed J1--J5 and `junction.py`. With a=h-s, d=F-q and delta=d-a,
after k>=1 R steps the rooted profile is
(a+q+(2^k-1)d,q,q,2^k d+q). Applying N then R gives
b-3t=2delta and, since t is odd,
gcd(b,t)=gcd(delta,2^k d+q). At the regular five-source the parameters
(a,q,d,delta)=(12,5,16,4) give the displayed primitive family.

The common-automorphism lemma is exact on aligned named graphs: any common
automorphism preserves their agreement arcs, which are precisely the
directed Hamilton spine, and hence fixes every vertex. It does not say
that the new graph has no new automorphisms.

An earlier draft confused satisfaction of the renewal arithmetic guard
with membership in the old eight-state projective controller. The audit
independently recovered the hostile k=2 ratio10 after C3. The current J4
has the correct bounded entry bridge: after C3 apply C5 unless the ratio
is7; in that exceptional case apply C3 then C5. The exceptional ratio
becomes3, and every successful C5 lands at0. The current receipt covers
all ten input k-classes. This repair has been checked; the initial direct
entry sentence must not be inherited. The parent records the correction
lineage separately.

The subsequent O(X^(2/3)) image bound sums the inherited balanced-family
bound over geometrically growing b_k. Its declared output language and
exclusion of later fresh wrappers/native choices are retained.

## Cross-audit of the axiom lane

The parent reviewed `AXIOM.md`: inherited finite core budget and decidability,
the SCC argument making prime/7p cores strong, conditional native compiler
clock, and the sigma/rho address map passed. The gluing lane independently
reviewed `COUNTER_CARRIER.md`: strongness, opposite-polarity increments,
the full eighteen-polynomial collision, zero-test separation, and the
native arc-transport budget passed. The three-fraction transfer program
and its clock are direct integer identities and have literal controls.

All conclusions remain within the named operation languages. Neither
the coverage axiom, machine simulations, modular closure, nor these
finite certificates establishes the open spectrum classification.
