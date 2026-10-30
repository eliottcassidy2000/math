# What the Pisano merger adds to the symmetry program

This session's unifying object is a finite module equipped with **an action,
an observation, and an allowed attachment map**. A symmetry of one component
need not preserve the other two. The useful theorems specify all three.
This is a proposed organizing framework, not a classification of all mathematics.

## 1. Three distinct losses, three exact repairs

1. **Observation loss.** In the Fibonacci algebra, two adjacent traces have
   determinant-five pairing. At ramification, a five-element radical is
   invisible. A marked coefficient digit restores a snapshot. Globally another
   CRT component can already restore the cyclic phase, even though the whole
   algebra still has a trace radical.
2. **Attachment loss.** Five child orbit signs can vary independently above one
   parent sign. A parent summary is exact only for the coherent subfamily, or
   after retaining the four relative signs. Kernel differences add four further
   bits. This is an explicit arithmetic family in which joint port data is an
   actual coordinate, not a metaphor.
3. **Translation-gauge loss.** The cokernel of I−Q^r records the failure to
   remove a carry. The carry after the first linear return determines affine
   period inflation. These two records are not equivalent: a nonzero
   obstruction can have zero return carry.

All are elementary kernel/image/cokernel phenomena, but identifying their
cardinalities does not identify their operations. At stride four the factor
2Q−I and discriminant five give an explicit common algebraic origin. At stride
five the translation obstruction is instead eleven; the trace pairing remains
controlled by five. This is a useful hostile to a single-number explanation.

## 2. An exact recursive tournament family, and its composition boundary

For mu^2−a mu−1=0 and an odd prime p with v_p(a^2+4)=1, multiplication by
mu^4 on R_k=Z[mu]/p^k has an orientation bank of B_k=p^k−1 independent bits.
Hence the number T_k of labelled invariant Cayley tournaments obeys

    B_k=p B_(k−1)+(p−1),      B_0=0,
    T_k=2^(p−1) T_(k−1)^p,   T_0=1.

The formal ordinary generating function for the **bit counts**, not T_k, is

    sum_(k>=0) B_k z^k = (p−1)z / ((1−z)(1−pz)).

The recurrence is supported by an actual attachment map. Each nonzero orbit
below has p orbit children, while the kernel contributes p−1 sign bits.
Choosing a reference child separates one reference sign from p−1 relative
signs. The choice is not canonical; the underlying p-child incidence is.

Coherence says all children take their parent's sign. For a fixed lower
orientation exactly 2^(p−1) choices remain. At p=5, these coherent lifts are
ordinary lexicographic products T_lower[H_kernel], where H_kernel is one of
the 16 invariant Cayley tournaments on 25 vertices. Any section of reduction
provides a product labeling; it need not split the additive group. Cross-fibre
orientation depends only on the reduced difference. Within each fibre,
translation identifies the same kernel tournament. Thus the fibres are modules.

The same proof applies to the general ramified setting with its established
kernel identification and orbit count. Fully coherent towers through level k
have 2^(k(p−1)) labelled choices. Generic top-level banks have 2^(p^k−1).
Their difference is precisely the joint child-sign freedom lost by a parse
that stores only parent tournaments. At level2,p=5 the counts are256 coherent
versus16,777,216 unrestricted. The independent audit checks all390,000 ordered
distinct vertex pairs of a coherent level2 example against product adjacency.

This is a concrete bridge to prior tournament modular composition, not a new
formula for A000568. Different connection sets can define isomorphic tournaments;
counts here retain the fixed additive labeling and specified multiplier.
The free Klein observer quotient identifies only quadratic conjugation and
reversal, not arbitrary tournament isomorphism. Likewise a noncoherent lift
is not thereby proved prime: it only fails this particular module partition.

## 3. Why 3, 5, and 19 behave differently here

For Q the Fibonacci matrix, these prime roles follow from its discriminant5,
not from the earlier divisor-sum encoding:

| p | splitting | pi(p) | preserving power | bits on p^2 vertices |
|---:|---|---:|---|---:|
|3|inert|8|Q^8|4|
|5|ramified|20|Q^4|4|
|11|split|10|Q^2|12|
|19|split|18|Q^2|20|

Every number in the table is independently checked in `synthesis/check.py`.
For3, Q^4=−I reverses all16 connection sets in the identity-action bank.
For5, Q^2 reverses all16 Q^4-invariant connection sets, despite having order10.
For19 and11 no power of Q reverses a tournament: the two-adic part of their
Pisano periods is only2, so the half-order involution has determinant−1,
not the required scalar−I. External negation always reverses a Cayley tournament;
its absence from this particular cyclic action is the relevant distinction.

The earlier sigma-squareclass roles of these same primes remain separate:
sigma(3)=4, sigma(5)=6, sigma(19)=20. No bijection identifying those arithmetic
fibres with these Pisano banks is supplied or claimed. The present connection
is structural: specify the observation and its exact invisible subspace before
calling two encodings equivalent.

## 4. A five-state phase can be read with three serial binary observations

A five-valued snapshot is minimal for restoring a five-state trace fibre from
one additional contemporaneous observation. That is **not** a lower bound of
five on a sensor's alphabet when its time history may be used.

Suppose the full cyclic clock is m=5t and the trace-history clock is t.
In each trace fibre, advancing by t cycles the five hidden phases. Choose a
representative and assign binary outputs in cyclic order 00011. Three readings
at times0,t,2t are respectively

    000, 001, 011, 110, 100.

They distinguish all five phases. Two binary readings cannot distinguish five
possibilities. Therefore three readings are sharp for this fixed-step,
fixed-length binary scheme. `synthesis/check.py` verifies the abstract code
and its explicit realization on all20 powers of Q modulo5 (four trace fibres).

This is a small exact bridge between Pisano clocks, words/necklaces, and the
prior minimum-observer program. It does not get information for free: the
sensor was designed using named phase representatives, and observation time
has increased. It supplies no implementation from a quotient that has already
irreversibly discarded that phase. Snapshot storage, temporal readout, and
the cost of implementing a predicate are different resources.

## 5. The Collatz transfer is affine composition, not a convergence argument

The exact pair (linear part, carry) composes by (P,c)→(BP,Bc+d). This gives a
valid common language for finite Fibonacci transformations and evaluated legal
Collatz words. Reduction modulo an odd modulus preserves endpoint congruence,
but discards magnitude, exact valuation guards and original-source descent.
The hostile27→41→31 remains: a child return below41 does not pay descent below27.
Neither a Pisano return nor a symmetry quotient provides that missing payment.

Likewise the old Zeckendorf two-register reader must retain legality/marker
memory. The current affine Smith classification computes congruence obstruction
modules, not accepted-word counts or a parser for arbitrary hidden marks.

## 6. Next precise questions, not promises or solved conjectures

- Classify unrestricted ramified banks up to full tournament isomorphism;
  the fixed Klein quotient is a strict, computable intermediate quotient.
- Determine which noncoherent banks admit different modular decompositions.
- Extend the discriminant analysis to nonmaximal orders v_p(a^2+4)>1;
  the DVR proof here does not cover them.
- Study temporal observer cost with restricted implementable predicates,
  rather than assuming arbitrary phase tests are free.

No Wall–Sun–Sun nonexistence theorem, LRC closure, Collatz convergence proof,
or unrestricted A000568 enumeration follows from these scoped results.
