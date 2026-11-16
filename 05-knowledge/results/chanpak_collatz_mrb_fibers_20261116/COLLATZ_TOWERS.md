# Decimal Collatz coincidences, and why they are not a decoder

## A complete fixed-offset theorem

The only positive triples satisfying
`2^u = 10^v + 3·2^w`
are (4,1,1) and (10,3,3). In particular the optional equality with
the Collatz odd rule u=3v+1 introduces no further solutions.

If v≠w, the terms on the right have different 2-adic valuations,
so u=min(v,w); division by that power of two gives 1 as a sum
of two positive integers, impossible. If v=w, divide by 2^v to
get 2^s=5^v+3, with s=u-v>0. For s≤7 monotonicity gives v≤3
and the only values are v=1,3.

If s≥8, then 5^v≡−3 mod 256. Since 5^3+3=128, division by
5^3 gives 5^(v−3)≡129. Repeated squaring gives 5^16≡193,
5^32≡129, 5^64≡1 mod 256, so v≡35 mod 64.
Modulo 641, both 2^8 and 5^8 are 256, and squaring gives
2^16≡5^16≡154 and 2^32≡5^32≡−1. Thus 5^v+3≡−125+3=519.
But 519²≡141 and 519⁴≡10. In consequence
519^64≡10^16≡−1, whereas (2^s)^64≡1. No primality
test is needed for this calculation.

This closes the mathematical case left bounded in the earlier
decimal-rank packet. The modular ideas 256 and 641 are not new:
see https://aopswootblog.wordpress.com/2013/03/09/number-theory-4-using-appropriate-moduli-to-solve-exponential-diophantine-equations/
for the known solution. The last power collapses its residue check
into one test. The equation is in the fixed-gap Pillai family, not an
instance of the gap-one Catalan theorem. The original synchronized
claim alone is easier: 2·4^v=5^v+3 and (5/4)^4>2.

## Functional graph versus bounded coincidence

Encode each positive integer n by y=2^n. The standard Collatz map
is then y↦√y on integer squares and y↦2y³ on nonsquares, restricted
to the positive powers of two other than 1. This is an exact conjugate,
not a new dynamical theorem. The formula is total there and each vertex
has one outgoing edge. The decimal relation picks just two odd edges,
1→4 and 3→10; extending it by congruences would not recover the graph.

For comparison use the accelerated map T(n)=n/2 (even) or
(3n+1)/2 (odd). For a prescribed parity word b_0...b_(k−1),
let s=sum b_i and let
B=sum over j with b_j=1 of 2^j 3^(sum of b_i for i>j).
Its map is (3^s n+B)/2^k. Exactly one residue class of starting
integers modulo 2^k has that word: inductively a legal word of
length j has one class mod 2^j; replacing n by n+2^j changes
T^j(n) by an odd factor, so exactly one lift has the next parity bit.
The same class is precisely the integrality class of the affine formula
since 3^s is coprime to 2^k.

Asking for a prescribed positive endpoint N instead gives an exact
bounded preimage formula
n=(2^k N−B)/3^s.
A legal positive start exists iff 3^s divides 2^k N−B.
Indeed, if divisible the integer start has the prescribed parities by
the same integrality-class proof, now on all integers. Both branches
keep the nonpositive integers nonpositive, so a positive endpoint
forces a positive start. Thus the arithmetic condition *does* give a
complete positive endpoint class, but only because of this barrier.
Equivalently the natural positivity inequality B<2^k N follows
automatically in the congruence class. For the all-odd word,
B=3^k−2^k; positive endpoints are −1 mod 3^k and the least is 3^k−1,
where the source is 2^k−1. The barrier is not formal: for the 3n+5
variant, the one-odd-step prescription with N=1 passes the congruence
but gives source −1, not a physical positive source.

A separate sharp branch obstruction comes from modulus 3. For the
standard accelerated map, a target divisible by 3 has no odd predecessor,
and its even predecessor is again divisible by 3. Thus its whole backward
tree is the doubling ray. Conversely, for odd N coprime to 3 the odd-only
predecessors are (2^a N−1)/3 with a≥1 and 2^a N≡1 mod3. Among the three
choices of a modulo 6 with the right parity, exactly one makes
2^a N≡1 mod9 and therefore makes the predecessor divisible by 3.
Choosing either other residue produces another positive odd target
coprime to 3 and can be repeated indefinitely. This proves unbounded
*backward* branching for those targets, not a statement about forward
termination; in particular it does not prove Collatz.

## Tournament design criterion

A candidate family must emulate not just the affine formula but the
unique integrality-class lift. For each terminal-reader word w, keep
(1) the image congruence class, (2) its first physically realized
endpoint, (3) the condition that every later endpoint in that class is
constructible with a legal continuation. Covering each residue at
some modulus does not prove (3). Likewise, a minimum from an abstract
matrix orbit is not a physical baseline. The promising stronger strategy
is to prove a physical counterpart of the sign barrier: every integral
inverse with an admissible output preserves the source cone and port
states. That would make the absolute bound automatic, as in Collatz,
rather than a second inequality to establish for each target.

A useful different attack is to fix a target N and permit one of two
families whose lower endpoints have opposing trends. When both families
are proved to realize full one-sided progressions in their actual
common lattice, interval gluing is purely arithmetic: if their
lattice intervals overlap or are adjacent, the union has no interior
hole. The hard part is furnishing a physical residue/bound certificate,
not another asymptotic density estimate. The Collatz example shows why
a certificate can be small when the right coordinate is known, but
does not supply this missing tournament lemma or solve either problem.

### An exact test for the barrier

For positive odd coefficients a,c, use odd branch (an+c)/2 and the
same even branch. At every fixed word all integrality-class statements
above remain valid, with 3 replaced by a and each summand in B multiplied
by c. The statement “fixed-word congruence alone certifies a positive
endpoint has a positive source” holds for all words exactly when c≤a.
Sufficiency follows by checking the largest nonpositive integer in
each branch: 0 on the even branch, −1 on the odd. Necessity follows
from −1 mapping to (c−a)/2>0 in a single odd step. This shows the
barrier should be checked on a residue-conditioned boundary, not on
the real cone without the lattice (where the +c term always leaks).

The same proof works in higher dimension: let each legal digit act
by a rational affine injection on a specified lattice class, and let
C be a proposed physical cone inside the state lattice. If on each
digit domain the complement of C maps into the complement of C, any
integral inverse satisfying the digit conditions and ending in C
must remain in C backwards. This is enough for “physical” only if
every state in C on the relevant class really has a tournament
witness and the moves act on the full retained interface. A cone
defined by inequalities alone gives no such lifting theorem.
