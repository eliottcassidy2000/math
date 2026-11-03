# A uniform equal-count switch with joint ports and F49 phase

**PROVED all t>=1.** The eighteen endpoint-cover polynomials and their
two-polarity birth operators are inherited from
[GRADED_OBSERVER](../source_selection_closure_20261103/GRADED_OBSERVER.md).
The contributions here are their extraction into the two-entry funnel,
the explicit shared family, and the uniform guarded native switch.

## 1 The joint incidence term is a directly extractable coordinate

A cover is an unordered collection of nonempty directed paths partitioning
the vertices. At two named ports record every component's ordered pair of
named endpoints, using -1 for an unnamed endpoint; completely unnamed
components contribute the exponent j. There are eighteen possible named
endpoint skeletons sigma. Write `Z_sigma(z)=sum_j c_(sigma,j) z^j`.
The degree is unbounded with order; eighteen polynomials are not eighteen
finite states. A port internal to a component is deliberately not its
endpoint, but still remains a physical vertex of the retained source.

For a term with c components, ordered-bin observations are obtained by
assigning distinct bins to components. With falling factorial notation,

    P(m) = sum c_sigma,j (m)_c,
    E_a(m) = sum c_sigma,j [some component ends at a] (m-1)_(c-1),
    S_a(m) = sum c_sigma,j [some component starts at a] (m-1)_(c-1).

For J comprising the two ports,

    C_J(m+1) = sum c_sigma,j N_J(sigma) (m-1)_(c-2),   (12)

where N_J counts ordered **distinct** components, the first ending in J
and the second starting in J. A single component joining the two ports
does not count twice. This is the exact joint correction in the inherited
multi-entry funnel formula. Forgetting which endpoints belong to the same
component loses(12); multiplication of endpoint marginals does not repair it.

## 2 Literal compiler with physical slots

The inherited command `(polarity,mask,discard)` adds a new vertex. Forgotten
old vertices beat it iff polarity=1; named old slot i beats it iff mask_i=1.
Discard0 or1 replaces that slot by the newborn; discard2 forgets the newborn.
The positive path-cover insertion operator preserves the complete joint
polynomial vector for every such command.

Start with physical ports `(u,v)` and J subset of them. Birth the funnel
root using polarity1 and the mask of ports outside J, replacing slot0.
Retain slot1 as the old v. For each tail bit epsilon append

    (1, 2+(1-epsilon), 2).

Every old core/tail vertex beats the new tail vertex, and its arc to the
root is exactly epsilon. To make a ready source, append `(0,1,1)`:
the newborn x is beaten only by the funnel root, and slots become `(r,x)`.
For a reset at selected slot j, let k=1-j and append

    (0, 2^j, k), then (0, 2^k, 2),

selecting slot k at the output. These are literal operations on the named
source. The executable compiler checks the full joint polynomials against
direct path-cover enumeration at small orders and separately checks the
scalar Hamilton count by subset DP.

This operation language has no free re-marking of a forgotten old vertex.
For the native programs below the full labelled source and path sidecar
are also retained. The incoming [counter zero-test
hostile](../coverage_machine_induction_20261103/axiom/COUNTER_CARRIER.md)
explains why the forward polynomials alone must not authorize deletion or
arbitrary selected-path complement.

## 3 Two infinite sources for the same integer

Take C3 on `0->1->2->0`. Use the funnels

    lane1: J={0},    epsilon = 1^t 101,
    lane2: J={0,1},  epsilon = 1^t 001,       t>=1,

then add the indegree-one newborn x at the funnel root r and select x.
Call the resulting strong sources A_t and B_t. Both have t+8 vertices.
Their readiness gaps are1 and2 by the parent start counts.

Here is a direct all-size scalar calculation using the preceding packet's
positive tail automaton. For C3, `P(m)=m^3+2m` and each port has
`E(m)=S(m)=m^2-m+1`. For J={0,1}, `C_J(m+1)=2m`; for J={0} it is zero.
Write u=2^t and v=3^t. The tail coefficients are:

|Lane|L(2)|D(1)|B(1)|L(3)|D(2)|B(2)|
|---|---:|---:|---:|---:|---:|---:|
|1|6u|0|3|21v|3|12u|
|2|5u|0|5|17v|7|16u|

These follow by running the fixed three-letter suffix after the repeated
1. In particular the ready-count formula from the funnel gives

    H(A_t)=15(6u)+3*3-3-1=90u+5,
    H(B_t)=18(5u)+5*3-2-6-2=90u+5.                   (13)

Thus the actual reset outputs are `270*2^t+17` and `270*2^t+19`.
This equality is proved from the word recurrence, not inferred from its
first three values. The finite broader search in research.py located
eleven shared affine formulas in its declared bank; (13) is the smallest
one by H at t=1, not an assertion of global minimality.

## 4 A guarded uniform native realization of the source switch

Use construction-order names: the triangle is0,1,2; r=3; the tail is
`v1,...,v_(t+3)`; x is the last newborn. Let w=v_(t+1), the first suffix
vertex. The two sources differ exactly by reversing `1->r` and `w->r`.
Every other arc and physical port agrees.

The inherited [one-auxiliary arc
compiler](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md#5-one-temporary-vertex-compiles-any-old-arc-reversal)
uses two actual complements for a traced arc u->v in a Hamilton path P,
when v is not last: birth z beating only u, first complement off Pz,
then complement off P with z inserted between u,v, and delete z. Every
old arc except uv returns to its original orientation.

For the first reversal choose the following explicit path in A_t:

    P1 = (1,r,x,2,0,v1,...,v_(t+3)).

For the second choose this path in the first reversal's target:

    P2 = (v1,...,w,r,x,0,1,2,v_(t+2),v_(t+3)).

Both words are legal by the displayed source construction, and r is not
their last vertex. The traced compiler therefore gives a uniform switch
`A_t -> B_t` using **two births, four native complements, two deletions**,
with peak order t+9. It is reversible by reversing the legal program.
All old names, both final physical ports, and the endpoint H are retained.
Each temporary birth has one exceptional out-port. The source row relation,
temporary cuts, two path pairs and deletion names are part of the certificate.

This is a stabilized operation: the four complements are not asserted to
form a path in the fixed-order t+8 metagraph. Intermediate H and strongness
are not asserted constant. It provides an actual common source choice for
the two ternary commands, not a return of their outputs to (13).

## 5 Full golden coordinates and a free phase alignment on every third clock

At selected x, extract the inherited golden observer

    (h,X,Y) = (P(1), (S_x(2)-E_x(1))/2,
                    (S_x(2)+E_x(1))/2).

Both sources have E_x(1)=3. The exact formulas are

    A_t: X=147*3^t+132*2^t-15,      Y=X+3,
    B_t: X=238*3^t+136*2^t-54,      Y=X+3.            (14)

For example at t=1 their observers are `(185,690,693)` and
`(185,932,935)`. To derive(14), apply the newborn identity
`S_new(2)=P_parent(2)-S_parent(2)` to the coefficient table in section3.
It uses the actual joint correction C, including its m=2 value4.

Use the inherited field convention `theta^2=4theta+2` over F7 and

    alpha=4(h+2X-Y)+Y theta,       g=h-Y mod7.          (15)

This is an invertible encoding of (h,X,Y) mod7 into `(alpha,g)`. Under q2,
`alpha -> theta alpha`, g is fixed, and Frobenius sends theta to4-theta.
All five Q4 matrices q0,q1,q2,q4,q5 are retained individually in the receipt.

The difference in X and Y between the sources is

    91*3^t + 4*2^t -39 = 4(2^t-1) mod7.

Since2 has order3 modulo7, **the native switch preserves the complete
golden field interface iff 3 divides t**. This is an endpoint equality
proved on the literal source pair, not an assumed phase-neutral reset.
For t=3m, m>=1, both sources have `(h,X,Y)=(4,5,1) mod7` and
`(alpha,g)=(3+theta,3)`. Every Q4 continuation starts in the same phase,
including all q2 powers and their Frobenius transports.

The phase-fixed subfamily lands naturally in the7p lane:

    h=90*8^m+5,
    3h+2 = 7q_m,     q_m=270(8^m-1)/7+41,
    3h+4 = p_m,      p_m=270*8^m+19,
    q_(m+1)=8q_m-17,     p_(m+1)=8p_m-133.            (16)

At m=1, these outputs are2177=7*311 and2179;311 and2179 are prime by
trial division. No infinitude of primes in either exponential family is
claimed. Both sources at this first phase-fixed clock have order11 and
the reset outputs have order13. The recurrences enumerate actual families,
not all primes or all7p values.

## 6 Why the repeated-tail clock is controllable

The six-state tail automaton records the last before-root status
(empty/good/bad) and whether the after-root path has begun. At a fixed bit,
its m-bin transition is

    M_b(m)=(m-1)I + L_b + A_b,

where L changes only the first coordinate and A only the second. Both are
idempotents and commute. Therefore

    (M-(m-1)I)(M-mI)(M-(m+1)I)=0.                    (17)

This proves the repeated-word exponential formulas: at m=1 only0,1,2
occur, and at m=2 only1,2,3. It also gives period dividing6 modulo7 for
any fixed prefix/suffix golden observer once the zero-eigenvalue initial
term has vanished. In (14), coefficients of3^t vanish modulo7, sharpening
the phase clock to3. This is the mechanism behind the alignment; the
integer growth clock is not itself Frobenius.

The next induction obligation is still substantive. The standard reset
outputs have no ready root by READY's theorem. The shared family cannot
be closed by merely shifting its t clock: a ternary branch multiplies its
leading coefficient by3, whereas a clock shift multiplies it by a power
of2. A finite collection of positive-leading-coefficient dyadic families
cannot have a cycle of exact branch identities using only integer clock
shifts, since such a cycle would force `3^L=2^d`. This is a scoped
obstruction, not a restriction on variable cores or more general sources.
