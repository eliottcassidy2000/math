# Center deficit, exact padding, and the meaning of D(z)=z^4

**INHERITED PROVED:** `D(z)=z^4`, from the concurrently published
[complete center classification](../deficit_two_centers_20261020/README.md).
**PROVED:** its divisor bookkeeping extension below. **FINITE-EXACT:**
independent center-deficit controls through order nine. General black-core
connectivity and a complete positive-deficit center alphabet remain OPEN.

## 1. Definition and inheritance

For an ordinary tournament type let h(T)=H(T)/|Aut(T)|, the number of
Hamiltonian-path orbits. For a self-converse type let b(T) be its blue
path-orbit count. Use **center deficit** e(T)=h(T)-b(T). This is different
from both feedback deficiency `2k+1-heavy_count` and symmetry slack
`3f(T)-moved_support`.

The inherited mechanism is the center factorization in
[THM-2454, pure-blue alphabet and center law](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md),
Section 2. Its older open pure-blue classification is superseded by the
[uniform prefix classification](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md).
The unique ordered SCC decomposition is the required sidecar. The hostile
is R5 -> R5^op: each half is pure-blue, but its orbit count is three,
giving h=9,b=3,e=6. A nonrigid pure-blue half is not harmless padding.

The concurrent [paired-black leaf reduction](../paired_black_leaf_reduction_20261020/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md)
already derives the deficit-two reduction and the sharp degree-three
nonrigid-half boundary. Those results are inherited here, not claimed as
new discoveries. The general divisor expression retains the same mechanism
at arbitrary positive deficit.

The blue involution acts on the h path orbits. Its fixed orbits are counted
by b, so e is even. The inherited odd Hamiltonian-path count and odd
automorphism order imply h is odd. These are path-orbit quantities;
dividing labelled path counts without their automorphism action would
change the proposed identity.

## 2. Exact divisor restriction

Every self-converse type has the unique SCC form

    T = A -> C -> A^op,

where A is an arbitrary ordered SCC word and C is empty or one strong
self-converse type. An empty center has h(C)=b(C)=1. Put a=h(A),
c=h(C), d=b(C). The inherited factorization gives

    h(T)=a^2 c,     b(T)=a d,
    e(T)=a(a c-d).                                      (1)

For q=e(T)>0, a is an odd divisor of q. If a>1, then d<=c implies

    q >= a(a-1)c >= a(a-1).                             (2)

Consequently only finitely many half-orbit counts a can occur at each
fixed q, and every nonrigid half pays at least six units of deficit.
For a>1 the center mass is also bounded:

    c <= floor(q/[a(a-1)]),     d=a c-q/a.               (3)

These restrictions are exact necessary conditions, not realizability
claims for every integer solution.

In particular, at q=2 or q=4 the half must have a=1. The proved rigid
classification then says A is a word over singleton and directed-triangle
SCCs. All the deficit is carried by the single center. At q=6 there are
exactly two algebraic alternatives:

* a=1 and center deficit six;
* a=3 and c=d=1, so the center is empty, a singleton, or C3.

The latter alternatives follow from (2)-(3) and the rigid strong alphabet.
The class of ordinary halves with h=3 is retained as a class series;
it is not silently replaced by one presumed strong atom.

## 3. A small exact recurrence, and its remaining input

Let H_a(z) count ordinary tournament types with h=a, and let R_q(z) count
strong self-converse types of center deficit q. For positive q let F_q(z)
count all self-converse types of deficit q. The unique SCC decomposition
and (1) give

    F_2(z) = R_2(z)/(1-z^2-z^6),
    F_4(z) = R_4(z)/(1-z^2-z^6),                         (4)
    F_6(z) = R_6(z)/(1-z^2-z^6)
             + H_3(z^2)(1+z+z^3).                      (5)

In (4), R_4 is a generating series indexed by deficit four; it is not
the name of the four-vertex tournament.

Set D(z)=R_2(z)=sum d_n z^n and F_2(z)=sum f_n z^n, with negative
indices zero. The exact recurrence is

    f_n = f_(n-2)+f_(n-6)+d_n.                          (6)

The concurrently proved `D(z)=z^4` gives unconditionally

    F_2(z)=z^4/(1-z^2-z^6),

and classifies all self-converse fundamental-black-degree-one nodes as
rigid palindromic padding of the strong four-vertex tournament. The
all-order classification is supplied by the separate center proof, not by
re-expanding (6) or extrapolating our independent order-nine controls.

More generally, for a>1 define C_(a,q)(z) to count empty/strong centers
with `(c,d)` satisfying (3). Then

    F_q(z)=H_1(z^2)R_q(z)
              + sum_(a>1 odd, a|q, a(a-1)<=q) H_a(z^2) C_(a,q)(z).

This is an exact divisor decomposition for the blue deficit layers. It
does not compute unrestricted A000568 without the indicated strong-type
data, and it does not assert black connectivity from degree counts.

## 4. Exact controls and incoming proof audit

Run from the repository root:

    python 05-knowledge/results/deficiency_three_20261020/centers_check.py
    python -O 05-knowledge/results/deficiency_three_20261020/centers_check.py

The checker reads every class in the authenticated order-three through
order-nine censuses, filters self-converse types, and computes strong
connectivity directly from each adjacency matrix. This includes 191,536
ordinary classes and 2,752 self-converse classes at order nine. Every
deficit-two survivor is rechecked by literal Hamiltonian paths and a
separate automorphism counter. The only strong survivor is order four,
with (h,b)=(5,3). Every other deficit-two survivor is its predicted padding.
The independent padding tests cover 89 constructed cases and the
nonrigid-half hostile (9,3).

The low positive strong-center deficits at orders 4,...,9 begin

    2, 4, 6, 6, 12, 20.

This is a finite census statistic, not a monotonicity theorem. The incoming
proof supplies the missing all-order argument: three half-supported
four-cycle-prefix reroutings force deficit at least four; this bounds
vertex connectivity by four. The cycle-cut rotation quotient then forces
an effective order cutoff at 28, leaving seven orders above nine. Every
one of the 2,121,984 necessary symmetry-normalized presentations has an
explicit prefix obstruction. The independent audit reconstructs the
entire affine presentation spaces and checks all certificates. It
correctly leaves the q=2 spaces free of an unjustified reflection condition.

The root read the cutoff, hand configuration supplement, complete census
proof, and independent full-census audit. Our local checker supplies
separate small-order and padding controls; it does not claim to rerun
all 2,121,984 certificates. The final status is taken from the incoming
packet's audited README, which explicitly supersedes its frozen
manuscripts' historical pending-audit lines.

An abandoned draft inference here claimed that a strong degree-one
vertex's descending neighbor must be nonstrong. That is invalid: it
transfers the degree-one hypothesis to the neighbor. A descending neighbor
could itself have another descending edge. The false inference was
removed before this packet's publication. Descent alone does not prove
the center classification or the connectivity of all nonleaf states.

Connection ledger: the map is unique SCC decomposition, preserving h,
b, and center deficit via (1); forgetting the half mass destroys (2), as
R5->R5^op shows. The preserved scalar does not retain a neighbor or path;
endpoint-changing detours still require their explicit Hamiltonian words.
