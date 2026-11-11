# Address arithmetic, a suffix obstruction, and a finite exact carry lift

**PROVED IN THIS NOTE:** a dyadic-cylinder obstruction to uniform ternary
returns inside a bounded fixed-core bank, and an explicit nonnegative
matrix lift for every fixed odd address multiplier. **FINITE-EXACT:** the
bounded controls in `check.json`. **OPEN:** exact target coverage and the
Chan–Pak conjecture. No external novelty or canon-promotion claim.

## Inheritance and the operation being studied

The closest proved mechanism is the arbitrary-cut two-flank compiler in
[the preceding H5 proof](../../chanpak_55_105_golden_20261110/h5/PROOF.md),
also independently obtained in
[Session251](../../session251_creative_completion_20261111/ternary/PROOF.md).
Its seven-coordinate reduction and integral lattice are already in
[Session253](../../session253_cross_problem_frontiers_20261111/integral/PROOF.md).
The incoming [H13 ternary rebase](../../h13_paired_descent_20261111/square_graceful/RESIDUE_PHASE.md)
already supplies the multiplier-three carry system and an actual paired
ray. The present carry formula generalizes that system to arbitrary fixed
odd multipliers; it does not rediscover its m=3 application.

The canonical hostile is the distinction between a scalar target and an
actual source: H105 has disjoint317/319 source families. The corrected
near miss is using congruence reachability as exact positive target
coverage. The least-used sidecar here is the integer address height. Our
board is fixed core / named cut / binary address / parent-count coefficient /
carry / actual graph reconstruction.

Let K be a nonempty named tournament core, O a named cut, and d=H(K).
The actual graph is an apex attachment to TT_a -> K -> TT_b, with binary
flank integers A,B in the inherited convention. Its count is

    Phi(A,B)=d(g(A)+g(B)-2)+qAB+eA+sB+t+u(A mod2)+v(B mod2),   (1)

where every coefficient is supplied by the named core/cut. Here t is the
count on K plus the apex; d is the deletion-core count, NOT the count of
the whole source graph Phi(A,B). The inherited function is

    g(0)=g(1)=1,
    g(2n+1)=g(n),
    g(2n)=g(n)+n+(n mod2).                                (2)

The finite five-core bank has deletion counts1,3,5,9. Its large finite
coverage is inherited and is not rerun here. We ask what an unlimited
*uniform address-rewrite proof* could do inside such a bank.

## C1. A four-point detector for the deletion core

The recurrence gives, for every a>=0,

    g(4a)   =g(a)+3a+(a mod2),
    g(4a+1) =g(a)+a+(a mod2),
    g(4a+2) =g(a)+2a+2,
    g(4a+3) =g(a).

Thus the signed stencil (1,-1,-1,1) gives exactly -2. It annihilates
constants, linear functions of the address, and address parity. Applying
it to (1) at any fixed B proves

    Phi(4a,B)-Phi(4a+1,B)-Phi(4a+2,B)+Phi(4a+3,B)=-2d.       (3)

This is a convenient coordinate-extraction form of the inherited minimal
compiler, not a new path-counting formula. In particular the four values
recover d modulo every odd modulus, prime OR composite, since2 is a unit.

For every fixed k>=0 and 0<=r<2^k, induction in the appended digits gives
integers c_1,c_2,c_3 with

    g(2^k a+r)=g(a)+c_1 a+c_2(a mod2)+c_3.                (4)

The coefficient of g(a) is exactly one. The induction is literal: an
appended1 leaves g unchanged, whereas an appended0 adds the old address
and its parity. Leading zero digits are allowed with their fixed width.
No numerical search or asymptotic approximation enters (4).

Consequently (3) still extracts d after restricting either address to a
binary cylinder, A=2^k a+r, B=2^l b+s. The other terms in (1) remain
bilinear, linear, parity or constant, so the free-a stencil annihilates
them. It also survives further fixed suffix appends and an interchange
of the two flank roles.

## C2. Uniform suffix rules cannot close a bounded core bank

**Theorem.** Suppose a rule replaces an actual coloured core/cut by another
one, appends fixed finite binary suffixes to the two addresses (possibly
interchanging the flanks), and satisfies the identity

    Phi_new(A_new,B_new)=lambda Phi_old(A,B)+kappa          (5)

for every A,B in a full binary cylinder. Then

    H(K_new)=lambda H(K_old).                             (6)

The same conclusion holds if (5) is required only for all sufficiently
large free addresses in that cylinder. The constants lambda,kappa may be
chosen separately for each rule, but cannot depend on the free addresses.

**Proof.** Restrict to the cylinder and apply the translated four-point
stencil at free addresses4a,...,4a+3, with the other free address fixed.
By (3)-(4), the two sides of (5) give -2H(K_new) and
-2lambda H(K_old). Large enough a handles the cofinite version. QED.

**Bounded low-bit policies.** A policy may choose among finitely many
rules using finitely many low bits of A and B. Refine to a cylinder on
which that choice is constant and apply the theorem there. For any fixed
number of successive rules, refine the initial cylinder further to fix
all their choices. This is possible because suffix appends and flank
interchange pull bounded low-bit predicates back to bounded low-bit
predicates. The composed addresses are again suffix appends, possibly
interchanged. If each step realizes3Phi+2 or3Phi+4 uniformly, then after j
steps the deletion-core count is exactly3^j d.

Therefore no finite bank of positive deletion-core counts supports
unlimited such uniform ternary returns on all address cylinders. For
the bank{1,3,5,9}, the only possible count transitions are1->3 and3->9;
from5 or9 a uniform next step would need15 or27. This is an obstruction
to the proposed proof language, NOT to the corresponding numerical
Hamilton counts or to individual one-vertex returns.

The quantifier boundary matters. This theorem does not exclude rules
guarded by unbounded height or by equations involving g, arbitrary
finite automata reading the full address, variable-length suffix choices,
or non-dyadic address maps. The existing H13 repeated-1001 rebase is
guarded by a specific defect equation, and its ternary rebase uses3B.
Neither is contradicted. These are exactly the sorts of operations a
bounded-core construction must add if it hopes to prove global closure.

## C3. An exact carry lift for every odd multiplier

Fix an odd integer m>=1 and retain

    V_m(n)=(G_0(n),...,G_(m-1)(n), n, epsilon_0(n),epsilon_1(n)),
    G_r(n)=g(mn+r), epsilon_i(n)=[n mod2=i].               (7)

This is m+3 integer coordinates, not m+3 finite states. The address height
n is retained, and epsilon_0+epsilon_1=1. The initial vector is
(g(0),...,g(m-1),0,1,0).

For an input digit a in{0,1} and each output residue r, write

    ma+r=2j+b,      0<=j<m, b in{0,1}.

The exact binary update is

    G_r(2n+a)=G_j(n)                              if b=1,
    G_r(2n+a)=G_j(n)+mn+j+epsilon_(1-(j mod2))(n) if b=0.   (8)

The remaining coordinates update by

    n'=2n+a(epsilon_0+epsilon_1),
    epsilon'_a=epsilon_0+epsilon_1, epsilon'_(1-a)=0.        (9)

**Proof.** Since m(2n+a)+r=2(mn+j)+b, apply (2). In the even case
(mn+j) mod2 equals(n+j) mod2 because m is odd; this is the indicated
one-hot parity coordinate. Equations(9) are ordinary binary reading.
QED.

Equations(8)-(9) are two explicit nonnegative integer matrices of order
m+3. Constants j are j(epsilon_0+epsilon_1), so no affine coordinate is
missing. They calculate all g(mn+r) from the binary digits of n in
O(m log(n+1)) elementary coordinate updates. This is a finite-dimensional
exact representation of an infinite address family, with shared height
and parity instead of m duplicate clocks. No minimal-dimension theorem
or bounded bit-complexity claim is asserted.

The construction works unchanged for odd composite m, including9 and15.
Its information-preserving step is the quotient/remainder pair(j,b), not
field division. Dropping b or the integer height is a genuine projection;
whether m is prime does not repair that lost data.

## C4. The lawful Collatz connection acts on addresses

For an odd address x=2n+1, one division by2 after the odd affine step gives

    (3x+1)/2=3n+2,       (3x-1)/2=3n+1.                   (10)

Thus the m=3 lift computes g at both destinations in its channels G_2 and
G_1. This is one halving, not the fully accelerated map that removes every
factor of2. Further acceleration needs its actual valuation/word data.

For any fixed named K,O and fixed left address A, substitute B=mn+r in
(1). Since g(B)=G_r(n), B=mn+r and its parity is epsilon_{1-(r mod2)},
the resulting Hamilton count is a fixed linear readout of V_m(n), with a
constant absorbed by epsilon_0+epsilon_1. The word of mn+r reconstructs
the actual right flank and all apex incidences. Every output therefore
has a real tournament certificate.

The preserved data are the named core/cut and the exact new address and
count. Arbitrary old right-flank names are not automatically retained by
this address rebase. This is not yet a native reset with two old ports,
nor a closed17-coordinate transition. If ports must survive, their
embeddings and full attachment responses remain a separate obligation.

A small hostile prevents confusing (10) with a scalar Hamilton recurrence.
For K=C3, cut O={0}, and A=1, (1) is

    Phi(1,B)=3g(B)+15B+10+(B mod2).

The odd address1 moves to2 under the plus branch in (10), but the actual
counts move29->49, not29->89 or91. Arithmetic on the address supplies a
construction language; it does not by itself prove3H+2/3H+4 returns.

The remaining positive problem can now be posed without an opaque
integer operation: find a legal chart switch and a linear equality among
these carry coordinates that attains the required exact target, and prove
that its zero fibre can always be reached with controlled height. The
inherited m=3 H13 ray supplies one such fibre. The present finite lift
makes all odd multipliers available for this question without claiming
that all those fibres are nonempty.

## Reproduction and controls

`check.py` checks all fixed suffix words through length5 at128 input
addresses; translated stencils; core stencils on four named small cores;
the m=1,3,5,7,9,15 matrices at256 input addresses; leading-zero
invariance; and actual graph certificates for both channels in (10).
Direct subset path counting, independent of the imported arithmetic
compiler, recounts the graph controls. The normal and optimized receipts
agree and contain source/script hashes. The declared finite checks
corroborate the all-order proofs; they do not supply the quantifiers.

    python3 -B 05-knowledge/results/central_modular_chanpak_20261110/compiler/check.py
    python3 -B -O 05-knowledge/results/central_modular_chanpak_20261110/compiler/check.py
