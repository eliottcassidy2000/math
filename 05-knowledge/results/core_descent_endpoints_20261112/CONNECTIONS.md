# Exact bridges for12,13,24,127 and the retained orientation

**PROVED algebraic maps and finite controls; inherited statements credited.
No prime infinitude, Monster action on tournaments, or Chan–Pak proof.**

## 1. The supplied formula, its cube, and its marked derivative

Let U=(127+5sqrt(645))/2. Its conjugate is U^-1, since
127^2-25*645=4. It is the least positive unit greater than1 in the ring
of integers of Q(sqrt(645)): norm-1 is impossible modulo3, and direct
checking y=1,2,3,4 in x^2-645y^2=4 leaves y=5,x=127 as the first choice.

Write a_k=binom(2k,k)binom(3k,k)binom(6k,3k) and

    b_k=sum_(j=0)^k a_j binom(k+j,k-j)(-432)^(k-j).

For A(w)=sum a_k w^k and B(z)=sum b_k z^k, formal summation proves

    B(z)=A(z/(1+432z)^2)/(1+432z).                     (1)

Thus the half-density sqrt(z)B(z) equals sqrt(w)A(w) with coherent
square roots, where w=z/(1+432z)^2. The derivative transports by

    Theta_z=(1-432z)/(1+432z) Theta_w.

Set X=1/z=-432 U^3. The two exact calculations are

    U^3+U^-3=127^3-3*127=2048002,
    1/w=X+864+432^2/X=-960^3,
    320sqrt(645)*(U^3+1)/(U^3-1)=8127.                 (2)

Applying (1) and its derivative to the classical discriminant43 series
gives precisely

    1/pi=24i sum_(k>=0) b_k
       [-3669+320sqrt(645)(k+1/2)]/(-432 U^3)^(k+1/2).

Indeed the transformed numerator is
`-3669+8127(k+1/2)=(3/2)(5418k+263)`, giving the usual coefficient36i
over (-960^3)^(k+1/2). Take sqrt(X)=i sqrt(-X), and sqrt(-960^3)=
i sqrt(960^3); the square roots of the reciprocals are their reciprocals,
not independently chosen principal branches.

The classical input and companion recurrence are **CITED**, from
[H.H.Chan–S.Cooper, Theorem2.1 and Table3](https://www.researchgate.net/publication/232016510_Rational_analogues_of_Ramanujan%27s_series_for_1p),
[publisher DOI](https://doi.org/10.1017/S0305004112000254).
The substitution algebra above is explicit. This is a different Chan from
Swee Hong Chan in the Chan–Pak paper.

With this conventional definition of U, the user's denominator is missing
the cube. The uncubed version converges to45.45380698224419..., whereas
the cubed version agrees with1/pi at the120-digit working precision.
The companion coefficients are independently checked by their recurrence.
The estimate |b_k|<=(4/3)6912^k gives a convergent geometric tail bound;
the receipt separates this analytic bound from numerical roundoff checking.
If U had instead been defined to mean the cube of the fundamental unit,
the notation would change, rather than the formula's mathematical value.

During this session, concurrent
[Session290](../session290_core_descent_20261113/moonshine/PROOF.md)
independently supplied this same pullback and derivative transport. We
integrated it and credit that convergence; this is not claimed as a new
Ramanujan–Sato identity. The label1B here denotes a companion-series
presentation; it does not create a second identity class of the Monster.

## 2. A literal127–13–7 bridge to a tournament

This bridge is an elementary calculation in the specified quadratic ring.
Reduce U modulo13. Its minimal polynomial becomes

    x^2-127x+1 = x^2+3x+1 over F_13.

Its discriminant5 is not a square modulo13, so the quotient is F_169.
Let u be the residue class of x. In coordinates a+b*u its powers0..7 are

    (1,0),(0,1),(12,10),(3,8),(5,5),(8,3),(10,12),(1,0).

Therefore u has order7. Frobenius is conjugation: u^13=u^-1. The
norm-one subgroup of F_169^* has order14; this particular unit lands
in its subgroup of order7, rather than generating the whole group.

Now choose u as generator, and take vertices u^a, a in Z/7Z. Orient
u^a -> u^b exactly when b-a belongs to{1,2,4}, the nonzero squares
modulo7. This is a tournament because exactly one of d,-d is a square.
It is the Paley tournament P7. Its symmetries are the maps

    a -> r*a+t, r in{1,2,4}, t in Z/7Z,

forming7:3. A complete permutation check gives21 automorphisms,
21 orientation-reversing automorphisms, and189 Hamiltonian paths.
With the generator/orientation gauge held fixed, Frobenius sends a to-a
and reverses every arc. Adjoining these reversing maps gives7:6.

**Map and loss.** The map is reduction of the integral unit, followed by
its cyclic subgroup and the written Cayley orientation. It preserves the
order7 multiplication action. It loses the infinite unit exponent modulo7,
the real embedding, and the ambient field addition when only the subgroup
is retained. A chosen generator and orientation gauge are necessary. If
the branch U versus U^-1 is forgotten, the two orientations are identified.
Transporting the gauge along with conjugation would instead give an
isomorphism of the transported structures; the anti-automorphism statement
is specifically about the fixed gauge.

**General family mechanism.** Let p be an odd prime, and let a norm-one
element of F_(p^2) have prime order ell congruent3 modulo4. Its cyclic
subgroup admits the same Paley orientation through a chosen generator.
Frobenius inverts the norm-one group, hence reverses this orientation.
The affine subgroup of order ell(ell-1)/2 acts by automorphisms. Its action
on directed Hamilton paths is free, since an automorphism fixing an ordered
spanning path fixes every vertex. Therefore ell(ell-1)/2 divides H.
This includes the189=21*9 check above. No infinitude of primes giving a
specified order for this particular integer unit is asserted.

The symmetry therefore creates a coverage obstruction as well as a
construction: for ell>=7 these Paley sources cannot have prime H. To reach
arbitrary odd prime counts one needs less symmetric sources. This is
consistent with the inherited
[vertex-transitive prime-count boundary](../fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md),
and explains why retaining local endpoint structure while changing global
symmetry is a more useful goal than maximizing symmetry everywhere.

This makes the orientation loss in the reciprocal invariant
`X+864+432^2/X` concrete. It does **not** turn the analytic1/pi identity
into a count identity or a compiler for an arbitrary odd c. The useful
structural parallel is to E2's forgotten selected gate: equality after
forgetting a coordinate need not support the same next operation.

The older [Monster7-normalizer discussion](../cyclic_attachment_lattices_20261023/MONSTER.md)
retains a shared square/nonsquare sign in its7:6 factor. The explicit
cyclic action above exhibits that same abstract sign mechanism. No embedding
of this entire field construction into the Monster is proved here. In
particular a nonabelian simple group of even order, such as the Monster or
Baby Monster, cannot act nontrivially by automorphisms on a finite tournament:
every tournament automorphism group has odd order. An involution would swap
some pair and reverse its edge; Cauchy's theorem gives the odd-order claim.

## 3. The different roles of12,13 and24 in prior work

The [older quadratic trace carrier](../cursor_cm_phase_20261028/TRACE_CARRIER.md)
used v^2=5v+1 modulo13, a norm-minus-one unit. It has order12. For
x=a+bv, one trace is `2a+5b` and has13 fibres of size13. Adding
`Tr(vx)=5a+27b` recovers all169 states: the two-coordinate determinant
is29, nonzero modulo13. This is a different clock from U645, whose
mod13 order is7. The distinction is a testable norm/discriminant difference.
More precisely, the old discriminant29 is a square modulo13, so its algebra
splits as F_13 x F_13, with roots11,7 and unit order dividing12. The new
discriminant5 modulo13 is nonsquare, and its norm-one torus has order14.
The orders12 and7 are thus consequences of different split types and of
the selected unit, not interchangeable occurrences of the same clock.

The requested set, forgetting its ordering, is exactly

    {4,2,5,3,7,13}={d+1:d divides12}.

Removing composite4 gives{2,3,5,7,13}; their product2730 is the denominator
of B_12=-691/2730. The script verifies this from the defining Bernoulli
recurrence. The repository's
[Borwein/divisor atlas](../borwein_ladder_owner_numbers_vsc_twelve_20261028.md)
already distinguishes this fact from the independent divisor-height cycle
`4->2->5->3->7->4`, with `13->2`, on the specified Monster/Baby-Monster
valuation core. The order of the labels is extra data; the set identity
does not prove that cycle or a tournament recursion.

In the present seed transition,24=(13-5)*3 is the constant required by
transporting the normalized anchor3. In the modular calculation,24 comes
from transport of the derivative coefficient. In the older moonshine
[eta-product calculation](../cursor_cm_phase_20261028/RAMANUJAN.md),24
occurs in a product exponent and a Fricke normalization. These are precise
mechanisms with a common lesson about transporting the observable. They
are not an identification of their underlying objects.

## 4. First25 integer nodes and why endpoints matter

The square-sum graph here has vertices1,...,n and undirected edges x+y
a square. It is not the graph on the first25 primes. The old
[small-graph theorem THM-4505](../../../01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md)
and the recovered
[Session211 port proof](../session257_permutohedron_squares_20261111/history/sources/session211__square-sum-port-bridge__PROOF.md)
already establish the pertinent endpoint mechanism.

At12 the graph has path components of sizes3,4,5. New vertex13 attaches
to3 and12; vertex3 is internal to its old path. Component sizes alone
therefore do not certify a legal concatenation. At25 the new vertex has
neighbors11,24, and all20 directed Hamilton paths use it internally.
Deleting it exposes ten unoriented spanning two-path covers of Q24 with
these prescribed endpoints in different components. Q24 itself has no
Hamilton path. The replay retains all birth endpoint/internal counts
through25; the unoriented totals at15,16,17,23,25 are1,1,1,3,10.

The map is delete/insert the degree-two vertex; its preserved predicate is
the exact port-cover condition. It loses that condition if only the scalar
number of paths or component sizes are kept. This is the useful analogy
for a seed extension. The tempting map H=2p(Q_n)+1 is impossible: Q23,Q25
would give7,21. Neither these square-sum graphs nor their symmetric
digraphs are tournaments.

## 5. Enumeration coordinates and the three requested sequences

The new rational generating function in D2 counts **guarded construction
states**. It makes this infinite family manageable without constructing
every graph. The endpoint quotient has36 classes on the78 original
presentations, but E2 shows why it cannot replace the marked seed census.
Unmarked isomorphism, selected-gate equivalence and ordinary continuation
equivalence are different quotients and should have separate sequences.

The requested OEIS objects also have different universes:
[A000568](https://oeis.org/A000568) counts unlabelled tournaments;
[A135244](https://oeis.org/A135244) maximizes an aliquot predecessor;
[A002496](https://oeis.org/A002496) lists primes of the form k^2+1.
No term extension or faster general enumeration for these three sequences
is claimed here. A cyclic unit subgroup, a formal-state generating function
and a marked endpoint quotient are useful additional recorded aspects,
but none is a substitute for the relevant equivalence or primality test.

## 6. Scope of the conjecture

The target is the spectrum statement in
[Chan–Pak, section6.5, Conjecture6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf):
all positive odd Hamilton-path counts except7,21. That paper credits the
conjecture to bof and Gordon Royle. The new five-base grammar supplies
physical decreasing moves on a sparse domain. Establishing a legal
predecessor or a lower-head repair for every remaining target is still OPEN.
