# The scale and capacity of a terminal absorber

**PROVED, scoped counting and response obstructions; FINITE-EXACT literal
controls. Chan--Pak remains OPEN.** This packet explains why a small relative
error is not yet a useful approximate prime source. It supplies no new global
Hamiltonian-count absence claim.

## Inheritance and the retained object

The anchor is the actual full count H(T). The nearest positive mechanism is
the shared-source core replacement in
[CORE_REPLACEMENTS](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md):
for any nonempty module M with h=H(M), p=P_M(2), its three frames give

    S(M)=15p-h,       B2(M)=45p-h,       B4(M)=45p+h.             (1)

Thus B2=3S+2h and B4=3S+4h. The arbitrary paired-module mirror in
[MIRROR_MODULE](../../sixty_three_bridge_20261113/MIRROR_MODULE.md) gives
another actual response, 12p^2+8hp-h^2, and a quadratic-norm exit. Its unit
defect becomes h^2 for a nontransitive input. These exact constructions are
inherited, not discoveries here.

The canonical hostile is a primitive response ideal with no desired actual
continuation: AP7 of
[all-prime tournament geometry](../../all_prime_tournament_geometry_20261103/PROOFS.md).
The all-prime uniform-module content and boundary-divisor theorems there,
and G3--G5 of
[the cover-content packet](../../coverage_machine_induction_20261103/gluing/NOTE.md),
are already proved. We do not claim them anew. The new obstructions below
also apply to modules whose cover content is one.

The corrected near miss is to infer target placement from a sublinear
residual, or prime supply from an unlimited number of reordered ordinal
blocks. The least-used sidecar is how many *different scalar values* the
whole architecture can place below a given height, before correction.

Concept board: full path-cover profile / scalar approximation gap / correction
capacity / fixed versus growing frame / commuting versus marked block order.
The anchor asks for exact critical targets, the niche studies cover-profile
commutation, and the wildcard tests whether an apparently rich reorderable
bank has enough distinct scalar outputs to support absorption.

## T1. A fixed frame and atom bank has only polylogarithmic scalar capacity

Let A_1,...,A_d be a fixed finite bank of nontrivial strong tournaments, and
let S be the one-vertex tournament. A module is any nonempty ordinal word in
these atoms and S, in any order and with arbitrary multiplicities. Let F be
a fixed finite bank of tournament frames, each with at most t marked vertices
to be replaced by such modules. Every replacement is uniform: all vertices
of a module inherit the incidences of its marked frame vertex. Other frame
vertices remain fixed singletons. Repeated modules, converse constraints,
and additional restrictions only shrink the output set.

**Theorem.** The number of distinct full Hamiltonian counts at most X supplied
by this entire architecture is

    O((1+log X)^D),       D=t(d+1).                           (2)

Constants may depend on the finite frame and atom banks. This statement
allows cyclic interactions between marked slots and arbitrarily many visits
to them. It is not limited to two-run frames.

**Proof, first mechanism: all cover profiles commute under ordinal sums.**
Let K_j(M) count ordered covers by j nonempty directed paths. Let P_M(q)
count ordered q-slot covers, allowing empty slots. Then

    P_M(q)=sum_j binom(q,j) K_j(M),
    P_(A->B)(q)=P_A(q) P_B(q).                              (3)

For the second identity, each slot concatenates its path in A and its path
in B; this is a bijection, including empty paths. Binomial inversion recovers
every K_j from P_M(q). Hence the complete cover profile of an ordinal word
depends only on atom multiplicities, not their order. This is stronger than
retaining only h or p.

The Hamilton paths of any uniformly inflated frame split uniquely into their
ordered runs through each module. For each fixed tuple of run numbers, their
number is a frame-word coefficient times the product of the corresponding
K_j's. For fixed module orders the sum is finite. Consequently the full H of
the inflated frame depends only on the atom-multiplicity vectors, even when
the number of possible visits grows. No independence of active switches is
being asserted.

**Second mechanism: cyclic module size is logarithmic.** If a marked frame
vertex lies in a nonsingleton SCC of the frame, every vertex of its inflated
module lies in a nonsingleton SCC of the result. Indeed a frame cycle through
that vertex lifts to a cycle through any chosen vertex of the module, using
one representative in the other nonempty modules. Thus its entire module
order is at most

    c(T)<=K(H(T))<=3 log_3 H(T),                            (4)

by the inherited cyclic-core theorem. All d+1 multiplicities at this slot
are O(log X).

If a marked frame vertex is a singleton SCC, its module is an ordinal-sum
factor of the result. Its contribution is just

    H(M)=product_i H(A_i)^(r_i).

The singleton multiplicity is invisible. Since each H(A_i)>=3 and all other
SCC factors are positive, each r_i<=log_3 X. All fixed frame SCCs and slots
therefore supply at most the bound in (2). QED.

The input for (4) is
[the exact singleton-SCC kernel/core budget](../../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md),
also used in [the previous growth packet](../../strategy_fibres_absorption_20261113/GROWTH.md).
No total-order bound is applied to neutral singleton padding.

For a single marked slot there is also an elementary alternative to (4).
If the frame permits j>=2 visits and |M|>=j, then
K_j(M)>=2^(|M|-j): fix one vertex in each of the j ordered bins, and assign
each remaining vertex to either of two selected bins. Every induced bin has
at least one Hamilton path. Thus any active multi-run term forces |M|=O(log H).
If only one visit is permitted, the response is a fixed multiple of H(M).

## T2. Sublinear terminal repair cannot supply the prime lanes from that bank

Let B be the set of base counts in T1. Give an absorber every possible
integer within C b^alpha of every b in B, with fixed C>0 and 0<=alpha<1.
This deliberately grants more correction capability than any actual
tournament construction has been proved to possess.

**Theorem.** At most

    O(X^alpha (1+log X)^D)                                (5)

positive targets up to X can be supplied in this way. In particular, the
proportion of primes up to X that can be supplied tends to zero; the same
holds for targets 7p with p prime. Neither critical lane is covered.

For sufficiently large b, C b^alpha<b/2. A corrected target <=X therefore
has b<=2X, apart from a fixed finite initial set. There are O(log^D X) such
base values and each supplies at most 2C(2X)^alpha+1 integers. This proves
(5).

No unproved prime-distribution conjecture is needed. For completeness,
binom(2n,n)>=4^n/(2n+1), whereas the power of any prime dividing this binomial
is at most 2n. The latter follows by writing its valuation as
sum_a (floor(2n/p^a)-2 floor(n/p^a)); each summand is zero or one. Hence

    pi(2n)>= (n log 4-log(2n+1))/log(2n),

which gives pi(X) >= c X/log X for large X. Dividing (5) by this bound gives
O(X^(alpha-1) log^(D+1) X), tending to zero. Replacing X by X/7 proves the
seven-prime assertion. Removing the finitely many excluded small cases has
no effect. QED.

More generally, if the allowed symmetric radius is a nondecreasing function
R(b)=o(b), coverage of every prime up to X requires

    R(2X) >= c' X/(1+log X)^(D+1)                          (6)

for all sufficiently large X, with a constant depending on the architecture.
This is a capacity lower bound, not an existence theorem for absorbers of
that radius. It concerns the magnitude of the count correction; it does not
bound the number of vertices or arcs used to perform it.

**Exact scope.** A fixed bank may cover infinitely many primes, and nothing
here disproves that. It cannot cover all primes with the stated sublinear
budget. A finite alphabet of recursively applied graph operations is not
automatically a fixed frame bank: its accumulated frame can grow with the
word. Nonuniform incidences to successive atoms also leave T1's hypothesis.
Those distinctions are necessary, not cosmetic loopholes.

## T3. Primitive modules still obstruct bounded fixed-frame corrections

Let Q=(2,4,9,3), the strong four-vertex tournament, and let M_r be its r-fold
ordinal sum, r>=1. Its slot polynomial and first two responses are

    P_Q(q)=q^4+4q^2,
    P_(M_r)(q)=(q^4+4q^2)^r,
    h_r=5^r,                 p_r=32^r.                       (7)

The full ordered-cover content is one: already
gcd(K_1,K_2)=gcd(5^r,32^r-2*5^r)=1. Thus a content divisor cannot explain
the following failure.

**Theorem.** For any fixed finite collection of single-slot tournament
frames and any fixed integer linear comparison of their full counts on M_r,
the comparison is either identically zero for all r>=1, or its absolute
value tends to infinity exponentially. In particular, a fixed finite menu
of such comparisons cannot equal a fixed nonzero bounded repair digit for
arbitrarily large r.

To prove this, write each response as sum_j a_j K_j(M_r), with bounded j
for each fixed frame. Binomial inversion in (3) gives

    K_j(M_r)=sum_(q=1)^j (-1)^(j-q) binom(j,q) lambda_q^r,
    lambda_q=q^4+4q^2.                                    (8)

The lambda_q are distinct strictly increasing integers, all at least five.
Combine coefficients in the finite linear comparison and take its largest
nonzero base. That term eventually dominates all smaller bases. If there
is no such base, the comparison is identically zero. QED.

This works for any fixed nontransitive base module whose H is greater than
one: its P(q) is strictly increasing for positive integers q, and P(1)>1.
The Q example records that gcd one does not remove the obstruction. The
root session's stronger two-run order-floor argument applies beyond ordinal
powers; T3 instead covers arbitrary fixed visit rank on this explicit family.

## T4. A real sublinear debt, and why it still does not place targets

Apply the inherited frames (1) to M_r. They are actual strong tournaments:

    s_r=15*32^r-5^r,
    b2_r=45*32^r-5^r=3s_r+2*5^r,
    b4_r=45*32^r+5^r=3s_r+4*5^r.                        (9)

Their orders are 4r+5, 4r+7, and 4r+8. With

    alpha=log(5)/log(32) < 1,

the residuals are O(s_r^alpha). This is a proved family of approximate ternary
transitions with controlled full counts and primitive cover response. The
ratios between successive source scales nevertheless tend to 32, so these
sources do not approximate every prescribed target. T2 quantifies the gap
even after finitely many such frame/module charts are combined.

There is a stronger positive calibration. For any prescribed epsilon>0,
choose a fixed integer a>=0 so that

    log(5)/((5+a)log(2)) <= epsilon.

Use M_r=(Q->TT_a)^r. Slotwise multiplication gives h=5^r and
p=2^((5+a)r), with the same primitive cover content one. The three literal
frames then give strong tournaments of orders (4+a)r+5, (4+a)r+7 and
(4+a)r+8, and

    B2=3S+2*5^r,       B4=3S+4*5^r,
    2*5^r,4*5^r = O(S^epsilon).                          (9a)

The cyclic cores grow unboundedly, since all three outputs are strong.
Thus even arbitrarily small power-law debt and primitive growing cores
can coexist with the capacity obstruction. This is a proved actual family,
not a claim that its scales approximate arbitrary primes. Varying a and r
still uses the same finite atom bank {Q,S}, so T1--T2 continue to apply.

The equal-two-slot substitution Q^j -> TT_(5j) does preserve p=32^j while
changing h=5^j to 1. More generally, on a fixed two-slot level p=32^r the
literal modules

    M_(r,j)=Q^j -> TT_(5(r-j)),       0<=j<=r,

give h=5^j. Thus (1) supplies the exact cluster

    45*32^r-5^j,  45*32^r+5^j,       0<=j<=r.             (10)

The cluster contains 2(r+1) values. Its largest adjacent gap for r>=1 is
4*5^(r-1). This honest full-count knob is not an interval absorber: equal
p does not make the h values additive or independently selectable. T1
predicts precisely this loss of scalar capacity.

## What should change in the next construction

The follow-up [two-run cone theorem and C5](TWO_RUN_CONE.md) strengthen the
specific Q^r probe: for r>=3, no uniform one-slot frame of any order or visit
rank can repair S(Q^r) to 3S+2 or 3S+4. The third cover response overshoots
the target, forcing two runs, and the all-order cone then excludes the
required coefficient. Thus merely growing the correction exterior does
not escape on this primitive family; the module itself must change.

The bounded-jump terminal-absorption lemma in
[the previous packet, A1](../../strategy_fibres_absorption_20261113/absorption/PROOF.md)
is still valid. The stronger obligation exposed here is to give an approximate
source architecture enough distinct outputs at the target scale before
asking a sublinear residual absorber to finish them.

Three precise escape tests are worth retaining:

1. Permit the external frame or the atom bank to grow with the target, and
   quantify how many distinct *full counts* this adds. Merely permuting an
   ordinal module's fixed atoms adds none under any uniform exterior.
2. Mark atom positions by nonuniform exterior incidences. Test order
   sensitivity on actual H or the full required marked response. Such a
   construction can escape cover-profile commutation; formal word entropy
   alone is insufficient.
3. If staying in a fixed commutative bank, budget an almost-linear correction
   radius as in (6), rather than extrapolating from constant digits or the
   sublinear debt in (9). This remains an unconstructed possibility.

These are requirements on possible proof mechanisms, not new assertions
that one of them succeeds. They clarify why the successful incoming core
replacements and primitive mirror sources do not yet combine into a uniform
prime selector.

## Reproduction

    python -B 05-knowledge/results/prime_terminal_20261113_1610/absorption/check.py
    python -B -O 05-knowledge/results/prime_terminal_20261113_1610/absorption/check.py

The checker verifies the Q slot polynomial, full cover content on three ordinal
powers, and actual frame inflations using Hamilton subset DP independently
of (9). It also checks finite equal-p clusters and the precise spectral
coefficient transform. The universal statements are proved above, not
extrapolated from the controls. No broad tournament census is repeated.
