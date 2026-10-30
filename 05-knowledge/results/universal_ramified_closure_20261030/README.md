# Universal ramified symmetry: exact closure, selectors and enumeration

2026-10-30. **PROVED scoped / FINITE-EXACT / CITED background / OPEN
individual-bank classification. Same-agent proof review.**

The previous residue-fibre subgroup is now the COMPLETE universal group:
it is exactly the intersection of the automorphism groups of every
independent-sign bank in the tower. This closes a universal symmetry
question without asserting that every individual tournament has that group.

## A complete all-size mechanism

For p>=5, k>=1 and v_p(D)=1, use the literal ramified ring
R=(Z/p^k)[delta]/(delta^2-D), N=p^(2k), with norm-one clock
g=(1+delta)/(1-delta). Let A consist of translations and clock powers.
The full group preserving every directed difference orbit is

    H=A^(2)=F A,
    F={x -> x+eta*h(rho(x)) : h:F_p -> F_p},
    eta=g^(p^(k-1))-1,  |H|=p^(3k+p-2).

The norm first forces any extra symmetry into its last-layer radical.
Directed differences of odd valuation then force that freedom to be
constant across each residue fibre. This proves completeness, not just
the existence of more automorphisms. See [U1](PROOFS.md).

The result includes the prior metallic/Fibonacci family, and also gives
actual tournament banks for p=7,11 and other primes>=5. The200 sqrt(pN)
spectral construction is not imported into the additional sources.
The prime3 is excluded for a reason: a retained mod9 clock has order3
instead of the required9.

## The smallest point-selector base, and information recovered for free

Exactly **p named vertices** are necessary and sufficient to kill H:
one must meet every residue fibre. If a finite pin set meets j>=2
fibres, its stabilizer has order p^(p-j). If it meets only one, the
exact order also depends on the common-ancestor depth in the actual
delta-adic coset tree. [U2](PROOFS.md) covers every finite pin set.

Minimal pins do not imply that direct pair observations resolve every
vertex. At625, five distinct points share their complete pair-color
signature to all five integer pins. A second stage reconstructs the
points j+delta from the existing pins; their observations then resolve
every vertex. **No new chosen selector is needed.** The collision and
the two-stage reconstruction are both retained in [U3](PROOFS.md).
These statements require the fully colored relation, not a single
unlabelled tournament or an arbitrary later complement.

## Exact counts for every number of jointly named vertices

For ordered r-tuples, allowing repetition, Burnside now reduces to
p+k-1 positive exponential terms. At fixed p,k their orbit counts have
a **minimal rational recurrence of order p+k-1**. The proof also gives
an independent count by residue occupancy and common-ancestor depth.
Distinct named vertices use falling factorials instead. See [U4](PROOFS.md).

For the pinned625 tournament, whose full group was already certified
to be H, the first four counts are

    1, 49, 12145, 3778225.

Its recurrence polynomial is

    (X-5) product_(j=1..5)(X-125j).

These count joint vertex markings under the full group of that parent;
they are not counts of tournaments on r vertices or new A000568 terms.
The previous full Q-color attachment recurrence is inherited unchanged.

## Joining the vertex symmetry to the incoming bank recurrence

The incoming Session111 recurrence counts actions ON bank coordinates.
Our H fixes those bank states and is the missing vertex-symmetry kernel.
For the additive affine normalizer M there is an exact extension

    1 -> H -> H M -> M/A -> 1.

Adjoining H leaves bank-state orbits unchanged while multiplying each
known stabilizer order by p^(p-2). At625, the inherited odd action on
all2^24 bank masks has665856,26112,256 orbits with effective stabilizer
orders1,5,25. The corresponding lifted vertex subgroups have orders
5^9,5^10,5^11. These include nonprime and unbalanced banks; full graph
groups may still be larger. [U6](PROOFS.md) states the exact scope.

An ordered packet of at most1+ceil(log2(p^k-1)) actual bank tournaments
can encode every directed color and therefore realize H as its exact
common group. One chosen tournament can be retained as its first layer.
This bit bound is optimal for color-faithful coding, not for the weaker
task of realizing H alone; one625 tournament already does the latter.

## Evidence, sources and remaining work

The checker exhausts the universal two-point stabilizers at25,49 and625;
enumerates the full15625-element group at25; checks all its singleton,
pair and triple pin sets; and tests every anchored pair in six larger
ring universes. It compares the two joint-count formulas and recurrence
in16 parameter families, retains the direct-signature collision and the
prime3 failure, and replays the81/27 control. Normal and optimized runs
are compared. [AUDIT.md](AUDIT.md) records precise universes and hashes.

The standard literature framework is the cyclotomic association scheme
and two-closure; the general nonlinear impurity mechanism is already
in [Evdokimov--Ponomarenko](https://arxiv.org/abs/math/0608469). The exact
quadratic upper bound and selector/counting deductions here have their
own proofs. No literature-priority claim is made.

From the repository root:

    python -X utf8 -B 05-knowledge/results/universal_ramified_closure_20261030/check.py
    python -X utf8 -O -B 05-knowledge/results/universal_ramified_closure_20261030/check.py --output /path/to/separate-receipt.json
    python -X utf8 -B 05-knowledge/results/universal_ramified_closure_20261030/seal.py

The next classification boundary is when a SINGLE bank has extra
symmetries beyond H, and whether its transporters normalize this
universal kernel. General complement closure and unrestricted A000568
reconstruction remain OPEN.
