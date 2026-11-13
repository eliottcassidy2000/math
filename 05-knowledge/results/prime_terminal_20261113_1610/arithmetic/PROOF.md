# Selective ternary proof forests and actual source reselection

**PROVED, scoped arithmetic compiler and source-selection criterion;
FINITE-EXACT literal controls; CONDITIONAL physical interpretation of
unfilled gates. Chan--Pak remains OPEN.** No prime infinitude is claimed for
any displayed exponential family.

## Inheritance and the question that changes

[Session193, prime induction §§2--3](../../session193_triangulated_proof_20261105/prime_induction/PROOF.md)
already proves the exact relative atoms and the prime-only ternary induction
with seeds 1,3,5,23,49,63,67,343. It explicitly permits unrelated witnesses for
different targets and does not require a common ready class. These are
inherited, not results first obtained here. Its §4 already proves that a
fixed finite bank of fixed-complexity clock families, even with arbitrary
ordinal joins, misses infinitely many targets in both critical lanes.

The older
[six-root forest](../../triple_reentry_interfaces_20261103/RETURN_FOREST.md)
already gives bijective ternary addressing for 3h,3h+2,3h+4. The
[Session304 arithmetic tree](../../session304_collatz_peeling_20261113/ARITHMETIC.md)
already proves the mixed-sign Collatz descent and a different conditional
six-seed ternary compiler. None is re-presented here as a new descent theorem.

The new bookkeeping is a finite *factorization forest*: digit gates are used
only for critical atoms, multiplication constructs their smaller source
counts, and the source used by the physical gate may be reselected. We prove
a total gate budget, an exact adequacy criterion for a declared gate language,
and a cost-optimal divisor recursion. A literal H1919 comparison makes source
reselection concrete, using the incoming
[core replacements](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md).

Inheritance pass: closest proved mechanism, Session193's prime induction;
hostile, the same-H1919 reducible/strong source pair below; corrected near
miss, the incoming C3 module changes digits 2,4 into 6,12; least-used sidecar,
the identity of the physical source selected at a terminal gate. The concept
board is relative atoms / factorization choice / one-shot source certificate /
literal gate interface / division versus reconstruction.

## 1. An explicit conditional proof forest

Let M be the positive odd integers excluding 7 and21. Use the inherited atoms

    A = {odd primes except7}
        union {7p: p prime, p>=5, p!=7}
        union {49,63,343},

and the fixed seed bank B={1,3,5,23,49,63,67,343}. All these seeds have actual
tournament witnesses, independently recounted by the companion.

For each n in M form a formal expression recursively:

1. If n is a seed, stop at its seed leaf.
2. If n is not an atom, factor it into atoms and recurse on the factors.
3. Otherwise set d=2 when n=2 mod3, and d=4 when n=1 mod3. Recurse on
   h=(n-d)/3, then attach a gate labelled “produce 3h+d from an applicable
   actual source of count h”.

This is the constructive form of the inherited induction, not proof that
the labelled gate is physically available. The inverse source is a smaller
allowed odd integer. The boundary h=1 gives5 or7; h=7 gives23 or25; h=21 gives65
or67. Only5,23,67 are unexcluded atoms there, and all are seeded. The carry
atoms49,63,343 are seeded separately. Every remaining gate has h>=3 and h in M.

Every digit target is a prime or seven-prime atom. A source may be composite,
and distinct branches of its factorization forest may use different physical
construction families. A product node is an actual ordinal join once its
children have actual witnesses. A gate remains conditional until a suitable
physical source and operation are certified.

## 2. A total gate budget, not merely a depth bound

**Theorem A.** If a forest for n has g digit gates and its seed leaves have
values b1,...,br, counted with multiplicity, then

    3^g product_i b_i <= n,                            (1)

strictly if g>0. In particular g<=floor(log_3 n), including all branches of
the forest, not just its longest path. With the stated seed convention,
if there are l nonunit leaves then the stronger combined bound is

    g+l <= floor(log_3 n).                             (2)

**Proof.** At a leaf the assertion is equality. A product node multiplies
the inequalities of its children and adds their gate counts. A digit gate
has value3h+d>3h, so it contributes one further factor3 and makes the
inequality strict. Every nonunit seed is at least3, proving the leaf bound.
Equivalently, factorization preserves the product of all pending target
counts; replacing one digit target by its source divides that product by
more than3. QED.

This measures *symbolic proof obligations*. It is not a runtime estimate
for integer factorization, gate search or path counting. It is not a bound
on tournament order: a source selector may replace the current witness by a
much larger graph. An order bound requires an additional cost theorem for
the actual source selector and physical gates. Nor does it bound temporary
Hamiltonian counts occurring inside an unexamined physical gate.

For example the forest of5761 is

    5761 = 3*(19*101)+4,
    19   = 3*5+4,
    101  = 3*(3*11)+2,
    11   = 3*3+2.

It has four gates, depth three and leaf product5*3*3=45. The check
3^4*45=3645<5761 illustrates the total-cost budget. These formulas alone
do not supply the four physical operations.

## 3. Exact adequacy of a declared one-shot gate language

For d=2,4 let D_d be a set of allowed source counts. Membership is certified
by an actual source tournament, its required marks, an allowed finite repair
or construction, and an actual output of count3h+d. The source certificate
must retain the operation's real interface, not just its scalar h. Only
critical outputs in A minus the carry atoms are used as gates.

Let C(D) be the least scalar language generated by the seed bank, products,
and such certified gates whose source count has already been generated.

**Theorem B (exact source-selector criterion, this language only).**

    C(D)=M

if and only if every unseeded atom a has its unique inverse source
h=(a-d)/3 in D_d.

**Proof.** For sufficiency use the strictly smaller induction above: source
h is generated first, and the stipulated gate then supplies a. Nonatoms
follow by products. For necessity, an unseeded atom cannot first enter the
language through a product of nonunit allowed values. Identity factors can
be discarded from its derivation. Therefore its first entry is a digit gate.
There is only one odd inverse source and one digit d in{2,4}, so h must
belong to D_d. QED.

This iff concerns the declared gate/product language. Global Chan--Pak
completeness does not imply that a particular proposed gate family is
adequate. Direct construction families can be added as authenticated leaves;
then their atoms are removed from the residual gate obligation.

No common renewable ready class is required. At a gate, the current product
witness T0 for h and the applicable source T1 may be different. The proof
must either give T1 independently, or give an actual count-preserving repair
from T0 to T1; equality H(T0)=H(T1) by itself certifies neither a marked
interface nor an allowed repair.

There is a useful precise meaning of *uniform* here. If valid finite gate
certificates are recursively enumerable with a decidable checker, a proved
existence assertion for every residual source yields a uniform terminating
selector: enumerate certificates until the source count and digit match.
This gives no uniform size or time bound and no fixed local move. It is
different from assuming that a count match magically transports marks.

## 4. Factorization choice is part of the selector

The inherited atom decomposition may choose any eligible prime to pair with
one remaining factor7. A greedy choice can waste gates or require an
unavailable subtarget even when another factorization is usable.

For the fixed seed bank above, let c(n) be the minimum number of formal
digit gates among these proof forests. It obeys the exact finite recursion

    c(b)=0                                      for seeds b,
    c(a)=1+c((a-d)/3)                          for unseeded atoms a,
    c(n)=min_{n=uv, u,v in M, u,v>1}(c(u)+c(v)) for nonatoms n.       (2)

All arguments on the right are smaller than n, so this is a terminating
divisor dynamic program. Necessity follows from the top operation of a
forest; sufficiency follows by joining minimizing child forests. This is
an optimality theorem for the stated formal language and cost, not for all
actual tournament constructions.

**Exact greedy hostile.** The prime1367 has source455=7*5*13. Pairing7 with
the smallest prime gives35*13, whose cost is2+1=3 before the final gate.
Choosing91*5 instead gives2+0=2. Thus

    greedy cost(1367)=4,       optimal cost(1367)=3.

These are the only relative-atom factorizations of455, so this small optimum
does not depend on an unbounded search. The incoming m=1 +4 frame gives an
actual strong H91 witness, and joining it to the H5 seed gives an actual
H455 source. The final H455->H1367 gate is not supplied here. This example
shows why the source factorization remains selectable instead of hard-coded.

## 5. A literal seven-prime step needing a different source witness

The incoming core-replacement theorem supplies, for every m>=1,

    S[TT_m]: H=15*2^m-1,
    B2[TT_m]: H=3H(S[TT_m])+2,
    B4[TT_m]: H=3H(S[TT_m])+4.                         (3)

The full source frame is S=(2,20,9,35,45,7), module slot0. The +4 frame is
W(10010110), module slot4. The incoming proof certifies the all-m identities
by module response vectors; our finite controls independently construct and
count the literal inflated tournaments.

At m=7 the actual step is

    1919=19*101  ->  5761=7*823.                      (4)

The cofactor823 is prime. Two source witnesses are available:

* a canonical ordinal join of an H19 witness and an H101 witness, of order15;
* the applicable strong source S[TT7], of order12.

Both have H1919. They cannot be isomorphic: the first has two nontrivial
strong components and the second is strong. Thus the product proof of the
source count does not itself put that supplied witness into the incoming
frame. The explicit source family nevertheless supplies an alternative
witness, and its +4 operation produces an actual strong order15 H5761 target.
This is a one-shot source-selection success, not a return theorem.

The literal H101 factor comes from the older
[five-core near-square construction](../../chanpak_five_core_bridges_20261110/history/NEAR_SQUARE_AND_ADJACENT.md):

    (250,508,505,240,224,192,128,256,121).

The H19 factor is W(10100). The receipt retains both factors, their join,
the eligible strong source and the target. Every count is independently
checked by subset path DP. Neutral births alone cannot turn the two-component
join into a strong tournament of the same H: a newborn that joins the old
cyclic components is cyclic and hence raises H by the inherited deletion
criterion. Reassembly or old-arc changes are a genuinely different resource.

There is an infinite arithmetic slice of this source-selection example.
For m=7+18k, k>=0, one has

    19 divides 15*2^m-1,       7 divides 45*2^m+1.       (5)

Indeed2^18=1 modulo both19 and7, while2^7=14 modulo19 and2 modulo7.
Write h_m=19q_m and y_m=7p_m. The actual strong source and actual +4 target
exist for every k by (3). Under the smaller-count induction hypothesis,
q_m>=101 is attained, so a reducible ordinal source19*q_m is also available.
Whenever p_m is prime, this is a critical seven-prime gate with explicit
alternative-source selection. No infinitude of prime p_m is asserted.

## 6. A retained module is not a unit digit

For an arbitrary nonempty module M, the incoming frame identities are

    H(B2[M])=3H(S[M])+2H(M),
    H(B4[M])=3H(S[M])+4H(M).                           (6)

With M=C3, the actual counts are177,537,543. The digits are6 and12,
not2 and4. This physical hostile is independently recounted. Replacing
the transitive module by a recursively generated source therefore requires
checking its H coordinate; the universal module formula does not preserve
unit digits under arbitrary substitution.

This is why a formal forest cannot be evaluated by blindly nesting the same
frame. One-shot source reselection is permitted by Theorem B, but each new
applicable source still needs its own literal or universal certificate.

## 7. How division by three differs from Collatz and from rebuilding

The two inverse parent formulas divide n-2 or n-4 by3. They are integral
because the atom residue determines the digit. Ordinary Collatz instead
multiplies an odd integer by3 and adds1, then divides by powers of2; the
mixed-sign result of Session304 permits both signs. Neither gives any
tournament operation or a source certificate for the forest above.

An actual graph operation dividing H by3 is a separate assertion. The
incoming
[activation/localization analysis](../../sixty_three_bridge_20261113/semigroup/PROOF.md)
already gives the decisive hostile63/3=21:63 is attained,21 is not. It also
constructs abstract multiplicative countermodels in which all permissible
3-divisible values are present but infinitely many prime atoms remain absent.
Consequently all-count source supply at3p does not justify cancellation.

Within strong induction, one can instead *reconstruct* an attained multiple
3h as a triangle joined to a separately constructed smaller witness h. That
does not peel a triangle from the original graph. In our arithmetic forest,
63 is an indivisible relative atom and a seed; it is never factored as3*21.
Among allowed multiples of3,63 is the unique one whose arithmetic quotient
by3 is one of the two forbidden values. This guard is necessary but does
not prove physical divisibility for the other values.

## 8. Exact controls and remaining obligation

Run from the repository root:

    python -B 05-knowledge/results/prime_terminal_20261113_1610/arithmetic/check.py
    python -O -B 05-knowledge/results/prime_terminal_20261113_1610/arithmetic/check.py

The finite arithmetic universe is all49999 allowed odd targets through100001.
Their forests, relative-atom products, total gate budgets and optimal divisor
costs are checked exactly. The greedy factorization is improved on3453
targets; the largest saving in this universe is five gates at34895. The
largest greedy gate count and depth are nine, first at96787. These statistics
describe the fixed formal cost model, not physical construction complexity.

Physical controls recount all eight inherited seed graphs; use independent
permutation counts on the five nontrivial seeds of order at most seven;
recount the incoming source and both children at m=1,2,3,7; test the C3 module
hostile; and verify the actual H455, both H1919 sources, and H5761. Source
files are pinned by LF-normalized SHA256. All checks use explicit runtime
exceptions, including under Python optimization.

The remaining uniform theorem is exactly the source obligation in Theorem B:
for every unseeded critical target, provide an applicable source at its
unique smaller count and a correct one-shot gate. The forest removes the
unnecessary renewal demand and bounds how many such obligations an individual
target uses. It does not discharge their universal physical availability.
