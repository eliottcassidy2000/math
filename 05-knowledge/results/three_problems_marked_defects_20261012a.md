# Marked defects across plane colouring, planar Jacobian, and Collatz

2026-10-12. **PROVED scoped results, FINITE-EXACT controls, and OPEN research
hypotheses.** None of the three main problems is settled by this packet.

## The productive common question

A small diagram may repeat exactly while the information needed for a proof
changes outside that diagram. The useful question is: **which actual defect
does one recursive step remove, and what certifies that removal at the
original source?** The four companion packets make that question concrete.

| Lane | New proved mechanism | Remaining obligation |
|---|---|---|
| Hadwiger–Nelson | The fixed Stechkin trihexagonal partition requires eight whole-interior colours, without any period assumption; a specified local split escapes its finite obstruction | A compatible refinement or geometric gadget with the required complete colour support |
| Planar Jacobian | Every nonlinear polynomial cotangent lift has a zero-differential obstruction for all its polynomial observables; quadratic iterates have an exact critical-value budget | A source-changing construction outside this excluded family, retaining every sheet and standard-area condition |
| Collatz | A variable-depth residue75 rule pays a smaller child; a direct infinite ROOT family pays its own positive measurement floor uniformly under refinement | Ground every remaining source, rather than nearby completed sources |
| Quadratic dynamics | Exact return thickness separates -7/4 and -29/16; repetition never improves that thickness for these cycles | A new operation that changes the retained defect, not merely the finite permutation |

The closest inherited mechanisms were complete marked boundary relations,
the separated rational Keller observable ring, and signed-charge receipts.
The hostile examples were a surviving separating colouring, the critical
zero-momentum source, and an uncovered source with only changing rooted
neighbours. The corrected near miss was the same in form, though not in
content: a compressed invariant did not retain the predicate being proved.
The least-used sidecars selected here were within-tile position, local
ramification order, the signed charge at the original source, and return
multiplicity. These are the session's live concept board.

Incoming commits b4eb117d0 and ddadfd166 supplied
[fixed-source shadow charge](collatz_shadow_charge_20261011.md),
[signed charge receipts](collatz_charge_receipts_20261011.md),
[rhombus field barriers](rhombus_overlap_elliptic_20261011.md), and
[eleventh-root port layers](eleven_port_layers_20261011.md).
The source75 obligation comes directly from that integration. The geometric
side avoids repeating a generator family already carrying an all-word
three-colouring. LRC(14) remains the repository's open anchor, but this
session's assigned anchor is the requested three-problem portfolio.

## 1. Plane colouring: change the cell, not just the period

[The Stechkin conflict packet](stechkin_refined_conflict_graph_20261012a.md)
uses the exact hexagons and two triangle orientations of the inherited
mixed-distance construction. Two regions are adjacent when some pair of
their interior points has distance one. Pairwise region conflicts need not
be realized by one simultaneous selection of points; this is a constraint
graph for constant-on-region colourings.

The full infinite region graph has chromatic number eight. The lower
certificate has only14 regions. Three form a universal K3; the remaining
11 need five colours, as shown by an exact independence-capacity argument.
The upper certificate uses four hexagon parity colours and four triangle
parity colours, with an explicit boundary ownership rule. Neither bound
depends on a chosen colouring period.

This closes a broader route than the earlier40-square fixed-period result.
Increasing the period of the unchanged trihexagonal cells cannot produce
six or seven colours. On the positive side, splitting one specified hexagon
in that14-region witness gives a15-piece graph with an explicit seven-colour
witness. Position inside the old cell is therefore a useful new coordinate.
A finite number of splits cannot fix the whole plane, because an untouched
translate of the obstruction remains.

The next cheap global probe also has a decisive answer. Splitting precisely
the even-even hexagons repeats that local repair with period two and hits
every translate of the old14-region witness. Its13-region period quotient
nevertheless has63 edges and contains a K8. Thus the new pattern still
requires eight colours in that period. This quotient clique is another
region-colour constraint, not eight plane points at mutual unit distance.
The local repair was real; the failure is a new compatibility obstruction.

**OPEN geometric search hypothesis.** A recursively refined configuration
can accumulate useful constraints only when it retains either within-cell
position or a complete multiport relation that changes at a strict step.
This is a strategy hypothesis, not a claim that all such refinements yield
six colours. A concrete six-port conditional mechanism is proved in
[marked return thickness, section4](marked_return_thickness_20261012a.md):
there are151 possible separating partitions and52 equal-port partitions.
Preserving one equal partition while authentically eliminating at least one
separating partition per step would finish in at most151 strict steps.
The missing input is an admissible geometric eliminating fragment at every
stage. At two ports, the complete separating support is a single colour
permutation orbit, so declining positive extension counts cannot certify
partial support elimination.

The next decisive geometric object is a compatible periodic or recursive
split, with all new unit-distance conflicts and boundary ownership retained,
or an exact module that kills one joint separating pattern while preserving
a supplied equal pattern. The two approaches should share the same complete
constraint checker, not just compare vertex counts or field degrees.

## 2. Planar Jacobian: a source sheet becomes a singular observation

[The polynomial-observable packet](quadratic_clock_return_20261012a.md)
extends [THM-3567, separated rational Keller maximal observables](../../01-canon/theorems/THM-3567-separated-rational-keller-maximal-observable-nodal-completion-no-go.md)
to repeated critical values. For F_n=(x^2+c)^{circ n}, v_2(c)<0, set

    v_i=f_c^i(0),              Delta_n(X)=product_(i=1)^n(X-v_i).

The 2^n-1 simple critical points have only n distinct critical values.
Inside the fixed target ring, exactly the observables

    H(X,Y)=sum_j a_j(X)Y^j,
    Delta_n(X)^ceil(j/2) divides a_j(X)

pull back polynomially under (F_n(x),p/F_n'(x)). The condition is necessary
and sufficient in every degree. A nonzero momentum-j coefficient therefore
costs at least n ceil(j/2) base degree. These target rings strictly descend;
their intersection over all n is k[X]. Arbitrarily deep finite polynomial
regularity is possible, while a fixed nonzero-momentum polynomial cannot
survive every depth.

At depth one the source generators are particularly transparent:

    P=x^2+c,          A=xp/2,          B=p^2/4,
    A^2=(P-c)B.

They forget the simultaneous sign (x,p)->(-x,-p). The two source sheets
and the node of the observation surface are the same explicit quotient.
All three differentials vanish at(0,0). More generally, for **any** nonlinear
polynomial F in characteristic zero, every polynomial element of
k(F,p/F') has zero differential at each critical source(alpha,0). Hence
every pair of such observables has Jacobian zero there. A different choice
of rational target formulas cannot turn this entire family into a polynomial
Keller pair on the full affine source.

This is a proved exclusion, not a conjecture about all rational symplectic
maps. The next candidate must genuinely change the source mechanism. Its
first tests should be the full coordinate pullback, generic number of sheets,
deleted divisor, and standard area form. A periodic fibre or a polynomial
observer alone cannot pass these tests.

## 3. Collatz: a paid residue rule and its independently grounded family

[The charged terminal packet](collatz_charged_terminal_family_20261012a.md)
targets the next incoming uncovered source75. With the odd map
U(n)=(3n+1)/2^v_2(3n+1), every n=75 mod128 has actual prefix

    (1,2,a),   a=v_2((27n+23)/8)>=4,
    h=(27n+23)/2^(a+3),      4(h+1)<=n+1.

Its odd-valuation charge is1+(a mod2)<=2, paid by the loss of at least two
binary rank bits. A known certificate for the smaller h transfers to n.
This is an all-height guarded rule; it does not assume all its children
converge.

A directly grounded subfamily removes the child premise entirely:

    n_j=(2^(18j+11)-23)/27,
    n_j --(1,2,18j+8)--> 1,                  j>=0.          (1)

It includes75 and has charge exactly one. Its normalized original-source
score is

    Z(n_j)=1/2-7/2^(18j+11) >= 1017/2048.                 (2)

The packet proves a general fixed-head completion consumer: the exact
ternary terminal phase, original-source decoder, positive score, and a
uniform measurement refinement depth remain tied to the same family
member. This is stronger than completing arbitrary nearby integers. It
still does not imply that an arbitrary supplied integer belongs to the
constructed family; the decoder rejects that leap.

The next source123 yields a second paid all-height rule:

    n=123 mod512,   a=v_2((243n+319)/64)>=3,
    n --(1,2,1,2,a)--> h=(243n+319)/2^(6+a),
    2(h+1)<=n+1,    odd charge=2+(a mod2)<=3.

The proper prefix endpoints all exceed n, and the half-band estimate
follows from26n>=1150, valid throughout this residue. Thus it is another
actual first descent with its own child premise.

For the stated policy with15 supplied seeds and inverse depth at most8,
adding the coarse75 rule increases accepted odd sources through8191 from
1018 to1193; adding the123 rule increases this to1278. The first remaining
source is127. These are relative gains
in an intentionally restricted proof grammar, not new numerical Collatz
verification records.

**OPEN quantitative hypothesis: paid dyadic-band exits.** Put
E(n)=floor(log2(n+1)). For every odd n>1, consider the first actual prefix
ending at h with2(h+1)<=n+1. The proposed inequality is

    odd_charge(prefix) <= 3 [E(n)^2-E(h)^2].              (3)

This is a specific hypothesis to attack, not an assumption used in the
proved packets. If every such first exit exists and(3) holds, induction
on the original positive integer gives a strict ROOT certificate and

    o(n) <= 3[E(n)^2-1].                                 (4)

The squared rank telescopes. Existence of the exits is part of the open
input, not licensed by the inequality alone. An independent algorithm or
guard grammar proving both would pay the incoming source-specific measure
floors. This strengthens mere eventual convergence by a quantitative charge
bound, so it must be challenged as such.

> **Update 2026-10-12 (opus-2026-10-12-S24): (3) is REFUTED at coefficient 3.**
> The least counterexample is `n = 56,924,955` (198 odd steps to the first
> half-band exit `h = 25,168,801`, charge 151 against `3(25^2-24^2) = 147`);
> there are 37 counterexamples below `2^33` and nine among the glide-record
> holders to `2.6e18`, all at drop one; coefficient 4 fails at
> `1,008,932,249,296,231` (ratio 415/97). The random model puts the threshold
> coefficient at `c* = 4.64`. See [the band-budget note](collatz_band_budget_hostile_20261012.md)
> and [HYP-9248](../hypotheses/HYP-9248-first-band-exit-charge-threshold-coefficient.md).
> The finite probe through `10^6` below was simply below the first large deviation.

There is a useful stronger consequence and a finite test interface. If a
strict nonroot trajectory spends at most B odd letters, every such letter
multiplies n-1 by at most2 and every other letter by at most3/4. After
B+t steps, with(n-1)2^B(3/4)^t<2, it must have reached1 or exceeded B.
With B=3(E(n)^2-1), that deadline is O(E(n)^2). The paid-band hypothesis
would consequently give a quadratic-logarithmic bound on the number of
odd steps, as well as charge. The incoming valuation invoice bounds the
total halving count by log2(3n+1)+(L-1)log2(16/5), so the full trajectory
certificate also has O(E(n)^2) binary-step length. These are **CONDITIONAL**
consequences, not measured complexity claims.

The exported finite tester can reject a failed invoice without first assuming
that the source converges: stop at a band exit or at charge B+1, using the
deadline above. What remains open is that this total tester accepts every
positive odd source. This is a concrete acceptance obligation with a source
and a finite bound, rather than an unspecified instruction to iterate longer.

The declared finite control through100001 has136 failures of coefficient1,
first at27, and maximum exact ratio26/9 at55->23. A subsequent hostile probe
through1000001 retains that maximum and finds no coefficient3 failure.
Thus coefficient3 survives that finite population while coefficient1 is
refuted. These bounds
retrospectively measured on these completed trajectories do not prove(3).

The next practical obligations are therefore two different tests: develop a
parameter-wide paid rule for the residual127 first-reset-2 branch, and seek
a source-only proof or explicit counterexample to the paid-band mechanism.
The first improves coverage; the second can change the universal strategy.
For127 itself, the exact first-band word(1^6,2,4,3,3) already reaches the
supplied seed29 with charge8. The missing object is a family rule with paid
guards, not an unknown trajectory for127. Its identity127=2*63+1 identifies
the Mersenne recursion without transferring the primitive-prime exception
from63 to a stopping-time theorem.

## 4. The special values select different levels of information

[The marked-return packet](marked_return_thickness_20261012a.md) proves
the exact distinction:

| Object | Small visible recurrence | Retained defect |
|---|---|---|
| Chosen three-cycle at c=-7/4 | Returns after three steps | Identity modulo P^2, but not P^3, at every repeated return |
| Rational three-cycle at c=-29/16 | Returns after three steps | Multiplier(35/8)^m prevents even first-order identity |
| Critical orbit at c=-7/4 | Third numerator repeats7 | Denominator grows as4^(2^(n-1)); the orbit does not close |
| u->2u+1, sixth term63 | No new prime in3^2*7 | The old prime3 acquires the order-six modulus9 |

The exact bridge63/64=-c(c+1)^2 records a unit cancellation. The affine
three-cycle reversal records a bijection of marked points. Neither one
preserves all the neighbourhood, denominator, or source data in this table.
The old no-primitive-divisor example now functions as a warning against
counting only *distinct* defect labels: multiplicity can carry the new
obligation even when the set of labels stops growing.

## 5. A unifying hypothesis with explicit ways to fail

**OPEN research hypothesis: paid defect growth.** A useful adaptive proof
controller should satisfy an alternative at each unresolved marked source:
either produce an authenticated terminal/smaller-child step, or expose an
independent defect whose cost is bounded by a separately justified potential.
Repeating an old visible pattern is progress only if its hidden multiplicity,
source charge, or surviving boundary support changes in the requisite way.

This is a design hypothesis, not an asserted theorem linking the three
problems. Its three concrete realizations have different statuses:

* Colour support is a finite set. Strict elimination plus a protected
  witness gives the proved151-step conditional bound. Realizable eliminators
  are missing.
* Polynomial regularity pays one new critical value per quadratic depth.
  Fixed coefficient degree gives a proved exact stopping budget. This
  successfully excludes a candidate family, rather than proving JC.
* Collatz charge pays the fixed-source floor. The75 family supplies an
  independent constant budget; universal band-exit payment remains open.

The hostile to the common claim is a controller that repeatedly changes
representation without decreasing any complete support or paying a new
source bill. The quadratic cycles make that failure exact at every depth;
63 shows why prime-support novelty alone cannot detect it. The profitable
next operation is therefore chosen by the first unpaid defect, with a
checker that rejects loss of the source or of a surviving witness.

## Validation and handoff

Each companion has an explicit reproduction command, finite universe,
positive and hostile controls, and retained witness or polynomial receipt.
The scoped infinite statements have displayed proofs. Peer audits and the
final exact check counts are recorded with the completed companion outputs.
No scarcity-sensitive theorem IDs were reserved or promoted in this session.

| Packet | Always-active exact checks | Independent review |
|---|---:|---|
| Marked return thickness and residual127 |304|Family-algebra proof/API audit, including final residual control |
| Polynomial-observable return budgets |1017|Tiling proof/API audit |
| Charged terminal families and band invoice |523424|Family-algebra proof/API audit |
| Stechkin region graph and split contexts |7008|Collatz-memory and family-algebra proof/API audits |

Total:531753 checks, with ordinary/optimized/saved output agreement in each
packet. These counts include explicitly bounded controls; they are not a
substitute for the separate infinite proofs. The documentation validator,
whitespace check and all43 local links in the five new notes pass.
