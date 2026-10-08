# Guard coverage and terminal proof defects are different ledgers

**PROVED:** the labelled dyadic partition, finite-cut accounting, least
uncovered-parameter algorithm, and all-alternatives literal proof consumer.
**FINITE-EXACT:** the complete small ledger universe and the authenticated
429-row parameter bank. **CONDITIONAL:** a rule's common-future child must be
grounded before its receipt becomes a ROOT proof. **OPEN:** covering and
grounding every parameter of the target family.

## 1. Inheritance and the missing interface

[paper_reset2_transfers_20261007](paper_reset2_transfers_20261007.md) already
subtracts dyadic cylinders without enumerating a common modulus. Its ledger
counts asymptotic guard mass; generic finite height cuts are separate.
[collatz_marked_completion_20261007c](collatz_marked_completion_20261007c.md)
already authenticates literal ROOT words and common-future transports, traces
a selected proof plan, and exposes holes or dependency cycles.
[paper_local_global_grounding_20261007](paper_local_global_grounding_20261007.md)
already relates finite authenticated network connectivity to grounding.
These principles are inherited, not new claims of this package.

The new instrument retains **every matching rule identifier on each parameter
cell**, incorporates finite cuts pointwise, and passes supplied literal
alternatives to the existing strict first-hit consumer. A rule contained in
another rule's guard adds no guard mass, but can supply a different child or
the only grounded continuation. Thus a coverage antichain cannot by itself
serve as a proof-choice ledger.

The anchor is the complement of the new paid parameter bank; the niche is a
least-significant-bit decision tree with exact finite heads; the wildcard is
the distinction between union coverage and an OR-network of proof premises.
The live concepts are source parameter, guard cell, finite cutoff, alternative
child, ROOT seed, and defect cut. The canonical hostile below uses actual
words `(1)` and `(1,4)` at source 3. The corrected near miss is to identify
zero new guard mass with zero new proof information. The least-used sidecar
is the complete set of matching rule labels and their child obligations.

The connection sends a finite collection of guarded recipes to an exact
partition of its parameter domain. It preserves all rule choices, exact
point membership and residual cells. It does not supply a ROOT premise.
The literal consumer separately preserves the exact integer and both ordered
common-future words. A cell's representative is never substituted for its
other members.

## 2. A complete dyadic cover witness without a huge common modulus

For nonnegative integer parameter t, write

\[
 C(r,b)=\{t\ge0:t\equiv r\pmod{2^b}\},\qquad
 0\le r<2^b,\quad b\ge0.                              \tag{1}
\]

These cells are the same cylinders in the binary integers when discussing
Haar mass. The address reads low bits first. Two cells are either disjoint
or one contains the other; containment is exactly

\[
 C(r,b)\subseteq C(s,a)
 \quad\Longleftrightarrow\quad b\ge a\ \hbox{and}\ r=s\pmod{2^a}. \tag{2}
\]

A rule has a unique name, a cell, a lower cut t>=m, and an explicit terminal
obligation label. Names and labels are provenance, not evidence that a
symbolic rule is true. In the actual adapter, the arithmetic compiler
recomputes every head identity, native phase, and cutoff before it enters
the bank.

First ignore the cuts. Split the root cell C(0,0) only where some declared
guard continues to a deeper node. At each node retain all rule names whose
guards contain it. Stop where no guard has a deeper boundary. The resulting
leaves are disjoint, cover the whole domain, and carry exactly the matching
name set. They form the maximal cells on which that entire set is constant.
Different child labels are not discarded when two guards coincide or nest.

This construction terminates. Every internal node is a proper binary prefix
of at least one guard, so there are at most `sum b_i` internal nodes and
at most

\[
                            1+\sum_i b_i               \tag{3}
\]

leaves. It does not enumerate 2^(max b_i) residues. Every binary split
preserves mass, so the exact Kraft accounting is

\[
 \sum_{\rm all\ leaves}2^{-b}=1,
 \qquad \mu(\hbox{guard union})=
       \sum_{\rm nonempty\ label\ leaves}2^{-b}.         \tag{4}
\]

An empty label set is an exact uncovered cylinder for this bank, rather than
just a sampled failure. A nonempty label set means at least one authenticated
route is available when the labels come from the arithmetic adapter. It says
nothing yet about a child's ROOT proof.

### Finite cuts are part of pointwise coverage

Sort the distinct cuts together with zero, and divide the nonnegative
integers into the resulting half-open intervals, with the last unbounded.
On each interval the active rule set is fixed; apply the preceding partition
to those rules. For a finite interval `[L,U)`, the exact number of members
of a cell is

\[
 \#(C(r,b)\cap[L,U))=
 \left\lfloor\frac{U-1-r}{2^b}\right\rfloor-
 \left\lfloor\frac{L-1-r}{2^b}\right\rfloor.             \tag{5}
\]

Its first member at or above L is `L+(r-L) mod 2^b`. These formulas yield
both exact uncovered counts and the least uncovered parameter, with no
enumeration of the interval. `least_uncovered` returns `None` iff the
declared finite bank covers every nonnegative integer parameter, including
its finite heads.

The last stratum gives the natural density, since all earlier intervals are
finite. Full tail mass is therefore necessary but not sufficient for this
pointwise cover: the single guard `t>=1` has mass one and still misses t=0.
For a finite bank, full pointwise coverage is equivalent to an empty tail
residual and zero residual counts in every earlier stratum. That is an exact
finite test; it does not claim that the routes in a full cover are grounded.

This is the limited analogue of the finite-core discipline in
[factorion_completion_defects_20261007b](factorion_completion_defects_20261007b.md):
an independently proved reduction justifies checking a finite remaining
head. Here the reduction is merely the explicit finite list of guard cuts,
not an all-source Collatz height theorem.

### Two ways that mass can mislead

Repeated or contained guards add no union mass. Their additional mass must
be measured against the **union already present**, not added as if they were
disjoint. Their proof alternatives can still matter, as Section 3 shows.

Also, increasing finite covers of masses tending to one need not cover a
specified ordinary integer. The disjoint cells

\[
                      C(2^k,k+1),\quad k=0,1,\ldots     \tag{6}
\]

cover exactly the positive ordinary integers, leaving t=0. After k=0 through
d-1 their residual is C(0,d) and their total mass is 1-2^(-d). The infinite
union has mass one but retains that exact source defect. Countable tail
arguments for other rule grammars must be proved separately; finite exact
Kraft accounting cannot be extrapolated into pointwise completion.

## 3. Keep alternatives until their proof obligations are settled

At a supplied **literal** positive odd integer n, the inherited consumer
authenticates three types of node:

* `Root(n,w)`: the actual strict first-hit word to 1;
* `Hole(n,label)`: an unsupplied proof at exactly n;
* `Transport(n,h,a,b,z)`: actual strict prefixes from n and h to the same z>1.

The new finite network allows several transports or a ROOT certificate at
the same source. Every child is explicitly listed as a port. Begin with all
supplied `Root` nodes and repeatedly ground a source if **any** one of its
transport children is grounded. This is the least grounded set of this
supplied network. Breadth-first propagation selects a finite proof forest;
the inherited exact suffix splice then exports an actual first-hit ROOT word.
Every common-future receipt is used in both directions: swapping its source
and child also swaps the two actual words. The reverse is a general join,
not a new assertion of source-size payment. This automatically allows a
grounded source to discharge an already observed forward checkpoint.

The complement is the greatest defect cut with no ROOT seed such that
**both directions of every** supplied common future stay inside the cut.
The proof is elementary reachability: a path to a seed exits any root-free
closed cut, while a vertex outside the least grounded set cannot have an
edge into it. A smaller closed root-free cut is also a valid non-grounding
witness for its members. These are claims about the supplied network, not
about mathematical nonconvergence.

Thus full grounding is equivalent to the existence of a choice of one
decreasing dependency-rank edge at each nonroot port, ending in supplied
ROOT certificates. It is not necessary for every unused alternative to
decrease that rank. A cycle can coexist with a valid exit. Conversely an
unanchored cycle supplies no proof merely because its local identities hold.

The consumer rechecks exact types, source identity, ordered valuations,
matching endpoints and the first-hit boundary. `consume_instances` additionally
checks the same parameter and an active rule identifier for each supplied
instance. Its bank chart string is provenance: for a generic caller it does
not prove a source-family formula. The caller must retain that binding; the
literal ROOT conclusion itself is verified at the explicitly supplied n.
Only the arithmetic adapter supplies the particular all-height symbolic
source-family theorem used in Section 4. No label can create a ROOT node.

### A contained rule with zero new mass can supply the only proof

Consider N(t)=8t+3. Word `(1)` is actual for every t>=0, with endpoint
`12t+5`. Word `(1,4)` is actual precisely when t=0 modulo8, with endpoint
`1+9t/4`. Both statements follow directly by checking the two valuations.
Adding the latter rule to the former changes the guard-union mass by zero.

At t=0, however, the broad rule gives `3 -> 5` with an open child 5, whereas
the narrower rule supplies the actual ROOT word `(1,4)` for 3. Pruning the
narrower rule because its guard is contained loses this proof. At other
members of its cell its endpoint is different; the one ROOT observation is
not promoted to a whole-cell ROOT theorem.

A draft used modulo4, which omitted the final oddness bit. The exact hostile
t=4 gives `35 -> 53` followed by valuation5, not4. The repaired t=8 example
is `67 -> 101 -> 19` with word `(1,4)`. This correction does not change the
zero-added-mass comparison or the literal ROOT example at t=0.

The actual reverse implication `5 -> 3`, witnessed by
`U_()(5)=U_(1)(3)=5`, is locally paid since 3<5. Together with `3 -> 5` it
gives an ungrounded cycle. Supplying the authentic ROOT word `(4)` at 5,
or the authentic `(1,4)` at 3, grounds both through an available exit.
The network retains all alternatives, while the export uses a grounded
proof forest. An endpoint identity without its supplied child proof remains
a defect.

## 4. The actual complement bank and its exact residual

The authenticated data in
[collatz_parameter_cover_20261007e](collatz_parameter_cover_20261007e.md)
give 429 retained paid rules at the source family

\[
             N(t)=2^{924745897+2^{32}t}-1,\qquad t\ge0.  \tag{7}
\]

`parameter_bank` calls that package's compiler/auditor for every row. It
recomputes each full carry identity, the exact pulled-back dyadic guard and
the first-hit height cut. Every retained cut is zero. The old `t=0 mod2^47`
rule remains as a separate labelled child obligation even though it is
already covered by the new union. We do not rerun the proposal search here;
the exact retained data and the authenticated identities are the inputs.

The independent labelled partition has **27,254 leaves**, of which **26,825**
are uncovered and 429 are routed. Its maximum depth is 238 bits. Thus the
calculation never expands a table with 2^238 positions. Its routed mass is

\[
 \frac{87577585436686537861655933644679889301594718852088371532346146423787809}
 {441711766194596082395824375185729628956870974218904739530401550323154944}.
                                                               \tag{8}
\]

This is approximately 0.19826862705329032 in the parameter t. The gain beyond
the old rule is exactly (8)-2^(-47), rather than (8) plus that old mass.
No ROOT mass is inferred from this route mass. The finite sample has 431
uncovered parameters among 0,...,1023; this is a different statistic from (8).

The least missing parameter is **t=4**, and its entire residual cell is

\[
                             t=4\pmod{512}.              \tag{9}
\]

The coarsest residual cells are `t=36 mod128` and `t=116 mod128`. These are
reproducible exact targets for further guarded rules, not a claim that their
integers fail to converge. A rule added elsewhere, including a longer rule
inside `t=0 mod2^47`, does not resolve (9). The ledger can decide its exact
incremental coverage without mistaking depth along one branch for coverage
of the complement.

The arithmetic package also develops a separate infinite run grammar. The
fixed finite partition (8) does not include that countable extension; its
own shrinking-tail proof and finite height cuts must be retained when it is
combined with this ledger. The API can account pointwise for any finite
truncation using the same strata and label rules.

### A separately labelled signed extension

`add_signed` additionally authenticates the 136 bounded proposals and the
one separately declared deeper t=6 rule in
[collatz_signed_parameter_fill_20261007e](collatz_signed_parameter_fill_20261007e.md).
Their negative sibling gaps retain the extra terminal-valuation reserve.
All 137 cuts are zero. Combined with the preceding finite bank, the exact
labelled partition has **38,259 leaves**, including **37,693 uncovered cells**;
its largest address has 455 bits. Its routed density is approximately
0.199776235881211, with the full rational printed in the saved output.
This equals the independent signed package's union calculation.

The least missing parameter remains t=4, and **all of t=4 modulo512 still
lies in the residual**. The finite sample now has 294 missing parameters
among 0,...,1023. The point t=6 has gained an actual rule, but its own tiny
guard is not generalized to the rest of that sample or to (9).

The two infinite run grammars and the separate ternary entries are outside
these finite partition counts. Their parent packages prove their own union
and tail statements. This separation prevents counting a refined old cell
as new coverage or confusing a binary residual with a residual after all
other available rule types have been added.

All symbolic child labels remain OPEN here. The script additionally consumes
three actual modest ordinary-source receipts from retained rules, using the
independent strict first-hit checker. Those tests authenticate the interface,
not ROOT proofs for the enormous Mersenne sources in (7).

## 5. What would turn a guard cover into a completion theorem

A full finite guard cover is one obligation. A well-founded child proof
schema is another. For a family of sources, one must prove that every
covered recipe ends either in an actual ROOT word, in an independently
grounded external source, or in another family member with a strictly
decreasing well-founded rank. All finite-head exceptions must be supplied
as well. Induction would then complete that family.

Merely deleting bits produces a smaller integer, but can leave the chosen
parameter family. The current children `M_(924745897+2^32t-D)` generally
have a different exponent residue. They cannot be declared recursively
grounded by the same chart without an additional closure theorem. This is
why the ledger preserves the child formula and why its root-free defect
cuts are useful even when local guard coverage improves.

## 6. Reproduction and explicit finite universe

```text
python -B -X utf8 04-computation/experiments/collatz_guard_cover_defects_20261007e.py
python -B -O -X utf8 04-computation/experiments/collatz_guard_cover_defects_20261007e.py
```

The 21,067 exact controls include all 576 banks with at most three guards
chosen from the 15 cells of depths 0..3, with declared nonzero cuts; independent
pointwise checks at t=0..23; Kraft identities, finite interval counts and
least-missing checks; 64 shrinking-residual hostiles; the actual `3,5` proof
network and contained-rule example; all 429 arithmetic rows, the exact full
partition and 1,024 sample comparisons; the separate 137 signed entries and
all 566 combined labelled routed cells; three literal receipt exports; and
twelve malformed/source-mismatch controls. Normal and optimized executions
agree. The algorithms use integer and rational arithmetic throughout.

The package owns no new orbit-discovery routine, changes no supplied source,
and imports no presumed ROOT certificate. It provides exact next obligations
and an operational test of whether a proposed repair improves coverage,
grounding, both, or neither.
