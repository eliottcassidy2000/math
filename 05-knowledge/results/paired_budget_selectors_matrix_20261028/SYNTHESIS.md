# What closes, what fails under composition, and the next useful objects

**Status:** proved mechanisms are identified below; suggested research targets
are OPEN. This is a connection ledger, not a claim that different major open
problems have become equivalent.

## A finite weighted controller that actually proves a uniform statement

The ROOT proof is an instance of a useful zero-sum pumping lemma. Consider
two synchronized vector-valued words, with a finite choice of front block,
middle event, back block, and optional head. Suppose a period-p addition to
either block has zero total in **each** channel. The objective is a quadratic
boundary term plus positive parts of paired suffix increments, minus a fixed
multiple of the sum of source norms.

Adding a zero-sum block leaves every increment outside that block unchanged.
Its new contribution depends only on the two suffix sums and the phase at
the insertion point. If those data depend on block lengths modulo p, the
objective is affine in the two numbers of complete blocks. Thus a uniform
upper bound is equivalent, for this grammar, to nonpositive pump slopes and
a bound on its finitely many remainder states. A positive slope supplies an
explicit unbounded hostile family. Zero slopes identify possible sharp
asymptotic coefficients.

Here p=3, there are 702 remainder states, and the sharp one-third slope is
zero in an A-mode state. The actual source, cutoff and integer-clock formulas
are kept as a reconstruction sidecar. Pumping is an **induction on charge
words**, not a newly licensed move along a fixed Collatz orbit.

Taking the user's “turning machines” in the computational sense, this is a
finite weighted automaton with a certified nonpositive loop cost. Its success
depends on proving that the arithmetic inputs enter its accepted grammar.
No finite controller for arbitrary Turing-machine computations or all Collatz
depths is inferred. This matches the model discipline of THM-4048
(`01-canon/theorems/THM-4048-rule30-periodicity-balance-and-model-firewalls.md`):
a finite representation and an eventual global statement need an explicit
bridge between their quantifiers.

## Composition exposes the missing coordinate

The new bounds include the uniform one-step endpoint norm 13. Yet at two
absorbing native steps, `(N,u)=(57,1)` gives endpoint norm 39 and outer paired
debt 27 from original source energy 16. The exact ancestry is in
`depth_two.json`. At `(19,1)` the three child charges all point along
`-omega`, so a diagonal energy sum misses their coherence immediately.

The next controller must retain a joint phase or Gram state for the child
fibres. Merely replacing the child energies with original source mass does
not suffice, as the N=57 witness shows. A useful next question is:

> Can two-generation cutoff events be ordered into a finite collection of
> zero-sum block languages while carrying every child cross term and the
> original clock address?

The cheap test is to reproduce both explicit two-step hostiles before
accepting a proposed state quotient. A successful construction would still
need closure under further depth and a summable boundary allowance.

## A concrete representation-theory bridge

The aligned tournament lift is graded by its reversing involution. Its
skew-adjacency operator exchanges the plus/minus spaces, giving the rectangular
factor `K` and integer Gram matrix `LL^T+2rr^T`. A cyclic three-edge signed
core then decomposes into the trivial and Eisenstein characters, producing
the exact determinant alternatives 9 and -7.

Source: aligned core tournaments and a marked root. Target: signed integer
core, Gram matrix, marked kernel. Map: the explicit formulas in
`SIGNED_CORE_AND_GRAM.md`. Preserved: the full named lift when the core and
orientation bit are retained; characteristic polynomial under the Gram
projection. Lost: global root sign under that projection, and potentially
much more under a spectrum-only projection. Test: the transitive triple and
directed triangle have the same Gram value three and opposite root bits.

The connection with ROOT is a shared phase algebra and a shared failure of
diagonal summaries to retain directional data. There is currently no
trajectory map sending Collatz words to these tournament spectra. That
missing map is stated rather than replaced by matching occurrences of 3,7,9.

## The geometric version retains actual contact multiplicity

The five-interior theorem converts a colouring failure into an exact critical
quotient, then recovers the unit-circle centres that quotient edges demand.
For an induced pentagon, four required distinct triple-contact centres
conflict with the proved maximum three. The boundary pairs remain distinct
through the entire argument. A quotient that reused one geometric point for
several pairs would erase precisely the obstruction.

The next possible architecture has six private interiors. One concrete
abstract test is a quotient `K3 join H`, where H is a seven-vertex
four-critical graph, such as a spindle graph. The three universal boundary
pairs would need to cover all six remaining private vertices. This is a
**candidate contact problem**, not a realized unit-distance selector: the
deleted abstract vertex is represented by a pair, so its old rigid embedding
cannot be assumed. A useful first certificate is the exact matching number
of the external-port coverage relation for a proposed realization.

## Incoming work changes the practical enumeration target

Session89 already supplies the aligned-pair generator through core size five,
the explicit self-converse lift, the residual double-coset gluing rule, and
the exact odd-order orientation-anchor criterion. These are incorporated by
reference. Our extension audits the full order-nine image and gives the
signed-core/Gram refinement.

The further incoming packet
`../consecutive_flags_order10_20261028/README.md` constructs all **9,733,056**
order-ten representatives with a separate full-child audit and compact parent
attachment addresses. It supersedes any navigation statement that only
counts, rather than representatives, are available beyond order nine. This
session reads that result and its declared contract; it does not claim to
have performed that order-ten generation itself.

Session89's anchor theorem also explains why designating an oriented core
helps: its odd-order stabilizer admits an equivariant local matrix codec.
One named vertex alone generally leaves an even permutation group and cannot
remove the obstruction. For compositions, retain the concrete residual
stabilizer and relative map, as its double-coset construction requires.

The signed-core filtration suggests a separate OPEN enumeration refinement:
can cycle-cover statistics modulo successive powers of two be updated cheaply
along the existing parent-only generator, and used as exact rejection filters
without replacing the full transporter certificate? The four-probe zero
fibre `{0,30,45,51}` is the first required hostile: determinant mod 64 restores
one bit, but leaves two genuine possibilities. This is a finite filter
question with a clear baseline, not a proposed scalar formula for A000568.
