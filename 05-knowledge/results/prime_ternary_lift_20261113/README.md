# Prime-lane ternary constructions with changing cores

2026-11-13, Codex. **PROVED IN NOTES / FINITE-EXACT**, as separately scoped.
Chan--Pak and ordinary Collatz remain **OPEN**. This is a continuation of
the prime-lane strategy, not a claim of a complete selector or a new
attained-count record.

The shared-source result is [CORE_REPLACEMENTS.md](CORE_REPLACEMENTS.md):
for every m>=1, one strong source with H=15*2^m-1 has explicit strong
children with H'=3H+2 and H'=3H+4. The +4 child retains every old vertex
and every old arc touching the growing transitive module; it reverses
one fixed core arc and adds three vertices. Universal path-cover formulas
prove the all-m identities.

The smaller [mixed-mode +4 move](HIGHER_MODES.md) starts from
H=6*3^m+6*2^m-5. It reverses one core arc and adds **one vertex**, giving
25->79,73->223,205->619 and all later members. Another family realizes
3H+4 using four old core reversals and **no new vertices**. The complete
small-frame search supplies312 full-response transitions, of which84
use three or four visits to the growing module. These extend the
construction interfaces beyond the dyadic-only class.

Two structural results explain why such exits are not automatically
renewable. [FINITE_PHASE.md](FINITE_PHASE.md) rules out a cycle containing
a +4 step in any finite one-module frame grammar with fixed nonnegative
length increments. [CLONE_PARITY.md](CLONE_PARITY.md) proves that, in any
fixed frame with any number of transitive blocks, H modulo4 depends only
on which blocks are singletons. Pure block growth must consume a fresh
singleton for every +4 step; an unbounded supply or a core change is
essential in that grammar.

The [Collatz comparison](COLLATZ_AND_TARGETS.md) has a successful lift
and a precise failure. The mixed-mode +4 construction realizes a genuine
r->3r+1 branch with r=(H+1)/2 odd. Its next halved count family cannot be
realized by a fixed one-module frame, even on an arithmetic progression;
the whole family is also excluded for fixed finitely many steadily
growing transitive blocks by the modulo4 theorem. Individual counts are
not excluded. A residue-dependent or changing-module repair remains OPEN.

[SURPLUS.md](SURPLUS.md) determines the direct triangle-context digits
at all tournament orders for attained h<=73, using a proved padding
lemma and the cited Busch minimum. Dropping indegree one adds47 root
occurrences in the order-eight bank, but cannot repair every direct exit.

The existential prime/7p induction needs one applicable source for each
target, not one finite universally renewable gadget. The remaining gap
is a proved unbounded source selector, or an alternative family entering
at each required smaller parent. Neither the312 transitions nor residue
coverage proves that premise.

[BOARD.md](BOARD.md) retains the inheritance pass, portfolio and live
concepts. [AUDIT.md](AUDIT.md) gives exact universes, independent controls,
dependency pins, reproduction commands and limits. No new theorem ID
was reserved, and no new A000568 or other global sequence term is claimed.
