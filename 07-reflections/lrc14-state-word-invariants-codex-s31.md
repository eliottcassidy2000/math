# LRC14 State-Word Invariants

The recurring warning from the last several sessions is that every scalar is
arriving too late.  Relation rank, fold count, positive/negative mass, and even
the `p0,p1` plateau are useful only after the proof has already chosen what to
forget.  The more faithful object is the measured cyclic word of missed-sector
states on the wall arrangement.

That word is small enough to compute exactly and large enough to explain the
recent coincidences.  HYP-2647's balanced endpoint sector transport is not an
accident: the AP-to-defect map is balanced on a coarse address but not after the
state-pair valuation.  HYP-2643's fold move `3/8 -> 3/9` is another address on
the same word.  HYP-2646's reciprocal kernel says the dual version of this
lesson: retain finite mod-7 address before taking signed/absolute sums.

The S31 table made the hierarchy visible.  For AP9 to `(0..7,9)`, the common
state transport has neutral measure `4393/5880`, support `76` state-pairs, and
yet the signed residue is exactly the known `-887/158760`.  So the proof should
not start with the signed number.  It should start with the state-pair support
and prove that the dangerous rows have only a few possible couplings.

There is also a useful negative lesson.  AP9, the near-AP defect, and the
doubled-top spectrum row can share the visible fold count while differing in
state support, entropy, single-sector bias, and transport behavior.  Fold count
alone is only a label on the state word.  The label matters, but it is not the
object.

My current guess: the LRC14 structure is determined by a finite addressed
automaton.  Its base alphabet is the missed-sector subset.  Its addresses are
cyclic transition type, fold target, moving endpoint sector, and mod-7 coimage
phase.  Every known successful invariant is a valuation of that automaton; every
stuck route has tried to value too early.

The next attack should look for template rigidity.  If a row has high `L_y`, its
state word should be low-complexity in a very specific sense: few missed-set
states, few high-jump transitions, and sector-singleton bias concentrated in a
near-AP pattern.  Failure of that rigidity should either increase entropy enough
to drop below the cap or expose a Freiman/GAP/signed-tail certificate already in
the backlog.
