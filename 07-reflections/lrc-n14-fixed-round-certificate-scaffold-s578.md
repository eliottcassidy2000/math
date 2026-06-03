# LRC n=14 Fixed-Round Certificate Scaffold (S578)

The previous step said: do not label all 190 converse-merged nodes equally; start with the 64 fixed classes.  Oracle then made the next correction: even the 64 are an overcount for tightness.  The tight speed family looks like `{AP,V*}`, while the 64 classes are a boundary surface whose meaning is tangled by the 14-gon tie-wall.

S578 tries to keep both truths.  The 64 fixed classes are still worth scaffolding, but not as final certificates.  They are the place where labels have to be reattached.

Opus S570 landed while this was being closed: all `8191` gcd-1 transversals mod `27` are lonely, AP is the unique floor-tight transversal, and every transversal has an unblocked small pair.  That is the clean Boolean/flip-lattice side.  S578 is the complementary composite side: once the unit spine is fixed, V* is the non-transversal floor cousin, and it still goes through the same cheap pair `(1,13)`.

The scaffold result is neat.  Among `4096` labelled round `d`-vectors, `820` have a dihedral anti-witness, and these collapse to `64` fixed groups.  Of the 64, `63` are strong and one is transitive; `63` have a unique dihedral anti-witness and the regular class has `13`.  So the fixed boundary is not random mush.  It has a rigid almost-unique reflection axis, except at the regular AP-like point.

The speed-level result is the more useful one.  On the canonical n=14 unit spine, through slack `42`, there are `531` full D/U/N quotient covers.  None fall below `1/14`.  Exactly two are floor rows: AP slack `(3,6,9,12)` and V* slack `(3,6,9,24)`.  Both use the HYP-2095 cheap pair `(1,13)` at `1/14`.  The only two full-cover rows with no unblocked small pair are positive-measure controls.

That is a very good sign for the paired-or-anchored proof split.  In both normal forms now, "block all cheap pairs" does not create a hard residual; it creates positive measure or leaves the transversal regime.  The measure-zero rows go straight through the cheap pair.

This also explains why class-only reasoning keeps wobbling.  The unlabelled class preserves round/converse/fixed-boundary structure, but it forgets the pair `(1,13)`, the endpoint owner, the unit-shell owner, and whether the row is AP/V*.  A proof cannot stop at the class table.  It has to say how a fixed class lifts to speed owners, or why any bad lift exposes a cheap pair before it matters.

The next lemma I would try to prove is a bridge:

Every tight-boundary n=14 realization either normalizes to the unit-spine/four-slack scaffold while preserving cheap-pair status, or the failed normalization already gives an unblocked small pair or a positive-measure interval.

That sounds more tractable than "prove every one of 64 classes lonely."  It is not class-by-class heroics; it is a preservation theorem for labels.

