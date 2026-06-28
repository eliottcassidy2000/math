# LRC14 Component-Spine Endpoint Certificate

HYP-3425 made the covering-floor branch local on interval components, HYP-3426
removes branch ambiguity by the mirror involution `u -> 1-u`, HYP-3427 records
exact wall signatures, and HYP-3428 records the two-adic loss ledger.  The next
useful observation, recorded here as HYP-3429, is that the surviving components
appear to be much lower rank than the raw interval family suggests.

In the 150-row exact audit, every row has a survivor window whose active
endpoint labels have rank at most `2`.  Most best windows are E-only: their
endpoints are even-safe walls and the whole small window dodges the odd
two-color core.  Those are not boring; they are exactly the kind of local
finite-ruler certificate that a proof can try to propagate under two-adic
descent.

The mixed cases are the real obstruction shape.  A mixed `E/B0` or `E/B1`
window says an even-safe wall and one odd branch wall form the local escape
gate.  That is a more theorem-shaped object than total survivor measure,
because it names the active constraints and still remembers the branch that
produces `t`.

Two audited rows lack a mixed even/odd spine, but both are easier E-only
exceptions.  That suggests the lemma should not demand mixed walls everywhere.
The cleaner statement is: every primitive covering row has either an E-only
free component or a rank-2 mixed endpoint spine, with owner-current labels
reserved only for genuine exceptions.

This is still not the covering-floor proof.  The next step is to prove the
endpoint-spine lemma without sampling: show that the two-color odd bad core
cannot swallow every even-safe component unless it emits a named finite
exception.
