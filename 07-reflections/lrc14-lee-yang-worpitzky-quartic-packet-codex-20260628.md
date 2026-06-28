# LRC14 Lee-Yang/Worpitzky/Quartic Packet

HYP-3153 is a follow-on to two fresh mainline results rather than a competing
claim.  HYP-3151 says every quotient has to be judged against a target
function; HYP-3152 says the miss-PGF has a Lee-Yang radius coordinate
`q0=q6*R^6`, with the dip living off the circle.  The useful packet is their
intersection:

```text
root radius + root spread
pair-Pascal mass + cap dip
L_y bimodality margin
even biquadratic fold
odd Worpitzky/ear sidecar
factor-through legality
terminal exit or named debt.
```

The exact scout verifies the low-level identities.  For consec `k=8..13`,
`q0=q6*R^6` is exact by Vieta and the numeric root-radius spreads are
`1.1427..1.3629`.  The pair-mass dips are exactly `1081/76440`, `1/4004`,
and then zero.  For `k=8,9,10`, `L_y=q0+q6+q3/10` is below cap with exact
margins `683/29400`, `106901/2102100`, and `69/910`.

The k=8 identity is the hinge:

```text
10q0 + q3 + 10q6 = 10S0 - 10S1 + 10S2 - 9S3 + 6S4.
```

The even part is the biquadratic `u^4-5u^2+4`, which folds through `v=u^2`.
The odd part is the larger term: `-9S3` has magnitude about `3.15` times
`6S4`, matching the Worpitzky/edge-flip sidecar from HYP-3147/HYP-3151.

After rebasing over HYP-3160/S31ai, the even half has a cleaner name:
total covariance / excess co-emptiness `Sigma-kappa_2`.  The landed scans
verify that consec/even-AP co-maximize `Var(N)`, `mu4`, raw bimodality, and
`L_yK8`, but not entropy or excess kurtosis.  The cumulant tower also warns
that consec does not maximize the associator `Sigma-kappa_3`, the exact `1/7`
claim was a near-coincidence, and `anchor=odd-residue` is false.  So the next
lemma should be a pairwise covariance lemma for the even fold, with the odd
Worpitzky/associator channel retained as genuine sidecar debt.

HYP-3161 strengthens this again: the covariance maximum is exhaustive over all
`3432` bounded k=8 clusters, and the binding rows are genuinely
ferromagnetic in the sector-spin sense.  Consec k=8 and k=10 have all 15
pairwise covariances positive; a minimizer has mixed signs.  The next lemma is
therefore not "FKG proves the upper bound" but "the consecutive cluster is the
ferromagnetic ground state/extremal coupling that maximizes the ordered-state
mass."

The main proof lesson is negative but useful: raw scalar `p0`, raw root radius,
or raw cap numerology is not a legal proof object.  The packet has to carry the
off-circle coordinate and the ordered/odd ear sidecar until a real LRC
predicate-preserving map consumes them.

HYP-3199 adds the parallel n=4 caution.  The fixed-path `a,b,c` tournament
table is a cover with a high-multiplicity `S` fiber; the exact chart is the
partial-score `(0,1,1,2)` section with `x,y`.  The compression
`x=a OR c`, `y=b OR c`, `S=c OR (a and b)` is useful, but it forgets whether
`S` came from the canary edge `c` or from the live pair.  So the HYP-3153
packet should carry the n=4 minimality/deletion fields beside the
Lee-Yang/quartic fields, not just the class table.

HYP-3154 adds the geometric bridge I would keep as a sidecar in the next
scout: `w=z+R^2/z` turns the Lee-Yang circle into the real-rooted De Moivre
axis, with the 7th cyclotomic profile as the cap ideal and off-circle `Im(w)`
as the dip defect.  This reframes the missing inequality as a stability
question, but it is still only a carrier until a Grace-Walsh-Szego or Asano
argument preserves the LRC predicate.

Assumption challenged: tournament vertices are proof packets here, not
runners, arcs, isomorphism classes, or score sequences.  The tournament
fingerprint is transitive with one Hamiltonian path led by the HYP-3152
root-curve packet, then HYP-3151 function legality, then the k=8 bimodality
certificate.

Next proof target: bound the off-circle dip/lambda as

```text
even covariance/excess-co-emptiness contribution
+ odd Worpitzky/ear/associator contribution,
```

then attach the resulting certificate to a live HYP-3141/HYP-3142 row without
discarding the ordered sidecar.
