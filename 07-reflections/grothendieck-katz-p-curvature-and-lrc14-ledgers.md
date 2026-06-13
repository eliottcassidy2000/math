# Grothendieck-Katz p-curvature and the LRC14 ledger

The useful import is not the word "curvature."  It is the discipline of not
confusing a scalar mod-p shadow with the operator that produced it.

Grothendieck-Katz says, in the form this repo needs: a differential equation
whose reductions have vanishing p-curvature for almost all primes should have
finite monodromy.  The local tests are not just values modulo p.  They remember
how the derivative interacts with Frobenius.  That is the little hinge.

The LRC14 analogue is now fairly concrete.  A denominator `q` being blocked is
only a scalar shadow.  The live object is the marked blocker ledger:

```text
which unit twists are blocked,
which runners block each twist,
which Pisano/shell-27 class is being spent,
which 13-clock or divisor fiber is being used,
which Bprime or owner-private deletion route is opened.
```

This is exactly where HYP-2443, HYP-2444, and HYP-2445 now meet.

HYP-2443 says: do not count blocked denominators nakedly; mark support.
HYP-2444 says: in the one-stranger family, the missing support coordinate is
the shell-27 antipodal/Pisano class plus the 13-clock.  HYP-2445 says: in a
serious geometric problem, the scalar finite-field condition can pass while a
retained Frobenius support channel still obstructs the global object.

The p-curvature toy atlas makes the warning sharp.  On the local p-jet
`F_p[z]/(z^p)`, the rational connection `a=1/(1-z)` has zero operator
p-curvature for all tested primes even though the naive scalar `a^p` acts with
full rank.  The connection `a=z/(1-z)` does the opposite: the naive scalar
vanishes on the p-jet, while the operator curvature is full rank.  So scalar
tests can be wrong both ways.  That is the mathematical reason the LRC proof
route must not stop at "q blocked."

## The proof picture

Try to build a denominator curvature sheaf over the ladder:

```text
plain shells q <= 27
        |
        v
Q27 = {d*m : d|14, m<=27}
        |
        v
carry/fiber maps q -> p*q, q -> d*m, 27 -> 9 -> 3
```

A local section over `q` is not a true/false value.  It is the blocker
hypergraph plus the side coordinates.  A counterexample would need compatible
sections everywhere.  The conjectural theorem is that compatible sections are
finite-monodromy-like: after normalization they are just AP, Vstar, or 2AP.
Every incompatible seam should produce a nonzero-curvature certificate:
a witness denominator, a Bprime(any runner) opening, or an owner-private
deletion.

This also says what the next computation should look like.  We should not only
search for multi-stranger rows blocking `Q27`.  We should compute how the
blocker hypergraphs transform from `q` to the related fiber denominators.
The first useful invariant is probably a transition defect:

```text
defect(q -> q') =
  minimal blocker support at q'
  minus image/extension of minimal blocker support at q.
```

Positive defect is curvature.  Zero defect over too many maps is rigidity.

## Concrete next experiments

1. For each high-pressure row in the HYP-2443 atlas, compute support transport
   along `q -> 7q`, `q -> 2q`, and `m -> 27/m` where those maps make sense.
2. For one-stranger rows, prove by hand that the `q=91` rescue is the fibered
   image of the 13-clock plus the missing `+-10` shell-27 class.
3. For two-stranger probes, stop pairing old bad residues naively; instead
   require low-clock support to remain compatible while shell-27 classes are
   spent.
4. For AP/Vstar/2AP, treat the all-ladder-compatible support ledger as the
   finite monodromy model.  Then prove primitive loose rows cannot mimic that
   ledger without opening an owner-private derivative.

The wild version: LRC14 may not be a search for one clever time.  It may be a
finite-monodromy theorem for a discrete connection whose local p-curvatures are
these denominator blocker ledgers.  A lonely time is nonzero curvature becoming
visible at one finite denominator.
