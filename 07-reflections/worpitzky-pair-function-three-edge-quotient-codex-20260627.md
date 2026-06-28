# Worpitzky Pair Functions And The K3 Compression Test

**Session:** codex-2026-06-27-S274 continuation
**Hypothesis:** HYP-3144

The useful result is not that the three-vertex tournament quotient is
complicated.  It is useful because it is the smallest place where a quotient is
almost legal but still loses the coordinate that the next operation wants.

The exact K3 score-class kernel is:

```text
      to T  to C
from T   2     1
from C   3     0
```

It preserves class transition counts, stationary mass, and the symmetric
Eulerian-looking aggregate.  It forgets the identity of the transitive
source-sink edge whose flip enters the cyclic class.  That is the local version
of the LRC14 rule: a scalar is allowed only after the target predicate factors
through the quotient.

The prompt's pair functions are the cleanest vocabulary for the rule.
`a+b` and `a*b` factor through the unordered-pair quotient.  `a^b` and `b^a`
do not, except at accidental equalities.  The analogous LRC observables are
HYP-3140 fiber-PGF curves, HYP-3141 tip/tail witnesses, HYP-3139 reflection
pages, and HYP-3143/HYP-3145/HYP-3149 packet quotient tables.  Each has to say
whether it is sum/product-like or exponent-like before it is scalarized.

S71 makes the local toy relevant to the current proof frontier.  Score-class
compression is safe through n=4 and first fails at n=5; the k=8 bounded-core
dip turns on at the same information threshold.  Its correction splits into an
even biquadratic face and a larger odd Worpitzky face.  KPS lambda adds the
root-curve version: off-circle variance tracks the coverage gap, so the whole
PGF/root curve is also a sidecar, not decoration.

The post-rebase HYP-3151/HYP-3152 pair sharpens the warning.  HYP-3151 makes
the n=4 cover explicit as `x=a OR c, y=b OR c` and verifies that no affine
`GF(2)^3 -> GF(2)^2` map can replace that nonlinear compression.  HYP-3152
corrects the group story: the two-bit slice may be V4-like, but the full flip
action is a transformation monoid with an absorbing apex arc.  Thus "below the
Abel-Ruffini wall" means degree `<=4` and Galois `<=S4`, plus named sidecars
for destroyed order data, not a blanket V4 quotient.

The later HYP-3153/HYP-3160 pair turns the next target into a sharper split
rather than a larger metaphor.  HYP-3153 reserves the Lee-Yang/Worpitzky/quartic
packet; HYP-3160 says the even k=8 face is degree-two variance, equivalently
total covariance of empty-sector indicators, while entropy is the wrong
principle.  Its S31ai follow-up makes the honest correction: consec maximizes
`Sigma kappa_2`, but it does not maximize the associator `Sigma kappa_3`; the
near `1/7` associativity-defect ratio is a coincidence, and removing the anchor
does not symmetrize the odd coverage distribution.  The odd `-9S3` Worpitzky
face should therefore be treated as non-associative residual sidecar data that
cannot be reconstructed from pairwise covariance alone.
HYP-3200 then verifies this correction exactly on the bounded bank: consec has
`Sigma kappa_3/S3=407891843/2855269200`, not `1/7`, no bank row has exact ratio
`1/7`, while consec is rank `0/3431` for primitive `Sigma kappa_2`.  That turns
the proof target away from an apex-prime scalar and toward a covariance theorem
plus an odd associator sidecar.

HYP-3199 adds the n=4 chart lesson that should travel with this K3 warning.
The fixed-path `a,b,c` cube is only an abundant cover: `S` has five
representatives and class multiplication is ambiguous after quotienting.  The
partial-score `(0,1,1,2)` chart with live `x=a,y=b` is the exact Einheit
section, and the map `x=a OR c, y=b OR c` is a lossy monotone circuit.  So the
same rule becomes sharper: a quotient must record minimality/deletion status,
not just representation count.

The next proof move should be executable, not metaphorical.  Pick one active
frontier quotient and write it as a map `q:X->Y` with observable `f:X->Z`.
Then prove `f` factors through `q`, or record the missing coordinate as
`minority_edge_gate`, `ordered_pair_exponent_sidecar`,
`worpitzky_descent_word`, `fiber_pgf_order_loss_alarm`, or a more specific
named debt.
